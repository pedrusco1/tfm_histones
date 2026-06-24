//
// Modelo de distribucion log-normal multivariante
// Sobre coordenadas ILR
// Con precisión residual estructurada por una matriz epistática
// Coeficientes expresados también en CLR con suma cero por covariable

functions {

  // Log-densidad de una normal multivariante parametrizada por
  // la Cholesky de la matriz de precisión Q = LQ * LQ'
  real multi_normal_prec_chol_lpdf(vector z, vector mu, matrix LQ) {
    int d = rows(z);
    vector[d] v;

    v = LQ' * (z - mu);

    return -0.5 * d * log(2 * pi())
           + sum(log(diagonal(LQ)))
           - 0.5 * dot_self(v);
  }

}

data {
  int<lower=1> N;                         // número de muestras
  int<lower=2> K;                         // número de componentes composicionales
  int<lower=1> K_ilr;                     // K - 1
  int<lower=1> P;                         // número de covariables

  matrix[N, K_ilr] Z;                     // datos observados en coordenadas ILR
  matrix[N, P] X;                         // matriz de diseño

  matrix[K_ilr, K_ilr] L_shape;           // Laplaciano ILR reescalado (solo forma)
  matrix[K, K_ilr] V;                     // base ILR: CLR = V * ILR
  matrix[N,K] Y;                          // Datos observados en el simplex
}

parameters {
  matrix[K_ilr, P] B_ilr;                 // coeficientes libres en coordenadas ILR
  real<lower=0> tau;                      // intensidad global de regularización epistática
  vector[K_ilr] log_delta;                // escala log de la precisión residual diagonal
}

transformed parameters {
  vector<lower=0>[K_ilr] delta;
  matrix[K_ilr, K_ilr] Q;
  matrix[K_ilr, K_ilr] LQ;
  matrix[N, K_ilr] MU;
  matrix[K, P] B;                         // coeficientes en CLR; cada columna suma cero

  delta = exp(log_delta);

  // Proyección de los coeficientes ILR al espacio CLR.
  // Como las columnas de V pertenecen al subespacio CLR, cada columna de B
  // satisface sum(B[, p]) = 0.
  B = V * B_ilr;

  // Matriz de precisión total
  Q = tau * L_shape + diag_matrix(delta);
  for (k in 1:K_ilr) {
    Q[k, k] = Q[k, k] + 1e-9; // Pequeño empuje para estabilidad
  }
  
  // Cholesky de la precisión
  LQ = cholesky_decompose(Q);

  // Medias por muestra
  for (i in 1:N) {
    MU[i] = (B_ilr * to_vector(X[i]'))';
  }
}

model {
  // Priors débilmente informativos
  to_vector(B_ilr) ~ normal(0, 2);
  tau ~ gamma(2, 1);
  log_delta ~ normal(0, 1);

  // Likelihood en parametrización por precisión
  for (i in 1:N) {
    target += multi_normal_prec_chol_lpdf(
      to_vector(Z[i]') | to_vector(MU[i]'), LQ
    );
  }
}

generated quantities {
  // -------------------------------------------------------------------------
  // LOG-VEROSIMILITUD EN EL SIMPLEX
  // -------------------------------------------------------------------------
  // El modelo se ajusta en coordenadas ILR:
  //
  //   Z_i = ilr(Y_i),        Z_i ~ Normal(MU_i, Sigma)
  //
  // Esta es la formulación natural de Aitchison/Mateu-Figueras: una normal
  // multivariante en coordenadas ortonormales del simplex.
  //
  // Para LOO/WAIC comparables con modelos definidos directamente en el simplex
  // respecto a la medida de Lebesgue, como Dirichlet, hay que escribir la
  // densidad inducida para Y_i. Para una base ILR ortonormal:
  //
  //   p_Y(y_i) = p_Z(ilr(y_i)) / (sqrt(K) * prod_k y_{i,k})
  //
  // Por tanto:
  //
  //   log p_Y(y_i) = log p_Z(ilr(y_i))
  //                  - sum_k log(y_{i,k})
  //                  - 0.5 * log(K)
  //
  // Nota: sum_log_y contiene sum_k log(y_{i,k}); no es el log-jacobiano
  // completo. El signo correcto en log_lik es negativo.
  // -------------------------------------------------------------------------
  vector[N] log_lik;

  // log_lik_ilr: densidad normal de Z_i en coordenadas ILR.
  // Se guarda para diagnóstico interno; no incluye el cambio de medida al simplex.
  vector[N] log_lik_ilr;

  matrix[N, K_ilr] Z_rep;
  matrix[N, K] Y_rep;

  // Solo en generated quantities: construir Sigma para simulación PPC
  matrix[K_ilr, K_ilr] Sigma;
  matrix[K_ilr, K_ilr] L_Sigma;

  Sigma = inverse_spd(Q);
  L_Sigma = cholesky_decompose(Sigma);

  for (i in 1:N) {
    vector[K_ilr] mu_i;
    vector[K_ilr] z_rep_i;
    vector[K] clr_rep;
    vector[K] y_rep;
    real log_lik_ilr_i;
    real sum_log_y;

    mu_i = to_vector(MU[i]');

    // -------------------------------------------------------
    // 1. log-likelihood en espacio ILR (parametrizacion interna)
    // -------------------------------------------------------
    log_lik_ilr_i = multi_normal_cholesky_lpdf(
      to_vector(Z[i]') | mu_i, L_Sigma
    );
    log_lik_ilr[i] = log_lik_ilr_i;

    // -------------------------------------------------------
    // 2. Suma de log-componentes de la composición observada
    //    sum_log_y = sum_k log(Y[i, k])
    //    En la densidad logística-normal respecto a Lebesgue entra
    //    con signo negativo y con la constante -0.5 * log(K).
    // -------------------------------------------------------
    sum_log_y = 0;
    for (k in 1:K) {
      sum_log_y += log(Y[i, k]);
    }

    // -------------------------------------------------------
    // 3. log-likelihood en el simplex = log p_Z(z_i) - sum_log_y - 0.5 * log(K)
    //    Esta es la cantidad comparable con una Dirichlet log_lik
    // -------------------------------------------------------
    log_lik[i] = log_lik_ilr_i - sum_log_y - 0.5 * log(K);

    // -------------------------------------------------------
    // 4. Replica posterior predictiva en espacio ILR
    // -------------------------------------------------------
    z_rep_i = multi_normal_cholesky_rng(mu_i, L_Sigma);
    Z_rep[i] = z_rep_i';

    // ILR^{-1}: ILR -> CLR -> simplex
    clr_rep = V * z_rep_i;

    for (k in 1:K) {
      y_rep[k] = exp(clr_rep[k]);
    }
    y_rep = y_rep / sum(y_rep);

    Y_rep[i] = y_rep';
  }
}

