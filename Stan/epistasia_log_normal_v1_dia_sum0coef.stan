//
// Modelo de distribucion log-normal multivariante
// Sobre corrdenadas ILR
// Con regularizacion epistatica
// Coeficientes expresados también en CLR con suma cero por covariable
// Incluye efecto aleatorio de día centrado en coordenadas ILR

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

  // Efecto aleatorio de día
  int<lower=1> D;                         // número de días experimentales
  array[N] int<lower=1, upper=D> day;     // índice de día para cada muestra
}

parameters {
  matrix[K_ilr, P] B_ilr;                 // coeficientes libres en coordenadas ILR
  real<lower=0> tau;                      // intensidad global de regularización epistática
  vector[K_ilr] log_delta;                // escala log de la precisión residual diagonal

  // Efectos aleatorios de día en coordenadas ILR
  // day_raw[d, j] es el efecto no centrado del día d sobre la coordenada ILR j
  matrix[D, K_ilr] day_raw;
  vector<lower=0>[K_ilr] sigma_day;       // escala del efecto de día por coordenada ILR
}

transformed parameters {
  vector<lower=0>[K_ilr] delta;
  matrix[K_ilr, K_ilr] Q;
  matrix[K_ilr, K_ilr] LQ;
  matrix[N, K_ilr] MU;
  matrix[K, P] B;                         // coeficientes en CLR; cada columna suma cero
  matrix[D, K_ilr] day_eff;

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

  // Efectos de día centrados:
  // Para cada coordenada ILR, los efectos de los D días suman cero.
  // Esto evita que el efecto aleatorio de día absorba el intercepto global
  // y hace el modelo comparable a la parametrización usada en Dirichlet.
  for (j in 1:K_ilr) {
    vector[D] raw_j;
    raw_j = col(day_raw, j);
    for (d in 1:D) {
      day_eff[d, j] = sigma_day[j] * (raw_j[d] - mean(raw_j));
    }
  }

  // Medias por muestra: efecto fijo + efecto aleatorio del día experimental
  for (i in 1:N) {
    MU[i] = (B_ilr * to_vector(X[i]'))' + day_eff[day[i]];
  }
}

model {
  // Priors débilmente informativos
  to_vector(B_ilr) ~ normal(0, 2);
  tau ~ gamma(2, 1);
  log_delta ~ normal(0, 1);

  // Priors del efecto aleatorio de día
  to_vector(day_raw) ~ normal(0, 1);
  sigma_day ~ normal(0, 0.5);

  // Likelihood en parametrización por precisión
  for (i in 1:N) {
    target += multi_normal_prec_chol_lpdf(
      to_vector(Z[i]') | to_vector(MU[i]'), LQ
    );
  }
}

generated quantities {
  // log_lik EN EL ESPACIO DEL SIMPLEX (con jacobiano ILR -> simplex)
  // Permite comparacion directa con modelos de verosimilitud Dirichlet
  // via loo::loo() / loo::waic() en R.
  //
  // Formula: log p(y_i) = log p_Z(z_i) + sum_k log(y_{i,k})
  //
  // El termino sum_k log(y_{i,k}) es el log-jacobiano de la transformacion
  // ILR con base ortonormal: al pasar de densidad en R^{K-1} a densidad
  // en el simplex S^K, el factor de cambio de volumen es:
  //   |J_{ILR -> simplex}|^{-1}  =>  log|J| = sum_k log(y_k)
  vector[N] log_lik;

  // log_lik solo en espacio ILR (para diagnostico interno / comparacion
  // entre modelos que usen la misma parametrizacion ILR)
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
    real log_jacobian;

    mu_i = to_vector(MU[i]');

    // -------------------------------------------------------
    // 1. log-likelihood en espacio ILR (parametrizacion interna)
    // -------------------------------------------------------
    log_lik_ilr_i = multi_normal_cholesky_lpdf(
      to_vector(Z[i]') | mu_i, L_Sigma
    );
    log_lik_ilr[i] = log_lik_ilr_i;

    // -------------------------------------------------------
    // 2. log-jacobiano de ILR -> simplex
    //    log|J| = sum_{k=1}^{K} log(y_{i,k})
    //    (valido para cualquier base ILR ortonormal, incluyendo
    //     la base de Egozcue implementada via la matriz V)
    // -------------------------------------------------------
    log_jacobian = 0;
    for (k in 1:K) {
      log_jacobian += log(Y[i, k]);
    }

    // -------------------------------------------------------
    // 3. log-likelihood en el simplex = log p_Z(z_i) + log|J|
    //    Esta es la cantidad comparable con una Dirichlet log_lik
    // -------------------------------------------------------
    log_lik[i] = log_lik_ilr_i + log_jacobian;

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
