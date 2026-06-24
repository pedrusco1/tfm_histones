// Modelo logístico-normal en coordenadas ILR
// Precisión residual: Q = tau * L_shape + diag(delta)
// Coeficientes CLR con suma cero inducida por una base ILR ortonormal
// Variante flexible: precisiones diagonales con contracción jerárquica.

functions {
  real multi_normal_prec_chol_lpdf(vector z, vector mu, matrix LQ) {
    int d = rows(z);
    vector[d] v = LQ' * (z - mu);
    return -0.5 * d * log(2 * pi())
           + sum(log(diagonal(LQ)))
           - 0.5 * dot_self(v);
  }
}

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> K_ilr;
  int<lower=1> P;
  matrix[N, K_ilr] Z;
  matrix[N, P] X;
  matrix[K_ilr, K_ilr] L_shape;
  matrix[K, K_ilr] V;
  matrix[N, K] Y;
}

transformed data {
  if (K_ilr != K - 1) reject("K_ilr debe ser igual a K - 1");

  for (i in 1:N) {
    real suma_y = 0;
    for (k in 1:K) {
      if (Y[i, k] <= 0) reject("Todas las componentes de Y deben ser positivas");
      suma_y += Y[i, k];
    }
    if (fabs(suma_y - 1.0) > 1e-6)
      reject("Cada fila de Y debe sumar 1; fila ", i, " suma ", suma_y);
  }
}

parameters {
  matrix[K_ilr, P] B_ilr;
  real log_tau;
  real mu_log_delta;
  real<lower=0> sigma_log_delta;
  vector[K_ilr] z_log_delta;
}

transformed parameters {
  real<lower=0> tau;
  vector[K_ilr] log_delta;
  vector<lower=0>[K_ilr] delta;
  matrix[K_ilr, K_ilr] Q;
  matrix[K_ilr, K_ilr] LQ;
  matrix[N, K_ilr] MU;
  matrix[K, P] B;

  B = V * B_ilr;

  tau = exp(log_tau);
  log_delta = mu_log_delta + sigma_log_delta * z_log_delta;
  delta = exp(log_delta);
  Q = tau * L_shape + diag_matrix(delta);
  for (k in 1:K_ilr) Q[k, k] += 1e-9;
  LQ = cholesky_decompose(Q);

  for (i in 1:N) {
    MU[i] = (B_ilr * to_vector(X[i]'))';
  }
}

model {
  // En escala ILR, normal(0, 3) es un prior amplio para interceptos y contrastes.
  to_vector(B_ilr) ~ normal(0, 3);
  // Media global amplia y desviaciones regularizadas hacia una precisión común.
  log_tau ~ normal(0, 1.5);
  mu_log_delta ~ normal(0, 1.5);
  sigma_log_delta ~ normal(0, 0.5);
  z_log_delta ~ std_normal();

  for (i in 1:N) {
    target += multi_normal_prec_chol_lpdf(
      to_vector(Z[i]') | to_vector(MU[i]'), LQ
    );
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] log_lik_ilr;
  matrix[N, K_ilr] Z_rep;
  matrix[N, K] Y_rep;
  matrix[K_ilr, K_ilr] Sigma;
  matrix[K_ilr, K_ilr] L_Sigma;

  Sigma = inverse_spd(Q);
  L_Sigma = cholesky_decompose(Sigma);

  for (i in 1:N) {
    vector[K_ilr] mu_i = to_vector(MU[i]');
    vector[K_ilr] z_rep_i;
    vector[K] clr_rep;
    real sum_log_y = 0;

    log_lik_ilr[i] = multi_normal_prec_chol_lpdf(
      to_vector(Z[i]') | mu_i, LQ
    );

    for (k in 1:K) sum_log_y += log(Y[i, k]);
    log_lik[i] = log_lik_ilr[i] - sum_log_y - 0.5 * log(K);

    z_rep_i = multi_normal_cholesky_rng(mu_i, L_Sigma);
    Z_rep[i] = z_rep_i';
    clr_rep = V * z_rep_i;
    Y_rep[i] = softmax(clr_rep)';
  }
}
