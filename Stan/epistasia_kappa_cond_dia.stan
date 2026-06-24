// Modelo de epistasia con kappa por condición + efecto aleatorio de día
// Versión comparable a Dirichlet en el CRITERIO DE CENTRADO de los efectos aleatorios
//
// Nota importante:
//   A diferencia de los modelos Dirichlet de regresión, aquí NO se fija la categoría 1
//   como referencia. En epistasia el campo energético phi se mantiene para las K
//   componentes y la identificabilidad frente a traslaciones constantes se controla
//   centrando los vectores del campo y del efecto de día.
//
// Criterio de centrado aplicado al efecto aleatorio de día:
//   Para cada componente k, los efectos de día se centran a través de los D días:
//       day_eff[k, d] = sigma_day[k] * (z_day[k, d] - mean_d z_day[k, d])
//   Por tanto, para cada componente k:
//       sum_d day_eff[k, d] = 0
//
// Este es el análogo directo del criterio usado en los modelos Dirichlet con día,
// donde cada componente/categoría tiene desviaciones de día centradas sobre días.

data {
  int<lower=1> K;                         // Número de estados PTM / componentes
  int<lower=1> N;                         // Número de observaciones
  int<lower=1> D;                         // Número de días experimentales

  matrix[N, K] Y;                         // Composiciones observadas
  array[N] int<lower=1, upper=2> cond;    // Condición: 1=async, 2=mitosis
  array[N] int<lower=1, upper=D> day;     // Día experimental de cada observación

  matrix[K, K] Delta;                     // Operador de mínima epistasia
  int<lower=0> rank_Delta;                // Rango efectivo de Delta
  real<lower=0> s;                        // Escala del prior epistático

  int<lower=0, upper=1> prior_only;       // 1 = prior predictive; 0 = usar verosimilitud
}

parameters {
  matrix[K, 2] phi_raw;                   // Campo energético bruto por condición
  vector<lower=0>[2] a;                   // Suavidad epistática por condición
  vector[2] log_kappa;                    // Log-concentración Dirichlet por condición

  matrix[K, D] z_day;                     // Efectos aleatorios estandarizados de día
  vector<lower=0>[K] sigma_day;           // Desviación típica del efecto de día por componente
}

transformed parameters {
  vector<lower=0>[2] kappa = exp(log_kappa);

  matrix[K, 2] phi;                       // Campos centrados por condición
  matrix[K, 2] Q_cond;                    // Distribuciones esperadas sin efecto de día
  matrix[K, D] day_eff;                   // Efectos de día centrados sobre días para cada componente

  // Centrado del campo por condición: softmax(-phi) es invariante a constantes aditivas.
  for (c in 1:2) {
    vector[K] phi_c = phi_raw[, c] - mean(phi_raw[, c]);
    phi[, c] = phi_c;
    Q_cond[, c] = softmax(-phi_c);
  }

  // Centrado comparable a Dirichlet: para cada componente k, se centra a través de los días.
  // No se usa categoría de referencia; se conservan las K componentes epistáticas.
  for (k in 1:K) {
    vector[D] z_k;
    for (d in 1:D) {
      z_k[d] = z_day[k, d];
    }

    for (d in 1:D) {
      day_eff[k, d] = sigma_day[k] * (z_k[d] - mean(z_k));
    }
  }
}

model {
  // Priors originales
  to_vector(phi_raw) ~ normal(0, 10);
  a ~ gamma(2, 0.1);
  log_kappa ~ normal(0, 3);

  // Priors del efecto aleatorio de día, análogos a Dirichlet pero para K componentes.
  to_vector(z_day) ~ normal(0, 1);
  sigma_day ~ exponential(2);

  // Prior epistático de mínima epistasia por condición
  for (c in 1:2) {
    vector[K] phi_c = phi[, c];
    target += 0.5 * rank_Delta * log(a[c])
              - (a[c] / (2 * s)) * quad_form(Delta, phi_c);
  }

  // Verosimilitud Dirichlet con media modificada por día
  if (prior_only == 0) {
    for (n in 1:N) {
      int c = cond[n];
      vector[K] phi_n = phi[, c] + day_eff[, day[n]];
      vector[K] Q_n = softmax(-phi_n);
      vector[K] alpha = kappa[c] * Q_n + 1e-6;

      target += dirichlet_lpdf(to_vector(Y[n, ]') | alpha);
    }
  }
}

generated quantities {
  matrix[N, K] Y_rep;
  vector[N] log_lik;

  for (n in 1:N) {
    int c = cond[n];
    vector[K] phi_n = phi[, c] + day_eff[, day[n]];
    vector[K] Q_n = softmax(-phi_n);
    vector[K] alpha = kappa[c] * Q_n + 1e-6;

    Y_rep[n] = to_row_vector(dirichlet_rng(alpha));
    log_lik[n] = dirichlet_lpdf(to_vector(Y[n, ]') | alpha);
  }
}
