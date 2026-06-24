// -----------------------------------------------------------------------------
// DOCUMENTACIÓN
// Modelo de regresión Dirichlet con precisión común y efecto aleatorio de día.
// La media composicional se modela con softmax(eta).
// Para las categorías 2..K, eta incluye efectos fijos X*beta y una desviación
// aleatoria u_day asociada al día experimental.
// Los efectos de día se construyen con parametrización no centrada y suma cero.
// alpha[n] = mu[n] * phi.
// -----------------------------------------------------------------------------

// Modelo de regresión Dirichlet
// Precisión (phi) común + efecto aleatorio de día
//
// Este modelo extiende dirichreg_phi_comun.stan añadiendo un efecto aleatorio
// del día experimental sobre la media composicional. La precisión Dirichlet
// sigue siendo común para todas las observaciones.
//
// Estructura:
//   Y[n] ~ Dirichlet(alpha[n])
//   alpha[n] = phi * mu[n]
//   mu[n] = softmax(eta[n])
//
// donde:
//   eta[n, 1] = 0                                             categoría de referencia
//   eta[n, k] = X[n, ] * beta_c[k-1, ]' + u_day[k-1, day[n]]  k = 2,...,K
//
// El efecto de día es un intercepto aleatorio específico de cada componente
// composicional, excepto la categoría de referencia.

data {
  int<lower=1> N;                       // Número de observaciones
  int<lower=2> K;                       // Número de partes/categorías composicionales
  int<lower=1> P;                       // Número de columnas de la matriz de diseño X
  int<lower=1> D;                       // Número de días experimentales
  array[N] int<lower=1, upper=D> day;   // Día experimental de cada observación

  matrix[N, K] Y;                       // Matriz respuesta composicional; filas positivas y suma 1
  matrix[N, P] X;                       // Matriz de diseño para la media composicional

  int<lower=0, upper=1> prior_only;     // 1 = prior predictive; 0 = usar verosimilitud
}

parameters {
  matrix[K - 1, P] beta_c;              // Coeficientes para K-1 categorías; categoría 1 es referencia
  //real logphi;                          // Log-precisión común
  real<lower=0> phi;
  matrix[K - 1, D] z_day;               // Efectos aleatorios estandarizados de día
  vector<lower=0>[K - 1] sigma_day;     // Desviación típica del efecto de día por categoría
}

transformed parameters {
  //real<lower=0> phi = exp(logphi);      // Precisión común Dirichlet

  matrix[K - 1, D] u_day;               // Efectos aleatorios de día en escala del predictor lineal
  matrix[K, P] beta;                    // Matriz completa de coeficientes, incluyendo referencia

  array[N] simplex[K] mu;               // Media composicional esperada
  array[N] vector[K] alpha;             // Parámetros Dirichlet

// Calcular efectos aleatorios forzando la suma cero desde el vector base (z_day)
  for (k in 1:(K - 1)) {
    vector[D] z_day_raw;
    for (d in 1:D) {
      z_day_raw[d] = z_day[k, d];
    }
    
    // Restamos la media al vector estandarizado primero
    real mean_z = mean(z_day_raw);
    for (d in 1:D) {
      // Ahora el efecto real es estrictamente identificable con sigma_day
      u_day[k, d] = sigma_day[k] * (z_day_raw[d] - mean_z); 
    }
  }

  // Categoría de referencia: coeficientes fijados a cero.
  for (p in 1:P) {
    beta[1, p] = 0;
  }
  for (k in 2:K) {
    beta[k, ] = beta_c[k - 1, ];
  }

  for (n in 1:N) {
    vector[K] eta;

    eta[1] = 0;
    for (k in 2:K) {
      eta[k] = X[n, ] * beta[k, ]' + u_day[k - 1, day[n]];
    }

    mu[n] = softmax(eta);
    alpha[n] = mu[n] * phi;
  }
}

model {
  // Priors de efectos fijos
  to_vector(beta_c) ~ normal(0, 2.5);
  //logphi ~ normal(0, 10);
  phi ~ gamma(2,0.001);

  // Priors del efecto aleatorio de día
  to_vector(z_day) ~ normal(0, 1);
  sigma_day ~ exponential(2);

  // Verosimilitud
  if (prior_only == 0) {
    for (n in 1:N) {
      target += dirichlet_lpdf(to_vector(Y[n, ]') | alpha[n]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  matrix[N, K] Y_rep;

  for (n in 1:N) {
    log_lik[n] = dirichlet_lpdf(to_vector(Y[n, ]') | alpha[n]);
    Y_rep[n] = dirichlet_rng(alpha[n])';
  }
}
