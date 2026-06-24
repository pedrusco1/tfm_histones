// -----------------------------------------------------------------------------
// DOCUMENTACIÓN
// Modelo de regresión Dirichlet con precisión por condición y efecto aleatorio
// de día sobre la media composicional.
// Combina phi[condition[n]] para la dispersión con u_day[k, day[n]] para capturar
// variación diaria en el predictor lineal de las categorías 2..K.
// ADVERTENCIA: el bloque transformed parameters contiene una sección duplicada
// que recalcula beta, mu y alpha. La segunda versión sobrescribe la primera.
// -----------------------------------------------------------------------------

// Modelo de regresión Dirichlet
// Precisión (phi) por condición + efecto aleatorio de día
//
// Este modelo extiende dirichreg_phi_cond.stan añadiendo un efecto aleatorio
// del día experimental sobre la media composicional. La precisión Dirichlet
// sigue dependiendo de la condición experimental.
//
// Estructura:
//   Y[n] ~ Dirichlet(alpha[n])
//   alpha[n] = phi[condition[n]] * mu[n]
//   mu[n] = softmax(eta[n])
//
// donde:
//   eta[n, 1] = 0
//   eta[n, k] = X[n, ] * beta_c[k-1, ]' + u_day[k-1, day[n]]

data {
  int<lower=1> N;                         // Número de observaciones
  int<lower=2> K;                         // Número de partes/categorías composicionales
  int<lower=1> P;                         // Número de predictores para la media
  int<lower=1> J;                         // Número de condiciones experimentales
  int<lower=1> D;                         // Número de días experimentales

  array[N] int<lower=1, upper=J> condition; // Condición de cada observación
  array[N] int<lower=1, upper=D> day;       // Día experimental de cada observación

  matrix[N, K] Y;                         // Datos composicionales
  matrix[N, P] X;                         // Matriz de diseño

  int<lower=0, upper=1> prior_only;       // 1 = prior predictive; 0 = usar verosimilitud
}

parameters {
  matrix[K - 1, P] beta_c;                // Coeficientes para K-1 categorías; categoría 1 referencia
  //vector[J] gamma;                        // Log-precisión por condición
  vector<lower=0>[J] phi;
  matrix[K - 1, D] z_day;                 // Efectos aleatorios estandarizados de día
  vector<lower=0>[K - 1] sigma_day;       // Desviación típica del efecto de día por categoría
}

transformed parameters {
  matrix[K, P] beta;
  matrix[K - 1, D] u_day;                 // Efectos aleatorios de día CENTRADOS
  array[N] simplex[K] mu;                 // Media composicional
  array[N] vector[K] alpha;               // Parámetros Dirichlet

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

  // 2. Definir matriz beta con la categoría de referencia
  for (p in 1:P) {
    beta[1, p] = 0;
  }
  for (k in 2:K) {
    beta[k, ] = beta_c[k - 1, ];
  }

  // 3. Bucle de observaciones limpio y ultra-eficiente para los gradientes
  for (n in 1:N) {
    vector[K] eta;
    eta[1] = 0;
    for (k in 2:K) {
      // u_day ya viene estrictamente centrado desde arriba
      eta[k] = X[n, ] * beta[k, ]' + u_day[k - 1, day[n]];
    }

    mu[n] = softmax(eta);
    alpha[n] = mu[n] * phi[condition[n]];
  }
}

model {
  // Priors de efectos fijos y precisión
  //to_vector(beta_c) ~ normal(0, 10);
  to_vector(beta_c) ~ normal(0, 2.5);
  //gamma ~ normal(0, 10);
  phi ~ gamma(2, 0.001);

  // Priors del efecto aleatorio de día
  to_vector(z_day) ~ normal(0, 1);
  //sigma_day ~ normal(0, 1);
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

