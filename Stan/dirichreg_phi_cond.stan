// -----------------------------------------------------------------------------
// DOCUMENTACIÓN
// Modelo de regresión Dirichlet con precisión dependiente de la condición.
// La media composicional se obtiene mediante softmax del predictor lineal.
// La categoría 1 es la referencia y no tiene coeficientes libres.
// alpha[n] = mu[n] * phi[condition[n]], permitiendo distinta dispersión por
// condición experimental.
// -----------------------------------------------------------------------------

// Modelo regresión de Dirichlet
// Precisión (phi) por condición experimental

data {
  int <lower=1> N;             // Número de observaciones
  int <lower=2> K;             // Número de categorías
  int <lower=1> P;             // Número de predictores para mu
  int <lower=1> J;             // Número de condiciones (ej. 2)
  int <lower=1, upper=J> condition[N]; // Índice de condición para cada observación
  matrix[N, K] Y;              // Datos composicionales
  matrix[N, P] X;              // Matriz de diseño
  int <lower=0, upper=1> prior_only;
}

parameters {
  matrix[K-1, P] beta_c;       // Coeficientes para las categorías (menos la de referencia)
  //vector[J] gamma;             // NUEVO: Parámetro de precisión por cada condición
  vector<lower=0>[J] phi;      //Vector de precisión por condición (tamaño J)
}

transformed parameters {
  matrix[K, P] beta;
  simplex[K] mu[N];            // Proporciones medias
  //vector<lower=0>[J] phi; // Vector de precisión por condición (tamaño J)
  vector[K] alpha[N];          // Parámetros de la distribución Dirichlet

  // Definir phi para cada una de las J condiciones
  //for (j in 1:J) {
    //phi[j] = exp(gamma[j]); 
  //}

  // Construcción de la matriz beta con la primera fila como referencia (0)
  for (i in 1:P) beta[1, i] = 0; 
  for (k in 2:K) beta[k, ] = beta_c[k-1, ];

  for (n in 1:N) {
    vector[K] eta;
    for (k in 1:K) {
      eta[k] = X[n, ] * beta[k, ]';
    }
    mu[n] = softmax(eta);
    
    // Multiplicamos mu por el phi correspondiente a la condición de la observación n
    alpha[n] = mu[n] * phi[condition[n]];
  }
}

model {
  // Priors
  to_vector(beta_c) ~ normal(0, 2.5); 
  //gamma ~ normal(0, 10);
  phi ~ gamma(2,0.001);

  // Likelihood
  if (prior_only == 0) {
    for (n in 1:N) {
      Y[n, ] ~ dirichlet(alpha[n]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  matrix[N, K] Y_rep;
  
  for (n in 1:N) {
    log_lik[n] = dirichlet_lpdf(Y[n, ] | alpha[n]);
    // Generar una observación simulada (Posterior Predictive Check)
    Y_rep[n] = dirichlet_rng(alpha[n])';
  }
}

