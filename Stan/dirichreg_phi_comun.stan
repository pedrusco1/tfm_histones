// -----------------------------------------------------------------------------
// DOCUMENTACIÓN
// Modelo de regresión Dirichlet con precisión común.
// La media composicional se obtiene con softmax aplicado al predictor lineal.
// La categoría 1 actúa como referencia y sus coeficientes se fijan a cero.
// alpha[n] = mu[n] * phi, donde phi es compartido por todas las observaciones.
// prior_only = 1 permite simular desde la prior sin usar la verosimilitud.
// -----------------------------------------------------------------------------

//Modelo regresión Dirichlet
// Precisión (phi) comun

data {
  int <lower=1> N;       //Numero de observaciones
  int <lower=2> K;       //Numero de partes
  int <lower=2> P;       //Numero de columnas matriz de diseño
  matrix[N,K] Y;         //Matriz respuesta
  matrix[N,P] X;         //Matriz de diseño
  int <lower=0, upper=1> prior_only; // 1=Hacer verificacion preditiva previa 
}

parameters {
  matrix[K-1,P] beta_c; // Coeficientes para K-1 categorías
  //real logphi;          // Log-precisión única para todo el modelo
  real <lower=0> phi;
}

transformed parameters {
 // real <lower=0> phi = exp(logphi);
  matrix[K,P] beta;
  simplex[K] mu[N];     // Proporciones medias estimadas
  vector[K] alpha[N];   // Parámetros de concentración final

  // Definición de la categoría de referencia (primera categoría = 0)
  for (i in 1:P) {
    beta[1,i] = 0;
  }
  for (k in 2:K) {
    beta[k, ] = beta_c[k-1, ];
  }

  for (n in 1:N) {
    vector[K] eta;
    for (k in 1:K) {
      eta[k] = X[n, ] * beta[k, ]';
    }
    mu[n] = softmax(eta);   // La media es el softmax del predictor lineal
    alpha[n] = mu[n] * phi; // Concentración = proporción * precisión
  }
}

model {
  to_vector(beta_c) ~ normal(0, 2.5);
  // logphi ~ normal(0, 10);
  phi ~ gamma(2, 0.001);

  // Likelihood SOLO si no estamos en prior predictive
  if (prior_only == 0) {
   for (n in 1:N) {
    Y[n, ] ~ dirichlet(alpha[n]);
   }
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = dirichlet_lpdf(Y[n, ] | alpha[n]);
  }
  matrix[N, K] Y_rep;
  for (n in 1:N) {
    // Se añade el apostrofe (') para transponer el vector a row_vector
    Y_rep[n] = dirichlet_rng(alpha[n])';
  }
}


