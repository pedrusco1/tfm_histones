// -----------------------------------------------------------------------------
// Modelo: epistasia_kappa_comun
// Propósito:
//   Ajusta un modelo Dirichlet para composiciones Y con dos condiciones
//   experimentales. Cada condición tiene su propio campo energético phi y su
//   propia suavidad epistática a, pero ambas comparten una concentración
//   Dirichlet global kappa.
//
// Estructura estadística:
//   phi_c  = campo energético centrado para la condición c
//   Q_c    = softmax(-phi_c), media composicional esperada por condición
//   alpha  = kappa * Q_c + 1e-6
//   Y[n]   ~ Dirichlet(alpha), si prior_only == 0
//
// Prior epistático:
//   target += 0.5 * rank_Delta * log(a[c])
//             - (a[c] / (2*s)) * phi_c' Delta phi_c
//   Este término penaliza campos energéticos poco suaves según Delta.
//
// Uso de prior_only:
//   prior_only = 1: ignora la verosimilitud y permite chequeos predictivos
//                 desde el prior.
//   prior_only = 0: ajusta el modelo usando los datos observados.
// -----------------------------------------------------------------------------
//Modelo epistasia kappa común

  data {
    int<lower=1> K;              // nº de PTM states
    int<lower=1> N;              // nº de observaciones
    matrix[N, K] Y;              // composiciones observadas
    array[N] int<lower=1, upper=2> cond; // condición de cada observación: 1=async, 2=mitosis

    matrix[K, K] Delta;          // matriz de Laplaciano (mínima epistasia)
    int rank_Delta;              // rango efectivo de Delta
    real s;                      // escala (como en tu modelo)
    int<lower=0, upper=1> prior_only;// 1 = prior predictive; 0 = usar verosimilitud
  }

  parameters {
    // Campo por condición
    matrix[K, 2] phi_raw;        // columnas: condición 1 y 2

    // Parámetros de suavidad por condición
    vector<lower=0>[2] a;

    // Concentración Dirichlet por condición
    real log_kappa;
  } 

  transformed parameters {
    real<lower=0> kappa = exp(log_kappa);

    // Campos centrados
    matrix[K, 2] phi;
    // Distribuciones de equilibrio por condición
    matrix[K, 2] Q;

    for (c in 1:2) {
        vector[K] phi_c_raw = phi_raw[, c];
        vector[K] phi_c = phi_c_raw - mean(phi_c_raw);
        phi[, c] = phi_c;
        Q[, c] = softmax(-phi_c);
      }
  }

  model {
    // Priors
    to_vector(phi_raw) ~ normal(0, 10);
    a ~ gamma(2, 0.1);
    log_kappa ~ normal(0, 3);

    // Prior de campo (mínima epistasia) por condición
    for (c in 1:2) {
      vector[K] phi_c = phi[, c];
      target += 0.5 * rank_Delta * log(a[c])
              - (a[c] / (2 * s)) * quad_form(Delta, phi_c);
    }

    // Likelihood Dirichlet
    if (prior_only == 0) {
    	for (n in 1:N) {
      	int c = cond[n];
      	vector[K] alpha = kappa * Q[, c] + 1e-6;
      	target += dirichlet_lpdf(to_vector(Y[n, ]') | alpha);
    	}
    }
  }

  generated quantities {
    matrix[N, K] Y_rep;
    vector[N] log_lik;

    for (n in 1:N) {
      int c = cond[n];
      vector[K] alpha = kappa * Q[, c] + 1e-6;
      Y_rep[n] = to_row_vector(dirichlet_rng(alpha));
      log_lik[n] = dirichlet_lpdf(to_vector(Y[n, ]') | alpha);
    }
  }
  
