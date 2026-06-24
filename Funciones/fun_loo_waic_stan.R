library(rstan)
library(loo)

# Tu función original optimizada para asegurar compatibilidad
compute_loo <- function(fit) {
  log_lik <- rstan::extract(fit, "log_lik")[[1]]
  # relative_eff ayuda a corregir las cadenas en el muestreo de HMC
  rel_eff <- loo::relative_eff(exp(log_lik), chain_id = rep(1:fit@sim$chains, each = fit@sim$iter - fit@sim$warmup))
  loo::loo(log_lik, r_eff = rel_eff)
}

# Nueva función para calcular el WAIC
compute_waic <- function(fit) {
  log_lik <- rstan::extract(fit, "log_lik")[[1]]
  loo::waic(log_lik)
}
