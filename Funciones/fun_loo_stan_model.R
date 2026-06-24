#############################
#     Funcion hacer loo     #
#############################

compute_loo <- function(fit) {
  
  log_lik <- rstan::extract(fit, "log_lik")[[1]]
  loo::loo(log_lik)
}
