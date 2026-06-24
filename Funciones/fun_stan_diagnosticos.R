#############################################
#      Diagnosticos muestreos de Stan       #
#############################################

compute_diagnostics <- function(fit) {
  
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  divergences <- sum(sapply(sampler_params,
                            function(x) sum(x[, "divergent__"])))
  
  max_rhat <- max(summary(fit)$summary[, "Rhat"], na.rm = T)
  min_neff <- min(summary(fit)$summary[, "n_eff"], na.rm = T)
  
  list(
    divergences = divergences,
    max_rhat = max_rhat,
    min_neff = min_neff
  )
}
