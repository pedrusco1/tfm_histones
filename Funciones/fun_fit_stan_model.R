###########################################
#    Funcion ajustar mododelo stan        #
###########################################

fit_stan_model <- function(stan_file,
                           data_list,
                           iter = 8000,
                           warmup = 4000,
                           chains = 4,
                           adapt_delta = 0.99,
                           max_treedepth = 15) {
  
  rstan::stan(
    file = stan_file,
    data = data_list,
    iter = iter,
    seed = 123,
    warmup = warmup,
    chains = chains,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
}
