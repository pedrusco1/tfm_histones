###########################################
#    Función ajustar modelo cmdstanr      #
###########################################

fit_cmdstan_model <- function(stan_file,
                              data_list,
                              iter_sampling = 4000,
                              iter_warmup = 4000,
                              chains = 4,
                              parallel_chains = 4,
                              adapt_delta = 0.99,
                              max_treedepth = 15) {
  
  # 1. Verificación de seguridad
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("El paquete 'cmdstanr' es necesario. Instálalo con: 
         remotes::install_github('stan-dev/cmdstanr')")
  }
  
  # 2. Compilar el modelo
  # Esto crea un objeto ejecutable a partir del archivo .stan
  mod <- cmdstanr::cmdstan_model(stan_file)
  
  # 3. Ajustar el modelo (Muestreo)
  fit <- mod$sample(
    data = data_list,
    seed = 123,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  # 4. Conversión a objeto 'stanfit' (clásico de rstan)
  # Leemos los archivos CSV de salida que generó cmdstanr
  #fit_rstan <- rstan::read_stan_csv(fit$output_files())
  
  return(fit)
}
