############################################################
# LOO Y WAIC: 4 MODELOS LOG-NORMAL V1 REVISADOS
############################################################

library(here)
source(here::here("ScriptsR", "40_stan_config.R"))
source(here::here("Funciones", "fun_loo_waic_stan.R"))

input_file <- file.path(data_dir, "stan_lognormal_sum0_fits_revisados.Rdata")
if (!file.exists(input_file)) stop("No se encuentra: ", input_file)

# Cargar en entorno aislado evita mezclar modelos antiguos y revisados.
env_revisado <- new.env(parent = globalenv())
load(input_file, envir = env_revisado)

model_names <- sort(ls(
  envir = env_revisado,
  pattern = "^fit_lognor_v1_.*_revisado$"
))

if (length(model_names) != 4) {
  stop("Se esperaban 4 modelos revisados y se encontraron ", length(model_names),
       ": ", paste(model_names, collapse = ", "))
}

list_of_fits <- mget(model_names, envir = env_revisado)
loo_results_revisados <- lapply(list_of_fits, compute_loo)
waic_results_revisados <- lapply(list_of_fits, compute_waic)

comp_loo_revisados <- loo::loo_compare(loo_results_revisados)
comp_waic_revisados <- loo::loo_compare(waic_results_revisados)

print(comp_loo_revisados, simplify = FALSE)
print(comp_waic_revisados, simplify = FALSE)

resultados_loo_waic_lognormal_revisados <- list(
  modelos = model_names,
  loo = loo_results_revisados,
  waic = waic_results_revisados,
  comparacion_loo = comp_loo_revisados,
  comparacion_waic = comp_waic_revisados
)

rds_file <- file.path(data_dir, "loo_waic_lognormal_sum0_revisados.rds")
saveRDS(resultados_loo_waic_lognormal_revisados, rds_file)
message("Resultados revisados guardados en: ", rds_file)
