######################################################
# MODELOS LOG-NORMAL V1 (SUMA CERO): LOO-WAIC        #
######################################################

# 1. Configuración de rutas y funciones
source("ScriptsR/40_stan_config.R") # Asegura la disponibilidad de data_dir
source(here::here("Funciones", "fun_loo_waic_stan.R"))

# 2. Cargar el archivo consolidado de los 4 modelos
out_file <- file.path(data_dir, "stan_lognormal_sum0_fits.Rdata")

if (!file.exists(out_file)) {
  stop("No se encuentra el archivo: ", out_file, "\nAsegúrate de ejecutar primero el script de ajuste de los 4 modelos.")
}

cat("Cargando ajustes y diagnósticos desde:", out_file, "\n")
load(out_file)

# 3. Obtener la lista de nombres de los modelos de forma dinámica
# El patrón "^fit_lognor_v1_" capturará exactamente tus 4 nuevos modelos
model_names <- ls(pattern = "^fit_lognor_v1_")

if (length(model_names) == 0) {
  stop("No se encontraron objetos de tipo stanfit con el patrón '^fit_lognor_v1_' en el entorno global.")
}

cat("Modelos detectados para la comparación:\n", 
    paste("  -", model_names, collapse = "\n"), "\n\n")

# 4. Recuperar los objetos stanfit del entorno
list_of_fits <- mget(model_names)

# 5. Calcular LOO y WAIC dinámicamente usando las funciones auxiliares
cat("Calculando LOO para todos los modelos...\n")
loo_results <- lapply(list_of_fits, compute_loo)
names(loo_results) <- model_names  # Renombramos la lista para identificar el output

cat("Calculando WAIC para todos los modelos...\n")
waic_results <- lapply(list_of_fits, compute_waic)
names(waic_results) <- model_names

# 6. Comparación multimodelo (Criterios de información relativos al simplex)
cat("\n======================================================\n")
cat("--- COMPARACIÓN LOO ---")
cat("\n======================================================\n")
comp_loo <- loo::loo_compare(loo_results)
print(comp_loo, simplify = FALSE) # simplify = FALSE muestra los SE individuales

cat("\n======================================================\n")
cat("--- COMPARACIÓN WAIC ---")
cat("\n======================================================\n")
comp_waic <- loo::loo_compare(waic_results)
print(comp_waic, simplify = FALSE)

# 7. Guardar todos los resultados en una lista RDS
resultados_loo_waic <- list(
  modelos = model_names,
  loo = loo_results,
  waic = waic_results,
  comparacion_loo = comp_loo,
  comparacion_waic = comp_waic
)

rds_file <- file.path(data_dir, "loo_waic_lognormal_sum0_resultados.rds")

saveRDS(resultados_loo_waic, file = rds_file)

cat("\nResultados LOO/WAIC guardados en:\n", rds_file, "\n")
