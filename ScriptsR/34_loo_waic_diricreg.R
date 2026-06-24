##############################
# MODELOS DIRICLET: LOO-WAIC #
##############################

source(here::here("Funciones", "fun_loo_waic_stan.R"))
load(here::here("DatosProcesados", "stan_dirichlet_fits.Rdata"))
# 1. Obtener la lista de nombres de los modelos
model_names <- ls(pattern = "^fit_dirich_")

if (length(model_names) == 0) {
  stop("No se encontraron objetos con nombres que empiecen por 'fit_dirich_'.")
}

# 2. Calcular LOO y WAIC para todos los modelos dinámicamente
# Usamos mget() para recuperar el objeto stanfit a partir de su nombre de texto
list_of_fits <- mget(model_names)

cat("Calculando LOO para todos los modelos...\n")
loo_results <- lapply(list_of_fits, compute_loo)
names(loo_results) <- model_names

cat("Calculando WAIC para todos los modelos...\n")
waic_results <- lapply(list_of_fits, compute_waic)
names(waic_results) <- model_names

# 3. Comparación multimodelo
# El paquete 'loo' usa loo_compare() tanto para objetos loo como waic
cat("\n--- COMPARACIÓN LOO ---\n")
comp_loo <- loo::loo_compare(loo_results)
print(comp_loo, simplify = FALSE)

cat("\n--- COMPARACIÓN WAIC ---\n")
comp_waic <- loo::loo_compare(waic_results)
print(comp_waic, simplify = FALSE)

# 4. Guardar una lista con dos tablas: una para LOO y otra para WAIC
loo_waic_dirichreg <- list(
  loo = as.data.frame(comp_loo),
  waic = as.data.frame(comp_waic)
)

# Añadir el nombre del modelo como columna para que no se pierda al exportar/leer
loo_waic_dirichreg$loo$modelo <- rownames(loo_waic_dirichreg$loo)
loo_waic_dirichreg$waic$modelo <- rownames(loo_waic_dirichreg$waic)

# Reordenar columnas para dejar 'modelo' al principio
loo_waic_dirichreg$loo <- loo_waic_dirichreg$loo[, c("modelo", setdiff(names(loo_waic_dirichreg$loo), "modelo"))]
loo_waic_dirichreg$waic <- loo_waic_dirichreg$waic[, c("modelo", setdiff(names(loo_waic_dirichreg$waic), "modelo"))]

# Guardar en disco
saveRDS(
  loo_waic_dirichreg,
  file = here::here("DatosProcesados", "loo_waic_dirichreg.rds")
)

cat("\nLista guardada en: ", here::here("Resultados", "loo_waic_dirichreg.rds"), "\n", sep = "")
