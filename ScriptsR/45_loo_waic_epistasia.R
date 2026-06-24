##################################################
# MODELOS EPISTASIA: LOO-WAIC Y COMPARACIÓN     #
##################################################

# 1. Cargar configuración y funciones necesarias
source("ScriptsR/40_stan_config.R")
# Cargamos el script que contiene compute_loo y compute_waic (ajusta la ruta si es necesario)
source("Funciones/fun_loo_waic_stan.R") 

# Cargar los modelos ajustados de epistasia
load(file.path(data_dir, "stan_epistasia_fits.Rdata"))

# 2. Obtener la lista de nombres de los modelos de forma dinámica
# Buscamos los objetos que empiecen por 'fit_kappa_'
model_names <- ls(pattern = "^fit_kappa_")

if (length(model_names) == 0) {
  stop("No se encontraron objetos con nombres que empiecen por 'fit_kappa_'.")
}

# Recuperamos los objetos stanfit a partir de su nombre en texto
list_of_fits <- mget(model_names)

# 3. Calcular LOO y WAIC para todos los modelos
cat("Calculando LOO para los modelos de epistasia...\n")
loo_results <- lapply(list_of_fits, compute_loo)
names(loo_results) <- model_names

cat("Calculando WAIC para los modelos de epistasia...\n")
waic_results <- lapply(list_of_fits, compute_waic)
names(waic_results) <- model_names

# 4. Comparación multimodelo usando el paquete 'loo'
cat("\n--- COMPARACIÓN LOO (EPISTASIA) ---\n")
comp_loo <- loo::loo_compare(loo_results)
print(comp_loo, simplify = FALSE)

cat("\n--- COMPARACIÓN WAIC (EPISTASIA) ---\n")
comp_waic <- loo::loo_compare(waic_results)
print(comp_waic, simplify = FALSE)

# 5. Estructurar los resultados en data.frames limpios
loo_waic_epistasia <- list(
  loo = as.data.frame(comp_loo),
  waic = as.data.frame(comp_waic)
)

# Añadir el nombre del modelo como columna principal
loo_waic_epistasia$loo$modelo <- rownames(loo_waic_epistasia$loo)
loo_waic_epistasia$waic$modelo <- rownames(loo_waic_epistasia$waic)

# Reordenar columnas para dejar 'modelo' al principio
loo_waic_epistasia$loo <- loo_waic_epistasia$loo[, c("modelo", setdiff(names(loo_waic_epistasia$loo), "modelo"))]
loo_waic_epistasia$waic <- loo_waic_epistasia$waic[, c("modelo", setdiff(names(loo_waic_epistasia$waic), "modelo"))]

# 6. Guardar en el directorio de datos procesados
output_file <- file.path(data_dir, "loo_waic_epistasia.rds")
saveRDS(loo_waic_epistasia, file = output_file)

cat("\nResultados guardados con éxito en: ", output_file, "\n", sep = "")