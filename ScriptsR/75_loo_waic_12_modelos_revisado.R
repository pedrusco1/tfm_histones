############################################################
# COMPARACIÓN GLOBAL: 8 MODELOS EXISTENTES + 4 LOG-NORMALES REVISADOS       #
#                                                          #
# Incluye:                                                 #
#   - 4 modelos Dirichlet regression: fit_dirich_*         #
#   - 4 modelos de epistasia:          fit_kappa_*         #
#   - 4 modelos log-normal revisados:   fit_lognor_v1_*_revisado     #
############################################################
rm(list = ls())
# ----------------------------------------------------------
# 0. Configuración
# ----------------------------------------------------------
source(here::here("ScriptsR","40_stan_config.R"))
source(here::here("Funciones", "fun_loo_waic_stan.R"))

if (!requireNamespace("loo", quietly = TRUE)) {
  stop("Hace falta instalar el paquete 'loo'.")
}

# Directorio de salida. Usamos data_dir si existe en 40_stan_config.R;
# si no, usamos DatosProcesados.
if (!exists("data_dir")) {
  data_dir <- here::here("DatosProcesados")
}

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# ----------------------------------------------------------
# 1. Especificación de familias de modelos
# ----------------------------------------------------------
# Cada familia se carga en su propio entorno para evitar mezclar objetos
# de distintos ficheros y para que los patrones detecten solo lo que toca.

familias_modelos <- list(
  dirichlet = list(
    archivo = here::here("DatosProcesados", "stan_dirichlet_fits.Rdata"),
    patron  = "^fit_dirich_"
  ),
  epistasia = list(
    archivo = file.path(data_dir, "stan_epistasia_fits.Rdata"),
    patron  = "^fit_kappa_"
  ),
  lognormal_sum0 = list(
    archivo = file.path(data_dir, "stan_lognormal_sum0_fits_revisados.Rdata"),
    patron  = "^fit_lognor_v1_.*_revisado$"
  )
)

# ----------------------------------------------------------
# 2. Funciones auxiliares
# ----------------------------------------------------------

comparar_familia <- function(nombre_familia, spec) {
  cat("\n============================================================\n")
  cat("Familia:", nombre_familia, "\n")
  cat("Archivo:", spec$archivo, "\n")
  cat("============================================================\n")

  if (!file.exists(spec$archivo)) {
    stop("No se encuentra el archivo de ajustes para la familia '",
         nombre_familia, "': ", spec$archivo)
  }

  env_familia <- new.env(parent = globalenv())
  load(spec$archivo, envir = env_familia)

  model_names <- ls(envir = env_familia, pattern = spec$patron)

  if (length(model_names) == 0) {
    stop("No se encontraron modelos para la familia '", nombre_familia,
         "' con el patrón ", spec$patron)
  }

  model_names <- sort(model_names)
  cat("Modelos detectados (", length(model_names), "):\n", sep = "")
  cat(paste("  -", model_names), sep = "\n")
  cat("\n")

  fits <- mget(model_names, envir = env_familia)

  cat("Calculando LOO...\n")
  loo_results <- lapply(fits, compute_loo)
  names(loo_results) <- paste(nombre_familia, model_names, sep = "__")

  cat("Calculando WAIC...\n")
  waic_results <- lapply(fits, compute_waic)
  names(waic_results) <- paste(nombre_familia, model_names, sep = "__")

  list(
    familia = nombre_familia,
    modelos = model_names,
    loo = loo_results,
    waic = waic_results
  )
}

comparacion_a_tabla <- function(comp, criterio) {
  tab <- as.data.frame(comp)
  tab$modelo_global <- rownames(tab)
  rownames(tab) <- NULL

  partes <- strsplit(tab$modelo_global, "__", fixed = TRUE)
  tab$familia <- vapply(partes, `[`, character(1), 1)
  tab$modelo  <- vapply(partes, `[`, character(1), 2)
  tab$criterio <- criterio

  tab <- tab[, c(
    "criterio", "familia", "modelo", "modelo_global",
    setdiff(names(tab), c("criterio", "familia", "modelo", "modelo_global"))
  )]

  tab
}

# ----------------------------------------------------------
# 3. Calcular LOO y WAIC por familia
# ----------------------------------------------------------
resultados_por_familia <- Map(
  f = comparar_familia,
  nombre_familia = names(familias_modelos),
  spec = familias_modelos
)

names(resultados_por_familia) <- names(familias_modelos)

# Unificar los objetos loo y waic de las 3 familias en listas globales
loo_todos <- unlist(
  lapply(resultados_por_familia, `[[`, "loo"),
  recursive = FALSE
)

waic_todos <- unlist(
  lapply(resultados_por_familia, `[[`, "waic"),
  recursive = FALSE
)

cat("\n============================================================\n")
cat("Número total de modelos detectados:", length(loo_todos), "\n")
cat("============================================================\n")

if (length(loo_todos) != 12) {
  warning("Se esperaban 12 modelos, pero se han detectado ", length(loo_todos),
          ". Revisa los patrones o los objetos guardados en los .Rdata.")
}

# ----------------------------------------------------------
# 4. Comparación global entre los 12 modelos
# ----------------------------------------------------------
cat("\n============================================================\n")
cat("COMPARACIÓN GLOBAL LOO-IC\n")
cat("============================================================\n")
comp_loo_global <- loo::loo_compare(loo_todos)
print(comp_loo_global, simplify = FALSE)

cat("\n============================================================\n")
cat("COMPARACIÓN GLOBAL WAIC\n")
cat("============================================================\n")
comp_waic_global <- loo::loo_compare(waic_todos)
print(comp_waic_global, simplify = FALSE)

# Tablas limpias para guardar/exportar
tabla_loo_global  <- comparacion_a_tabla(comp_loo_global,  "loo")
tabla_waic_global <- comparacion_a_tabla(comp_waic_global, "waic")
unir_tablas_modelos <- function(...) {
  tablas <- list(...)
  todas_columnas <- unique(unlist(lapply(tablas, names)))
  
  tablas_alineadas <- lapply(tablas, function(tab) {
    faltan <- setdiff(todas_columnas, names(tab))
    for (col in faltan) {
      tab[[col]] <- NA
    }
    tab[, todas_columnas, drop = FALSE]
  })
  
  do.call(rbind, tablas_alineadas)
}

tabla_comparacion_global <- unir_tablas_modelos(
  tabla_loo_global,
  tabla_waic_global
)
# ----------------------------------------------------------
# 5. Comparaciones dentro de cada familia
# ----------------------------------------------------------
comparaciones_por_familia <- lapply(resultados_por_familia, function(res) {
  comp_loo  <- loo::loo_compare(res$loo)
  comp_waic <- loo::loo_compare(res$waic)

  list(
    loo = comp_loo,
    waic = comp_waic,
    tabla_loo = comparacion_a_tabla(comp_loo, "loo"),
    tabla_waic = comparacion_a_tabla(comp_waic, "waic")
  )
})

# ----------------------------------------------------------
# 6. Guardar resultados
# ----------------------------------------------------------
resultados_loo_waic_12_modelos_revisados <- list(
  modelos = names(loo_todos),
  resultados_por_familia = resultados_por_familia,
  loo_todos = loo_todos,
  waic_todos = waic_todos,
  comparacion_loo_global = comp_loo_global,
  comparacion_waic_global = comp_waic_global,
  tabla_loo_global = tabla_loo_global,
  tabla_waic_global = tabla_waic_global,
  tabla_comparacion_global = tabla_comparacion_global,
  comparaciones_por_familia = comparaciones_por_familia
)

rds_file <- file.path(data_dir, "loo_waic_12_modelos_revisados_resultados.rds")
csv_loo_file <- file.path(data_dir, "loo_12_modelos_revisados_comparacion.csv")
csv_waic_file <- file.path(data_dir, "waic_12_modelos_revisados_comparacion.csv")
csv_global_file <- file.path(data_dir, "loo_waic_12_modelos_revisados_comparacion.csv")

saveRDS(resultados_loo_waic_12_modelos_revisados, file = rds_file)
write.csv(tabla_loo_global, file = csv_loo_file, row.names = FALSE)
write.csv(tabla_waic_global, file = csv_waic_file, row.names = FALSE)
write.csv(tabla_comparacion_global, file = csv_global_file, row.names = FALSE)

cat("\n============================================================\n")
cat("Resultados guardados\n")
cat("============================================================\n")
cat("RDS completo:  ", rds_file, "\n", sep = "")
cat("CSV LOO:       ", csv_loo_file, "\n", sep = "")
cat("CSV WAIC:      ", csv_waic_file, "\n", sep = "")
cat("CSV combinado: ", csv_global_file, "\n", sep = "")

# ----------------------------------------------------------
# 7. Resumen rápido en consola
# ----------------------------------------------------------
cat("\nMejor modelo por LOO-IC:\n")
print(tabla_loo_global[1, c("familia", "modelo", "elpd_diff", "se_diff")])

cat("\nMejor modelo por WAIC:\n")
print(tabla_waic_global[1, c("familia", "modelo", "elpd_diff", "se_diff")])

