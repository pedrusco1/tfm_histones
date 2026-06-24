# ==============================================================================
# SCRIPT DE COMPARACIÓN DE MODELOS (LOO & WAIC) - COMPLETO
# Incluye modelos Dirichlet, Epistasia y los modelos LogNormal v1
#
# Versión modificada:
#   - Mantiene compatibilidad con los 2 LogNormal v1 previos.
#   - Añade los 6 últimos modelos LogNormal v1 ajustados:
#       * sum0coef
#       * dia_sum0coef
#       * delta_comun
#       * dia_delta_comun
#       * sum0coef_delta_comun
#       * dia_sum0coef_delta_comun
#   - También detecta la lista `fits_lognor_v1_6modelos` si fue guardada por
#     ScriptsR/56_fit_lognor_v1_6_modelos_h3917.R.
# ==============================================================================

# 1. Cargar librerías y entornos compartidos ------------------------------------
library(rstan)
library(loo)

source("ScriptsR/40_stan_config.R")      # Carga rutas como 'data_dir' y 'stan_dir'
source("Funciones/fun_loo_waic_stan.R")  # Carga compute_loo() y compute_waic()

# 2. Cargar outputs de los scripts de ajuste ------------------------------------
message("Cargando archivos de resultados...")

load_if_exists <- function(path) {
  if (file.exists(path)) {
    message("  - Cargando: ", path)
    load(path, envir = .GlobalEnv)
    return(TRUE)
  } else {
    message("  - No encontrado: ", path)
    return(FALSE)
  }
}

load_if_exists(file.path(data_dir, "stan_dirichlet_fits.Rdata"))
load_if_exists(file.path(data_dir, "stan_epistasia_fits.Rdata"))

# Ajustes lognormal v1.
# Se prueban varios nombres habituales para que el script sea robusto frente
# a pequeñas diferencias en cómo se guardaron los fits.
lognormal_files <- c(
  # Salidas antiguas / alternativas
  file.path(data_dir, "stan_lognormal_fits.Rdata"),
  file.path(data_dir, "stan_lognor_fits.Rdata"),
  file.path(data_dir, "stan_lognormal_v1_fits.Rdata"),
  file.path(data_dir, "stan_lognor_v1_fits.Rdata"),
  file.path(data_dir, "stan_lognormal_v1_y_dia_fits.Rdata"),
  file.path(data_dir, "stan_lognor_v1_y_dia_fits.Rdata"),

  # Salida de los dos modelos sum0coef
  file.path(data_dir, "stan_lognormal_v1_sum0coef_fits_h3917.Rdata"),

  # Salida de los 6 últimos modelos lognormal v1
  file.path(data_dir, "stan_lognormal_v1_6modelos_fits.Rdata")
)

invisible(lapply(unique(lognormal_files), load_if_exists))

# 3. Construir la lista de modelos de forma dinámica -----------------------------
lista_modelos <- list()

add_model_if_exists <- function(lista, nombre_salida, nombre_objeto) {
  if (exists(nombre_objeto, envir = .GlobalEnv)) {
    lista[[nombre_salida]] <- get(nombre_objeto, envir = .GlobalEnv)
  }
  lista
}

add_model_from_list <- function(lista, nombre_salida, nombre_lista, nombre_elemento) {
  if (exists(nombre_lista, envir = .GlobalEnv)) {
    x <- get(nombre_lista, envir = .GlobalEnv)
    if (is.list(x) && nombre_elemento %in% names(x)) {
      lista[[nombre_salida]] <- x[[nombre_elemento]]
    }
  }
  lista
}

# --- Bloque 1: Modelos Dirichlet ---
lista_modelos <- add_model_if_exists(lista_modelos, "Dirich_Phi_Comun",     "fit_dirich_phi_comun")
lista_modelos <- add_model_if_exists(lista_modelos, "Dirich_Phi_Cond",      "fit_dirich_phi_cond")
lista_modelos <- add_model_if_exists(lista_modelos, "Dirich_Phi_Comun_Dia", "fit_dirich_phi_comun_dia")
lista_modelos <- add_model_if_exists(lista_modelos, "Dirich_Phi_Cond_Dia",  "fit_dirich_phi_cond_dia")

# --- Bloque 2: Modelos Epistasia ---
lista_modelos <- add_model_if_exists(lista_modelos, "Epist_Kappa_Comun",     "fit_kappa_comun")
lista_modelos <- add_model_if_exists(lista_modelos, "Epist_Kappa_Cond",      "fit_kappa_cond")
lista_modelos <- add_model_if_exists(lista_modelos, "Epist_Kappa_Comun_Dia", "fit_kappa_comun_dia")
lista_modelos <- add_model_if_exists(lista_modelos, "Epist_Kappa_Cond_Dia",  "fit_kappa_cond_dia")

# --- Bloque 3: Modelos LogNormal v1 antiguos / base ---
# Nombres esperados si se usó el script de dos modelos v1:
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1",     "fit_lognor_v1")
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Dia", "fit_lognor_v1_dia")

# Nombres alternativos usados en scripts previos:
if (!"LogNormal_v1" %in% names(lista_modelos)) {
  lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1", "fit_lognor")
}
if (!"LogNormal_v1" %in% names(lista_modelos)) {
  lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1", "fit_lognormal")
}
if (!"LogNormal_v1_Dia" %in% names(lista_modelos)) {
  lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Dia", "fit_lognor_dia")
}
if (!"LogNormal_v1_Dia" %in% names(lista_modelos)) {
  lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Dia", "fit_lognormal_dia")
}

# --- Bloque 4: 6 últimos modelos LogNormal v1 ---
# Primero se intentan tomar desde la lista guardada por 56_fit_lognor_v1_6_modelos_h3917.R.
lista_modelos <- add_model_from_list(
  lista_modelos,
  "LogNormal_v1_Sum0Coef",
  "fits_lognor_v1_6modelos",
  "fit_lognor_v1_sum0coef"
)

lista_modelos <- add_model_from_list(
  lista_modelos,
  "LogNormal_v1_Dia_Sum0Coef",
  "fits_lognor_v1_6modelos",
  "fit_lognor_v1_dia_sum0coef"
)

lista_modelos <- add_model_from_list(
  lista_modelos,
  "LogNormal_v1_Delta_Comun",
  "fits_lognor_v1_6modelos",
  "fit_lognor_v1_delta_comun"
)

lista_modelos <- add_model_from_list(
  lista_modelos,
  "LogNormal_v1_Dia_Delta_Comun",
  "fits_lognor_v1_6modelos",
  "fit_lognor_v1_dia_delta_comun"
)

lista_modelos <- add_model_from_list(
  lista_modelos,
  "LogNormal_v1_Sum0Coef_Delta_Comun",
  "fits_lognor_v1_6modelos",
  "fit_lognor_v1_sum0coef_delta_comun"
)

lista_modelos <- add_model_from_list(
  lista_modelos,
  "LogNormal_v1_Dia_Sum0Coef_Delta_Comun",
  "fits_lognor_v1_6modelos",
  "fit_lognor_v1_dia_sum0coef_delta_comun"
)

# Si no estaba la lista, se intentan los objetos individuales guardados en el .Rdata.
# Esto también cubre el caso en que el script de ajuste hizo list2env().
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Sum0Coef",                  "fit_lognor_v1_sum0coef")
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Dia_Sum0Coef",              "fit_lognor_v1_dia_sum0coef")
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Delta_Comun",               "fit_lognor_v1_delta_comun")
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Dia_Delta_Comun",           "fit_lognor_v1_dia_delta_comun")
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Sum0Coef_Delta_Comun",      "fit_lognor_v1_sum0coef_delta_comun")
lista_modelos <- add_model_if_exists(lista_modelos, "LogNormal_v1_Dia_Sum0Coef_Delta_Comun",  "fit_lognor_v1_dia_sum0coef_delta_comun")

# Eliminar entradas duplicadas por nombre, manteniendo la última versión cargada.
# Esto es útil si un modelo aparece tanto en la lista como como objeto individual.
lista_modelos <- lista_modelos[!duplicated(names(lista_modelos), fromLast = TRUE)]

message("Se han detectado y cargado ", length(lista_modelos), " modelos para evaluar.")
message("Modelos detectados: ", paste(names(lista_modelos), collapse = ", "))

if (length(lista_modelos) < 2) {
  stop("Se necesitan al menos dos modelos cargados para comparar LOO/WAIC.")
}

lognormal_6_esperados <- c(
  "LogNormal_v1_Sum0Coef",
  "LogNormal_v1_Dia_Sum0Coef",
  "LogNormal_v1_Delta_Comun",
  "LogNormal_v1_Dia_Delta_Comun",
  "LogNormal_v1_Sum0Coef_Delta_Comun",
  "LogNormal_v1_Dia_Sum0Coef_Delta_Comun"
)

lognormal_6_faltantes <- setdiff(lognormal_6_esperados, names(lista_modelos))
if (length(lognormal_6_faltantes) > 0) {
  warning(
    "No se han detectado algunos de los 6 últimos modelos lognormal: ",
    paste(lognormal_6_faltantes, collapse = ", "),
    "\nComprueba que existe el archivo: ",
    file.path(data_dir, "stan_lognormal_v1_6modelos_fits.Rdata")
  )
}

# 4. Calcular LOO y WAIC ---------------------------------------------------------
message("Calculando LOO de manera secuencial...")
resultados_loo  <- lapply(lista_modelos, compute_loo)

message("Calculando WAIC de manera secuencial...")
resultados_waic <- lapply(lista_modelos, compute_waic)

# 5. Funciones auxiliares de comparación ----------------------------------------
compare_block <- function(resultados, patron, titulo) {
  idx <- grep(patron, names(resultados))
  if (length(idx) > 1) {
    cat("\n--- ", titulo, " ---\n", sep = "")
    comp <- loo::loo_compare(resultados[idx])
    print(comp, simplify = FALSE)
    out <- as.data.frame(comp)
    out$modelo <- rownames(out)
    out <- out[, c("modelo", setdiff(names(out), "modelo"))]
    return(out)
  } else if (length(idx) == 1) {
    cat("\n--- ", titulo, " ---\n", sep = "")
    cat("Solo hay un modelo en este bloque: ", names(resultados)[idx], "\n", sep = "")
  }
  return(NULL)
}

safe_global_compare <- function(resultados, etiqueta) {
  tryCatch({
    comp <- loo::loo_compare(resultados)
    print(comp, simplify = FALSE)
    out <- as.data.frame(comp)
    out$modelo <- rownames(out)
    out <- out[, c("modelo", setdiff(names(out), "modelo"))]
    out
  }, error = function(e) {
    cat("Aviso: No se pueden comparar todos los modelos juntos en una única matriz.\n")
    cat("Motivo: ", conditionMessage(e), "\n", sep = "")
    cat("Separando comparaciones por bloques independientes.\n")
    NULL
  })
}

# 6. Comparación global y por bloques -------------------------------------------
# loo_compare exige que los modelos compartan exactamente la misma variable
# respuesta observada y el mismo número de observaciones. Primero se intenta una
# comparación global; si falla, se comparan bloques independientes.

cat("\n==================================================\n")
cat("      COMPARACIÓN DE MODELOS BASADA EN LOO       \n")
cat("==================================================\n")

comp_loo_global <- safe_global_compare(resultados_loo, "LOO")

comp_loo_bloques <- list(
  Dirichlet = compare_block(resultados_loo,  "^Dirich_",    "Bloque Regresión Dirichlet (LOO)"),
  Epistasia = compare_block(resultados_loo,  "^Epist_",     "Bloque Campo Epistasia (LOO)"),
  LogNormal = compare_block(resultados_loo,  "^LogNormal_", "Bloque LogNormal v1 completo (LOO)")
)

cat("\n==================================================\n")
cat("      COMPARACIÓN DE MODELOS BASADA EN WAIC      \n")
cat("==================================================\n")

comp_waic_global <- safe_global_compare(resultados_waic, "WAIC")

comp_waic_bloques <- list(
  Dirichlet = compare_block(resultados_waic, "^Dirich_",    "Bloque Regresión Dirichlet (WAIC)"),
  Epistasia = compare_block(resultados_waic, "^Epist_",     "Bloque Campo Epistasia (WAIC)"),
  LogNormal = compare_block(resultados_waic, "^LogNormal_", "Bloque LogNormal v1 completo (WAIC)")
)

# 7. Guardar resultados ----------------------------------------------------------
loo_waic_todos_modelos <- list(
  modelos_evaluados = names(lista_modelos),
  loo = resultados_loo,
  waic = resultados_waic,
  comparacion_loo_global = comp_loo_global,
  comparacion_waic_global = comp_waic_global,
  comparacion_loo_bloques = comp_loo_bloques,
  comparacion_waic_bloques = comp_waic_bloques
)

out_file <- file.path(data_dir, "loo_waic_modelos_con_lognormal.rds")

saveRDS(
  loo_waic_todos_modelos,
  file = out_file
)

message("Resultados guardados en: ", out_file)
