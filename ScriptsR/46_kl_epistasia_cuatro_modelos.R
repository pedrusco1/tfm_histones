###############################################################
# Aplicar divergencia KL Dirichlet a los cuatro ajustes Stan
# H3917: Modelos Epistasia (Q y Q_cond) con/sin kappa por condición y con/sin día
###############################################################

#-------------------------------------------------------------
# 0. Configuración
#-------------------------------------------------------------

if (file.exists(here::here("ScriptsR", "40_stan_config.R"))) {
  source(here::here("ScriptsR", "40_stan_config.R"))
}

if (!exists("data_dir")) {
  data_dir <- here::here("DatosProcesados")
}

output_dir <- file.path(data_dir, "resultados_kl_epistasia")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "stan_epistasia_fits.Rdata")
if (!file.exists(input_file)) {
  input_file <- file.path(data_dir, "stan_episia_fits.Rdata")
}

# Las condiciones en estos modelos suelen estar codificadas como 1=async, 2=mitosis.
# Si son NULL, la función detectará los dos niveles del vector automáticamente.
control_level <- NULL
tratamiento_level <- NULL

# TRUE = 0.5 * KL(P||Q) + 0.5 * KL(Q||P)
# FALSE = KL(Control||Tratamiento)
symmetric <- TRUE

# Guardar también los draws completos por modelo
save_draws <- TRUE

#-------------------------------------------------------------
# 1. Funciones auxiliares (Adaptadas a Q y kappa)
#-------------------------------------------------------------

kl_dirichlet_unidirectional <- function(a1, a2) {
  if (any(!is.finite(a1)) || any(!is.finite(a2))) return(NA_real_)
  if (any(a1 <= 0) || any(a2 <= 0)) return(NA_real_)
  
  sum_a1 <- sum(a1)
  sum_a2 <- sum(a2)
  
  lgamma(sum_a1) - lgamma(sum_a2) -
    sum(lgamma(a1)) + sum(lgamma(a2)) +
    sum((a1 - a2) * (digamma(a1) - digamma(sum_a1)))
}

kl_dirichlet_draws <- function(Q_control,
                               kappa_control,
                               Q_tratamiento,
                               kappa_tratamiento,
                               symmetric = TRUE) {
  ndraws <- nrow(Q_control)
  out <- numeric(ndraws)
  
  for (i in seq_len(ndraws)) {
    # alpha = kappa * Q
    a1 <- as.numeric(Q_control[i, ]) * as.numeric(kappa_control[i])
    a2 <- as.numeric(Q_tratamiento[i, ]) * as.numeric(kappa_tratamiento[i])
    
    kl12 <- kl_dirichlet_unidirectional(a1, a2)
    
    if (isTRUE(symmetric)) {
      kl21 <- kl_dirichlet_unidirectional(a2, a1)
      out[i] <- 0.5 * (kl12 + kl21)
    } else {
      out[i] <- kl12
    }
  }
  out
}

summarise_draws <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q025 = unname(quantile(x, 0.025, na.rm = TRUE)),
    q975 = unname(quantile(x, 0.975, na.rm = TRUE)),
    n_draws = sum(!is.na(x))
  )
}

as_draws_df_safe <- function(fit) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Necesitas instalar el paquete 'posterior'")
  }
  out <- try(posterior::as_draws_df(fit), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(as.data.frame(out, check.names = FALSE))
  }
  if (inherits(fit, "CmdStanMCMC")) {
    out <- posterior::as_draws_df(fit$draws())
    return(as.data.frame(out, check.names = FALSE))
  }
  stop("No he podido convertir el fit a draws.")
}

# NUEVA FUNCIÓN: Extrae Q o Q_cond directamente de los transformed parameters
extract_Q_for_condition <- function(draws_df, cond_level, K) {
  ndraws <- nrow(draws_df)
  out <- matrix(NA_real_, nrow = ndraws, ncol = K)
  
  # Buscar si el modelo usa Q (sin día) o Q_cond (con día)
  base_name <- NULL
  if (paste0("Q[1,", cond_level, "]") %in% names(draws_df)) {
    base_name <- "Q"
  } else if (paste0("Q_cond[1,", cond_level, "]") %in% names(draws_df)) {
    base_name <- "Q_cond"
  } else {
    stop("No encuentro la matriz Q ni Q_cond para la condicion ", cond_level, ". Revisa si la condición existe en el modelo.")
  }
  
  for (k in seq_len(K)) {
    # CmdStanR/RStan vectorizan matrices en orden de columnas: parametro[fila, columna]
    col_name <- paste0(base_name, "[", k, ",", cond_level, "]")
    if (!col_name %in% names(draws_df)) {
      stop("Falta la columna esperada: ", col_name)
    }
    out[, k] <- as.numeric(draws_df[[col_name]])
  }
  colnames(out) <- paste0("comp_", seq_len(K))
  return(out)
}

extract_kappa_for_condition <- function(draws_df, cond_level) {
  # Modelo con kappa comun
  if ("kappa" %in% names(draws_df)) {
    return(as.numeric(draws_df[["kappa"]]))
  }
  
  # Modelo con kappa por condicion
  kappa_name <- paste0("kappa[", cond_level, "]")
  if (kappa_name %in% names(draws_df)) {
    return(as.numeric(draws_df[[kappa_name]]))
  }
  
  stop("No encuentro el parámetro kappa correspondiente a la condicion ", cond_level, ".")
}

calc_kl_for_epistasia_model <- function(fit,
                                        condition_vector,
                                        K,
                                        control_level,
                                        tratamiento_level,
                                        symmetric = TRUE) {
  
  idx_control <- which(condition_vector == control_level)
  idx_trat <- which(condition_vector == tratamiento_level)
  
  if (length(idx_control) == 0) stop("No hay observaciones para control_level = ", control_level)
  if (length(idx_trat) == 0) stop("No hay observaciones para tratamiento_level = ", tratamiento_level)
  
  draws_df <- as_draws_df_safe(fit)
  draws_df <- draws_df[, !grepl("^\\.", names(draws_df)), drop = FALSE]
  
  # Extraemos Q (o Q_cond) y kappa
  Q_control <- extract_Q_for_condition(draws_df, control_level, K)
  Q_trat <- extract_Q_for_condition(draws_df, tratamiento_level, K)
  
  kappa_control <- extract_kappa_for_condition(draws_df, control_level)
  kappa_trat <- extract_kappa_for_condition(draws_df, tratamiento_level)
  
  kl <- kl_dirichlet_draws(
    Q_control = Q_control,
    kappa_control = kappa_control,
    Q_tratamiento = Q_trat,
    kappa_tratamiento = kappa_trat,
    symmetric = symmetric
  )
  
  list(
    kl_draws = kl,
    Q_control = Q_control,
    Q_tratamiento = Q_trat,
    kappa_control = kappa_control,
    kappa_tratamiento = kappa_trat,
    n_control = length(idx_control),
    n_tratamiento = length(idx_trat)
  )
}

choose_two_levels <- function(cond, control_level = NULL, tratamiento_level = NULL) {
  lev <- sort(unique(cond))
  if (length(lev) != 2 && (is.null(control_level) || is.null(tratamiento_level))) {
    stop("Hay ", length(lev), " niveles de condicion. Define manualmente.")
  }
  if (is.null(control_level)) control_level <- lev[1]
  if (is.null(tratamiento_level)) tratamiento_level <- lev[2]
  list(control = control_level, tratamiento = tratamiento_level)
}

#-------------------------------------------------------------
# 2. Cargar datos base y resultados Stan
#-------------------------------------------------------------

# Cargamos los datos originales para obtener K y el vector cond
load(file.path(data_dir, "data_inputs_h3917.Rdata"))
K <- ncol(Y_mat) # Hay 20 proteoformas
condition_vector <- cond

if (!file.exists(input_file)) {
  stop("No encuentro el fichero de modelos de epistasia: ", input_file)
}
load(input_file)

levels_used <- choose_two_levels(condition_vector, control_level, tratamiento_level)
control_level <- levels_used$control
tratamiento_level <- levels_used$tratamiento

message("Comparando condiciones: ", control_level, " vs ", tratamiento_level)
message("Divergencia: ", ifelse(symmetric, "KL simétrica", "KL(Control||Tratamiento)"))

#-------------------------------------------------------------
# 3. Preparar lista de modelos a procesar
#-------------------------------------------------------------

fits_epistasia <- list()
if (exists("fit_kappa_comun")) fits_epistasia$kappa_comun <- fit_kappa_comun
if (exists("fit_kappa_comun_dia")) fits_epistasia$kappa_comun_dia <- fit_kappa_comun_dia
if (exists("fit_kappa_cond")) fits_epistasia$kappa_cond <- fit_kappa_cond
if (exists("fit_kappa_cond_dia")) fits_epistasia$kappa_cond_dia <- fit_kappa_cond_dia

if (length(fits_epistasia) == 0) {
  stop("No se han encontrado objetos fit_kappa_* en el archivo cargado.")
}

#-------------------------------------------------------------
# 4. Calcular KL en cada modelo
#-------------------------------------------------------------

kl_results <- list()
summary_rows <- list()

a <- 1
for (model_name in names(fits_epistasia)) {
  message("Procesando modelo: ", model_name)
  
  res <- calc_kl_for_epistasia_model(
    fit = fits_epistasia[[model_name]],
    condition_vector = condition_vector,
    K = K,
    control_level = control_level,
    tratamiento_level = tratamiento_level,
    symmetric = symmetric
  )
  
  kl_results[[model_name]] <- res
  
  s <- summarise_draws(res$kl_draws)
  summary_rows[[a]] <- data.frame(
    modelo = model_name,
    control_level = as.character(control_level),
    tratamiento_level = as.character(tratamiento_level),
    symmetric = symmetric,
    n_control = res$n_control,
    n_tratamiento = res$n_tratamiento,
    mean = unname(s["mean"]),
    sd = unname(s["sd"]),
    median = unname(s["median"]),
    q025 = unname(s["q025"]),
    q975 = unname(s["q975"]),
    n_draws = unname(s["n_draws"]),
    row.names = NULL
  )
  a <- a + 1
}

kl_summary <- do.call(rbind, summary_rows)

#-------------------------------------------------------------
# 5. Guardar resultados
#-------------------------------------------------------------

suffix <- ifelse(symmetric, "kl_simetrica", "kl_control_vs_tratamiento")
summary_file <- file.path(output_dir, paste0("resumen_", suffix, "_epistasia.csv"))
rdata_file <- file.path(output_dir, paste0("draws_", suffix, "_epistasia.Rdata"))

write.csv(kl_summary, summary_file, row.names = FALSE)

if (isTRUE(save_draws)) {
  save(
    kl_results,
    kl_summary,
    control_level,
    tratamiento_level,
    symmetric,
    file = rdata_file
  )
}

print(kl_summary)
message("Resumen guardado en: ", summary_file)
if (isTRUE(save_draws)) message("Draws guardados en: ", rdata_file)
