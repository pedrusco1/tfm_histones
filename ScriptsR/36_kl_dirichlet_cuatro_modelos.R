###############################################################
# Aplicar divergencia KL Dirichlet a los cuatro ajustes Stan
# H3917: modelos Dirichlet con/sin phi por condicion y con/sin dia
###############################################################

# Este script carga los resultados guardados por 33_fit_dirichreg_models.R:
#   file.path(data_dir, "stan_dirichlet_fits.Rdata")
# y calcula, para cada modelo ajustado, la divergencia entre las dos condiciones.
#
# La comparacion se hace en la escala generativa Dirichlet:
#   alpha = mu * phi
# usando los draws posteriores de mu y phi.
#
# Para cada draw:
#   1) mu_control      = media de mu[n, ] en las observaciones control
#   2) mu_tratamiento  = media de mu[n, ] en las observaciones tratamiento
#   3) phi_control / phi_tratamiento segun el modelo
#   4) KL(Dir(alpha_control) || Dir(alpha_tratamiento))
#      o KL simetrica si symmetric = TRUE

#-------------------------------------------------------------
# 0. Configuracion
#-------------------------------------------------------------

# Si tu script de configuracion define data_dir, se usa automaticamente.
if (file.exists(here::here("ScriptsR", "40_stan_config.R"))) {
  source(here::here("ScriptsR", "40_stan_config.R"))
}

if (!exists("data_dir")) {
  data_dir <- "DatosProcesados"
}

input_file <- file.path(data_dir, "stan_dirichlet_fits.Rdata")
output_dir <- file.path(data_dir, "resultados_kl_dirichlet")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Comparacion. Si son NULL se usan automaticamente las dos condiciones encontradas.
control_level <- NULL
tratamiento_level <- NULL

# TRUE = 0.5 * KL(P||Q) + 0.5 * KL(Q||P)
# FALSE = KL(Control||Tratamiento)
symmetric <- TRUE

# Guardar tambien los draws completos por modelo
save_draws <- TRUE

#-------------------------------------------------------------
# 1. Funciones auxiliares
#-------------------------------------------------------------

kl_dirichlet_unidirectional <- function(a1, a2) {
  # KL(Dir(a1) || Dir(a2))
  if (any(!is.finite(a1)) || any(!is.finite(a2))) return(NA_real_)
  if (any(a1 <= 0) || any(a2 <= 0)) return(NA_real_)

  sum_a1 <- sum(a1)
  sum_a2 <- sum(a2)

  lgamma(sum_a1) - lgamma(sum_a2) -
    sum(lgamma(a1)) + sum(lgamma(a2)) +
    sum((a1 - a2) * (digamma(a1) - digamma(sum_a1)))
}

kl_dirichlet_draws <- function(mu_control,
                               phi_control,
                               mu_tratamiento,
                               phi_tratamiento,
                               symmetric = TRUE) {
  ndraws <- nrow(mu_control)
  out <- numeric(ndraws)

  for (i in seq_len(ndraws)) {
    a1 <- as.numeric(mu_control[i, ]) * as.numeric(phi_control[i])
    a2 <- as.numeric(mu_tratamiento[i, ]) * as.numeric(phi_tratamiento[i])

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
    stop("Necesitas instalar el paquete 'posterior': install.packages('posterior')")
  }

  # CmdStanR y RStan suelen ser compatibles con posterior::as_draws_df().
  out <- try(posterior::as_draws_df(fit), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(as.data.frame(out, check.names = FALSE))
  }

  # Fallback para objetos CmdStanMCMC.
  if (inherits(fit, "CmdStanMCMC")) {
    out <- posterior::as_draws_df(fit$draws())
    return(as.data.frame(out, check.names = FALSE))
  }

  stop("No he podido convertir el fit a draws. Revisa si es objeto cmdstanr/rstan compatible.")
}

get_stan_dims <- function(stan_data) {
  if (!is.null(stan_data$N) && !is.null(stan_data$K)) {
    return(list(N = stan_data$N, K = stan_data$K))
  }
  if (!is.null(stan_data$Y)) {
    return(list(N = nrow(stan_data$Y), K = ncol(stan_data$Y)))
  }
  stop("No puedo inferir N y K desde stan_data.")
}

get_condition_vector <- function(stan_data_ref, n_expected) {
  cond <- NULL

  if (!is.null(stan_data_ref$condition)) cond <- stan_data_ref$condition
  if (is.null(cond) && !is.null(stan_data_ref$condicion)) cond <- stan_data_ref$condicion

  if (is.null(cond)) {
    stop("No encuentro vector de condicion en los objetos stan_data guardados.")
  }
  if (length(cond) != n_expected) {
    stop("El vector de condicion no tiene longitud N.")
  }

  cond
}

choose_two_levels <- function(cond, control_level = NULL, tratamiento_level = NULL) {
  lev <- sort(unique(cond))
  if (length(lev) != 2 && (is.null(control_level) || is.null(tratamiento_level))) {
    stop("Hay ", length(lev), " niveles de condicion. Define control_level y tratamiento_level manualmente.")
  }

  if (is.null(control_level)) control_level <- lev[1]
  if (is.null(tratamiento_level)) tratamiento_level <- lev[2]

  list(control = control_level, tratamiento = tratamiento_level)
}

mean_mu_by_condition <- function(draws_df, idx, K) {
  ndraws <- nrow(draws_df)
  out <- matrix(NA_real_, nrow = ndraws, ncol = K)

  for (k in seq_len(K)) {
    cols <- paste0("mu[", idx, ",", k, "]")
    missing_cols <- setdiff(cols, names(draws_df))
    if (length(missing_cols) > 0) {
      stop("Faltan columnas de mu, por ejemplo: ", missing_cols[1])
    }
    out[, k] <- rowMeans(as.matrix(draws_df[, cols, drop = FALSE]))
  }

  colnames(out) <- paste0("comp_", seq_len(K))
  out
}

extract_phi_for_condition <- function(draws_df, cond_level, common_ok = TRUE) {
  ndraws <- nrow(draws_df)

  # Modelo con phi comun: parametro 'phi'.
  if ("phi" %in% names(draws_df)) {
    return(as.numeric(draws_df[["phi"]]))
  }

  # Modelo con phi por condicion: parametros 'phi[1]', 'phi[2]', ...
  phi_name <- paste0("phi[", cond_level, "]")
  if (phi_name %in% names(draws_df)) {
    return(as.numeric(draws_df[[phi_name]]))
  }

  # Fallback por si cond_level es factor/character y las condiciones internas son 1:J.
  phi_cols <- grep("^phi\\[[0-9]+\\]$", names(draws_df), value = TRUE)
  if (length(phi_cols) > 0) {
    cond_numeric <- suppressWarnings(as.integer(as.character(cond_level)))
    if (!is.na(cond_numeric)) {
      phi_name <- paste0("phi[", cond_numeric, "]")
      if (phi_name %in% names(draws_df)) {
        return(as.numeric(draws_df[[phi_name]]))
      }
    }
  }

  stop("No encuentro el phi correspondiente a la condicion ", cond_level, ".")
}

calc_kl_for_model <- function(fit,
                              stan_data,
                              condition_vector,
                              control_level,
                              tratamiento_level,
                              symmetric = TRUE) {
  dims <- get_stan_dims(stan_data)
  N <- dims$N
  K <- dims$K

  if (length(condition_vector) != N) {
    stop("condition_vector no coincide con N para este modelo.")
  }

  idx_control <- which(condition_vector == control_level)
  idx_trat <- which(condition_vector == tratamiento_level)

  if (length(idx_control) == 0) stop("No hay observaciones para control_level = ", control_level)
  if (length(idx_trat) == 0) stop("No hay observaciones para tratamiento_level = ", tratamiento_level)

  draws_df <- as_draws_df_safe(fit)
  draws_df <- draws_df[, !grepl("^\\.", names(draws_df)), drop = FALSE]

  mu_control <- mean_mu_by_condition(draws_df, idx_control, K)
  mu_trat <- mean_mu_by_condition(draws_df, idx_trat, K)

  phi_control <- extract_phi_for_condition(draws_df, control_level)
  phi_trat <- extract_phi_for_condition(draws_df, tratamiento_level)

  kl <- kl_dirichlet_draws(
    mu_control = mu_control,
    phi_control = phi_control,
    mu_tratamiento = mu_trat,
    phi_tratamiento = phi_trat,
    symmetric = symmetric
  )

  list(
    kl_draws = kl,
    mu_control = mu_control,
    mu_tratamiento = mu_trat,
    phi_control = phi_control,
    phi_tratamiento = phi_trat,
    n_control = length(idx_control),
    n_tratamiento = length(idx_trat)
  )
}

#-------------------------------------------------------------
# 2. Cargar resultados
#-------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("No encuentro el fichero: ", input_file)
}

load(input_file)

if (!exists("fits_dirich")) {
  stop("El fichero cargado no contiene el objeto fits_dirich.")
}
if (!exists("modelos_dirich")) {
  stop("El fichero cargado no contiene el objeto modelos_dirich.")
}
if (!exists("stan_data_phi_cond")) {
  stop("Necesito stan_data_phi_cond para recuperar el vector de condicion.")
}

# Usamos stan_data_phi_cond como referencia de condicion para todos los modelos.
ref_dims <- get_stan_dims(stan_data_phi_cond)
condition_vector <- get_condition_vector(stan_data_phi_cond, ref_dims$N)
levels_used <- choose_two_levels(condition_vector, control_level, tratamiento_level)
control_level <- levels_used$control
tratamiento_level <- levels_used$tratamiento

message("Comparando condiciones: ", control_level, " vs ", tratamiento_level)
message("Divergencia: ", ifelse(symmetric, "KL simetrica", "KL(Control||Tratamiento)"))

#-------------------------------------------------------------
# 3. Asociar cada fit con su stan_data correspondiente
#-------------------------------------------------------------

stan_data_by_model <- list(
  phi_comun = stan_data_phi_comun,
  phi_cond = stan_data_phi_cond
)

if (exists("stan_data_phi_comun_dia")) {
  stan_data_by_model$phi_comun_dia <- stan_data_phi_comun_dia
}
if (exists("stan_data_phi_cond_dia")) {
  stan_data_by_model$phi_cond_dia <- stan_data_phi_cond_dia
}

model_names <- intersect(names(fits_dirich), names(stan_data_by_model))
if (length(model_names) == 0) {
  stop("No hay modelos comunes entre fits_dirich y stan_data_by_model.")
}

#-------------------------------------------------------------
# 4. Calcular KL en cada modelo
#-------------------------------------------------------------

kl_results <- list()
summary_rows <- list()

a <- 1
for (model_name in model_names) {
  message("Procesando modelo: ", model_name)

  res <- calc_kl_for_model(
    fit = fits_dirich[[model_name]],
    stan_data = stan_data_by_model[[model_name]],
    condition_vector = condition_vector,
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
summary_file <- file.path(output_dir, paste0("resumen_", suffix, "_dirichlet.csv"))
rdata_file <- file.path(output_dir, paste0("draws_", suffix, "_dirichlet.Rdata"))

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
