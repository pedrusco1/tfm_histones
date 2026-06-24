###############################################################
# Aplicar divergencia KL continua entre las dos distribuciones 
# normales multivariantes (Espacio ILR) para los cuatro ajustes Stan revisados
# H3917: Modelos Log-normal revisados con coeficientes CLR de suma cero
###############################################################

#-------------------------------------------------------------
# 0. Configuración
#-------------------------------------------------------------

library(rstan)
library(compositions) # Requerido para ilrBase y clrInv

if (file.exists(here::here("ScriptsR", "40_stan_config.R"))) {
  source(here::here("ScriptsR", "40_stan_config.R"))
}

if (!exists("data_dir")) {
  data_dir <- here::here("DatosProcesados")
}

output_dir <- file.path(data_dir, "resultados_kl_lognormal_sum0_revisados")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_dir, "stan_lognormal_sum0_fits_revisados.Rdata")

control_level <- NULL
tratamiento_level <- NULL
symmetric <- TRUE
save_draws <- TRUE

#-------------------------------------------------------------
# 1. Funciones auxiliares para KL Continua Multivariante
#-------------------------------------------------------------

# Calcula la KL unidimensional continua entre N(mu1, Sigma1) y N(mu2, Sigma2)
# utilizando directamente los factores de Cholesky L_S1 y L_S2 para máxima estabilidad.
kl_mvn_cholesky_unidirectional <- function(mu1, L_S1, mu2, L_S2) {
  d <- length(mu1)
  
  # Reconstruir matrices de covarianza completas a partir de Cholesky
  Sigma1 <- L_S1 %*% t(L_S1)
  Sigma2 <- L_S2 %*% t(L_S2)
  
  # Inversa de Sigma2
  inv_Sigma2 <- solve(Sigma2)
  
  # Componentes de la fórmula clásica de la KL para Normales Multivariantes
  term_trace <- sum(diag(inv_Sigma2 %*% Sigma1))
  term_mean  <- as.numeric(t(mu2 - mu1) %*% inv_Sigma2 %*% (mu2 - mu1))
  
  # El log-determinante calculado de forma estable usando la diagonal de Cholesky
  log_det1 <- 2 * sum(log(abs(diag(L_S1))))
  log_det2 <- 2 * sum(log(abs(diag(L_S2))))
  term_det  <- log_det2 - log_det1
  
  # Combinación final
  0.5 * (term_trace + term_mean - d + term_det)
}

calc_kl_for_lognormal_model <- function(fit,
                                        X_mat,
                                        condition_vector,
                                        control_level,
                                        tratamiento_level,
                                        symmetric = TRUE) {
  
  idx_control <- which(condition_vector == control_level)
  idx_trat <- which(condition_vector == tratamiento_level)
  
  if (length(idx_control) == 0) stop("No hay observaciones para control_level = ", control_level)
  if (length(idx_trat) == 0) stop("No hay observaciones para tratamiento_level = ", tratamiento_level)
  
  post <- rstan::extract(fit)
  if (!"B" %in% names(post)) {
    stop("No se encuentra el parámetro matricial 'B' (CLR) en el objeto stanfit.")
  }
  if (!"L_Sigma" %in% names(post)) {
    stop("No se encuentra el parámetro 'L_Sigma' (Factor Cholesky ILR) en el objeto stanfit.")
  }
  
  B_array       <- post$B
  L_Sigma_array <- post$L_Sigma
  
  ndraws <- dim(B_array)[1]
  K      <- dim(B_array)[2]
  P      <- dim(B_array)[3]
  K_ilr  <- K - 1
  
  # Recuperar exactamente la base ILR utilizada al ajustar este modelo.
  # Evita asumir que compositions::ilrBase() coincide con la base guardada en Stan.
  stan_data_fit <- fit@stan_args[[1]]$data
  if (is.null(stan_data_fit$V)) {
    stop("El objeto stanfit no contiene la matriz V usada en el ajuste.")
  }
  V <- as.matrix(stan_data_fit$V)
  if (!all(dim(V) == c(K, K_ilr))) {
    stop("Dimensiones incompatibles de V en el objeto stanfit.")
  }
  
  X_mean_control <- colMeans(X_mat[idx_control, , drop = FALSE])
  X_mean_trat    <- colMeans(X_mat[idx_trat, , drop = FALSE])
  
  kl_draws <- numeric(ndraws)
  
  # Conservamos las matrices de proporciones por compatibilidad y análisis posterior
  prop_control_mat <- matrix(NA_real_, nrow = ndraws, ncol = K)
  prop_trat_mat    <- matrix(NA_real_, nrow = ndraws, ncol = K)
  
  for (i in seq_len(ndraws)) {
    # 1. Extraer efectos y calcular centros en espacio abierto CLR
    B_i <- B_array[i, , , drop = FALSE]
    dim(B_i) <- c(K, P) 
    
    clr_control <- as.numeric(B_i %*% X_mean_control)
    clr_trat    <- as.numeric(B_i %*% X_mean_trat)
    
    # Guardar las proporciones en el simplex (clrInv) para históricos/gráficos
    prop_control_mat[i, ] <- as.numeric(compositions::clrInv(clr_control))
    prop_trat_mat[i, ]    <- as.numeric(compositions::clrInv(clr_trat))
    
    # 2. PROYECTAR MEDIAS AL ESPACIO ILR (no singular) usando la base V
    ilr_control <- as.numeric(t(V) %*% clr_control)
    ilr_trat    <- as.numeric(t(V) %*% clr_trat)
    
    # 3. Extraer matriz Cholesky de la covarianza (dimensión K_ilr x K_ilr)
    L_S_i <- L_Sigma_array[i, , , drop = FALSE]
    dim(L_S_i) <- c(K_ilr, K_ilr)
    
    # 4. Calcular Divergencia KL Continua Multivariante
    kl12 <- kl_mvn_cholesky_unidirectional(ilr_control, L_S_i, ilr_trat, L_S_i)
    
    if (isTRUE(symmetric)) {
      kl21 <- kl_mvn_cholesky_unidirectional(ilr_trat, L_S_i, ilr_control, L_S_i)
      kl_draws[i] <- 0.5 * (kl12 + kl21)
    } else {
      kl_draws[i] <- kl12
    }
  }
  
  list(
    kl_draws = kl_draws,
    prop_control = prop_control_mat,
    prop_tratamiento = prop_trat_mat,
    n_control = length(idx_control),
    n_tratamiento = length(idx_trat)
  )
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
# 2. Cargar datos base y resultados Stan Log-normal
#-------------------------------------------------------------

load(file.path(data_dir, "data_inputs_h3917.Rdata"))
condition_vector <- cond

if (!file.exists(input_file)) {
  stop("No encuentro el fichero de modelos log-normales de suma cero: ", input_file)
}
load(input_file)

# Usar la matriz de diseño guardada con los ajustes revisados.
if (exists("X_mat_revisado")) {
  X_mat <- X_mat_revisado
} else if (!exists("X_mat")) {
  message("-> No se detectó X_mat_revisado; reconstruyendo la matriz de diseño.")
  if (exists("df.ancho")) {
    X_mat <- model.matrix(~ Estado, data = df.ancho)
  } else {
    X_mat <- model.matrix(~ factor(condition_vector))
  }
}

levels_used <- choose_two_levels(condition_vector, control_level, tratamiento_level)
control_level <- levels_used$control
tratamiento_level <- levels_used$tratamiento

message("Comparando condiciones (Densidades continuas MVN): ", control_level, " vs ", tratamiento_level)

#-------------------------------------------------------------
# 3. Preparar lista de modelos lognormal_sum0 a procesar
#-------------------------------------------------------------

fits_lognormal <- list()
if (exists("fit_lognor_v1_sum0coef_revisado")) {
  fits_lognormal$lognor_v1_sum0coef_revisado <- fit_lognor_v1_sum0coef_revisado
}
if (exists("fit_lognor_v1_dia_sum0coef_revisado")) {
  fits_lognormal$lognor_v1_dia_sum0coef_revisado <- fit_lognor_v1_dia_sum0coef_revisado
}
if (exists("fit_lognor_v1_sum0coef_delta_comun_revisado")) {
  fits_lognormal$lognor_v1_sum0coef_delta_comun_revisado <- fit_lognor_v1_sum0coef_delta_comun_revisado
}
if (exists("fit_lognor_v1_dia_sum0coef_delta_comun_revisado")) {
  fits_lognormal$lognor_v1_dia_sum0coef_delta_comun_revisado <- fit_lognor_v1_dia_sum0coef_delta_comun_revisado
}

if (length(fits_lognormal) == 0) {
  stop("No se han encontrado los objetos fit_lognor_v1_* en el entorno.")
}

#-------------------------------------------------------------
# 4. Calcular KL continua en cada uno de los 4 modelos
#-------------------------------------------------------------

kl_results <- list()
summary_rows <- list()

a <- 1
for (model_name in names(fits_lognormal)) {
  message("Procesando modelo log-normal: ", model_name)
  
  res <- calc_kl_for_lognormal_model(
    fit = fits_lognormal[[model_name]],
    X_mat = X_mat,
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
# 5. Guardar resultados consolidados
#-------------------------------------------------------------

suffix <- ifelse(symmetric, "kl_mvn_simetrica", "kl_mvn_control_vs_tratamiento")
summary_file <- file.path(output_dir, paste0("resumen_", suffix, "_lognormal_sum0_revisados.csv"))
rdata_file <- file.path(output_dir, paste0("draws_", suffix, "_lognormal_sum0.Rdata"))

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
message("\nResumen de KL Continua Normal Multivariante guardado correctamente.")

