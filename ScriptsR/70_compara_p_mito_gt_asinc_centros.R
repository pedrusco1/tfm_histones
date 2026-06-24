###########################################################################
# SCRIPT DE INFERENCIA POSTERIOR PREDICTIVA:
# P(centro composicional predictivo mitótico >
#   centro composicional predictivo asíncrono)
#
# Usa Y_rep de generated quantities.
#
# Centro composicional:
#   centro(Y_1, ..., Y_m) = cierre(exp(media(log(Y_i))))
#
# Observaciones:
#   asíncrono: 1, 3, 5
#   mitosis:   2, 4, 6
#
# Procesamiento por bloques y limpieza estricta de RAM.
###########################################################################

rm(list = ls())
gc()

library(rstan)

# -------------------------------------------------------------------------
# 1. Configuración de directorios y carga de nombres de proteoformas
# -------------------------------------------------------------------------

data_dir <- "DatosProcesados"

load(file.path(data_dir, "df_h3917.Rdata"))

if (!exists("df.ancho")) {
  stop("Error: No se encuentra el objeto `df.ancho` en el archivo de entrada.")
}

proteoformas <- colnames(df.ancho)[3:22]
K_total <- length(proteoformas)

cat(">>> Inicializando análisis predictivo posterior para",
    K_total, "proteoformas.\n")

tabla_resultados <- data.frame(
  Proteoforma = proteoformas,
  stringsAsFactors = FALSE
)

# Observaciones según el orden de df.ancho:
#   asíncrono: 1, 3, 5
#   mitosis:   2, 4, 6
idx_asinc <- c(1, 3, 5)
idx_mitot <- c(2, 4, 6)

# -------------------------------------------------------------------------
# 2. Funciones auxiliares
# -------------------------------------------------------------------------

# Cierre composicional
closure_vec <- function(x) {
  sx <- sum(x)
  
  if (!is.finite(sx) || sx <= 0) {
    stop("No se puede cerrar un vector con suma no positiva o no finita.")
  }
  
  x / sx
}

# Centro composicional de una matriz de composiciones
# Filas = observaciones
# Columnas = partes/proteoformas
centro_composicional <- function(Y, eps = 1e-12) {
  
  # Evitar log(0), aunque Y_rep debería ser estrictamente positivo
  Y <- pmax(Y, eps)
  
  # Media geométrica por parte:
  # exp(media(log(parte)))
  g <- exp(colMeans(log(Y)))
  
  # Cierre
  closure_vec(g)
}

# -------------------------------------------------------------------------
# 3. Función principal:
#    P(centro composicional mitótico > centro composicional asíncrono)
#    desde Y_rep
# -------------------------------------------------------------------------

calc_prob_centros_yrep <- function(fit_obj,
                                   K,
                                   proteoformas,
                                   idx_asinc = c(1, 3, 5),
                                   idx_mitot = c(2, 4, 6),
                                   par_yrep = "Y_rep",
                                   eps = 1e-12) {
  
  if (!par_yrep %in% fit_obj@model_pars) {
    stop("El objeto stanfit no contiene el parámetro `", par_yrep, "`.")
  }
  
  message("   Extrayendo únicamente ", par_yrep, "...")
  
  # Extrae únicamente Y_rep.
  # Matriz plana:
  #   filas = draws posteriores
  #   columnas = Y_rep[n, k]
  Y_mat <- as.matrix(fit_obj, pars = par_yrep)
  
  n_draws <- nrow(Y_mat)
  
  centro_asinc <- matrix(NA_real_, nrow = n_draws, ncol = K)
  centro_mitot <- matrix(NA_real_, nrow = n_draws, ncol = K)
  
  colnames(centro_asinc) <- proteoformas
  colnames(centro_mitot) <- proteoformas
  
  # Preconstruir nombres de columnas para evitar repetir comprobaciones
  cols_asinc_por_k <- vector("list", K)
  cols_mitot_por_k <- vector("list", K)
  
  for (k in seq_len(K)) {
    cols_asinc <- paste0(par_yrep, "[", idx_asinc, ",", k, "]")
    cols_mitot <- paste0(par_yrep, "[", idx_mitot, ",", k, "]")
    
    faltan_asinc <- setdiff(cols_asinc, colnames(Y_mat))
    faltan_mitot <- setdiff(cols_mitot, colnames(Y_mat))
    
    if (length(faltan_asinc) > 0) {
      stop(
        "Faltan columnas asíncronas para la proteoforma ", k, ": ",
        paste(faltan_asinc, collapse = ", ")
      )
    }
    
    if (length(faltan_mitot) > 0) {
      stop(
        "Faltan columnas mitóticas para la proteoforma ", k, ": ",
        paste(faltan_mitot, collapse = ", ")
      )
    }
    
    cols_asinc_por_k[[k]] <- cols_asinc
    cols_mitot_por_k[[k]] <- cols_mitot
  }
  
  message("   Calculando centros composicionales draw a draw...")
  
  for (s in seq_len(n_draws)) {
    
    Y_asinc_s <- matrix(NA_real_, nrow = length(idx_asinc), ncol = K)
    Y_mitot_s <- matrix(NA_real_, nrow = length(idx_mitot), ncol = K)
    
    for (k in seq_len(K)) {
      Y_asinc_s[, k] <- as.numeric(Y_mat[s, cols_asinc_por_k[[k]]])
      Y_mitot_s[, k] <- as.numeric(Y_mat[s, cols_mitot_por_k[[k]]])
    }
    
    centro_asinc[s, ] <- centro_composicional(Y_asinc_s, eps = eps)
    centro_mitot[s, ] <- centro_composicional(Y_mitot_s, eps = eps)
  }
  
  # Probabilidad posterior predictiva:
  # P(centro composicional mitosis_k > centro composicional asinc_k)
  prob <- colMeans(centro_mitot > centro_asinc)
  names(prob) <- proteoformas
  
  # Limpieza local
  rm(
    Y_mat,
    centro_asinc,
    centro_mitot,
    cols_asinc_por_k,
    cols_mitot_por_k
  )
  gc()
  
  prob
}

# -------------------------------------------------------------------------
# 4. Función para procesar un bloque de modelos desde un fichero .Rdata
# -------------------------------------------------------------------------

procesar_bloque_yrep <- function(file_rdata,
                                 modelos,
                                 nombre_bloque,
                                 tabla_resultados,
                                 K,
                                 proteoformas,
                                 idx_asinc = c(1, 3, 5),
                                 idx_mitot = c(2, 4, 6)) {
  
  if (!file.exists(file_rdata)) {
    warning("No se localizó el fichero: ", file_rdata)
    
    for (col_name in names(modelos)) {
      tabla_resultados[[col_name]] <- NA_real_
    }
    
    return(tabla_resultados)
  }
  
  message("\n=========================================================")
  message(">>> PROCESANDO BLOQUE: ", nombre_bloque)
  message(">>> Usando centros composicionales de Y_rep")
  message("=========================================================")
  
  load(file_rdata)
  
  for (col_name in names(modelos)) {
    
    obj_name <- modelos[[col_name]]
    
    if (!exists(obj_name)) {
      warning("No se encontró el objeto ", obj_name, " en ", file_rdata)
      tabla_resultados[[col_name]] <- NA_real_
      next
    }
    
    message("-> Modelo: ", obj_name)
    
    fit_obj <- get(obj_name)
    
    tabla_resultados[[col_name]] <- calc_prob_centros_yrep(
      fit_obj = fit_obj,
      K = K,
      proteoformas = proteoformas,
      idx_asinc = idx_asinc,
      idx_mitot = idx_mitot,
      par_yrep = "Y_rep"
    )
    
    # Borrar inmediatamente el ajuste y objetos temporales del modelo
    rm(fit_obj)
    rm(list = obj_name)
    gc()
  }
  
  # Limpieza adicional de objetos que puedan venir en el .Rdata.
  # Conservamos solo lo necesario para seguir creando la tabla.
  objetos_a_conservar <- c(
    "tabla_resultados",
    "proteoformas",
    "K_total",
    "idx_asinc",
    "idx_mitot",
    "data_dir",
    "df.ancho",
    "closure_vec",
    "centro_composicional",
    "calc_prob_centros_yrep",
    "procesar_bloque_yrep"
  )
  
  rm(
    list = setdiff(ls(envir = .GlobalEnv), objetos_a_conservar),
    envir = .GlobalEnv
  )
  gc()
  
  tabla_resultados
}

# -------------------------------------------------------------------------
# 5. BLOQUE 1: Modelos Log-Normal
# -------------------------------------------------------------------------

file_lognor <- file.path(data_dir, "stan_lognormal_sum0_fits.Rdata")

modelos_lognor <- list(
  Lognormal_V1_SinDia             = "fit_lognor_v1_sum0coef",
  Lognormal_V1_ConDia             = "fit_lognor_v1_dia_sum0coef",
  Lognormal_DeltaComun_SinDia     = "fit_lognor_v1_sum0coef_delta_comun",
  Lognormal_DeltaComun_ConDia     = "fit_lognor_v1_dia_sum0coef_delta_comun"
)

tabla_resultados <- procesar_bloque_yrep(
  file_rdata = file_lognor,
  modelos = modelos_lognor,
  nombre_bloque = "LOG-NORMAL",
  tabla_resultados = tabla_resultados,
  K = K_total,
  proteoformas = proteoformas,
  idx_asinc = idx_asinc,
  idx_mitot = idx_mitot
)

# -------------------------------------------------------------------------
# 6. BLOQUE 2: Modelos de Epistasia
# -------------------------------------------------------------------------

file_epistasia <- file.path(data_dir, "stan_epistasia_fits.Rdata")

modelos_epistasia <- list(
  Epistasia_KappaComun_SinDia = "fit_kappa_comun",
  Epistasia_KappaComun_ConDia = "fit_kappa_comun_dia",
  Epistasia_KappaCond_SinDia  = "fit_kappa_cond",
  Epistasia_KappaCond_ConDia  = "fit_kappa_cond_dia"
)

tabla_resultados <- procesar_bloque_yrep(
  file_rdata = file_epistasia,
  modelos = modelos_epistasia,
  nombre_bloque = "EPISTASIA KAPPA",
  tabla_resultados = tabla_resultados,
  K = K_total,
  proteoformas = proteoformas,
  idx_asinc = idx_asinc,
  idx_mitot = idx_mitot
)

# -------------------------------------------------------------------------
# 7. BLOQUE 3: Modelos Dirichlet
# -------------------------------------------------------------------------

file_dirich <- file.path(data_dir, "stan_dirichlet_fits.Rdata")

modelos_dirich <- list(
  Dirichlet_Comun_SinDia = "fit_dirich_phi_comun",
  Dirichlet_Cond_SinDia  = "fit_dirich_phi_cond",
  Dirichlet_Comun_ConDia = "fit_dirich_phi_comun_dia",
  Dirichlet_Cond_ConDia  = "fit_dirich_phi_cond_dia"
)

tabla_resultados <- procesar_bloque_yrep(
  file_rdata = file_dirich,
  modelos = modelos_dirich,
  nombre_bloque = "REGRESIÓN DIRICHLET",
  tabla_resultados = tabla_resultados,
  K = K_total,
  proteoformas = proteoformas,
  idx_asinc = idx_asinc,
  idx_mitot = idx_mitot
)

# -------------------------------------------------------------------------
# 8. Guardado de resultados
# -------------------------------------------------------------------------

message("\n=========================================================")
message(">>> ANÁLISIS FINALIZADO")
message(">>> Probabilidades basadas en centros composicionales de Y_rep")
message("=========================================================")

print(head(tabla_resultados))

write.csv(
  tabla_resultados,
  file = file.path(data_dir, "tabla_probabilidades_centros_composicionales_yrep_tfm.csv"),
  row.names = FALSE
)

save(
  tabla_resultados,
  file = file.path(data_dir, "tabla_probabilidades_centros_composicionales_yrep_tfm.Rdata")
)

cat("\n[✓] Fichero generado: tabla_probabilidades_centros_composicionales_yrep_tfm.csv\n")
cat("[✓] Fichero generado: tabla_probabilidades_centros_composicionales_yrep_tfm.Rdata\n")
