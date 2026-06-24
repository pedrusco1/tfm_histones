#############################################################################
#              FUNCIONES CALCULO DIVERGENCIA JENSEN-TSALLIS                 #
#############################################################################

# ==========================================================================
# BLOQUE 1: Función Atómica (Cálculo puro)
# ==========================================================================
calc_jtd <- function(p, q, q_param = 0.5) {
  # Aseguramos que no haya ceros estrictos para evitar NaN
  eps <- 1e-10
  p <- (p + eps) / sum(p + eps)
  q <- (q + eps) / sum(q + eps)
  
  tsallis_s <- function(x, q_p) (1 - sum(x^q_p)) / (q_p - 1)
  
  m <- 0.5 * (p + q)
  val <- tsallis_s(m, q_param) - 0.5 * (tsallis_s(p, q_param) + tsallis_s(q, q_param))
  return(max(0, val))
}

# ==========================================================================
# BLOQUE 2: Función de Interfaz para PPC (Manejo de datos de Stan)
# ==========================================================================
apply_jtd_ppc <- function(y_obs, y_rep_draws, q_val = 0.5) {
  # y_obs: matriz [n_muestras_condicion, 20]
  # y_rep_draws: matriz [16000, 20] (ya filtrada por condición)
  
  # 1. Calculamos el centro observado de la condición
  centro_obs <- colMeans(y_obs)
  
  # 2. Aplicamos la función atómica a cada draw de la posterior
  # Usamos vapply para que sea más rápido y seguro en el tipo de dato
  t_rep <- vapply(1:nrow(y_rep_draws), function(s) {
    calc_jtd(centro_obs, y_rep_draws[s, ], q_param = q_val)
  }, FUN.VALUE = numeric(1))
  
  return(t_rep)
}
