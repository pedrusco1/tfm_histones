###############################################################
#   Preparar datos para regresión simetria compuesta          #
###############################################################

prepara_cosy_stan_data <- function(Y, X_cov, id_vector, V) {
  # Y: Matriz de composiciones (N x K)
  # X_cov: Data frame con covariables (ej. Estado)
  # id_vector: Vector que identifica cada muestra original
  # V: Matriz de base ILR (K x K-1) proporcionada por el usuario
  
  if(!requireNamespace("compositions", quietly = TRUE)) stop("Package 'compositions' required")
  
  # 1. Transformación a coordenadas ILR usando la V externa
  # Convertimos Y a objeto acomp para asegurar el cierre (sum 1)
  Y_acomp <- compositions::acomp(Y)
  
  # Calculamos Z: log(Y) %*% V (proyección en la base V)
  # Nota: usamos la multiplicación matricial manual para asegurar el uso de TU V
  # log_Y <- log(as.matrix(Y_acomp))
  # Z <- log_Y %*% V
  # Mas limpio es
  Z <- compositions::ilr(Y, V)
  
  colnames(Z) <- paste0("ILR", 1:ncol(Z))
  
  # 2. Estructurar datos en formato Long
  df_wide <- data.frame(id = id_vector, X_cov, Z)
  
  # Derretimos la matriz para tener una observación por fila (Requerido por cosy_log_normal.stan)
  df_long <- tidyr::pivot_longer(
    df_wide,
    cols = starts_with("ILR"),
    names_to = "Balance",
    values_to = "Valor"
  )
  
  # 3. Crear Matriz de Diseño (X) [cite: 3, 4]
  # El modelo cosy_log_normal.stan pide: mu = X * b + r_id [cite: 14]
  # Configuramos Balance e Interacción con Estado
  X_formula <- model.matrix(Valor ~ Balance + Balance:Estado - 1, data = df_long)
  
  # 4. Construir lista para Stan conforme al archivo cosy_log_normal.stan [cite: 1]
  stan_data <- list(
    N          = nrow(df_long),        # [cite: 2]
    Y          = df_long$Valor,        # [cite: 3]
    K          = ncol(X_formula),      # [cite: 3]
    X          = X_formula,            # [cite: 4]
    J_id       = length(unique(id_vector)), # [cite: 5]
    id         = as.integer(as.factor(df_long$id)) # [cite: 5]
  )
  
  return(stan_data)
}
