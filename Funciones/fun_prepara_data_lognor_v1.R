# ============================================================
# Preparar datos para modelos epistasia log-normal v1
# Compatible con:
#   1) epistasia_log_normal_v1.stan
#   2) epistasia_log_normal_v1_dia.stan
#
# Si se proporciona `day`, la lista incluye D y day para el modelo
# con efecto aleatorio de día. Si `day = NULL`, devuelve la lista
# original compatible con el modelo sin efecto día.
# ============================================================

make_stan_data <- function(Y, L, X, day = NULL) {
  
  if (!requireNamespace("compositions", quietly = TRUE)) {
    stop("Package 'compositions' required")
  }
  
  #========================
  # 1. Convertir a matrices
  #========================
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  L <- as.matrix(L)
  
  N <- nrow(Y)
  K <- ncol(Y)
  
  if (is.null(N) || is.null(K)) {
    stop("Y must be a matrix-like object with rows = samples and columns = compositional parts")
  }
  
  if (nrow(X) != N) {
    stop("X must have same number of rows as Y")
  }
  
  if (nrow(L) != K || ncol(L) != K) {
    stop("L must be K x K where K = ncol(Y)")
  }
  
  if (any(!is.finite(Y))) {
    stop("Y contains non-finite values")
  }
  
  if (any(Y <= 0)) {
    stop(
      "Y must contain strictly positive values. ",
      "The v1 log-normal model computes log(Y) in generated quantities, ",
      "so zeros must be handled before calling this function."
    )
  }
  
  if (any(rowSums(Y) <= 0)) {
    stop("All rows of Y must have positive sum")
  }
  
  P <- ncol(X)
  K_ilr <- K - 1
  
  #========================
  # 2. Cerrar composiciones
  #========================
  
  Y <- Y / rowSums(Y)
  
  #========================
  # 3. Base ILR
  #========================
  
  V <- compositions::ilrBase(compositions::acomp(rep(1, K)))
  V <- as.matrix(V)
  
  #========================
  # 4. Transformación ILR
  #========================
  
  Z <- compositions::ilr(compositions::acomp(Y))
  Z <- as.matrix(Z)
  
  #========================
  # 5. Laplaciano en ILR
  #========================
  
  L_ilr <- t(V) %*% L %*% V
  L_shape <- L_ilr / mean(diag(L_ilr))
  
  #========================
  # 6. Lista base para Stan
  #========================
  
  stan_data <- list(
    N = N,
    K = K,
    K_ilr = K_ilr,
    P = P,
    Z = Z,
    X = X,
    L_shape = L_shape,
    V = V,
    Y = Y
  )
  
  #========================
  # 7. Añadir efecto día, si procede
  #========================
  
  if (!is.null(day)) {
    if (length(day) != N) {
      stop("day must have length N = nrow(Y)")
    }
    
    if (any(is.na(day))) {
      stop("day contains NA values")
    }
    
    # Convertir días arbitrarios, por ejemplo 1/2/3, fechas o etiquetas,
    # a índices enteros consecutivos 1:D, como espera Stan.
    day_factor <- factor(day)
    day_index <- as.integer(day_factor)
    D <- nlevels(day_factor)
    
    stan_data$D <- D
    stan_data$day <- day_index
    
    # Guardar el mapeo como atributo para poder inspeccionarlo en R.
    # Stan ignorará los atributos al recibir la lista.
    attr(stan_data, "day_levels") <- levels(day_factor)
  }
  
  return(stan_data)
}
