#-------Funcion para el analisis de prior de a

analizar_formas_cuadraticas <- function(a_sim, Matriz_Delta, s) {
  library(MASS)
  
  # 1. Preparación inicial
  n_comp <- ncol(Matriz_Delta)
  mu <- rep(0, n_comp)
  cov_base <- ginv(Matriz_Delta)
  k <- qr(Matriz_Delta)$rank
  
  # 2. Generación de muestras (Iterando sobre el vector a_sim)
  lista_muestras <- lapply(a_sim, function(ai) {
    # Sigma específica para cada valor de 'a'
    Sigma_i <- (s / ai) * cov_base
    # Forzar simetría para evitar errores numéricos de eigen
    Sigma_i <- (Sigma_i + t(Sigma_i)) / 2
    
    return(mvrnorm(n = 1, mu = mu, Sigma = Sigma_i, tol = 1e-6))
  })
  
  # 3. Construcción de matriz y cálculo de formas cuadráticas
  matriz_final <- do.call(rbind, lista_muestras)
  formas_cuadraticas <- rowSums((matriz_final %*% Matriz_Delta) * matriz_final)
  
  # 4. Cálculos solicitados
  # Verificación teórica: (a/s) * x' * Delta * x debería tender al rango k
  verificacion <- (a_sim / s) * formas_cuadraticas
  promedio_verificacion <- mean(verificacion)
  
  # Valor teórico esperado de las formas cuadráticas
  valor_esperado_teorico <- (s / mean(a_sim)) * k
  
  # Probabilidad P(formas < valor_esperado)
  probabilidad_inferior <- mean(formas_cuadraticas < valor_esperado_teorico)
  
  # 5. Retornar resultados como lista
  return(list(
    formas_cuadraticas = formas_cuadraticas,
    promedio_verificacion = promedio_verificacion,
    valor_esperado = valor_esperado_teorico,
    probabilidad_menor_esperado = probabilidad_inferior,
    rango_matriz = k
  ))
}