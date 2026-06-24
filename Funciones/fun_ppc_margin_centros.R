###############################################
# FUNCION: PPC MARGINALES. EVALUACION GRAFICA #
###############################################

library(compositions)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

library(compositions)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

check_ppc_margin <- function(y_obs, y_rep, title = "PPC: Ajuste Marginal de Proteoformas") {
  
  # 0. Preparación básica
  epsilon <- 1e-9
  nombres_marcas <- colnames(y_obs)
  K <- ncol(y_obs)
  n_iter <- dim(y_rep)[1]
  
  # 1. CÁLCULO DEL CENTROIDE OBSERVADO
  centroide_obs <- as.numeric(mean(acomp(y_obs + epsilon)))
  
  # 2. PROCESAMIENTO DE LAS RÉPLICAS DEL MODELO
  centroides_rep <- matrix(NA, nrow = n_iter, ncol = K)
  for(i in 1:n_iter) {
    centroides_rep[i, ] <- as.numeric(mean(acomp(y_rep[i, , ] + epsilon)))
  }
  
  # 3. CÁLCULO DE p-values PPC POR PROTEOFORMA
  # p-PPC = Pr(Centroide_rep > Centroide_obs)
  p_ppc <- colMeans(centroides_rep > matrix(rep(centroide_obs, n_iter), 
                                            nrow = n_iter, byrow = TRUE))
  
  # 4. TABLA DE RESUMEN
  resumen_modelo <- data.frame(
    Marca = nombres_marcas,
    Obs = centroide_obs,
    Media_Post = apply(centroides_rep, 2, mean),
    Low = apply(centroides_rep, 2, quantile, probs = 0.025),
    High = apply(centroides_rep, 2, quantile, probs = 0.975),
    p_PPC = p_ppc
  )
  
  # 5. GENERACIÓN DEL GRÁFICO
  p <- ggplot(resumen_modelo, aes(x = reorder(Marca, Obs))) +
    geom_linerange(aes(ymin = Low, ymax = High), 
                   color = "#377EB8", linewidth = 1.2, alpha = 0.5) +
    geom_point(aes(y = Media_Post), color = "#377EB8", size = 2, shape = 1) +
    geom_point(aes(y = Obs), color = "red", size = 2.5) +
    scale_y_log10(labels = scales::label_percent(accuracy = 0.1)) +
    coord_flip() +
    labs(
      title = title,
      subtitle = "Centroide observado (rojo) vs Intervalo de Credibilidad 95% (azul)",
      x = "Marca de Histona",
      y = "Abundancia Relativa (Log)"
    ) +
    theme_minimal()
  
  # 6. SALIDA COMO LISTA
  return(list(
    grafico = p,
    stats = resumen_modelo
  ))
}

check_ppc_margin_var_clr <- function(
    y_obs,
    y_rep,
    title = "PPC: Variación marginal CLR por proteoforma"
) {
  
  epsilon <- 1e-9
  nombres_marcas <- colnames(y_obs)
  K <- ncol(y_obs)
  n_iter <- dim(y_rep)[1]
  
  # 1. Varianza CLR observada por parte
  y_obs_clr <- clr(acomp(y_obs + epsilon))
  var_obs <- apply(y_obs_clr, 2, var)
  
  # 2. Varianza CLR en cada réplica posterior
  var_rep <- matrix(NA, nrow = n_iter, ncol = K)
  
  for(i in seq_len(n_iter)) {
    y_rep_clr_i <- clr(acomp(y_rep[i, , ] + epsilon))
    var_rep[i, ] <- apply(y_rep_clr_i, 2, var)
  }
  
  # 3. p-values PPC por proteoforma
  p_ppc <- colMeans(
    var_rep > matrix(rep(var_obs, n_iter),
                     nrow = n_iter,
                     byrow = TRUE)
  )
  
  # 4. Tabla resumen
  resumen_modelo <- data.frame(
    Marca = nombres_marcas,
    Obs = var_obs,
    Media_Post = apply(var_rep, 2, mean),
    Low = apply(var_rep, 2, quantile, probs = 0.025),
    High = apply(var_rep, 2, quantile, probs = 0.975),
    p_PPC = p_ppc
  )
  
  # 5. Gráfico
  p <- ggplot(resumen_modelo, aes(x = reorder(Marca, Obs))) +
    geom_linerange(
      aes(ymin = Low, ymax = High),
      color = "#377EB8",
      linewidth = 1.2,
      alpha = 0.5
    ) +
    geom_point(
      aes(y = Media_Post),
      color = "#377EB8",
      size = 2,
      shape = 1
    ) +
    geom_point(
      aes(y = Obs),
      color = "red",
      size = 2.5
    ) +
    coord_flip() +
    labs(
      title = title,
      subtitle = "Varianza CLR observada (rojo) vs intervalo predictivo posterior 95% (azul)",
      x = "Marca de Histona",
      y = "Varianza de coordenada CLR"
    ) +
    theme_minimal()
  
  return(list(
    grafico = p,
    stats = resumen_modelo,
    var_rep = var_rep
  ))
}
