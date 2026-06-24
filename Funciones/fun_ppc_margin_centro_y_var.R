###############################################
# FUNCIONES: PPC MARGINALES
#   1. Centro composicional marginal
#   2. Varianza CLR marginal
###############################################

library(compositions)
library(ggplot2)
library(scales)

# ==========================================================
# Función auxiliar común para construir el gráfico marginal
# ==========================================================

plot_ppc_marginal <- function(resumen_modelo,
                              title,
                              subtitle,
                              y_label,
                              log_y = FALSE) {

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
      subtitle = subtitle,
      x = "Marca de Histona",
      y = y_label
    ) +
    theme_minimal()

  if(log_y) {
    p <- p + scale_y_log10(
      labels = scales::label_percent(accuracy = 0.1)
    )
  }

  return(p)
}

# ==========================================================
# 1. PPC marginal del centro composicional
# ==========================================================

check_ppc_margin <- function(
    y_obs,
    y_rep,
    title = "PPC: Ajuste marginal de proteoformas"
) {

  epsilon <- 1e-9
  nombres_marcas <- colnames(y_obs)
  K <- ncol(y_obs)
  n_iter <- dim(y_rep)[1]

  centroide_obs <- as.numeric(
    mean(acomp(y_obs + epsilon))
  )

  centroides_rep <- matrix(
    NA,
    nrow = n_iter,
    ncol = K
  )

  for(i in seq_len(n_iter)) {
    centroides_rep[i, ] <- as.numeric(
      mean(acomp(y_rep[i, , ] + epsilon))
    )
  }

  p_ppc <- colMeans(
    centroides_rep > matrix(
      rep(centroide_obs, n_iter),
      nrow = n_iter,
      byrow = TRUE
    )
  )

  resumen_modelo <- data.frame(
    Marca = nombres_marcas,
    Obs = centroide_obs,
    Media_Post = apply(centroides_rep, 2, mean),
    Low = apply(centroides_rep, 2, quantile, probs = 0.025),
    High = apply(centroides_rep, 2, quantile, probs = 0.975),
    p_PPC = p_ppc,
    stringsAsFactors = FALSE
  )

  p <- plot_ppc_marginal(
    resumen_modelo = resumen_modelo,
    title = title,
    subtitle = "Centroide observado (rojo) vs intervalo predictivo posterior 95% (azul)",
    y_label = "Abundancia relativa",
    log_y = TRUE
  )

  return(list(
    grafico = p,
    stats = resumen_modelo,
    stat_rep = centroides_rep
  ))
}

# ==========================================================
# 2. PPC marginal de la varianza CLR
# ==========================================================

check_ppc_margin_var_clr <- function(
    y_obs,
    y_rep,
    title = "PPC: Variación marginal CLR por proteoforma"
) {

  epsilon <- 1e-9
  nombres_marcas <- colnames(y_obs)
  K <- ncol(y_obs)
  n_iter <- dim(y_rep)[1]

  y_obs_clr <- clr(acomp(y_obs + epsilon))
  var_obs <- apply(y_obs_clr, 2, var)

  var_rep <- matrix(
    NA,
    nrow = n_iter,
    ncol = K
  )

  for(i in seq_len(n_iter)) {
    y_rep_clr_i <- clr(acomp(y_rep[i, , ] + epsilon))
    var_rep[i, ] <- apply(y_rep_clr_i, 2, var)
  }

  p_ppc <- colMeans(
    var_rep > matrix(
      rep(var_obs, n_iter),
      nrow = n_iter,
      byrow = TRUE
    )
  )

  resumen_modelo <- data.frame(
    Marca = nombres_marcas,
    Obs = var_obs,
    Media_Post = apply(var_rep, 2, mean),
    Low = apply(var_rep, 2, quantile, probs = 0.025),
    High = apply(var_rep, 2, quantile, probs = 0.975),
    p_PPC = p_ppc,
    stringsAsFactors = FALSE
  )

  p <- plot_ppc_marginal(
    resumen_modelo = resumen_modelo,
    title = title,
    subtitle = "Varianza CLR observada (rojo) vs intervalo predictivo posterior 95% (azul)",
    y_label = "Varianza de coordenada CLR",
    log_y = FALSE
  )

  return(list(
    grafico = p,
    stats = resumen_modelo,
    stat_rep = var_rep
  ))
}
