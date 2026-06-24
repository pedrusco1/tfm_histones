################################
# FUNCION PCC VALORES GLOBALES #
################################

library(compositions)
library(ggplot2)
library(gridExtra)
library(dplyr)

#' Función para PPC de modelos composicionales
#' @param y_obs Matriz o data.frame [N x K] con composiciones observadas.
#' @param y_rep Array [Iteraciones x N x K] con réplicas predictivas en el símplex.
#' @return Lista con p-valores y objeto ggplot.

check_compo_ppc <- function(y_obs, y_rep) {
  
  # --- Funciones internas ---
  calc_total_variation <- function(data_mat) {
    data_clr <- clr(acomp(data_mat + 1e-9))
    return(sum(apply(data_clr, 2, var)))
  }
  
  calc_aitchison_norm <- function(data_mat) {
    centroide <- mean(acomp(data_mat + 1e-9))
    centroide_clr <- as.numeric(clr(centroide))
    return(sqrt(sum(centroide_clr^2)))
  }
  
  # 1. ESTADÍSTICOS OBSERVADOS
  t_obs_norma <- calc_aitchison_norm(y_obs)
  t_obs_var   <- calc_total_variation(y_obs)
  
  # 2. PROCESAMIENTO DE RÉPLICAS
  n_iter <- dim(y_rep)[1]
  t_rep_norma <- numeric(n_iter)
  t_rep_var   <- numeric(n_iter)
  
  for(i in 1:n_iter) {
    y_rep_i <- y_rep[i, , ]
    t_rep_norma[i] <- calc_aitchison_norm(y_rep_i)
    t_rep_var[i]   <- calc_total_variation(y_rep_i)
  }
  
  # 3. CÁLCULO DE P-VALORES BAYESIANOS
  p_norma <- mean(t_rep_norma > t_obs_norma)
  p_var   <- mean(t_rep_var > t_obs_var)
  
  # 4. GENERACIÓN DE GRÁFICOS
  # Función interna para estética de histogramas
  plot_hist <- function(t_rep, t_obs, p_val, titulo, label_x, color_fill) {
    ggplot(data.frame(x = t_rep), aes(x)) +
      geom_histogram(fill = color_fill, color = "white", bins = 30, alpha = 0.7) +
      geom_vline(xintercept = t_obs, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = titulo,
           subtitle = paste("p-PPC:", round(p_val, 3)),
           x = label_x, y = "Frecuencia") +
      theme_minimal()
  }
  
  p1 <- plot_hist(t_rep_norma, t_obs_norma, p_norma, 
                  "Localización: Norma del Centroide", "||g(y)||clr", "#377EB8")
  
  p2 <- plot_hist(t_rep_var, t_obs_var, p_var, 
                  "Dispersión: Variación Total", "Variación Total CLR", "#4DAF4A")
  
  grafico_final <- gridExtra::arrangeGrob(p1, p2, ncol = 2)
  
  # 5. RETORNO DE RESULTADOS
  return(list(
    p_valor_centroide = p_norma,
    p_valor_variacion = p_var,
    plot = grafico_final
  ))
}
