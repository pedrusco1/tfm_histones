###########################################################################
# FUNCIÓN: PPC PARA BALANCES BIOLÓGICOS (Localización y Varianza)
###########################################################################

library(ggplot2)
library(gridExtra)
library(dplyr)

check_balances_ppc <- function(y_obs, y_rep) {
  
  n_obs <- nrow(y_obs)
  n_iter <- dim(y_rep)[1]
  
  # 1. DEFINICIÓN DE ÍNDICES BIOLÓGICOS -----------------------------------
  idx_S10ph   <- c(10:16)                   
  idx_noS10ph <- setdiff(1:20, idx_S10ph)
  idx_mod     <- c(2:20)                    
  idx_unmod   <- c(1)
  idx_ac      <- c(5, 8, 9, 12, 13, 16, 19) 
  idx_no_ac   <- setdiff(1:20, idx_ac)
  
  list_balances <- list(
    list(name = "S10ph vs No-S10ph", pos = idx_S10ph, neg = idx_noS10ph),
    list(name = "Mod vs Unmod", pos = idx_mod, neg = idx_unmod),
    list(name = "Ac vs No-Ac", pos = idx_ac, neg = idx_no_ac)
  )
  
  # 2. FUNCIÓN INTERNA: CÁLCULO DE BALANCE --------------------------------
  calc_bal <- function(Y, pos, neg) {
    # Media geométrica de las partes (espacio composicional)
    gm_pos <- apply(Y[, pos, drop=FALSE], 1, function(x) exp(mean(log(x + 1e-9))))
    gm_neg <- apply(Y[, neg, drop=FALSE], 1, function(x) exp(mean(log(x + 1e-9))))
    return(log(gm_pos / gm_neg))
  }
  
  # Contenedores de resultados
  p_vals_mean <- list()
  p_vals_var  <- list()
  plots_mean  <- list()
  plots_var   <- list()
  
  # 3. PROCESAMIENTO POR BALANCE ------------------------------------------
  for(b in list_balances) {
    # A. Datos Observados
    bal_obs <- calc_bal(y_obs, b$pos, b$neg)
    stat_obs_mean <- mean(bal_obs)
    stat_obs_var  <- var(bal_obs)
    
    # B. Datos Replicados (Simulaciones del modelo)
    stat_rep_mean <- numeric(n_iter)
    stat_rep_var  <- numeric(n_iter)
    
    for(i in 1:n_iter) {
      bal_rep_i <- calc_bal(y_rep[i,,], b$pos, b$neg)
      stat_rep_mean[i] <- mean(bal_rep_i)
      stat_rep_var[i]  <- var(bal_rep_i)
    }
    
    # C. P-valores Bayesiano
    p_mean <- mean(stat_rep_mean > stat_obs_mean)
    p_var  <- mean(stat_rep_var > stat_obs_var)
    
    p_vals_mean[[b$name]] <- p_mean
    p_vals_var[[b$name]]  <- p_var
    
    # D. Generación de Gráficos (Localización)
    plots_mean[[b$name]] <- ggplot(data.frame(x = stat_rep_mean), aes(x)) +
      geom_histogram(fill = "#377EB8", color = "white", alpha = 0.7, bins = 35) +
      geom_vline(xintercept = stat_obs_mean, color = "red", linetype = "dashed", size = 1) +
      labs(title = paste("Media:", b$name), subtitle = paste("p =", round(p_mean, 3)),
           x = "Log-Ratio", y = "Frecuencia") +
      theme_minimal(base_size = 9)
    
    # E. Generación de Gráficos (Dispersión/Varianza)
    plots_var[[b$name]] <- ggplot(data.frame(x = stat_rep_var), aes(x)) +
      geom_histogram(fill = "#4DAF4A", color = "white", alpha = 0.7, bins = 35) +
      geom_vline(xintercept = stat_obs_var, color = "red", linetype = "dashed", size = 1) +
      labs(title = paste("Var:", b$name), subtitle = paste("p =", round(p_var, 3)),
           x = "Varianza del Balance", y = "Frecuencia") +
      theme_minimal(base_size = 9)
  }
  
  # 4. ORGANIZACIÓN DE SALIDA ---------------------------------------------
  return(list(
    p_pcc_medias     = p_vals_mean,
    p_pcc_varianzas  = p_vals_var,
    plot_localizacion = gridExtra::arrangeGrob(grobs = plots_mean, ncol = 3),
    plot_dispersion   = gridExtra::arrangeGrob(grobs = plots_var, ncol = 3)
  ))
}
