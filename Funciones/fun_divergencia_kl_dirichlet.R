#################################################################
# Divergencia KL posterior entre distribuciones Dirichlet
#
# Pensado para modelos Dirichlet donde:
#   alpha = mu * phi
#
# Permite comparar dos condiciones/grupos usando draws posteriores
# de mu y phi.
#################################################################

# ---------------------------------------------------------------
# KL entre dos Dirichlet individuales: KL(Dir(a1) || Dir(a2))
# ---------------------------------------------------------------
kl_dirichlet_alpha <- function(a1, a2) {
  if (length(a1) != length(a2)) {
    stop("a1 y a2 deben tener la misma longitud.")
  }
  
  if (any(a1 <= 0) || any(a2 <= 0)) {
    stop("Todos los parámetros alpha deben ser positivos.")
  }
  
  sum_a1 <- sum(a1)
  sum_a2 <- sum(a2)
  
  kl <- lgamma(sum_a1) - sum(lgamma(a1)) -
    lgamma(sum_a2) + sum(lgamma(a2)) +
    sum((a1 - a2) * (digamma(a1) - digamma(sum_a1)))
  
  as.numeric(kl)
}


# ---------------------------------------------------------------
# KL posterior draw a draw entre dos distribuciones Dirichlet
# ---------------------------------------------------------------
# Argumentos:
#   mu_control: matriz draws x K con medias composicionales posteriores
#   phi_control: vector de draws de phi para control
#   mu_tratamiento: matriz draws x K con medias composicionales posteriores
#   phi_tratamiento: vector de draws de phi para tratamiento
#   symmetric: si TRUE calcula 0.5 * (KL(control||tratamiento) + KL(tratamiento||control))
#   probs: cuantiles a devolver
#
# Devuelve:
#   lista con draws de KL, media, mediana, sd e intervalo posterior
# ---------------------------------------------------------------

divergencia_kl_dirichlet <- function(mu_control,
                                     phi_control,
                                     mu_tratamiento,
                                     phi_tratamiento,
                                     symmetric = FALSE,
                                     probs = c(0.025, 0.975)) {
  
  if (!is.matrix(mu_control)) {
    mu_control <- as.matrix(mu_control)
  }
  
  if (!is.matrix(mu_tratamiento)) {
    mu_tratamiento <- as.matrix(mu_tratamiento)
  }
  
  if (nrow(mu_control) != nrow(mu_tratamiento)) {
    stop("mu_control y mu_tratamiento deben tener el mismo número de draws.")
  }
  
  if (ncol(mu_control) != ncol(mu_tratamiento)) {
    stop("mu_control y mu_tratamiento deben tener el mismo número de categorías/composiciones.")
  }
  
  ndraws <- nrow(mu_control)
  
  if (length(phi_control) != ndraws) {
    stop("phi_control debe tener la misma longitud que el número de draws de mu_control.")
  }
  
  if (length(phi_tratamiento) != ndraws) {
    stop("phi_tratamiento debe tener la misma longitud que el número de draws de mu_tratamiento.")
  }
  
  if (any(mu_control <= 0) || any(mu_tratamiento <= 0)) {
    stop("Todos los valores de mu deben ser positivos. Revisa posibles ceros.")
  }
  
  if (any(phi_control <= 0) || any(phi_tratamiento <= 0)) {
    stop("Todos los valores de phi deben ser positivos.")
  }
  
  kl_draws <- numeric(ndraws)
  
  for (i in seq_len(ndraws)) {
    a1 <- mu_control[i, ] * phi_control[i]
    a2 <- mu_tratamiento[i, ] * phi_tratamiento[i]
    
    kl12 <- kl_dirichlet_alpha(a1, a2)
    
    if (symmetric) {
      kl21 <- kl_dirichlet_alpha(a2, a1)
      kl_draws[i] <- 0.5 * (kl12 + kl21)
    } else {
      kl_draws[i] <- kl12
    }
  }
  
  list(
    kl_draws = kl_draws,
    mean = mean(kl_draws),
    median = median(kl_draws),
    sd = sd(kl_draws),
    ci = quantile(kl_draws, probs = probs),
    symmetric = symmetric,
    probs = probs
  )
}


# ---------------------------------------------------------------
# Ejemplo mínimo de uso
# ---------------------------------------------------------------
# res <- divergencia_kl_dirichlet(
#   mu_control = mu_post_control,          # matriz draws x K
#   phi_control = phi_post_control,        # vector draws
#   mu_tratamiento = mu_post_tratamiento,  # matriz draws x K
#   phi_tratamiento = phi_post_tratamiento,# vector draws
#   symmetric = FALSE
# )
#
# res$mean
# res$median
# res$ci
# hist(res$kl_draws)
