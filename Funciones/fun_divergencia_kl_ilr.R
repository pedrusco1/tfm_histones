####################################################################
# Divergencia basada en la diferencia de las medias en espacio ILR #
####################################################################

divergencia_kl_ilr <- function(mu_control, mu_tratamiento) {
  # mu_control y mu_tratamiento deben ser vectores de medias (o matrices de draws)
  # Si pasas las matrices de posteriores (16000x19):
  m1 <- colMeans(mu_control)
  m2 <- colMeans(mu_tratamiento)
  s1 <- cov(mu_control)
  s2 <- cov(mu_tratamiento)
  
  inv_s2 <- solve(s2)
  k <- length(m1)
  
  # KL(Control || Tratamiento) para normales multivariantes
  kl <- 0.5 * (sum(diag(inv_s2 %*% s1)) + 
                 t(m2 - m1) %*% inv_s2 %*% (m2 - m1) - 
                 k + log(det(s2)/det(s1)))
  
  return(as.numeric(kl))
}
