####################################################################
# Calcular la divergencia de Jenssen-Shannon sobre coordenadas ilr #
####################################################################
# Esta acotada entre 0 y 1
# simetrica

divergencia_jsd_ilr <- function(z1, z2) {
  # z1 y z2 son matrices [16000, 19]
  
  # Función auxiliar para KL entre dos normales
  kl_norm <- function(m1, s1, m2, s2) {
    k <- length(m1)
    inv_s2 <- solve(s2)
    0.5 * (sum(diag(inv_s2 %*% s1)) + t(m2-m1) %*% inv_s2 %*% (m2-m1) - k + log(det(s2)/det(s1)))
  }
  
  # Parámetros de las dos distribuciones
  mu1 <- colMeans(z1); sig1 <- cov(z1)
  mu2 <- colMeans(z2); sig2 <- cov(z2)
  
  # Distribución mezcla (M = 0.5*P + 0.5*Q) aproximada
  mu_m <- 0.5 * (mu1 + mu2)
  sig_m <- 0.5 * (sig1 + sig2) + 0.25 * (mu1 - mu2) %*% t(mu1 - mu2)
  
  # JSD = 0.5 * KL(P||M) + 0.5 * KL(Q||M)
  jsd <- 0.5 * kl_norm(mu1, sig1, mu_m, sig_m) + 0.5 * kl_norm(mu2, sig2, mu_m, sig_m)
  
  return(as.numeric(jsd))
}
