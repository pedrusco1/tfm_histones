###################################################
#         PPC TOTAL VARANCE LOG-RATIO             #
###################################################

# T(Y) = suma de varianzas por columna en CLR, equivalente a tr(cov(CLR(Y)))
T_clr <- function(Y, eps = 1e-8) {
  Y <- pmax(Y, eps)
  Z <- log(Y)
  Z <- Z - rowMeans(Z)                     # CLR por fila (composición)
  Zc <- scale(Z, center = TRUE, scale = FALSE)  # centrar entre réplicas
  sum(Zc^2) / (nrow(Z) - 1)                # = sum(Var(col)) con ddof=1
}

# Calcula T para cada draw posterior en un array S x n x k
posterior_T_clr <- function(post, eps = 1e-8, Y_obs = NULL) {
  d <- dim(post)
  stopifnot(length(d) == 3)
  S <- d[1]
  
  T_draws <- vapply(
    seq_len(S),
    function(s) T_clr(post[s, , , drop = FALSE][1, , ], eps = eps),
    numeric(1)
  )
  
  out <- list(T_draws = T_draws)
  
  if (!is.null(Y_obs)) {
    T_obs <- T_clr(Y_obs, eps = eps)
    out$T_obs <- T_obs
    out$p_bayes <- mean(T_draws >= T_obs)
  }
  
  out
}