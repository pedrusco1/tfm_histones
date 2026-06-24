library(rstan)
library(bayesplot)
library(posterior)

diagnostico_funnel_ligero <- function(fit,
                                      pars,
                                      n_draws = 1000,
                                      file = NULL,
                                      width = 1200,
                                      height = 1000,
                                      res = 150,
                                      seed = 123) {
  
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("Package 'bayesplot' required")
  }
  
  # Extrae solo los parámetros pedidos y aplana cadenas
  draws_mat <- as.matrix(fit, pars = pars)
  draws_df  <- as.data.frame(draws_mat)
  
  # Comprobar nombres exactos
  missing <- setdiff(pars, colnames(draws_df))
  if (length(missing) > 0) {
    stop(
      "Estos parámetros no existen exactamente: ",
      paste(missing, collapse = ", "),
      "\nParámetros disponibles parecidos:\n",
      paste(
        grep(paste(gsub("\\[|\\]", "\\\\[|\\\\]", missing), collapse = "|"),
             colnames(draws_df),
             value = TRUE),
        collapse = ", "
      )
    )
  }
  
  set.seed(seed)
  idx <- sample(seq_len(nrow(draws_df)), size = min(n_draws, nrow(draws_df)))
  draws_small <- draws_df[idx, pars, drop = FALSE]
  
  if (!is.null(file)) {
    png(file, width = width, height = height, res = res)
    on.exit({
      dev.off()
      gc()
    }, add = TRUE)
  }
  
  p <- bayesplot::mcmc_pairs(
    draws_small,
    pars = pars
  )
  
  print(p)
  invisible(p)
}
pp <- diagnostico_funnel_ligero(
  fit_lognor_v1_sum0coef_delta_comun,
  pars = c("tau", "B[11,1]", "B[19,2]"),
  n_draws = 800
)
