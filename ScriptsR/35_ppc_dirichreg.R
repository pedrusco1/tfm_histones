############################################################
# PPC PARA LOS 4 MODELOS DIRICHLET
# Compatible con objetos stanfit de rstan
# Todas las rutas usan here::here()
#
# Requiere:
#   - DatosProcesados/data_inputs_dirichlet.Rdata
#   - DatosProcesados/stan_dirichlet_fits.Rdata
#   - Funciones/fun_ppc_globales.R
#   - Funciones/fun_ppc_margin_centro_y_var.R
#   - Funciones/fun_ppc_balances.R
#   - Funciones/fun_ppc_varianza.R
#
# Salidas en DatosProcesados/PPC:
#   - Figuras PNG por modelo
#   - RDS completo por modelo
#   - CSV resumen completo por modelo
#   - CSV/RDS resumen conjunto para todos los modelos
############################################################

# ==========================================================
# 1. LIBRERÍAS Y FUNCIONES
# ==========================================================

library(rstan)
library(here)
library(ggplot2)

source(here("Funciones", "fun_ppc_globales.R"))
source(here("Funciones", "fun_ppc_margin_centro_y_var.R"))
source(here("Funciones", "fun_ppc_balances.R"))
source(here("Funciones", "fun_ppc_varianza.R"))

# ==========================================================
# 2. CARGAR DATOS Y AJUSTES STAN
# ==========================================================

load(here("DatosProcesados", "data_inputs_dirichlet.Rdata"))
load(here("DatosProcesados", "stan_dirichlet_fits.Rdata"))

# ==========================================================
# 3. RUTAS Y DATOS OBSERVADOS
# ==========================================================

out_dir <- here("DatosProcesados", "PPC")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

Y_obs <- as.matrix(datos_dirich[, response_cols])

# ==========================================================
# 4. FUNCIONES AUXILIARES
# ==========================================================

ruta_ppc <- function(...) {
  file.path(out_dir, ...)
}

guardar_grid_png <- function(plot,
                             filename,
                             width = 1800,
                             height = 900,
                             res = 150) {

  png(
    filename = ruta_ppc(filename),
    width = width,
    height = height,
    res = res
  )

  grid::grid.draw(plot)
  dev.off()
}

guardar_ggplot_png <- function(plot,
                               filename,
                               width = 10,
                               height = 7) {

  ggsave(
    filename = ruta_ppc(filename),
    plot = plot,
    width = width,
    height = height
  )
}

extraer_y_rep <- function(fit, model_name) {

  post <- rstan::extract(fit)
  y_rep <- post$Y_rep

  if(is.null(y_rep)) {
    stop(paste("No existe Y_rep en:", model_name))
  }

  return(y_rep)
}

crear_resumen_ppc <- function(model_name,
                              ppc_global,
                              ppc_margin_centro,
                              ppc_margin_var_clr,
                              ppc_bal,
                              ppc_var_global) {

  # --------------------------------------------------------
  # PPC globales
  # --------------------------------------------------------

  resumen_global <- data.frame(
    modelo = model_name,
    bloque_ppc = "global",
    estadistico = c(
      "centroide_aitchison",
      "variacion_total_clr"
    ),
    componente = "global",
    p_ppc = c(
      ppc_global$p_valor_centroide,
      ppc_global$p_valor_variacion
    ),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------
  # PPC marginal del centro composicional por proteoforma
  # --------------------------------------------------------

  resumen_marginal_centro <- data.frame(
    modelo = model_name,
    bloque_ppc = "marginal",
    estadistico = "centroide_por_proteoforma",
    componente = ppc_margin_centro$stats$Marca,
    p_ppc = ppc_margin_centro$stats$p_PPC,
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------
  # PPC marginal de la varianza CLR por proteoforma
  # IMPORTANTE: este bloque queda incluido en el rbind final.
  # --------------------------------------------------------

  resumen_marginal_var_clr <- data.frame(
    modelo = model_name,
    bloque_ppc = "marginal",
    estadistico = "varianza_clr_por_proteoforma",
    componente = ppc_margin_var_clr$stats$Marca,
    p_ppc = ppc_margin_var_clr$stats$p_PPC,
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------
  # PPC de balances: media y varianza
  # --------------------------------------------------------

  resumen_bal_medias <- data.frame(
    modelo = model_name,
    bloque_ppc = "balances",
    estadistico = "media_balance",
    componente = names(ppc_bal$p_pcc_medias),
    p_ppc = as.numeric(unlist(ppc_bal$p_pcc_medias)),
    stringsAsFactors = FALSE
  )

  resumen_bal_varianzas <- data.frame(
    modelo = model_name,
    bloque_ppc = "balances",
    estadistico = "varianza_balance",
    componente = names(ppc_bal$p_pcc_varianzas),
    p_ppc = as.numeric(unlist(ppc_bal$p_pcc_varianzas)),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------
  # PPC global adicional: T_clr
  # --------------------------------------------------------

  resumen_var_clr_global <- data.frame(
    modelo = model_name,
    bloque_ppc = "varianza_clr",
    estadistico = "T_clr",
    componente = "global",
    p_ppc = ppc_var_global$p_bayes,
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------
  # Resumen completo
  # --------------------------------------------------------

  resumen <- rbind(
    resumen_global,
    resumen_marginal_centro,
    resumen_marginal_var_clr,
    resumen_bal_medias,
    resumen_bal_varianzas,
    resumen_var_clr_global
  )

  return(resumen)
}

# ==========================================================
# 5. LOOP SOBRE LOS 4 MODELOS
# ==========================================================

resumen_todos_modelos <- list()

for(model_name in names(fits_dirich)) {

  cat("\n=====================================\n")
  cat("Procesando:", model_name, "\n")
  cat("=====================================\n")

  fit <- fits_dirich[[model_name]]
  y_rep <- extraer_y_rep(fit, model_name)

  # --------------------------------------------------------
  # 1. PPC globales
  # --------------------------------------------------------

  ppc_global <- check_compo_ppc(
    y_obs = Y_obs,
    y_rep = y_rep
  )

  guardar_grid_png(
    plot = ppc_global$plot,
    filename = paste0(model_name, "_ppc_global.png"),
    width = 1800,
    height = 900
  )

  # --------------------------------------------------------
  # 2. PPC marginal: centro composicional por proteoforma
  # --------------------------------------------------------

  ppc_margin_centro <- check_ppc_margin(
    y_obs = Y_obs,
    y_rep = y_rep,
    title = paste("PPC Marginal - Centroide -", model_name)
  )

  guardar_ggplot_png(
    plot = ppc_margin_centro$grafico,
    filename = paste0(model_name, "_ppc_marginal_centroide.png"),
    width = 10,
    height = 7
  )

  # --------------------------------------------------------
  # 3. PPC marginal: varianza CLR por proteoforma
  # --------------------------------------------------------

  ppc_margin_var_clr <- check_ppc_margin_var_clr(
    y_obs = Y_obs,
    y_rep = y_rep,
    title = paste("PPC Marginal - Varianza CLR -", model_name)
  )

  guardar_ggplot_png(
    plot = ppc_margin_var_clr$grafico,
    filename = paste0(model_name, "_ppc_marginal_varianza_clr.png"),
    width = 10,
    height = 7
  )

  # --------------------------------------------------------
  # 4. PPC balances biológicos
  # --------------------------------------------------------

  ppc_bal <- check_balances_ppc(
    y_obs = Y_obs,
    y_rep = y_rep
  )

  guardar_grid_png(
    plot = ppc_bal$plot_localizacion,
    filename = paste0(model_name, "_balances_media.png"),
    width = 1800,
    height = 700
  )

  guardar_grid_png(
    plot = ppc_bal$plot_dispersion,
    filename = paste0(model_name, "_balances_varianza.png"),
    width = 1800,
    height = 700
  )

  # --------------------------------------------------------
  # 5. PPC global: variación total CLR
  # --------------------------------------------------------

  ppc_var_global <- posterior_T_clr(
    post = y_rep,
    Y_obs = Y_obs
  )

  # --------------------------------------------------------
  # 6. Guardar resultados completos del modelo
  # --------------------------------------------------------

  saveRDS(
    list(
      ppc_global = ppc_global,
      ppc_margin_centro = ppc_margin_centro,
      ppc_margin_var_clr = ppc_margin_var_clr,
      ppc_balances = ppc_bal,
      ppc_varianza_global = ppc_var_global
    ),
    file = ruta_ppc(paste0(model_name, "_ppc_results.rds"))
  )

  # --------------------------------------------------------
  # 7. Tabla resumen completa del modelo
  # --------------------------------------------------------

  resumen_modelo <- crear_resumen_ppc(
    model_name = model_name,
    ppc_global = ppc_global,
    ppc_margin_centro = ppc_margin_centro,
    ppc_margin_var_clr = ppc_margin_var_clr,
    ppc_bal = ppc_bal,
    ppc_var_global = ppc_var_global
  )

  write.csv(
    resumen_modelo,
    file = ruta_ppc(paste0(model_name, "_resumen_ppc_completo.csv")),
    row.names = FALSE
  )

  resumen_todos_modelos[[model_name]] <- resumen_modelo

  cat("Finalizado:", model_name, "\n")
}

# ==========================================================
# 6. TABLA RESUMEN CONJUNTA PARA TODOS LOS MODELOS
# ==========================================================

resumen_ppc_todos_modelos <- do.call(
  rbind,
  resumen_todos_modelos
)

write.csv(
  resumen_ppc_todos_modelos,
  file = ruta_ppc("resumen_ppc_todos_modelos_completo.csv"),
  row.names = FALSE
)

saveRDS(
  resumen_ppc_todos_modelos,
  file = ruta_ppc("resumen_ppc_todos_modelos_completo.rds")
)

cat("\n=====================================\n")
cat("PPC COMPLETADO PARA TODOS LOS MODELOS\n")
cat("Resumen conjunto guardado en DatosProcesados/PPC\n")
cat("=====================================\n")
