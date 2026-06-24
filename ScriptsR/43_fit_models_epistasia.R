###############################################
#     Ajustar los modelos campo epistasia     #
#       4 modelos epistasia_kappa revisados   #
###############################################

source("ScriptsR/40_stan_config.R")
source("Funciones/fun_prepara_epistasia_data.R")
source("Funciones/fun_fit_stan_model.R")
source("Funciones/fun_stan_diagnosticos.R")
source("ScriptsR/41_cargar_datos_h3917.R")

load(file.path(data_dir, "data_inputs_h3917.Rdata"))
load(file.path(data_dir, "delta_obj.Rdata"))

# ------------------------------------------------------------
# Datos Stan sin efecto día: para epistasia_kappa_comun.stan
# y epistasia_kappa_cond.stan
# ------------------------------------------------------------
stan_data <- prepara_stan_data(
  Y_mat = Y_mat,
  cond = cond,
  Delta = delta_obj$Delta,
  rank_Delta = delta_obj$rank_Delta,
  prior_only = 0
)

# ------------------------------------------------------------
# Datos Stan con efecto día: para modelos *_dia.stan
#
# La función prepara_stan_data() añade D y day cuando se pasa dia.
# Se intenta usar primero un objeto llamado dia; si no existe,
# un objeto llamado Dia; y, finalmente, df.ancho$Dia.
# ------------------------------------------------------------
if (exists("dia")) {
  dia_vec <- dia
} else if (exists("Dia")) {
  dia_vec <- Dia
} else if (exists("df.ancho") && "Dia" %in% names(df.ancho)) {
  dia_vec <- df.ancho$Dia
} else {
  stop(
    "No encuentro el vector de día experimental. ",
    "Debe existir un objeto `dia`, `Dia` o `df.ancho$Dia` ",
    "en data_inputs_h3917.Rdata o en ScriptsR/41_cargar_datos_h3917.R."
  )
}

stan_data_dia <- prepara_stan_data(
  Y_mat = Y_mat,
  cond = cond,
  Delta = delta_obj$Delta,
  rank_Delta = delta_obj$rank_Delta,
  prior_only = 0,
  dia = dia_vec
)

# ------------------------------------------------------------
# Ajuste de los cuatro modelos epistasia_kappa actuales
# ------------------------------------------------------------
fit_kappa_comun <- fit_stan_model(
  file.path(stan_dir, "epistasia_kappa_comun.stan"),
  stan_data
)

fit_kappa_comun_dia <- fit_stan_model(
  file.path(stan_dir, "epistasia_kappa_comun_dia.stan"),
  stan_data_dia
)

fit_kappa_cond <- fit_stan_model(
  file.path(stan_dir, "epistasia_kappa_cond.stan"),
  stan_data
)

fit_kappa_cond_dia <- fit_stan_model(
  file.path(stan_dir, "epistasia_kappa_cond_dia.stan"),
  stan_data_dia
)

# ------------------------------------------------------------
# Diagnósticos
# ------------------------------------------------------------
diag_comun <- compute_diagnostics(fit_kappa_comun)
diag_comun_dia <- compute_diagnostics(fit_kappa_comun_dia)
diag_cond <- compute_diagnostics(fit_kappa_cond)
diag_cond_dia <- compute_diagnostics(fit_kappa_cond_dia)

# ------------------------------------------------------------
# Guardar resultados
# ------------------------------------------------------------
save(
  fit_kappa_comun,
  fit_kappa_comun_dia,
  fit_kappa_cond,
  fit_kappa_cond_dia,
  diag_comun,
  diag_comun_dia,
  diag_cond,
  diag_cond_dia,
  file = file.path(data_dir, "stan_epistasia_fits.Rdata")
)
