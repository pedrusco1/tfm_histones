######################################################
#   Ajustar modelos de regresión Dirichlet           #
#   H3917: modelos sin/con phi por condición y día   #
######################################################

source("ScriptsR/40_stan_config.R")
source("Funciones/fun_prepara_dirich_data.R")
source("Funciones/fun_fit_stan_model.R")
source("Funciones/fun_stan_diagnosticos.R")
source("ScriptsR/31_cargar_datos_dirichreg_h3917.R")

#-----------------------------------------------------
# 1. Cargar datos
#-----------------------------------------------------
load(file.path(data_dir, "data_inputs_dirichlet.Rdata"))

# Variables auxiliares comunes
# - condicion_var entra en la matriz X como predictor de la media composicional.
# - condicion_vector se pasa aparte a Stan cuando phi depende de la condición.
datos_dirich$condicion_var <- datos_dirich$Estado
condicion_vector <- datos_dirich$Estado

# `dia_vector` queda definido en 31_cargar_datos_dirichreg_h3917.R.
# Si no existe columna de día, los modelos *_dia no se ajustan.
ajustar_modelos_dia <- !is.null(dia_vector)

#-----------------------------------------------------
# 2. Preparar datos Stan
#-----------------------------------------------------

# Modelos sin phi por condición y sin día
stan_data_phi_comun <- prepara_dirich_data(
  data = datos_dirich,
  formula = ~ condicion_var,
  response_cols = response_cols,
  prior_only = 0
)

# Modelos con phi por condición y sin día
stan_data_phi_cond <- prepara_dirich_data(
  data = datos_dirich,
  formula = ~ condicion_var,
  response_cols = response_cols,
  condicion = condicion_vector,
  prior_only = 0
)

# Modelos con efecto aleatorio de día
if (ajustar_modelos_dia) {
  stan_data_phi_comun_dia <- prepara_dirich_data(
    data = datos_dirich,
    formula = ~ condicion_var,
    response_cols = response_cols,
    dia = dia_vector,
    prior_only = 0
  )

  stan_data_phi_cond_dia <- prepara_dirich_data(
    data = datos_dirich,
    formula = ~ condicion_var,
    response_cols = response_cols,
    condicion = condicion_vector,
    dia = dia_vector,
    prior_only = 0
  )
}

#-----------------------------------------------------
# 3. Definir catálogo de modelos a ajustar
#-----------------------------------------------------

modelos_dirich <- list(
  phi_comun = list(
    stan_file = file.path(stan_dir, "dirichreg_phi_comun.stan"),
    data_list = stan_data_phi_comun
  ),
  phi_cond = list(
    stan_file = file.path(stan_dir, "dirichreg_phi_cond.stan"),
    data_list = stan_data_phi_cond
  )
)

if (ajustar_modelos_dia) {
  modelos_dirich$phi_comun_dia <- list(
    stan_file = file.path(stan_dir, "dirichreg_phi_comun_dia.stan"),
    data_list = stan_data_phi_comun_dia
  )

  modelos_dirich$phi_cond_dia <- list(
    stan_file = file.path(stan_dir, "dirichreg_phi_cond_dia.stan"),
    data_list = stan_data_phi_cond_dia
    #adapt_delta = 0.99,
    #max_treedepth = 18
  )
}

# Comprobar que los ficheros Stan existen antes de lanzar ajustes largos.
stan_files_exist <- vapply(
  modelos_dirich,
  function(x) file.exists(x$stan_file),
  logical(1)
)

if (any(!stan_files_exist)) {
  stop(
    "Faltan ficheros Stan: ",
    paste(vapply(modelos_dirich[!stan_files_exist], `[[`, character(1), "stan_file"),
          collapse = ", ")
  )
}

#-----------------------------------------------------
# 4. Ajustar modelos Dirichlet
#-----------------------------------------------------

fits_dirich <- lapply(modelos_dirich, function(m) {
  fit_stan_model(
    stan_file = m$stan_file,
    data_list = m$data_list
  )
})


# Mantener objetos individuales por compatibilidad con scripts previos
fit_dirich_phi_comun <- fits_dirich$phi_comun
fit_dirich_phi_cond  <- fits_dirich$phi_cond

if (ajustar_modelos_dia) {
  fit_dirich_phi_comun_dia <- fits_dirich$phi_comun_dia
  fit_dirich_phi_cond_dia  <- fits_dirich$phi_cond_dia
}

#-----------------------------------------------------
# 5. Diagnósticos
#-----------------------------------------------------

diag_dirich <- lapply(fits_dirich, compute_diagnostics)

diag_dirich_phi_comun <- diag_dirich$phi_comun
diag_dirich_phi_cond  <- diag_dirich$phi_cond

if (ajustar_modelos_dia) {
  diag_dirich_phi_comun_dia <- diag_dirich$phi_comun_dia
  diag_dirich_phi_cond_dia  <- diag_dirich$phi_cond_dia
}

#-----------------------------------------------------
# 6. Guardar resultados de ajuste
#-----------------------------------------------------

save(
  fits_dirich,
  diag_dirich,
  modelos_dirich,
  stan_data_phi_comun,
  stan_data_phi_cond,
  fit_dirich_phi_comun,
  fit_dirich_phi_cond,
  diag_dirich_phi_comun,
  diag_dirich_phi_cond,
  list = if (ajustar_modelos_dia) {
    c(
      "stan_data_phi_comun_dia",
      "stan_data_phi_cond_dia",
      "fit_dirich_phi_comun_dia",
      "fit_dirich_phi_cond_dia",
      "diag_dirich_phi_comun_dia",
      "diag_dirich_phi_cond_dia"
    )
  } else {
    character(0)
  },
  file = file.path(data_dir, "stan_dirichlet_fits.Rdata")
)





