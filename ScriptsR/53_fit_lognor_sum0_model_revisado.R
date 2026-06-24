###########################################################################
# AJUSTE DE 4 MODELOS LOG-NORMAL V1 REVISADOS (SUMA CERO EN CLR)
# H3917
#
# Genera objetos nuevos, sin sobrescribir los ajustes log-normales antiguos.
###########################################################################

library(rstan)
library(here)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# -------------------------------------------------------------------------
# 1. Configuración y funciones auxiliares
# -------------------------------------------------------------------------
source(here::here("ScriptsR", "40_stan_config.R"))
source(here::here("Funciones", "fun_prepara_data_lognor_v1.R"))
source(here::here("ScriptsR", "41_cargar_datos_h3917.R"))
source(here::here("Funciones", "fun_stan_diagnosticos.R"))

# -------------------------------------------------------------------------
# 2. Cargar datos
# -------------------------------------------------------------------------
load(file.path(data_dir, "delta_obj.Rdata"))
load(file.path(data_dir, "data_inputs_h3917.Rdata"))

objetos_necesarios <- c("Y_mat", "df.ancho", "delta_obj")
faltan <- objetos_necesarios[!vapply(objetos_necesarios, exists, logical(1), inherits = TRUE)]
if (length(faltan) > 0) stop("Faltan objetos necesarios: ", paste(faltan, collapse = ", "))
if (!"Delta" %in% names(delta_obj)) stop("No se encuentra delta_obj$Delta")
if (!all(c("Estado", "Dia") %in% names(df.ancho))) stop("df.ancho debe contener Estado y Dia")
if (nrow(df.ancho) != nrow(Y_mat)) stop("df.ancho e Y_mat deben tener el mismo número de filas")
if (!"day" %in% names(formals(make_stan_data))) {
  stop("make_stan_data() debe aceptar el argumento day para los modelos revisados con día")
}

# -------------------------------------------------------------------------
# 3. Archivos Stan revisados
# -------------------------------------------------------------------------
stan_files_revisados <- c(
  lognor_v1_sum0coef_revisado = file.path(
    stan_dir, "epistasia_log_normal_v1_sum0coef_revisado.stan"
  ),
  lognor_v1_dia_sum0coef_revisado = file.path(
    stan_dir, "epistasia_log_normal_v1_dia_sum0coef_revisado.stan"
  ),
  lognor_v1_sum0coef_delta_comun_revisado = file.path(
    stan_dir, "epistasia_log_normal_v1_sum0coef_delta_comun_revisado.stan"
  ),
  lognor_v1_dia_sum0coef_delta_comun_revisado = file.path(
    stan_dir, "epistasia_log_normal_v1_dia_sum0coef_delta_comun_revisado.stan"
  )
)

missing_stan_files <- stan_files_revisados[!file.exists(stan_files_revisados)]
if (length(missing_stan_files) > 0) {
  stop(
    "Faltan los siguientes archivos Stan revisados:\n",
    paste(names(missing_stan_files), missing_stan_files, sep = ": ", collapse = "\n")
  )
}

# -------------------------------------------------------------------------
# 4. Matriz de diseño y datos Stan
# -------------------------------------------------------------------------
X_mat_revisado <- model.matrix(~ Estado, data = df.ancho)

stan_data_lognor_v1_sin_dia_revisado <- make_stan_data(
  Y = Y_mat,
  L = delta_obj$Delta,
  X = X_mat_revisado
)

stan_data_lognor_v1_con_dia_revisado <- make_stan_data(
  Y = Y_mat,
  L = delta_obj$Delta,
  X = X_mat_revisado,
  day = df.ancho$Dia
)

message(
  "Niveles de día usados en Stan: ",
  paste(attr(stan_data_lognor_v1_con_dia_revisado, "day_levels"), collapse = ", ")
)

# -------------------------------------------------------------------------
# 5. Configuración de muestreo
# -------------------------------------------------------------------------
iter <- 8000
warmup <- 4000
chains <- 4
cores <- min(chains, parallel::detectCores())
adapt_delta <- 0.999
max_treedepth <- 15

model_specs_revisados <- list(
  fit_lognor_v1_sum0coef_revisado = list(
    stan_file = stan_files_revisados[["lognor_v1_sum0coef_revisado"]],
    data = stan_data_lognor_v1_sin_dia_revisado,
    descripcion = "log-normal v1 revisado, delta jerárquico, sin día"
  ),
  fit_lognor_v1_dia_sum0coef_revisado = list(
    stan_file = stan_files_revisados[["lognor_v1_dia_sum0coef_revisado"]],
    data = stan_data_lognor_v1_con_dia_revisado,
    descripcion = "log-normal v1 revisado, delta jerárquico, con día"
  ),
  fit_lognor_v1_sum0coef_delta_comun_revisado = list(
    stan_file = stan_files_revisados[["lognor_v1_sum0coef_delta_comun_revisado"]],
    data = stan_data_lognor_v1_sin_dia_revisado,
    descripcion = "log-normal v1 revisado, delta común, sin día"
  ),
  fit_lognor_v1_dia_sum0coef_delta_comun_revisado = list(
    stan_file = stan_files_revisados[["lognor_v1_dia_sum0coef_delta_comun_revisado"]],
    data = stan_data_lognor_v1_con_dia_revisado,
    descripcion = "log-normal v1 revisado, delta común, con día"
  )
)

ajusta_modelo_lognor <- function(spec, nombre) {
  message("\n============================================================")
  message("Ajustando ", nombre, ": ", spec$descripcion)
  message("Archivo Stan: ", spec$stan_file)
  message("============================================================\n")

  rstan::stan(
    file = spec$stan_file,
    data = spec$data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )
}

check_B_sum0 <- function(fit, tol = 1e-6) {
  s <- rstan::summary(fit, pars = "B")$summary
  rn <- rownames(s)
  indices <- regmatches(rn, gregexpr("[0-9]+", rn))
  k_idx <- as.integer(vapply(indices, `[`, character(1), 1))
  p_idx <- as.integer(vapply(indices, `[`, character(1), 2))
  B_mean <- matrix(NA_real_, nrow = max(k_idx), ncol = max(p_idx))
  for (i in seq_along(rn)) B_mean[k_idx[i], p_idx[i]] <- s[i, "mean"]
  out <- colSums(B_mean)
  if (any(abs(out) > tol)) warning("Alguna columna de B no suma aproximadamente cero")
  out
}

# -------------------------------------------------------------------------
# 6. Ajustar y crear objetos nuevos
# -------------------------------------------------------------------------
fits_lognor_v1_4modelos_revisados <- list()
inicio <- Sys.time()

for (nm in names(model_specs_revisados)) {
  fits_lognor_v1_4modelos_revisados[[nm]] <- ajusta_modelo_lognor(
    model_specs_revisados[[nm]], nm
  )
  gc()
}

list2env(fits_lognor_v1_4modelos_revisados, envir = .GlobalEnv)
message("Tiempo total: ", round(difftime(Sys.time(), inicio, units = "mins"), 2), " minutos")

# -------------------------------------------------------------------------
# 7. Diagnósticos y comprobación de suma cero
# -------------------------------------------------------------------------
diagnosticos_lognor_v1_4modelos_revisados <- lapply(
  fits_lognor_v1_4modelos_revisados,
  compute_diagnostics
)

names_diag <- sub("^fit_", "diag_", names(diagnosticos_lognor_v1_4modelos_revisados))
for (i in seq_along(diagnosticos_lognor_v1_4modelos_revisados)) {
  assign(
    names_diag[i],
    diagnosticos_lognor_v1_4modelos_revisados[[i]],
    envir = .GlobalEnv
  )
}

B_sum0_checks_revisados <- lapply(
  fits_lognor_v1_4modelos_revisados,
  check_B_sum0
)

# -------------------------------------------------------------------------
# 8. Registro y guardado separado de los modelos antiguos
# -------------------------------------------------------------------------
model_registry_revisado <- data.frame(
  objeto = names(model_specs_revisados),
  descripcion = vapply(model_specs_revisados, `[[`, character(1), "descripcion"),
  stan_file = vapply(model_specs_revisados, `[[`, character(1), "stan_file"),
  usa_dia = grepl("_dia_", names(model_specs_revisados)),
  delta_comun = grepl("delta_comun", names(model_specs_revisados)),
  delta_jerarquico = !grepl("delta_comun", names(model_specs_revisados)),
  revisado = TRUE,
  stringsAsFactors = FALSE
)

out_file_revisado <- file.path(data_dir, "stan_lognormal_sum0_fits_revisados.Rdata")

save(
  fits_lognor_v1_4modelos_revisados,
  diagnosticos_lognor_v1_4modelos_revisados,
  B_sum0_checks_revisados,
  X_mat_revisado,
  stan_files_revisados,
  model_registry_revisado,
  stan_data_lognor_v1_sin_dia_revisado,
  stan_data_lognor_v1_con_dia_revisado,

  fit_lognor_v1_sum0coef_revisado,
  fit_lognor_v1_dia_sum0coef_revisado,
  fit_lognor_v1_sum0coef_delta_comun_revisado,
  fit_lognor_v1_dia_sum0coef_delta_comun_revisado,

  diag_lognor_v1_sum0coef_revisado,
  diag_lognor_v1_dia_sum0coef_revisado,
  diag_lognor_v1_sum0coef_delta_comun_revisado,
  diag_lognor_v1_dia_sum0coef_delta_comun_revisado,
  file = out_file_revisado
)

message("Modelos revisados guardados en: ", out_file_revisado)
