###########################################
#    Cargar datos para modelos Dirichlet  #
###########################################

source(here::here("ScriptsR", "40_stan_config.R"))

# Cargar data frame procesado
load(here::here("DatosProcesados", "df_h3917.Rdata"))

#--------------------------------------------------
# 1. Definir data frame base
#--------------------------------------------------

datos_dirich <- df.ancho

#--------------------------------------------------
# 2. Definir columnas de la respuesta composicional
#--------------------------------------------------

response_cols <- colnames(datos_dirich)[3:22]

#--------------------------------------------------
# 3. Limpiar / recodificar predictor categórico
#--------------------------------------------------

datos_dirich$Estado <- factor(
  datos_dirich$Estado,
  levels = c("asinc", "mitot")
)

if (anyNA(datos_dirich$Estado)) {
  stop(
    "`Estado` contiene niveles distintos de 'asinc' y 'mitot' o valores perdidos. ",
    "Revisa la recodificación antes de ajustar los modelos."
  )
}

#--------------------------------------------------
# 4. Identificar variable de día, si existe
#--------------------------------------------------
# Los modelos dirichreg_phi_*_dia.stan necesitan un índice de día.
# Ajusta `dia_col` si tu columna tiene otro nombre.

posibles_cols_dia <- c("Dia", "dia", "Día", "day", "Day")
dia_col <- intersect(posibles_cols_dia, colnames(datos_dirich))

if (length(dia_col) > 1) {
  stop(
    "Se encontraron varias columnas candidatas para día: ",
    paste(dia_col, collapse = ", "),
    ". Deja solo una o define explícitamente `dia_col`."
  )
}

if (length(dia_col) == 1) {
  dia_col <- dia_col[1]
  datos_dirich$dia_var <- factor(datos_dirich[[dia_col]])

  if (anyNA(datos_dirich$dia_var)) {
    stop("La variable de día `", dia_col, "` contiene NA.")
  }

  dia_vector <- datos_dirich$dia_var
} else {
  dia_col <- NULL
  dia_vector <- NULL
  warning(
    "No se encontró columna de día. ",
    "Se podrán ajustar los modelos sin día, pero no los modelos *_dia."
  )
}

#--------------------------------------------------
# 5. Comprobar que la respuesta es válida
#--------------------------------------------------

Y_mat <- as.matrix(datos_dirich[, response_cols, drop = FALSE])

# Comprobar no negatividad
if (any(Y_mat < 0, na.rm = TRUE)) {
  stop("Hay valores negativos en las columnas de respuesta.")
}

# Comprobar filas con suma positiva
if (any(rowSums(Y_mat, na.rm = TRUE) <= 0)) {
  stop("Hay filas de la respuesta con suma <= 0.")
}

# Comprobar valores perdidos
if (anyNA(Y_mat)) {
  stop("Hay NA en las columnas de respuesta.")
}

#--------------------------------------------------
# 6. Guardar inputs
#--------------------------------------------------

save(
  datos_dirich,
  response_cols,
  dia_col,
  dia_vector,
  file = file.path(data_dir, "data_inputs_dirichlet.Rdata")
)
