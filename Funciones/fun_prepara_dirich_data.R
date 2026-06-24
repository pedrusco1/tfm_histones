######################################################
#   Preparar datos para regresión Dirichlet          #
#   Con soporte opcional para condición y día        #
######################################################

#' Preparar datos para modelos Stan de regresión Dirichlet
#'
#' Construye la lista de datos que requieren los modelos Stan de regresión
#' Dirichlet. La función prepara la matriz de diseño `X`, la matriz de
#' composiciones `Y`, y opcionalmente los índices de condición (`J`,
#' `condition`) y de día (`D`, `day`).
#'
#' Los argumentos `condicion` y `dia` son opcionales porque no todos los
#' modelos Stan usan los mismos bloques de datos:
#'
#' - modelos comunes: no necesitan `J` ni `condition`;
#' - modelos condicionales: necesitan `J` y `condition`;
#' - modelos con efecto aleatorio de día: necesitan `D` y `day`.
#'
#' @param data `data.frame` con las variables predictoras, respuestas y,
#'   opcionalmente, condición y día.
#' @param formula fórmula de R usada para construir la matriz de diseño `X`
#'   mediante `model.matrix()`.
#' @param response_cols vector de nombres o índices de las columnas que forman
#'   la composición respuesta. Cada fila debe contener valores no negativos y
#'   con suma positiva.
#' @param condicion vector opcional con la condición experimental de cada fila.
#'   Si se proporciona, se convierte a índices enteros 1, ..., J y se añaden
#'   `J` y `condition` al objeto devuelto.
#' @param dia vector opcional con el día experimental de cada fila, por ejemplo
#'   `df.ancho$Dia`. Si se proporciona, se convierte a índices enteros
#'   1, ..., D y se añaden `D` y `day` al objeto devuelto.
#' @param prior_only entero 0/1. Si es 1, el modelo Stan debe simular del prior
#'   sin evaluar la likelihood; si es 0, ajusta usando los datos observados.
#'
#' @return Lista preparada para pasar al argumento `data` de `rstan::stan()`.
#'
#' @examples
#' stan_data <- prepara_dirich_data(
#'   data = df.ancho,
#'   formula = ~ Condicion,
#'   response_cols = cols_ptm,
#'   condicion = df.ancho$Condicion,
#'   dia = df.ancho$Dia,
#'   prior_only = 0
#' )
prepara_dirich_data <- function(data,
                                formula,
                                response_cols,
                                condicion = NULL,
                                dia = NULL,
                                prior_only = 0) {

  # -----------------------------
  # Validaciones básicas
  # -----------------------------
  if (!is.data.frame(data)) {
    stop("`data` debe ser un data.frame.")
  }

  if (!prior_only %in% c(0, 1)) {
    stop("`prior_only` debe ser 0 o 1.")
  }

  # -----------------------------
  # Matriz de diseño y respuesta
  # -----------------------------
  X <- model.matrix(formula, data = data)
  Y <- as.matrix(data[, response_cols, drop = FALSE])

  if (nrow(X) != nrow(Y)) {
    stop("`X` e `Y` deben tener el mismo número de filas.")
  }

  if (anyNA(Y)) {
    stop("`Y` contiene NA. Stan no acepta valores perdidos.")
  }

  if (any(Y < 0)) {
    stop("Las composiciones en `Y` no pueden contener valores negativos.")
  }

  row_sums <- rowSums(Y)
  if (any(row_sums <= 0)) {
    stop("Todas las filas de `Y` deben tener suma positiva.")
  }

  # Normalización para evitar pequeños errores de redondeo.
  Y <- sweep(Y, 1, row_sums, "/")

  # Los modelos Dirichlet requieren componentes estrictamente positivos.
  if (any(Y <= 0)) {
    stop(
      "La Dirichlet requiere composiciones estrictamente positivas. ",
      "Hay ceros en `Y`; considera aplicar un ajuste/pseudoconteo antes."
    )
  }

  # -----------------------------
  # Lista mínima común
  # -----------------------------
  stan_data <- list(
    N = nrow(X),
    K = ncol(Y),
    P = ncol(X),
    Y = Y,
    X = X,
    prior_only = as.integer(prior_only)
  )

  # -----------------------------
  # Condición experimental opcional
  # -----------------------------
  if (!is.null(condicion)) {
    if (length(condicion) != nrow(Y)) {
      stop("`condicion` debe tener longitud N.")
    }
    if (anyNA(condicion)) {
      stop("`condicion` contiene NA.")
    }

    cond_vec <- as.integer(as.factor(condicion))
    stan_data$J <- max(cond_vec)
    stan_data$condition <- cond_vec
  }

  # -----------------------------
  # Día experimental opcional
  # -----------------------------
  if (!is.null(dia)) {
    if (length(dia) != nrow(Y)) {
      stop("`dia` debe tener longitud N.")
    }
    if (anyNA(dia)) {
      stop("`dia` contiene NA.")
    }

    day_vec <- as.integer(as.factor(dia))
    stan_data$D <- max(day_vec)
    stan_data$day <- day_vec
  }

  return(stan_data)
}
