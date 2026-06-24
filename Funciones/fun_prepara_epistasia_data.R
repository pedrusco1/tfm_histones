###################################################
#      Preparar datos para modelos de epistasia   #
#      Con soporte opcional para día experimental #
###################################################

#' Preparar datos para modelos Stan de epistasia
#'
#' Construye la lista de datos necesaria para los modelos Stan de epistasia.
#' Incluye la composición observada `Y`, el índice de condición `cond`, el
#' operador epistático `Delta`, su rango efectivo `rank_Delta`, la escala `s`
#' y el indicador `prior_only`.
#'
#' Si se proporciona `dia`, la función añade además `D` y `day`, necesarios
#' para las versiones del modelo con efecto aleatorio de día.
#'
#' @param Y_mat matriz N x K con composiciones observadas. Cada fila debe sumar
#'   1, contener valores no negativos y no contener ceros si se usa likelihood
#'   Dirichlet.
#' @param cond vector de condición experimental de longitud N. Se convierte a
#'   índices enteros 1, ..., J si no viene ya codificado.
#' @param Delta matriz K x K que define la estructura epistática.
#' @param rank_Delta rango efectivo de `Delta` usado en el prior epistático.
#' @param s escala positiva del prior epistático.
#' @param prior_only entero 0/1. Si es 1, el modelo Stan simula del prior; si
#'   es 0, evalúa la likelihood.
#' @param dia vector opcional con el día experimental de cada observación, por
#'   ejemplo `df.ancho$Dia`. Si se proporciona, se convierte a índices enteros
#'   1, ..., D y se añaden `D` y `day`.
#'
#' @return Lista preparada para pasar al argumento `data` de `rstan::stan()`.
#'
#' @examples
#' stan_data <- prepara_stan_data(
#'   Y_mat = Y,
#'   cond = df.ancho$Condicion,
#'   Delta = Delta,
#'   rank_Delta = rank_Delta,
#'   s = 45,
#'   prior_only = 0,
#'   dia = df.ancho$Dia
#' )
prepara_stan_data <- function(Y_mat,
                              cond,
                              Delta,
                              rank_Delta,
                              s = 45,
                              prior_only = 0,
                              dia = NULL) {

  # -----------------------------
  # Validaciones básicas
  # -----------------------------
  Y_mat <- as.matrix(Y_mat)

  if (!prior_only %in% c(0, 1)) {
    stop("`prior_only` debe ser 0 o 1.")
  }

  if (anyNA(Y_mat)) {
    stop("`Y_mat` contiene NA. Stan no acepta valores perdidos.")
  }

  if (any(Y_mat < 0)) {
    stop("`Y_mat` no puede contener valores negativos.")
  }

  row_sums <- rowSums(Y_mat)
  if (any(row_sums <= 0)) {
    stop("Todas las filas de `Y_mat` deben tener suma positiva.")
  }

  # Normalización suave por robustez frente a redondeo.
  Y_mat <- sweep(Y_mat, 1, row_sums, "/")

  if (any(abs(rowSums(Y_mat) - 1) > 1e-6)) {
    stop("Las filas de `Y_mat` deben sumar 1.")
  }

  if (any(Y_mat <= 0)) {
    stop(
      "La likelihood Dirichlet requiere composiciones estrictamente positivas. ",
      "Hay ceros en `Y_mat`; considera aplicar un ajuste/pseudoconteo antes."
    )
  }

  if (length(cond) != nrow(Y_mat)) {
    stop("`cond` debe tener longitud N.")
  }

  if (anyNA(cond)) {
    stop("`cond` contiene NA.")
  }

  if (!is.matrix(Delta) && !is.data.frame(Delta)) {
    stop("`Delta` debe ser una matriz K x K.")
  }
  Delta <- as.matrix(Delta)

  if (nrow(Delta) != ncol(Y_mat) || ncol(Delta) != ncol(Y_mat)) {
    stop("`Delta` debe tener dimensión K x K, donde K = ncol(Y_mat).")
  }

  if (!is.numeric(rank_Delta) || length(rank_Delta) != 1 || rank_Delta < 0) {
    stop("`rank_Delta` debe ser un número no negativo.")
  }

  if (!is.numeric(s) || length(s) != 1 || s <= 0) {
    stop("`s` debe ser un escalar positivo.")
  }

  # Convertimos condición a índices enteros 1, ..., J.
  cond_vec <- as.integer(as.factor(cond))

  stan_data <- list(
    K = ncol(Y_mat),
    N = nrow(Y_mat),
    Y = Y_mat,
    cond = cond_vec,
    Delta = Delta,
    rank_Delta = rank_Delta,
    s = s,
    prior_only = as.integer(prior_only)
  )

  # Algunos modelos epistáticos condicionales usan explícitamente J.
  # Incluirlo aquí no perjudica a los modelos que no lo declaran en Stan,
  # siempre que el backend usado acepte datos extra. Si se usa rstan/cmdstanr
  # con validación estricta, elimina J para modelos que no lo declaren.
  stan_data$J <- max(cond_vec)

  # -----------------------------
  # Día experimental opcional
  # -----------------------------
  if (!is.null(dia)) {
    if (length(dia) != nrow(Y_mat)) {
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
