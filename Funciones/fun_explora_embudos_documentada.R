# ==============================================================================
# Funciones ligeras para explorar posibles embudos en modelos Stan
# ==============================================================================
#
# Esta función está pensada para diagnosticar geometrías problemáticas de la
# posterior, especialmente embudos asociados a parámetros de escala y efectos
# latentes, sin saturar la memoria de RStudio.
#
# En lugar de usar mcmc_pairs() sobre todos los draws y muchos parámetros, la
# función submuestrea un número limitado de draws y grafica únicamente los
# parámetros indicados por el usuario.
#
# Ejemplos típicos para modelos con efectos aleatorios:
#   - sigma_day[1] frente a day_raw[1,1]
#   - sigma_day[2] frente a day_eff[1,2]
#   - tau frente a log_delta
#
# La idea es buscar patrones tipo embudo: regiones donde un parámetro de escala
# pequeño fuerza a los efectos latentes a concentrarse fuertemente cerca de cero,
# generando una geometría difícil para HMC/NUTS.
#

#' Diagnóstico ligero de embudos con bayesplot::mcmc_pairs()
#'
#' Genera un gráfico de pares (`mcmc_pairs`) usando solo un subconjunto aleatorio
#' de draws posteriores. Está diseñada para explorar posibles embudos en modelos
#' ajustados con Stan sin consumir tanta memoria como un `mcmc_pairs()` aplicado
#' directamente al objeto completo.
#'
#' La función es especialmente útil para modelos jerárquicos o con efectos
#' aleatorios, donde los embudos suelen aparecer entre un parámetro de escala
#' —por ejemplo `sigma_day[1]`— y efectos latentes asociados —por ejemplo
#' `day_raw[1,1]` o `day_eff[1,1]`—.
#'
#' @param fit Objeto de ajuste de Stan. Normalmente un objeto `stanfit` generado
#'   por `rstan::stan()`. También puede funcionar con otros objetos compatibles
#'   con `posterior::as_draws_df()`.
#'
#' @param pars Vector de caracteres con los nombres exactos de los parámetros
#'   que se quieren representar. Los nombres deben coincidir con los nombres de
#'   columna devueltos por `posterior::as_draws_df(fit)`, por ejemplo
#'   `c("sigma_day[1]", "day_raw[1,1]")`.
#'
#' @param n_draws Número máximo de draws posteriores que se usarán en el gráfico.
#'   Por defecto usa `1000`. Si el ajuste contiene menos draws, se usan todos.
#'   Reducir este valor disminuye el consumo de memoria.
#'
#' @param file Ruta opcional a un archivo `.png`. Si se proporciona, el gráfico se
#'   guarda directamente en disco en lugar de depender del panel gráfico de
#'   RStudio. Esto ayuda a evitar problemas de memoria. Si es `NULL`, el gráfico
#'   se dibuja en el dispositivo gráfico activo.
#'
#' @param width Anchura del archivo PNG en píxeles si `file` no es `NULL`.
#'   Por defecto `1200`.
#'
#' @param height Altura del archivo PNG en píxeles si `file` no es `NULL`.
#'   Por defecto `1000`.
#'
#' @param res Resolución del archivo PNG si `file` no es `NULL`. Por defecto `150`.
#'
#' @return Devuelve invisiblemente el objeto gráfico generado por
#'   `bayesplot::mcmc_pairs()`. La función se usa principalmente por su efecto
#'   secundario: mostrar o guardar el gráfico.
#'
#' @details
#' Para diagnosticar embudos, no conviene graficar todos los parámetros del
#' modelo. Es mejor seleccionar pares pequeños de parámetros con interpretación
#' jerárquica clara, por ejemplo una escala y su efecto aleatorio asociado.
#'
#' En modelos lognormales con efecto de día, pares útiles podrían ser:
#'
#' ```r
#' c("sigma_day[1]", "day_raw[1,1]")
#' c("sigma_day[2]", "day_raw[1,2]")
#' c("sigma_day[1]", "day_eff[1,1]")
#' c("tau", "log_delta")
#' ```
#'
#' Si aparecen divergencias en el modelo, conviene complementar estos gráficos
#' con `bayesplot::mcmc_nuts_divergence()` y `bayesplot::mcmc_nuts_energy()`.
#'
#' @examples
#' \dontrun{
#' diagnostico_funnel_ligero(
#'   fit = fit_lognor_v1_dia,
#'   pars = c("sigma_day[1]", "day_raw[1,1]"),
#'   n_draws = 800,
#'   file = "Resultados/pairs_sigma_day1_dayraw11.png"
#' )
#'
#' diagnostico_funnel_ligero(
#'   fit = fit_lognor_v1,
#'   pars = c("tau", "log_delta"),
#'   n_draws = 800
#' )
#' }
#'
#' @export

diagnostico_funnel_ligero <- function(fit,
                                      pars,
                                      n_draws = 1000,
                                      file = NULL,
                                      width = 1200,
                                      height = 1000,
                                      res = 150) {
  
  # Comprobar que los paquetes necesarios están instalados.
  # Se usa requireNamespace() para evitar cargar paquetes globalmente dentro de
  # la función y para que el mensaje de error sea más claro.
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("El paquete 'posterior' es necesario para usar esta función.")
  }
  
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("El paquete 'bayesplot' es necesario para usar esta función.")
  }
  
  # Validaciones básicas de argumentos.
  if (!is.character(pars) || length(pars) < 2) {
    stop("'pars' debe ser un vector de caracteres con al menos dos parámetros.")
  }
  
  if (!is.numeric(n_draws) || length(n_draws) != 1 || n_draws <= 0) {
    stop("'n_draws' debe ser un número positivo.")
  }
  
  n_draws <- as.integer(n_draws)
  
  # Convertir el ajuste de Stan a un data frame de draws.
  # Las columnas serán parámetros como 'tau', 'sigma_day[1]', 'day_raw[1,1]', etc.
  draws <- posterior::as_draws_df(fit)
  
  # Submuestrear draws para reducir el consumo de memoria.
  # Esto es clave porque mcmc_pairs() puede ser muy pesado si se usa sobre todos
  # los draws o sobre muchos parámetros simultáneamente.
  set.seed(123)
  idx <- sample(seq_len(nrow(draws)), size = min(n_draws, nrow(draws)))
  draws_small <- draws[idx, ]
  
  # Comprobar que todos los parámetros solicitados existen realmente en el ajuste.
  # Esto evita errores poco informativos de bayesplot.
  if (!all(pars %in% names(draws_small))) {
    missing <- setdiff(pars, names(draws_small))
    stop("Estos parámetros no existen en el ajuste: ",
         paste(missing, collapse = ", "))
  }
  
  # Si el usuario proporciona una ruta de archivo, abrir un dispositivo PNG.
  # Esto reduce el riesgo de que RStudio se quede sin memoria manteniendo gráficos
  # grandes en el panel de plots.
  if (!is.null(file)) {
    png(filename = file, width = width, height = height, res = res)
    on.exit({
      dev.off()
      gc()
    }, add = TRUE)
  }
  
  # Construir el gráfico de pares con bayesplot.
  # Se recomienda usar pocos parámetros, idealmente 2-4.
  p <- bayesplot::mcmc_pairs(
    draws_small,
    pars = pars
  )
  
  # Mostrar o enviar el gráfico al dispositivo activo.
  print(p)
  
  # Devolver el gráfico de forma invisible por si el usuario quiere guardarlo o
  # modificarlo posteriormente.
  invisible(p)
}
