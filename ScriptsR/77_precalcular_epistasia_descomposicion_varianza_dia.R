# 77_precalcular_epistasia_descomposicion_varianza_dia.R
# -------------------------------------------------------------------------
# Objetivo
# -------------------------------------------------------------------------
# Este script precalcula los objetos necesarios para analizar la
# descomposición aditiva/epistática media por condición a partir del
# modelo epistático con efecto aleatorio de día.
#
# La idea es separar los cálculos pesados del renderizado del documento.
# El .qmd debería limitarse a:
#   1) cargar el .rds generado aquí,
#   2) construir tablas,
#   3) construir gráficos,
#   4) interpretar los resultados.
#
# Salida principal:
#
#   DatosProcesados/epistasia_descomp_dia/77_resultados_epistasia_descomposicion_varianza_dia.rds
#
# El objeto guardado es una lista con tablas y data.frames ya preparados para
# ser usados directamente en Quarto.
#
# -------------------------------------------------------------------------
# Entradas esperadas
# -------------------------------------------------------------------------
# Se asume que existen estos archivos:
#
#   DatosProcesados/<archivo que contiene fit_kappa_comun_dia>.Rdata
#   DatosProcesados/df_h3917.Rdata
#   DatosProcesados/delta_obj_h3917.Rdata
#
# Y que contienen, respectivamente:
#
#   - fit_kappa_comun_dia:
#       ajuste Stan del modelo epistático con kappa común y efecto
#       aleatorio de día.
#
# Nota importante:
#   Este análisis usa directamente `phi`, no `phi + day_eff`, porque el
#   objetivo es estudiar la epistasia media por condición ajustada por día,
#   no la epistasia específica de cada día.
#
#   - df.ancho:
#       data.frame ancho con columnas:
#         Estado, Dia, y después las 20 proteoformas de H3_9_17.
#
#   - delta_obj:
#       lista con la matriz de contrastes locales B y sus etiquetas B_label.
#
# -------------------------------------------------------------------------
# Paquetes
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(rstan)
  library(tidyverse)
  library(here)
})

# -------------------------------------------------------------------------
# Configuración de rutas
# -------------------------------------------------------------------------

dir_datos <- here::here("DatosProcesados")

# Carpeta de salida específica para el modelo con efecto aleatorio de día.
# Así evitamos sobrescribir los resultados generados a partir del otro modelo.
dir_salida <- file.path(dir_datos, "epistasia_descomp_dia")

# Primero comprobamos/creamos las carpetas necesarias.
dir.create(dir_datos, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_salida, showWarnings = FALSE, recursive = TRUE)

# Archivo que contiene el objeto stanfit `fit_kappa_comun_dia`.
# Si tu archivo tiene otro nombre, cambia solo esta línea.
archivo_fit <- file.path(dir_datos, "stan_epis_kappa_comun_dia_h3917.Rdata")

archivo_datos     <- file.path(dir_datos, "df_h3917.Rdata")
archivo_delta_obj <- file.path(dir_datos, "delta_obj_h3917.Rdata")

archivo_salida <- file.path(
  dir_salida,
  "77_resultados_epistasia_descomposicion_varianza_dia.rds"
)

# -------------------------------------------------------------------------
# Funciones auxiliares generales
# -------------------------------------------------------------------------

comprobar_archivo <- function(path) {
  if (!file.exists(path)) {
    stop("No se encuentra el archivo esperado: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

comprobar_objeto <- function(nombre_objeto, entorno = parent.frame()) {
  if (!exists(nombre_objeto, envir = entorno, inherits = FALSE)) {
    stop("No se encuentra el objeto esperado: `", nombre_objeto, "`.", call. = FALSE)
  }
  invisible(TRUE)
}

localizar_archivo_fit <- function(archivo_preferido, dir_datos, objeto_fit) {
  # Devuelve el archivo que contiene el objeto stanfit esperado.
  # Primero intenta usar el archivo preferido. Si no existe, busca entre los
  # .Rdata de DatosProcesados hasta encontrar uno que contenga `objeto_fit`.
  if (file.exists(archivo_preferido)) {
    return(archivo_preferido)
  }

  archivos_rdata <- list.files(
    dir_datos,
    pattern = "\\.Rdata$",
    full.names = TRUE
  )

  for (archivo in archivos_rdata) {
    env_tmp <- new.env(parent = emptyenv())
    obj_names <- load(archivo, envir = env_tmp)
    if (objeto_fit %in% obj_names) {
      message("Archivo preferido no encontrado. Usando archivo detectado: ", archivo)
      return(archivo)
    }
    rm(env_tmp)
    gc(verbose = FALSE)
  }

  stop(
    "No se encontró ningún archivo .Rdata en ", dir_datos,
    " que contenga el objeto `", objeto_fit, "`. ",
    "Edita `archivo_fit` para apuntar al archivo correcto.",
    call. = FALSE
  )
}

row_vars <- function(x) {
  # Varianza por fila.
  #
  # Equivale a:
  #   apply(x, 1, var)
  #
  # pero evita el coste de apply() en matrices grandes.
  n <- ncol(x)
  mu <- rowMeans(x)
  rowSums((x - mu)^2) / (n - 1)
}

col_vars <- function(x) {
  # Varianza por columna.
  #
  # Equivale a:
  #   apply(x, 2, var)
  #
  # pero evita el coste de apply() en matrices grandes.
  n <- nrow(x)
  mu <- colMeans(x)
  colSums((x - rep(mu, each = n))^2) / (n - 1)
}

# -------------------------------------------------------------------------
# Diseño aditivo de proteoformas
# -------------------------------------------------------------------------

parsear_diseno_h3917 <- function(nombres) {
  # Construye un data.frame con los estados marginales de cada proteoforma.
  #
  # Para cada proteoforma de H3_9_17 se codifican tres factores:
  #
  #   K9  : unmod, me1, me2, me3, ac
  #   S10 : unmod, ph
  #   K14 : unmod, ac
  #
  # Esta codificación define el modelo aditivo:
  #
  #   phi ~ K9 + S10 + K14
  #
  # frente al que se mide el residuo epistático.
  tibble::tibble(
    Estado = nombres,
    K9 = factor(
      dplyr::case_when(
        stringr::str_detect(nombres, "K9me1") ~ "me1",
        stringr::str_detect(nombres, "K9me2") ~ "me2",
        stringr::str_detect(nombres, "K9me3") ~ "me3",
        stringr::str_detect(nombres, "K9ac")  ~ "ac",
        TRUE                                  ~ "unmod"
      ),
      levels = c("unmod", "me1", "me2", "me3", "ac")
    ),
    S10 = factor(
      dplyr::if_else(stringr::str_detect(nombres, "S10ph"), "ph", "unmod"),
      levels = c("unmod", "ph")
    ),
    K14 = factor(
      dplyr::if_else(stringr::str_detect(nombres, "K14ac"), "ac", "unmod"),
      levels = c("unmod", "ac")
    )
  )
}

construir_diseno_aditivo <- function(df_diseno) {
  # Construye la matriz del modelo aditivo.
  #
  # Se usan contrastes de suma cero para que la descomposición sea interpretable
  # como una proyección sobre efectos marginales:
  #
  #   intercepto + efecto K9 + efecto S10 + efecto K14
  #
  # La matriz P permite calcular los coeficientes por mínimos cuadrados:
  #
  #   beta = P %*% phi
  #
  # donde:
  #
  #   P = (X'X)^(-1) X'
  old_contrasts <- options("contrasts")
  on.exit(options(old_contrasts), add = TRUE)

  options(contrasts = c("contr.sum", "contr.poly"))

  X <- stats::model.matrix(~ K9 + S10 + K14, data = df_diseno)
  P <- solve(crossprod(X), t(X))

  list(
    X = X,
    P = P,
    columnas = colnames(X)
  )
}

# -------------------------------------------------------------------------
# Descomposición global de varianza
# -------------------------------------------------------------------------

resumir_varianza_global <- function(var_comp,
                                    cond_labels = c("Asíncrona", "Mitosis")) {
  # Resume la fracción de varianza explicada por:
  #
  #   K9, S10, K14 y componente epistático.
  #
  # var_comp tiene dimensiones:
  #
  #   iteraciones x componentes x condiciones
  #
  # donde los componentes son:
  #
  #   V_Total, V_K9, V_S10, V_K14, V_Epi
  frac_var <- sweep(
    var_comp[, 2:5, , drop = FALSE],
    c(1, 3),
    var_comp[, 1, ],
    "/"
  )

  tibble::tibble(
    Componente = c("K9", "S10", "K14", "Epistasia"),
    Asincrono = apply(frac_var[, , 1], 2, mean) * 100,
    SD_Asincrono = apply(frac_var[, , 1], 2, sd) * 100,
    Mitosis = apply(frac_var[, , 2], 2, mean) * 100,
    SD_Mitosis = apply(frac_var[, , 2], 2, sd) * 100
  )
}

descomponer_phi_vectorizado <- function(phi_samples,
                                        X,
                                        P,
                                        nombres,
                                        cond_labels = c("Asíncrona", "Mitosis"),
                                        guardar_draws = FALSE) {
  # Descompone los campos phi posteriores en:
  #
  #   phi_total = phi_aditivo + phi_epistasia
  #
  # donde phi_aditivo es la proyección de phi sobre el modelo aditivo:
  #
  #   K9 + S10 + K14
  #
  # y phi_epistasia es el residuo:
  #
  #   phi_epistasia = phi - phi_aditivo
  #
  # Esta función está vectorizada por condición:
  #
  #   Betas  = Phi %*% t(P)
  #   Phi_ad = Betas %*% t(X)
  #
  # evitando bucles draw-a-draw.
  #
  # Entrada:
  #   phi_samples:
  #     array con dimensiones:
  #       iteraciones x proteoformas x condiciones
  #
  #   X:
  #     matriz de diseño aditivo.
  #
  #   P:
  #     matriz de proyección para obtener coeficientes aditivos.
  #
  #   nombres:
  #     nombres de las proteoformas.
  #
  # Salida:
  #   lista con:
  #     var_comp,
  #     resumen_varianza,
  #     resumen_estados,
  #     analisis_por_estado,
  #     df_probs,
  #     df_delta,
  #     opcionalmente phi_aditivo y phi_epistasia.
  n_iter  <- dim(phi_samples)[1]
  n_prot  <- dim(phi_samples)[2]
  n_cond  <- dim(phi_samples)[3]

  stopifnot(ncol(X) == nrow(P))
  stopifnot(nrow(X) == n_prot)
  stopifnot(length(nombres) == n_prot)
  stopifnot(n_cond >= 2)

  var_comp <- array(
    NA_real_,
    dim = c(n_iter, 5, n_cond),
    dimnames = list(
      NULL,
      c("V_Total", "V_K9", "V_S10", "V_K14", "V_Epi"),
      cond_labels[seq_len(n_cond)]
    )
  )

  analisis_por_estado <- vector("list", n_cond)
  names(analisis_por_estado) <- cond_labels[seq_len(n_cond)]

  if (guardar_draws) {
    phi_aditivo   <- array(NA_real_, dim = dim(phi_samples))
    phi_epistasia <- array(NA_real_, dim = dim(phi_samples))
  } else {
    phi_aditivo   <- NULL
    phi_epistasia <- NULL
  }

  epi_mean <- epi_sd <- prob_pos <- prob_neg <- matrix(
    NA_real_,
    nrow = n_prot,
    ncol = n_cond,
    dimnames = list(nombres, cond_labels[seq_len(n_cond)])
  )

  for (cc in seq_len(n_cond)) {
    Phi <- phi_samples[, , cc, drop = TRUE]  # iteraciones x proteoformas

    # Coeficientes del modelo aditivo para todos los draws a la vez.
    Betas <- Phi %*% t(P)                    # iteraciones x columnas_X

    # Componentes marginales reconstruidas en el espacio de proteoformas.
    # La columna 1 de X es el intercepto.
    comp_K9 <- Betas[, 2:5, drop = FALSE] %*%
      t(X[, 2:5, drop = FALSE])

    comp_S10 <- Betas[, 6, drop = FALSE] %*%
      t(X[, 6, drop = FALSE])

    comp_K14 <- Betas[, 7, drop = FALSE] %*%
      t(X[, 7, drop = FALSE])

    Phi_ad <- Betas %*% t(X)
    Phi_epi <- Phi - Phi_ad

    if (guardar_draws) {
      phi_aditivo[, , cc] <- Phi_ad
      phi_epistasia[, , cc] <- Phi_epi
    }

    # Varianza global por draw.
    var_comp[, "V_Total", cc] <- row_vars(Phi)
    var_comp[, "V_K9", cc]    <- row_vars(comp_K9)
    var_comp[, "V_S10", cc]   <- row_vars(comp_S10)
    var_comp[, "V_K14", cc]   <- row_vars(comp_K14)
    var_comp[, "V_Epi", cc]   <- row_vars(Phi_epi)

    # Varianza atribuible a cada componente por proteoforma.
    analisis_por_estado[[cc]] <- tibble::tibble(
      Estado = nombres,
      V_K9 = col_vars(comp_K9),
      V_S10 = col_vars(comp_S10),
      V_K14 = col_vars(comp_K14),
      V_Epi = col_vars(Phi_epi),
      V_Total = col_vars(Phi)
    )

    # Resumen posterior del residuo epistático por proteoforma.
    epi_mean[, cc] <- colMeans(Phi_epi)
    epi_sd[, cc]   <- apply(Phi_epi, 2, sd)
    prob_pos[, cc] <- colMeans(Phi_epi > 0)
    prob_neg[, cc] <- colMeans(Phi_epi < 0)
  }

  resumen_estados <- tidyr::expand_grid(
    Marcas = nombres,
    Condicion = cond_labels[seq_len(n_cond)]
  ) %>%
    dplyr::mutate(
      epi_mean = as.vector(epi_mean),
      epi_sd = as.vector(epi_sd),
      prob_pos = as.vector(prob_pos),
      prob_neg = as.vector(prob_neg)
    )

  # Tabla ligera para los gráficos de probabilidad de epistasia positiva.
  # Se mantienen nombres compatibles con el qmd original.
  df_probs <- tibble::tibble(
    Marcas = nombres,
    prop_gt_0_asinc = prob_pos[, 1],
    prob_gt_0_mitot = prob_pos[, 2],
    prop_lt_0_asinc = prob_neg[, 1],
    prop_lt_0_mitot = prob_neg[, 2]
  )

  # Tabla de cambio de probabilidad de epistasia positiva entre condiciones.
  df_delta <- df_probs %>%
    dplyr::mutate(
      delta_prob = prob_gt_0_mitot - prop_gt_0_asinc,
      abs_delta = abs(delta_prob),
      log2ratio = log2((prob_gt_0_mitot + 1e-6) / (prop_gt_0_asinc + 1e-6)),
      cambio = dplyr::case_when(
        delta_prob > 0.4  ~ "Gana Sinergia (Mitosis)",
        delta_prob < -0.4 ~ "Pierde Sinergia (Hacia Aditividad)",
        TRUE ~ "Estable"
      ),
      twofold = dplyr::case_when(
        log2ratio > 1  ~ "Gana sinergia (Mitosis)",
        log2ratio < -1 ~ "Hacia aditividad (Mitosis)",
        TRUE ~ "Estable"
      )
    )

  list(
    var_comp = var_comp,
    resumen_varianza = resumir_varianza_global(var_comp, cond_labels),
    resumen_estados = resumen_estados,
    analisis_por_estado = analisis_por_estado,
    df_probs = df_probs,
    df_delta = df_delta,
    phi_aditivo = phi_aditivo,
    phi_epistasia = phi_epistasia
  )
}

# -------------------------------------------------------------------------
# Índice de complejidad epistática por condición
# -------------------------------------------------------------------------

calcular_complejidad_local <- function(analisis_por_estado,
                                       nombres,
                                       cond_labels = c("Asíncrona", "Mitosis")) {
  # Calcula el índice de complejidad epistática LOCAL agregado.
  #
  # Este índice se calcula primero por proteoforma:
  #
  #   Prop_Epi_Local_j = Var_post(Phi_epi[, j]) / Var_post(Phi[, j])
  #
  # y después se agrega sumando sobre proteoformas:
  #
  #   Indice_Local = sum_j Prop_Epi_Local_j
  #
  # IMPORTANTE:
  #   - Es un índice local agregado, no un porcentaje global de varianza.
  #   - Puede ser mayor que 1, porque se calcula a partir de varianzas
  #     posteriores por proteoforma, no de la descomposición ortogonal global
  #     del paisaje energético.
  #   - Por eso se mantiene en escala de índice y no se multiplica por 100.

  stopifnot(length(analisis_por_estado) >= 2)
  stopifnot(length(nombres) == nrow(analisis_por_estado[[1]]))

  detalle_async <- analisis_por_estado[[1]] %>%
    dplyr::mutate(
      Prop_Epi_Local = V_Epi / V_Total,
      Prop_Adit_Local = (V_K9 + V_S10 + V_K14) / V_Total
    )

  detalle_mitot <- analisis_por_estado[[2]] %>%
    dplyr::mutate(
      Prop_Epi_Local = V_Epi / V_Total,
      Prop_Adit_Local = (V_K9 + V_S10 + V_K14) / V_Total
    )

  indice_local_async <- sum(detalle_async$Prop_Epi_Local, na.rm = TRUE)
  indice_local_mitot <- sum(detalle_mitot$Prop_Epi_Local, na.rm = TRUE)

  tabla_detalle_prop <- tibble::tibble(
    Estado = nombres,
    Prop_Epi_Local_Async = detalle_async$Prop_Epi_Local,
    Prop_Epi_Local_Mitot = detalle_mitot$Prop_Epi_Local,
    Cambio_Complejidad_Local = Prop_Epi_Local_Mitot - Prop_Epi_Local_Async,
    Prop_Adit_Local_Async = detalle_async$Prop_Adit_Local,
    Prop_Adit_Local_Mitot = detalle_mitot$Prop_Adit_Local
  )

  list(
    indice_local = c(
      Asincrona = indice_local_async,
      Mitosis = indice_local_mitot
    ),
    tabla_detalle_prop = tabla_detalle_prop
  )
}

calcular_complejidad_global <- function(var_comp,
                                        cond_labels = c("Asíncrona", "Mitosis")) {
  # Calcula métricas globales usando la descomposición ortogonal por draw.
  #
  # Estas métricas sí son porcentajes globales interpretables porque se calculan
  # dentro de cada draw posterior como:
  #
  #   Peso_epistatico_global = V_Epi / V_Total
  #   Aditividad_global      = (V_K9 + V_S10 + V_K14) / V_Total
  #
  # con varianzas medidas sobre la superficie energética completa.

  componentes_necesarios <- c("V_Total", "V_K9", "V_S10", "V_K14", "V_Epi")
  if (!all(componentes_necesarios %in% dimnames(var_comp)[[2]])) {
    stop(
      "`var_comp` debe contener los componentes: ",
      paste(componentes_necesarios, collapse = ", "),
      call. = FALSE
    )
  }

  frac_epi <- var_comp[, "V_Epi", ] / var_comp[, "V_Total", ]

  frac_adit <- (
    var_comp[, "V_K9", ] +
      var_comp[, "V_S10", ] +
      var_comp[, "V_K14", ]
  ) / var_comp[, "V_Total", ]

  # Comprobación de cierre de la descomposición global.
  # En una proyección ortogonal, aditividad + epistasia debería cerrar a 1.
  error_cierre <- max(abs(frac_epi + frac_adit - 1), na.rm = TRUE)
  if (is.finite(error_cierre) && error_cierre > 1e-6) {
    warning(
      "La descomposición aditiva + epistática no cierra exactamente a 1. ",
      "Máximo error de cierre: ", signif(error_cierre, 3),
      call. = FALSE
    )
  }

  tabla_complejidad_draws <- tibble::tibble(
    Iteracion = rep(seq_len(dim(var_comp)[1]), times = length(cond_labels)),
    Condicion = rep(cond_labels, each = dim(var_comp)[1]),
    Peso_Epistatico_Global = as.vector(frac_epi),
    Grado_Aditividad_Global = as.vector(frac_adit)
  )

  list(
    peso_epi_global = c(
      Asincrona = mean(frac_epi[, 1], na.rm = TRUE) * 100,
      Mitosis = mean(frac_epi[, 2], na.rm = TRUE) * 100
    ),
    aditividad_global = c(
      Asincrona = mean(frac_adit[, 1], na.rm = TRUE) * 100,
      Mitosis = mean(frac_adit[, 2], na.rm = TRUE) * 100
    ),
    tabla_complejidad_draws = tabla_complejidad_draws,
    error_cierre_descomposicion = error_cierre
  )
}

resumir_complejidad <- function(var_comp,
                                analisis_por_estado,
                                nombres,
                                cond_labels = c("Asíncrona", "Mitosis")) {
  # Construye la tabla sintética final para el informe.
  #
  # La Tabla 4 combina dos escalas complementarias:
  #
  #   1) Índice de Complejidad Local:
  #      suma de Prop_Epi_Local por proteoforma. Es un índice agregado, no un %.
  #
  #   2) Peso epistático global (%) y grado de aditividad global (%):
  #      calculados desde var_comp, es decir, desde la descomposición global
  #      ortogonal de la varianza del paisaje energético.

  complejidad_local <- calcular_complejidad_local(
    analisis_por_estado = analisis_por_estado,
    nombres = nombres,
    cond_labels = cond_labels
  )

  complejidad_global <- calcular_complejidad_global(
    var_comp = var_comp,
    cond_labels = cond_labels
  )

  tabla_complejidad_corregida <- tibble::tibble(
    Metrica = c(
      "Índice de Complejidad Local (sum Prop_Epi por proteoforma)",
      "Peso epistático global (%)",
      "Grado de aditividad global (%)"
    ),
    Asincrona = c(
      complejidad_local$indice_local[["Asincrona"]],
      complejidad_global$peso_epi_global[["Asincrona"]],
      complejidad_global$aditividad_global[["Asincrona"]]
    ),
    Mitosis = c(
      complejidad_local$indice_local[["Mitosis"]],
      complejidad_global$peso_epi_global[["Mitosis"]],
      complejidad_global$aditividad_global[["Mitosis"]]
    )
  ) %>%
    dplyr::mutate(
      Delta_Absoluto = Mitosis - Asincrona
    )

  list(
    tabla_complejidad_corregida = tabla_complejidad_corregida,
    tabla_detalle_prop = complejidad_local$tabla_detalle_prop,
    tabla_complejidad_draws = complejidad_global$tabla_complejidad_draws,
    error_cierre_descomposicion = complejidad_global$error_cierre_descomposicion
  )
}

# -------------------------------------------------------------------------
# Epistasia local por cuadrados
# -------------------------------------------------------------------------

identificar_estados <- function(fila_B, nombres_phi) {
  # Traduce una fila de la matriz B en una fórmula legible.
  #
  # Si una fila de B contiene coeficientes:
  #
  #   +1, -1, -1, +1
  #
  # se representa como:
  #
  #   (estado_a + estado_d) - (estado_b + estado_c)
  pos_positivas <- nombres_phi[which(fila_B == 1)]
  pos_negativas <- nombres_phi[which(fila_B == -1)]

  paste0(
    "(", paste(pos_positivas, collapse = " + "), ") - (",
    paste(pos_negativas, collapse = " + "), ")"
  )
}

resumir_epistasia_local <- function(phi_samples,
                                    delta_obj,
                                    cond_labels = c("Asíncrona", "Mitosis")) {
  # Calcula los contrastes locales de epistasia:
  #
  #   epsilon = B %*% phi
  #
  # para cada draw posterior y cada condición.
  #
  # B codifica los cuadrados epistáticos del grafo de proteoformas.
  # La salida es una tabla resumida por cuadrado y condición.
  B <- delta_obj$B
  B_label <- delta_obj$B_label

  if (is.null(B) || is.null(B_label)) {
    stop("`delta_obj` debe contener los elementos `B` y `B_label`.", call. = FALSE)
  }

  nombres_cuadrados <- rownames(B)
  nombres_unicos <- make.unique(nombres_cuadrados)
  nombres_estados_phi <- colnames(B)

  mapa_cuadrados <- tibble::tibble(
    Cuadrado = nombres_unicos,
    Formula = apply(B, 1, identificar_estados, nombres_phi = nombres_estados_phi),
    Etiqueta = B_label
  )

  df_epistasia <- purrr::map_dfr(seq_along(cond_labels), function(cc) {
    eps <- phi_samples[, , cc, drop = TRUE] %*% t(B)

    tibble::as_tibble(eps, .name_repair = ~ nombres_unicos) %>%
      dplyr::mutate(
        Iteracion = dplyr::row_number(),
        Condicion = cond_labels[[cc]]
      ) %>%
      tidyr::pivot_longer(
        cols = -c(Iteracion, Condicion),
        names_to = "Cuadrado",
        values_to = "Epsilon"
      )
  })

  df_epistasia %>%
    dplyr::left_join(mapa_cuadrados, by = "Cuadrado") %>%
    dplyr::mutate(
      Modulo_Constante = stringr::str_extract(Cuadrado, "(?<=\\|)[^.]+"),
      Tipo_Interaccion = dplyr::case_when(
        stringr::str_detect(Modulo_Constante, "K14") ~ "Eje K9-S10",
        stringr::str_detect(Modulo_Constante, "S10") ~ "Eje K9-K14",
        stringr::str_detect(Modulo_Constante, "K9")  ~ "Eje S10-K14",
        TRUE ~ "Otra"
      )
    ) %>%
    dplyr::group_by(
      Condicion,
      Cuadrado,
      Etiqueta,
      Formula,
      Tipo_Interaccion
    ) %>%
    dplyr::summarise(
      Media = mean(Epsilon),
      IC_Lower = stats::quantile(Epsilon, 0.025),
      IC_Upper = stats::quantile(Epsilon, 0.975),
      Prob_Positiva = mean(Epsilon > 0),
      .groups = "drop"
    )
}

# -------------------------------------------------------------------------
# Ejecución principal
# -------------------------------------------------------------------------

message("== Precalculando descomposición de varianza epistática: modelo con día ==")

archivo_fit <- localizar_archivo_fit(
  archivo_preferido = archivo_fit,
  dir_datos = dir_datos,
  objeto_fit = "fit_kappa_comun_dia"
)

comprobar_archivo(archivo_fit)
comprobar_archivo(archivo_datos)
comprobar_archivo(archivo_delta_obj)

load(archivo_fit)
load(archivo_datos)
load(archivo_delta_obj)

comprobar_objeto("fit_kappa_comun_dia")
comprobar_objeto("df.ancho")
comprobar_objeto("delta_obj")

# -------------------------------------------------------------------------
# Extracción posterior
# -------------------------------------------------------------------------

message("Extrayendo draws posteriores de phi...")

fit <- fit_kappa_comun_dia

# Para epistasia media por condición usamos `phi` directamente.
# No sumamos `day_eff`, porque ese término representa desviaciones
# experimentales día-a-día y no la arquitectura epistática media.
phi_samples <- rstan::extract(fit, pars = "phi")$phi

if (length(dim(phi_samples)) != 3) {
  stop(
    "`phi_samples` debería tener dimensiones: iteraciones x proteoformas x condiciones.",
    call. = FALSE
  )
}

# En df.ancho se asume:
#   columna 1: Estado
#   columna 2: Dia
#   columnas 3:22: composiciones/proteoformas
nombres <- names(df.ancho[, 3:22])

if (length(nombres) != dim(phi_samples)[2]) {
  stop(
    "El número de proteoformas en df.ancho no coincide con la dimensión de phi.",
    call. = FALSE
  )
}

# -------------------------------------------------------------------------
# Construcción de diseño aditivo
# -------------------------------------------------------------------------

message("Construyendo matriz de diseño aditivo K9 + S10 + K14...")

df_diseno <- parsear_diseno_h3917(nombres)
diseno <- construir_diseno_aditivo(df_diseno)

# -------------------------------------------------------------------------
# Descomposición de phi
# -------------------------------------------------------------------------

message("Calculando descomposición aditiva/epistática...")

descomp <- descomponer_phi_vectorizado(
  phi_samples = phi_samples,
  X = diseno$X,
  P = diseno$P,
  nombres = nombres,
  cond_labels = c("Asíncrona", "Mitosis"),
  guardar_draws = FALSE
)

# -------------------------------------------------------------------------
# Tablas derivadas
# -------------------------------------------------------------------------

message("Calculando índices globales de complejidad epistática...")

complejidad <- resumir_complejidad(
  var_comp = descomp$var_comp,
  analisis_por_estado = descomp$analisis_por_estado,
  nombres = nombres,
  cond_labels = c("Asíncrona", "Mitosis")
)

message("Calculando epistasia local por cuadrados...")

df_resumen_epistasia_local <- resumir_epistasia_local(
  phi_samples = phi_samples,
  delta_obj = delta_obj,
  cond_labels = c("Asíncrona", "Mitosis")
)

# -------------------------------------------------------------------------
# Objeto final
# -------------------------------------------------------------------------

resultados_epistasia_varianza <- c(
  list(
    metadata = list(
      script = "77_precalcular_epistasia_descomposicion_varianza_dia.R",
      modelo = "fit_kappa_comun_dia",
      fecha = as.character(Sys.time()),
      archivo_fit = archivo_fit,
      archivo_datos = archivo_datos,
      archivo_delta_obj = archivo_delta_obj,
      archivo_salida = archivo_salida,
      dim_phi_samples = dim(phi_samples)
    ),
    nombres = nombres,
    df_diseno = df_diseno,
    X = diseno$X,
    columnas_X = diseno$columnas
  ),
  descomp[setdiff(names(descomp), c("phi_aditivo", "phi_epistasia"))],
  complejidad,
  list(
    df_resumen_epistasia_local = df_resumen_epistasia_local
  )
)

# -------------------------------------------------------------------------
# Guardado
# -------------------------------------------------------------------------

saveRDS(
  resultados_epistasia_varianza,
  file = archivo_salida
)

message("Guardado correctamente en: ", archivo_salida)
message("Objetos disponibles en el .rds:")

print(names(resultados_epistasia_varianza))

message("== Fin del precálculo ==")
