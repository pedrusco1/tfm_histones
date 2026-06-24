# ==============================================================================
# 32_expresion_diferencial.R
# ------------------------------------------------------------------------------
# Precalcula los análisis diferenciales de PTMs del fragmento H3(9-17)
# usados en 32_resultados_bioinfo.qmd.
#
# Objetivo:
#   - Sacar del QMD los cálculos pesados/repetidos.
#   - Guardar un único objeto .rds con las tablas y datos necesarios para
#     renderizar rápidamente tablas y figuras.
#
# Entrada esperada:
#   DatosProcesados/df_h3917.Rdata
#     Debe contener el objeto df.ancho con columnas:
#       Estado, Dia, y las 20 proteoformas H3(9-17) en las columnas 3:22.
#
# Salida:
#   DatosProcesados/32_resultados_bioinfo.rds
# ==============================================================================

# --- 1. Librerías --------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(edgeR)
  library(compositions)
  library(dplyr)
  library(tibble)
  library(purrr)
})

# --- 2. Rutas y carga de datos -------------------------------------------------

archivo_datos <- here::here("DatosProcesados", "df_h3917.Rdata")
archivo_salida <- here::here("DatosProcesados", "32_resultados_bioinfo.rds")

if (!file.exists(archivo_datos)) {
  stop(
    "No se encuentra el archivo de datos: ", archivo_datos, "\n",
    "Comprueba que existe DatosProcesados/df_h3917.Rdata."
  )
}

load(archivo_datos)

if (!exists("df.ancho")) {
  stop("El archivo ", archivo_datos, " no contiene el objeto `df.ancho`.")
}

if (!all(c("Estado", "Dia") %in% names(df.ancho))) {
  stop("`df.ancho` debe contener al menos las columnas `Estado` y `Dia`.")
}

if (ncol(df.ancho) < 22) {
  stop("`df.ancho` debe contener las 20 proteoformas en las columnas 3:22.")
}

# --- 3. Preparación común ------------------------------------------------------

# En el documento original se usan las columnas 3:22 como proteoformas.
marcas <- colnames(df.ancho)[3:22]
data_prop <- df.ancho[, marcas]
matriz_prop <- as.matrix(data_prop)
storage.mode(matriz_prop) <- "double"

if (any(is.na(matriz_prop))) {
  stop("Hay valores NA en las proporciones de proteoformas.")
}

if (any(matriz_prop < 0)) {
  stop("Hay proporciones negativas en las columnas de proteoformas.")
}

grupo <- factor(df.ancho$Estado)
condicion <- grupo

# Matriz de diseño única para modelos con intercepto + contraste mitótico.
design_mat <- model.matrix(~ condicion)
colnames(design_mat) <- c("Intercepto", "MitoticoVsAsinc")

# --- 4. Modelo 1: limma sobre coordenadas CLR ---------------------------------

# Se añade pseudocuenta pequeña para evitar log(0) antes de la transformación CLR,
# manteniendo el criterio usado en el QMD original.
data_clr <- as.matrix(compositions::clr(matriz_prop + 1e-5))
rownames(data_clr) <- seq_len(nrow(data_clr))
colnames(data_clr) <- marcas

fit_limma <- limma::lmFit(t(data_clr), design_mat)
fit_limma <- limma::eBayes(fit_limma, trend = TRUE, robust = TRUE)

resultados_limma <- limma::topTable(
  fit_limma,
  coef = "MitoticoVsAsinc",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Versión estandarizada con la proteoforma como columna explícita.
resultados_limma_tbl <- resultados_limma %>%
  tibble::rownames_to_column(var = "PTM")

volcano_limma <- resultados_limma_tbl %>%
  mutate(
    neglog10FDR = -log10(adj.P.Val),
    significativa = adj.P.Val < 0.05
  )

# --- 5. Modelo 2: edgeR sobre pseudo-conteos ----------------------------------

# El QMD original construye pseudo-conteos multiplicando proporciones por 1e6.
# Esto permite aplicar edgeR/voom de forma comparativa, pero debe interpretarse
# como análisis sobre pseudo-conteos derivados de abundancias relativas.
conteos <- t(matriz_prop) * 1e6
rownames(conteos) <- marcas

obj_edgeR <- edgeR::DGEList(counts = conteos, group = condicion)
obj_edgeR <- edgeR::calcNormFactors(obj_edgeR, method = "TMM")
obj_edgeR <- edgeR::estimateDisp(obj_edgeR, design_mat)

fit_edgeR_qlf <- edgeR::glmQLFit(obj_edgeR, design_mat, robust = TRUE)
qlf <- edgeR::glmQLFTest(fit_edgeR_qlf, coef = 2)

resultados_edgeR <- edgeR::topTags(qlf, n = Inf)$table

resultados_edgeR_tbl <- resultados_edgeR %>%
  tibble::rownames_to_column(var = "PTM")

volcano_edgeR <- resultados_edgeR_tbl %>%
  mutate(
    neglog10FDR = -log10(FDR),
    significativa = FDR < 0.05
  )

# --- 6. Modelo 3: voom-limma ---------------------------------------------------

voom_obj <- limma::voom(obj_edgeR, design_mat, plot = FALSE)
fit_voom <- limma::lmFit(voom_obj, design_mat)
fit_voom <- limma::eBayes(fit_voom)

resultados_voom <- limma::topTable(
  fit_voom,
  coef = 2,
  number = Inf,
  sort.by = "P"
) %>%
  tibble::rownames_to_column(var = "PTM")

volcano_voom <- resultados_voom %>%
  mutate(
    neglog10FDR = -log10(adj.P.Val),
    significativa = adj.P.Val < 0.05
  )

# --- 7. Comparación centralizada entre métodos --------------------------------

tabla_comparativa <- resultados_limma_tbl %>%
  select(PTM, logFC_limma = logFC, FDR_limma = adj.P.Val) %>%
  inner_join(
    resultados_edgeR_tbl %>%
      select(PTM, logFC_edgeR = logFC, FDR_edgeR = FDR),
    by = "PTM"
  ) %>%
  inner_join(
    resultados_voom %>%
      select(PTM, logFC_voom = logFC, FDR_voom = adj.P.Val),
    by = "PTM"
  )

lista_significativas <- list(
  `limma-CLR`  = tabla_comparativa$PTM[tabla_comparativa$FDR_limma < 0.05],
  `voom-limma` = tabla_comparativa$PTM[tabla_comparativa$FDR_voom < 0.05],
  `edgeR-NB`   = tabla_comparativa$PTM[tabla_comparativa$FDR_edgeR < 0.05]
)

calcular_jaccard <- function(set1, set2) {
  union_set <- union(set1, set2)
  if (length(union_set) == 0) return(NA_real_)
  length(intersect(set1, set2)) / length(union_set)
}

df_jaccard <- tibble::tibble(
  `Comparación de Modelos` = c(
    "limma-CLR vs limma-voom",
    "limma-CLR vs edgeR-NB",
    "limma-voom vs edgeR-NB"
  ),
  `Índice de Jaccard` = c(
    calcular_jaccard(lista_significativas$`limma-CLR`, lista_significativas$`voom-limma`),
    calcular_jaccard(lista_significativas$`limma-CLR`, lista_significativas$`edgeR-NB`),
    calcular_jaccard(lista_significativas$`voom-limma`, lista_significativas$`edgeR-NB`)
  )
)

resumen_significativas <- tibble::tibble(
  metodo = names(lista_significativas),
  n_significativas = purrr::map_int(lista_significativas, length),
  PTMs_significativas = purrr::map_chr(lista_significativas, ~ paste(.x, collapse = ", "))
)

# --- 8. Objeto final -----------------------------------------------------------

res_bioinfo <- list(
  metadata = list(
    script = "32_expresion_diferencial.R",
    fecha = Sys.time(),
    entrada = archivo_datos,
    salida = archivo_salida,
    fragmento = "H3_9_17",
    n_muestras = nrow(df.ancho),
    n_proteoformas = length(marcas),
    proteoformas = marcas,
    pseudocuenta_clr = 1e-5,
    factor_pseudoconteos = 1e6
  ),
  datos = list(
    grupo = grupo,
    condicion = condicion,
    marcas = marcas,
    matriz_prop = matriz_prop,
    data_clr = data_clr,
    design_mat = design_mat,
    conteos = conteos
  ),
  ajustes = list(
    fit_limma = fit_limma,
    fit_edgeR_qlf = fit_edgeR_qlf,
    qlf = qlf,
    voom_obj = voom_obj,
    fit_voom = fit_voom
  ),
  tablas = list(
    resultados_limma = resultados_limma,
    resultados_limma_tbl = resultados_limma_tbl,
    resultados_edgeR = resultados_edgeR,
    resultados_edgeR_tbl = resultados_edgeR_tbl,
    resultados_voom = resultados_voom,
    tabla_comparativa = tabla_comparativa,
    lista_significativas = lista_significativas,
    resumen_significativas = resumen_significativas,
    df_jaccard = df_jaccard
  ),
  graficas = list(
    volcano_limma = volcano_limma,
    volcano_edgeR = volcano_edgeR,
    volcano_voom = volcano_voom
  )
)

# --- 9. Guardado ---------------------------------------------------------------

dir.create(dirname(archivo_salida), recursive = TRUE, showWarnings = FALSE)
saveRDS(res_bioinfo, archivo_salida)

message("Archivo guardado correctamente en: ", archivo_salida)
message("Objetos principales disponibles en res_bioinfo$tablas y res_bioinfo$graficas.")
