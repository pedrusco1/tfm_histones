#----------------------------------------------------------#
# Estudio descriptivo composicional del fragmento H3, 9-17 #
#----------------------------------------------------------#
#
# Este script prepara los resultados descriptivos usados en el informe
# 31_resultados_descri_v1.qmd. La lógica general es:
#   1. Cargar los datos del fragmento H3 9-17.
#   2. Calcular descriptores composicionales básicos.
#   3. Generar gráficos descriptivos.
#   4. Calcular PCA en coordenadas ILR.
#   5. Retroproyectar las cargas del PCA a coordenadas CLR para interpretarlas
#      en términos de proteoformas.
#   6. Guardar todos los resultados en un único objeto RDS.
#
# Nota: se elimina deliberadamente el bloque de diferencias CLR con IC clásico,
# porque con n = 3 por condición es poco informativo y será sustituido por la
# inferencia bayesiana posterior.

#----------------------------------------------------------#
# 0. Librerías                                             #
#----------------------------------------------------------#

library(here)
library(compositions)
library(tidyverse)
library(grid)
library(gridExtra)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggtern)

#----------------------------------------------------------#
# 1. Carga de datos y comprobaciones                       #
#----------------------------------------------------------#

load(here::here("DatosProcesados", "df_h3917.Rdata"))

stopifnot(exists("df.ancho"))
stopifnot(all(c("Estado", "Dia") %in% names(df.ancho)))
stopifnot(ncol(df.ancho) >= 22)

# Variables composicionales: columnas 3:22
ptm_vars <- 3:22
ptm_names <- colnames(df.ancho[, ptm_vars])

# Fijar orden de condiciones
cond <- factor(df.ancho$Estado, levels = c("asinc", "mitot"))
df.ancho$Estado <- cond

# Objeto composicional principal
compo <- acomp(df.ancho[, ptm_vars])

#----------------------------------------------------------#
# 2. Transformación CLR y boxplot por proteoforma           #
#----------------------------------------------------------#

# Transformar composición a coordenadas CLR
compo_clr <- clr(compo)

# Media composicional global, media de Aitchison
media_global <- mean(compo)

# Pasar CLR a formato largo para ggplot
clr_df <- as.data.frame(compo_clr)
clr_df$cond <- cond

clr_df_largo <- clr_df %>%
  pivot_longer(
    cols = -cond,
    names_to = "ptm",
    values_to = "clr"
  )

# Boxplot de coordenadas CLR por proteoforma
box_plot_clr <- clr_df_largo %>%
  ggplot(aes(x = ptm, y = clr)) +
  geom_boxplot(
    fill = "lightblue",
    outlier.shape = NA,
    alpha = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  geom_point(
    aes(colour = cond),
    position = position_jitter(width = 0.1, height = 0),
    size = 2
  ) +
  labs(
    x = "Proteoforma",
    y = "CLR",
    colour = "Condición"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#----------------------------------------------------------#
# 3. Matriz de variación composicional global               #
#----------------------------------------------------------#

# La matriz de variación contiene varianzas de log-cocientes:
# T_ij = var(log(x_i / x_j))
var_matrix_global <- variation(compo)

mapa_var_global <- pheatmap(
  var_matrix_global,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Matriz de variación - conjunto de datos",
  silent = TRUE
)

#----------------------------------------------------------#
# 4. Medias composicionales por condición                   #
#----------------------------------------------------------#

# Medias composicionales por condición. Se evita tapply() sobre objetos acomp
# multivariantes y se calculan las medias sobre subconjuntos de filas.
compo_medias_cond <- lapply(
  split(seq_len(nrow(df.ancho)), cond),
  function(ii) {
    mean(acomp(df.ancho[ii, ptm_vars]))
  }
)

medias_cond <- tibble(
  ptm = ptm_names,
  media_asinc = as.numeric(compo_medias_cond$asinc),
  media_mitot = as.numeric(compo_medias_cond$mitot),
  media_global = as.numeric(media_global)
)

#----------------------------------------------------------#
# 5. Matrices de variación composicional por condición      #
#----------------------------------------------------------#

compo_asinc <- acomp(df.ancho[df.ancho$Estado == "asinc", ptm_vars])
compo_mitot <- acomp(df.ancho[df.ancho$Estado == "mitot", ptm_vars])

var_matrix_asinc <- variation(compo_asinc)
var_matrix_mitot <- variation(compo_mitot)

mapa_var_asinc <- pheatmap(
  var_matrix_asinc,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Matriz de variación - estado asíncrono",
  silent = TRUE
)

mapa_var_mitot <- pheatmap(
  var_matrix_mitot,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Matriz de variación - estado mitótico",
  silent = TRUE
)

#----------------------------------------------------------#
# 6. Agregación ternaria: noS10, S10_only, S10_K14          #
#----------------------------------------------------------#
#
# noS10:     proteoformas sin fosforilación en S10.
# S10_only: proteoformas con S10ph y sin K14ac.
#            Puede incluir K9acS10ph; por tanto, "only" se refiere a ausencia
#            de K14ac, no a ausencia total de acetilación.
# S10_K14:   proteoformas con S10ph y K14ac.

columnas_necesarias_ternario <- c(
  "unmod", "K9me1", "K9me2", "K9me3", "K9ac", "K14ac",
  "K9me1K14ac", "K9me2K14ac", "K9me3K14ac", "K9acK14ac",
  "S10ph", "K9me1S10ph", "K9me2S10ph", "K9me3S10ph", "K9acS10ph",
  "S10phK14ac", "K9me1S10phK14ac", "K9me2S10phK14ac",
  "K9me3S10phK14ac", "K9acS10phK14ac"
)

stopifnot(all(columnas_necesarias_ternario %in% names(df.ancho)))

df_grupos <- df.ancho %>%
  transmute(
    Estado = Estado,
    Dia = Dia,
    noS10 = unmod + K9me1 + K9me2 + K9me3 + K9ac + K14ac +
      K9me1K14ac + K9me2K14ac + K9me3K14ac + K9acK14ac,
    S10_only = S10ph + K9me1S10ph + K9me2S10ph + K9me3S10ph + K9acS10ph,
    S10_K14 = S10phK14ac + K9me1S10phK14ac + K9me2S10phK14ac +
      K9me3S10phK14ac + K9acS10phK14ac
  )

Y_grupos <- acomp(df_grupos[, c("noS10", "S10_only", "S10_K14")])

df_ternario <- bind_cols(
  df_grupos[, c("Estado", "Dia")],
  as.data.frame(Y_grupos)
)

g_ternario <- ggtern(
  df_ternario,
  aes(x = noS10, z = S10_only, y = S10_K14, colour = Estado)
) +
  geom_point(size = 3) +
  labs(
    title = "Agrupación ternaria de proteoformas H3 9-17",
    x = "Sin S10ph",
    y = "S10ph + K14ac",
    z = "S10ph sin K14ac",
    colour = "Condición"
  ) +
  theme_bw()

#----------------------------------------------------------#
# 7. PCA composicional en coordenadas ILR                  #
#----------------------------------------------------------#

Y_comp <- acomp(df.ancho[, ptm_vars])
Y_ilr <- ilr(Y_comp)

# El PCA se calcula en coordenadas ILR. No escalamos las coordenadas porque
# ILR ya proporciona coordenadas ortonormales en el espacio de Aitchison.
res_pca_ilr <- PCA(
  as.data.frame(Y_ilr),
  scale.unit = FALSE,
  graph = FALSE
)

pca_ilr_eig <- as.data.frame(res_pca_ilr$eig)
pca_ilr_eig <- pca_ilr_eig %>%
  rownames_to_column("componente")

# Coordenadas de individuos en PC1 y PC2
pca_ind_df <- as.data.frame(res_pca_ilr$ind$coord[, 1:2]) %>%
  mutate(
    Estado = df.ancho$Estado,
    Dia = df.ancho$Dia,
    muestra = row_number()
  )

g_pca_ilr_ind <- ggplot(
  pca_ind_df,
  aes(x = Dim.1, y = Dim.2, colour = Estado, label = Dia)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(size = 3) +
  geom_text(
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA composicional en coordenadas ILR",
    subtitle = "Números: día experimental",
    x = paste0("PC1 (", round(res_pca_ilr$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res_pca_ilr$eig[2, 2], 1), "%)"),
    colour = "Condición"
  ) +
  theme_minimal()

#----------------------------------------------------------#
# 8. Cargas PCA retroproyectadas a coordenadas CLR          #
#----------------------------------------------------------#

# Loadings del PCA en coordenadas ILR
loadings_ilr <- as.matrix(res_pca_ilr$var$coord)

# Base ILR usada por compositions. Tiene dimensión K x (K - 1).
V_ilr <- ilrBase(D = ncol(Y_comp))

# Proyección de direcciones principales desde ILR a CLR.
# Cada fila corresponde a una proteoforma; cada columna a una dimensión del PCA.
loadings_clr <- V_ilr %*% loadings_ilr
rownames(loadings_clr) <- ptm_names
colnames(loadings_clr) <- colnames(loadings_ilr)

loadings_clr_df <- as.data.frame(loadings_clr[, 1:2, drop = FALSE]) %>%
  rownames_to_column("PTM") %>%
  rename(
    PC1 = Dim.1,
    PC2 = Dim.2
  )

loadings_clr_long <- loadings_clr_df %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "loading_clr"
  )

# Visualización como barras: contribución de cada proteoforma a PC1 y PC2
g_pca_clr_loadings_bar <- loadings_clr_long %>%
  ggplot(aes(x = reorder(PTM, loading_clr), y = loading_clr)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  labs(
    title = "Direcciones principales expresadas en coordenadas CLR",
    subtitle = "PCA calculado en ILR; cargas proyectadas al espacio CLR",
    x = "Proteoforma",
    y = "Carga CLR"
  ) +
  theme_minimal()

# Visualización como plano PC1-PC2: posición de cada PTM según sus cargas
g_pca_clr_loadings_scatter <- ggplot(
  loadings_clr_df,
  aes(x = PC1, y = PC2, label = PTM)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point() +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    title = "Cargas de proteoformas en PC1 y PC2",
    subtitle = "PCA calculado en ILR; cargas retroproyectadas a CLR",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

# Gráfico principal recomendado para el informe
g_pca_clr_loadings <- g_pca_clr_loadings_scatter

#----------------------------------------------------------#
# 9. Alias temporales para compatibilidad con el QMD actual #
#----------------------------------------------------------#
# Estos nombres permiten que partes del QMD antiguo sigan funcionando mientras
# se actualizan los chunks. En la versión final del QMD conviene usar solo los
# nombres nuevos con guion bajo.

medias.con <- medias_cond
box.plot.clr <- box_plot_clr
mapa.var <- mapa_var_global
mapa.asinc <- mapa_var_asinc
mapa.mitot <- mapa_var_mitot
res_pca <- res_pca_ilr
ptm_loadings <- loadings_clr
df_loadings <- loadings_clr_df

#----------------------------------------------------------#
# 10. Guardar resultados                                   #
#----------------------------------------------------------#

resultados_descri_h3917 <- list(
  # Datos básicos
  df_ancho = df.ancho,
  ptm_vars = ptm_vars,
  ptm_names = ptm_names,
  cond = cond,

  # Objetos composicionales
  compo = compo,
  compo_clr = compo_clr,
  clr_df = clr_df,
  clr_df_largo = clr_df_largo,
  media_global = media_global,
  compo_medias_cond = compo_medias_cond,
  medias_cond = medias_cond,

  # Matrices de variación
  var_matrix_global = var_matrix_global,
  var_matrix_asinc = var_matrix_asinc,
  var_matrix_mitot = var_matrix_mitot,

  # Figuras descriptivas
  box_plot_clr = box_plot_clr,
  mapa_var_global = mapa_var_global,
  mapa_var_asinc = mapa_var_asinc,
  mapa_var_mitot = mapa_var_mitot,

  # Ternario
  df_grupos = df_grupos,
  Y_grupos = Y_grupos,
  df_ternario = df_ternario,
  g_ternario = g_ternario,

  # PCA-ILR
  Y_comp = Y_comp,
  Y_ilr = Y_ilr,
  res_pca_ilr = res_pca_ilr,
  pca_ilr_eig = pca_ilr_eig,
  pca_ind_df = pca_ind_df,
  g_pca_ilr_ind = g_pca_ilr_ind,

  # Cargas PCA expresadas en CLR
  loadings_ilr = loadings_ilr,
  V_ilr = V_ilr,
  loadings_clr = loadings_clr,
  loadings_clr_df = loadings_clr_df,
  loadings_clr_long = loadings_clr_long,
  g_pca_clr_loadings_bar = g_pca_clr_loadings_bar,
  g_pca_clr_loadings_scatter = g_pca_clr_loadings_scatter,
  g_pca_clr_loadings = g_pca_clr_loadings,

  # Alias temporales
  medias.con = medias.con,
  box.plot.clr = box.plot.clr,
  mapa.var = mapa.var,
  mapa.asinc = mapa.asinc,
  mapa.mitot = mapa.mitot,
  res_pca = res_pca,
  ptm_loadings = ptm_loadings,
  df_loadings = df_loadings
)

saveRDS(
  resultados_descri_h3917,
  here::here("DatosProcesados", "31_resultados_descri_h3917.rds")
)

message(
  "Resultados guardados en: ",
  here::here("DatosProcesados", "31_resultados_descri_h3917.rds")
)
