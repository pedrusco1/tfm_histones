##############################################################
#     CONSTRUIR MATRIZ B PARA EL FRAGMENTO H3_9_17           #
##############################################################

### Cargar las librerias
library(tidyverse)
library(Matrix)

### Cargar datos
load(here::here("DatosProcesados", "df_h3917.Rdata"))

### Cargar funcion para matriz B
source(here::here("Funciones", "fun_hacer_matriz_B_OK.R"))

# 5. CALCULAR DELTA (Matriz Laplaciana de Epistasia)
# Delta = B^T * B
B_list <-  construir_matriz_B(
  nombres_estados = names(df.ancho[,3:22]),
  niveles = list(
    K9 = c("un", "me1", "me2", "me3", "ac" ),
    S10 = c("un", "ph"),
    K14 = c("un", "ac")
  )
  )

B <- B_list$B
B_label <- B_list$info_B$Etiqueta

Matriz_Delta <- t(B) %*% B

# 6. Verificación Final. Esto es muy importante para interpretar los resultados
# El orden de Delta DEBE ser el mismo que el de tus columnas de datos
#all(colnames(Delta) == names(df.ancho[ ,3:22])) # Devuelve TRUE, todo ok de momento


# Matriz_Delta <- Delta

# Calculamos el rango de la matriz Delta (necesario para la Ec. 9)
# Generalmente rank = K - (número de componentes conectadas)
# Obtener los autovalores (eigenvalues) de la matriz Delta
# Definir una tolerancia (threshold)
# Los valores menores a esto se consideran 0 debido a la precisión de la máquina
ev <- eigen(Matriz_Delta, only.values = TRUE)$values
tol <- max(1e-12 * ev[1], 1e-12)
rango_Delta <- qr(Matriz_Delta, tol = 1e-7)$rank 
s <- length(rownames(B)) # 's: numero ecuaciones epistasia local' 

delta_obj <- list(
  Delta = Matriz_Delta,
  rank_Delta = rango_Delta,
  B = B,
  B_label = B_label
  )
save(delta_obj, file = file.path(data_dir, "delta_obj.Rdata"))
