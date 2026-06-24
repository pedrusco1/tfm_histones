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
 
Delta <- t(B) %*% B
nombres_estados = names(df.ancho[,3:22])

# 6. Verificación Final. Esto es muy importante para interpretar los resultados
# El orden de Delta DEBE ser el mismo que el de tus columnas de datos
all(colnames(Delta) == nombres_estados) # Devuelve TRUE, todo ok de momento


Matriz_Delta <- Delta

# Calculamos el rango de la matriz Delta (necesario para la Ec. 9)
# Generalmente rank = K - (número de componentes conectadas)
rango_Delta <- qr(Matriz_Delta, tol = 1e-7)$rank 
num_coeficientes <- length(rownames(B)) # 's: numero cuefecientes epistasia local' 

# 1. Obtener los autovalores (eigenvalues) de la matriz Delta
ev <- eigen(Matriz_Delta, only.values = TRUE)$values

# 2. Definir una tolerancia (threshold)
# Los valores menores a esto se consideran 0 debido a la precisión de la máquina
tol <- max(1e-12 * ev[1], 1e-12) 

# 3. Contar autovalores positivos (Este es el rango real según el documento)
rango_Delta <- sum(ev > tol)

# 4. Cálculo del s-pseudo-determinante (opcional, pero mencionado en Eq. 8)
# log_det_pseudo <- sum(log(ev[ev > tol]))

print(paste("El rango calculado es:", rango_Delta)) # coicide con el calculado por descomposición QR

