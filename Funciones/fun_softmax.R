###################
# Funcion soxtmax #
###################

softmax <- function(z) {
  # Restamos el máximo para estabilidad numérica
  e_z <- exp(z - max(z))
  return(e_z / sum(e_z))
}
