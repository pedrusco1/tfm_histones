construir_matriz_B <- function(nombres_estados, niveles) {
  
  # 1. Función interna: parsear_estados (Se mantiene igual)
  parsear_estados <- function(nombres_estados, niveles) {
    posiciones <- names(niveles)
    df <- data.frame(nombre_estado = nombres_estados, stringsAsFactors = FALSE)
    for (pos in posiciones) {
      niveles_pos <- niveles[[pos]]
      patron <- paste0(pos, "(", paste(niveles_pos, collapse = "|"), ")")
      extraido <- stringr::str_extract(nombres_estados, patron)
      nivel <- stringr::str_remove(extraido, pos)
      nivel[is.na(nivel)] <- "un"
      df[[pos]] <- nivel
    }
    df
  }
  
  # 2. Parsear estados
  nodos <- parsear_estados(nombres_estados, niveles)
  posiciones <- names(niveles)
  n_pos <- length(posiciones)
  pares_pos <- combn(posiciones, 2, simplify = FALSE)
  
  # 3. Calcular número total de filas de B (Ajustado para n >= 2)
  total_filas <- 0
  for (par in pares_pos) {
    p1 <- par[1]; p2 <- par[2]
    p3_names <- setdiff(posiciones, par)
    
    # Si hay 3 o más factores, multiplicamos por los niveles de los "otros"
    # Si solo hay 2 factores, el multiplicador es 1
    multiplicador <- if(length(p3_names) > 0) length(niveles[[p3_names[1]]]) else 1
    
    total_filas <- total_filas +
      multiplicador *
      choose(length(niveles[[p1]]), 2) *
      choose(length(niveles[[p2]]), 2)
  }
  
  # 4. Inicializar matriz B
  B <- matrix(0, nrow = total_filas, ncol = length(nombres_estados))
  colnames(B) <- nombres_estados
  nombres_epistasia <- character(total_filas)
  fila <- 1
  
  # 5. Construir matriz B con lógica condicional
  for (par in pares_pos) {
    p1 <- par[1]; p2 <- par[2]
    p3_names <- setdiff(posiciones, par)
    
    comb_p1 <- combn(niveles[[p1]], 2, simplify = FALSE)
    comb_p2 <- combn(niveles[[p2]], 2, simplify = FALSE)
    
    # Definir los niveles del tercer factor (si no hay, usamos un placeholder)
    niveles_p3 <- if(length(p3_names) > 0) niveles[[p3_names[1]]] else "NO_P3"
    
    for (s3 in niveles_p3) {
      for (c1 in comb_p1) {
        for (c2 in comb_p2) {
          
          # Lógica de filtrado dinámico
          filtro <- nodos[[p1]] == c1[2] & nodos[[p2]] == c2[2]
          filtro10 <- nodos[[p1]] == c1[2] & nodos[[p2]] == c2[1]
          filtro01 <- nodos[[p1]] == c1[1] & nodos[[p2]] == c2[2]
          filtro00 <- nodos[[p1]] == c1[1] & nodos[[p2]] == c2[1]
          
          # Si existe p3, añadimos esa restricción al filtro
          if(s3 != "NO_P3") {
            p3 <- p3_names[1]
            filtro <- filtro & nodos[[p3]] == s3
            filtro10 <- filtro10 & nodos[[p3]] == s3
            filtro01 <- filtro01 & nodos[[p3]] == s3
            filtro00 <- filtro00 & nodos[[p3]] == s3
          }
          
          idx11 <- which(filtro)
          idx10 <- which(filtro10)
          idx01 <- which(filtro01)
          idx00 <- which(filtro00)
          
          B[fila, idx11] <-  1
          B[fila, idx10] <- -1
          B[fila, idx01] <- -1
          B[fila, idx00] <-  1
          
          # Nombre de la fila
          suffix <- if(s3 != "NO_P3") paste0("|", p3_names[1], s3) else ""
          nombres_epistasia[fila] <- paste0(p1, c1[2], ":", p2, c2[2], suffix)
          
          fila <- fila + 1
        }
      }
    }
  }
  
  rownames(B) <- nombres_epistasia
  return(B)
}
