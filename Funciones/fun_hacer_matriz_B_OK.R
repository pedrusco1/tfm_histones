# Construye la funcion de coeficientes de epistasia local
# saca las etiquetas de los cuadrados con los cambios en la 
# posicion inferior y en la superior con sus contextos

construir_matriz_B <- function(nombres_estados, niveles) {
  # niveles: lista con elementos K9, S10, K14, etc.
  # nombres_estados: vector con los nombres de las proteoformas (columnas de df.ancho)
  
  # 1) Parsear estados a niveles por locus -------------------------------------
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
  
  nodos      <- parsear_estados(nombres_estados, niveles)
  posiciones <- names(niveles)
  pares_pos  <- combn(posiciones, 2, simplify = FALSE)
  
  # 2) Calcular número máximo de filas de B -----------------------------------
  total_filas <- 0
  for (par in pares_pos) {
    p1 <- par[1]; p2 <- par[2]
    p3_names <- setdiff(posiciones, par)
    
    multiplicador <- if (length(p3_names) > 0) {
      length(niveles[[p3_names[1]]])
    } else {
      1
    }
    
    total_filas <- total_filas +
      multiplicador *
      choose(length(niveles[[p1]]), 2) *
      choose(length(niveles[[p2]]), 2)
  }
  
  # 3) Inicializar B e info_B --------------------------------------------------
  B <- matrix(0, nrow = total_filas, ncol = length(nombres_estados))
  colnames(B) <- nombres_estados
  
  info_B <- data.frame(
    fila           = integer(total_filas),
    p1             = character(total_filas),
    p2             = character(total_filas),
    p3             = character(total_filas),
    nivel_p3       = character(total_filas),
    nivel_p1_low   = character(total_filas),
    nivel_p1_high  = character(total_filas),
    nivel_p2_low   = character(total_filas),
    nivel_p2_high  = character(total_filas),
    Cuadrado       = character(total_filas),
    Etiqueta       = character(total_filas),
    stringsAsFactors = FALSE
  )
  
  fila <- 1
  
  # 4) Construir B e info_B ----------------------------------------------------
  for (par in pares_pos) {
    p1 <- par[1]; p2 <- par[2]
    p3_names <- setdiff(posiciones, par)
    
    comb_p1 <- combn(niveles[[p1]], 2, simplify = FALSE)
    comb_p2 <- combn(niveles[[p2]], 2, simplify = FALSE)
    
    niveles_p3 <- if (length(p3_names) > 0) {
      niveles[[p3_names[1]]]
    } else {
      "NO_P3"
    }
    
    for (s3 in niveles_p3) {
      for (c1 in comb_p1) {
        for (c2 in comb_p2) {
          
          # filtros 11,10,01,00
          filtro   <- nodos[[p1]] == c1[2] & nodos[[p2]] == c2[2]
          filtro10 <- nodos[[p1]] == c1[2] & nodos[[p2]] == c2[1]
          filtro01 <- nodos[[p1]] == c1[1] & nodos[[p2]] == c2[2]
          filtro00 <- nodos[[p1]] == c1[1] & nodos[[p2]] == c2[1]
          
          if (s3 != "NO_P3") {
            p3 <- p3_names[1]
            filtro   <- filtro   & nodos[[p3]] == s3
            filtro10 <- filtro10 & nodos[[p3]] == s3
            filtro01 <- filtro01 & nodos[[p3]] == s3
            filtro00 <- filtro00 & nodos[[p3]] == s3
          } else {
            p3 <- NA_character_
          }
          
          idx11 <- which(filtro)
          idx10 <- which(filtro10)
          idx01 <- which(filtro01)
          idx00 <- which(filtro00)
          
          # Requerimos exactamente un estado para cada vértice
          if (length(idx11) != 1 ||
              length(idx10) != 1 ||
              length(idx01) != 1 ||
              length(idx00) != 1) {
            next
          }
          
          # Asignar pesos del cuadrado: 11 - 10 - 01 + 00
          B[fila, idx11] <-  1
          B[fila, idx10] <- -1
          B[fila, idx01] <- -1
          B[fila, idx00] <-  1
          
          # Construir nombres ricos (con niveles bajos y altos)
          nivel_p1_low  <- c1[1]
          nivel_p1_high <- c1[2]
          nivel_p2_low  <- c2[1]
          nivel_p2_high <- c2[2]
          
          # nombre "clásico" (solo niveles altos + fondo)
          suffix_fondo <- if (!is.na(p3) && s3 != "NO_P3") {
            paste0("|", p3, s3)
          } else {
            ""
          }
          nombre_clasico <- paste0(p1, nivel_p1_high, ":", p2, nivel_p2_high, suffix_fondo)
          
          # nombre detallado con niveles bajos y altos
          detalle <- paste0(
            p1, nivel_p1_low, "→", p1, nivel_p1_high, " : ",
            p2, nivel_p2_low, "→", p2, nivel_p2_high
          )
          if (!is.na(p3) && s3 != "NO_P3") {
            detalle <- paste0(detalle, " | ", p3, s3)
          }
          
          info_B$fila[fila]          <- fila
          info_B$p1[fila]            <- p1
          info_B$p2[fila]            <- p2
          info_B$p3[fila]            <- ifelse(is.na(p3), "", p3)
          info_B$nivel_p3[fila]      <- ifelse(is.na(p3), "", s3)
          info_B$nivel_p1_low[fila]  <- nivel_p1_low
          info_B$nivel_p1_high[fila] <- nivel_p1_high
          info_B$nivel_p2_low[fila]  <- nivel_p2_low
          info_B$nivel_p2_high[fila] <- nivel_p2_high
          info_B$Cuadrado[fila]      <- nombre_clasico
          info_B$Etiqueta[fila]      <- detalle
          
          fila <- fila + 1
        }
      }
    }
  }
  
  # 5) Recortar si se han saltado combinaciones -------------------------------
  if (fila <= total_filas) {
    B      <- B[1:(fila - 1), , drop = FALSE]
    info_B <- info_B[1:(fila - 1), , drop = FALSE]
  }
  
  rownames(B) <- info_B$Cuadrado
  
  # salida: lista con B e info_B
  return (list(B = B, info_B = info_B))
}

