#!/usr/bin/env Rscript

# 78_precalcular_resultados_proporciones_epistasia.R
# ------------------------------------------------------------
# Precalcula resúmenes posteriores para el subcapítulo de resultados
# del modelo Dirichlet epistático con kappa común y efecto aleatorio de día.
#
# Entradas esperadas, en DatosProcesados/:
#   - df_h3917.Rdata: contiene df.ancho o df_ancho
#   - stan_epistasia_fits.Rdata: contiene fit_kappa_comun_dia
#
# Salidas, en DatosProcesados/:
#   - 78_resultados_proporciones_epistasia.rds
#   - 78_tabla_beta_logratio_epistasia.csv
#   - 78_tabla_Q_cond_epistasia.csv
#   - 78_tabla_balances_acetilacion_epistasia.csv
#   - 78_draws_Q_long_epistasia.rds
#   - 78_draws_balances_acetilacion_epistasia.rds
#
# Balances ILR incluidos (versión z1-z4):
#   z1_S10ph_vs_NoS10ph:
#      Balance de fosforilación global. Proteoformas con S10ph frente
#      a proteoformas sin S10ph.
#   z2_Ac_vs_NoAc_given_NoS10ph:
#      Balance de acetilación en ausencia de fosforilación. Dentro del
#      subconjunto S10ph-, compara proteoformas acetiladas frente a
#      no acetiladas.
#   z3_Ac_vs_NoAc_given_S10ph:
#      Balance de acetilación en presencia de fosforilación. Dentro del
#      subconjunto S10ph+, compara proteoformas acetiladas frente a
#      no acetiladas.
#   z4_MonoAc_vs_DiAc:
#      Balance de grado de acetilación. Compara variantes monoacetiladas
#      frente a doblemente acetiladas.
#
# Nota: además de las tablas resumen, este script guarda los draws de Q_cond
# como array (iter x K x condición) en el objeto final, para evitar recalcular
# balances costosos dentro del documento Quarto.
#
# Este script está pensado para ejecutarse antes de renderizar el book:
#   source(here::here("ScriptsR", "78_precalcular_resultados_proporciones_epistasia.R"))

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(readr)
  library(rstan)
})

# ---------- Configuración ----------

dir_proc <- here::here("DatosProcesados")
file_datos <- file.path(dir_proc, "df_h3917.Rdata")
file_fits  <- file.path(dir_proc, "stan_epistasia_fits.Rdata")

out_rds          <- file.path(dir_proc, "78_resultados_proporciones_epistasia.rds")
out_beta_csv     <- file.path(dir_proc, "78_tabla_beta_logratio_epistasia.csv")
out_Q_csv        <- file.path(dir_proc, "78_tabla_Q_cond_epistasia.csv")
out_bal_csv      <- file.path(dir_proc, "78_tabla_balances_acetilacion_epistasia.csv")
out_Q_long_rds   <- file.path(dir_proc, "78_draws_Q_long_epistasia.rds")
out_bal_long_rds <- file.path(dir_proc, "78_draws_balances_acetilacion_epistasia.rds")

if (!file.exists(file_datos)) {
  stop("No encuentro el archivo de datos: ", file_datos)
}
if (!file.exists(file_fits)) {
  stop("No encuentro el archivo de ajustes: ", file_fits)
}

# ---------- Carga de objetos ----------

load(file_datos)
load(file_fits)

# El objeto de datos aparece en tus scripts unas veces como df.ancho y otras como df_ancho.
if (exists("df_ancho")) {
  df_h3917 <- df_ancho
} else if (exists("df.ancho")) {
  df_h3917 <- df.ancho
} else {
  stop("No encuentro `df_ancho` ni `df.ancho` dentro de ", file_datos)
}

if (!exists("fit_kappa_comun_dia")) {
  stop("No encuentro `fit_kappa_comun_dia` dentro de ", file_fits)
}

ptm_names <- colnames(df_h3917)[3:22]
K <- length(ptm_names)

# ---------- Funciones auxiliares ----------

summ_draws <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble(
    media   = mean(x, na.rm = TRUE),
    mediana = median(x, na.rm = TRUE),
    q025    = unname(quantile(x, probs[1], na.rm = TRUE)),
    q975    = unname(quantile(x, probs[3], na.rm = TRUE)),
    p_gt_0  = mean(x > 0, na.rm = TRUE),
    p_lt_0  = mean(x < 0, na.rm = TRUE)
  )
}

softmax_minus_phi <- function(phi_mat) {
  # phi_mat: draws x K x C
  # Devuelve Q = softmax(-phi) con la misma dimensión.
  dims <- dim(phi_mat)
  if (length(dims) != 3L) stop("`phi_mat` debe tener dimensión iter x K x C")
  out <- array(NA_real_, dim = dims)
  for (cc in seq_len(dims[3])) {
    eta <- -phi_mat[, , cc, drop = FALSE][, , 1]
    eta <- eta - apply(eta, 1, max)
    exp_eta <- exp(eta)
    out[, , cc] <- exp_eta / rowSums(exp_eta)
  }
  out
}

balance_two_groups <- function(x, idx_g1, idx_g2) {
  # Balance ILR entre dos grupos de partes:
  # sqrt(r*s/(r+s)) * (mean(log x_g1) - mean(log x_g2))
  r <- length(idx_g1)
  s <- length(idx_g2)
  if (r == 0L || s == 0L) return(NA_real_)
  sqrt(r * s / (r + s)) * (mean(log(x[idx_g1])) - mean(log(x[idx_g2])))
}

# ---------- Extracción posterior ----------

pars_fit <- names(rstan::extract(fit_kappa_comun_dia, permuted = TRUE, inc_warmup = FALSE))

if (!("phi" %in% pars_fit)) {
  stop("El ajuste no contiene el parámetro/generated quantity `phi`.")
}

phi_draws <- rstan::extract(fit_kappa_comun_dia, pars = "phi", permuted = TRUE)$phi
# Esperado: iter x K x 2, condición 1 = asincronía, condición 2 = mitosis.
if (length(dim(phi_draws)) != 3L) {
  stop("`phi` no tiene dimensión iter x K x condición. Dimensiones encontradas: ",
       paste(dim(phi_draws), collapse = " x "))
}
if (dim(phi_draws)[2] != K) {
  stop("El número de PTMs en `phi` (", dim(phi_draws)[2],
       ") no coincide con df_h3917[,3:22] (", K, ").")
}
if (dim(phi_draws)[3] < 2L) {
  stop("`phi` debe contener al menos dos condiciones: asincronía y mitosis.")
}

n_draws <- dim(phi_draws)[1]

if ("Q_cond" %in% pars_fit) {
  Q_draws <- rstan::extract(fit_kappa_comun_dia, pars = "Q_cond", permuted = TRUE)$Q_cond
} else {
  message("No encuentro `Q_cond`; calculo Q_cond = softmax(-phi) a partir de `phi`.")
  Q_draws <- softmax_minus_phi(phi_draws)
}

if (length(dim(Q_draws)) != 3L) {
  stop("`Q_cond` no tiene dimensión iter x K x condición. Dimensiones encontradas: ",
       paste(dim(Q_draws), collapse = " x "))
}
if (dim(Q_draws)[2] != K) {
  stop("El número de PTMs en `Q_cond` (", dim(Q_draws)[2],
       ") no coincide con df_h3917[,3:22] (", K, ").")
}

# ---------- 1) Efecto en campo energético beta_k y log(Q_mit/Q_asinc) ----------

beta_tbl <- map_dfr(seq_len(K), function(k) {
  beta_k <- phi_draws[, k, 2] - phi_draws[, k, 1]
  log_ratio_k <- log(Q_draws[, k, 2] / Q_draws[, k, 1])

  bind_cols(
    tibble(PTM = ptm_names[k]),
    summ_draws(beta_k) %>% rename_with(~ paste0("beta_", .x)),
    summ_draws(log_ratio_k) %>% rename_with(~ paste0("log_ratio_Q_", .x))
  ) %>%
    mutate(
      direccion_beta = case_when(
        beta_q025 > 0 ~ "beta > 0: aumenta energía en mitosis",
        beta_q975 < 0 ~ "beta < 0: disminuye energía en mitosis",
        TRUE ~ "IC95 incluye 0"
      ),
      direccion_Q = case_when(
        log_ratio_Q_q025 > 0 ~ "Q aumenta en mitosis",
        log_ratio_Q_q975 < 0 ~ "Q disminuye en mitosis",
        TRUE ~ "IC95 incluye 0"
      )
    )
})

# ---------- 2) Tabla de proporciones medias posteriores Q por condición ----------

Q_cond_tbl <- map_dfr(seq_len(K), function(k) {
  q_async <- Q_draws[, k, 1]
  q_mito  <- Q_draws[, k, 2]
  diff_abs <- q_mito - q_async
  ratio_rel <- q_mito / q_async
  log_ratio <- log(ratio_rel)

  bind_cols(
    tibble(PTM = ptm_names[k]),
    summ_draws(q_async) %>% select(-p_gt_0, -p_lt_0) %>% rename_with(~ paste0("Q_asinc_", .x)),
    summ_draws(q_mito)  %>% select(-p_gt_0, -p_lt_0) %>% rename_with(~ paste0("Q_mitosis_", .x)),
    summ_draws(diff_abs) %>% rename_with(~ paste0("diff_abs_", .x)),
    summ_draws(ratio_rel) %>% select(-p_gt_0, -p_lt_0) %>% rename_with(~ paste0("ratio_", .x)),
    summ_draws(log_ratio) %>% rename_with(~ paste0("log_ratio_", .x))
  ) %>%
    mutate(
      direccion_Q = case_when(
        diff_abs_q025 > 0 ~ "aumenta en mitosis",
        diff_abs_q975 < 0 ~ "disminuye en mitosis",
        TRUE ~ "IC95 incluye 0"
      )
    )
})

# Draws largos para gráficos de densidades por PTM.
Q_long <- map_dfr(seq_len(K), function(k) {
  tibble(
    draw = seq_len(n_draws),
    PTM = ptm_names[k],
    asincronia = Q_draws[, k, 1],
    mitosis = Q_draws[, k, 2]
  )
}) %>%
  pivot_longer(
    cols = c("asincronia", "mitosis"),
    names_to = "condicion",
    values_to = "Q"
  ) %>%
  mutate(condicion = factor(condicion, levels = c("asincronia", "mitosis")))

# ---------- 3) Balances ILR z1-z4 desde Q_cond ----------

# Clasificación de proteoformas.
# Usamos búsqueda literal en los nombres extraídos de df_h3917[, 3:22].
idx_s10ph    <- which(str_detect(ptm_names, fixed("S10ph")))
idx_no_s10ph <- setdiff(seq_len(K), idx_s10ph)

idx_ac    <- which(str_detect(ptm_names, fixed("ac")))
idx_no_ac <- setdiff(seq_len(K), idx_ac)

# Acetilación en ausencia/presencia de fosforilación.
idx_ac_no_s10ph    <- intersect(idx_ac, idx_no_s10ph)
idx_noac_no_s10ph <- intersect(idx_no_ac, idx_no_s10ph)
idx_ac_s10ph       <- intersect(idx_ac, idx_s10ph)
idx_noac_s10ph     <- intersect(idx_no_ac, idx_s10ph)

# Grado de acetilación: monoacetiladas frente a doblemente acetiladas.
tiene_K9ac  <- str_detect(ptm_names, fixed("K9ac"))
tiene_K14ac <- str_detect(ptm_names, fixed("K14ac"))
n_acetilos <- as.integer(tiene_K9ac) + as.integer(tiene_K14ac)
idx_monoac <- which(n_acetilos == 1L)
idx_diac   <- which(n_acetilos == 2L)

# Comprobaciones explícitas de todos los grupos.
check_group <- function(idx, nombre) {
  if (length(idx) == 0L) {
    stop("El grupo `", nombre, "` está vacío. Revisa los nombres de PTMs: ",
         paste(ptm_names, collapse = ", "))
  }
}

check_group(idx_s10ph, "S10ph")
check_group(idx_no_s10ph, "No-S10ph")
check_group(idx_ac_no_s10ph, "Ac & No-S10ph")
check_group(idx_noac_no_s10ph, "NoAc & No-S10ph")
check_group(idx_ac_s10ph, "Ac & S10ph")
check_group(idx_noac_s10ph, "NoAc & S10ph")
check_group(idx_monoac, "monoacetiladas")
check_group(idx_diac, "diacetiladas")

# Definición final de los cuatro balances biológicos.
# Convención: balance = grupo_positivo / grupo_negativo.
# Valores mayores implican mayor peso relativo del grupo_positivo.
grupos_balances <- list(
  z1_S10ph_vs_NoS10ph = list(
    etiqueta = "z1: S10ph / No-S10ph",
    grupo_positivo = idx_s10ph,
    grupo_negativo = idx_no_s10ph,
    interpretacion = paste(
      "Balance de la fosforilación global:",
      "proteoformas con S10ph frente a proteoformas sin S10ph."
    )
  ),
  z2_Ac_vs_NoAc_given_NoS10ph = list(
    etiqueta = "z2: Ac / No-Ac | S10ph-",
    grupo_positivo = idx_ac_no_s10ph,
    grupo_negativo = idx_noac_no_s10ph,
    interpretacion = paste(
      "Balance de la acetilación en ausencia de fosforilación:",
      "dentro del sustrato S10ph-, compara formas acetiladas frente a formas no acetiladas."
    )
  ),
  z3_Ac_vs_NoAc_given_S10ph = list(
    etiqueta = "z3: Ac / No-Ac | S10ph+",
    grupo_positivo = idx_ac_s10ph,
    grupo_negativo = idx_noac_s10ph,
    interpretacion = paste(
      "Balance de la acetilación en presencia de fosforilación:",
      "dentro del sustrato S10ph+, compara formas acetiladas frente a formas no acetiladas."
    )
  ),
  z4_MonoAc_vs_DiAc = list(
    etiqueta = "z4: monoacetiladas / diacetiladas",
    grupo_positivo = idx_monoac,
    grupo_negativo = idx_diac,
    interpretacion = paste(
      "Balance de grado de acetilación:",
      "compara variantes monoacetiladas frente a variantes doblemente acetiladas."
    )
  )
)

# Tabla diagnóstica de grupos. Se guarda en metadata y permite verificar
# qué proteoformas entran en el numerador y el denominador de cada balance.
tabla_grupos_balances <- imap_dfr(grupos_balances, function(g, nombre_balance) {
  tibble(
    balance = nombre_balance,
    etiqueta = g$etiqueta,
    grupo = c("Numerador", "Denominador"),
    n = c(length(g$grupo_positivo), length(g$grupo_negativo)),
    proteoformas = c(
      paste(ptm_names[g$grupo_positivo], collapse = ", "),
      paste(ptm_names[g$grupo_negativo], collapse = ", ")
    )
  )
})

# Cálculo vectorizado de un balance para todos los draws de una condición.
balance_two_groups_matrix <- function(Q_mat, idx_g1, idx_g2) {
  # Q_mat: draws x K
  r <- length(idx_g1)
  s <- length(idx_g2)
  if (r == 0L || s == 0L) return(rep(NA_real_, nrow(Q_mat)))
  sqrt(r * s / (r + s)) * (
    rowMeans(log(Q_mat[, idx_g1, drop = FALSE])) -
      rowMeans(log(Q_mat[, idx_g2, drop = FALSE]))
  )
}

balances_long <- imap_dfr(grupos_balances, function(g, nombre_balance) {
  bal_async <- balance_two_groups_matrix(
    Q_mat = Q_draws[, , 1, drop = FALSE][, , 1],
    idx_g1 = g$grupo_positivo,
    idx_g2 = g$grupo_negativo
  )
  bal_mito <- balance_two_groups_matrix(
    Q_mat = Q_draws[, , 2, drop = FALSE][, , 1],
    idx_g1 = g$grupo_positivo,
    idx_g2 = g$grupo_negativo
  )

  tibble(
    draw = seq_len(n_draws),
    balance = nombre_balance,
    etiqueta = g$etiqueta,
    interpretacion = g$interpretacion,
    asincronia = bal_async,
    mitosis = bal_mito,
    diferencia_mitosis_menos_asincronia = bal_mito - bal_async
  ) %>%
    pivot_longer(
      cols = c("asincronia", "mitosis"),
      names_to = "condicion",
      values_to = "valor"
    ) %>%
    mutate(condicion = factor(condicion, levels = c("asincronia", "mitosis")))
})

balances_tbl_cond <- balances_long %>%
  group_by(balance, etiqueta, interpretacion, condicion) %>%
  summarise(
    media = mean(valor, na.rm = TRUE),
    mediana = median(valor, na.rm = TRUE),
    q025 = unname(quantile(valor, 0.025, na.rm = TRUE)),
    q975 = unname(quantile(valor, 0.975, na.rm = TRUE)),
    .groups = "drop"
  )

balances_tbl_diff <- balances_long %>%
  distinct(draw, balance, etiqueta, interpretacion, diferencia_mitosis_menos_asincronia) %>%
  group_by(balance, etiqueta, interpretacion) %>%
  summarise(
    diff_media = mean(diferencia_mitosis_menos_asincronia, na.rm = TRUE),
    diff_mediana = median(diferencia_mitosis_menos_asincronia, na.rm = TRUE),
    diff_q025 = unname(quantile(diferencia_mitosis_menos_asincronia, 0.025, na.rm = TRUE)),
    diff_q975 = unname(quantile(diferencia_mitosis_menos_asincronia, 0.975, na.rm = TRUE)),
    diff_p_gt_0 = mean(diferencia_mitosis_menos_asincronia > 0, na.rm = TRUE),
    diff_p_lt_0 = mean(diferencia_mitosis_menos_asincronia < 0, na.rm = TRUE),
    .groups = "drop"
  )

balances_tbl <- balances_tbl_cond %>%
  select(balance, etiqueta, interpretacion, condicion, media, mediana, q025, q975) %>%
  pivot_wider(
    names_from = condicion,
    values_from = c(media, mediana, q025, q975),
    names_glue = "{.value}_{condicion}"
  ) %>%
  left_join(balances_tbl_diff, by = c("balance", "etiqueta", "interpretacion"))

# ---------- Guardado ----------

resultados <- list(
  metadata = list(
    fecha = Sys.time(),
    modelo = "fit_kappa_comun_dia",
    archivo_datos = file_datos,
    archivo_fits = file_fits,
    n_draws = n_draws,
    K = K,
    ptm_names = ptm_names,
    grupos_balances = list(
      S10ph = ptm_names[idx_s10ph],
      no_S10ph = ptm_names[idx_no_s10ph],
      ac_no_S10ph = ptm_names[idx_ac_no_s10ph],
      no_ac_no_S10ph = ptm_names[idx_noac_no_s10ph],
      ac_S10ph = ptm_names[idx_ac_s10ph],
      no_ac_S10ph = ptm_names[idx_noac_s10ph],
      monoacetiladas = ptm_names[idx_monoac],
      diacetiladas = ptm_names[idx_diac]
    ),
    tabla_grupos_balances = tabla_grupos_balances
  ),
  beta_logratio = beta_tbl,
  Q_cond = Q_cond_tbl,
  Q_draws_array = Q_draws,
  Q_long = Q_long,
  balances = balances_tbl,
  balances_long = balances_long
)

saveRDS(resultados, out_rds)
write_csv(beta_tbl, out_beta_csv)
write_csv(Q_cond_tbl, out_Q_csv)
write_csv(balances_tbl, out_bal_csv)
saveRDS(Q_long, out_Q_long_rds)
saveRDS(balances_long, out_bal_long_rds)

message("Resultados guardados en:")
message("  ", out_rds)
message("  ", out_beta_csv)
message("  ", out_Q_csv)
message("  ", out_bal_csv)
message("  ", out_Q_long_rds)
message("  ", out_bal_long_rds)
