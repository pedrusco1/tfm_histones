######################################################
# Posterior predictive check: modelo epistasa H3_9-17#
######################################################

library(bayesplot)
library(ggplot2)
library(rstan)
library(gridExtra)

# Cargar los datos y funciones
source(file.path("ScriptsR", "40_stan_config.R"))
source(file.path("ScriptsR", "41_cargar_datos_h3917.R"))
load(file.path(data_dir, "stan_epis_h3917.Rdata"))
load(file.path(data_dir, "df_h3917.Rdata"))
source(file.path("Funciones", "fun_ppc_varianza.R"))

# Extraer las matrices de réplicas (Y_rep) y observaciones (Y_mat)
y_rep_comun <- rstan::extract(fit_kappa_comun)$Y_rep
y_rep_cond <- rstan::extract(fit_kappa_cond)$Y_rep
y_obs <- Y_mat
idx_asin <- c(1,3,5)
idx_mito <- c(2,4,6)
samples <- sample(1:dim(y_rep_comun)[1], 100)
# Definir nombres de los estados de PTM para los gráficos 
colnames(y_obs) <- names(df.ancho[ ,3:22]) 

# -------------------------------------------------
# PPC para la media de unmmod
# -------------------------------------------------
# PPC para el modelo de Kappa Común
p1 <- ppc_dens_overlay(y_obs[idx_asin,1], y_rep_comun[samples, idx_asin, 1]) + 
  ggtitle("Modelo Kappa Común: Densidad PTM 1") +
  theme_bw()

# PPC para el modelo de Kappa Condicional
p2 <- ppc_dens_overlay(y_obs[idx_asin,1], y_rep_cond[samples, idx_asin , 1]) + 
  ggtitle("Modelo Kappa Condicional: Densidad PTM 1")+
  theme_bw()

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

# Comparación de la media de las proporciones
ppc_stat(y = y_obs[,1], yrep = y_rep_cond[, , 1], stat = "mean")

# Comparación de la dispersión (relevante para el parámetro kappa)
ppc_stat(y = y_obs[,1], yrep = y_rep_cond[, , 1], stat = "sd")
ppc_stat(y = y_obs[,1], yrep = y_rep_comun[, , 1], stat = "sd")
# Creamos un vector con las etiquetas de las condiciones
condiciones <- ifelse(df.ancho$Estado == "asinc", "Async", "Mitosis")

# Visualizar la distribución de una PTM específica agrupada por condición
ppc_means <- list()
# 2. Ejecutar el bucle y guardar cada gráfico
for(i in 1:20) {
  ppc_means[[i]] <- ppc_stat_grouped(
    y = y_obs[, i], 
    yrep = y_rep_cond[, , i], 
    group = condiciones, 
    stat = "mean"
  ) +
    ggtitle(paste0("Marca", colnames(y_obs)[i]))
}
for (i in 1:20){
  print(ppc_means[[i]])
}
# Para visualizar todos juntos (paquete gridExtra'):
library(gridExtra)
do.call(grid.arrange, ppc_means)

# -------------------------------------------------
# P-valores bayesianos por PTM y condición (media)
# -------------------------------------------------

# Índices de réplicas por condición (ya definidos antes)
# idx_asin <- c(1, 3, 5)
# idx_mito <- c(2, 4, 6)

K <- ncol(y_obs)                    # número de PTM (20)
S_comun <- dim(y_rep_comun)[1]
S_cond  <- dim(y_rep_cond)[1]

# Matrices para guardar p-valores:
# filas = PTM, columnas = condición
p_bayes_mean_comun_async   <- numeric(K)
p_bayes_mean_comun_mitosis <- numeric(K)
p_bayes_mean_cond_async    <- numeric(K)
p_bayes_mean_cond_mitosis  <- numeric(K)

for (k in 1:K) {
  # Estadísticos observados: media por condición para la PTM k
  T_obs_async   <- mean(y_obs[idx_asin, k])
  T_obs_mitosis <- mean(y_obs[idx_mito, k])
  
  ## Modelo kappa común
  T_draws_async_comun   <- numeric(S_comun)
  T_draws_mitosis_comun <- numeric(S_comun)
  
  for (s in 1:S_comun) {
    # draw s, PTM k, condición Async / Mitosis
    T_draws_async_comun[s]   <- mean(y_rep_comun[s, idx_asin,  k])
    T_draws_mitosis_comun[s] <- mean(y_rep_comun[s, idx_mito,  k])
  }
  
  p_bayes_mean_comun_async[k]   <- mean(T_draws_async_comun   >= T_obs_async)
  p_bayes_mean_comun_mitosis[k] <- mean(T_draws_mitosis_comun >= T_obs_mitosis)
  
  ## Modelo kappa condicional
  T_draws_async_cond   <- numeric(S_cond)
  T_draws_mitosis_cond <- numeric(S_cond)
  
  for (s in 1:S_cond) {
    T_draws_async_cond[s]   <- mean(y_rep_cond[s, idx_asin,  k])
    T_draws_mitosis_cond[s] <- mean(y_rep_cond[s, idx_mito,  k])
  }
  
  p_bayes_mean_cond_async[k]   <- mean(T_draws_async_cond   >= T_obs_async)
  p_bayes_mean_cond_mitosis[k] <- mean(T_draws_mitosis_cond >= T_obs_mitosis)
}

# Resumen en data.frame con nombres de PTM
ptm_names <- colnames(y_obs)

pvals_ptm <- data.frame(
  PTM                        = ptm_names,
  p_mean_comun_async         = p_bayes_mean_comun_async,
  p_mean_comun_mitosis       = p_bayes_mean_comun_mitosis,
  p_mean_cond_async          = p_bayes_mean_cond_async,
  p_mean_cond_mitosis        = p_bayes_mean_cond_mitosis,
  row.names                  = NULL
)

print(pvals_ptm)
# Hit-mat p.bayesiano para la media por marca
library(tidyr)
library(dplyr)
library(ggplot2)

# Partimos de pvals_ptm
pvals_long <- pvals_ptm |>
  pivot_longer(
    cols = starts_with("p_mean_"),
    names_to  = "modelo_cond",
    values_to = "p_bayes"
  ) |>
  mutate(
    Modelo    = case_when(
      grepl("comun", modelo_cond) ~ "kappa común",
      grepl("cond",  modelo_cond) ~ "kappa condicional",
      TRUE                        ~ "otro"
    ),
    Condicion = case_when(
      grepl("async",   modelo_cond) ~ "Async",
      grepl("mitosis", modelo_cond) ~ "Mitosis",
      TRUE                          ~ "?"
    )
  )

# Construimos niveles únicos y en orden inverso
ptm_levels <- rev(unique(pvals_long$PTM))

pvals_long <- pvals_long |>
  mutate(
    PTM = factor(PTM, levels = ptm_levels),
    ModeloCond = interaction(Modelo, Condicion, sep = " - ")
  )

ggplot(pvals_long, aes(x = ModeloCond, y = PTM, fill = p_bayes)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    limits   = c(0, 1),
    midpoint = 0.5,
    low  = "blue",
    mid  = "white",
    high = "red",
    name = "p_Bayes"
  ) +
  labs(
    x = "Modelo × Condición",
    y = "PTM / Marca epigenética",
    title = "Heatmap de p-valores bayesianos por PTM, modelo y condición"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )



# -------------------------------------------------
# PPC para la variación composicional total
# -------------------------------------------------
ppc_var_comun_a <- posterior_T_clr(y_rep_comun[,idx_asin, ], Y_obs = y_obs[idx_asin,])
ppc_var_comun_m <- posterior_T_clr(y_rep_comun[,idx_mito, ], Y_obs = y_obs[idx_mito,])
ppc_var_cond_a  <- posterior_T_clr(y_rep_cond[,idx_asin,],  Y_obs = y_obs[idx_asin,])
ppc_var_cond_m <- posterior_T_clr(y_rep_cond[,idx_mito, ], Y_obs = y_obs[idx_mito,])
ppc_var_comun_a$p_bayes
ppc_var_comun_m$p_bayes
ppc_var_cond_a$p_bayes
ppc_var_cond_m$p_bayes

# Visualización con bayesplot del estadístico T_clr
T_mat_comun <- cbind(
  `Kappa común-Asinc`       = ppc_var_comun_a$T_draws,
  `Kappa comun-Mitosis`     = ppc_var_comun_m$T_draws
)

p_T_comun <- mcmc_areas(
  as.data.frame(T_mat_comun),
  prob = 0.8
) +
  vline_at(ppc_var_comun_a$T_obs, color = "red") +
  ggtitle("PPC varianza total CLR: epistasia-comun") +
  xlab("T_clr (suma de varianzas CLR)")

print(p_T_comun)


T_mat_cond <- cbind(
  `Kappa cond-Asinc`       = ppc_var_cond_a$T_draws,
  `Kappa cond-Mitosis`     = ppc_var_cond_m$T_draws
)

p_T_cond <- mcmc_areas(
  as.data.frame(T_mat_cond),
  prob = 0.8
) +
  vline_at(ppc_var_cond_m$T_obs, color = "red") +
  ggtitle("PPC varianza total CLR: epistasia-cond") +
  xlab("T_clr (suma de varianzas CLR)")

print(p_T_cond)
