###########################################
#       Cargar datos H3_9_17              #
###########################################
source(here::here("ScriptsR", "40_stan_config.R"))

load(here::here("DatosProcesados", "df_h3917.Rdata"))

Y_mat <- as.matrix(df.ancho[, 3:22])

cond <- as.integer(factor(df.ancho$Estado,
                          levels = c("asinc", "mitot")))

save(Y_mat, cond, file = file.path(data_dir, "data_inputs_h3917.Rdata"))
