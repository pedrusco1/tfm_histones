#################################################
#                Configurar Stan                #
#################################################

# Cargar librerías
library(rstan)
library(loo)
library(Matrix)
library(stringr)
library(dplyr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(123)

# Rutas
root <- getwd()
stan_dir <- file.path(root, "Stan")
data_dir <- file.path(root, "DatosProcesados")
