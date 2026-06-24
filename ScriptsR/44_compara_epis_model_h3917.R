#############################################
#     Comparar modelos epistasia H3 9-17    #
#############################################

source("ScriptsR/40_stan_config.R")
source("Funciones/fun_loo_stan_model.R")

load(file.path(data_dir, "stan_epis_h3917.Rdata"))

loo_cond  <- compute_loo(fit_kappa_cond)
loo_cond_hier <- compute_loo(fit_kappa_cond_hier)
loo_comun <- compute_loo(fit_kappa_comun)

print(loo::loo_compare(loo_cond, loo_cond_hier, loo_comun))

save(loo_cond, loo_comun, loo_cond_hier,
     file = file.path(data_dir, "loo_results_h3917.Rdata"))
