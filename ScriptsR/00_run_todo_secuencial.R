# 00_run_todo_secuencial.R

scripts <- c(
  "ScriptsR/31_DesCompoH3917_corregido-2.R",
  "ScriptsR/32_expresion_diferencial-3.R",
  "ScriptsR/31_cargar_datos_dirichreg_h3917.R",
  "ScriptsR/33_fit_dirichreg_models_h3917-4.R",
  "ScriptsR/34_loo_waic_diricreg-5.R",
  "ScriptsR/35_ppc_dirichreg-6.R",
  "ScriptsR/36_kl_dirichlet_cuatro_modelos-7.R",
  "ScriptsR/40_epistasia_matriz_B_h3917-8.R",
  "ScriptsR/40_stan_config-9.R",
  "ScriptsR/41_cargar_datos_h3917-10.R",
  "ScriptsR/42_matriz_delta_h3917-7.R",
  "ScriptsR/43_fit_models_epistasia-8.R",
  "ScriptsR/44_compara_epis_model_h3917-9.R",
  "ScriptsR/44_epistasia_ppc_h3917-10.R",
  "ScriptsR/45_loo_waic_epistasia.R",
  "ScriptsR/46_kl_epistasia_cuatro_modelos-2.R",
  "ScriptsR/53_fit_lognor_sum0_model_revisado-3.R",
  "ScriptsR/54_loo_waic_lognor_sum0_revisado-4.R",
  "ScriptsR/55_ppc_lognor_sum0_revisado-5.R",
  "ScriptsR/57_kl_lognormal_sum0_cuatro_modelos_revisado_corregido.R",
  "ScriptsR/70_compara_p_mito_gt_asinc_centros-2.R",
  "ScriptsR/75_loo_waic_12_modelos_revisado-3.R",
  "ScriptsR/77_grafica_epistasia_acetilacion-4.R",
  "ScriptsR/77_precalcular_epistasia_descomposicion_varianza_dia-5.R",
  "ScriptsR/78_precalcular_resultados_proporciones_epistasia_v4-6.R"
)

for (s in scripts) {
  cat("\n==============================\n")
  cat("Ejecutando:", s, "\n")
  source(s)
}
