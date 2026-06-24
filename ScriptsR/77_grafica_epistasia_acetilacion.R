resultados_epistasia_dia <- readRDS(here::here("DatosProcesados", 
                   "epistasia_descomp_dia",
                   "77_resultados_epistasia_descomposicion_varianza_dia.rds"))

library(dplyr)
library(stringr)

df_filtrado <- resultados_epistasia_dia$df_resumen_epistasia_local %>%
  filter(str_detect(Cuadrado, "K9ac:K14ac"))
df_filtrado %>% select(Condicion,Cuadrado, Etiqueta)
df_epsitasia_ac_no_ph <- df_filtrado %>% 
  filter(str_detect(Cuadrado, "S10un"))
df_epsitasia_ac_si_ph <- df_filtrado %>% 
  filter(str_detect(Cuadrado, "S10ph"))
df_filtrado %>% 
  select(Etiqueta)
df_filtrado <- df_filtrado %>%
  mutate(K9 = str_extract(Etiqueta, "K9[^→ ]+"))
df_filtrado %>% select(Condicion, K9, Etiqueta)

library(ggplot2)

df_plot <- df_filtrado %>%
  mutate(
    S10 = case_when(
      str_detect(Etiqueta, "S10un") ~ "S10un",
      str_detect(Etiqueta, "S10ph") ~ "S10ph",
      TRUE ~ NA_character_
    )
  )



ggplot(df_plot, aes(x = K9, y = Media, color = Condicion, group = Condicion)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = IC_Lower, ymax = IC_Upper), width = 0.1) +
  facet_wrap(~ S10) +
  theme_bw(base_size = 14) +
  labs(
    x = "Estado inicial de K9",
    y = "Media (con IC)",
    color = "Condición",
    title = "Efectos epistáticos K9ac:K14ac separados por estado de S10"
  )




