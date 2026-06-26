# RProyecto_composicional

Repositorio del Trabajo Fin de Máster **“Modelos bayesianos de respuesta composicional para el análisis de proteoformas del fragmento 9–17 de la histona H3 en la transición a mitosis”**.

Este proyecto estudia la reorganización composicional de 20 proteoformas del fragmento H3(9–17), cuantificadas mediante LC-MS/MS en células madre embrionarias humanas en estado asíncrono y en mitosis. El análisis combina estadística bayesiana, análisis de datos composicionales y modelización de interacciones epistáticas entre modificaciones postraduccionales (PTMs).

## Objetivos

El objetivo general es caracterizar cómo cambian las proporciones relativas de las proteoformas del fragmento H3(9–17) entre el estado asíncrono y la mitosis, prestando especial atención al papel de la fosforilación en S10, las acetilaciones y su interacción con los estados de metilación en K9.

Desde el punto de vista metodológico, se realiza:

-  análisis exploratorio composicional;

- Comparación dre varias familias de modelos:

  -  modelos lineales moderados sobre transformaciones CLR/ILR;

  - regresión de Dirichlet;

  -  modelos bayesianos de Dirichlet con campo epistático;

  - modelos logístico‑normales con estructura epistática y regularización.

## Estructura del repositroio

``` text
RProyecto_composicional/
├── DatosProcesados/
├── Funciones/
├── Informes/
│   └── Proyecto/
│       └── .quarto/
│       ├── _book/
│       └── Imagenes/
├── ScriptsR/
└── Stan/
```

-  \`DatosProcesados/\`: Datos proceasdos preparados para el análisis

-  \`Funciones/\`: funciones auxiliares usadas por los scripts principales.

-  \`Informes/Proyecto/\`: documentos del manuscrito y materiales derivados.

-  \`ScriptsR/\`: scripts de carga de datos, ajuste de modelos, comparación y generación de resultados.

-  \`Stan/\`: modelos probabilísticos implementados en Stan.

## Reproducibilidad

- Para facilitar la ejecución de los código se ha incluido un script de R que se llama `00_run_todo_secuencial.R`. Este escript llama a los restantes y se generan las carpetas y archivos inermedios que se utlilizarán para compilar el informe.

- Para ejercutarlo, desde la carpeta raíz del proyecto, hacer en Rstudio:

```{r}
source("ScriptsR/00_run_todo_secuencial.R")
```

- Los scripts están organizados en bloques que cubren:

  1.   análisis exploratorio y expresión diferencial;

  2.   modelos Dirichlet;

  3.   modelos epistáticos;

  4.   modelos logístico‑normales;

  5.  comparación entre modelos y resúmenes finales.

## Documento de referencia

La descripción completa del contexto biológico, los métodos estadísticos, los resultados y la discusión se encuentra en el manuscrito del TFM, incluido en la ruta, \~/Informes/Proyecto/\_book “Modelos bayesianos de respuesta composicional para el análisis de proteoformas del fragmento 9–17 de la histona H3 en la transición a mitosis”.

# Resumen breve de los resultados

Considerado en su conjunto, este trabajo pone de manifiesto que el paso a la mitosis lleva consigo un aumento relativo de las proteoformas que presentan fosforilación en S10, una reorganización de la acetilación supeditada a ciertas condiciones y y una menor influencia de los efectos epistásicos, lo cual concuerda con un sistema de carácter más aditivo durante la división celular.
