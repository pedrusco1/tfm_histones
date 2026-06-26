---

editor_options: 
  markdown: 
    wrap: 72
---

# RProyecto_composicional

Repositorio del Trabajo Fin de Máster **“Modelos bayesianos de respuesta composicional para el análisis de proteoformas del fragmento 9–17 de la histona H3 en la transición a mitosis”**.

Este proyecto estudia la reorganización composicional de 20 proteoformas del fragmento H3(9–17), cuantificadas mediante LC-MS/MS en células madre embrionarias humanas en estado asíncrono y en mitosis. El análisis combina estadística bayesiana, análisis de datos composicionales y modelización de interacciones epistáticas entre modificaciones postraduccionales (PTMs).


```text
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

- Para facilitar la ejecución de los código se ha incluido un script de R que se llama 00_run_todo_secuencial.R. Este escript llama a los restantes y se generan las carpetas y archivos inermedios que se utlilizarán para compilar el informe.
- Para ejercutarlo, desde la carpeta raiz del proyecto, hacer en Rstudio

source("ScriptsR/00_run_todo_secuencial.R")
