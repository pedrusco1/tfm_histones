# RProyecto_composicional

Repositorio del Trabajo Fin de Máster **“Modelos bayesianos de respuesta composicional para el análisis de proteoformas del fragmento 9–17 de la histona H3 en la transición a mitosis”**.

Este proyecto estudia la reorganización composicional de 20 proteoformas del fragmento H3(9–17), cuantificadas mediante LC-MS/MS en células madre embrionarias humanas en estado asíncrono y en mitosis. El análisis combina estadística bayesiana, análisis de datos composicionales y modelización de interacciones epistáticas entre modificaciones postraduccionales (PTMs).

## Objetivos

El objetivo general es caracterizar cómo cambian las proporciones relativas de las proteoformas del fragmento H3(9–17) entre el estado asíncrono y la mitosis, prestando especial atención al papel de la fosforilación en S10, las acetilaciones y su interacción con los estados de metilación en K9.

Desde el punto de vista metodológico, se realiza:

- análisis exploratorio composicional;

- Comparación dre varias familias de modelos:

  - modelos lineales moderados sobre transformaciones CLR/ILR;

  - regresión de Dirichlet;

  - modelos bayesianos de Dirichlet con campo epistático;

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

- \`DatosProcesados/\`: Datos proceasdos preparados para el análisis

- \`Funciones/\`: funciones auxiliares usadas por los scripts principales.

- \`Informes/Proyecto/\`: documentos del manuscrito y materiales derivados.

- \`ScriptsR/\`: scripts de carga de datos, ajuste de modelos, comparación y generación de resultados.

- \`Stan/\`: modelos probabilísticos implementados en Stan.

## Reproducibilidad

Para reproducir este proyecto, seguir las siguientes instrucciones:

1.  **Clonar el repositorio**, en una terminal:

    ```{bash}
    git clone https://github.com/usuario/tfm_histones.git cd tfm_histones
    ```

2.  **Abrir el proyecto de RStudio**

    - Abrir RStudio.

    - Ir a `File` → `Open Project...`.

    - Navegar hasta la carpeta `tfm_histones/` recién clonada.

    - Seleccionar el archivo `RProyecto_composicional.Rproj`.

    > Nota: El nombre del proyecto de RStudio (`RProyecto_composicional.Rproj`) no coincide con el nombre de la carpeta (`tfm_histones`), pero esto no afecta al funcionamiento. Lo importante es abrir ese archivo `.Rproj` dentro de la carpeta clonada.

3.  **Restaurar el entorno de paquetes con `renv`** En la consola de R:

    ```{r}
    install.packages("renv")  # solo si no está disponible renv::restore()
    install.packages("BiocManager") # necesario para paquetes de Bioconductor 
    renv::restore()
    ```

    Esto reinstala las versiones de paquetes registradas en `renv.lock` para que el entorno sea el mismo que el usado en el análisis.

4.  **Ejecutar el pipeline de análisis** En la consola de R, ya dentro del proyecto:

    ```{r}
    source("ScriptsR/00_run_todo_secuencial.R")
    ```

    Este script ejecuta de forma secuencial todos los scripts de `ScriptsR/` y genera los resultados necesarios en las carpetas de datos procesados.

- Los scripts están organizados en bloques que cubren:

  1.  análisis exploratorio y expresión diferencial;

  2.  modelos Dirichlet;

  3.  modelos epistáticos;

  4.  modelos logístico‑normales;

  5.  comparación entre modelos y resúmenes finales.

6.  **Generar el PDF de la memoria.** Desde la terminal de RStudio:

    ```{bash}
    cd Informes/Proyecto
    quarto render
    ```

    Esto compila el proyecto Quarto y genera el PDF final en `Informes/Proyecto/_book/`.

## Documento de referencia

La descripción completa del contexto biológico, los métodos estadísticos, los resultados y la discusión se encuentra en el manuscrito del TFM, generado en la ruta, \~/Informes/Proyecto/\_book y que se llama “Modelos bayesianos de respuesta composicional para el análisis de proteoformas del fragmento 9–17 de la histona H3 en la transición a mitosis”.

# Resumen breve de los resultados

Considerado en su conjunto, este trabajo pone de manifiesto que el paso a la mitosis lleva consigo un aumento relativo de las proteoformas que presentan fosforilación en S10, una reorganización de la acetilación supeditada a ciertas condiciones y y una menor influencia de los efectos epistásicos, lo cual concuerda con un sistema de carácter más aditivo durante la división celular.
