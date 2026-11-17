# miRNA en vesículas extracelulares de suero en glioblastoma y resistencia a temozolomida

Trabajo Fin de Máster – Máster en Bioinformática

Autora: Diana Campoverde

## Descripción
Análisis in silico para identificar miRNA diferencialmente expresados en
vesículas extracelulares de suero de pacientes con glioblastoma y evaluar si
sus genes diana predichos se solapan con una firma transcriptómica de
resistencia a la temozolomida (TMZ). Es un trabajo exploratorio, orientado a
generar hipótesis.

## Datos
Todos los datos son públicos y se descargan de GEO (no se incluyen en el repositorio):

- **GSE122488**: miRNA en vesículas extracelulares de suero (controles sanos,
  glioblastoma y glioma de grado II-III), valores en CPM (Conteos Por Millón).
- **GSE151680**, **GSE327029** y **GSE68029**: células de glioma sensibles y resistentes
  a TMZ, usadas para construir la firma de resistencia.

Los archivos de GSE122488, GSE151680 y GSE327029 se guardan en una carpeta `data/`.
GSE68029 se descarga desde R con el paquete GEOquery.

## Estructura
- `scripts/` – código en R, numerado por módulos.
- `results/` – tablas de resultados en formato CSV.
- `results/figuras/` – figuras en formato PNG.

## Módulos
1. `01_cargar_y_QC.R`: carga de datos, control de calidad y filtrado.
2. `02_expresion_diferencial.R`: miRNA diferencialmente expresados (glioblastoma vs controles).
3. `03_dianas.R` y `03b_filtrar_dianas.R`: dianas predichas (multiMiR) y filtro por consenso (>= 3 bases de datos).
4. `04_firma_RNAseq.R`, `04a_cargar_GSE68029.R` y `04b_firma_GSE68029.R`: genes diferenciales resistente vs sensible en cada dataset.
5. `05_solapamiento.R`: firma de consenso y test hipergeométrico.
6. `06_sensibilidad.R`: análisis exploratorio del solapamiento por dataset.
7. `07a_dianas_nulo.R` y `07b_modelo_nulo.R`: modelo nulo con miRNA aleatorios.
8. `08_figuras.R`: figuras de resultados.

## Requisitos
R 4.6.2. Paquetes: BiocManager, GEOquery, limma, multiMiR, ggplot2.
Ejecutar los scripts en orden, abriendo el proyecto desde la carpeta principal.

## Estado
Trabajo en curso. Los resultados son exploratorios y deben interpretarse con cautela.
