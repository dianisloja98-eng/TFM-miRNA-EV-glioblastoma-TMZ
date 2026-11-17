# miRNA en vesículas extracelulares de suero en glioblastoma y resistencia a temozolomida

Trabajo Fin de Máster – Máster en Bioinformática 
Autora: Diana Campoverde

## Descripción
Análisis in silico para identificar miRNA diferencialmente expresados en
vesículas extracelulares de suero de pacientes con glioblastoma y evaluar si
sus genes diana predichos se solapan con una firma transcriptómica de
resistencia a la temozolomida.

## Datos
- Dataset principal: GEO GSE122488 (matriz normalizada de miRNA, CPM).
  Descargar desde GEO y guardar en una carpeta `data/`.

## Estructura
- `scripts/` – código en R, numerado por módulos.
- `results/` – tablas de resultados.

## Módulos
1. `01_cargar_y_QC.R` – carga de datos, control de calidad y filtrado.
2. `02_expresion_diferencial.R` – miRNA diferencialmente expresados (GBM vs HC).

## Requisitos
R 4.6.2. Paquetes: BiocManager, GEOquery, limma, multiMiR, ggplot2, dplyr.
