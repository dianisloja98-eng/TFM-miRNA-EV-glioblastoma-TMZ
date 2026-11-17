# 04a_cargar_GSE68029.R
# Objetivo: descargar GSE68029 (microarray de células madre) y revisar su estructura

library(GEOquery)

# 1. Descargar la serie
gse <- getGEO("GSE68029", GSEMatrix = TRUE)
gse <- gse[[1]]

# 2. Matriz de expresión: filas = sondas, columnas = muestras
expr <- exprs(gse)
dim(expr)                       

# 3. Escala de los valores
summary(as.vector(expr))

# 4. Nombres de las muestras y columnas de anotación de las sondas
pData(gse)[, "title"]
colnames(fData(gse))