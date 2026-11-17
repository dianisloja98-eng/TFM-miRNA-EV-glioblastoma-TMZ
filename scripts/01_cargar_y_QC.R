# 1. Instalar paquetes
install.packages("BiocManager")          # gestor de paquetes de bioinformática
BiocManager::install(c("GEOquery",       # descarga datos de GEO
                       "limma",          # expresión diferencial
                       "multiMiR"))      # dianas de miRNA
install.packages(c("ggplot2", "dplyr"))

# 2. 01_cargar_y_QC.R (quality control)
# Objetivo: cargar la matriz de GSE122488, comprobar que está bien leída y filtrar

# 1. Ruta del archivo.
ruta <- "data/GSE122488_normalized_microRNA_counts.txt"

# 2. Leer la tabla. skip = 3 salta las 3 líneas de texto; la cuarta es el encabezado
datos <- read.delim(ruta, skip = 3, header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)

# 3. Poner el nombre del miRNA como nombre de fila
rownames(datos) <- datos[[1]]
datos <- datos[, -1]
colnames(datos) <- trimws(colnames(datos))      # quita espacios sobrantes
datos <- datos[, colnames(datos) != ""]         # quita columnas vacías si las hubiera

# 4. Asignar cada muestra a su grupo según su nombre (HC (healthy controls), GBM (glioblastoma) o GII-III (glioblastomas de grado II y III))
grupo <- ifelse(grepl("^HC", colnames(datos)), "HC",
                ifelse(grepl("^GBM", colnames(datos)), "GBM", "GII_III"))

# 5. Comprobaciones
cat("miRNA x muestras:", dim(datos), "\n")
print(table(grupo))
print(summary(colSums(datos)))   # cada muestra suma ~1,000,000 (son CPM)

# 6. Filtro: conservar miRNA con CPM (counts per million) >= 1 en al menos el 25 % de las muestras
filtro  <- rowSums(datos >= 1) >= ceiling(0.25 * ncol(datos))
datos_f <- datos[filtro, ]
cat("miRNA tras el filtro:", nrow(datos_f), "\n") 

# 7. Guardar el resultado 
dir.create("results", showWarnings = FALSE)
saveRDS(list(datos = datos_f, grupo = grupo), "results/gse122488_filtrado.rds")