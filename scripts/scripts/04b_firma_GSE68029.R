# 04b_firma_GSE68029.R
# Objetivo: genes diferenciales resistente vs parental en GSE68029 (microarray)

library(GEOquery)
library(limma)

# 1. Descargar y extraer la matriz
gse  <- getGEO("GSE68029", GSEMatrix = TRUE)[[1]]
expr <- exprs(gse)
cat("dimensiones:", dim(expr), "\n")
cat("rango de valores:", range(expr, na.rm = TRUE), "\n")
if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)  

# 2. Símbolo de gen de cada sonda 
anot <- fData(gse)
col_simbolo <- grep("symbol", colnames(anot), ignore.case = TRUE, value = TRUE)[1]
cat("columna de símbolo usada:", col_simbolo, "\n")
simbolo <- trimws(sub(" ///.*", "", anot[[col_simbolo]]))   
ok <- !is.na(simbolo) & simbolo != "" & simbolo != "---"
expr <- expr[ok, ]; simbolo <- simbolo[ok]

# 3. Una sola fila por gen: la sonda con mayor expresión media
orden  <- order(rowMeans(expr), decreasing = TRUE)
expr   <- expr[orden, ]; simbolo <- simbolo[orden]
quedan <- !duplicated(simbolo)
expr   <- expr[quedan, ]
rownames(expr) <- simbolo[quedan]
cat("genes tras colapsar sondas:", nrow(expr), "\n")

# 4. Grupos: resistente = clones que sobrevivieron a TMZ; línea = D431, S496 o E445
titulo      <- pData(gse)[, "title"]
resistencia <- factor(ifelse(grepl("TMZ", titulo), "resistente", "sensible"),
                      levels = c("sensible", "resistente"))
linea       <- substr(titulo, 1, 4)
print(table(resistencia, linea))        

# 5. limma ajustando por línea celular
diseno <- model.matrix(~ linea + resistencia)
ajuste <- eBayes(lmFit(expr, diseno), trend = TRUE, robust = TRUE)
tt3 <- topTable(ajuste, coef = "resistenciaresistente", number = Inf, sort.by = "none")
tt3$gen <- rownames(tt3)

# 6. Resumen y guardado
cat("significativos GSE68029 (FDR<0.05 y |log2FC|>=1):",
    sum(tt3$adj.P.Val < 0.05 & abs(tt3$logFC) >= 1), "\n")
saveRDS(tt3, "results/DE_resistencia_GSE68029.rds")