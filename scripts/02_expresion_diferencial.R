# 02_expresion_diferencial.R
# Objetivo: encontrar miRNA distintos entre GBM y controles sanos (HC)

# 1. Cargar lo que guardó el script 01
objeto  <- readRDS("results/gse122488_filtrado.rds")
datos_f <- objeto$datos
grupo   <- objeto$grupo

# 2. Pasar a escala log2
logcpm <- log2(datos_f + 1)

# 3. Quedarnos solo con las muestras GBM y HC (dejamos fuera GII_III)
sel <- grupo %in% c("GBM", "HC")
x   <- logcpm[, sel]
g   <- factor(grupo[sel], levels = c("HC", "GBM"))
table(g)  

# 4. Prueba de Wilcoxon para cada miRNA (cada fila de x)
p <- apply(x, 1, function(v) {
  wilcox.test(v[g == "GBM"], v[g == "HC"], exact = FALSE)$p.value
})

# 5. Cambio en log2: media en GBM menos media en HC
log2fc <- rowMeans(x[, g == "GBM"]) - rowMeans(x[, g == "HC"])

# 6. Corregir por pruebas múltiples (Benjamini-Hochberg)
fdr <- p.adjust(p, method = "BH")

# 7. Juntar todo en una tabla y ordenarla
resultados <- data.frame(miRNA = rownames(x), log2FC = log2fc,
                         p = p, FDR = fdr)
resultados$direccion <- ifelse(resultados$FDR >= 0.05, "no significativo",
                               ifelse(resultados$log2FC > 0, "sube en GBM", "baja en GBM"))
resultados <- resultados[order(resultados$FDR), ]

# 8. Resumen
print(table(resultados$direccion))   
head(resultados, 10)

# 9. Guardar
write.csv(resultados, "results/DE_GBM_vs_HC.csv", row.names = FALSE)
saveRDS(resultados, "results/DE_GBM_vs_HC.rds")