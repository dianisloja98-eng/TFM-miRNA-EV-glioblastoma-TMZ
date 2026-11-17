# 03_dianas.R
# Objetivo: obtener los genes diana predichos de los miRNA diferenciales

setwd("C:/Documentos/MAESTRIA/TFM/TRABAJO_FIN_DE_MASTER")
library(multiMiR)

# 1. Leer los resultados del módulo 2 y quedarnos con los significativos
res  <- readRDS("results/DE_GBM_vs_HC.rds")
sig  <- res[res$direccion != "no significativo", ]
cat("miRNA significativos:", nrow(sig), "\n")  

# 2. Consultar dianas predichas por lotes de 10 miRNA
lotes <- split(sig$miRNA, ceiling(seq_along(sig$miRNA) / 10))
lista <- list()
for (i in seq_along(lotes)) {
  cat("Lote", i, "de", length(lotes), "\n")
  r <- get_multimir(org = "hsa", mirna = lotes[[i]],
                    table = "predicted", summary = TRUE,
                    predicted.cutoff.type = "p", predicted.cutoff = 20)
  lista[[i]] <- r@summary
  saveRDS(lista, "results/dianas_parcial.rds")   
}

# 3. Unificar columnas
cols <- Reduce(union, lapply(lista, colnames))
lista2 <- lapply(lista, function(x) {
  for (f in setdiff(cols, colnames(x))) x[[f]] <- 0   
  x[, cols]
})

# 4. Unir todos los lotes y revisar
tabla <- do.call(rbind, lista2)
dim(tabla)
print(colnames(tabla))
head(tabla)

# 5. Guardar
saveRDS(tabla, "results/dianas_predichas.rds")