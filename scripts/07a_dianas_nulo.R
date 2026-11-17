# 07a_dianas_nulo.R
# Objetivo: obtener dianas de miRNA NO significativos para el modelo nulo

library(multiMiR)
set.seed(123)                 

res  <- readRDS("results/DE_GBM_vs_HC.rds")
pool <- res$miRNA[res$direccion == "no significativo"]
elegidos <- sample(pool, 150)                    # 150 miRNA al azar

lotes <- split(elegidos, ceiling(seq_along(elegidos) / 10))
lista_n <- list()
for (i in seq_along(lotes)) {
  cat("Lote", i, "de", length(lotes), "\n")
  r <- get_multimir(org = "hsa", mirna = lotes[[i]],
                    table = "predicted", summary = TRUE,
                    predicted.cutoff.type = "p", predicted.cutoff = 20)
  lista_n[[i]] <- r@summary
  saveRDS(lista_n, "results/dianas_nulo_parcial.rds")
}

# Unificar columnas, unir y aplicar el mismo consenso (>= 3 bases)
cols <- Reduce(union, lapply(lista_n, colnames))
lista_n2 <- lapply(lista_n, function(x) {
  for (f in setdiff(cols, colnames(x))) x[[f]] <- 0
  x[, cols]
})
tabla_n <- do.call(rbind, lista_n2)
cons_n <- tabla_n[tabla_n$predicted.sum >= 3 & !is.na(tabla_n$target_symbol) &
                    tabla_n$target_symbol != "", ]
saveRDS(cons_n, "results/dianas_nulo_consenso.rds")
cat("pares:", nrow(cons_n), "| miRNA con dianas:", length(unique(cons_n$mature_mirna_id)), "\n")