# 03b_filtrar_dianas.R
# Objetivo: quedarnos con dianas de consenso y separarlas por dirección del miRNA

tabla <- readRDS("results/dianas_predichas.rds")
res   <- readRDS("results/DE_GBM_vs_HC.rds")

# 1. Cuántos pares quedan si exigimos 1, 2, 3, 4 o 5 bases de datos
sapply(1:5, function(k) sum(tabla$predicted.sum >= k))

# 2. Consenso: pares predichos por al menos 3 bases
cons <- tabla[tabla$predicted.sum >= 3 & !is.na(tabla$target_symbol) &
                tabla$target_symbol != "", ]

# 3. Añadir la dirección del miRNA (sube / baja en GBM)
cons$direccion <- res$direccion[match(cons$mature_mirna_id, res$miRNA)]
print(table(cons$direccion, useNA = "ifany"))

# 4. Genes diana únicos por dirección
dianas_up   <- unique(cons$target_symbol[cons$direccion == "sube en GBM"])
dianas_down <- unique(cons$target_symbol[cons$direccion == "baja en GBM"])
cat("genes diana de miRNA que suben:", length(dianas_up), "\n")
cat("genes diana de miRNA que bajan:", length(dianas_down), "\n")

# 5. Cuántos genes diana tiene cada miRNA 
n_dianas <- table(cons$mature_mirna_id)
print(summary(as.numeric(n_dianas)))

# 6. Guardar
saveRDS(list(consenso = cons, up = dianas_up, down = dianas_down),
        "results/dianas_consenso.rds")