# 05_solapamiento.R
# Objetivo: firma de resistencia a TMZ y solapamiento con las dianas de los miRNA

# 1. Cargar resultados de los módulos 3b y 4
tt1 <- readRDS("results/DE_resistencia_GSE151680.rds")
tt2 <- readRDS("results/DE_resistencia_GSE327029.rds")
tt3 <- readRDS("results/DE_resistencia_GSE68029.rds")
d   <- readRDS("results/dianas_consenso.rds")   

# 2. Universo: genes analizados en al menos 2 de los 3 datasets
cuenta_test <- table(c(tt1$gen, tt2$gen, tt3$gen))
universo <- names(cuenta_test)[cuenta_test >= 2]
cat("genes en el universo:", length(universo), "\n")

# 3. Firma: significativo (FDR<0.05 y |log2FC|>=1) en >= 2 datasets, mismo sentido
marca_up   <- function(tt) universo %in% tt$gen[which(tt$adj.P.Val < 0.05 & tt$logFC >=  1)]
marca_down <- function(tt) universo %in% tt$gen[which(tt$adj.P.Val < 0.05 & tt$logFC <= -1)]
n_up   <- marca_up(tt1)   + marca_up(tt2)   + marca_up(tt3)
n_down <- marca_down(tt1) + marca_down(tt2) + marca_down(tt3)
firma_up   <- universo[n_up   >= 2]    # suben en la resistencia
firma_down <- universo[n_down >= 2]    # bajan en la resistencia
cat("firma: suben =", length(firma_up), "| bajan =", length(firma_down), "\n")

# 4. Coherencia entre datasets (correlación de Spearman de los cambios log2FC)
coh <- function(a, b) {
  comun <- intersect(a$gen, b$gen)
  cor(a$logFC[match(comun, a$gen)], b$logFC[match(comun, b$gen)], method = "spearman")
}
cat("correlación 151680 vs 327029:", round(coh(tt1, tt2), 3), "\n")
cat("correlación 151680 vs 68029 :", round(coh(tt1, tt3), 3), "\n")
cat("correlación 327029 vs 68029 :", round(coh(tt2, tt3), 3), "\n")

# 5. Test hipergeométrico (una cola: "más solapamiento que el esperado")
test_hiper <- function(dianas, firma) {
  dianas <- intersect(dianas, universo)
  firma  <- intersect(firma,  universo)
  k <- length(intersect(dianas, firma))          # genes que están en ambos
  N <- length(universo); K <- length(firma); n <- length(dianas)
  esperado <- n * K / N
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(dianas = n, firma = K, solapan = k,
             esperado = round(esperado, 1),
             razon = round(k / esperado, 2), p = p)
}

# 6. 4 comparaciones: 2 según la hipótesis y 2 de control
tabla <- rbind(
  cbind(comparacion = "miRNA suben -> firma baja (hipótesis)", test_hiper(d$up,   firma_down)),
  cbind(comparacion = "miRNA bajan -> firma sube (hipótesis)", test_hiper(d$down, firma_up)),
  cbind(comparacion = "miRNA suben -> firma sube (control)",   test_hiper(d$up,   firma_up)),
  cbind(comparacion = "miRNA bajan -> firma baja (control)",   test_hiper(d$down, firma_down))
)
tabla$p_BH <- p.adjust(tabla$p, method = "BH")   # corrección por 4 pruebas
print(tabla)

# 7. Guardar tabla y lista de genes que solapan
dir.create("results", showWarnings = FALSE)
write.csv(tabla, "results/solapamiento_hipergeometrico.csv", row.names = FALSE)
solapan <- list(suben_miRNA_baja_firma = intersect(intersect(d$up, universo), firma_down),
                bajan_miRNA_sube_firma = intersect(intersect(d$down, universo), firma_up))
saveRDS(list(firma_up = firma_up, firma_down = firma_down, solapan = solapan,
             universo = universo), "results/firma_y_solapamiento.rds")
print(lapply(solapan, head, 20))