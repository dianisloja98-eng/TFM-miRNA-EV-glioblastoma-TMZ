# 06_sensibilidad.R
# Análisis EXPLORATORIO: solapamiento de las dianas con la firma de cada dataset por separado

tt1 <- readRDS("results/DE_resistencia_GSE151680.rds")
tt2 <- readRDS("results/DE_resistencia_GSE327029.rds")
tt3 <- readRDS("results/DE_resistencia_GSE68029.rds")
d   <- readRDS("results/dianas_consenso.rds")

test_hiper <- function(dianas, firma, universo) {
  dianas <- intersect(dianas, universo)
  firma  <- intersect(firma,  universo)
  k <- length(intersect(dianas, firma))
  N <- length(universo); K <- length(firma); n <- length(dianas)
  esperado <- n * K / N
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(dianas = n, firma = K, solapan = k,
             esperado = round(esperado, 1), razon = round(k / esperado, 2), p = p)
}

lista_tt <- list(GSE151680 = tt1, GSE327029 = tt2, GSE68029 = tt3)
filas <- list()
for (nombre in names(lista_tt)) {
  tt     <- lista_tt[[nombre]]
  univ   <- tt$gen                                          # universo = genes de ese dataset
  f_up   <- tt$gen[which(tt$adj.P.Val < 0.05 & tt$logFC >=  1)]
  f_down <- tt$gen[which(tt$adj.P.Val < 0.05 & tt$logFC <= -1)]
  filas[[nombre]] <- rbind(
    cbind(dataset = nombre, comparacion = "miRNA suben -> firma baja", test_hiper(d$up,   f_down, univ)),
    cbind(dataset = nombre, comparacion = "miRNA bajan -> firma sube", test_hiper(d$down, f_up,   univ))
  )
}
tabla_sens <- do.call(rbind, filas)
tabla_sens$p_BH <- p.adjust(tabla_sens$p, method = "BH")    # corrección por 6 pruebas
print(tabla_sens)
write.csv(tabla_sens, "results/solapamiento_sensibilidad.csv", row.names = FALSE)