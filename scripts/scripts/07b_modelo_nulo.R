# 07b_modelo_nulo.R
# Objetivo: comparar el solapamiento observado con el de miRNA aleatorios

set.seed(2026)
d   <- readRDS("results/dianas_consenso.rds")
cn  <- readRDS("results/dianas_nulo_consenso.rds")
tt1 <- readRDS("results/DE_resistencia_GSE151680.rds")
tt2 <- readRDS("results/DE_resistencia_GSE327029.rds")

# miRNA reales con al menos 1 diana de consenso, por dirección
obs <- d$consenso
mir_up   <- unique(obs$mature_mirna_id[obs$direccion == "sube en GBM"])
mir_down <- unique(obs$mature_mirna_id[obs$direccion == "baja en GBM"])

# Fondo del nulo: miRNA -> sus genes diana
dianas_pool <- lapply(split(cn$target_symbol, cn$mature_mirna_id), unique)
cat("miRNA en el fondo del nulo:", length(dianas_pool), "\n")

nulo <- function(tt, firma, n_mir, dianas_obs, B = 20000, tol = 0.20) {
  univ  <- tt$gen
  firma <- intersect(firma, univ)
  N <- length(univ); K <- length(firma)
  dianas_obs <- intersect(dianas_obs, univ)
  n_obs <- length(dianas_obs)
  k_obs <- length(intersect(dianas_obs, firma))
  razon_obs <- k_obs / (n_obs * K / N)
  razones <- numeric(0)
  for (b in 1:B) {
    ids <- sample(names(dianas_pool), n_mir, replace = TRUE)
    g <- intersect(unique(unlist(dianas_pool[ids])), univ)
    n <- length(g)
    if (abs(n - n_obs) <= tol * n_obs) {            # mismo nº de dianas (±20 %)
      razones <- c(razones, length(intersect(g, firma)) / (n * K / N))
    }
  }
  data.frame(razon_obs = round(razon_obs, 2), sorteos_validos = length(razones),
             media_nulo = round(mean(razones), 2),
             p_empirico = (1 + sum(razones >= razon_obs)) / (1 + length(razones)))
}

f1_down <- tt1$gen[which(tt1$adj.P.Val < 0.05 & tt1$logFC <= -1)]
f2_up   <- tt2$gen[which(tt2$adj.P.Val < 0.05 & tt2$logFC >=  1)]

r1 <- nulo(tt1, f1_down, length(mir_up),   d$up)
r2 <- nulo(tt2, f2_up,   length(mir_down), d$down)
res_nulo <- rbind(cbind(comparacion = "GSE151680: miRNA suben -> firma baja", r1),
                  cbind(comparacion = "GSE327029: miRNA bajan -> firma sube", r2))
print(res_nulo)
write.csv(res_nulo, "results/modelo_nulo.csv", row.names = FALSE)