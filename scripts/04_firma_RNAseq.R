# 04_firma_RNAseq.R
# Objetivo: genes diferenciales resistente vs sensible en GSE151680 y GSE327029

library(limma)

# Función 1: si un símbolo de gen se repite, quedarse con la fila de mayor expresión media
colapsar <- function(m, genes) {
  orden <- order(rowMeans(m), decreasing = TRUE)
  m <- m[orden, ]; genes <- genes[orden]
  quedan <- !duplicated(genes)
  m <- m[quedan, ]
  rownames(m) <- genes[quedan]
  m
}

# --- Función 2: análisis con limma, ajustando por línea celular
analizar <- function(expr, resistencia, linea) {
  expr <- expr[rowSums(expr >= 1) >= 3, ]          # filtro de expresión mínima
  logx <- log2(expr + 1)                           # escala log2
  resistencia <- factor(resistencia, levels = c("sensible", "resistente"))
  diseno <- model.matrix(~ linea + resistencia)
  ajuste <- eBayes(lmFit(logx, diseno), trend = TRUE, robust = TRUE)
  tt <- topTable(ajuste, coef = "resistenciaresistente", number = Inf, sort.by = "none")
  tt$gen <- rownames(tt)
  tt
}

# ================= GSE151680 =================
a <- read.delim("data/GSE151680_gene.expression.symbol_matrix.txt", check.names = FALSE)
m1 <- colapsar(as.matrix(a[, -1]), a[[1]])
res1 <- grepl("TR", colnames(m1))                   # U251TR y U87TR = resistentes
lin1 <- ifelse(grepl("^U251", colnames(m1)), "U251", "U87")
tt1 <- analizar(m1, ifelse(res1, "resistente", "sensible"), lin1)

# ================= GSE327029 =================
b <- read.delim("data/GSE327029_mRNA-12.TPM.genes.tsv", check.names = FALSE)
m2 <- colapsar(as.matrix(b[, -(1:2)]), b$gene_name)
res2 <- grepl("R-[0-9]$", colnames(m2))             # SF126R y U118MGR = resistentes
lin2 <- ifelse(grepl("^SF126", colnames(m2)), "SF126", "U118MG")
tt2 <- analizar(m2, ifelse(res2, "resistente", "sensible"), lin2)

# ================= Comprobaciones =================
print(table(res1, lin1)); print(table(res2, lin2))   # deben ser 3 muestras en cada casilla
cat("genes analizados GSE151680:", nrow(tt1), "\n")
cat("genes analizados GSE327029:", nrow(tt2), "\n")
cat("significativos GSE151680 (FDR<0.05 y |log2FC|>=1):",
    sum(tt1$adj.P.Val < 0.05 & abs(tt1$logFC) >= 1), "\n")
cat("significativos GSE327029 (FDR<0.05 y |log2FC|>=1):",
    sum(tt2$adj.P.Val < 0.05 & abs(tt2$logFC) >= 1), "\n")

# ================= Guardar =================
saveRDS(tt1, "results/DE_resistencia_GSE151680.rds")
saveRDS(tt2, "results/DE_resistencia_GSE327029.rds")