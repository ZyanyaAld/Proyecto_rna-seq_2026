#!/usr/bin/env Rscript

# 03_EDA.R
# Objetivo:
# 1) Analisis exploratorio (EDA)
# 2) Analisis de expresion diferencial (DEA)
# Comparacion:
# left atrial appendage vs right atrial appendage

# 1. CARGA DE LIBRERÍAS

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)


# 2. CARGAR DATOS FILTRADOS

load("results/rse_gene_SRP188219_filtered.RData")
dim(rse_gene_SRP188219)


# 3. NORMALIZACION edgeR

dge <- DGEList(
  counts = assay(rse_gene_SRP188219, "counts")
)

dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE)


# 4. EDA

# -------- PCA --------

pca <- prcomp(t(logCPM))

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Tissue = rse_gene_SRP188219$sra_attribute.source_name
)

p <- ggplot(pca_df, aes(PC1, PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw(base_size = 15) +
  ggtitle("PCA - RNA-seq SRP188219")

ggsave(
  "results/PCA_plot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


# -------- MDS (corregido) --------

png("results/MDS_plot.png", width = 1200, height = 1000, res = 200)

group <- factor(rse_gene_SRP188219$sra_attribute.source_name)
colors <- c("#1B9E77", "#D95F02")

plotMDS(dge, col = colors[group], pch = 16, cex = 1.5)

legend("topright", legend = levels(group), col = colors, pch = 16, cex = 1.2)

dev.off()


# -------- Distribucion global --------

png("results/Density_logCPM.png", width = 1000, height = 800, res = 150)

plotDensities(logCPM, legend = FALSE)

dev.off()


# 5. MODELO LINEAL (DEA)

model <- model.matrix(
  ~sra_attribute.source_name,
  data = colData(rse_gene_SRP188219)
)


# -------- VOOM --------

png("results/Voom_mean_variance.png", width = 1000, height = 800, res = 150)

v <- voom(dge, model, plot = TRUE)

dev.off()


# -------- LIMMA --------

fit <- lmFit(v, model)
fit <- eBayes(fit)


# 6. RESULTADOS

de_results <- topTable(
  fit,
  coef = 2,
  number = Inf
)

write.csv(
  de_results,
  file = "results/DE_results_SRP188219.csv",
  row.names = TRUE
)

# Conteo correcto de significativos
sig_genes <- sum(de_results$adj.P.Val < 0.05)

cat("Genes significativos (FDR < 0.05):", sig_genes, "\n")


# 7. MA PLOT

png("results/MA_plot.png", width = 1000, height = 800, res = 150)

plotMA(fit, coef = 2)

dev.off()


# 8. VOLCANO PLOT (mejorado)

de_results$Significant <- de_results$adj.P.Val < 0.05 &
  abs(de_results$logFC) > 1

# Crear clasificación biológica
de_results$Regulation <- "No significativo"

de_results$Regulation[
  de_results$adj.P.Val < 0.05 & de_results$logFC > 0
] <- "Sobreexpresado"

de_results$Regulation[
  de_results$adj.P.Val < 0.05 & de_results$logFC < 0
] <- "Subexpresado"

# Convertir a data.frame
volcano_df <- de_results

volcano_df$negLogFDR <- -log10(volcano_df$adj.P.Val)

# Seleccionar top genes para etiquetar
top_up <- volcano_df[
  volcano_df$Regulation == "Sobreexpresado",
][order(volcano_df$adj.P.Val), ][1:2, ]

top_down <- volcano_df[
  volcano_df$Regulation == "Subexpresado",
][order(volcano_df$adj.P.Val), ][1:2, ]

top_genes <- rbind(top_up, top_down)

png("results/Volcano_plot.png", width = 1200, height = 1000, res = 200)

ggplot(volcano_df, aes(x = logFC, y = negLogFDR, color = Regulation)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Sobreexpresado" = "#D73027",
      "Subexpresado" = "#4575B4",
      "No significativo" = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 4,
    max.overlaps = 20
  ) +
  theme_bw(base_size = 15) +
  labs(
    title = "Volcano Plot – Expresión Diferencial",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

dev.off()

# 9. HEATMAP TOP 50

top50 <- head(order(de_results$adj.P.Val), 50)

exprs_heatmap <- v$E[top50, ]

annotation_df <- data.frame(
  Tissue = rse_gene_SRP188219$sra_attribute.source_name
)

rownames(annotation_df) <- colnames(exprs_heatmap)

png("results/Heatmap_top50.png", width = 1200, height = 1000, res = 150)

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_df,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = hcl.colors(50, "RdPu")
)

dev.off()


# ============================================================

cat("Analisis EDA y DEA completado correctamente.\n")
cat("Todos los resultados fueron guardados en la carpeta 'results'.\n")
