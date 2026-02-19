#!/usr/bin/env Rscript

# 02-Procesamiento de Datos.R
# Script procesar y limpiar los datos

# 1. CARGA DE LIBRERÍAS

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# 2. CARGA DE DATOS PROCESADOS

rse_gene_SRP188219 <- readRDS("processed-data/rse_gene_SRP188219")

# Exploramos dimensiones
dim(rse_gene_SRP188219)

# 3. PREPARACIÓN DE VARIABLES

# Convertimos el tejido a factor
colData(rse_gene_SRP188219)$sra_attribute.source_name <-
  factor(colData(rse_gene_SRP188219)$sra_attribute.source_name)

# Verificamos niveles
levels(rse_gene_SRP188219$sra_attribute.source_name)
# "left atrial appendage" "right atrial appendage"

table(colData(rse_gene_SRP188219)$sra_attribute.source_name)
# left atrial appendage right atrial appendage
#                  200                    200

# 4. CONTROL DE CALIDAD
# Proporción de lecturas asignadas a genes
rse_gene_SRP188219$assigned_gene_prop <-
  rse_gene_SRP188219$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP188219$recount_qc.gene_fc_count_all.total

summary(rse_gene_SRP188219$assigned_gene_prop)

# Visualización
# Crear carpeta results si no existe
if (!dir.exists("results")) {
  dir.create("results")
}

# Guardar histograma
png(
  "results/assigned_gene_prop_hist.png",
  width = 1000,
  height = 800,
  res = 150
)

hist(
  rse_gene_SRP188219$assigned_gene_prop,
  breaks = 20,
  col = "#8DA0CB",
  border = "white",
  main = "Proportion of Reads Assigned to Genes",
  xlab = "Assigned proportion",
  ylab = "Frequency"
)

dev.off()

png(
  "results/assigned_gene_prop_boxplot.png",
  width = 1000,
  height = 800,
  res = 150
)

boxplot(
  rse_gene_SRP188219$assigned_gene_prop,
  col = "#E78AC3",
  border = "grey30",
  main = "Assigned Gene Proportion",
  ylab = "Proportion"
)
dev.off()

# Se puede filtrar por baja calidad, sin embargo nestro minimo fue de 49, asi que no vamos a quitar nada al filtrar
#Filtrar muestras de baja calidad
#rse_gene_SRP188219_filtrado <-
#  rse_gene_SRP188219[, rse_gene_SRP188219$assigned_gene_prop > 0.3]

# 5. FILTRADO DE GENES POCO EXPRESADOS

gene_means <- rowMeans(assay(rse_gene_SRP188219, "counts"))

# Guardamos copia sin filtrar
rse_gene_SRP188219_unfiltered <- rse_gene_SRP188219

# Creamos objeto NUEVO filtrado (no sobreescribimos el original)
rse_gene_SRP188219_filtered <-
  rse_gene_SRP188219[gene_means > 0.1, ]

# Porcentaje retenido
round(
  nrow(rse_gene_SRP188219_filtered) /
    nrow(rse_gene_SRP188219_unfiltered) *
    100,
  2
)

# Guardamos objeto filtrado
saveRDS(
  rse_gene_SRP188219_filtered,
  file = "processed-data/rse_gene_SRP188219_filtered.rds"
)
