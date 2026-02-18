#!/usr/bin/env Rscript

# 01-descarga-datos.R
# Script para descargar los datos del proyecto SRP060355 desde recount3

# Cargar librerías
library(recount3)
library(SummarizedExperiment)
library(edgeR)
library(limma)

# Obtener proyecto
human_projects <- available_projects(organism = "human")

# Descargar datos del proyecto
project_info <- subset(
  human_projects,
  project == "SRP188219" & project_type == "data_sources"
)

rse_gene_SRP188219 <- create_rse(project_info)

rse_gene_SRP188219

assay(rse_gene_SRP188219, "counts") <- compute_read_counts(rse_gene_SRP188219)

rse_gene_SRP188219 <- expand_sra_attributes(rse_gene_SRP188219)

# Guardar datos crudos
saveRDS(create_rse(project_info), file = "raw-data/raw_rse_gene_SRP188219") # Raw RSE

saveRDS(rse_gene_SRP188219, file = "processed-data/rse_gene_SRP188219") # Processed RSE

unique(rse_gene_SRP188219$sra_attribute.source_name)
#[1] "right atrial appendage" "left atrial appendage"

colData(rse_gene_SRP188219)[,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP188219)))
]

#
#DataFrame with 400 rows and 3 columns
#           sra_attribute.Sex sra_attribute.source_name sra_attribute.tissue
#                 <character>               <character>          <character>
#SRR8715000              male    right atrial appendage     atrial appendage
#SRR8715001              male    right atrial appendage     atrial appendage
#SRR8715002              male    right atrial appendage     atrial appendage
#SRR8715003              male    right atrial appendage     atrial appendage
#SRR8715004              male    right atrial appendage     atrial appendage
#...                      ...                       ...                  ...
#SRR8714981              male     left atrial appendage     atrial appendage
#SRR8714983              male     left atrial appendage     atrial appendage
#SRR8714984              male     left atrial appendage     atrial appendage
#SRR8714985              male     left atrial appendage     atrial appendage
#SRR8714986              male     left atrial appendage     atrial appendage
#
