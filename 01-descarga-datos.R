#!/usr/bin/env Rscript

# 01-descarga-datos.R
# Script para descargar los datos del proyecto SRP060355 desde recount3

# Cargar librerías
library("recount3")

# Crear carpeta para datos crudos si no existe
if (!dir.exists("raw-data")) {
  dir.create("raw-data")
}

# Descargar datos del proyecto
cat("Buscando proyecto SRP060355 en recount3...\n")
human_projects <- available_projects(organism = "human")

project_info <- subset(
  human_projects,
  project == "SRP060355" & project_type == "data_sources"
)

if (nrow(project_info) == 0) {
  stop("No se encontró el proyecto SRP060355")
}

cat("Proyecto encontrado. Descargando datos...\n")
rse_gene_SRP060355 <- create_rse(project_info)

# Guardar datos crudos
saveRDS(rse_gene_SRP060355, file = "raw-data/rse_gene_SRP060355_raw.rds")

cat("Datos guardados en raw-data/rse_gene_SRP060355_raw.rds\n")
cat("Dimensiones del objeto:", dim(rse_gene_SRP060355), "\n")
