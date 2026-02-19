# Proyecto RNA-seq 2026

***Zyanya Valentina Velazquez Aldrete***

## Análisis de expresión diferencial – SRP188219

Este proyecto realiza un análisis completo de RNA-seq comparando:

- Apéndice auricular izquierdo
- Apéndice auricular derecho

El objetivo fue identificar genes diferencialmente expresados entre ambos tejidos cardíacos.

El flujo de trabajo incluye:

- Control de calidad
- Filtrado de genes poco expresados
- Normalización (edgeR)
- Análisis exploratorio (PCA, MDS)
- Expresión diferencial (limma + voom)
- Visualizaciones (MA plot, Volcano, Heatmap)

## Estructura del proyecto

Code/
Scripts de análisis

results/
Resultados y figuras generadas

reporte_RNAseq_SRP188219.md
Reporte detallado del análisis

## Tecnologías utilizadas

- R
- edgeR
- limma
- ggplot2
- pheatmap
- recount3

---