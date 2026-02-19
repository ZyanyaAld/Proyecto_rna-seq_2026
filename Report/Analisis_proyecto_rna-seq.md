# Análisis de Expresión Diferencial – SRP188219
**Zyanya Valentina Velazquez Aldrete**

## Descripción del proyecto

En este trabajo se realizó un análisis de RNA-seq para comparar la expresión génica entre:

- Apéndice auricular izquierdo
- Apéndice auricular derecho

Dataset: SRP188219  
Número total de muestras: 400 (200 por tejido)

El objetivo fue identificar genes diferencialmente expresados entre ambos tejidos cardíacos.

---

# 1. Control de calidad

Se calculó la proporción de lecturas asignadas a genes:

assigned_gene_prop = lecturas_asignadas / lecturas_totales

Resumen:

- Mínimo: 0.49
- Mediana: 0.56
- Máximo: 0.63

Todas las muestras presentaron buena calidad (> 0.3), por lo que no se eliminaron muestras.

Se generaron:
- Histograma de proporción de lecturas asignadas
- Boxplot de calidad

---

# 2. Filtrado de genes poco expresados

Se eliminaron genes con media de expresión ≤ 0.1.

Resultados:

- Genes iniciales: 63,856
- Genes retenidos: 35,699
- Porcentaje retenido: 55.91%
- Genes eliminados: 44.09%

Este filtrado reduce ruido y mejora el poder estadístico.

---

# 3. Normalización

Se utilizó:

- edgeR (normalización TMM)
- Transformación a logCPM

---

# 4. Análisis Exploratorio (EDA)

Se realizaron:

- PCA
- MDS
- Gráficas de densidad de logCPM

Resultados:

El PCA y el MDS mostraron separación entre tejidos, indicando diferencias transcriptómicas claras entre el apéndice auricular izquierdo y derecho.

---

# 5. Análisis de Expresión Diferencial

Se utilizó:

- Modelo lineal con limma
- Corrección de varianza con voom
- Moderación Bayesiana (eBayes)

Modelo:

Expresión ~ Tejido

Criterio de significancia:

FDR < 0.05

Los resultados completos se guardaron en:

results/DE_results_SRP188219.csv

---

# 6. Visualizaciones generadas

- PCA
- MDS
- MA plot
- Volcano plot
- Heatmap de los 50 genes más significativos
- Tendencia media-varianza (voom)

---

# Interpretacion Biológica

El análisis de expresión diferencial identificó a BMP10 como el gen más significativamente sobreexpresado en la aurícula derecha respecto a la izquierda (logFC = 9.00, FDR < 0.001). BMP10 es un gen fundamental en el desarrollo y mantenimiento del tejido cardíaco, particularmente involucrado en la proliferación y diferenciación de cardiomiocitos. Su mayor expresión en la aurícula derecha podría estar relacionada con diferencias funcionales y estructurales entre ambas aurículas, incluyendo variaciones en carga hemodinámica y regulación contráctil. Estos resultados sugieren que existen programas moleculares específicos que distinguen ambos compartimentos auriculares.

En cuanto al análisis de componentes principales (PCA), se mostró una tendencia de separación entre muestras de aurícula izquierda y aurícula derecha, indicando que existen diferencias globales en los perfiles de expresión génica entre ambos tejidos. Sin embargo, también se observa cierta superposición entre grupos, lo que sugiere que, aunque comparten un programa transcripcional cardíaco común, presentan diferencias específicas asociadas a su función anatómica y fisiológica.

# Conclusión

El análisis identificó un conjunto importante de genes diferencialmente expresados entre el apéndice auricular izquierdo y derecho. El análisis exploratorio confirmó que los perfiles transcriptómicos difieren entre tejidos, apoyando la relevancia biológica de los resultados obtenidos.

El flujo de trabajo sigue buenas prácticas en análisis de RNA-seq utilizando edgeR y limma.
