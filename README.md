# Análisis de expresión diferencial - Proyecto SRP060355

## Descripción
Análisis de expresión diferencial de genes farmacogenéticos en tejido hepático y renal usando datos de recount3 (SRP060355).

## Autora

*Zyanya Valentina Velazquez Aldrete*

## Datos
- Proyecto: SRP060355
- Título: Transcriptomic variation of pharmacogenes in multiple human tissues
- Muestras: 94 (24 hígado, 20 riñón, 25 corazón, 25 adiposo)

## Estructura del proyecto

```bash
.
├── 01-descarga-datos.R        # Script para descargar datos
├── 02-procesamiento-datos.R   # Script para procesar y filtrar
├── 03-analisis-expresion.R    # Script para análisis DE
├── 04-graficas.R              # Script para visualizaciones
├── analisis.Rmd               # Reporte final en R Markdown
├── raw-data/                  # Datos crudos (no versionados)
├── processed-data/            # Datos procesados
├── output/                    # Resultados (tablas, figuras)
└── figs/                      # Gráficas para el reporte
```
