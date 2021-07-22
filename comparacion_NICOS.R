### SCRIPT PARA VER SI HAY NICOS EN LOS CLUSTERS SIGNIFICATIVOS
library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(pheatmap)
library(readr)


## 1. HACER HEATMAP DE LOS NICOS

NICO <- read_csv("~/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/NICO/rlog_NICO_rescalado.csv")
NICO <- column_to_rownames(NICO, var="X1")

pheatmap(as.matrix(NICO))

nombresNico <- as.data.frame(rownames(NICO))
names(nombresNico) [1] <- "X1"

ClusJuntos <- read_csv("/Users/janire/Desktop/TFM/TFM_workflow/007_visualizacion_resultados/heatmap_2_normalized counts/ClusSignificativos_normalizados.csv")

nicos_significativos <- merge(ClusJuntos, nombresNico, by="X1")





library(DESeq2)
counts_normalizadas_NICO <- read_csv("~/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/NICOS/counts_normalizadas_NICO.csv")
counts_normalizadas_NICO <- column_to_rownames(counts_normalizadas_NICO, var="X1")

df <- counts_normalizadas_NICO %>% t() %>% scale() %>% t()
df <- as.data.frame(rescale(df, to=c(1,30)))
pheatmap(as.matrix(df))





