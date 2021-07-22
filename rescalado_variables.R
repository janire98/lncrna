library(tidyverse)
library(readr)
library(scales)

RutaArchivo <- "/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/NICOS/counts_normalizadas_NICO.csv"
filtrados <- read_csv(RutaArchivo)
lncrna <- column_to_rownames(filtrados, var="X1")
genes <- lncrna %>% t() %>% scale() %>% t()
df <- rescale(genes, to=c(1,30))
df<- as.data.frame(df)

RutaGuardado <-"/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/NICO"
write.csv(df, file=paste0(RutaGuardado, "/NICO_counts_rescalado.csv"))








