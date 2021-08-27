library(tidyverse)
library(readr)
library(tibble)
library(dplyr)

## LECTURA DEL ARCHIVO DE SUBREAD Y HTSEQ UNA VEZ RESCALADO
tcga_counts_rescalado <- read_csv("/Users/janire/Desktop/TFM/tcga/cuentas_rescaladas/tcga_vst_rescalado.csv")
subread_counts_rescalado <- read_csv("/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv")
comunes <- tcga_counts_rescalado [tcga_counts_rescalado$...1 %in% subread_counts_rescalado$...1, ] # comparten 1427 genes una vez pasado el filtro

## ahora voy a cargar los datasets que me han salido significativos en subread y veo cuántos de esos genes están

# pongo el directorio donde están los clusters significativos, para así extraer sus nombres y seleccionarlos en el montón de clusters
dirClusSignificativos <- "/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/rlog"
setwd(dirClusSignificativos)

archivosSig <- list.files()
nombreArchivosSig <- gsub("\\_.*","",archivosSig)
nombreArchivosSig <- mixedsort(nombreArchivosSig)

# directorio que contiene todos los clusters. Voy a buscar los que tienen mi patrón y los abro
dirTodosClusters <- "/Users/janire/Desktop/TFM/TFM_workflow/005_clusterings/rlog_50"
setwd(dirTodosClusters)

lecturaArchivos <- function(archivo) {
  df <- list.files(pattern=archivo)
  df <- read.csv(df)
  df <- df [ , 1:31]
  df <- column_to_rownames(df, var="X")
  return(df)
}

for (archivo in nombreArchivosSig) {
  assign(archivo, lecturaArchivos(archivo))
}


## Ahora hago las posibles combinaciones de mis genes para cada cluster (no repetidas). 
# Si fuesen repetidas, saldrían muchas más


Pattern1 <- mixedsort(grep("cluster", names(.GlobalEnv), value = TRUE))
Pattern1_list <- mget(Pattern1)

# junto en una tabla todos los genes significativos
genes_significativos_juntos <- do.call("rbind", unname(Pattern1_list))
genes_significativos_juntos <- rownames_to_column(genes_significativos_juntos, var="...1")

genes_significativos_juntos <- column_to_rownames(genes_significativos_juntos, var="...1")
write.csv(genes_significativos_juntos, file="/Users/janire/Desktop/TFM/TFM_workflow/007_visualizacion_resultados/rlog_genes_juntos/clusters_juntos.csv")

# Comparo la de tcga con los significativos

Significativos_comunes <- tcga_counts_rescalado [tcga_counts_rescalado$...1 %in% genes_significativos_juntos$...1, ]

# de los 164 genes significativos de subread, 45 se encuentran en TCGA

## voy a cargar el archivo con la transformación svt
vst_LNCRNA_filtrados <- read_csv("~/Desktop/TFM/TFM_desarrollo/003_genes_filtered/tcga/vst_LNCRNA_filtrados.csv")

vst_comunes <- tcga_counts_rescalado [tcga_counts_rescalado$...1 %in% genes_significativos_juntos$...1, ]
vst_comunes <- column_to_rownames(vst_comunes, var="...1")

library(pheatmap)

z_score <- as.data.frame(t(scale(t(vst_comunes))))

pheatmap(as.matrix(z_score), show_colnames = FALSE)

seleccion <- order(rowMeans(vst_comunes), decreasing=TRUE)

prueba <- vst_comunes [seleccion, ]
prueba2 <- as.data.frame(t(scale(t(prueba))))

pheatmap(prueba2, show_colnames = FALSE)







