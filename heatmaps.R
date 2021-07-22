#### SCRIPT QUE ESCOJA LOS CLUSTERS SIGNIFICATIVOS, LOS LEA Y HAGA UN HEATMAP PARA VER SI ESTÁN DOWNREGULADOS O UPREGULADOS
library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(pheatmap)
library(ComplexHeatmap)
library(gtools)


dirClusSignificativos <- "/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/Subread_normalized_counts"
setwd(dirClusSignificativos)

files <- list.files(pattern=".pdf")
nombreClusSignificativos <- gsub("\\_.*","",files)

dirTodosClust <- "/Users/janire/Desktop/TFM/TFM_workflow/005_clusterings/Subread_50_counts_normal"
setwd(dirTodosClust)

selClust <- vector("character")

for (i in seq_along(files)) {
  selClust [i] <- list.files(pattern=nombreClusSignificativos[i])
  df <-read.csv(selClust [i])
  df <- column_to_rownames(df, var="X")
  df <- df [ , 1:30]
  assign(nombreClusSignificativos[i], df )
}


nombresHeatmaps <- mixedsort(paste0("heatmap_", nombreClusSignificativos))
Pattern1 <- mixedsort(grep("cluster",names(.GlobalEnv),value=TRUE))
Pattern1_list <- mget(Pattern1)

for (i in seq_along(Pattern1)) {
  df <- as.matrix(Pattern1_list [[i]])
  assign(nombresHeatmaps [i], pheatmap(df))
}

Pattern2 <- mixedsort(grep("heatmap_",names(.GlobalEnv),value=TRUE))
Pattern2_list <- mget(Pattern2)

setwd("/Users/janire/Desktop/TFM/TFM_workflow/007_visualizacion_resultados/heatmap_2")
for(i in seq_along(Pattern2_list)){
  p <-draw(Pattern2_list [[i]])
  nombrearchivo <- names(Pattern2_list)
  archivo <- paste0(nombrearchivo[i],"_heatmap.pdf")
  
  pdf(archivo)
  print(p, newpage = TRUE)
  dev.off()
}

ClusJuntos <- do.call("rbind", unname(Pattern1_list))
write.csv(ClusJuntos, file="/Users/janire/Desktop/TFM/TFM_workflow/007_visualizacion_resultados/ClusSignificativos_normalizados.csv")
ptodos <-pheatmap(as.matrix(ClusJuntos))
d <- draw(ptodos)
archivo <- paste0("ClusJuntos_heatmap.pdf")
pdf(archivo)
print(d, newpage=FALSE)
dev.off()

#############################################################################

# VISUALIZACION HEATMAP CON EL RLOG

# 1. Cargo los Rlog sin rescalar

rlog_LNCRNA_filtrados <- read_csv("~/Desktop/TFM/TFM_workflow/003_genes_filtered/Subread/rlog_LNCRNA_filtrados.csv")
rlog_LNCRNA_filtrados <- column_to_rownames(rlog_LNCRNA_filtrados, var="X1")
rlog_LNCRNA_filtrados <- as.data.frame(rlog_LNCRNA_filtrados %>% t() %>% scale() %>% t())
rlog_LNCRNA_filtrados <- rownames_to_column(rlog_LNCRNA_filtrados, var="X1")
setwd("/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/Subread_normalized_counts")


rlog_LNCRNA_filtrados <- read_csv("~/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv")
setwd("/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/Subread_normalized_counts")



files <- list.files(pattern=".pdf")
nombreClusSignificativos <- gsub("\\_.*","",files)
nombreClusSignificativos <- mixedsort(nombreClusSignificativos)

dirGrupos <- "/Users/janire/Desktop/TFM/TFM_workflow/005_clusterings/Subread_50_counts_normal"
setwd(dirGrupos)

selClust <- vector("character")
for (i in seq_along(files)) {
  selClust [i] <- list.files(pattern=nombreClusSignificativos[i])
  df <-read.csv(selClust [i])
  df <- column_to_rownames(df, var="X")
  df <- df [ , 1:30]
  assign(nombreClusSignificativos[i], df )
}


# Ahora me quedo con los genes que estén en los clusters significativos por aparte
Pattern1 <- mixedsort(grep("cluster",names(.GlobalEnv),value=TRUE))
Pattern1_list <- mget(Pattern1)


for (i in seq_along (Pattern1_list)){
  df <- as.data.frame(Pattern1_list [[i]])
  df <- rownames_to_column (df, var="X1")
  df <- as.data.frame(df [ ,1])
  names(df) [1] <- "X1"
  assign(paste0("rlog_",nombreClusSignificativos [i]), inner_join(rlog_LNCRNA_filtrados, df, by="X1"))
}

Pattern3 <- mixedsort(grep("rlog_clus",names(.GlobalEnv),value=TRUE))
Pattern3_list <- mget(Pattern3)

for (i in seq_along(Pattern3)) {
  df <- as.data.frame(Pattern3_list [[i]])
  df <- column_to_rownames(df, var="X1")
  assign(paste0("heatmap_",Pattern3[[i]] ), pheatmap(as.matrix(df)))
}


Pattern4 <- mixedsort(grep("heatmap_rlog_cl", names(.GlobalEnv), value=TRUE))
Pattern4_list <- mget(Pattern4)



setwd("/Users/janire/Desktop/TFM/TFM_workflow/007_visualizacion_resultados/heatmap_2_normalized counts/rlog_heatmap/Rescalado")
for(i in seq_along(Pattern4_list)){
  p <-draw(Pattern4_list [[i]])
  nombrearchivo <- names(Pattern4_list)
  archivo <- paste0(nombrearchivo[i],"_heatmap.pdf")
  
  pdf(archivo)
  print(p, newpage = TRUE)
  dev.off()
}







