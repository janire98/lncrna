### SCRIPT QUE HACE COMBINACIONES DE GENES

library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(gtools)

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

caminoGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/008_combinacion_genes/rlog"

for (j in seq_along(Pattern1_list)){
  df <- as.data.frame(Pattern1_list [[j]])
  combinaciones <- as.data.frame(combinat::combn(rownames(df), 2))
  cluster <- names(Pattern1_list [j])
  for (i in seq_along(combinaciones)) {
    CombCluster <- df [combinaciones [ ,i], ]
    write.csv(CombCluster, file=paste0(caminoGuardado,"/",cluster,"_", i, ".csv"))
    print(paste0("Se ha guardado el archivo", cluster, ",que pertenece a la combinación", i ))
  }
}

# nrow(unique(t(combinaciones))) == nrow(t(combinaciones))
# combinaciones <- t(combinaciones)











