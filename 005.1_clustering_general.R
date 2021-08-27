### SCRIPT DE CLUSTERING GENERAL -- seleccionar el rango de clusters a los que se quiere forzar a la máquina

library(tidyverse)
library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(factoextra)
library(cluster)

## Seleccionar en rutaArchivoApertura el archivo que se quiere abrir. Se preparan los datos
rutaArchivoApertura <- "/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv"
lncrna <- read_csv(rutaArchivoApertura)
lncrna <- column_to_rownames(lncrna, var="X1")

## Seleccionar el rango de clusters que se quiere hacer
DesdeClust <-2
HastaClust <- 50
cuantosClusters <- DesdeClust:HastaClust

## Función para realizar un clustering jerárquico al data frame seleccionado. El método de distancia escogido es el euclideo,
# y el de clustering, el average. El resultado es una lista que contiene data frames separados por cluster.
FunClusterizacion <- function(df) {
  df <- as.data.frame(df)
  d <- dist(df, method="euclidean")
  hc1 <- hclust(d, method="average")
  grupo <- cutree(hc1, k=i)
  df <- df %>% mutate(cluster=grupo)
  lista <- split(df, df$cluster)
  return(lista)
}

# Se carga la ruta a la que se quiere ir para crear las carpetas y guardar los clusters
RutaClusters <- "/Users/janire/Desktop/TFM/TFM_workflow/005_clusterings/cluster_general"
setwd(RutaClusters)

# Bucle que crea las carpetas que van a contener un número i de clusters, se realiza la función de clusterización
# y se guarda cada uno de los clusters en formato .csv dentro de la carpeta correspondiente.
for (i in cuantosClusters) {
  dir.create(paste0(RutaClusters, "/numClusters_",i))
  setwd(paste0(RutaClusters, "/numClusters_",i))
  lista <- FunClusterizacion(lncrna)
  
  for(j in 1:i){
    cluster <- as.data.frame(lista [[j]])
    cluster <- cluster [ ,1:30]
    write.csv(cluster, file=paste0("cluster",j, ".csv"))
  }
  
}





