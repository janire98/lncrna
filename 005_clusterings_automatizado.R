library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(cluster)
library(factoextra)
library(gtools)
library(scales)
library(gtools)

rutaDF <- "/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv"
df <- read_csv(rutaDF)
df <- column_to_rownames(df, var="...1")

## FUNCIÓN QUE COGE UN DATA FRAME Y LE HACE UN CLUSTERING JERÁRQUICO, DIVIDIENDO EN DOS GRUPOS Y CREANDO UNA LISTA CON
# CADA UNO DE LOS CLÚSTER
grupos <- function(df) {
  df <- as.data.frame(df)
  d <- dist(df, method="euclidean")
  hc1 <- hclust(d, method="average")
  grupo <- cutree(hc1, k=2)
  df <- df %>% mutate(cluster=grupo)
  cluster1 <- subset(df, cluster==1)
  cluster2 <- subset(df, cluster==2)
  output <- list(cluster1, cluster2)
  return(output)
}

# HAGO LA PRIMERA RONDA MANUAL PARA TENER YA UNA LISTA QUE INTRODUCIR Y CON LA QUE TRABAJAR
lista_1 <- grupos(df)
intentos <- paste0("lista_", 2:500) #nombres para mis listas
# introduzco los contadores que creo necesarios para las numerizaciones
#contador <- 1
contaCluster <-1
contaResultados <- 2
i <-1

## ESCOJO EL NÚMERO DE GENES QUE QUIERO EN CADA CLUSTER
numGenes <- as.numeric(50)
## ESCOJO EL PATH EN EL QUE GUARDAR EL ARCHIVO CLUSTER
camino <- "/Users/janire/Desktop/TFM/TFM_workflow/005_clusterings/rlog_50"


# BUCLE QUE INTERACTÚA CON LA LISTA GENERADA PARA SEGUIR HACIENDO LOS CLUSTERINGS JERÁRQUICOS Y CREAR NUEVAS LISTAS,
# YA QUE NO SE CONOCE CUÁNTAS ITERACIONES SE NECESITAN. EL BUCLE SE ROMPE UNA VEZ LA LARGURA DE LA LISTA CON LA QUE SE HA TRABAJADO ES 0
# EL BUCLE INCLUYE UNA PARTE EN LA QUE BORRA LA LISTA CREADA PARA SUSTITUIRLA POR UNA MODIFICADA QUE SOLO CONTIENE LOS DF
# QUE TIENEN QUE SEGUIR CLUSTERIZÁNDOSE
while(TRUE) {
  # me extrae del environment las variables con el patrón lista y pido que me los ordene según número
  Pattern1<-mixedsort(grep("lista",names(.GlobalEnv),value=TRUE))
  #Pattern1 <- mixedsort(Pattern1)
  Pattern1_list<-mget(Pattern1)
  assign(intentos [i], mapply(Pattern1_list [[i]], FUN=grupos)) #aquí realizo el clustering
  # vuelvo a solicitar una actualización del patrón lista
  Pattern2<-mixedsort(grep("lista",names(.GlobalEnv),value=TRUE))
  #Pattern2 <- mixedsort(Pattern2)
  Pattern2_list<-mget(Pattern2)
  #pido que vaya recorriendo cada data frame para que escriba un csv de aquellos que tienen menos de 50 genes. Esto se puede modificar
  for (j in seq_along(Pattern2_list [[contaResultados]])) {
    dfInt <- Pattern2_list [[contaResultados]] [[j]]
    if (dim(dfInt) [1] <=numGenes){
      write.csv(dfInt, file=paste0(camino,"/cluster", contaCluster, ".csv"))
      contaCluster <- contaCluster +1
    }
  }
  # pido que me filtre y que la lista del patrón se quede solo con los df de más de 50 genes
  Pattern2_list [[contaResultados]] <-Filter(function(x) dim(x)[1] >= numGenes, Pattern2_list[[contaResultados]])
  # preparo para quedarme solo con la lista con la que acabo de trabajar, para así borrarla después
  varname <- as.character(intentos[i])
  tablaVars <- ls()
  tabla_to_remove <- tablaVars[ tablaVars %in% varname]
  
  #elimino la lista que contiene todos los dataframes y la sustituyo por una con el mismo nombre pero que solo tenga los df
  # con más de 50 genes, para que así no se clustericen los menores a 50 y la líe con clusters que no deberían estar
  rm(list=tabla_to_remove)
  assign(intentos[i], Pattern2_list [[contaResultados]])
  
  
  # condición para que finalice el bucle. Si la lista que se ha trabajado ahora tiene longitud de 0,
  # es decir, no tiene data frames porque todos son menores a 50 genes, se detiene el bucle while.
  if(length(Pattern2_list [[contaResultados]])==0)  {
    break
  }
  contaResultados <- contaResultados +1
  i <- i+1
  
}
