library(tidyverse)
library(igraph)
library(ggraph)
library(stringr)
library(readr)


archivos <- list.files(path="/Users/janire/Desktop/TFM/TFM_workflow/010_resultados_combinaciones/rlog001", pattern="cluster")

for (archivo in archivos) {
  setwd("/Users/janire/Desktop/TFM/TFM_workflow/010_resultados_combinaciones/rlog001")
  nombrearchivo <- gsub("\\..*", "", archivo)
  #combinaciones$Cluster <- gsub("\\_.*","",combinaciones$Cluster)
  combinaciones <- read_csv(archivo)
  data <- combinaciones [ ,-1]
  
  g <- graph.data.frame(data, directed = FALSE) 
  
  #plot(g)
  
  #plot(g, layout = layout.grid(g))
  
  jpeg(paste0(nombrearchivo, "_grafo.jpeg"), quality = 75)
  
  
  g%>%
    ggraph() +
    geom_edge_link() +
    geom_node_label(aes(label = name)) +
    theme_graph()
  
  dev.off()
  
  
  ### Quiero ver cuántas veces se encuentra cada gen
  
  # apilo los datos
  apilamiento <- stack(data)
  apilamiento <- apilamiento [ ,1]
  
  # obtengo los nombres de los genes involucrados
  nombres_unicos <- unique(apilamiento)
  i <-1
  frecuencia_nombre <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("gen", "num_relaciones"))
  
  
  for (nombre in nombres_unicos) {
    frecuencia_nombre [i , 1] <- paste0(nombre)
    frecuencia_nombre [i , 2]  <- sum(str_count(data, pattern=nombre))
    i <- i+1
  }
  
  orden <- order(frecuencia_nombre$num_relaciones, decreasing=TRUE)
  
  frecuencia_nombre <- frecuencia_nombre [ orden, ]
  
  write.csv(frecuencia_nombre, file=paste0("/Users/janire/Desktop/TFM/TFM_workflow/010_resultados_combinaciones/rlog001/numero_relaciones_gen_",nombrearchivo, ".csv"))
  
  
}



