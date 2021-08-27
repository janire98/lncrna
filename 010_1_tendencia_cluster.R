library(clustertend)
library(tidyverse)
library(readr)
library(stringr)
library(tibble)
library(sleuth)

## Cargo mi dataset de cuentas normalizadas y escaladas
camino <- "/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv"



  df <- read_csv(camino)
  df <- column_to_rownames(df, var="...1")
  df <- t(df) # la traspongo porque necesita que las muestras estén como features
  # cálculo del estadístico de Hopkins para los datos reales
  set.seed(1)
  estadistico1 <-hopkins(df, n = nrow(df)-1)
  # genero datos aleatorios y calculo su estadístico
  random_df <- apply(df, 2,
                     function(x){runif(length(x), min(x), (max(x)))})
  random_df <- as.data.frame(random_df)
  
  set.seed(1)
  estadistico2 <-hopkins(scale(random_df), n = nrow(random_df)-1)
  estadisticofinal <- data.frame(datos_reales=estadistico1, datos_aleatorios=estadistico2)
  names(estadisticofinal) [1] <- "datos_reales"
  names(estadisticofinal) [2] <- "datos_aleatorios"
  
  

setwd("/Users/janire/Desktop/TFM/TFM_workflow/011_pruebas_clusterizaciones")

 
  jpeg("HCC_Data.jpeg")
  fviz_dist(dist(df), show_labels = FALSE)+ labs(title = "HCC data")
  dev.off()

  
  jpeg("random_Data.jpeg")
  fviz_dist(dist(scale(random_df)), show_labels = FALSE)+ labs(title = "random data")
  dev.off()
  
  









