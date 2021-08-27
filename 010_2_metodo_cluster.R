library(tidyverse)
library(tibble)
library(readr)
library(cluster)

df <- "/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv"

Lectura_archivo <- function(df) {
  df <- read_csv(df)
  df <- column_to_rownames(df, var="...1")
  df <- df [ ,-31]
  return(df)
}

df <- Lectura_archivo(df)



# voy a medir la adecuación de la métrica y el método de unión
metricas <- c("euclidean", "manhattan",  "canberra", "binary", "minkowski")
LinkMethod <- c("single", "complete", "average", "ward.D2", "centroid", "median")


EstadisticaCoph <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("metrica_distancia", "metodo_linkage", "cophenetic"))
i <- 1


for (metrica in metricas) {
  d <- dist(df, method=metrica)
  for (metodo in LinkMethod) {
    hc <- hclust(d, method=metodo)
    estadistico1 <- cophenetic(hc)
    EstadisticaCoph [i, ] <- data.frame(metrica_distancia=metrica, metodo_linkage=metodo,cophenetic= cor(d, estadistico1) )
    i <- i+1
  }
}

EstadisticaCoph <- na.omit(EstadisticaCoph)
EstadisticoOrdenado <- EstadisticaCoph [order(EstadisticaCoph$cophenetic, decreasing=TRUE), ]
# el de elección es la distancia euclidea con el linkageMethod average.
write.csv(EstadisticoOrdenado, file="/Users/janire/Desktop/TFM/TFM_workflow/011_pruebas_clusterizaciones/coefs_copheneticos.csv")

# elección del método adecuado
library(clValid)
clmethods <- c("hierarchical","pam", "kmeans", "clara", "diana")

intern <- clValid(df,  nClust = 2:6, clMethods = clmethods,
                  validation = "internal", maxitems = 3000,
                  metric = "euclidean", method = "average")
summary(intern)
plot(intern)


setwd("/Users/janire/Desktop/TFM/TFM_workflow/011_pruebas_clusterizaciones")
sink("sumario_validacion_interna.txt")
print(summary(intern))
sink() 


### Más pruebas para la validación de los clusters
library(factoextra)
library(cluster)
# enhanced clustering para ver luego aquellos valores que están mal asociados a un cluster
res.hc <- eclust(df, "hclust", k = 2, hc_metric = "euclidean",
                 method = "average", graph = FALSE) 
fviz_silhouette(res.hc)


jpeg(filename="Silhouette_coef.jpeg")
fviz_silhouette(res.hc)
dev.off()




















# mclust for Gaussian mixture models -- not use for now
library(mclust)

mb <- Mclust(df [1:20, 1:4]  , 2)
mb$modelName

# optimal number of cluster
mb$G

head(mb$z)

# get probabilities, means, variances
summary(mb, parameters = TRUE)
plot(mb, what=c("classification"))
plot(mb, "density")

install.packages("fpc")
library("fpc")

cs <- cluster.stats(dist(df), mb$classification)
cs[c("within.cluster.ss","avg.silwidth")]


library(pvclust)
df <- t(df)
fit <- pvclust(df, method.hclust="average", method.dist="euclidean")





