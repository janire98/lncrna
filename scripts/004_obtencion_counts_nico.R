library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(readr)

# ESCOJO LA RUTA DE LOS COUNTS NICO PORQUE TENGO AQUÍ LOS NOMBRES DE LOS NICO

ruta <- "/Users/janire/Desktop/TFM/TFM_workflow/001_Files/Subread_NICO"
setwd(ruta)

nombres_NICO <- read_delim("nombres_NICO.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
nombres_NICO <- nombres_NICO [ ,c(1:3,12)]

nombres_NICO$ENSEMBLEID <- gsub(nombres_NICO$ENSEMBLEID [25],nombres_NICO$X12 [25],nombres_NICO$ENSEMBLEID)
nombres_NICO$ENSEMBLEID <- gsub(nombres_NICO$ENSEMBLEID [90],nombres_NICO$X12 [90],nombres_NICO$ENSEMBLEID)
nombres_NICO <- nombres_NICO [ ,1:3]


## AHORA VOY A CARGAR LA CARPETA DE LOS COUNTS NORMALIZADOS CON LA TRANSFORMACIÓN RLOG

rutaRLOG <- "/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/subread/rlog_transformacion.csv"
rlog_todos_los_genes <- read_csv(rutaRLOG)


# LES VOY A QUITAR LOS PUNTOS PORQUE SI NO, NO ME COINCIDE NI UN NOMBRE
rlog_todos_los_genes$X1 <- gsub("\\..*","",rlog_todos_los_genes$X1)
names(rlog_todos_los_genes) [1] <- "ENSEMBLEID"

# JUNTO LAS DOS TABLAS
NICO <- inner_join(rlog_todos_los_genes, nombres_NICO, by="ENSEMBLEID")

NICO <- column_to_rownames(NICO, var="Gene name")
NICO <- NICO [ ,2:31]

# SE ME HAN QUEDADO 95 GENES, QUE TIENE SENTIDO DADO QUE SI MIRO LA TABLA EXCEL, CINCO DE ELLOS NO ESTABAN EN KALLISTO

rutaGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/NICOS"
write.csv(NICO, file=paste0(rutaGuardado, "/rlog_NICO.csv"))



### COUNTS NORMALIZADAS SIN MÁS

counts_normalizadas_todos_los_genes <- read_csv("/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/subread/counts_normalizadas.csv")
counts_normalizadas_todos_los_genes$X1 <- gsub("\\..*","",counts_normalizadas_todos_los_genes$X1)
names(counts_normalizadas_todos_los_genes) [1] <- "ENSEMBLEID"

NICO <- inner_join(counts_normalizadas_todos_los_genes, nombres_NICO, by="ENSEMBLEID")
NICO <- column_to_rownames(NICO, var="Gene name")
NICO <- NICO [ ,2:31]

write.csv(NICO, file=paste0(rutaGuardado, "/counts_normalizadas_NICO.csv"))




