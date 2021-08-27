library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
setwd("/Users/janire/Desktop/TFM/TFM_desarrollo")

# descargo datos cínicos
query <- GDCquery(project = "TCGA-LIHC", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")


# descargo datos de muestras
query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  #data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts"
                  )

GDCdownload(query)
data <- GDCprepare(query)
# obtengo la información de las muestras
informacion_muestra <-as.data.frame(colData(data[ , !data$shortLetterCode=="NT"]))

informacion <- as.data.frame(colData(data))
write.csv(informacion_muestra, file="/Users/janire/Desktop/TFM/tcga/datos_clinicos/informacion_muestra.csv")

counts <- as.data.frame(assay(data [ ,!data$shortLetterCode=="NT"], "HTSeq - Counts"))

datos_summarizedexperiment <-data[ , !data$shortLetterCode=="NT"]


library(DESeq2)
ddsSE <- DESeqDataSet(datos_summarizedexperiment, design = ~1)
ddsSE <- estimateSizeFactors(ddsSE)
counts_normalized <- as.data.frame(counts(ddsSE, normalized=TRUE))
write.csv(counts_normalized, file="/Users/janire/Desktop/TFM/tcga/cuentas_normalizadas/counts_normalizadas.csv")

vst_transformacion <- vst(ddsSE, blind=FALSE)
vst_counts <- as.data.frame(assay(vst_transformacion))
write.csv(vst_counts, file="/Users/janire/Desktop/TFM/tcga/cuentas_normalizadas/vst_transformacion.csv")


#### --------------------------------------------- FILTRADO DE GENES

dirGTF <- '/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Metadata/gencode.v38.long_noncoding_RNAs.gtf'
gtf <- rtracklayer::import(dirGTF)
gtf.df=as.data.frame(gtf)
gtf.df <- gtf.df [ ,c(10,12)]
ids_completos <- data.frame(gtf.df, stringsAsFactors = FALSE)

## FILTRO POR NOMBRE DE GEN PORQUE LOS ENSG TIENEN ALGUNOS GENES IGUALES, SOLO QUE PONE ADEMÁS ("_PAR_Y"), POR EJEMPLO.
ids_unicos <- ids_completos %>% distinct(gene_name, .keep_all=TRUE)
ids_unicos <- drop_na(ids_unicos)
ids_unicos$gene_id <- gsub("\\..*","",ids_unicos$gene_id)

counts_normalizadas <- read_csv("~/Desktop/TFM/tcga/cuentas_normalizadas/counts_normalizadas.csv")
names(counts_normalizadas) [1] <- "gene_id"

lncrna <- merge(counts_normalizadas, ids_unicos, by="gene_id")
lncrna <- column_to_rownames(lncrna, var="gene_name")
lncrna <- lncrna [ ,-1]

lncrnaFiltro <- as.data.frame(lncrna [rowSums(lncrna>=10) >=62, ])
lncrnaFiltro$max <- apply(lncrnaFiltro [ ,1:374], 1 , FUN=max)
lncrnaFiltro$sd <- apply(lncrnaFiltro [ ,1:374], 1 , FUN=sd)

filtro_final <- lncrnaFiltro %>% filter(max>=50 & sd>=10)

lncrna_filtrados <- filtro_final [ , -c(375,376)]
write.csv(lncrna_filtrados, file="/Users/janire/Desktop/TFM/tcga/cuentas_filtradas/counts_filtradas.csv")

vst_transformacion <- read_csv("~/Desktop/TFM/tcga/cuentas_normalizadas/vst_transformacion.csv")
names(vst_transformacion) [1] <- "gene_id"
vst_LNCRNA <- merge(vst_transformacion, ids_unicos, var="gene_id")
vst_LNCRNA <- vst_LNCRNA[ ,-1]
KeepGenes <- data.frame(gene_name=rownames(lncrna_filtrados))

vst_filtrados <- inner_join(vst_LNCRNA, KeepGenes, by="gene_name")
vst_filtrados <- column_to_rownames(vst_filtrados, var="gene_name")
write.csv(vst_filtrados, file="/Users/janire/Desktop/TFM/tcga/cuentas_filtradas/vst_filtrados.csv")

##### -----------------------------------

#### ESCALADO DE LAS VARIABLES
library(scales)
RutaArchivo <- "/Users/janire/Desktop/TFM/tcga/cuentas_filtradas/vst_filtrados.csv"
filtrados <- read_csv(RutaArchivo)
lncrna <- column_to_rownames(filtrados, var="...1")
genes <- lncrna %>% t() %>% scale() %>% t()
df <- as.data.frame(rescale(genes, to=c(1,30)))

RutaGuardado <-"/Users/janire/Desktop/TFM/tcga/cuentas_rescaladas"
write.csv(df, file=paste0(RutaGuardado, "/30_vst.csv"))

################### CLUSTERINGS
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

rutaDF <- "/Users/janire/Desktop/TFM/tcga/cuentas_rescaladas/30_vst.csv"
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
camino <- "/Users/janire/Desktop/TFM/tcga/30rango_clusterings"


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


####### ANÁLISIS DE SUPERVIVENCIA

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(dplyr)
library(readr)


tcga_counts_rescalado <- read_csv("/Users/janire/Desktop/TFM/tcga/cuentas_rescaladas/30_vst.csv")
tcga_counts_rescalado <- column_to_rownames(tcga_counts_rescalado, var="...1")
colnames(tcga_counts_rescalado) <- substr(colnames(tcga_counts_rescalado), 1, 12)

clinical <- read_csv("/Users/janire/Desktop/TFM/tcga/datos_clinicos/clinical.csv")
clinical <- clinical [ ,-1]


datos_clinicos <- clinical [ ,c(1,23,24,25,26)]
informacion_filtrada <- informacion_muestra %>% filter(!shortLetterCode=="NT")

datos_clinicos <- datos_clinicos [datos_clinicos$bcr_patient_barcode %in% informacion_filtrada$patient, ]
datos_clinicos <- datos_clinicos %>% distinct(bcr_patient_barcode, .keep_all=TRUE)


informacion_filtrada <- informacion_filtrada %>% distinct(patient, .keep_all=TRUE)

datos_clinicos$OS <- as.integer(
  ifelse( is.na(datos_clinicos$days_to_death),
          datos_clinicos$days_to_last_followup,
          datos_clinicos$days_to_death))

write.csv(datos_clinicos, file="/Users/janire/Desktop/TFM/tcga/datos_clinicos/datos_clinicos_final_30.csv")

surv_object <- Surv(time = datos_clinicos$OS)

rutaArchivo <- "/Users/janire/Desktop/TFM/tcga/30rango_clusterings"
files <- list.files(rutaArchivo, pattern=".csv", recursive=TRUE, full.names=TRUE)

filtroMinimoPacientes <- as.numeric(74)

rutaGuardado <- "/Users/janire/Desktop/TFM/tcga/analisis_supervivencia_30rango"


# Import t_rog dataset rescaled and make the hierarchical clustering
for (file in files) {
  tryCatch ({
    datos_rescalado <- read_csv(file)
  datos_rescalado <- datos_rescalado [ ,-376]
  datos_rescalado <- tibble::column_to_rownames(datos_rescalado, var="...1")
  #colnames(datos_rescalado) <- substr(colnames(datos_rescalado), 1, 12)
  datos_final <- as.data.frame(t(datos_rescalado ))
  datos_final <- rownames_to_column(datos_final, var="...1")
  #datos_final$...1 <- gsub('\\.', '-', datos_final$...1)
  for (i in seq_along(datos_clinicos$bcr_patient_barcode)) {
    #datos_final$...1 <-str_replace(datos_final$...1, i, i)
    indice <-as.numeric(grep(datos_clinicos$bcr_patient_barcode [i], datos_final$...1))
    
    datos_final$...1 [indice] <- datos_clinicos$bcr_patient_barcode [i]
  }
  datos_final <- datos_final %>% distinct(...1, .keep_all=TRUE)
  datos_final <- column_to_rownames(datos_final, var="...1")
  
  #rlog <- as.data.frame(rlog [ -24, ])
  
  res.dist <- dist(datos_final, method = "euclidean")
  hc1 <- hclust(res.dist, method = "average")
  clust <- cutree(hc1, k = 2)
  
  
  fit1 <- survfit(surv_object ~ clust, data=datos_clinicos)
  numPacientesMinimo <- as.numeric(min(fit1 [['strata']]))
  
  
  
  
  pvalor <-surv_pvalue(fit1)
  
  if (pvalor [2]<=0.05 & numPacientesMinimo>=filtroMinimoPacientes) {
    name <- gsub(".*/", "", file)
    cluster <- gsub("\\..*","",name)
    survp <-ggsurvplot(fit1,
                       pval = TRUE, conf.int = FALSE,
                       title =name,
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       linetype = "strata", # Change line type by groups
                       surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_bw(), # Change ggplot2 theme
    )
    
    setwd(rutaGuardado)
    archivo <- paste0(cluster,"_surv.pdf")
    
    pdf(archivo)
    print(survp, newpage = FALSE)
    dev.off()
    
  } else{
    print(paste0("El cluster del archivo ", file, " no tiene un p-valor significativo, por lo que estadísticamente no se puede relacionar el cluster con una mayor supervivencia"))
  }
  }, error=function(e){})
}




informacion_muestra
view(datos_clinicos)



prueba <- datos_final %>% distinct(...1, .keep_all=TRUE)













