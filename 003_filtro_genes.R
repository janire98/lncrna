## FILTRO DE GENES

library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(readr)


### OBTENCIÓN DE LOS LNCRNA. LECTURA DEL GTF DE LOS LNCRNA. ME QUEDO CON NOMBRES ÚNICOS Y OMITO NAS
dirGTF <- '/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Metadata/gencode.v38.long_noncoding_RNAs.gtf'
gtf <- rtracklayer::import(dirGTF)
gtf.df=as.data.frame(gtf)
gtf.df <- gtf.df [ ,c(10,12)]
ids_completos <- data.frame(gtf.df, stringsAsFactors = FALSE)

## FILTRO POR NOMBRE DE GEN PORQUE LOS ENSG TIENEN ALGUNOS GENES IGUALES, SOLO QUE PONE ADEMÁS ("_PAR_Y"), POR EJEMPLO.
ids_unicos <- ids_completos %>% distinct(gene_name, .keep_all=TRUE)
ids_unicos <- drop_na(ids_unicos)

# Cargo el archivo de las counts de subread para filtrar por lncRNA
caminoCountsNormalizadas <- "~/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/subread/counts_normalizadas.csv"
counts_normalizadas <- read_csv(caminoCountsNormalizadas)
names(counts_normalizadas) [1] <- "gene_id"

lncrna <- merge(counts_normalizadas, ids_unicos, by="gene_id")
lncrna <- column_to_rownames(lncrna, var="gene_name")
lncrna <- lncrna [ ,2:31]

#### ahora ya tengo los lncrna de las counts normalizadas


################  1. FILTRO DE LAS NORMALIZED COUNTS PARA SABER QUÉ GENES ME QUEDO  #############

# FILTRADO DE LNCRNA NORMALIZED COUNTS
# AHORA QUE TENGO LOS LNCRNA, VOY A FILTRARLOS POR TRES CRITERIOS:
  #   1 ME QUEDO CON LOS GENES QUE TENGAN AL MENOS 10 COUNTS EN UN MÍNIMO DE 5 PACIENTES
  #   2 ME QUEDO CON LOS GENES CON UN MÁXIMO IGUAL O SUPERIOR A 50
  #   3 ME QUEDO CON LOS GENES CON UNA DESVIACIÓN ESTÁNDAR IGUAL O SUPERIOR A 10
lncrnaFiltro <- lncrna [rowSums(lncrna>=10) >=5, ]
# cuenta el número de muestras que tengan 10 counts o más. PRUEBA DE QUE VA BIEN LO QUE HE HECHO
#contar_mayores_10 <- apply(lncrna,1, function(x) sum(x >= 10))
# quédate solo con los genes que tengan en al menos cinco muestras esas 10 counts
#filtro1 <- lncrna[contar_mayores_10 >=5,]

lncrnaFiltro <- as.data.frame(lncrnaFiltro)
lncrnaFiltro$max <- apply(lncrnaFiltro [ ,1:30], 1 , FUN=max)
lncrnaFiltro$sd <- apply(lncrnaFiltro [ ,1:30], 1 , FUN=sd)

filtro_final <- lncrnaFiltro %>% filter(max>=50 & sd>=10)

lncrna_filtrados <- filtro_final [ , 1:30]

camGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/003_genes_filtered/Subread"
write.csv(lncrna_filtrados, file=paste0(camGuardado, "/LNCRNA_filtrados.csv"))


### GENES FILTRADOS PARA UTILIZAR EN EL RLOG TRANSFORMATION
KeepGenes <- data.frame(gene_name=rownames(lncrna_filtrados))


## CARGO EL RLOG DE TODOS LOS GENES Y HAGO UN MERGE O INNER JOIN PARA QUEDARME EN LOS RLOG SOLO LOS MENCIONADOS. 
# ANTES TENGO QUE PONER LOS GENE_NAMES, YA QUE SON DIFERENTES EN EL ARCHIVO YA FILTRADO Y EL QUE VOY A CARGAR

camRlogGenes <- "/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/subread/rlog_transformacion.csv"
rlog_LNCRNA <- read_csv(camRlogGenes)
names(rlog_LNCRNA) [1] <- "gene_id"
rlog_LNCRNA <- merge(rlog_LNCRNA, ids_unicos, var="gene_id")
rlog_LNCRNA <- rlog_LNCRNA[ ,2:32]


rlog_filtrados <- inner_join(rlog_LNCRNA, KeepGenes, by="gene_name")
rlog_filtrados <- column_to_rownames(rlog_filtrados, var="gene_name")

rlogGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/003_genes_filtered/Subread"
write.csv(rlog_filtrados, file= paste0(rlogGuardado, "/rlog_LNCRNA_filtrados.csv"))













