library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(DESeq2)
library(SummarizedExperiment)
library("tximport")
library("readr")
library("tximportData")

#### INICIO: ESTABLEZCO EL PATH PARA LLEGAR A LOS ARCHIVOS

camino <- "/Users/janire/Desktop/TFM/TFM_workflow/001_Files/Subread_NASIR"
setwd(camino)


## GENERO LA LISTA DE ARCHIVOS Y LES ELIMINO TODA LA SOBREINFORMACIÓN DE LOS NOMBRES A LAS MUESTRAS
temp <- list.files()
nombre_muestra <- gsub("^[^_]*_","",temp)
nombre_muestra <- gsub("\\_.*","",nombre_muestra)

# GENERO UN BUCLE PARA LEER TODOS LOS ARCHIVOS Y PONER SOLO LOS NOMBRES DE LAS MUESTRAS. ES MÁS MANEJABLE
for (i in 1:length(temp)){
  assign(nombre_muestra[i], read_delim(temp[i], "\t", escape_double = FALSE, trim_ws = TRUE, 
                             skip = 1))
} 

## ESTA ES IGUAL DE GUAY PARA IMPORTAR, SOLO QUE LO DE LOS NOMBRES DE LAS MUESTRAS ESTÁ MÁS COMPLICADO
prueba <- lapply(temp, function (x) {read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE, 
                                  skip = 1)})


### TENGO QUE GENERAR UNA SOLA TABLA CON TODAS LAS COUNTS DE LAS MUESTRAS

solo_counts_df <- function (df) {
  df <- df [ ,c(1,7)]
  return(df)
}

Pattern1<-grep("RNAT",names(.GlobalEnv),value=TRUE)
Pattern1_list<-do.call("list",mget(Pattern1))

for (i in seq_along(Pattern1_list)) {
  Pattern1_list [[i]] <- solo_counts_df(Pattern1_list [[i]])
  names(Pattern1_list [[i]]) [1] <- "Gene_id"
  names(Pattern1_list [[i]]) [2] <- Pattern1 [[i]]
}

### UNO TODOS LOS DATA FRAMES DE LA LISTA SEGÚN EL GENE_ID
df_final <- Reduce(merge, Pattern1_list)
# no utilizo ahora la función para quitar los puntos de después del nombre porque hay unos genes que tienen _PAR_Y
#df_final$Gene_id <- gsub("\\..*","",df_final$Gene_id)
df_final <- column_to_rownames(df_final, var="Gene_id")
# Para comprobar que los nombres son únicos
prueba <- df_final %>% distinct(Gene_id, .keep_all=TRUE)

camGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/subread"
write.csv(df_final, file=paste0(camGuardado, "/counts_brutas.csv"))


###### -----------------------------------------------#### EMPIEZO DESEQ

sample <- dir(camino)
sample <- gsub("^[^_]*_","",sample)
sample <- gsub("\\_.*","",sample)
kal_dirs <- file.path(camino, sample)
# ---------------------------------------------------------------
###### GENERACIÓN DEL DATA FRAME CON LA METADATA
condition <- rep("tumor",30)
proyecto <- rep("NASIR", 30)
s2c <- as.data.frame(sample)
files <- file.path(kal_dirs)
s2c <- mutate(s2c, condition=condition, path=files, proyecto=proyecto)


# CARGO EL ARCHIVO DE MIS COUNTS
camCountsRaw <- paste0(camGuardado, "/counts_brutas.csv")
cts <- read_csv(camCountsRaw)
cts <- column_to_rownames(cts, var="X1")
rownames(s2c) <- s2c$sample
all(rownames(s2c) %in% colnames(cts))

all(rownames(s2c) == colnames(cts))
cts <- cts[,rownames(s2c)]
all(rownames(s2c) == colnames(cts))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = s2c,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_normalized <- counts(dds, normalized=TRUE)
write.csv(counts_normalized, file=paste0(camGuardado, "/counts_normalizadas.csv"))
dds.rlog <- rlog(dds, blind = FALSE)
rlog.norm.counts <- assay(dds.rlog)
write.csv(rlog.norm.counts, file=paste0(camGuardado, "/rlog_transformacion.csv"))

##############################################################################3

