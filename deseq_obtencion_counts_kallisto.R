library(tidyverse)
library(dplyr)
library(stringr)
library(tibble)
library(DESeq2)
library(SummarizedExperiment)
library("tximport")
library("readr")
library("tximportData")
## PARA LOS LNCRNA

directorio <- "/Users/janire/Desktop/TFM/TFM_workflow/001_Files/Kallisto_NASIR"
setwd(directorio)

sample <- dir(directorio)
kal_dirs <- file.path(directorio, sample)
# ---------------------------------------------------------------
###### GENERACIÓN DEL DATA FRAME CON LA METADATA
condition <- rep("tumor",30)
proyecto <- rep("NASIR", 30)
s2c <- as.data.frame(sample)
files <- file.path(kal_dirs, "abundance.tsv")
s2c <- mutate(s2c, condition=condition, path=files, proyecto=proyecto)

######### METADATA TRÁNSCRITO -- GEN
dirGTF <- '/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Metadata/gencode.v38.long_noncoding_RNAs.gtf'
gtf <- rtracklayer::import(dirGTF)
gtf.df=as.data.frame(gtf)
gtf.df <- gtf.df [ ,c(10,12,17)]
ids_completos <- data.frame(gtf.df, stringsAsFactors = FALSE)


#### Los nombres duplicados pueden generar problemas, por lo que hay que tratarlos.
tx2 <-ids_completos %>% distinct(transcript_id, .keep_all = TRUE)
tx2 <- drop_na(tx2)
names(tx2) [3] <- "TXNAME"
names(tx2) [1] <- "GENEID"

tx2 <- data.frame(TXNAME=tx2$TXNAME, GENEID=tx2$GENEID, SYMBOL=tx2$gene_name)


######### DESEQ   ##################################

###### OBTENCIÓN DE ARCHIVOS: EL .TSV
names(files) <- sample
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2, ignoreAfterBar = TRUE)
### compruebo que todo coincide
all(colnames(txi.kallisto.tsv) %in% rownames(tx2))
all(colnames(txi.kallisto.tsv) == rownames(tx2))

## DESeqDataSet
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData=s2c, design = ~1)

dds <- estimateSizeFactors(dds)
counts_normalized <- counts(dds, normalized=TRUE)
dirGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Kallisto"
write.csv(counts_normalized, file=paste0(dirGuardado, "/counts_normalizadas_LNCRNA.csv"))
dds.rlog <- rlog(dds, blind = FALSE)
rlog.norm.counts <- assay(dds.rlog)
write.csv(rlog.norm.counts, file=paste0(dirGuardado, "/rlog_LNCRNA.csv"))

counts_raw <- counts(dds, normalized=FALSE)
write.csv(counts_raw, file=paste0(dirGuardado, "/LNCRNA_raw.csv"))





































##### PARA TODOS LOS GENES

#### INICIO: ESTABLEZCO EL PATH PARA LLEGAR A LOS ARCHIVOS
setwd("/Users/janire/Desktop/TFM/TFM_workflow/001_Kallisto_files/NASIR")

sample <- dir("/Users/janire/Desktop/TFM/TFM_workflow/001_Kallisto_files/NASIR")
kal_dirs <- file.path("/Users/janire/Desktop/TFM/TFM_workflow/001_Kallisto_files/NASIR", sample)
# ---------------------------------------------------------------
###### GENERACIÓN DEL DATA FRAME CON LA METADATA
condition <- rep("tumor",30)
proyecto <- rep("NASIR", 30)
s2c <- as.data.frame(sample)
files <- file.path(kal_dirs, "abundance.tsv")
s2c <- mutate(s2c, condition=condition, path=files, proyecto=proyecto)

######### METADATA TRÁNSCRITO -- GEN
gtf <- rtracklayer::import('/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Metadata/gencode.v38.annotation.gtf')
gtf_df=as.data.frame(gtf)
gtf.df <- gtf_df [ ,c(10,12,16)]
ids_completos <- data.frame(gtf.df, stringsAsFactors = FALSE)


#### Los nombres duplicados pueden generar problemas, por lo que hay que tratarlos.
tx2 <-ids_completos %>% distinct(transcript_id, .keep_all = TRUE)
tx2 <- drop_na(tx2)
names(tx2) [3] <- "TXNAME"
names(tx2) [1] <- "GENEID"

tx2 <- data.frame(TXNAME=tx2$TXNAME, GENEID=tx2$GENEID, SYMBOL=tx2$gene_name)


######### DESEQ   ##################################

###### OBTENCIÓN DE ARCHIVOS: EL .TSV
names(files) <- sample
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2, ignoreAfterBar = TRUE)
### compruebo que todo coincide
all(colnames(txi.kallisto.tsv) %in% rownames(tx2))
all(colnames(txi.kallisto.tsv) == rownames(tx2))

## DESeqDataSet
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData=s2c, design = ~1)

dds <- estimateSizeFactors(dds)
counts_normalized <- counts(dds, normalized=TRUE)
write.csv(counts_normalized, file="/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Counts_normalizadas/counts_normalizadas_todos_los_genes.csv")
dds.rlog <- rlog(dds, blind = FALSE)
rlog.norm.counts <- assay(dds.rlog)
write.csv(rlog.norm.counts, file="/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Rlog_transformacion/rlog_todos_los_genes.csv")
counts_raw <- counts(dds, normalized=FALSE)
write.csv(counts_raw, file="/Users/janire/Desktop/TFM/TFM_workflow/002_gene_counts_obtention(deseq)/Counts_sin_normalizar/todos_los_genes_raw.csv")








