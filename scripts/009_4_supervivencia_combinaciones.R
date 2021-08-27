### Script de análisis de supervivencia de la combinación de genes

library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(dplyr)
library(readr)

datos_supervivencia <- read_delim("/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/datos_supervivencia.csv", 
                                  ";", escape_double = FALSE, trim_ws = TRUE)

datos <- datos_supervivencia [ ,c(10,12:14)]

surv_object <- Surv(time = datos$OS_days_68)

rutaArchivo <- "/Users/janire/Desktop/TFM/TFM_workflow/008_combinacion_genes/rlog"
files <- list.files(rutaArchivo, pattern=".csv", recursive=TRUE, full.names=TRUE)

filtroMinimoPacientes <- as.numeric(6)

rutaGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/009_supervivencia_combinaciones/rlog_001"

df_comb_genes <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("cluster", "gen1", "gen2", "p-valor", "num_pacientes1", "num_pacientes2"))
i <-1

# Import t_rog dataset rescaled and make the hierarchical clustering
for (file in files) {
  Comb <- read_csv(file)
  Comb <- tibble::column_to_rownames(Comb, var="...1")
  tComb <- as.data.frame(t(Comb))
  
  #rlog <- as.data.frame(rlog [ -24, ])
  tComb <- data.frame (tComb [ c(1, 22, 26:30, 2:21, 23:25), ])
  res.dist <- dist(tComb, method = "euclidean")
  hc1 <- hclust(res.dist, method = "average")
  clust <- cutree(hc1, k = 2)
  
  fit1 <- survfit(surv_object ~ clust, data=datos)
  numPacientesMinimo <- as.numeric(min(fit1 [['strata']]))
  numPacientesMaximo <- as.numeric(max(fit1 [['strata']]))
  
  
  
  pvalor <-surv_pvalue(fit1)
  
  if (pvalor [2]<=0.01 & numPacientesMinimo>=filtroMinimoPacientes) {
    name <- gsub(".*/", "", file)
    cluster <- gsub("\\..*","",name)
    nombGenes <- colnames(tComb)
    nombCluster <- gsub("\\_.*","",cluster)
    
    #survp <-ggsurvplot(fit1,
                       #pval = TRUE, conf.int = FALSE,
                       #title =paste0(name,"  ->  ",nombGenes[1],"  &  ", nombGenes[2]),
                       #risk.table = TRUE, # Add risk table
                       #risk.table.col = "strata", # Change risk table color by groups
                       #linetype = "strata", # Change line type by groups
                       #surv.median.line = "hv", # Specify median survival
                       #ggtheme = theme_bw(), # Change ggplot2 theme
    #)
    
    #setwd(rutaGuardado)
    #archivo <- paste0(cluster,"_surv.pdf")
    
    #pdf(archivo)
    #print(survp, newpage = FALSE)
    #dev.off()
    
    df_comb_genes [i ,1] <- nombCluster
    df_comb_genes [i, 2] <- nombGenes [1]
    df_comb_genes [i, 3] <- nombGenes [2]
    df_comb_genes [i, 4] <- pvalor [2]
    df_comb_genes [i, 5] <- numPacientesMinimo
    df_comb_genes [i, 6] <- numPacientesMaximo
    i <- i+1
    
  } else{
    print(paste0("El cluster del archivo ", file, " no tiene un p-valor significativo, por lo que estadísticamente no se puede relacionar el cluster con una mayor supervivencia"))
  }
}


write.csv(df_comb_genes, file="/Users/janire/Desktop/TFM/TFM_workflow/010_resultados_combinaciones/rlog001/genes_comb.csv")
df_separados <- split(df_comb_genes, df_comb_genes$cluster)

for (grupo in seq_along(df_separados)) {
  nombre <- names(df_separados [grupo])
  df <-as.data.frame(df_separados [grupo])
  df <- df [ ,-1]
  write.csv(df, file=paste0("/Users/janire/Desktop/TFM/TFM_workflow/010_resultados_combinaciones/rlog001/", nombre, ".csv"))
}


