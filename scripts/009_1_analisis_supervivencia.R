# ANÁLISIS DE SUPERVIVENCIA

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

rutaArchivo <- "/Users/janire/Desktop/TFM/TFM_workflow/005_clusterings/rlog_50"
files <- list.files(rutaArchivo, pattern=".csv", recursive=TRUE, full.names=TRUE)

filtroMinimoPacientes <- as.numeric(6)

rutaGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/rlog"


# Import t_rog dataset rescaled and make the hierarchical clustering
for (file in files) {
rlog_rescalado <- read_csv(file)
rlog_rescalado <- tibble::column_to_rownames(rlog_rescalado, var="...1")
rlog <- t(rlog_rescalado )

#rlog <- as.data.frame(rlog [ -24, ])
rlog <- data.frame (rlog [ c(1, 22, 26:30, 2:21, 23:25), ])
res.dist <- dist(rlog, method = "euclidean")
hc1 <- hclust(res.dist, method = "average")
clust <- cutree(hc1, k = 2)

fit1 <- survfit(surv_object ~ clust, data=datos)
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
}











