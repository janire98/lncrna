library(tidyverse)
library(dplyr)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(scales)

datos_supervivencia <- read_delim("/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/datos_supervivencia.csv", 
                                  ";", escape_double = FALSE, trim_ws = TRUE)

datos <- datos_supervivencia [ ,c(10,12:14)]
names(datos) [1] <- "X1"
datos$X1 <-paste0("RNA", datos$X1)
surv_object <- Surv(time = datos$OS_days_68)


rutaArchivo <- "/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/NICO/rlog_NICO_rescalado.csv"
rlog_NICO <- read_csv(rutaArchivo)


df <- column_to_rownames(rlog_NICO, var="X1")

rutaGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/NICO"


for (i in seq_along(as.data.frame(t(df)))) {
  gen <- df [i, ]
  nombre <- rownames(gen)
  t_gen <- as.data.frame(t(gen))
  t_gen <- rownames_to_column(t_gen, var="X1")
  names(t_gen) [2] <- "value"
  cuantiles <- sapply(t_gen[-1], function(x) x[x <= quantile(x, 0.25) | x >= quantile(x, 0.75)])
  cuantiles <- as.data.frame(cuantiles)
  names(cuantiles) [1] <- "value"
  filtrados <- inner_join(t_gen, cuantiles, by= "value")
  filtrados$X1 <- gsub(".*_", "", filtrados$X1)
  
  datos_gen <- inner_join(datos, filtrados, by="X1")
  datos_gen <- datos_gen [ ,-5]
  
  surv_object_gen <- Surv(time = datos_gen$OS_days_68)
  
  filtrados <- column_to_rownames(filtrados, var="X1")

  res.dist <- dist(filtrados, method = "euclidean")
  hc1 <- hclust(res.dist, method = "average")
  clust <- cutree(hc1, k = 2)
  
  fit1 <- survfit(surv_object_gen ~ clust, data=datos_gen)
  pvalor <-surv_pvalue(fit1)
  
  if (pvalor [2]<=0.05) {
    
    survp <-ggsurvplot(fit1,
                       pval = TRUE, conf.int = FALSE,
                       title =nombre,
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       linetype = "strata", # Change line type by groups
                       surv.median.line = "hv", # Specify median survival
                       ggtheme = theme_bw(), # Change ggplot2 theme
    )
    
    setwd(rutaGuardado)
    archivo <- paste0(nombre,"_surv.pdf")
    
    pdf(archivo)
    print(survp, newpage = FALSE)
    dev.off()
    
  } else{
    print(paste0("El gen ", nombre, " no tiene un p-valor significativo,", pvalor,  " por lo que estadísticamente no se puede relacionar el cluster con una mayor supervivencia"))
  }
}













