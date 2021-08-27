## Análisis de supervivencia gen a gen de los clusters significativos

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
genes <- read_csv(rutaArchivo)


df <- column_to_rownames(genes, var="...1")

rutaGuardado <- "/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/subread_gen_a_gen"

listado_genes <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("nombre_gen"))
contador <- 1

for (i in seq_along(as.data.frame(t(df)))) {
  gen <- df [i, ]
  nombre <- rownames(gen)
  t_gen <- as.data.frame(t(gen))
  t_gen <- rownames_to_column(t_gen, var="X1")
  names(t_gen) [2] <- "value"
  
  t_gen$X1 <- gsub('\\.', '-', t_gen$X1)
  
  datos_gen <- inner_join(datos, t_gen, by="X1")
  datos_gen <- datos_gen [ ,-5]
  
  surv_object_gen <- Surv(time = datos_gen$OS_days_68)
  
  t_gen <- column_to_rownames(t_gen, var="X1")
  
  res.dist <- dist(t_gen, method = "euclidean")
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
    
    listado_genes [contador, 1] <- nombre
    contador <- contador +1
    
  } else{
    print(paste0("El gen ", nombre, " no tiene un p-valor significativo",pvalor [2], "por lo que estadísticamente no se puede relacionar el cluster con una mayor supervivencia"))
  }
}

write.csv(listado_genes, file="/Users/janire/Desktop/TFM/TFM_workflow/006_survival_analysis/Subread/subread_gen_a_gen/listado_genes_significativos.csv")







