library(tidyverse)
library(readr)



LNCRNA_filtrados <- read_csv("~/Desktop/TFM/TFM_workflow/003_genes_filtered/Subread/LNCRNA_filtrados.csv")
LNCRNA_filtrados <- column_to_rownames(LNCRNA_filtrados, var="...1")
orden <- gtools::mixedsort(colnames(LNCRNA_filtrados))

LNCRNA_filtrados <- (LNCRNA_filtrados) [ , orden ]

rlog_LNCRNA_rescalado <- read_csv("~/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv")
rlog_LNCRNA_rescalado <- column_to_rownames(rlog_LNCRNA_rescalado, var="...1")
rlog_LNCRNA_rescalado <- rlog_LNCRNA_rescalado [ ,orden]

counts_rescalado <- read_csv("~/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/subread_counts_rescalado.csv")
counts_rescalado <- column_to_rownames(counts_rescalado, var="...1")
counts_rescalado <- counts_rescalado [ ,orden]


library(moments)

listado_skewness <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("muestra", "counts_filtradas","counts_rescalados" ,"rlog_rescalado"))

for (i in seq_along(counts_rescalado)) {
  muestra <- names(counts_rescalado) [i]
  listado_skewness [i,1] <- muestra
  listado_skewness [i ,2] <- skewness(LNCRNA_filtrados [ ,i])
  listado_skewness [i, 3] <- skewness(counts_rescalado [ ,i])
  listado_skewness [i, 4] <- skewness(rlog_LNCRNA_rescalado [ ,i])
}

write.csv(listado_skewness, file="/Users/janire/Desktop/TFM/TFM_workflow/012_visualizacion_datos/asimetrias.csv")





