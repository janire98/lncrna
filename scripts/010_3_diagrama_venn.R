tcga_filtradas <- read_csv("/Users/janire/Desktop/TFM/tcga/cuentas_rescaladas/tcga_vst_rescalado.csv")
subread_counts <- read_csv("/Users/janire/Desktop/TFM/TFM_workflow/004_genes_rescaled(1,30)/Subread/rlog_LNCRNA_rescalado.csv")
ClusSignificativos <- read_csv("/Users/janire/Desktop/TFM/TFM_workflow/007_visualizacion_resultados/rlog_genes_juntos/clusters_juntos.csv")

library(VennDiagram)


# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(tcga_filtradas$...1, ClusSignificativos$...1, subread_counts$...1),
  category.names = c("TCGA" , "genes significativos ", "NASIR" ),
  filename = 'venn_diagramm3.png',
  output=TRUE)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
setwd("/Users/janire/Desktop/TFM/tcga")

venn.diagram(
  x = list(tcga_filtradas$...1, ClusSignificativos$...1, subread_counts$...1),
  category.names = c("TCGA" , "Significativos ", "NASIR"),
  filename = '#diagrama_tcga.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)






