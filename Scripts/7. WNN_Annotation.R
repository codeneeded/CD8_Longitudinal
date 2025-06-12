#Load Required Libraries
library(Seurat)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scCustomize)
library(Polychrome)
library(viridis)

### Read Files in

setwd('~/Documents/CD8_Longitudinal/Annotation')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"



### TARA
load(paste0(load.path,'TARA_TCR_Combined.RData'))


######## Rejoin layers ####
TARA_ALL[["RNA"]] <- JoinLayers(TARA_ALL[["RNA"]])
TARA_ALL[["ADT"]] <- JoinLayers(TARA_ALL[["ADT"]])

rm(TARA_ALL_TRB_0, TARA_ALL_TRB_1, TARA_ALL_TRB_2)

#################################### RNA and Protein Features of interest ##############################################

# create multimodal heatmap 
rna.features <-  c('CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
                   'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
                   'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
                   'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
                   'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
                   'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
                   'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
                   'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67')

prots <- rownames(TARA_ALL@assays$ADT)

Idents(TARA_ALL) <- 'snn.louvianmlr_1'
Idents(EARTH) <- 'snn.louvianmlr_1'

########################################## Feature Plots and VLN Plots ###############################################
### ADT

DefaultAssay(TARA_ALL) <-'ADT'

# VLN Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Violin_Plot/ADT')

for (i in prots) {
  
  vln.pl <-VlnPlot_scCustom(TARA_ALL,features = i)
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 14, vln.pl)
}

# Feature Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Feature_Plot/ADT')

for (i in prots) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(TARA_ALL, reduction = 'wnn.umap', features = i, 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_Magma.png'),dpi=500, width = 8, fea.pl)
}

### RNA

DefaultAssay(TARA_ALL) <-'RNA'

# VLN Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Violin_Plot/RNA')

for (i in rna.features) {
  
  vln.pl <-VlnPlot_scCustom(TARA_ALL,features = i)
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 14, vln.pl)
}

# Feature Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Feature_Plot/RNA')

for (i in rna.features) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(TARA_ALL, reduction = 'wnn.umap', features = i, 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_Magma.png'),dpi=500, width = 8, fea.pl)
}

