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
library(readxl)
library(SeuratExtend)
### Read Files in

setwd('~/Documents/CD8_Longitudinal/Annotation')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"



### TARA
load(paste0(load.path,'TARA_TCR_Combined.RData'))
load(paste0(load.path,'EARTH_TCR_Combined.RData'))


######## Rejoin layers ####
TARA_ALL[["RNA"]] <- JoinLayers(TARA_ALL[["RNA"]])
TARA_ALL[["ADT"]] <- JoinLayers(TARA_ALL[["ADT"]])

rm(TARA_ALL_TRB_0, TARA_ALL_TRB_1, TARA_ALL_TRB_2)
gc()


EARTH[["RNA"]] <- JoinLayers(EARTH[["RNA"]])
EARTH[["ADT"]] <- JoinLayers(EARTH[["ADT"]])

rm(EARTH_TRB_0, EARTH_TRB_1, EARTH_TRB_2)
gc()

# Remove CP012_1m from Tara_ALL
TARA_ALL <- subset(TARA_ALL, subset = orig.ident != "CP013_1m")

#################################### RNA and Protein Features of interest ##############################################library(dplyr)


rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B',
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18',
  'TNFRSF4','TNFRSF9','TOX','TBX21','TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1',
  'XCL1','XCL2','ZBTB16','ZEB2'
)

prots <- rownames(TARA_ALL@assays$ADT)

Idents(TARA_ALL) <- 'snn.louvianmlr_1'
Idents(EARTH) <- 'snn.louvianmlr_1'

#Remove clusters with <50 cells
# Identify clusters with fewer than 20 cells
cluster_sizes <- table(Idents(TARA_ALL))
small_clusters <- names(cluster_sizes[cluster_sizes < 100])
large_clusters <- names(cluster_sizes[cluster_sizes > 20])
TARA_ALL <- subset(TARA_ALL, idents = large_clusters)

cluster_sizes <- table(Idents(EARTH))
small_clusters <- names(cluster_sizes[cluster_sizes < 100])
large_clusters <- names(cluster_sizes[cluster_sizes > 20])
EARTH <- subset(EARTH, idents = large_clusters)

############# Plots ###############################
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Cluster_Plot')
p1 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

ggsave(
  filename = "TARA_WNN_lmlr1.png",
  plot=p1,
  width = 10,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300
)

p2 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

ggsave(
  filename = "EARTH_WNN_lmlr1.png",
  plot=p2,
  width = 10,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300 
)


############################################################ Annotated Plots ####################################################
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL')

TARA_ALL$Manual_Annotation <- NULL
# Set cluster identity to Louvain clustering
Idents(TARA_ALL) <- 'snn.louvianmlr_1'

# Load and process annotation file
TARA_annotation_df <- read_excel("TARA_ALL_Annotation.xlsx") %>%
  rename(cluster = `Cluster Number`, annotation = `Cell Type`) %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(cluster %in% as.character(0:35)) %>%
  mutate(label = paste0(cluster, ": ", annotation))

# Get cluster order from Seurat object
cluster_levels <- levels(Idents(TARA_ALL))

# Create named vector for mapping cluster -> label
TARA_cluster_map <- setNames(TARA_annotation_df$label, TARA_annotation_df$cluster)

# Map annotations to cells based on cluster identity
cluster_ids <- as.character(Idents(TARA_ALL))
TARA_ALL$Manual_Annotation <- unname(TARA_cluster_map[cluster_ids])

# Ensure factor levels follow cluster order
TARA_ALL$Manual_Annotation <- factor(TARA_ALL$Manual_Annotation, levels = TARA_cluster_map[cluster_levels])

# Set identity class to Manual_Annotation
Idents(TARA_ALL) <- "Manual_Annotation"

# UMAP Plot 1: pt.size = 1
p1 <- DimPlot2(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "Manual_Annotation",
  label = TRUE,
  box = TRUE,
  repel = TRUE,
  label.color = "black",
  cols = 'default',
  label.size = 3.5,
  pt.size = 1
) + ggtitle("TARA_ALL: Manual_Annotation")

ggsave(
  filename = "TARA_Annotated.png",
  plot = p1,
  width = 13,
  height = 8,
  dpi = 300,
  bg = "white"
)

p2 <- ClusterDistrBar(TARA_ALL$orig.ident, TARA_ALL$Manual_Annotation, cols = "default", flip = FALSE, border = "black") +
  theme(axis.title.x = element_blank())

ggsave(
  filename = "TARA_Cluster_Distribution.png",
  plot = p2,
  width = 25,
  height = 14,
  dpi = 300,
  bg = "white"
)

 saveRDS(TARA_ALL, file = file.path(load.path, "TARA_ALL_post_annotation.rds"))


#####################
#### Average Expression ##############


avg_TARA <- AverageExpression(TARA_ALL, return.seurat = FALSE)
write.csv(avg_TARA$RNA, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/avg_rna_TARA.csv")
write.csv(avg_TARA$ADT, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/avg_adt_TARA.csv")

avg_EARTH <- AverageExpression(EARTH, return.seurat = FALSE)
write.csv(avg_EARTH$RNA, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/EARTH/avg_rna_EARTH.csv")
write.csv(avg_EARTH$ADT, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/EARTH/avg_adt_EARTH.csv")


########################################## Feature Plots and VLN Plots ###############################################
### TARA
### ADT
rna.features <- intersect(rna.features, rownames(TARA_ALL[["RNA"]]))

DefaultAssay(TARA_ALL) <- "ADT"

# VLN Plots (SeuratExtend)
setwd("~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Violin_Plot/ADT")

for (i in prots) {
  if (i %in% rownames(TARA_ALL[["ADT"]])) {
    vln.pl <- VlnPlot2(
      TARA_ALL,
      features = i,
      cols='default',
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("ADT |", i))
    ggsave(paste0(i, "_VLNplot.png"), dpi = 500, width = 14, height = 8, plot = vln.pl, bg = "white")
  }
}

# Feature Plots (keep original)
setwd("~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Feature_Plot/ADT")

for (i in prots) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(
    TARA_ALL, reduction = "wnn.umap", features = i,
    colors_use = pal, order = TRUE
  )
  ggsave(paste0(i, "_Featureplot_Magma.png"), dpi = 500, width = 8, plot = fea.pl)
}

### RNA

DefaultAssay(TARA_ALL) <- "RNA"

# VLN Plots (SeuratExtend)
setwd("~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Violin_Plot/RNA")

for (i in rna.features) {
  if (i %in% rownames(TARA_ALL[["RNA"]])) {
    vln.pl <- VlnPlot2(
      TARA_ALL,
      features = i,
      cols='default',
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("RNA |", i))
    ggsave(paste0(i, "_VLNplot.png"), dpi = 500, width = 14, height = 8, plot = vln.pl, bg = "white")
  }
}

# Feature Plots (keep original)
setwd("~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Feature_Plot/RNA")

for (i in rna.features) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(
    TARA_ALL, reduction = "wnn.umap", features = i,
    colors_use = pal, order = TRUE
  )
  ggsave(paste0(i, "_Featureplot_Magma.png"), dpi = 500, width = 8, plot = fea.pl)
}

### EARTH
### ADT

DefaultAssay(EARTH) <- "ADT"

# VLN Plots (SeuratExtend)
setwd("~/Documents/CD8_Longitudinal/Annotation/EARTH/Violin_Plot/ADT")

for (i in prots) {
  if (i %in% rownames(EARTH[["ADT"]])) {
    vln.pl <- VlnPlot2(
      EARTH,
      features = i,
      cols='default',
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("ADT |", i))
    ggsave(paste0(i, "_VLNplot.png"), dpi = 500, width = 14, height = 8, plot = vln.pl, bg = "white")
  }
}

# Feature Plots (keep original)
setwd("~/Documents/CD8_Longitudinal/Annotation/EARTH/Feature_Plot/ADT")

for (i in prots) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(
    EARTH, reduction = "wnn.umap", features = i,
    colors_use = pal, order = TRUE
  )
  ggsave(paste0(i, "_Featureplot_Magma.png"), dpi = 500, width = 8, plot = fea.pl)
}

### RNA

DefaultAssay(EARTH) <- "RNA"

# VLN Plots (SeuratExtend)
setwd("~/Documents/CD8_Longitudinal/Annotation/EARTH/Violin_Plot/RNA")

for (i in rna.features) {
  if (i %in% rownames(EARTH[["RNA"]])) {
    vln.pl <- VlnPlot2(
      EARTH,
      features = i,
      cols='default',
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("RNA |", i))
    ggsave(paste0(i, "_VLNplot.png"), dpi = 500, width = 14, height = 8, plot = vln.pl, bg = "white")
  }
}

# Feature Plots (keep original)
setwd("~/Documents/CD8_Longitudinal/Annotation/EARTH/Feature_Plot/RNA")

for (i in rna.features) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(
    EARTH, reduction = "wnn.umap", features = i,
    colors_use = pal, order = TRUE
  )
  ggsave(paste0(i, "_Featureplot_Magma.png"), dpi = 500, width = 8, plot = fea.pl)
}

