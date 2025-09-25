library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(SeuratExtend)
library(scCustomize)
# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "TARA_ALL_post_annotation.rds")

TARA_ALL <- readRDS(rds_in)
DefaultAssay(TARA_ALL) <- 'RNA'
assay_use   <- DefaultAssay(TARA_ALL)




################ Cluster Comparison Heatmaps ##################################

### Function
make_heatmap <- function(seu, clusters, outfile, n_top_genes = 10) {
  # Subset
  seu_sub <- subset(seu, subset = Manual_Annotation %in% clusters)
  
  # Get variable features
  genes <- VariableFeatures(seu_sub)
  
  # Calculate z-scores
  toplot <- CalcStats(seu_sub, features = genes, method = "zscore", order = "p", n = n_top_genes)
  
  # Heatmap
  p <- Heatmap(toplot, lab_fill = "zscore")
  
  # Save
  ggsave(outfile, p, width = 10, height = 6)
}

make_heatmap_adt <- function(seu, clusters, outfile, n_top_genes = 10) {
  # Subset
  seu_sub <- subset(seu, subset = Manual_Annotation %in% clusters)
  
  # Use ADT assay
  DefaultAssay(seu_sub) <- "ADT"
  
  # Use all proteins
  genes <- rownames(seu_sub[["ADT"]])
  
  # Calculate z-scores
  toplot <- CalcStats(
    seu_sub,
    features = genes,
    method   = "zscore",
    order    = "p",
    n        = n_top_genes
  )
  
  # Heatmap
  p <- Heatmap(toplot, lab_fill = "zscore")
  
  # Save
  ggsave(outfile, p, width = 10, height = 6)
}
### RNA

DefaultAssay(TARA_ALL) <-'RNA'

setwd ('/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Cluster_Comparison_Heatmaps/RNA')

# DN T Cells
make_heatmap(
  TARA_ALL,
  clusters = c("18: DN T cell_1","23: DN T cell_2","26: DN T cell_3","28: DN T cell_4"),
  outfile = "DN_T_Cells_Heatmap_TARA_ALL.png"
)

# CD8 
make_heatmap(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_All_Heatmap_TARA_ALL.png"
)

make_heatmap(
  TARA_ALL,
  clusters = c("8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_Only_Heatmap_TARA_ALL.png"
)

make_heatmap(
  TARA_ALL,
  clusters = c("6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_and_Naive_Heatmap_TARA_ALL.png"
)

make_heatmap(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell"),
  outfile = "CD8_MemoryvsNaive_Heatmap_TARA_ALL.png"
)

# NK ALL
make_heatmap(
  TARA_ALL,
  clusters = c("3: IL2RB+ NK cell","15: NK cell_1","17: NK cell_2","25: CD56bright NK"),
  outfile = "NK_All_Heatmap_TARA_ALL.png"
)

# B cells
make_heatmap(
  TARA_ALL,
  clusters = c("4: IgD+IgM+ B cell","11: Naïve B cell","13: Transitional B cell","22: CD69+ B cell","31: B cell","32: Plasmablast","35: TNF+ B cell"),
  outfile = "B_Cells_All_Heatmap_TARA_ALL.png"
)



### ADT
DefaultAssay(TARA_ALL) <-'ADT'


setwd ('/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Cluster_Comparison_Heatmaps/ADT')

# DN T Cells
make_heatmap_adt(
  TARA_ALL,
  clusters = c("18: DN T cell_1","23: DN T cell_2","26: DN T cell_3","28: DN T cell_4"),
  outfile = "DN_T_Cells_Heatmap_TARA_ALL.png"
)

# CD8 
make_heatmap_adt(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_All_Heatmap_TARA_ALL.png"
)

make_heatmap_adt(
  TARA_ALL,
  clusters = c("8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_Only_Heatmap_TARA_ALL.png"
)

make_heatmap_adt(
  TARA_ALL,
  clusters = c("6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_and_Naive_Heatmap_TARA_ALL.png"
)

make_heatmap_adt(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell"),
  outfile = "CD8_MemoryvsNaive_Heatmap_TARA_ALL.png"
)

# NK ALL
make_heatmap_adt(
  TARA_ALL,
  clusters = c("3: IL2RB+ NK cell","15: NK cell_1","17: NK cell_2","25: CD56bright NK"),
  outfile = "NK_All_Heatmap_TARA_ALL.png"
)

# B cells
make_heatmap_adt(
  TARA_ALL,
  clusters = c("4: IgD+IgM+ B cell","11: Naïve B cell","13: Transitional B cell","22: CD69+ B cell","31: B cell","32: Plasmablast","35: TNF+ B cell"),
  outfile = "B_Cells_All_Heatmap_TARA_ALL.png"
)

######################### NK Cell FasL Analysis #################################




####################### Find FasL+ and - Cells within the cluster and then deferentiate them ###########
DimPlot2(TARA_ALL, features = c("FAS", "FASLG"),reduction = 'wnn.umap')
FeaturePlot3(TARA_ALL, color = "ryb", feature.1 = "FAS", feature.2 = "FASLG", pt.size = 0.5,reduction = 'wnn.umap')
VlnPlot2(TARA_ALL,'FASLG',show.mean = T)


### Set viral load status
TARA_entry <- subset(TARA_ALL, subset = Age %in% c(1,2) & Condition == "HEI")
TARA_entry$Viral_Load_Category <- ifelse(TARA_entry$Viral_Load >= 100000, "High", "Low")
VlnPlot2(TARA_entry,'FAS',show.mean = T,split.by = 'Viral_Load_Category',
         stat.method = "wilcox.test")
VlnPlot2(TARA_entry,'FASLG',show.mean = T,split.by = 'Viral_Load_Category',
         stat.method = "wilcox.test")

ClusterDistrPlot(
  origin = TARA_ALL$orig.ident,
  cluster = TARA_ALL$Manual_Annotation,
  condition = TARA_ALL$Condition
)

ClusterDistrPlot(
  origin = TARA_entry$orig.ident,
  cluster = TARA_entry$Manual_Annotation,
  condition = TARA_entry$Viral_Load_Category
)

ClusterDistrBar(origin = TARA_entry$Viral_Load_Category, cluster = TARA_entry$orig.ident, flip = FALSE, reverse_order = FALSE)

levels(TARA_ALL$Manual_Annotation)




# Select some variable features
# Subset to clusters 8, 9, and 27


genes <- VariableFeatures(TARA_subset)[1:20]
DotPlot2(TARA_ALL, features = genes)

# Create grouped features
grouped_features <- list(
  "B_cell_markers" = c("MS4A1", "CD79A"),
  "T_cell_markers" = c("CD3D", "CD8A", "IL7R"),
  "Myeloid_markers" = c("CD14", "FCGR3A", "S100A8")
)

DotPlot2(TARA_ALL, features = grouped_features)
# Using colors instead of borders for split groups
DotPlot2(TARA_ALL, 
         features = genes, 
         group.by = "Manual_Annotation", 
         split.by = "Condition", 
         split.by.method = "color", 
         show_grid = FALSE)
