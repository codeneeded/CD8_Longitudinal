#Load Required Libraries
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(batchelor)
library(harmony)
library(reticulate)
##set path to load data


setwd('~/Documents/CD8_Longitudinal/Integration')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"

load(paste0(load.path,'Seuratv5_isotype_Assay3.RData'))


################ Subset Data as needed ##############
seurat_isotype[["RNA"]] <- JoinLayers(seurat_isotype[["RNA"]])

# Subset for Cohort = TARA and Condition = HEI
TARA_HEI <- subset(seurat_isotype, subset = Cohort == "TARA" & Condition == "HEI")

# Subset for Cohort = TARA (all conditions)
TARA_ALL <- subset(seurat_isotype, subset = Cohort == "TARA")

# Subset for Cohort = EARTH (all conditions)
EARTH <- subset(seurat_isotype, subset = Cohort == "EARTH")

rm(seurat_isotype)
gc()


################################################# Integrate RNA #############################################

### Split Layers for RNA for each batch
TARA_ALL[["RNA"]] <- split(TARA_ALL[["RNA"]], f = TARA_ALL$orig.ident)
TARA_HEI[["RNA"]] <- split(TARA_HEI[["RNA"]], f = TARA_HEI$orig.ident)
EARTH[["RNA"]] <- split(EARTH[["RNA"]], f = EARTH$orig.ident)


# Standard Processing TARA_ALL
TARA_ALL <- NormalizeData(TARA_ALL)
TARA_ALL <- FindVariableFeatures(TARA_ALL)
TARA_ALL <- ScaleData(TARA_ALL)
TARA_ALL <- RunPCA(TARA_ALL)
TARA_ALL <- FindNeighbors(TARA_ALL, dims = 1:30, reduction = "pca")
TARA_ALL <- FindClusters(TARA_ALL, resolution = 2, cluster.name = "unintegrated_clusters")
TARA_ALL <- RunUMAP(TARA_ALL, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Standard Processing TARA_HEI
TARA_HEI <- NormalizeData(TARA_HEI)
TARA_HEI <- FindVariableFeatures(TARA_HEI)
TARA_HEI <- ScaleData(TARA_HEI)
TARA_HEI <- RunPCA(TARA_HEI)
TARA_HEI <- FindNeighbors(TARA_HEI, dims = 1:30, reduction = "pca")
TARA_HEI <- FindClusters(TARA_HEI, resolution = 2, cluster.name = "unintegrated_clusters")
TARA_HEI <- RunUMAP(TARA_HEI, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Standard Processing EARTH
EARTH <- NormalizeData(EARTH)
EARTH <- FindVariableFeatures(EARTH)
EARTH <- ScaleData(EARTH)
EARTH <- RunPCA(EARTH)
EARTH <- FindNeighbors(EARTH, dims = 1:30, reduction = "pca")
EARTH <- FindClusters(EARTH, resolution = 2, cluster.name = "unintegrated_clusters")
EARTH <- RunUMAP(EARTH, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

### TARA ALL


TARA_ALL <- IntegrateLayers(
  TARA_ALL,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca.rna"
)

TARA_ALL <- IntegrateLayers(
  TARA_ALL,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn.rna"
)

### TARA HEI
TARA_HEI <- IntegrateLayers(
  TARA_HEI,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn.rna"
)

TARA_HEI <- IntegrateLayers(
  TARA_HEI,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca.rna"
)

### EARTH
EARTH <- IntegrateLayers(
  EARTH,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn.rna"
)


EARTH <- IntegrateLayers(
  EARTH,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca.rna"
)


#### UMAP and Clustering #########
### TARA All
TARA_ALL <- FindNeighbors(TARA_ALL, reduction = "integrated.cca.rna", dims = 1:30)
TARA_ALL <- FindClusters(TARA_ALL, resolution = 2, cluster.name = "cca_clusters_rna")
TARA_ALL <- RunUMAP(TARA_ALL, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna")
TARA_ALL <- FindNeighbors(TARA_ALL, reduction = "integrated.mnn.rna", dims = 1:30)
TARA_ALL <- FindClusters(TARA_ALL, resolution = 2, cluster.name = "mnn_clusters_rna")
TARA_ALL <- RunUMAP(TARA_ALL, reduction = "integrated.mnn.rna", dims = 1:30, reduction.name = "umap.mnn.rna")

### TARA HEI
TARA_HEI <- FindNeighbors(TARA_HEI, reduction = "integrated.cca.rna", dims = 1:30)
TARA_HEI <- FindClusters(TARA_HEI, resolution = 2, cluster.name = "cca_clusters_rna")
TARA_HEI <- RunUMAP(TARA_HEI, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna")
TARA_HEI <- FindNeighbors(TARA_HEI, reduction = "integrated.mnn.rna", dims = 1:30)
TARA_HEI <- FindClusters(TARA_HEI, resolution = 2, cluster.name = "mnn_clusters_rna")
TARA_HEI <- RunUMAP(TARA_HEI, reduction = "integrated.mnn.rna", dims = 1:30, reduction.name = "umap.mnn.rna")

### EARTH
EARTH <- FindNeighbors(EARTH, reduction = "integrated.cca.rna", dims = 1:30)
EARTH <- FindClusters(EARTH, resolution = 2, cluster.name = "cca_clusters_rna")
EARTH <- RunUMAP(EARTH, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna")
EARTH <- FindNeighbors(EARTH, reduction = "integrated.mnn.rna", dims = 1:30)
EARTH <- FindClusters(EARTH, resolution = 2, cluster.name = "mnn_clusters_rna")
EARTH <- RunUMAP(EARTH, reduction = "integrated.mnn.rna", dims = 1:30, reduction.name = "umap.mnn.rna")

### PLOT

### TARA ALL
p1 <- DimPlot(
  TARA_ALL,
  reduction = "umap.cca.rna",
  group.by = c("predicted.celltype.l2", "cca_clusters_rna"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  TARA_ALL,
  reduction = "umap.mnn.rna",
  group.by = c("predicted.celltype.l2", "mnn_clusters_rna"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "TARA_ALL_RNA_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

### TARA HEI

p1 <- DimPlot(
  TARA_HEI,
  reduction = "umap.cca.rna",
  group.by = c("predicted.celltype.l2", "cca_clusters_rna"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  TARA_HEI,
  reduction = "umap.mnn.rna",
  group.by = c("predicted.celltype.l2", "mnn_clusters_rna"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "TARA_HEI_RNA_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

### EARTH
p1 <- DimPlot(
  EARTH,
  reduction = "umap.cca.rna",
  group.by = c("predicted.celltype.l2", "cca_clusters_rna"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  EARTH,
  reduction = "umap.mnn.rna",
  group.by = c("predicted.celltype.l2", "mnn_clusters_rna"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "EARTH_RNA_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

############ ADT Integration #####################

DefaultAssay(TARA_ALL) <-'ADT'
DefaultAssay(TARA_HEI) <-'ADT'
DefaultAssay(EARTH) <-'ADT'

### Convert assay to v5
TARA_ALL[["ADT"]] <- as(object = TARA_ALL[["ADT"]], Class = "Assay5")
TARA_HEI[["ADT"]] <- as(object = TARA_HEI[["ADT"]], Class = "Assay5")
EARTH[["ADT"]] <- as(object = EARTH[["ADT"]], Class = "Assay5")

### Split Layers for ADT for each batch
TARA_ALL[["ADT"]] <- split(TARA_ALL[["ADT"]], f = TARA_ALL$orig.ident)
TARA_HEI[["ADT"]] <- split(TARA_HEI[["ADT"]], f = TARA_HEI$orig.ident)
EARTH[["ADT"]] <- split(EARTH[["ADT"]], f = EARTH$orig.ident)


# define proteins to use in clustering (non-isptype controls)
prots <- rownames(TARA_ALL[["ADT"]]$data)
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b', 'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')
prots <- setdiff(prots, isotype_genes)
VariableFeatures(TARA_ALL) <- prots
VariableFeatures(TARA_HEI) <- prots
VariableFeatures(EARTH) <- prots


# Standard Processing TARA_ALL
TARA_ALL <- ScaleData(TARA_ALL)
TARA_ALL <- RunPCA(TARA_ALL,reduction.name = "apca")
TARA_ALL <- FindNeighbors(TARA_ALL, dims = 1:30, reduction = "apca")
TARA_ALL <- FindClusters(TARA_ALL, resolution = 2, cluster.name = "unintegrated_clusters_adt")
TARA_ALL <- RunUMAP(TARA_ALL, dims = 1:30, reduction = "apca", reduction.name = "umap.unintegrated.adt")

# Standard Processing TARA_HEI
TARA_HEI <- ScaleData(TARA_HEI)
TARA_HEI <- RunPCA(TARA_HEI, reduction.name = "apca")
TARA_HEI <- FindNeighbors(TARA_HEI, dims = 1:30, reduction = "apca")
TARA_HEI <- FindClusters(TARA_HEI, resolution = 2, cluster.name = "unintegrated_clusters_adt")
TARA_HEI <- RunUMAP(TARA_HEI, dims = 1:30, reduction = "apca", reduction.name = "umap.unintegrated.adt")


# Standard Processing EARTH
EARTH <- ScaleData(EARTH)
EARTH <- RunPCA(EARTH, reduction.name = "apca")
EARTH <- FindNeighbors(EARTH, dims = 1:30, reduction = "apca")
EARTH <- FindClusters(EARTH, resolution = 2, cluster.name = "unintegrated_clusters_adt")
EARTH <- RunUMAP(EARTH, dims = 1:30, reduction = "apca", reduction.name = "umap.unintegrated.adt")

### TARA ALL
TARA_ALL <- IntegrateLayers(
  TARA_ALL,
  orig.reduction = 'apca',
  features = prots,
  method= FastMNNIntegration,
  assay = 'ADT',
  new.reduction = "integrated.mnn.adt"
)

TARA_ALL <- IntegrateLayers(
  TARA_ALL,
  method= CCAIntegration,
  orig.reduction = "apca",
  features = prots,
  assay = 'ADT',
  new.reduction = "integrated.cca.adt"
)

### TARA HEI
TARA_HEI <- IntegrateLayers(
  TARA_HEI,
  orig.reduction = 'apca',
  features = prots,
  method= FastMNNIntegration,
  assay = 'ADT',
  new.reduction = "integrated.mnn.adt"
)

TARA_HEI<- IntegrateLayers(
  TARA_HEI,
  method= CCAIntegration,
  orig.reduction = "apca",
  features = prots,
  assay = 'ADT',
  new.reduction = "integrated.cca.adt"
)

### EARTH
EARTH <- IntegrateLayers(
  EARTH,
  orig.reduction = 'apca',
  features = prots,
  method= FastMNNIntegration,
  assay = 'ADT',
  new.reduction = "integrated.mnn.adt"
)

EARTH <- IntegrateLayers(
  EARTH,
  method= CCAIntegration,
  orig.reduction = "apca",
  features = prots,
  assay = 'ADT',
  new.reduction = "integrated.cca.adt"
)

#### UMAP and Clustering #########
### TARA All
TARA_ALL <- FindNeighbors(TARA_ALL, reduction = "integrated.cca.adt", dims = 1:30)
TARA_ALL <- FindClusters(TARA_ALL, resolution = 2, cluster.name = "cca_clusters_adt")
TARA_ALL <- RunUMAP(TARA_ALL, reduction = "integrated.cca.adt", dims = 1:30, reduction.name = "umap.cca.adt")
TARA_ALL <- FindNeighbors(TARA_ALL, reduction = "integrated.mnn.adt", dims = 1:30)
TARA_ALL <- FindClusters(TARA_ALL, resolution = 2, cluster.name = "mnn_clusters_adt")
TARA_ALL <- RunUMAP(TARA_ALL, reduction = "integrated.mnn.adt", dims = 1:30, reduction.name = "umap.mnn.adt")

### TARA HEI
TARA_HEI <- FindNeighbors(TARA_HEI, reduction = "integrated.cca.adt", dims = 1:30)
TARA_HEI <- FindClusters(TARA_HEI, resolution = 2, cluster.name = "cca_clusters_adt")
TARA_HEI <- RunUMAP(TARA_HEI, reduction = "integrated.cca.adt", dims = 1:30, reduction.name = "umap.cca.adt")
TARA_HEI <- FindNeighbors(TARA_HEI, reduction = "integrated.mnn.adt", dims = 1:30)
TARA_HEI <- FindClusters(TARA_HEI, resolution = 2, cluster.name = "mnn_clusters_adt")
TARA_HEI <- RunUMAP(TARA_HEI, reduction = "integrated.mnn.adt", dims = 1:30, reduction.name = "umap.mnn.adt")

### EARTH
EARTH <- FindNeighbors(EARTH, reduction = "integrated.cca.adt", dims = 1:30)
EARTH <- FindClusters(EARTH, resolution = 2, cluster.name = "cca_clusters_adt")
EARTH <- RunUMAP(EARTH, reduction = "integrated.cca.adt", dims = 1:30, reduction.name = "umap.cca.adt")
EARTH <- FindNeighbors(EARTH, reduction = "integrated.mnn.adt", dims = 1:30)
EARTH <- FindClusters(EARTH, resolution = 2, cluster.name = "mnn_clusters_adt")
EARTH <- RunUMAP(EARTH, reduction = "integrated.mnn.adt", dims = 1:30, reduction.name = "umap.mnn.adt")


### PLOT

### TARA ALL
p1 <- DimPlot(
  TARA_ALL,
  reduction = "umap.cca.adt",
  group.by = c("predicted.celltype.l2", "cca_clusters_adt"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  TARA_ALL,
  reduction = "umap.mnn.adt",
  group.by = c("predicted.celltype.l2", "mnn_clusters_adt"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "TARA_ALL_ADT_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

### TARA HEI

p1 <- DimPlot(
  TARA_HEI,
  reduction = "umap.cca.adt",
  group.by = c("predicted.celltype.l2", "cca_clusters_adt"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  TARA_HEI,
  reduction = "umap.mnn.adt",
  group.by = c("predicted.celltype.l2", "mnn_clusters_adt"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "TARA_HEI_adt_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

### EARTH
p1 <- DimPlot(
  EARTH,
  reduction = "umap.cca.adt",
  group.by = c("predicted.celltype.l2", "cca_clusters_adt"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  EARTH,
  reduction = "umap.mnn.adt",
  group.by = c("predicted.celltype.l2", "mnn_clusters_adt"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "EARTH_adt_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

#### Save to load path

save(TARA_ALL, file = paste0(load.path, "TARA_ALL_RNA_ADT.Rdata"))
save(TARA_HEI, file = paste0(load.path, "TARA_HEI_RNA_ADT.Rdata"))
save(EARTH, file = paste0(load.path, "EARTH_RNA_ADT.Rdata"))

############# WSNN #######################################

TARA_ALL <- FindMultiModalNeighbors(
  TARA_ALL, reduction.list = list("integrated.mnn.rna", "integrated.cca.adt"), 
  dims.list = list(1:30, 1:14), modality.weight.name = "wnn.weight"
)

TARA_HEI <- FindMultiModalNeighbors(
  TARA_HEI, reduction.list = list("integrated.mnn.rna", "integrated.cca.adt"), 
  dims.list = list(1:30, 1:14), modality.weight.name = "wnn.weight"
)

EARTH <- FindMultiModalNeighbors(
  EARTH, reduction.list = list("integrated.mnn.rna", "integrated.cca.adt"), 
  dims.list = list(1:30, 1:14), modality.weight.name = "wnn.weight"
)

################ UMAP ##########################################

### TARA ALL

TARA_ALL <- RunUMAP(TARA_ALL, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

TARA_ALL <- FindClusters(TARA_ALL, graph.name = "wsnn", algorithm = 2, resolution = 0.5, cluster.name = 'snn.louvianmlr_0.5')
TARA_ALL <- FindClusters(TARA_ALL, graph.name = "wsnn", algorithm = 2, resolution = 1, cluster.name = 'snn.louvianmlr_1')
TARA_ALL <- FindClusters(TARA_ALL, graph.name = "wsnn", algorithm = 2, resolution = 1.5, cluster.name = 'snn.louvianmlr_1.5')

TARA_ALL <- FindClusters(TARA_ALL, graph.name = "wsnn", algorithm = 3, resolution = 0.5, cluster.name = 'snn.slm_0.5')
TARA_ALL <- FindClusters(TARA_ALL, graph.name = "wsnn", algorithm = 3, resolution = 1, cluster.name = 'snn.slm_1')
TARA_ALL <- FindClusters(TARA_ALL, graph.name = "wsnn", algorithm = 3, resolution = 1.5, cluster.name = 'snn.slm_1.5')

### TARA ALL

TARA_HEI <- RunUMAP(TARA_HEI, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

TARA_HEI <- FindClusters(TARA_HEI, graph.name = "wsnn", algorithm = 2, resolution = 0.5, cluster.name = 'snn.louvianmlr_0.5')
TARA_HEI <- FindClusters(TARA_HEI, graph.name = "wsnn", algorithm = 2, resolution = 1, cluster.name = 'snn.louvianmlr_1')
TARA_HEI <- FindClusters(TARA_HEI, graph.name = "wsnn", algorithm = 2, resolution = 1.5, cluster.name = 'snn.louvianmlr_1.5')

TARA_HEI <- FindClusters(TARA_HEI, graph.name = "wsnn", algorithm = 3, resolution = 0.5, cluster.name = 'snn.slm_0.5')
TARA_HEI <- FindClusters(TARA_HEI, graph.name = "wsnn", algorithm = 3, resolution = 1, cluster.name = 'snn.slm_1')
TARA_HEI <- FindClusters(TARA_HEI, graph.name = "wsnn", algorithm = 3, resolution = 1.5, cluster.name = 'snn.slm_1.5')

### TARA ALL

EARTH <- RunUMAP(EARTH, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

EARTH <- FindClusters(EARTH, graph.name = "wsnn", algorithm = 2, resolution = 0.5, cluster.name = 'snn.louvianmlr_0.5')
EARTH <- FindClusters(EARTH, graph.name = "wsnn", algorithm = 2, resolution = 1, cluster.name = 'snn.louvianmlr_1')
EARTH <- FindClusters(EARTH, graph.name = "wsnn", algorithm = 2, resolution = 1.5, cluster.name = 'snn.louvianmlr_1.5')

EARTH <- FindClusters(EARTH, graph.name = "wsnn", algorithm = 3, resolution = 0.5, cluster.name = 'snn.slm_0.5')
EARTH <- FindClusters(EARTH, graph.name = "wsnn", algorithm = 3, resolution = 1, cluster.name = 'snn.slm_1')
EARTH <- FindClusters(EARTH, graph.name = "wsnn", algorithm = 3, resolution = 1.5, cluster.name = 'snn.slm_1.5')

### Convert cluster to numeric/factors
# Convert cluster identities to factors with numeric sorting
# Convert cluster identities to factors with numeric sorting
for (cluster_col in c("snn.louvianmlr_0.5", "snn.louvianmlr_1", "snn.louvianmlr_1.5",
                      "snn.slm_0.5", "snn.slm_1", "snn.slm_1.5")) {
  TARA_ALL[[cluster_col]][, 1] <- factor(
    TARA_ALL[[cluster_col]][, 1],
    levels = sort(as.numeric(levels(TARA_ALL[[cluster_col]][, 1])))
  )
}
for (cluster_col in c("snn.louvianmlr_0.5", "snn.louvianmlr_1", "snn.louvianmlr_1.5",
                      "snn.slm_0.5", "snn.slm_1", "snn.slm_1.5")) {
  TARA_HEI[[cluster_col]][, 1] <- factor(
    TARA_HEI[[cluster_col]][, 1],
    levels = sort(as.numeric(levels(TARA_HEI[[cluster_col]][, 1])))
  )
}
for (cluster_col in c("snn.louvianmlr_0.5", "snn.louvianmlr_1", "snn.louvianmlr_1.5",
                      "snn.slm_0.5", "snn.slm_1", "snn.slm_1.5")) {
  EARTH[[cluster_col]][, 1] <- factor(
    EARTH[[cluster_col]][, 1],
    levels = sort(as.numeric(levels(EARTH[[cluster_col]][, 1])))
  )
}
####### Plot comparing cluster
# Load required library
library(patchwork)

# TARA ALL

# Azimuth plots (first row)
p1 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L1")

p2 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L2")

p3 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l3",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L3")

# SLM plots (second row)
p4 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.slm_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_0.5")

p5 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.slm_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1")

p6 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.slm_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1.5")

# LouvainMLR plots (third row)
p7 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_0.5")

p8 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

p9 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1.5")

# Combine all plots into a 3x3 grid
combined_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)

# Display the combined plot
print(combined_plot)


DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")
# Save the combined plot to a file
ggsave(
  filename = "TARA_ALL_WNN_Clusters.png",
  plot = combined_plot,
  width = 24,  # Adjust width as needed
  height = 17,  # Adjust height as needed
  dpi = 300
)


# TARA HEI

# Azimuth plots (first row)
p1 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L1")

p2 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L2")

p3 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l3",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L3")

# SLM plots (second row)
p4 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "snn.slm_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_0.5")

p5 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "snn.slm_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1")

p6 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "snn.slm_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1.5")

# LouvainMLR plots (third row)
p7 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_0.5")

p8 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

p9 <- DimPlot_scCustom(
  TARA_HEI,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1.5")

# Combine all plots into a 3x3 grid
combined_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)


# Save the combined plot to a file
ggsave(
  filename = "TARA_HEI_WNN_Clusters.png",
  plot = combined_plot,
  width = 24,  # Adjust width as needed
  height = 17,  # Adjust height as needed
  dpi = 300
)


# EARTH

# Azimuth plots (first row)
p1 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L1")

p2 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L2")

p3 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l3",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L3")

# SLM plots (second row)
p4 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.slm_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_0.5")

p5 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.slm_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1")

p6 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.slm_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1.5")

# LouvainMLR plots (third row)
p7 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_0.5")

p8 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

p9 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1.5")

# Combine all plots into a 3x3 grid
combined_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)

# Save the combined plot to a file
ggsave(
  filename = "EARTH_WNN_Clusters.png",
  plot = combined_plot,
  width = 24,  # Adjust width as needed
  height = 17,  # Adjust height as needed
  dpi = 300
)

#### Save to load path

save(TARA_ALL, file = paste0(load.path, "TARA_ALL_WNN.Rdata"))
save(TARA_HEI, file = paste0(load.path, "TARA_HEI_WNN.Rdata"))
save(EARTH, file = paste0(load.path, "EARTH_WNN.Rdata"))

