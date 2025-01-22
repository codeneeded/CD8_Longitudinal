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
  new.reduction = "integrated.cca"
)

TARA_ALL <- IntegrateLayers(
  TARA_ALL,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn"
)

### TARA HEI
TARA_HEI <- IntegrateLayers(
  TARA_HEI,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn"
)

TARA_HEI <- IntegrateLayers(
  TARA_HEI,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca"
)

### EARTH
EARTH <- IntegrateLayers(
  EARTH,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn"
)


EARTH <- IntegrateLayers(
  EARTH,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca"
)



#### Save to load path

saveRDS(TARA_ALL, file = paste0(load.path, "TARA_ALL.rds"))
saveRDS(TARA_HEI, file = paste0(load.path, "TARA_HEI.rds"))
saveRDS(EARTH, file = paste0(load.path, "EARTH.rds"))

