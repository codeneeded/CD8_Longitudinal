library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(SeuratExtend)
library(scCustomize)
library(tidyr)
library(clustree)

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "TARA_ALL_post_annotation.rds")

TARA_ALL <- readRDS(rds_in)
DefaultAssay(TARA_ALL) <- 'RNA'
assay_use   <- DefaultAssay(TARA_ALL)

# Subset by Manual_Annotation
tara_cdnk <- subset(
  TARA_ALL,
  subset = Manual_Annotation %in% c(
    "1: Memory CD8 T cell",
    "16: Gamma Delta 1 T cells",
    "8: CTL-like",
    "9: TRDV1+ CTL-like",
    "27: GZMK+ CD8 T cell",
    "6: Naïve CD8 T cell",
    "21: Gamma Delta 2 T cells",
    "3: IL2RB+ NK cell",
    "17: NK cell_2",
    "15: NK cell_1",
    "25: CD56bright NK"
  )
)


# Recalculate modality reductions
tara_cdnk <- tara_cdnk %>%
  NormalizeData(assay = "RNA") %>%
  FindVariableFeatures(assay = "RNA", selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(assay = "RNA") %>%
  RunPCA(assay = "RNA", reduction.name = "pca")


# --- ADT (DSB-normalized): DO NOT NormalizeData; just scale + PCA ---
DefaultAssay(tara_cdnk) <- "ADT"

# Use all proteins (or provide a curated panel vector instead)
VariableFeatures(tara_cdnk) <- rownames(tara_cdnk[["ADT"]])

# Center/scale works fine on DSB values (they can be < 0)
tara_cdnk <- ScaleData(tara_cdnk, assay = "ADT", features = VariableFeatures(tara_cdnk), verbose = FALSE)

# PCA on ADT (DSB)
tara_cdnk <- RunPCA(tara_cdnk, assay = "ADT", features = VariableFeatures(tara_cdnk),
                    reduction.name = "apca")

# --- Rebuild WNN (always redo after subsetting multimodal data) ---
tara_cdnk <- FindMultiModalNeighbors(
  tara_cdnk,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30, 1:20)  # tweak as needed
)

# --- Cluster + UMAP on the new WNN graph ---
tara_cdnk <- RunUMAP(tara_cdnk, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


# Define the range of resolutions you want to test
resolutions <- seq(0.2, 2.0, by = 0.2)

for (res in resolutions) {
  tara_cdnk <- FindClusters(
    tara_cdnk,
    graph.name = "wsnn",
    resolution = res,
    algorithm = 3,  # Leiden (recommended)
    verbose = FALSE
  )
}

setwd('/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets')
clustree(tara_cdnk, prefix = "wsnn_res.")
ggsave('clustree.png',  width = 15,  # Adjust width as needed
       height = 9,  # Adjust height as needed
       dpi = 300)

DimPlot2(
  tara_cdnk,
  reduction = "wnn.umap",
  group.by = "wsnn_res.0.4",
  cols = 'light',
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
)

ggsave('cd8nk_subcluster_res0.4.png',
       dpi = 300)
