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

for (col in grep("^wsnn_res\\.", colnames(tara_cdnk@meta.data), value = TRUE)) {
  lvls <- as.character(sort(as.numeric(levels(tara_cdnk[[col]][, 1]))))
  tara_cdnk[[col]] <- factor(tara_cdnk[[col]][, 1], levels = lvls)
}

setwd('/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets')



clustree(tara_cdnk, prefix = "wsnn_res.")


ggsave('clustree.png',  width = 15,  # Adjust width as needed
       height = 9,  # Adjust height as needed
       dpi = 300,
       bg='white')


DimPlot2(
  tara_cdnk,
  reduction = "wnn.umap",
  group.by = "wsnn_res.0.4",
  cols = 'light',
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
)

ggsave('cd8nk_subcluster_res0.4.png',
       dpi = 300,
       bg='white')


saveRDS(
  tara_cdnk,
  file = "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data/tara_cdnk.rds"
)

############################### ANNOTATION PLOTTING ##############################################

# Base output directory
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets"

# Create subfolders
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "RNA_VLN"), showWarnings = FALSE)
dir.create(file.path(out_dir, "ADT_VLN"), showWarnings = FALSE)
dir.create(file.path(out_dir, "RNA_FeaturePlot"), showWarnings = FALSE)
dir.create(file.path(out_dir, "ADT_FeaturePlot"), showWarnings = FALSE)

rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B',
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FASLG','FAS','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18',
  'TNFRSF4','TNFRSF9','TOX','TBX21','TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1',
  'XCL1','XCL2','ZBTB16','ZEB2'
)

prot.features <- rownames(tara_cdnk@assays$ADT)


### RNA VLN PLOT
DefaultAssay(tara_cdnk) <- 'RNA'
for (gene in rna.features) {
  if (gene %in% rownames(tara_cdnk[["RNA"]])) {
    vln.pl <- VlnPlot2(
      tara_cdnk,
      features = gene,
      group.by = 'wsnn_res.0.4',
      cols = 'light',
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("RNA |", gene))
    
    ggsave(
      filename = file.path(out_dir, "RNA_VLN", paste0(gene, "_VLN.png")),
      plot = vln.pl,
      dpi = 500, width = 14, height = 8, bg = "white"
    )
  }
}

### RNA Feature Plot
for (gene in rna.features) {
  if (gene %in% rownames(tara_cdnk[["RNA"]])) {
    fp <- DimPlot2(
      tara_cdnk,
      features = gene,
      reduction = "wnn.umap"
    ) + ggtitle(paste("RNA |", gene))
    
    ggsave(
      filename = file.path(out_dir, "RNA_FeaturePlot", paste0(gene, "_FeaturePlot.png")),
      plot = fp,
      dpi = 500, width = 8, height = 6, bg = "white"
    )
  }
}

### ADT VLN Plot
DefaultAssay(tara_cdnk) <- 'ADT'

for (prot in prot.features) {
  if (prot %in% rownames(tara_cdnk[["ADT"]])) {
    vln.pl <- VlnPlot2(
      tara_cdnk,
      features = prot,
      group.by = 'wsnn_res.0.4',
      assay = "ADT",
      cols = 'light',
      show.mean = TRUE,
      mean_colors = c("red", "blue")
    ) + ggtitle(paste("ADT |", prot))
    
    ggsave(
      filename = file.path(out_dir, "ADT_VLN", paste0(prot, "_VLN.png")),
      plot = vln.pl,
      dpi = 500, width = 14, height = 8, bg = "white"
    )
  }
}

for (prot in prot.features) {
  if (prot %in% rownames(tara_cdnk[["ADT"]])) {
    fp <- DimPlot2(
      tara_cdnk,
      features = prot,
      assay = "ADT",
      reduction = "wnn.umap"
    ) + ggtitle(paste("ADT |", prot))
    
    ggsave(
      filename = file.path(out_dir, "ADT_FeaturePlot", paste0(prot, "_FeaturePlot.png")),
      plot = fp,
      dpi = 500, width = 8, height = 6, bg = "white"
    )
  }
}


?DimPlot2






##########3
TARA_ALL$Viral_Load
TARA_entry <- subset(TARA_ALL, subset = Age %in% c(1, 2) & Condition == "HEI")  # VL is defined only for HEI
TARA_entry$Viral_Load_Category <- ifelse(TARA_entry$Viral_Load >= 100000, "High", "Low")

VlnPlot2(
  TARA_entry,
  features = 'HLA-C',
  group.by = 'Manual_Annotation',
  split.by = 'Viral_Load_Category',
  assay = "RNA",
  cols = 'light',
  stat.method = 'wilcox.test',
  show.mean = TRUE,
  mean_colors = c("red", "blue")
) 
