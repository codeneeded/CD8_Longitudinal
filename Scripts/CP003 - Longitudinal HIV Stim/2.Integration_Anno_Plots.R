############################################################
# CP003 RNA Integration (Seurat v5) from qs2 singlets object
# - Loads: CP003_postQC_cellcycle_scDblFinder_singlets.qs2
# - Unintegrated: Normalize/Variable/Scale/PCA/Neighbors/UMAP
# - Integrations: CCA + FastMNN via IntegrateLayers
# - Clustering: clustree across resolutions for CCA + MNN
# - Plotting: DimPlot2 everywhere + ggsave
# - Saves integrated object back to saved_R_data (qs2 only)
############################################################


library(qs2)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)
library(SeuratExtend)   # for DimPlot2 / VlnPlot2


# -----------------------------
# 0) Paths
# -----------------------------
save_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data"
in_qs2   <- file.path(save_dir, "CP003_postQC_cellcycle_scDblFinder_singlets.qs2")
out_qs2  <- file.path(save_dir, "CP003_RNA_integrated_CCA_MNN.qs2")

plot_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/Integration"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Load object
# -----------------------------
seu <- qs_read(in_qs2)
DefaultAssay(seu) <- "RNA"


# Your "batch" for integration (4 libraries)
seu$batch <- seu$Sample_Subfolder

# -----------------------------
# 2) Unintegrated processing (RNA)
# IMPORTANT: split RNA layers by batch like your example
# -----------------------------
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$batch)

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca", verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, cluster.name = "unintegrated_clusters", verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated", verbose = FALSE)

# Unintegrated quick plots
p_unint_1 <- DimPlot2(seu, reduction = "umap.unintegrated", group.by = "batch")
p_unint_2 <- DimPlot2(seu, reduction = "umap.unintegrated", group.by = "Timepoint")
p_unint_3 <- DimPlot2(seu, reduction = "umap.unintegrated", group.by = "CellType_Sort")
p_unint_4 <- if ("predicted.celltype.l2" %in% colnames(seu@meta.data)) {
  DimPlot2(seu, reduction = "umap.unintegrated", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE, box = TRUE)
} else {
  DimPlot2(seu, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", label = TRUE, repel = TRUE, box = TRUE)
}

ggsave(file.path(plot_dir, "CP003_RNA_Unintegrated_UMAPs.png"),
       plot = wrap_plots(p_unint_1, p_unint_2, p_unint_3, p_unint_4, ncol = 2),
       width = 14, height = 10, dpi = 300)

# -----------------------------
# 3) Integration: CCA
# -----------------------------
seu <- IntegrateLayers(
  object = seu,
  method = CCAIntegration,
  orig.reduction = "pca",
  assay = "RNA",
  new.reduction = "integrated.cca.rna",
  verbose = FALSE
)

# -----------------------------
# 4) Integration: FastMNN
# -----------------------------
seu <- IntegrateLayers(
  object = seu,
  method = FastMNNIntegration,
  assay = "RNA",
  new.reduction = "integrated.mnn.rna",
  verbose = FALSE
)

# -----------------------------
# 5) Clustering + clustree (CCA and MNN)
#   - run a resolution grid and save clustree plots
# -----------------------------
res_grid <- seq(0.2, 3.0, by = 0.2)

# --- CCA graph + clusters across resolutions
seu <- FindNeighbors(seu, reduction = "integrated.cca.rna", dims = 1:30, graph.name = "cca_snn", verbose = FALSE)

for (r in res_grid) {
  seu <- FindClusters(
    seu,
    graph.name = "cca_snn",
    resolution = r,
    algorithm = 1,
    verbose = FALSE
  )
  # This creates meta cols like: cca_snn_res.0.2, cca_snn_res.0.4, ...
}

p_tree_cca <- clustree(seu, prefix = "cca_snn_res.")
ggsave(file.path(plot_dir, "CP003_clustree_CCA.png"),
       plot = p_tree_cca, width = 14, height = 10, dpi = 300)

# pick one resolution to finalize (you will decide after viewing clustree)
# default pick:
cca_res_final <- 0.4
cca_col <- paste0("cca_snn_res.", cca_res_final)
seu$cca_clusters_rna <- seu[[cca_col]][,1]

seu <- RunUMAP(seu, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna", verbose = FALSE)

# --- MNN graph + clusters across resolutions
seu <- FindNeighbors(seu, reduction = "integrated.mnn.rna", dims = 1:30, graph.name = "mnn_snn", verbose = FALSE)

for (r in res_grid) {
  seu <- FindClusters(
    seu,
    graph.name = "mnn_snn",
    resolution = r,
    algorithm = 1,
    verbose = FALSE
  )
}

p_tree_mnn <- clustree(seu, prefix = "mnn_snn_res.")
ggsave(file.path(plot_dir, "CP003_clustree_MNN.png"),
       plot = p_tree_mnn, width = 14, height = 10, dpi = 300)

# default pick:
mnn_res_final <- 0.6
mnn_col <- paste0("mnn_snn_res.", mnn_res_final)
if (!mnn_col %in% colnames(seu@meta.data)) stop("Chosen MNN resolution column not found: ", mnn_col)
seu$mnn_clusters_rna <- seu[[mnn_col]][,1]

seu <- RunUMAP(seu, reduction = "integrated.mnn.rna", dims = 1:30, reduction.name = "umap.mnn.rna", verbose = FALSE)

# -----------------------------
# 6) Plot integrated embeddings (DimPlot2 only)
# -----------------------------
# CCA panels
p1 <- if ("predicted.celltype.l2" %in% colnames(seu@meta.data)) {
  DimPlot2(seu, reduction = "umap.cca.rna", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE, box = TRUE)
} else {
  DimPlot2(seu, reduction = "umap.cca.rna", group.by = "cca_clusters_rna", label = TRUE, repel = TRUE, box = TRUE)
}
p2 <- DimPlot2(seu, reduction = "umap.cca.rna", group.by = "cca_clusters_rna", label = TRUE, repel = TRUE, box = TRUE)

# MNN panels
p3 <- if ("predicted.celltype.l2" %in% colnames(seu@meta.data)) {
  DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE, box = TRUE)
} else {
  DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "mnn_clusters_rna", label = TRUE, repel = TRUE, box = TRUE)
}
p4 <- DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "mnn_clusters_rna", label = TRUE, repel = TRUE, box = TRUE)

combined_plot <- wrap_plots(p1, p2, p3, p4, ncol = 2, byrow = TRUE)

ggsave(
  filename = file.path(plot_dir, "CP003_RNA_Integration_CCA_MNN.png"),
  plot = combined_plot,
  width = 18,
  height = 12,
  dpi = 300
)

# Also save batch/timepoint checks (post integration)
p_chk_cca <- wrap_plots(
  DimPlot2(seu, reduction = "umap.cca.rna", group.by = "batch"),
  DimPlot2(seu, reduction = "umap.cca.rna", group.by = "Timepoint"),
  DimPlot2(seu, reduction = "umap.cca.rna", group.by = "CellType_Sort"),
  ncol = 3
)
ggsave(file.path(plot_dir, "CP003_CCA_UMAP_checks.png"), p_chk_cca, width = 18, height = 6, dpi = 300)

p_chk_mnn <- wrap_plots(
  DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "batch"),
  DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "Timepoint"),
  DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "CellType_Sort"),
  ncol = 3
)
ggsave(file.path(plot_dir, "CP003_MNN_UMAP_checks.png"), p_chk_mnn, width = 18, height = 6, dpi = 300)

p <- DimPlot2(seu, reduction = "umap.mnn.rna", split.by = "CellType_Sort")
ggsave(file.path(plot_dir, "CP003_MNN_UMAP_sample.png"), p, width = 8, height = 6, dpi = 300)

# -----------------------------
# 7) Save integrated object (qs2 only)
# -----------------------------
qs_save(seu, file = out_qs2)

################
############################################################
# CP003 annotation helper plots (MNN integrated; mnn_res_final = 0.6)
# - Loads integrated qs2 object
# - Makes VlnPlot2 + FeaturePlot2 for RNA features
# - Makes ClusterDistrBar for cluster composition
# - Saves to /home/akshay-iyer/Documents/CD8_Longitudinal/CP003/Annotation_Plots/...
############################################################


# -----------------------------

# Use MNN clusters @ res 0.6 (as you decided)
# If you already created seu$mnn_clusters_rna earlier, keep it.
# Otherwise, pull from the stored clustering column.
mnn_res_final <- 0.6
mnn_col <- paste0("mnn_snn_res.", mnn_res_final)

if (!"mnn_clusters_rna" %in% colnames(seu@meta.data)) {
  if (!mnn_col %in% colnames(seu@meta.data)) stop("Cannot find MNN clustering column: ", mnn_col)
  seu$mnn_clusters_rna <- seu[[mnn_col]][,1]
}
Idents(seu) <- "mnn_clusters_rna"

# -----------------------------
# 1) Output folders
# -----------------------------
out_root <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/Annotation_Plots"
dir_vln  <- file.path(out_root, "VlnPlot2_RNA")
dir_ftr  <- file.path(out_root, "FeaturePlot2_RNA")
dir_dist <- file.path(out_root, "Cluster_Distribution")

dir.create(dir_vln,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_ftr,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_dist, recursive = TRUE, showWarnings = FALSE)
seu$mnn_clusters_rna <- factor(
  seu$mnn_snn_res.0.6,
  levels = as.character(sort(as.numeric(unique(as.character(seu$mnn_snn_res.0.6)))))
)
Idents(seu) <- "mnn_clusters_rna"
levels(Idents(seu))
# -----------------------------
# 2) Features
# -----------------------------
rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B',
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FAS','FASLG','FABP5','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18',
  'TNFRSF4','TNFRSF9','TOX','TBX21','TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1',
  'XCL1','XCL2','ZBTB16','ZEB2'
)

# Keep only genes that exist (avoid errors)
rna.features.use <- intersect(rna.features, rownames(seu))
missing_genes <- setdiff(rna.features, rna.features.use)
if (length(missing_genes) > 0) message("Missing (not in object): ", paste(missing_genes, collapse = ", "))

# -----------------------------
# 3) VlnPlot2 (per gene) across MNN clusters
# -----------------------------
for (g in rna.features.use) {
  
  vln.pl <- VlnPlot2(
    seu,
    features = g,
    cols = "default",
    show.mean = TRUE,
    mean_colors = c("red", "blue")
  ) + ggtitle(paste("RNA |", g))
  
  ggsave(
    filename = file.path(dir_vln, paste0(g, "_VLNplot.png")),
    plot = vln.pl,
    dpi = 500,
    width = 14,
    height = 8,
    bg = "white"
  )
}

# -----------------------------
# 4) FeaturePlot2 (per gene) on MNN UMAP
# -----------------------------

for (g in rna.features.use) {
  
  fp <- FeaturePlot_scCustom(
    seu,
    features = g,
    reduction = "umap.mnn.rna"
  ) + ggtitle(paste("RNA |", g))
  
  ggsave(
    filename = file.path(dir_ftr, paste0(g, "_FeaturePlot.png")),
    plot = fp,
    dpi = 500,
    width = 10,
    height = 8,
    bg = "white"
  )
}

# -----------------------------
# 5) Cluster distribution barplot
# -----------------------------
# Cluster distribution (counts)
p_dist_counts <- ClusterDistrBar(
  origin  = seu$Sample_Subfolder,
  cluster = seu$mnn_clusters_rna,
  cols    = "default",
  flip    = FALSE,
  border  = "black",
  percent = FALSE
) +
  theme(axis.title.x = element_blank()) +
  ggtitle("Cluster distribution | mnn_clusters_rna (counts)")

ggsave(
  filename = file.path(dir_dist, "CP003_ClusterDistrBar_counts.png"),
  plot = p_dist_counts,
  dpi = 400,
  width = 14,
  height = 7,
  bg = "white"
)

# Cluster distribution (percent)
p_dist_percent <- ClusterDistrBar(
  origin  = seu$Sample_Subfolder,
  cluster = seu$mnn_clusters_rna,
  cols    = "default",
  flip    = FALSE,
  border  = "black",
  percent = TRUE
) +
  theme(axis.title.x = element_blank()) +
  ggtitle("Cluster distribution | mnn_clusters_rna (percent)")

ggsave(
  filename = file.path(dir_dist, "CP003_ClusterDistrBar_percent.png"),
  plot = p_dist_percent,
  dpi = 400,
  width = 14,
  height = 7,
  bg = "white"
)
