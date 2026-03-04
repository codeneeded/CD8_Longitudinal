############################################################
# CP003 longitudinal RNA: Read 10x H5 -> Seurat -> Metadata
############################################################

library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(qs2)
library(SingleCellExperiment)
library(scDblFinder)
library(SeuratExtend)
library(Azimuth)

# -----------------------------
# 0) Paths
# -----------------------------
base_dir <- "/home/akshay-iyer/CP003_multi_out/outs/per_sample_outs"

subfolders <- c(
  "CP003_2m_001A",
  "CP003_2m_001B",
  "CP003_101m_003",
  "CP003_101m_004"
)

# -----------------------------
# 1) Sample annotation map
# -----------------------------
sample_map <- tibble::tibble(
  Sample_Subfolder = subfolders,
  PID = "CP003",
  Months = as.integer(str_match(Sample_Subfolder, "CP003_(\\d+)m_")[,2]),
  OCM_Barcode = str_match(Sample_Subfolder, "CP003_\\d+m_(.+)$")[,2],
  CellType_Sort = dplyr::case_when(
    OCM_Barcode %in% c("001A","001B","003") ~ "CD8+",
    OCM_Barcode %in% c("004")              ~ "CD8-",
    TRUE                                   ~ NA_character_
  ),
  Timepoint_Label = paste0(Months, "m")
)

# -----------------------------
# 2) Read each H5 and make Seurat
# -----------------------------
objs <- list()

for (i in seq_len(nrow(sample_map))) {
  samp <- sample_map[i, ]
  
  h5_path <- file.path(base_dir, samp$Sample_Subfolder, "count", "sample_filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_path)) {
    stop("Missing H5: ", h5_path)
  }
  
  # Read 10x H5
  x <- Read10X_h5(h5_path)
  
  # Some H5s return a list (RNA/ADT/HTO). We only want RNA.
  if (is.list(x)) {
    if (!"Gene Expression" %in% names(x) && !"RNA" %in% names(x)) {
      stop("H5 returned a list but couldn't find 'Gene Expression' or 'RNA' in: ", h5_path,
           "\nNames found: ", paste(names(x), collapse = ", "))
    }
    if ("Gene Expression" %in% names(x)) {
      x_rna <- x[["Gene Expression"]]
    } else {
      x_rna <- x[["RNA"]]
    }
  } else {
    x_rna <- x
  }
  
  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = x_rna,
    project = "CP003",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
  
  # Add metadata
  obj$PID            <- samp$PID
  obj$Months         <- samp$Months
  obj$Timepoint      <- samp$Timepoint_Label
  obj$OCM_Barcode    <- samp$OCM_Barcode
  obj$CellType_Sort  <- samp$CellType_Sort
  obj$Sample_Subfolder <- samp$Sample_Subfolder
  
  # Make cell names unique after merge by prefixing with subfolder
  obj <- RenameCells(obj, add.cell.id = samp$Sample_Subfolder)
  
  objs[[samp$Sample_Subfolder]] <- obj
  
  message("Loaded: ", samp$Sample_Subfolder,
          " | cells=", ncol(obj),
          " | features=", nrow(obj))
}

# -----------------------------
# 3) Merge into one object
# -----------------------------
seu_CP003 <- Reduce(function(x, y) merge(x, y), objs)

############################################################
# CP003 longitudinal RNA QC
# - Uses your existing QC style, adapted to CP003 object
# - Saves outputs to: /home/akshay-iyer/Documents/CD8_Longitudinal/CP003/QC/Pre-QC + Post-QC
############################################################


DefaultAssay(seu_CP003) <- "RNA"
############################################################
# CP003 longitudinal RNA QC (split by Timepoint + CellType_Sort)
# Saves to: /home/akshay-iyer/Documents/CD8_Longitudinal/CP003/QC
############################################################


DefaultAssay(seu_CP003) <- "RNA"

# Convenience grouping column
seu_CP003$Sample <- seu_CP003$Sample_Subfolder

# -----------------------------
# 1) Output folders
# -----------------------------
out_root <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003"
qc_root  <- file.path(out_root, "QC")
pre_root <- file.path(qc_root, "Pre-QC")
post_root<- file.path(qc_root, "Post-QC")

dir.create(pre_root,  recursive = TRUE, showWarnings = FALSE)
dir.create(post_root, recursive = TRUE, showWarnings = FALSE)

save_png <- function(filename, plot_obj, width = 1800, height = 1200, res = 200) {
  png(filename = filename, width = width, height = height, res = res)
  print(plot_obj)
  dev.off()
}

# -----------------------------
# 2) QC metrics (adapted exactly from your template)
# -----------------------------
seu_CP003$log10GenesPerUMI <- log10(seu_CP003$nFeature_RNA) / log10(seu_CP003$nCount_RNA)

seu_CP003 <- PercentageFeatureSet(seu_CP003, pattern = "^MT-",      col.name = "percent_mito")
seu_CP003 <- PercentageFeatureSet(seu_CP003, pattern = "^RP[SL]",   col.name = "percent_ribo")
seu_CP003 <- PercentageFeatureSet(seu_CP003, pattern = "^HB[^(P)]", col.name = "percent_hb")
seu_CP003 <- PercentageFeatureSet(seu_CP003, pattern = "PECAM1|PF4",col.name = "percent_plat")

feats_qc <- c("nCount_RNA", "nFeature_RNA",
              "percent_mito", "percent_ribo", "percent_hb", "percent_plat",
              "log10GenesPerUMI")

# -----------------------------
# 3) PRE-QC plotting function (runs per grouping variable)
# -----------------------------
run_pre_qc_plots <- function(seu, group_col, out_dir) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- seu@meta.data
  
  # Cells per group
  p_cells <- meta %>%
    ggplot(aes(x = .data[[group_col]], fill = .data[[group_col]])) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(paste0("NCells by ", group_col))
  save_png(file.path(out_dir, "Cells_per_group.png"), p_cells)
  
  # Violin QC features grouped
  p_vln <- VlnPlot(seu, group.by = group_col, features = feats_qc, pt.size = 0.1, ncol = 2) + NoLegend()
  save_png(file.path(out_dir, "QC_features_grouped.png"), p_vln)
  
  # Density plots (your style)
  p_umi <- meta %>%
    ggplot(aes(color = .data[[group_col]], x = nCount_RNA, fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  save_png(file.path(out_dir, "UMI_Count.png"), p_umi)
  
  p_ngenes <- meta %>%
    ggplot(aes(color = .data[[group_col]], x = nFeature_RNA, fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 600)
  save_png(file.path(out_dir, "nGenes.png"), p_ngenes)
  
  p_complex <- meta %>%
    ggplot(aes(x = log10GenesPerUMI, color = .data[[group_col]], fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  save_png(file.path(out_dir, "Complexity_Score.png"), p_complex)
  
  p_mito <- meta %>%
    ggplot(aes(color = .data[[group_col]], x = percent_mito, fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 15)
  save_png(file.path(out_dir, "Mito_Ratio.png"), p_mito)
  
  p_ribo <- meta %>%
    ggplot(aes(color = .data[[group_col]], x = percent_ribo, fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 5)
  save_png(file.path(out_dir, "Ribo_Ratio.png"), p_ribo)
  
  p_hb <- meta %>%
    ggplot(aes(color = .data[[group_col]], x = percent_hb, fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 20)
  save_png(file.path(out_dir, "Heme_Ratio.png"), p_hb)
  
  p_plat <- meta %>%
    ggplot(aes(color = .data[[group_col]], x = percent_plat, fill = .data[[group_col]])) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 2)
  save_png(file.path(out_dir, "Platelet_Ratio.png"), p_plat)
  
  # scCustomize grouped cutoff panels
  p1 <- QC_Plots_Genes(seurat_object = seu, low_cutoff = 600, high_cutoff = 5500, group.by = group_col)
  p2 <- QC_Plots_UMIs(seurat_object = seu, low_cutoff = 1200, high_cutoff = 45000, group.by = group_col)
  p3 <- QC_Plots_Mito(seurat_object = seu, high_cutoff = 20, group.by = group_col)
  p4 <- QC_Plots_Complexity(seurat_object = seu, high_cutoff = 0.8, group.by = group_col)
  
  save_png(
    file.path(out_dir, "Grouped_Cutoff.png"),
    wrap_plots(p1, p2, p3, p4, ncol = 4),
    width = 3600, height = 1200, res = 200
  )
  
  # Scatter QC
  save_png(
    file.path(out_dir, "UMIvsGene.png"),
    QC_Plot_UMIvsGene(
      seurat_object = seu,
      low_cutoff_gene = 600, high_cutoff_gene = 5500,
      low_cutoff_UMI = 500, high_cutoff_UMI = 50000,
      group.by = group_col
    )
  )
  
  save_png(
    file.path(out_dir, "MitovsGene.png"),
    QC_Plot_GenevsFeature(
      seurat_object = seu,
      feature1 = "percent_mito",
      low_cutoff_gene = 600, high_cutoff_gene = 5500,
      high_cutoff_feature = 20,
      group.by = group_col
    )
  )
  
  save_png(
    file.path(out_dir, "MitovsGene_gradient.png"),
    QC_Plot_UMIvsGene(
      seurat_object = seu,
      meta_gradient_name = "percent_mito",
      low_cutoff_gene = 600, high_cutoff_gene = 5500,
      high_cutoff_UMI = 45000
    )
  )
}

# -----------------------------
# 4) Run PRE-QC for each grouping you care about
# -----------------------------
group_vars <- c("Sample_Subfolder", "Timepoint", "CellType_Sort")

for (g in group_vars) {
  run_pre_qc_plots(
    seu     = seu_CP003,
    group_col = g,
    out_dir = file.path(pre_root, g)
  )
}

# -----------------------------
# 5) FILTERING (your thresholds)
# -----------------------------
# JoinLayers if Seurat v5
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  seu_CP003 <- JoinLayers(seu_CP003)
}
DefaultAssay(seu_CP003) <- "RNA"

filtered_CP003 <- subset(
  x = seu_CP003,
  subset =
    (nCount_RNA >= 500) &
    (nFeature_RNA >= 600) &
    (log10GenesPerUMI > 0.80) &
    (percent_mito < 15) &
    (percent_ribo > 5) &
    (percent_hb < 20) &
    (percent_plat < 2)
)

# Remove genes expressed in >=10 cells
counts_mat <- NULL
if (!is.null(filtered_CP003@assays$RNA@counts) && nrow(filtered_CP003@assays$RNA@counts) > 0) {
  counts_mat <- filtered_CP003@assays$RNA@counts
} else if (!is.null(filtered_CP003@assays$RNA@layers$counts)) {
  counts_mat <- filtered_CP003@assays$RNA@layers$counts
} else {
  stop("Could not find RNA counts matrix (counts slot or layers$counts).")
}

genes_in_10_cells <- Matrix::rowSums(counts_mat > 0) >= 10
filtered_CP003 <- subset(filtered_CP003, features = names(genes_in_10_cells[genes_in_10_cells]))

# -----------------------------
# 6) POST-QC plots (same loop)
# -----------------------------
for (g in group_vars) {
  run_pre_qc_plots(
    seu     = filtered_CP003,
    group_col = g,
    out_dir = file.path(post_root, g)
  )
}

############################################################
# CP003: Cell cycle scoring + scDblFinder doublet removal
# - Input: filtered_CP003  (post-QC singlets not yet removed)
# - Output: CP003_postQC_cellcycle_scDblFinder_singlets.qs2 (saved_R_data only)
# - Also writes QC plots under: /home/akshay-iyer/Documents/CD8_Longitudinal/CP003/QC
############################################################




# -----------------------------
# 0) Paths
# -----------------------------
qc_root   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/QC"
cc_dir    <- file.path(qc_root, "Cell_Cycle")
dbl_dir   <- file.path(qc_root, "Doublets")
dir.create(cc_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(dbl_dir, recursive = TRUE, showWarnings = FALSE)

save_dir  <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

out_qs2   <- file.path(save_dir, "CP003_postQC_cellcycle_scDblFinder_singlets.qs2")
filtered_CP003 <- JoinLayers(filtered_CP003)
# -----------------------------
# 1) Input object (post-QC)
# -----------------------------
# Expecting you already created this in the QC script:
#   filtered_CP003
stopifnot(exists("filtered_CP003"))

seu <- filtered_CP003
DefaultAssay(seu) <- "RNA"

# Use your library/sample identifier for splitting (instead of orig.ident)
seu$Sample <- seu$Sample_Subfolder

# -----------------------------
# 2) Normalize + HVGs + PCA/UMAP (lightweight; enough to visualize cell cycle & doublets)
# -----------------------------
# If already normalized in your pipeline, this will just overwrite standard steps cleanly.
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap.pre", verbose = FALSE)

# -----------------------------
# 3) Cell cycle scoring (built-in Seurat gene sets)
# -----------------------------
# Seurat provides cc.genes.updated.2019
s.genes   <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# Make sure genes exist in your object (important for custom references / gene naming)
s.genes.use   <- intersect(s.genes, rownames(seu))
g2m.genes.use <- intersect(g2m.genes, rownames(seu))


seu <- CellCycleScoring(
  seu,
  s.features   = s.genes.use,
  g2m.features = g2m.genes.use,
  set.ident = T
)

############################################################
# Cell-cycle visualizations (DimPlot2 / VlnPlot2 + ggsave)
############################################################

# UMAP by Cell Cycle Phase
p_cc_umap <- DimPlot2(
  seu,
  reduction = "umap.pre",
  group.by = "Phase",
  pt.size = 0.4
)

ggsave(
  filename = file.path(cc_dir, "UMAP_by_CellCyclePhase.png"),
  plot = p_cc_umap,
  width = 10,
  height = 6,
  dpi = 600
)

# ----------------------------------------------------------

# Violin: Cell cycle scores by Sample
p_cc_sample <- VlnPlot2(
  seu,
  features = c("S.Score", "G2M.Score"),
  group.by = "Sample",
  pt.size = 0.1
)

ggsave(
  filename = file.path(cc_dir, "CellCycle_Score_Vln_by_Sample.png"),
  plot = p_cc_sample,
  width = 10,
  height = 6,
  dpi = 600
)

# ----------------------------------------------------------

# Violin: Cell cycle scores by Timepoint
p_cc_timepoint <- VlnPlot2(
  seu,
  features = c("S.Score", "G2M.Score"),
  group.by = "Timepoint",
  pt.size = 0.1
)

ggsave(
  filename = file.path(cc_dir, "CellCycle_Score_Vln_by_Timepoint.png"),
  plot = p_cc_timepoint,
  width = 10,
  height = 6,
  dpi = 600
)

# ----------------------------------------------------------

# Violin: Cell cycle scores by CellType_Sort
p_cc_sort <- VlnPlot2(
  seu,
  features = c("S.Score", "G2M.Score"),
  group.by = "CellType_Sort",
  pt.size = 0.1
)

ggsave(
  filename = file.path(cc_dir, "CellCycle_Score_Vln_by_CellType_Sort.png"),
  plot = p_cc_sort,
  width = 10,
  height = 6,
  dpi = 600
)

# Optional ridge markers (only if present)
ridge_feats <- intersect(c("PCNA","TOP2A","MCM6","MKI67"), rownames(seu))
if (length(ridge_feats) > 0) {
  png(file.path(cc_dir, "Ridge_CellCycle_Markers.png"), width = 2100, height = 1200, res = 200)
  print(RidgePlot(seu, features = ridge_feats, ncol = 2))
  dev.off()
}

# Azimuth predicted cell types
seu <- RunAzimuth(seu, reference = "pbmcref")

p_azimuth <- DimPlot2(
  seu,
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  box = TRUE
)

ggsave(
  filename = file.path(cc_dir, "UMAP_Azimuth_predicted_celltype_l2.png"),
  plot = p_azimuth,
  width = 10,
  height = 8,
  dpi = 600,
  bg='white'
)

# Split by CD8 sort
p_azimuth_sort <- DimPlot2(
  seu,
  group.by = "predicted.celltype.l2",
  split.by = "CellType_Sort",
  label = TRUE,
  repel = TRUE,
  box = TRUE
)

ggsave(
  filename = file.path(cc_dir, "UMAP_Azimuth_predicted_celltype_l2_by_Sort.png"),
  plot = p_azimuth_sort,
  width = 14,
  height = 8,
  dpi = 600,
  bg='white'
)

# Split by Group sort
p_azimuth_sample <- DimPlot2(
  seu,
  group.by = "predicted.celltype.l2",
  split.by = "Sample",
  label = TRUE,
  repel = TRUE,
  box = TRUE
)

ggsave(
  filename = file.path(cc_dir, "UMAP_Azimuth_predicted_celltype_l2_by_Sample.png"),
  plot = p_azimuth_sample,
  width = 14,
  height = 8,
  dpi = 600,
  bg='white'
)
# -----------------------------
# 4) Doublet detection with scDblFinder (split per library)
# -----------------------------
# Split by Sample_Subfolder (your true “library” unit)
split_seu <- SplitObject(seu, split.by = "Sample")

samples <- names(split_seu)

# Doublet calling + plots
for (i in samples) {
  
  sce <- scDblFinder(GetAssayData(split_seu[[i]], assay = "RNA", layer = "counts"))
  
  split_seu[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seu[[i]]$scDblFinder.class <- sce$scDblFinder.class
  
  p <- DimPlot2(
    split_seu[[i]],
    reduction = "umap.pre",
    group.by = "scDblFinder.class"
  ) + ggtitle(paste0(i, " Doublets"))
  
  ggsave(
    filename = file.path(dbl_dir, paste0(i, "_doublets.png")),
    plot = p,
    device = "png",
    width = 6,
    height = 6,
    dpi = 600
  )
}

# Remove doublets (keep singlets)
for (i in samples) {
  split_seu[[i]] <- subset(split_seu[[i]], subset = scDblFinder.class == "singlet")
}
# Recombine (merge) back into one Seurat object
seu_singlets <- Merge_Seurat_List(split_seu,merge.data = TRUE)
seu_singlets <- JoinLayers(seu_singlets)

# (optional but nice) reset ident
Idents(seu_singlets) <- "Sample"   # or "Timepoint" / "CellType_Sort"

# Save (qs2 only, in saved_R_data only)
qs_save(seu_singlets, file = out_qs2)
