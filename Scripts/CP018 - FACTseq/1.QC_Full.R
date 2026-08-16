############################################################
# CP018 FACT-Seq — QC, cell cycle, doublet removal
#
# Adapted from "CP003 - Longitudinal HIV Stim/1. QC_Full.R".
# Differences from CP003:
#   1. Subfolders are CP018_2m_A / _B, CP018_90m_A / _B (not CP003_2m_001A etc.)
#   2. No CD8+/CD8- split — all four libraries are the same sort, so
#      CellType_Sort is constant and is NOT used as a grouping variable.
#   3. Replicate is a real grouping variable here (A/B at BOTH timepoints),
#      which CP003 could not support.
#
# Input : CR_OUT/<subfolder>/count/sample_filtered_feature_bc_matrix.h5
# Output: QS_SINGLETS + QC plots under QC_ROOT
############################################################

# Works whether run interactively (RStudio, from the .Rproj) or via Rscript.
source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(scCustomize)
library(qs2)
library(SingleCellExperiment)
library(scDblFinder)
library(SeuratExtend)
library(Azimuth)

pre_root  <- file.path(QC_ROOT, "Pre-QC")
post_root <- file.path(QC_ROOT, "Post-QC")
cc_dir    <- file.path(QC_ROOT, "Cell_Cycle")
dbl_dir   <- file.path(QC_ROOT, "Doublets")
mkdirs(pre_root, post_root, cc_dir, dbl_dir, SAVED_DIR)

# -----------------------------
# 1) Read each H5 -> Seurat
# -----------------------------
objs <- list()

for (i in seq_len(nrow(CP018_SAMPLES))) {
  samp <- CP018_SAMPLES[i, ]

  # Cell Ranger 10.x puts the H5 directly in the sample dir; 9.x used count/.
  h5_candidates <- c(
    file.path(CR_OUT, samp$Sample_Subfolder, "sample_filtered_feature_bc_matrix.h5"),
    file.path(CR_OUT, samp$Sample_Subfolder, "count", "sample_filtered_feature_bc_matrix.h5")
  )
  h5_path <- h5_candidates[file.exists(h5_candidates)][1]
  if (is.na(h5_path)) {
    stop("No filtered matrix H5 for ", samp$Sample_Subfolder, ". Looked in:\n  ",
         paste(h5_candidates, collapse = "\n  "), call. = FALSE)
  }

  x <- Read10X_h5(h5_path)

  # OCM GEX H5 may return a list; take Gene Expression
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      x_rna <- x[["Gene Expression"]]
    } else if ("RNA" %in% names(x)) {
      x_rna <- x[["RNA"]]
    } else {
      stop("H5 list has no 'Gene Expression'/'RNA': ", h5_path,
           "\nNames: ", paste(names(x), collapse = ", "))
    }
  } else {
    x_rna <- x
  }

  obj <- CreateSeuratObject(
    counts = x_rna, project = "CP018", assay = "RNA",
    min.cells = 3, min.features = 200
  )

  obj$PID              <- samp$PID
  obj$Months           <- samp$Months
  obj$Timepoint        <- samp$Timepoint
  obj$Replicate        <- samp$Replicate
  obj$OCM_Barcode      <- samp$OCM_Barcode
  obj$Tube_Label       <- samp$Tube_Label
  obj$CellType_Sort    <- samp$CellType_Sort
  obj$Sample_Subfolder <- samp$Sample_Subfolder
  # Timepoint_Rep is the useful CP018-specific grouping (2m_A, 2m_B, ...)
  obj$Timepoint_Rep    <- paste0(samp$Timepoint, "_", samp$Replicate)

  obj <- RenameCells(obj, add.cell.id = samp$Sample_Subfolder)
  objs[[samp$Sample_Subfolder]] <- obj

  message("Loaded: ", samp$Sample_Subfolder,
          " | cells=", ncol(obj), " | features=", nrow(obj))
}

# -----------------------------
# 2) Merge
# -----------------------------
seu_CP018 <- Reduce(function(x, y) merge(x, y), objs)
DefaultAssay(seu_CP018) <- "RNA"
seu_CP018$Sample <- seu_CP018$Sample_Subfolder

# Record raw per-library cell counts (needed for concordance in script 2)
raw_counts <- as.data.frame(table(seu_CP018$Sample_Subfolder))
names(raw_counts) <- c("Sample_Subfolder", "n_cells_raw")
write.csv(raw_counts, file.path(QC_ROOT, "CP018_cells_per_library_raw.csv"),
          row.names = FALSE)
print(raw_counts)

# -----------------------------
# 3) QC metrics
# -----------------------------
seu_CP018$log10GenesPerUMI <- log10(seu_CP018$nFeature_RNA) / log10(seu_CP018$nCount_RNA)

seu_CP018 <- PercentageFeatureSet(seu_CP018, pattern = "^MT-",       col.name = "percent_mito")
seu_CP018 <- PercentageFeatureSet(seu_CP018, pattern = "^RP[SL]",    col.name = "percent_ribo")
seu_CP018 <- PercentageFeatureSet(seu_CP018, pattern = "^HB[^(P)]",  col.name = "percent_hb")
seu_CP018 <- PercentageFeatureSet(seu_CP018, pattern = "PECAM1|PF4", col.name = "percent_plat")

feats_qc <- c("nCount_RNA", "nFeature_RNA",
              "percent_mito", "percent_ribo", "percent_hb", "percent_plat",
              "log10GenesPerUMI")

# -----------------------------
# 4) QC plotting function (unchanged in style from CP003)
# -----------------------------
run_qc_plots <- function(seu, group_col, out_dir) {

  mkdirs(out_dir)
  meta <- seu@meta.data

  p_cells <- meta %>%
    ggplot(aes(x = .data[[group_col]], fill = .data[[group_col]])) +
    geom_bar(colour = "black", linewidth = 0.3, width = 0.7) +
    geom_text(stat = "count", aes(label = after_stat(count)),
              vjust = -0.4, size = 5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                       labels = scales::label_number(big.mark = ",")) +
    labs(x = NULL, y = "Cells", title = paste0("Cells by ", group_col)) +
    theme_cp018() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
          legend.position = "none")
  save_png(file.path(out_dir, "Cells_per_group.png"), p_cells,
           width = 8, height = 5.5)

  # pt.size = 0: with thousands of cells the jittered points form a solid mass
  # that hides the violin. Distribution shape is the point of this plot.
  p_vln <- VlnPlot(seu, group.by = group_col, features = feats_qc,
                   pt.size = 0, ncol = 2, log = FALSE) &
    theme_cp018(base_size = 15) &
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                                labels = scales::label_number(big.mark = ",",
                                                              drop0trailing = TRUE)) &
    ggplot2::labs(x = NULL) &
    ggplot2::theme(legend.position = "none")
  ps <- panel_size(length(feats_qc), ncol = 2, panel_w = 5.2, panel_h = 3.4)
  save_png(file.path(out_dir, "QC_features_grouped.png"), p_vln,
           width = ps[["width"]], height = ps[["height"]])

  dens <- function(var, vline, logx = TRUE) {
    p <- meta %>%
      ggplot(aes(colour = .data[[group_col]], x = .data[[var]],
                 fill = .data[[group_col]])) +
      geom_density(alpha = 0.25, linewidth = 0.7) +
      geom_vline(xintercept = vline, linetype = "dashed",
                 colour = "grey25", linewidth = 0.7) +
      labs(x = var, y = "Cell density", colour = NULL, fill = NULL,
           title = var, subtitle = paste0("dashed line = QC threshold (", vline, ")")) +
      theme_cp018()
    if (logx) {
      p <- p + scale_x_log10(labels = scales::label_number(big.mark = ",",
                                                           drop0trailing = TRUE))
    } else {
      p <- p + scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE))
    }
    p
  }

  save_png(file.path(out_dir, "UMI_Count.png"),        dens("nCount_RNA", QC_MIN_COUNT),   width = 8, height = 5.5)
  save_png(file.path(out_dir, "nGenes.png"), dens("nFeature_RNA", QC_MIN_FEATURE), width = 8, height = 5.5)
  save_png(file.path(out_dir, "Complexity_Score.png"), dens("log10GenesPerUMI", QC_MIN_COMPLEX, logx = FALSE), width = 8, height = 5.5)
  save_png(file.path(out_dir, "Mito_Ratio.png"), dens("percent_mito", QC_MAX_MITO), width = 8, height = 5.5)
  save_png(file.path(out_dir, "Ribo_Ratio.png"), dens("percent_ribo", QC_MIN_RIBO), width = 8, height = 5.5)
  save_png(file.path(out_dir, "Heme_Ratio.png"), dens("percent_hb", QC_MAX_HB), width = 8, height = 5.5)
  save_png(file.path(out_dir, "Platelet_Ratio.png"), dens("percent_plat", QC_MAX_PLAT), width = 8, height = 5.5)

  p1 <- QC_Plots_Genes(seurat_object = seu, low_cutoff = QC_MIN_FEATURE,
                       high_cutoff = 5500, group.by = group_col)
  p2 <- QC_Plots_UMIs(seurat_object = seu, low_cutoff = 1200,
                      high_cutoff = 45000, group.by = group_col)
  p3 <- QC_Plots_Mito(seurat_object = seu, high_cutoff = 20, group.by = group_col)
  p4 <- QC_Plots_Complexity(seurat_object = seu, high_cutoff = QC_MIN_COMPLEX,
                            group.by = group_col)
  save_png(file.path(out_dir, "Grouped_Cutoff.png"),
           wrap_plots(p1, p2, p3, p4, ncol = 2) & theme_cp018(base_size = 14),
           width = 12, height = 9)

  save_png(file.path(out_dir, "UMIvsGene.png"), width = 8.5, height = 6.5,
           plot_obj = QC_Plot_UMIvsGene(seurat_object = seu,
                             low_cutoff_gene = QC_MIN_FEATURE, high_cutoff_gene = 5500,
                             low_cutoff_UMI = 500, high_cutoff_UMI = 50000,
                             group.by = group_col))

  save_png(file.path(out_dir, "MitovsGene.png"), width = 8.5, height = 6.5,
           plot_obj = QC_Plot_GenevsFeature(seurat_object = seu, feature1 = "percent_mito",
                                 low_cutoff_gene = QC_MIN_FEATURE, high_cutoff_gene = 5500,
                                 high_cutoff_feature = 20, group.by = group_col))
}

# CP018 grouping variables. CellType_Sort is dropped (constant); Replicate and
# Timepoint_Rep added — the whole point of the duplicate design.
group_vars <- c("Sample_Subfolder", "Timepoint", "Replicate", "Timepoint_Rep")

for (g in group_vars) run_qc_plots(seu_CP018, g, file.path(pre_root, g))

# -----------------------------
# 5) Filtering (thresholds from config; identical to CP003)
# -----------------------------
if ("JoinLayers" %in% getNamespaceExports("Seurat")) seu_CP018 <- JoinLayers(seu_CP018)
DefaultAssay(seu_CP018) <- "RNA"

filtered_CP018 <- subset(
  x = seu_CP018,
  subset = (nCount_RNA        >= QC_MIN_COUNT)   &
           (nFeature_RNA      >= QC_MIN_FEATURE) &
           (log10GenesPerUMI  >  QC_MIN_COMPLEX) &
           (percent_mito      <  QC_MAX_MITO)    &
           (percent_ribo      >  QC_MIN_RIBO)    &
           (percent_hb        <  QC_MAX_HB)      &
           (percent_plat      <  QC_MAX_PLAT)
)

# Keep genes detected in >= QC_MIN_CELLS_PER_GENE cells
counts_mat <- tryCatch(
  LayerData(filtered_CP018, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(filtered_CP018, assay = "RNA", slot = "counts")
)
if (is.null(counts_mat) || nrow(counts_mat) == 0) stop("Could not retrieve RNA counts matrix.")

genes_keep <- Matrix::rowSums(counts_mat > 0) >= QC_MIN_CELLS_PER_GENE
filtered_CP018 <- subset(filtered_CP018, features = names(genes_keep[genes_keep]))

# Filtering summary — how many cells each library lost
filt_summary <- merge(
  raw_counts,
  setNames(as.data.frame(table(filtered_CP018$Sample_Subfolder)),
           c("Sample_Subfolder", "n_cells_postQC")),
  by = "Sample_Subfolder", all = TRUE
)
filt_summary$pct_retained <- round(100 * filt_summary$n_cells_postQC /
                                     filt_summary$n_cells_raw, 1)
write.csv(filt_summary, file.path(QC_ROOT, "CP018_QC_filtering_summary.csv"),
          row.names = FALSE)
print(filt_summary)

for (g in group_vars) run_qc_plots(filtered_CP018, g, file.path(post_root, g))

# -----------------------------
# 6) Cell cycle scoring
# -----------------------------
filtered_CP018 <- JoinLayers(filtered_CP018)
seu <- filtered_CP018
DefaultAssay(seu) <- "RNA"
seu$Sample <- seu$Sample_Subfolder

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
               reduction.name = "umap.pre", verbose = FALSE)

s.genes.use   <- intersect(Seurat::cc.genes.updated.2019$s.genes,   rownames(seu))
g2m.genes.use <- intersect(Seurat::cc.genes.updated.2019$g2m.genes, rownames(seu))

seu <- CellCycleScoring(seu, s.features = s.genes.use,
                        g2m.features = g2m.genes.use, set.ident = TRUE)

ggsave(file.path(cc_dir, "UMAP_by_CellCyclePhase.png"),
       DimPlot2(seu, reduction = "umap.pre", group.by = "Phase", pt.size = 0.4),
       width = 10, height = 6, dpi = 600, bg = "white")

for (g in c("Sample", "Timepoint", "Replicate", "Timepoint_Rep")) {
  ggsave(file.path(cc_dir, paste0("CellCycle_Score_Vln_by_", g, ".png")),
         VlnPlot2(seu, features = c("S.Score", "G2M.Score"),
                  group.by = g, pt.size = 0.1),
         width = 10, height = 6, dpi = 600, bg = "white")
}

# Phase composition table — FACT-Seq sorts proliferating cells, so phase
# distribution is a biological readout, not just a QC metric.
phase_tab <- as.data.frame.matrix(table(seu$Sample_Subfolder, seu$Phase))
phase_tab$Sample_Subfolder <- rownames(phase_tab)
write.csv(phase_tab, file.path(cc_dir, "CP018_phase_by_library.csv"), row.names = FALSE)
print(phase_tab)

ridge_feats <- intersect(c("PCNA", "TOP2A", "MCM6", "MKI67"), rownames(seu))
if (length(ridge_feats) > 0) {
  save_png(file.path(cc_dir, "Ridge_CellCycle_Markers.png"),
           RidgePlot(seu, features = ridge_feats, ncol = 2) &
             theme_cp018(base_size = 14),
           width = 10, height = 7)
}

# -----------------------------
# 7) Azimuth reference annotation
# -----------------------------
seu <- RunAzimuth(seu, reference = "pbmcref")

ggsave(file.path(cc_dir, "UMAP_Azimuth_predicted_celltype_l2.png"),
       DimPlot2(seu, group.by = "predicted.celltype.l2",
                label = TRUE, repel = TRUE, box = TRUE),
       width = 10, height = 8, dpi = 600, bg = "white")

ggsave(file.path(cc_dir, "UMAP_Azimuth_predicted_celltype_l2_by_Sample.png"),
       DimPlot2(seu, group.by = "predicted.celltype.l2", split.by = "Sample",
                label = TRUE, repel = TRUE, box = TRUE),
       width = 16, height = 8, dpi = 600, bg = "white")

ggsave(file.path(cc_dir, "UMAP_Azimuth_predicted_celltype_l2_by_Timepoint.png"),
       DimPlot2(seu, group.by = "predicted.celltype.l2", split.by = "Timepoint",
                label = TRUE, repel = TRUE, box = TRUE),
       width = 14, height = 8, dpi = 600, bg = "white")

# -----------------------------
# 8) Doublet detection (per library)
# -----------------------------
split_seu <- SplitObject(seu, split.by = "Sample")
samples   <- names(split_seu)

dbl_rows <- list()
for (i in samples) {
  sce <- scDblFinder(GetAssayData(split_seu[[i]], assay = "RNA", layer = "counts"))
  split_seu[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seu[[i]]$scDblFinder.class <- sce$scDblFinder.class

  ggsave(file.path(dbl_dir, paste0(i, "_doublets.png")),
         DimPlot2(split_seu[[i]], reduction = "umap.pre",
                  group.by = "scDblFinder.class") + ggtitle(paste0(i, " Doublets")),
         device = "png", width = 6, height = 6, dpi = 600, bg = "white")

  dbl_rows[[i]] <- data.frame(
    Sample_Subfolder = i,
    n_cells   = ncol(split_seu[[i]]),
    n_doublet = sum(split_seu[[i]]$scDblFinder.class == "doublet"),
    pct_doublet = round(100 * mean(split_seu[[i]]$scDblFinder.class == "doublet"), 2)
  )
}
dbl_summary <- do.call(rbind, dbl_rows)
write.csv(dbl_summary, file.path(dbl_dir, "CP018_doublet_summary.csv"), row.names = FALSE)
print(dbl_summary)

for (i in samples) {
  split_seu[[i]] <- subset(split_seu[[i]], subset = scDblFinder.class == "singlet")
}

seu_singlets <- Merge_Seurat_List(split_seu, merge.data = TRUE)
seu_singlets <- JoinLayers(seu_singlets)
Idents(seu_singlets) <- "Sample"

cat("\nFinal singlet counts per library:\n")
print(table(seu_singlets$Sample_Subfolder))
cat("Total:", ncol(seu_singlets), "cells\n")

qs_save(seu_singlets, file = QS_SINGLETS)
message("Saved: ", QS_SINGLETS)
