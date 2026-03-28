################################################################################
# FIGURE 4-5 PREP: CD8 T Cell Sub-clustering, Annotation & All Analyses
#   — TARA Cohort, HEI only
#
# Workflow:
#   1. Load TARA_ALL → filter to HEI → subset CD8 + rescue
#   2. Reprocess: RNA → ADT → WNN → UMAP → clustering  [SAVE: _integrated.qs2]
#   3. Explore clusters (avg expression CSVs + heatmaps)
#   4. Annotate (fill in after inspection)
#   5. Save annotated object                            [SAVE: _annotated.qs2]
#   8A. DGE/DPE: naïve cluster comparisons
#   8B. DGE/DPE: expanding clones (Sup vs Unsup)
#   8C. Cluster proportions
#   8D. Clonal expansion plots
#   8E. Trajectory (Slingshot + Monocle3)               [SAVE: _with_trajectory.qs2]
#   8F. Expansion vs exhaustion control
#   8G. Viral load correlations
#   8H. Module scores (LAST — can edit/rerun independently)
#   9. Save final object                                [SAVE: _final.qs2]
#
# Checkpoints: To rerun module scores only, load _with_trajectory.qs2
#              and run from section 8H.
#
# Input:  TARA_ALL_sorted_refined.qs2
# Output: analysis/ with numbered subfolders + saved objects
# Figures: Figure3_Plots.R
#
# ANNOTATION UPDATE (v2):
#   Naïve CD8 (innate-like)  → Naïve CD8 2
#   Naïve CD8 (post-ART)     → Naïve CD8 3
#   Effector CD8              → TEM CD8
#   Terminal TEMRA CD8        → TEMRA CD8
#   TRDV1+ γδ T cell          → γδ1 T cell
#   Naïve-like TRDV1+ γδ      → Naïve γδ1 T cell
#   Vγ9Vδ2 γδ T cell          → γδ2 T cell
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(SeuratExtend)
library(scCustomize)
library(tidyr)
library(clustree)
library(cowplot)
library(patchwork)
library(pheatmap)
library(grid)
library(qs2)
library(RColorBrewer)
library(batchelor)
library(SeuratWrappers)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir    <- "~/Documents/CD8_Longitudinal"
saved_dir   <- file.path(base_dir, "saved_R_data")
out_dir     <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 4-5/"
analysis_dir <- file.path(out_dir, "analysis")

# Numbered analysis subfolders — clean organization
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

dir_01_integration <- file.path(analysis_dir, "01_integration_diagnostics")
dir_02_clustering  <- file.path(analysis_dir, "02_clustering_exploration")
dir_03_dge         <- file.path(analysis_dir, "03_DGE")
dir_04_dpe         <- file.path(analysis_dir, "04_DPE")
dir_05_proportions <- file.path(analysis_dir, "05_proportions")
dir_06_clonal      <- file.path(analysis_dir, "06_clonal_expansion")
dir_07_trajectory  <- file.path(analysis_dir, "07_trajectory")
dir_08_expansion   <- file.path(analysis_dir, "08_expansion_vs_exhaustion")
dir_09_viralload   <- file.path(analysis_dir, "09_viral_load")
dir_10_modules     <- file.path(analysis_dir, "10_module_scores")

for (d in ls(pattern = "^dir_\\d")) {
  dir.create(get(d), recursive = TRUE, showWarnings = FALSE)
}

cat("Analysis output structure created.\n")

# ── Load refined object ──────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted_refined.qs2"))

# ── Filter to HEI only ──────────────────────────────────────────────────────
# Remove HEU and HUU — focus on HIV-exposed infected infants across ART stages
TARA_ALL <- subset(TARA_ALL, Condition == "HEI")
cat("Filtered to HEI only:", ncol(TARA_ALL), "cells\n")
cat("Timepoint_Group distribution:\n")
print(table(TARA_ALL$Timepoint_Group))

################################################################################
# 1. SUBSET: CD8 clusters + rescue CD8+ cells from other clusters
################################################################################

# ── Define the 5 known CD8 clusters ──────────────────────────────────────────
cd8_cluster_labels <- c(
  "Tcm/Tscm CD8",
  "Naïve CD8",
  "TEMRA/CTL",
  "TRDV1+ γδ",
  "Tpex CD8"
)

# Flag cells from known CD8 clusters
in_cd8_cluster <- TARA_ALL$Manual_Annotation_refined %in% cd8_cluster_labels

# ── Diagnostic: CD8 ADT distribution in CD8 vs non-CD8 clusters ─────────────
# This helps pick a rescue threshold
DefaultAssay(TARA_ALL) <- "ADT"

# CD8 ADT feature name
cd8_adt_name <- "CD8A"

cd8_adt_vals  <- FetchData(TARA_ALL, vars = cd8_adt_name)[, 1]
diag_df <- data.frame(
  CD8_ADT   = cd8_adt_vals,
  is_cd8    = ifelse(in_cd8_cluster, "CD8 Clusters", "Other Clusters"),
  cluster   = as.character(TARA_ALL$Manual_Annotation_refined)
)

# Plot distribution
p_diag <- ggplot(diag_df, aes(x = CD8_ADT, fill = is_cd8)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "CD8 ADT Expression: CD8 clusters vs Others",
       subtitle = "Red dashed line = proposed rescue threshold (adjust as needed)",
       x = paste0(cd8_adt_name, " (DSB-normalized)"),
       y = "Density", fill = NULL) +
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(dir_01_integration, "Diagnostic_CD8_ADT_distribution.png"),
       plot = p_diag, width = 10, height = 6, dpi = 300, bg = "white")

cat("Diagnostic plot saved. Check the distribution and adjust threshold below.\n")

# ── Rescue: tight multi-gate ──────────────────────────────────────────────────
cd8a_threshold <- 2.0   # stringent CD8A positivity
cd3d_threshold <- 1.0   # must be a T cell
cd4_max        <- 1.0   # exclude CD4+ T cells
cd14_max       <- 0.5   # exclude monocytes
cd19_max       <- 0.5   # exclude B cells

# Fetch all gating markers at once
gate_data <- FetchData(TARA_ALL, vars = c("CD8A", "CD3D", "CD4", "CD14", "CD19"),
                       layer = "data")

# Rescue mask: NOT already in CD8 clusters, AND passes all gates
rescue_mask <- (!in_cd8_cluster) &
  (gate_data$CD8A > cd8a_threshold) &
  (gate_data$CD3D > cd3d_threshold) &
  (gate_data$CD4  < cd4_max) &
  (gate_data$CD14 < cd14_max) &
  (gate_data$CD19 < cd19_max)

cat("Cells in CD8 clusters:", sum(in_cd8_cluster), "\n")
cat("Cells rescued from other clusters:", sum(rescue_mask), "\n")
cat("  Rescued from:\n")
print(table(TARA_ALL$Manual_Annotation_refined[rescue_mask]))

# ── Diagnostic: scatter plots of gating markers for rescued cells ────────────
gate_data$is_rescued <- rescue_mask
gate_data$is_cd8     <- in_cd8_cluster
gate_data$group      <- case_when(
  in_cd8_cluster ~ "CD8 Clusters",
  rescue_mask    ~ "Rescued",
  TRUE           ~ "Excluded"
)

p_gate <- ggplot(gate_data %>% filter(group != "Excluded"),
                 aes(x = CD8A, y = CD3D, color = group)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_vline(xintercept = cd8a_threshold, linetype = "dashed", color = "red") +
  geom_hline(yintercept = cd3d_threshold, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("CD8 Clusters" = "grey70", "Rescued" = "#E41A1C")) +
  labs(title = "Rescue gate: CD8A vs CD3D (ADT)",
       subtitle = paste0("Rescued: ", sum(rescue_mask), " cells"),
       x = "CD8A (DSB)", y = "CD3D (DSB)", color = NULL) +
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(dir_01_integration, "Diagnostic_rescue_gate_CD8A_vs_CD3D.png"),
       plot = p_gate, width = 9, height = 7, dpi = 300, bg = "white")

cat("Rescue gate diagnostic plot saved.\n")

# ── Combine and subset ───────────────────────────────────────────────────────
keep_cells <- in_cd8_cluster | rescue_mask
TARA_cd8 <- subset(TARA_ALL, cells = colnames(TARA_ALL)[keep_cells])

cat("\nTotal CD8 subset:", ncol(TARA_cd8), "cells\n")

################################################################################
# 2. REPROCESS + INTEGRATE: RNA → ADT → WNN → UMAP → Clustering
#    Uses FastMNN integration for both RNA and ADT to correct batch effects
################################################################################

# ── Reproducibility: set global seed ─────────────────────────────────────────
SEED <- 42
set.seed(SEED)
cat("Global seed set to", SEED, "for reproducible UMAP + clustering\n")

# ── RNA: Split by sample → Normalize → Variable Features → Scale → PCA ──────
DefaultAssay(TARA_cd8) <- "RNA"

# Split RNA layers by sample (batch) for integration
TARA_cd8[["RNA"]] <- split(TARA_cd8[["RNA"]], f = TARA_cd8$orig.ident)

TARA_cd8 <- NormalizeData(TARA_cd8, assay = "RNA")
TARA_cd8 <- FindVariableFeatures(TARA_cd8, assay = "RNA",
                                 selection.method = "vst", nfeatures = 3000)
TARA_cd8 <- ScaleData(TARA_cd8, assay = "RNA")
TARA_cd8 <- RunPCA(TARA_cd8, assay = "RNA", reduction.name = "pca")

# ── RNA Integration: FastMNN ─────────────────────────────────────────────────
TARA_cd8 <- IntegrateLayers(
  TARA_cd8,
  method        = FastMNNIntegration,
  assay         = "RNA",
  new.reduction = "integrated.mnn.rna"
)

# Rejoin RNA layers after integration (needed for downstream operations)
TARA_cd8[["RNA"]] <- JoinLayers(TARA_cd8[["RNA"]])

# ── Diagnostic: RNA integration — orig.ident + Timepoint_Group ────────────────
TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "pca", dims = 1:30, seed.use = SEED,
                    reduction.name = "umap.unintegrated", reduction.key = "unintUMAP_")

TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "integrated.mnn.rna", dims = 1:30, seed.use = SEED,
                    reduction.name = "umap.mnn.rna", reduction.key = "mnnUMAP_")

# By orig.ident
p_rna_unint_sample <- DimPlot2(TARA_cd8, reduction = "umap.unintegrated",
                               group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("RNA Unintegrated — by Sample")

p_rna_int_sample <- DimPlot2(TARA_cd8, reduction = "umap.mnn.rna",
                             group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("RNA FastMNN Integrated — by Sample")

ggsave(file.path(dir_01_integration, "CD8_RNA_integration_by_sample.png"),
       plot = p_rna_unint_sample | p_rna_int_sample,
       width = 22, height = 8, dpi = 300, bg = "white")

# By Timepoint_Group
p_rna_unint_tp <- DimPlot2(TARA_cd8, reduction = "umap.unintegrated",
                           group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("RNA Unintegrated — by Timepoint Group")

p_rna_int_tp <- DimPlot2(TARA_cd8, reduction = "umap.mnn.rna",
                         group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("RNA FastMNN Integrated — by Timepoint Group")

ggsave(file.path(dir_01_integration, "CD8_RNA_integration_by_timepoint_group.png"),
       plot = p_rna_unint_tp | p_rna_int_tp,
       width = 22, height = 8, dpi = 300, bg = "white")

cat("RNA integration comparison saved (by sample + by timepoint group).\n")

# ── ADT: Split by sample → Scale → PCA → Integration (FastMNN) ──────────────
DefaultAssay(TARA_cd8) <- "ADT"

# Split ADT layers by sample for integration
TARA_cd8[["ADT"]] <- split(TARA_cd8[["ADT"]], f = TARA_cd8$orig.ident)

VariableFeatures(TARA_cd8) <- rownames(TARA_cd8[["ADT"]])
TARA_cd8 <- ScaleData(TARA_cd8, assay = "ADT",
                      features = VariableFeatures(TARA_cd8), verbose = FALSE)
TARA_cd8 <- RunPCA(TARA_cd8, assay = "ADT",
                   features = VariableFeatures(TARA_cd8),
                   reduction.name = "apca")

# ADT Integration: FastMNN
TARA_cd8 <- IntegrateLayers(
  TARA_cd8,
  method        = FastMNNIntegration,
  assay         = "ADT",
  new.reduction = "integrated.mnn.adt"
)

# Rejoin ADT layers after integration
TARA_cd8[["ADT"]] <- JoinLayers(TARA_cd8[["ADT"]])

# ── Diagnostic: ADT integration — orig.ident + Timepoint_Group ────────────────
TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "apca", dims = 1:20, seed.use = SEED,
                    reduction.name = "umap.adt.unint", reduction.key = "adtunintUMAP_")

TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "integrated.mnn.adt", dims = 1:20, seed.use = SEED,
                    reduction.name = "umap.mnn.adt", reduction.key = "mnnADTUMAP_")

# By orig.ident
p_adt_unint_sample <- DimPlot2(TARA_cd8, reduction = "umap.adt.unint",
                               group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("ADT Unintegrated — by Sample")

p_adt_int_sample <- DimPlot2(TARA_cd8, reduction = "umap.mnn.adt",
                             group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("ADT FastMNN Integrated — by Sample")

ggsave(file.path(dir_01_integration, "CD8_ADT_integration_by_sample.png"),
       plot = p_adt_unint_sample | p_adt_int_sample,
       width = 22, height = 8, dpi = 300, bg = "white")

# By Timepoint_Group
p_adt_unint_tp <- DimPlot2(TARA_cd8, reduction = "umap.adt.unint",
                           group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("ADT Unintegrated — by Timepoint Group")

p_adt_int_tp <- DimPlot2(TARA_cd8, reduction = "umap.mnn.adt",
                         group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("ADT FastMNN Integrated — by Timepoint Group")

ggsave(file.path(dir_01_integration, "CD8_ADT_integration_by_timepoint_group.png"),
       plot = p_adt_unint_tp | p_adt_int_tp,
       width = 22, height = 8, dpi = 300, bg = "white")

cat("ADT integration comparison saved (by sample + by timepoint group).\n")

# ── WNN: Combine integrated RNA + integrated ADT ────────────────────────────
TARA_cd8 <- FindMultiModalNeighbors(
  TARA_cd8,
  reduction.list = list("integrated.mnn.rna", "integrated.mnn.adt"),
  dims.list      = list(1:30, 1:20)
)

# ── UMAP on WNN graph ───────────────────────────────────────────────────────
TARA_cd8 <- RunUMAP(TARA_cd8,
                    nn.name        = "weighted.nn",
                    reduction.name = "wnn.umap",
                    reduction.key  = "wnnUMAP_",
                    seed.use       = SEED)

# ── Multi-resolution clustering for exploration ──────────────────────────────
resolutions <- seq(0.2, 2.0, by = 0.2)

for (res in resolutions) {
  TARA_cd8 <- FindClusters(
    TARA_cd8,
    graph.name  = "wsnn",
    resolution  = res,
    algorithm   = 3,   # Leiden
    random.seed = SEED,
    verbose     = FALSE
  )
}

# Fix factor levels so they sort numerically
for (col in grep("^wsnn_res\\.", colnames(TARA_cd8@meta.data), value = TRUE)) {
  lvls <- as.character(sort(as.numeric(levels(TARA_cd8[[col]][, 1]))))
  TARA_cd8[[col]] <- factor(TARA_cd8[[col]][, 1], levels = lvls)
}

# ── Clustree: explore resolution stability ───────────────────────────────────
p_clustree <- clustree(TARA_cd8, prefix = "wsnn_res.")

ggsave(file.path(dir_02_clustering, "CD8_clustree.png"),
       plot = p_clustree, width = 15, height = 9, dpi = 300, bg = "white")

p_clustree_stab <- clustree(TARA_cd8, prefix = "wsnn_res.",
                            node_colour = "sc3_stability")

ggsave(file.path(dir_02_clustering, "CD8_clustree_stability.png"),
       plot = p_clustree_stab, width = 15, height = 9, dpi = 300, bg = "white")

cat("Clustree saved (standard + stability-colored).\n")

# ── Pick resolution ──────────────────────────────────────────────────────────
chosen_res <- 0.4
TARA_cd8$seurat_clusters <- TARA_cd8[[paste0("wsnn_res.", chosen_res)]][, 1]
Idents(TARA_cd8) <- "seurat_clusters"

cat("Chosen resolution:", chosen_res, "\n")
cat("Number of clusters:", length(levels(Idents(TARA_cd8))), "\n")
print(table(Idents(TARA_cd8)))

# ── UMAP at each resolution for comparison ───────────────────────────────────
for (res in c(0.4, 0.6, 0.8, 1.0)) {
  p <- DimPlot2(
    TARA_cd8,
    reduction = "wnn.umap",
    group.by  = paste0("wsnn_res.", res),
    cols      = "light",
    label     = TRUE,
    box       = TRUE,
    repel     = TRUE,
    label.color = "black",
    pt.size   = 0.6
  ) + ggtitle(paste0("Resolution: ", res))
  
  ggsave(file.path(dir_02_clustering, paste0("CD8_UMAP_res", res, ".png")),
         plot = p, width = 8, height = 7, dpi = 300, bg = "white")
}

# ── CHECKPOINT: Save after integration ───────────────────────────────────────
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_integrated.qs2"))
cat("CHECKPOINT: Integrated object saved as TARA_cd8_HEI_integrated.qs2\n")

################################################################################
# 3. ANNOTATION EXPLORATION: Average expression per cluster
################################################################################

# ── Key markers for CD8 sub-population identification ────────────────────────
annotation_markers_rna <- c(
  # Naïve / Stem-like
  "CCR7", "SELL", "TCF7", "LEF1", "IL7R",
  # Central memory
  "CD27", "CD28", "BCL2", "BACH2",
  # Effector memory
  "GZMK", "EOMES", "CXCR3",
  # TEMRA / CTL
  "GZMB", "GNLY", "PRF1", "NKG7", "TBX21", "CX3CR1",
  # Exhaustion / Tpex
  "TOX", "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CXCL13",
  # Proliferation
  "MKI67", "TOP2A", "STMN1",
  # Tissue residency
  "ITGAE", "CXCR6", "CD69",
  # γδ TCR
  "TRDV1", "TRGV9", "TRDC",
  # NK-like
  "TYROBP", "KLRD1", "FCGR3A", "NCAM1",
  # Activation
  "HLA-DRA", "CD38"
)

annotation_markers_adt <- c(
  # Naïve
  "CD45RA", "SELL", "CD7",
  # Memory
  "CD45RO", "IL7R", "CD27", "CD28",
  # Effector
  "B3GAT1", "KLRG1", "KIR3DL1",
  # Exhaustion
  "TIGIT", "PDCD1", "LAG3",
  # Activation
  "CD38", "ENTPD1", "NT5E",
  # NK-like
  "NCAM1", "FCGR3A", "SIGLEC7",
  # TCR
  "TCR-AB", "TCR-vA7.2", "TCR-vD2"
)

# ── Compute average expression per cluster ───────────────────────────────────
DefaultAssay(TARA_cd8) <- "RNA"
avg_rna <- AverageExpression(
  TARA_cd8,
  assays   = "RNA",
  features = annotation_markers_rna,
  group.by = "seurat_clusters",
  slot     = "data"
)$RNA

DefaultAssay(TARA_cd8) <- "ADT"
avg_adt <- AverageExpression(
  TARA_cd8,
  assays   = "ADT",
  features = annotation_markers_adt,
  group.by = "seurat_clusters",
  slot     = "data"
)$ADT

# ── Quick heatmap for exploration (not final figure) ─────────────────────────
scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)))
}

avg_rna_scaled <- scale_01(avg_rna)
avg_adt_scaled <- scale_01(log10(avg_adt + 1))

colnames(avg_rna_scaled) <- gsub("^g\\s*", "", colnames(avg_rna_scaled))
colnames(avg_adt_scaled) <- gsub("^g\\s*", "", colnames(avg_adt_scaled))

col_order <- as.character(sort(as.numeric(colnames(avg_rna_scaled))))
avg_rna_scaled <- avg_rna_scaled[, col_order]
avg_adt_scaled <- avg_adt_scaled[, col_order]

heatmap_colors <- colorRampPalette(c(
  "#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"
))(100)

png(file.path(dir_02_clustering, "CD8_AvgExpression_RNA_exploration.png"),
    width = 14, height = 16, units = "in", res = 300, bg = "white")
pheatmap(
  avg_rna_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color        = heatmap_colors,
  border_color = "white",
  cellwidth    = 40,
  cellheight   = 14,
  fontsize     = 10,
  main         = "RNA Average Expression (scaled) per CD8 Sub-cluster"
)
dev.off()

png(file.path(dir_02_clustering, "CD8_AvgExpression_ADT_exploration.png"),
    width = 14, height = 12, units = "in", res = 300, bg = "white")
pheatmap(
  avg_adt_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color        = heatmap_colors,
  border_color = "white",
  cellwidth    = 40,
  cellheight   = 14,
  fontsize     = 10,
  main         = "ADT Average Expression (scaled) per CD8 Sub-cluster"
)
dev.off()

cat("\n=== RNA Average Expression (top markers per cluster) ===\n")
print(round(avg_rna_scaled, 2))

cat("\n=== ADT Average Expression (top markers per cluster) ===\n")
print(round(avg_adt_scaled, 2))

################################################################################
# 4. ANNOTATION: Inspect clusters, then annotate
################################################################################

cat("\n=== Timepoint_Group distribution per cluster ===\n")
print(round(prop.table(table(TARA_cd8$seurat_clusters, TARA_cd8$Timepoint_Group), margin = 1), 3))

cat("\n=== Cluster sizes ===\n")
print(table(TARA_cd8$seurat_clusters))

cat("\n=== has_TCR per cluster ===\n")
print(round(prop.table(table(TARA_cd8$seurat_clusters, TARA_cd8$has_TCR), margin = 1), 3))

# ── Export avg expression CSVs for annotation ────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"
avg_adt_explore <- AverageExpression(TARA_cd8, assays = "ADT",
                                     features = rownames(TARA_cd8[["ADT"]]),
                                     group.by = "seurat_clusters", slot = "data")$ADT
colnames(avg_adt_explore) <- gsub("^g\\s*", "", colnames(avg_adt_explore))
write.csv(round(avg_adt_explore, 3),
          file.path(dir_02_clustering, "CD8_HEI_avg_ADT_ALL_by_cluster.csv"))

DefaultAssay(TARA_cd8) <- "RNA"
rna_markers_expanded <- c(
  "CCR7", "SELL", "TCF7", "LEF1", "IL7R", "FOXP1", "KLF2", "S1PR1",
  "CD27", "CD28", "BCL2", "BACH2", "CD44", "ID3",
  "FAS", "IL2RB", "CXCR3", "CXCR4",
  "GZMK", "EOMES", "GZMB", "GNLY", "PRF1", "NKG7", "TBX21", "CX3CR1",
  "FGFBP2", "GZMA", "GZMH", "GZMM",
  "TOX", "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "ENTPD1",
  "MKI67", "TOP2A", "STMN1",
  "ITGAE", "CXCR6", "CD69", "ZNF683", "ITGA1",
  "TRDV1", "TRDV2", "TRGV9", "TRDC",
  "TYROBP", "KLRD1", "KLRB1", "KLRC1", "NCAM1", "FCGR3A",
  "KIR2DL3", "KIR3DL1", "KIR3DL2",
  "HLA-DRA", "CD38", "IFNG", "TNF",
  "SLC4A10", "ZBTB16", "RORC", "NCR3",
  "RUNX3", "ZEB2", "PRDM1", "ID2", "S1PR5",
  "FOXP3", "IL2RA"
)
rna_markers_present <- rna_markers_expanded[rna_markers_expanded %in% rownames(TARA_cd8[["RNA"]])]

avg_rna_explore <- AverageExpression(TARA_cd8, assays = "RNA",
                                     features = rna_markers_present,
                                     group.by = "seurat_clusters", slot = "data")$RNA
colnames(avg_rna_explore) <- gsub("^g\\s*", "", colnames(avg_rna_explore))
write.csv(round(avg_rna_explore, 3),
          file.path(dir_02_clustering, "CD8_HEI_avg_RNA_expanded_by_cluster.csv"))

cat("\nExploration CSVs saved. Upload these for annotation help.\n")

################################################################################
# 4. ANNOTATION: HEI-only CD8, Resolution 0.4 (13 clusters: X0–X12)
#    UPDATED ANNOTATIONS (v2)
################################################################################

# ── Step 1: Remove contaminants ──────────────────────────────────────────────
clusters_to_remove <- c("11", "12")  # X11=APC/DC, X12=tiny suspect cluster
TARA_cd8 <- subset(TARA_cd8, seurat_clusters %in% clusters_to_remove, invert = TRUE)
cat("Removed clusters:", paste(clusters_to_remove, collapse = ", "), "\n")
cat("Remaining cells:", ncol(TARA_cd8), "\n")

# ── Step 2: Initial annotation (UPDATED v2) ─────────────────────────────────
cd8_annotations <- c(
  "0"  = "Naïve CD8",
  "1"  = "TEM CD8",
  "2"  = "Naïve CD8 2",
  "3"  = "γδ1 T cell",
  "4"  = "TEMRA CD8",
  "5"  = "Naïve CD8 3",
  "6"  = "Transitional Tem CD8",
  "7"  = "Naïve γδ1 T cell",          # will be split
  "8"  = "KIR+ innate-like CD8",
  "9"  = "γδ2 T cell",
  "10" = "MAIT-like Trm"
)

TARA_cd8$seurat_clusters <- droplevels(TARA_cd8$seurat_clusters)
cat("Remaining cluster levels:", paste(levels(TARA_cd8$seurat_clusters), collapse = ", "), "\n")

annot_vec <- cd8_annotations[as.character(TARA_cd8$seurat_clusters)]
names(annot_vec) <- colnames(TARA_cd8)
TARA_cd8 <- AddMetaData(TARA_cd8, metadata = annot_vec, col.name = "CD8_Annotation")

# ── Step 3: TCR-based splits ─────────────────────────────────────────────────

# X3 (γδ1 T cell): move αβ TCR contaminants → TEM CD8
c3_ab <- TARA_cd8$seurat_clusters == "3" & TARA_cd8$has_TCR == TRUE
cat("\nX3 clean: moving", sum(c3_ab), "αβ TCR contaminants to TEM CD8\n")
TARA_cd8$CD8_Annotation[c3_ab] <- "TEM CD8"

# X7 (Naïve γδ1 T cell): αβ TCR cells → Naïve CD8, rest stays γδ
c7_ab <- TARA_cd8$seurat_clusters == "7" & TARA_cd8$has_TCR == TRUE
cat("X7 split: moving", sum(c7_ab), "αβ TCR cells to Naïve CD8\n")
TARA_cd8$CD8_Annotation[c7_ab] <- "Naïve CD8"

# ── Step 4: Finalize ─────────────────────────────────────────────────────────
TARA_cd8$CD8_Annotation <- factor(TARA_cd8$CD8_Annotation)
Idents(TARA_cd8) <- "CD8_Annotation"

cat("\n=== Final CD8 sub-cluster sizes (HEI only) ===\n")
print(table(TARA_cd8$CD8_Annotation))

cat("\n=== αβ TCR proportions after cleanup ===\n")
print(round(prop.table(table(TARA_cd8$CD8_Annotation, TARA_cd8$has_TCR), margin = 1), 3))

cat("\n=== Timepoint_Group per annotation ===\n")
print(round(prop.table(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 1), 3))

# ── Canonical cluster order (biological progression) — UPDATED ───────────────
col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
  "Transitional Tem CD8", "TEM CD8", "TEMRA CD8",
  "KIR+ innate-like CD8", "MAIT-like Trm",
  "γδ1 T cell", "Naïve γδ1 T cell", "γδ2 T cell"
)

################################################################################
# 5. SAVE ANNOTATED OBJECT
################################################################################
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_annotated.qs2"))
cat("Annotated HEI CD8 object saved.\n")

################################################################################
# 8. DOWNSTREAM ANALYSES — organized into numbered subfolders
################################################################################

dge_dir   <- dir_03_dge
dpe_dir   <- dir_04_dpe
prop_dir  <- dir_05_proportions
clone_dir <- dir_06_clonal

# ── Shared aesthetics (used across multiple sections) ────────────────────────
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

# UPDATED: effector cluster names
effector_clusters <- c("TEM CD8", "TEMRA CD8", "Transitional Tem CD8")

################################################################################
# 8A. DGE/DPE: Pairwise naïve cluster comparisons (UPDATED names)
################################################################################
cat("\n=== Running DGE/DPE (MAST) for naïve cluster comparisons ===\n")

Idents(TARA_cd8) <- "CD8_Annotation"

# UPDATED naïve cluster names
naive_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3")

naive_pairs <- combn(naive_clusters, 2, simplify = FALSE)
for (pair in naive_pairs) {
  safe_name <- paste0(gsub("[^A-Za-z0-9]", "", pair[1]), "_vs_",
                      gsub("[^A-Za-z0-9]", "", pair[2]))
  cat("  DGE:", pair[1], "vs", pair[2], "...\n")
  
  tryCatch({
    DefaultAssay(TARA_cd8) <- "RNA"
    dge <- FindMarkers(TARA_cd8, ident.1 = pair[1], ident.2 = pair[2],
                       test.use = "MAST", logfc.threshold = 0.25,
                       min.pct = 0.1, verbose = FALSE)
    write.csv(dge, file.path(dge_dir, paste0("DGE_MAST_RNA_", safe_name, ".csv")))
    
    DefaultAssay(TARA_cd8) <- "ADT"
    dpe <- FindMarkers(TARA_cd8, ident.1 = pair[1], ident.2 = pair[2],
                       test.use = "MAST", logfc.threshold = 0.1,
                       min.pct = 0.05, verbose = FALSE)
    write.csv(dpe, file.path(dpe_dir, paste0("DPE_MAST_ADT_", safe_name, ".csv")))
  }, error = function(e) cat("    ERROR:", e$message, "\n"))
}

# ── TEM CD8 vs TEMRA CD8 (UPDATED names) ────────────────────────────────────
cat("  DGE: TEM CD8 vs TEMRA CD8...\n")
tryCatch({
  DefaultAssay(TARA_cd8) <- "RNA"
  dge_eff <- FindMarkers(TARA_cd8, ident.1 = "TEM CD8",
                         ident.2 = "TEMRA CD8",
                         test.use = "MAST", logfc.threshold = 0.25,
                         min.pct = 0.1, verbose = FALSE)
  write.csv(dge_eff, file.path(dge_dir, "DGE_MAST_RNA_TEMCD8_vs_TEMRACD8.csv"))
  
  DefaultAssay(TARA_cd8) <- "ADT"
  dpe_eff <- FindMarkers(TARA_cd8, ident.1 = "TEM CD8",
                         ident.2 = "TEMRA CD8",
                         test.use = "MAST", logfc.threshold = 0.1,
                         min.pct = 0.05, verbose = FALSE)
  write.csv(dpe_eff, file.path(dpe_dir, "DPE_MAST_ADT_TEMCD8_vs_TEMRACD8.csv"))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

cat("  Pairwise DGE/DPE saved.\n")

################################################################################
# 8B. DGE/DPE: Expanding clones — PostART_Suppressed vs PostART_Unsuppressed
################################################################################
cat("\n=== Running DGE/DPE: Expanding clones Suppressed vs Unsuppressed ===\n")

# ── Helper function for expanding clone DGE/DPE ─────────────────────────────
run_expand_dge <- function(obj, cluster_name, label, dge_dir, dpe_dir) {
  if (!is.null(cluster_name)) {
    cells <- subset(obj,
                    CD8_Annotation == cluster_name &
                      !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)" &
                      Timepoint_Group %in% c("PostART_Suppressed", "PostART_Unsuppressed")
    )
  } else {
    cells <- subset(obj,
                    !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)" &
                      Timepoint_Group %in% c("PostART_Suppressed", "PostART_Unsuppressed")
    )
  }
  Idents(cells) <- "Timepoint_Group"
  
  n_sup   <- sum(Idents(cells) == "PostART_Suppressed")
  n_unsup <- sum(Idents(cells) == "PostART_Unsuppressed")
  cat(sprintf("    %s: Suppressed=%d, Unsuppressed=%d\n", label, n_sup, n_unsup))
  
  if (min(n_sup, n_unsup) < 10) {
    cat("    Too few cells, skipping.\n")
    return(invisible(NULL))
  }
  
  safe_label <- gsub("[^A-Za-z0-9]", "", label)
  
  tryCatch({
    DefaultAssay(cells) <- "RNA"
    dge <- FindMarkers(cells, ident.1 = "PostART_Suppressed",
                       ident.2 = "PostART_Unsuppressed", test.use = "MAST",
                       logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
    write.csv(dge, file.path(dge_dir, paste0("DGE_MAST_Expanding_", safe_label, "_Sup_vs_Unsup.csv")))
    
    DefaultAssay(cells) <- "ADT"
    dpe <- FindMarkers(cells, ident.1 = "PostART_Suppressed",
                       ident.2 = "PostART_Unsuppressed", test.use = "MAST",
                       logfc.threshold = 0.1, min.pct = 0.05, verbose = FALSE)
    write.csv(dpe, file.path(dpe_dir, paste0("DPE_MAST_ADT_Expanding_", safe_label, "_Sup_vs_Unsup.csv")))
    cat("    Saved.\n")
  }, error = function(e) cat("    ERROR:", e$message, "\n"))
}

# All expanding clones combined
cat("  All expanding clones...\n")
run_expand_dge(TARA_cd8, NULL, "AllClusters", dge_dir, dpe_dir)

# UPDATED cluster names
cat("  Within TEM CD8...\n")
run_expand_dge(TARA_cd8, "TEM CD8", "TEMCD8", dge_dir, dpe_dir)

cat("  Within TEMRA CD8...\n")
run_expand_dge(TARA_cd8, "TEMRA CD8", "TEMRACD8", dge_dir, dpe_dir)

cat("  Within Transitional Tem CD8...\n")
run_expand_dge(TARA_cd8, "Transitional Tem CD8", "TransitionalTem", dge_dir, dpe_dir)

cat("  Expanding clone DGE/DPE saved.\n")

################################################################################
# 8C. CLUSTER PROPORTIONS: by Timepoint_Group + by Patient (HEI only)
################################################################################
cat("\n=== Computing cluster proportions ===\n")

prop_tp <- as.data.frame(prop.table(
  table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 2
))
colnames(prop_tp) <- c("Cluster", "Timepoint_Group", "Proportion")
write.csv(prop_tp, file.path(prop_dir, "ClusterProportions_by_TimepointGroup.csv"),
          row.names = FALSE)

p_prop_tp <- ggplot(prop_tp, aes(x = Timepoint_Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion", fill = "CD8 Sub-cluster",
       title = "CD8 Sub-cluster Composition by ART Status (HEI only)") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x      = element_text(size = 11, angle = 35, hjust = 1),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(prop_dir, "ClusterProportions_by_TimepointGroup.png"),
       plot = p_prop_tp, width = 12, height = 8, dpi = 300, bg = "white")

# ── By Patient ───────────────────────────────────────────────────────────────
pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
names(pid_vec) <- colnames(TARA_cd8)
TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")

prop_pid <- as.data.frame(prop.table(
  table(TARA_cd8$CD8_Annotation, TARA_cd8$PID), margin = 2
))
colnames(prop_pid) <- c("Cluster", "PID", "Proportion")
write.csv(prop_pid, file.path(prop_dir, "ClusterProportions_by_Patient.csv"),
          row.names = FALSE)

p_prop_pid <- ggplot(prop_pid, aes(x = PID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion", fill = "CD8 Sub-cluster",
       title = "CD8 Sub-cluster Composition by Patient (HEI)") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x      = element_text(size = 11, angle = 35, hjust = 1),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(prop_dir, "ClusterProportions_by_Patient.png"),
       plot = p_prop_pid, width = 12, height = 8, dpi = 300, bg = "white")

write.csv(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group),
          file.path(prop_dir, "ClusterCounts_by_TimepointGroup.csv"))
write.csv(table(TARA_cd8$CD8_Annotation, TARA_cd8$PID),
          file.path(prop_dir, "ClusterCounts_by_Patient.csv"))

cat("  Proportion plots and CSVs saved to:", prop_dir, "\n")

################################################################################
# 8D. CLONAL EXPANSION: UMAP overlay + per-cluster cloneSize distribution
################################################################################
cat("\n=== Plotting clonal expansion ===\n")

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

clone_size_levels <- c("Single (0 < X <= 1)", "Small (1 < X <= 5)",
                       "Medium (5 < X <= 20)", "Large (20 < X <= 100)",
                       "Hyperexpanded (100 < X <= 500)")
clone_size_colors <- setNames(c(colorblind_vector[c(1, 3, 4, 5, 7)]), clone_size_levels)

p_clone_umap <- DimPlot_scCustom(
  TARA_cd8,
  group.by   = "cloneSize",
  reduction  = "wnn.umap",
  pt.size    = 0.5
) +
  scale_color_manual(values = clone_size_colors) +
  ggtitle("Clonal Expansion — CD8 Sub-clusters (HEI)") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(clone_dir, "CD8_UMAP_cloneSize.png"),
       plot = p_clone_umap, width = 10, height = 8, dpi = 300, bg = "white")

p_clone_split <- DimPlot_scCustom(
  TARA_cd8,
  group.by     = "cloneSize",
  reduction    = "wnn.umap",
  split.by     = "Timepoint_Group",
  split_seurat = TRUE,
  pt.size      = 0.4
) +
  scale_color_manual(values = clone_size_colors) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(clone_dir, "CD8_UMAP_cloneSize_byTimepointGroup.png"),
       plot = p_clone_split, width = 18, height = 6, dpi = 300, bg = "white")

hyper_label <- ifelse(
  !is.na(TARA_cd8$cloneSize) & TARA_cd8$cloneSize == "Hyperexpanded (100 < X <= 500)",
  "Hyperexpanded", "Other"
)
names(hyper_label) <- colnames(TARA_cd8)
TARA_cd8 <- AddMetaData(TARA_cd8, metadata = hyper_label, col.name = "is_hyperexpanded")

p_hyper <- DimPlot_scCustom(
  TARA_cd8,
  group.by  = "is_hyperexpanded",
  reduction = "wnn.umap",
  pt.size   = 0.5,
  colors_use = c("Other" = "grey85", "Hyperexpanded" = "#E41A1C"),
  order      = c("Hyperexpanded", "Other")
) +
  ggtitle("Hyperexpanded Clones — CD8 Sub-clusters (HEI)") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(clone_dir, "CD8_UMAP_hyperexpanded.png"),
       plot = p_hyper, width = 10, height = 8, dpi = 300, bg = "white")

# ── CloneSize distribution per sub-cluster ───────────────────────────────────
clone_meta <- TARA_cd8@meta.data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(CD8_Annotation, cloneSize) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CD8_Annotation) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

clone_size_order <- c("Single (0 < X <= 1)", "Small (1 < X <= 5)",
                      "Medium (5 < X <= 20)", "Large (20 < X <= 100)",
                      "Hyperexpanded (100 < X <= 500)")
clone_colors <- setNames(c(colorblind_vector[c(1, 3, 4, 5, 7)]), clone_size_order)
clone_meta$cloneSize <- factor(clone_meta$cloneSize, levels = clone_size_order)

# UPDATED cluster order
clone_meta$CD8_Annotation <- factor(clone_meta$CD8_Annotation, levels = col_order_cd8)

p_clone_bar <- ggplot(clone_meta, aes(x = CD8_Annotation, y = n, fill = cloneSize)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.75, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = clone_colors, name = "Clone Size", drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of Cells") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x       = element_text(size = 10, angle = 40, hjust = 1),
    legend.text       = element_text(size = 9),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(clone_dir, "CD8_CloneSize_per_cluster.png"),
       plot = p_clone_bar, width = 14, height = 8, dpi = 300, bg = "white")

clone_tp_meta <- TARA_cd8@meta.data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(CD8_Annotation, Timepoint_Group, cloneSize) %>%
  summarise(n = n(), .groups = "drop")

write.csv(clone_tp_meta, file.path(clone_dir, "CD8_CloneSize_per_cluster_per_timepoint.csv"),
          row.names = FALSE)

write.csv(clone_meta, file.path(clone_dir, "CD8_CloneSize_per_cluster.csv"),
          row.names = FALSE)

cat("  Clonal expansion plots saved to:", clone_dir, "\n")

################################################################################
# 8E. TRAJECTORY ANALYSIS: Slingshot pseudotime on effector CD8 sub-clusters
################################################################################
cat("\n=== Running trajectory analysis (Slingshot) ===\n")

if (!requireNamespace("slingshot", quietly = TRUE)) BiocManager::install("slingshot")
if (!requireNamespace("tradeSeq", quietly = TRUE)) BiocManager::install("tradeSeq")
library(slingshot)

traj_dir <- dir_07_trajectory

# ── Subset to αβ CD8 clusters (exclude γδ, MAIT — different lineage) — UPDATED
ab_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
                 "Transitional Tem CD8", "TEM CD8", "TEMRA CD8",
                 "KIR+ innate-like CD8")
traj_cells <- subset(TARA_cd8, CD8_Annotation %in% ab_clusters)
cat("  Trajectory cells:", ncol(traj_cells), "\n")
cat("  Clusters included:", paste(ab_clusters, collapse = ", "), "\n")

# ── Extract WNN UMAP coordinates and cluster labels ──────────────────────────
umap_coords <- Embeddings(traj_cells, reduction = "wnn.umap")
cluster_labels <- traj_cells$CD8_Annotation

# ── Run Slingshot ────────────────────────────────────────────────────────────
cat("  Running Slingshot (start = Naïve CD8)...\n")
sling <- slingshot(
  data       = umap_coords,
  clusterLabels = cluster_labels,
  start.clus = "Naïve CD8",
  stretch    = 0
)

pseudotime_mat <- slingPseudotime(sling)
cat("  Lineages found:", ncol(pseudotime_mat), "\n")
for (i in seq_len(ncol(pseudotime_mat))) {
  cat("    Lineage", i, ":", paste(slingLineages(sling)[[i]], collapse = " → "), "\n")
}

pt_min <- apply(pseudotime_mat, 1, min, na.rm = TRUE)
pt_min[is.infinite(pt_min)] <- NA

pt_vec <- pt_min
names(pt_vec) <- colnames(traj_cells)
traj_cells <- AddMetaData(traj_cells, metadata = pt_vec, col.name = "pseudotime")

# ── PLOT 1: UMAP colored by pseudotime ───────────────────────────────────────
umap_df <- data.frame(
  UMAP_1      = umap_coords[, 1],
  UMAP_2      = umap_coords[, 2],
  pseudotime  = traj_cells$pseudotime,
  Cluster     = traj_cells$CD8_Annotation,
  Timepoint   = traj_cells$Timepoint_Group
)

curves_list <- slingCurves(sling)
curve_dfs <- lapply(seq_along(curves_list), function(i) {
  crv <- curves_list[[i]]$s[curves_list[[i]]$ord, ]
  data.frame(UMAP_1 = crv[, 1], UMAP_2 = crv[, 2], Lineage = paste0("Lineage ", i))
})
curves_df <- do.call(rbind, curve_dfs)

p_pt_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.3, alpha = 0.6) +
  geom_path(data = curves_df, aes(x = UMAP_1, y = UMAP_2, group = Lineage),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_viridis_c(option = "inferno", name = "Pseudotime", na.value = "grey80") +
  labs(title = "Slingshot pseudotime — αβ CD8 trajectory") +
  theme_cowplot(font_size = 12) +
  NoAxes() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(traj_dir, "Trajectory_UMAP_pseudotime.png"),
       plot = p_pt_umap, width = 10, height = 8, dpi = 300, bg = "white")

# ── PLOT 2: UMAP colored by ART status with trajectory curves ────────────────
p_art_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Timepoint)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_path(data = curves_df, aes(x = UMAP_1, y = UMAP_2, group = Lineage),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_manual(values = art_colors, name = "ART Status") +
  labs(title = "ART status on CD8 trajectory") +
  theme_cowplot(font_size = 12) +
  NoAxes() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(traj_dir, "Trajectory_UMAP_ART_status.png"),
       plot = p_art_umap, width = 10, height = 8, dpi = 300, bg = "white")

# ── PLOT 3: Pseudotime density by ART status ─────────────────────────────────
pt_density_df <- data.frame(
  pseudotime = traj_cells$pseudotime,
  Timepoint  = traj_cells$Timepoint_Group
) %>% filter(!is.na(pseudotime))

p_pt_density <- ggplot(pt_density_df, aes(x = pseudotime, fill = Timepoint)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  scale_fill_manual(values = art_colors, name = "ART Status") +
  labs(x = "Pseudotime", y = "Density",
       title = "Pseudotime distribution by ART status") +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position  = "right",
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(traj_dir, "Trajectory_pseudotime_density_byART.png"),
       plot = p_pt_density, width = 10, height = 5, dpi = 300, bg = "white")

# ── PLOT 4: Pseudotime by cluster (box plot) ─────────────────────────────────
pt_cluster_df <- data.frame(
  pseudotime = traj_cells$pseudotime,
  Cluster    = traj_cells$CD8_Annotation,
  Timepoint  = traj_cells$Timepoint_Group
) %>% filter(!is.na(pseudotime))

pt_cluster_df$Cluster <- factor(pt_cluster_df$Cluster, levels = ab_clusters)

p_pt_cluster <- ggplot(pt_cluster_df, aes(x = Cluster, y = pseudotime, fill = Cluster)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, alpha = 0.7) +
  labs(x = NULL, y = "Pseudotime",
       title = "Pseudotime distribution by CD8 sub-cluster") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x      = element_text(size = 9, angle = 35, hjust = 1),
    legend.position  = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(traj_dir, "Trajectory_pseudotime_byCluster.png"),
       plot = p_pt_cluster, width = 10, height = 6, dpi = 300, bg = "white")

# ── PLOT 5: Pseudotime by ART status within effector clusters only ───────────
pt_eff_df <- pt_cluster_df %>%
  filter(Cluster %in% effector_clusters)

p_pt_art_eff <- ggplot(pt_eff_df, aes(x = Timepoint, y = pseudotime, fill = Timepoint)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, alpha = 0.7) +
  stat_compare_means(
    comparisons = list(
      c("PreART_Entry", "PostART_Suppressed"),
      c("PostART_Suppressed", "PostART_Unsuppressed"),
      c("PreART_Entry", "PostART_Unsuppressed")
    ),
    method     = "wilcox.test",
    label      = "p.signif",
    size       = 3.5,
    step.increase = 0.12,
    tip.length    = 0.01,
    hide.ns       = TRUE
  ) +
  facet_wrap(~ Cluster, nrow = 1) +
  scale_fill_manual(values = art_colors, name = "ART Status") +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "Pseudotime",
       title = "Pseudotime by ART status — effector clusters") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x      = element_text(size = 9, angle = 30, hjust = 1),
    strip.text        = element_text(size = 11),
    strip.background  = element_rect(fill = "grey95", color = NA),
    legend.position   = "none",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(traj_dir, "Trajectory_pseudotime_ART_effectorClusters.png"),
       plot = p_pt_art_eff, width = 12, height = 6, dpi = 300, bg = "white")

# ── Export pseudotime metadata ───────────────────────────────────────────────
write.csv(
  data.frame(
    cell          = colnames(traj_cells),
    CD8_Annotation = traj_cells$CD8_Annotation,
    Timepoint     = traj_cells$Timepoint_Group,
    pseudotime    = traj_cells$pseudotime,
    pseudotime_mat
  ),
  file.path(traj_dir, "Slingshot_pseudotime_per_cell.csv"),
  row.names = FALSE
)

cat("  Slingshot trajectory saved to:", traj_dir, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 8E (cont). MONOCLE3 TRAJECTORY ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Running Monocle3 trajectory analysis ===\n")

if (!requireNamespace("monocle3", quietly = TRUE)) BiocManager::install("monocle3")
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) remotes::install_github("satijalab/seurat-wrappers")
library(monocle3)
library(SeuratWrappers)

monocle_dir <- file.path(dir_07_trajectory, "monocle3")

DefaultAssay(traj_cells) <- "RNA"
cds <- as.cell_data_set(traj_cells)

reducedDims(cds)[["UMAP"]] <- Embeddings(traj_cells, "wnn.umap")

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# ── Order cells: root at naïve cluster ───────────────────────────────────────
naive_cells_ids <- colnames(traj_cells)[traj_cells$CD8_Annotation == "Naïve CD8"]

get_earliest_principal_node <- function(cds, cells) {
  cell_ids <- which(colnames(cds) %in% cells)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_ids <- closest_vertex[cell_ids, 1]
  vertex_counts <- table(closest_vertex_ids)
  most_common_vertex <- names(which.max(vertex_counts))
  graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  root_node <- graph_nodes[as.numeric(most_common_vertex)]
  root_node
}

root_node <- get_earliest_principal_node(cds, naive_cells_ids)
cat("  Root node:", root_node, "\n")
cds <- order_cells(cds, root_pr_nodes = root_node)

cat("  Monocle3 pseudotime computed. Root node:", root_node, "\n")

m3_pt <- pseudotime(cds)
m3_pt_vec <- m3_pt[colnames(traj_cells)]
names(m3_pt_vec) <- colnames(traj_cells)
traj_cells <- AddMetaData(traj_cells, metadata = m3_pt_vec, col.name = "monocle3_pseudotime")

# ── Monocle3 plots ───────────────────────────────────────────────────────────
p_m3_traj <- plot_cells(
  cds,
  color_cells_by     = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves       = FALSE,
  label_branch_points = TRUE,
  graph_label_size   = 4,
  cell_size          = 0.5
) +
  scale_color_viridis_c(option = "inferno", name = "Pseudotime") +
  ggtitle("Monocle3 trajectory — αβ CD8 sub-clusters") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(monocle_dir, "Monocle3_UMAP_pseudotime.png"),
       plot = p_m3_traj, width = 10, height = 8, dpi = 300, bg = "white")

p_m3_clust <- plot_cells(
  cds,
  color_cells_by     = "CD8_Annotation",
  label_groups_by_cluster = TRUE,
  label_leaves       = FALSE,
  label_branch_points = FALSE,
  group_label_size   = 3.5,
  cell_size          = 0.5
) +
  ggtitle("Monocle3 graph — CD8 sub-cluster annotations") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(monocle_dir, "Monocle3_UMAP_clusters.png"),
       plot = p_m3_clust, width = 12, height = 8, dpi = 300, bg = "white")

colData(cds)$Timepoint_Group <- traj_cells$Timepoint_Group

p_m3_art <- plot_cells(
  cds,
  color_cells_by     = "Timepoint_Group",
  label_groups_by_cluster = FALSE,
  label_leaves       = FALSE,
  label_branch_points = FALSE,
  cell_size          = 0.5
) +
  scale_color_manual(values = art_colors, name = "ART Status") +
  ggtitle("Monocle3 trajectory — ART status") +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(monocle_dir, "Monocle3_UMAP_ART_status.png"),
       plot = p_m3_art, width = 10, height = 8, dpi = 300, bg = "white")

m3_pt_df <- data.frame(
  pseudotime = traj_cells$monocle3_pseudotime,
  Timepoint  = traj_cells$Timepoint_Group
) %>% filter(!is.na(pseudotime) & is.finite(pseudotime))

p_m3_density <- ggplot(m3_pt_df, aes(x = pseudotime, fill = Timepoint)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  scale_fill_manual(values = art_colors, name = "ART Status") +
  labs(x = "Monocle3 Pseudotime", y = "Density",
       title = "Monocle3 pseudotime distribution by ART status") +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position  = "right",
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(monocle_dir, "Monocle3_pseudotime_density_byART.png"),
       plot = p_m3_density, width = 10, height = 5, dpi = 300, bg = "white")

# UPDATED effector cluster names
m3_eff_df <- data.frame(
  pseudotime = traj_cells$monocle3_pseudotime,
  Cluster    = traj_cells$CD8_Annotation,
  Timepoint  = traj_cells$Timepoint_Group
) %>%
  filter(!is.na(pseudotime) & is.finite(pseudotime)) %>%
  filter(Cluster %in% effector_clusters)

p_m3_eff <- ggplot(m3_eff_df, aes(x = Timepoint, y = pseudotime, fill = Timepoint)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, alpha = 0.7) +
  stat_compare_means(
    comparisons = list(
      c("PreART_Entry", "PostART_Suppressed"),
      c("PostART_Suppressed", "PostART_Unsuppressed"),
      c("PreART_Entry", "PostART_Unsuppressed")
    ),
    method     = "wilcox.test",
    label      = "p.signif",
    size       = 3.5,
    step.increase = 0.12,
    tip.length    = 0.01,
    hide.ns       = TRUE
  ) +
  facet_wrap(~ Cluster, nrow = 1) +
  scale_fill_manual(values = art_colors, name = "ART Status") +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "Monocle3 Pseudotime",
       title = "Monocle3 pseudotime by ART status — effector clusters") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x      = element_text(size = 9, angle = 30, hjust = 1),
    strip.text        = element_text(size = 11),
    strip.background  = element_rect(fill = "grey95", color = NA),
    legend.position   = "none",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(monocle_dir, "Monocle3_pseudotime_ART_effectorClusters.png"),
       plot = p_m3_eff, width = 12, height = 6, dpi = 300, bg = "white")

# ── graph_test ───────────────────────────────────────────────────────────────
cat("  Running Monocle3 graph_test (Moran's I)...\n")
tryCatch({
  graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  graph_test_sig <- graph_test_res %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)
  
  write.csv(graph_test_sig, file.path(monocle_dir, "Monocle3_graph_test_significant.csv"))
  cat("    Significant trajectory-varying genes:", nrow(graph_test_sig), "\n")
  cat("    Top 20:\n")
  print_cols <- intersect(c("gene_short_name", "morans_I", "q_value"), colnames(graph_test_sig))
  print(head(graph_test_sig[, print_cols, drop = FALSE], 20))
}, error = function(e) cat("    graph_test error:", e$message, "\n"))

# ── Genes along pseudotime panels ────────────────────────────────────────────
cat("  Generating gene expression along pseudotime panels...\n")

exhaustion_plot_genes <- c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1")
stemness_plot_genes   <- c("TCF7", "SELL", "CCR7", "LEF1", "BACH2", "BCL2", "IL7R", "S1PR1", "KLF2")
cytotoxic_plot_genes  <- c("GZMB", "GNLY", "PRF1", "NKG7", "GZMA", "FGFBP2")
terminal_plot_genes   <- c("ZEB2", "PRDM1", "TBX21", "CX3CR1", "ID2")
stress_plot_genes     <- c("HSPA1A", "HSPA1B", "IFIT1", "ISG15", "CCL3", "CCL4L2")

all_genes_avail <- rownames(traj_cells[["RNA"]])
filter_present <- function(genes) genes[genes %in% all_genes_avail]

exhaustion_present <- filter_present(exhaustion_plot_genes)
stemness_present   <- filter_present(stemness_plot_genes)
cytotoxic_present  <- filter_present(cytotoxic_plot_genes)
terminal_present   <- filter_present(terminal_plot_genes)
stress_present     <- filter_present(stress_plot_genes)

# ── Helper: gene × pseudotime plot ──────────────────────────────────────────
plot_genes_along_pt <- function(obj, genes, ncol_plot = 4, title_text = "",
                                filename, width, height) {
  DefaultAssay(obj) <- "RNA"
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  
  df <- data.frame(
    pseudotime = obj$monocle3_pseudotime,
    Timepoint  = obj$Timepoint_Group,
    t(as.matrix(expr_mat)),
    check.names = FALSE
  ) %>%
    filter(!is.na(pseudotime) & is.finite(pseudotime)) %>%
    tidyr::pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")
  
  df$Gene <- factor(df$Gene, levels = genes)
  df$Timepoint <- factor(df$Timepoint,
                         levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
  
  p <- ggplot(df, aes(x = pseudotime, y = Expression, color = Timepoint)) +
    geom_point(size = 0.1, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2, span = 0.5) +
    facet_wrap(~ Gene, ncol = ncol_plot, scales = "free_y") +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    labs(x = "Monocle3 Pseudotime", y = "Expression", title = title_text) +
    theme_cowplot(font_size = 11) +
    theme(
      strip.text        = element_text(face = "italic", size = 10),
      strip.background  = element_rect(fill = "grey95", color = NA),
      legend.position   = "bottom",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300, bg = "white")
  cat("    Saved:", basename(filename), "(", length(genes), "genes )\n")
}

# Top variable genes from graph_test
if (exists("graph_test_sig") && nrow(graph_test_sig) > 0) {
  top_variable <- graph_test_sig %>% arrange(q_value) %>% head(16)
  top_variable <- if ("gene_short_name" %in% colnames(top_variable)) top_variable$gene_short_name else rownames(top_variable)
  top_variable <- filter_present(top_variable)
  
  if (length(top_variable) >= 4) {
    plot_genes_along_pt(traj_cells, top_variable, ncol_plot = 4,
                        title_text = "Top trajectory-varying genes (graph_test)",
                        filename = file.path(monocle_dir, "Monocle3_TopVariableGenes_in_pseudotime.png"),
                        width = 16, height = ceiling(length(top_variable) / 4) * 3.5)
  }
}

if (length(exhaustion_present) >= 2) {
  plot_genes_along_pt(traj_cells, exhaustion_present, ncol_plot = 4,
                      title_text = "Exhaustion genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_ExhaustionGenes_in_pseudotime.png"),
                      width = 16, height = ceiling(length(exhaustion_present) / 4) * 3.5)
}

if (length(stemness_present) >= 2) {
  plot_genes_along_pt(traj_cells, stemness_present, ncol_plot = 4,
                      title_text = "Stemness genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_StemnessGenes_in_pseudotime.png"),
                      width = 16, height = ceiling(length(stemness_present) / 4) * 3.5)
}

if (length(cytotoxic_present) >= 2) {
  plot_genes_along_pt(traj_cells, cytotoxic_present, ncol_plot = 3,
                      title_text = "Cytotoxicity genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_CytotoxicityGenes_in_pseudotime.png"),
                      width = 12, height = ceiling(length(cytotoxic_present) / 3) * 3.5)
}

if (length(terminal_present) >= 2) {
  plot_genes_along_pt(traj_cells, terminal_present, ncol_plot = 3,
                      title_text = "Terminal differentiation genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_TerminalDiffGenes_in_pseudotime.png"),
                      width = 12, height = ceiling(length(terminal_present) / 3) * 3.5)
}

if (length(stress_present) >= 2) {
  plot_genes_along_pt(traj_cells, stress_present, ncol_plot = 3,
                      title_text = "Stress / IFN / chemokine genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_StressChemokineGenes_in_pseudotime.png"),
                      width = 12, height = ceiling(length(stress_present) / 3) * 3.5)
}

narrative_genes <- c(
  "TOX", "PDCD1", "HAVCR2", "TIGIT",
  "TCF7", "BACH2", "SELL", "IL7R",
  "GNLY", "GZMB", "PRF1",
  "HSPA1A", "CCL3", "IFIT1"
)
narrative_present <- filter_present(narrative_genes)

if (length(narrative_present) >= 4) {
  plot_genes_along_pt(traj_cells, narrative_present, ncol_plot = 4,
                      title_text = "Key narrative genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_NarrativeGenes_in_pseudotime.png"),
                      width = 16, height = ceiling(length(narrative_present) / 4) * 3.5)
}

# ── UMAP colored by gene expression (top graph_test genes) ───────────────
if (exists("graph_test_sig") && nrow(graph_test_sig) > 0) {
  top_umap_df <- graph_test_sig %>% arrange(q_value) %>% head(12)
  top_umap_genes <- if ("gene_short_name" %in% colnames(top_umap_df)) top_umap_df$gene_short_name else rownames(top_umap_df)
  top_umap_genes <- filter_present(top_umap_genes)
  
  if (length(top_umap_genes) >= 4) {
    if (!"gene_short_name" %in% colnames(rowData(cds))) {
      rowData(cds)$gene_short_name <- rownames(rowData(cds))
    }
    
    p_gene_umap <- plot_cells(
      cds,
      genes               = top_umap_genes,
      label_cell_groups   = FALSE,
      label_leaves        = FALSE,
      label_branch_points = FALSE
    ) + ggtitle("Top pseudotime-varying genes - UMAP expression")
    
    ggsave(file.path(monocle_dir, "Monocle3_TopGenes_UMAP_expression.png"),
           plot = p_gene_umap, width = 14, height = 12, dpi = 300, bg = "white")
    cat("    Top genes UMAP expression plot saved.\n")
  }
}

cat("  Gene expression along pseudotime panels done.\n")

# ── Module scores along Monocle3 pseudotime (deferred — requires 8H) ────────
if (exists("expand_eff")) {
  cat("  Plotting module scores along Monocle3 pseudotime...\n")
  
  common_cells_m3 <- intersect(colnames(traj_cells), colnames(expand_eff))
  available_modules <- score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)]
  
  if (length(common_cells_m3) > 100 & length(available_modules) > 0) {
    
    m3_module_df <- data.frame(
      pseudotime = traj_cells$monocle3_pseudotime[match(common_cells_m3, colnames(traj_cells))],
      Timepoint  = expand_eff@meta.data[common_cells_m3, "Timepoint_Group"],
      Cluster    = expand_eff@meta.data[common_cells_m3, "CD8_Annotation"]
    )
    for (mod in available_modules) {
      m3_module_df[[mod]] <- expand_eff@meta.data[common_cells_m3, mod]
    }
    
    m3_module_df <- m3_module_df %>%
      filter(!is.na(pseudotime) & is.finite(pseudotime))
    
    cat("    Cells with pseudotime + module scores:", nrow(m3_module_df), "\n")
    
    m3_exh_stem <- m3_module_df %>%
      select(pseudotime, Timepoint, Exhaustion, Stemness) %>%
      tidyr::pivot_longer(cols = c(Exhaustion, Stemness),
                          names_to = "Module", values_to = "Score")
    
    p_m3_exh_stem <- ggplot(m3_exh_stem, aes(x = pseudotime, y = Score, color = Timepoint)) +
      geom_point(size = 0.2, alpha = 0.15) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.4, alpha = 0.2, span = 0.4) +
      facet_wrap(~ Module, nrow = 1, scales = "free_y") +
      scale_color_manual(values = art_colors, name = "ART Status") +
      labs(x = "Monocle3 Pseudotime", y = "Module Score",
           title = "Exhaustion and stemness along Monocle3 trajectory — expanding clones") +
      theme_cowplot(font_size = 12) +
      theme(
        strip.text        = element_text(size = 12, face = "bold"),
        strip.background  = element_rect(fill = "grey95", color = NA),
        legend.position   = "bottom",
        plot.background   = element_rect(fill = "white", color = NA),
        panel.background  = element_rect(fill = "white", color = NA)
      ) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(monocle_dir, "Monocle3_ExhaustionStemness_along_pseudotime.png"),
           plot = p_m3_exh_stem, width = 14, height = 6, dpi = 300, bg = "white")
    
    m3_all_mods <- m3_module_df %>%
      select(pseudotime, Timepoint, all_of(available_modules)) %>%
      tidyr::pivot_longer(cols = all_of(available_modules),
                          names_to = "Module", values_to = "Score")
    
    m3_all_mods$Module_Label <- module_labels[m3_all_mods$Module]
    m3_all_mods$Module_Label <- factor(m3_all_mods$Module_Label, levels = module_labels)
    
    p_m3_all_mods <- ggplot(m3_all_mods, aes(x = pseudotime, y = Score, color = Timepoint)) +
      geom_point(size = 0.1, alpha = 0.1) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2, span = 0.4) +
      facet_wrap(~ Module_Label, ncol = 2, scales = "free_y") +
      scale_color_manual(values = art_colors, name = "ART Status") +
      labs(x = "Monocle3 Pseudotime", y = "Module Score",
           title = "All module scores along Monocle3 trajectory — expanding clones") +
      theme_cowplot(font_size = 11) +
      theme(
        strip.text        = element_text(size = 10),
        strip.background  = element_rect(fill = "grey95", color = NA),
        legend.position   = "bottom",
        panel.spacing     = unit(0.5, "lines"),
        plot.background   = element_rect(fill = "white", color = NA),
        panel.background  = element_rect(fill = "white", color = NA)
      ) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(monocle_dir, "Monocle3_AllModules_along_pseudotime.png"),
           plot = p_m3_all_mods, width = 12, height = 16, dpi = 300, bg = "white")
    
    for (mod in available_modules) {
      m3_mod_sub <- m3_module_df %>%
        select(pseudotime, Timepoint, Score = !!sym(mod))
      
      p_m3_ind <- ggplot(m3_mod_sub, aes(x = pseudotime, y = Score, color = Timepoint)) +
        geom_point(size = 0.2, alpha = 0.15) +
        geom_smooth(method = "loess", se = TRUE, linewidth = 1.4, alpha = 0.2, span = 0.4) +
        scale_color_manual(values = art_colors, name = "ART Status") +
        labs(x = "Monocle3 Pseudotime", y = "Module Score",
             title = paste0(module_labels[mod], " along Monocle3 trajectory")) +
        theme_cowplot(font_size = 12) +
        theme(
          legend.position  = "bottom",
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
      
      safe_mod <- gsub("[^A-Za-z0-9]", "_", mod)
      ggsave(file.path(monocle_dir, paste0("Monocle3_", safe_mod, "_along_pseudotime.png")),
             plot = p_m3_ind, width = 10, height = 6, dpi = 300, bg = "white")
    }
    
    # Faceted by cluster: Exhaustion along pseudotime per effector cluster
    m3_exh_cluster <- m3_module_df %>%
      filter(Cluster %in% effector_clusters) %>%
      select(pseudotime, Timepoint, Cluster, Exhaustion)
    
    p_m3_exh_cl <- ggplot(m3_exh_cluster, aes(x = pseudotime, y = Exhaustion, color = Timepoint)) +
      geom_point(size = 0.2, alpha = 0.15) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.4, alpha = 0.2, span = 0.5) +
      facet_wrap(~ Cluster, nrow = 1) +
      scale_color_manual(values = art_colors, name = "ART Status") +
      labs(x = "Monocle3 Pseudotime", y = "Exhaustion Module Score",
           title = "Exhaustion along trajectory — per effector cluster") +
      theme_cowplot(font_size = 12) +
      theme(
        strip.text        = element_text(size = 11, face = "bold"),
        strip.background  = element_rect(fill = "grey95", color = NA),
        legend.position   = "bottom",
        plot.background   = element_rect(fill = "white", color = NA),
        panel.background  = element_rect(fill = "white", color = NA)
      ) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(monocle_dir, "Monocle3_Exhaustion_per_cluster_along_pseudotime.png"),
           plot = p_m3_exh_cl, width = 14, height = 6, dpi = 300, bg = "white")
    
    cat("  Module score × pseudotime plots saved.\n")
  } else {
    cat("  Too few overlapping cells for module score × pseudotime plots.\n")
  }
} else {
  cat("  Module score × pseudotime plots deferred — expand_eff not yet created (runs in 8H).\n")
}

# ── Export Monocle3 pseudotime ───────────────────────────────────────────────
write.csv(
  data.frame(
    cell               = colnames(traj_cells),
    CD8_Annotation     = traj_cells$CD8_Annotation,
    Timepoint          = traj_cells$Timepoint_Group,
    monocle3_pseudotime = traj_cells$monocle3_pseudotime
  ),
  file.path(monocle_dir, "Monocle3_pseudotime_per_cell.csv"),
  row.names = FALSE
)

cat("  Monocle3 analysis saved to:", monocle_dir, "\n")

# ── CHECKPOINT: Save after trajectory analysis ───────────────────────────────
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_with_trajectory.qs2"))
cat("CHECKPOINT: Object with trajectory saved as TARA_cd8_HEI_with_trajectory.qs2\n")

################################################################################
# 8F. CLONAL EXPANSION vs EXHAUSTION: controlling for expansion level
################################################################################
cat("\n=== Clonal expansion vs exhaustion control analysis ===\n")

expansion_dir <- dir_08_expansion

if (!"PID" %in% colnames(TARA_cd8@meta.data)) {
  pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
  names(pid_vec) <- colnames(TARA_cd8)
  TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")
}

# ── PLOT G1: Expanding cells per sample, by ART status × cluster ──────────────
total_per_sample <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters) %>%
  group_by(orig.ident, PID, Timepoint_Group, CD8_Annotation) %>%
  summarise(total_cells = n(), .groups = "drop")

expand_per_sample <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters &
           !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") %>%
  group_by(orig.ident, PID, Timepoint_Group, CD8_Annotation) %>%
  summarise(n_expanding = n(), .groups = "drop")

sample_expansion <- total_per_sample %>%
  left_join(expand_per_sample, by = c("orig.ident", "PID", "Timepoint_Group", "CD8_Annotation")) %>%
  mutate(
    n_expanding = replace_na(n_expanding, 0),
    pct_expanding = n_expanding / total_cells * 100
  )

sample_expansion$Timepoint_Group <- factor(sample_expansion$Timepoint_Group,
                                           levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
sample_expansion$CD8_Annotation <- factor(sample_expansion$CD8_Annotation, levels = effector_clusters)

p_g1a <- ggplot(sample_expansion, aes(x = Timepoint_Group, y = n_expanding, fill = Timepoint_Group)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("PreART_Entry", "PostART_Suppressed"),
      c("PostART_Suppressed", "PostART_Unsuppressed"),
      c("PreART_Entry", "PostART_Unsuppressed")
    ),
    method = "wilcox.test", label = "p.signif",
    size = 3.5, step.increase = 0.12, tip.length = 0.01, hide.ns = TRUE
  ) +
  facet_wrap(~ CD8_Annotation, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = art_colors, name = NULL) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "Expanding cells per sample",
       title = "Clonal expansion level per sample by ART status") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x       = element_text(size = 9, angle = 30, hjust = 1),
    strip.text        = element_text(size = 11, face = "bold"),
    strip.background  = element_rect(fill = "grey95", color = NA),
    legend.position   = "none",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(expansion_dir, "G1a_ExpandingCells_perSample_byART.png"),
       plot = p_g1a, width = 14, height = 6, dpi = 300, bg = "white")

p_g1b <- ggplot(sample_expansion, aes(x = Timepoint_Group, y = pct_expanding, fill = Timepoint_Group)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("PreART_Entry", "PostART_Suppressed"),
      c("PostART_Suppressed", "PostART_Unsuppressed"),
      c("PreART_Entry", "PostART_Unsuppressed")
    ),
    method = "wilcox.test", label = "p.signif",
    size = 3.5, step.increase = 0.12, tip.length = 0.01, hide.ns = TRUE
  ) +
  facet_wrap(~ CD8_Annotation, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = art_colors, name = NULL) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "% cells that are clonally expanded",
       title = "Proportion of clonally expanded cells per sample by ART status") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x       = element_text(size = 9, angle = 30, hjust = 1),
    strip.text        = element_text(size = 11, face = "bold"),
    strip.background  = element_rect(fill = "grey95", color = NA),
    legend.position   = "none",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(expansion_dir, "G1b_PctExpanding_perSample_byART.png"),
       plot = p_g1b, width = 14, height = 6, dpi = 300, bg = "white")

write.csv(sample_expansion, file.path(expansion_dir, "PerSample_expansion_counts.csv"),
          row.names = FALSE)

# ── PLOT G2: Clone size distribution per ART status ──────────────────────────
clone_art_df <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters &
           !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") %>%
  mutate(
    Timepoint_Group = factor(Timepoint_Group,
                             levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed")),
    cloneSize = factor(cloneSize, levels = c(
      "Small (1 < X <= 5)", "Medium (5 < X <= 20)",
      "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 500)"))
  )

clone_art_summary <- clone_art_df %>%
  group_by(CD8_Annotation, Timepoint_Group, cloneSize) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CD8_Annotation, Timepoint_Group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

clone_colors_short <- setNames(
  hcl.colors(n = 7, palette = "plasma", fixup = TRUE)[c(3, 4, 5, 7)],
  levels(clone_art_summary$cloneSize)
)

p_g2 <- ggplot(clone_art_summary,
               aes(x = Timepoint_Group, y = prop, fill = cloneSize)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.65, color = "white", linewidth = 0.3) +
  facet_wrap(~ CD8_Annotation, nrow = 1) +
  scale_fill_manual(values = clone_colors_short, name = "Clone Size") +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of expanding cells",
       title = "Clone size distribution by ART status — effector clusters") +
  theme_cowplot(font_size = 12) +
  theme(
    axis.text.x       = element_text(size = 9, angle = 30, hjust = 1),
    strip.text        = element_text(size = 11, face = "bold"),
    strip.background  = element_rect(fill = "grey95", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 8),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(expansion_dir, "G2_CloneSizeDistribution_byART_byCluster.png"),
       plot = p_g2, width = 14, height = 6, dpi = 300, bg = "white")

cat("  G3/G3b (exhaustion/stemness by clone size) deferred to section 8H.\n")
cat("  Expansion plots saved to:", expansion_dir, "\n")

################################################################################
# 8G. VIRAL LOAD CORRELATION — deferred to 8H
################################################################################
cat("\n=== Viral load correlation analysis ===\n")
cat("  NOTE: Viral load correlations require module scores.\n")
cat("  Deferred to section 8H where module scores are computed.\n")

################################################################################
# 8H. MODULE SCORES: Functional signatures in expanding clones by ART status
################################################################################
cat("\n=== Computing module scores for expanding clones ===\n")

module_dir <- dir_10_modules

module_gene_lists <- list(
  Exhaustion = c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1"),
  Stemness   = c("TCF7", "SELL", "CCR7", "LEF1", "BACH2", "BCL2", "IL7R", "S1PR1", "KLF2"),
  Cytotoxicity = c("GZMB", "GNLY", "PRF1", "NKG7", "GZMA", "GZMH", "FGFBP2"),
  Stress_IFN = c("HSPA1A", "HSPA1B", "DNAJB1", "IFI27", "IFI44L"),
  TypeI_IFN_Memory = c("IFIT1", "IFIT3", "ISG15", "MX1"),
  Inflammatory_Chemokines = c("CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5"),
  Terminal_Differentiation = c("ZEB2", "PRDM1", "TBX21", "CX3CR1", "S1PR5", "ID2")
)

# UPDATED effector cluster names
expand_eff <- subset(TARA_cd8,
                     CD8_Annotation %in% effector_clusters &
                       !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)"
)

cat("  Expanding effector cells:", ncol(expand_eff), "\n")
cat("  Timepoint breakdown:\n")
print(table(expand_eff$CD8_Annotation, expand_eff$Timepoint_Group))

DefaultAssay(expand_eff) <- "RNA"

for (mod_name in names(module_gene_lists)) {
  genes <- module_gene_lists[[mod_name]]
  genes_present <- genes[genes %in% rownames(expand_eff[["RNA"]])]
  
  if (length(genes_present) < 2) {
    cat("  Skipping", mod_name, "— too few genes found\n")
    next
  }
  
  cat("  Computing", mod_name, "(", length(genes_present), "/", length(genes), "genes )...\n")
  expand_eff <- AddModuleScore(
    expand_eff,
    features = list(genes_present),
    name     = mod_name,
    ctrl     = min(100, length(genes_present) * 5)
  )
}

score_cols <- paste0(names(module_gene_lists), "1")
score_cols_clean <- names(module_gene_lists)
for (i in seq_along(score_cols)) {
  if (score_cols[i] %in% colnames(expand_eff@meta.data)) {
    colnames(expand_eff@meta.data)[colnames(expand_eff@meta.data) == score_cols[i]] <- score_cols_clean[i]
  }
}

# ── Build long-format data for plotting ──────────────────────────────────────
plot_data <- expand_eff@meta.data %>%
  select(CD8_Annotation, Timepoint_Group, all_of(score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)])) %>%
  tidyr::pivot_longer(
    cols      = all_of(score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)]),
    names_to  = "Module",
    values_to = "Score"
  )

plot_data$Timepoint_Group <- factor(plot_data$Timepoint_Group,
                                    levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
plot_data$CD8_Annotation <- factor(plot_data$CD8_Annotation, levels = effector_clusters)

module_labels <- c(
  "Exhaustion"               = "Exhaustion",
  "Stemness"                 = "Stemness / naïve",
  "Cytotoxicity"             = "Cytotoxicity",
  "Stress_IFN"               = "Stress / IFN (acute)",
  "TypeI_IFN_Memory"         = "Type I IFN memory",
  "Inflammatory_Chemokines"  = "Inflammatory chemokines",
  "Terminal_Differentiation"  = "Terminal differentiation"
)
plot_data$Module_Label <- module_labels[plot_data$Module]
plot_data$Module_Label <- factor(plot_data$Module_Label, levels = module_labels)

art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel", repos = "https://cloud.r-project.org")
library(ggpubr)
library(ggridges)
library(ggrepel)

comparisons_list <- list(
  c("PreART_Entry", "PostART_Suppressed"),
  c("PostART_Suppressed", "PostART_Unsuppressed"),
  c("PreART_Entry", "PostART_Unsuppressed")
)

# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 1: COHEN'S D BAR PLOTS
# ══════════════════════════════════════════════════════════════════════════════
cat("  Drawing Cohen's d bar plots...\n")

stat_results_fc <- list()
for (mod in score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)]) {
  for (cl in effector_clusters) {
    cl_data <- expand_eff@meta.data %>% filter(CD8_Annotation == cl)
    sup_vals   <- cl_data %>% filter(Timepoint_Group == "PostART_Suppressed") %>% pull(!!sym(mod))
    unsup_vals <- cl_data %>% filter(Timepoint_Group == "PostART_Unsuppressed") %>% pull(!!sym(mod))
    pre_vals   <- cl_data %>% filter(Timepoint_Group == "PreART_Entry") %>% pull(!!sym(mod))
    
    if (length(sup_vals) >= 5 & length(unsup_vals) >= 5) {
      wt <- wilcox.test(sup_vals, unsup_vals)
      stat_results_fc[[paste(mod, cl, "Sup_vs_Unsup")]] <- data.frame(
        Module = mod, CD8_Annotation = cl, Comparison = "Suppressed vs Unsuppressed",
        p_value = wt$p.value, stringsAsFactors = FALSE)
    }
    if (length(pre_vals) >= 5 & length(unsup_vals) >= 5) {
      wt2 <- wilcox.test(pre_vals, unsup_vals)
      stat_results_fc[[paste(mod, cl, "Pre_vs_Unsup")]] <- data.frame(
        Module = mod, CD8_Annotation = cl, Comparison = "Pre-ART vs Unsuppressed",
        p_value = wt2$p.value, stringsAsFactors = FALSE)
    }
    if (length(sup_vals) >= 5 & length(pre_vals) >= 5) {
      wt3 <- wilcox.test(sup_vals, pre_vals)
      stat_results_fc[[paste(mod, cl, "Sup_vs_Pre")]] <- data.frame(
        Module = mod, CD8_Annotation = cl, Comparison = "Suppressed vs Pre-ART",
        p_value = wt3$p.value, stringsAsFactors = FALSE)
    }
  }
}

stat_fc <- do.call(rbind, stat_results_fc)

stat_fc <- stat_fc %>%
  group_by(Comparison) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

stat_fc$star <- ifelse(stat_fc$p_adj < 0.001, "***",
                       ifelse(stat_fc$p_adj < 0.01, "**",
                              ifelse(stat_fc$p_adj < 0.05, "*", "")))

write.csv(stat_fc, file.path(module_dir, "ModuleScore_Wilcoxon_all_comparisons.csv"), row.names = FALSE)

cat("\n=== All module score Wilcoxon results ===\n")
print(stat_fc %>% arrange(Comparison, Module, CD8_Annotation) %>%
        mutate(p_value = formatC(p_value, format = "e", digits = 2),
               p_adj   = formatC(p_adj, format = "e", digits = 2)),
      n = 100)

# ── Helper: build FC data for a given comparison ─────────────────────────────
build_fc <- function(plot_data_raw, num_tp, denom_tp, comp_label, stat_fc_df,
                     module_labels, effector_clusters) {
  fc <- plot_data_raw %>%
    filter(Timepoint_Group %in% c(num_tp, denom_tp)) %>%
    group_by(Module, Module_Label, CD8_Annotation) %>%
    summarise(
      mean_num   = mean(Score[Timepoint_Group == num_tp], na.rm = TRUE),
      mean_denom = mean(Score[Timepoint_Group == denom_tp], na.rm = TRUE),
      sd_num     = sd(Score[Timepoint_Group == num_tp], na.rm = TRUE),
      sd_denom   = sd(Score[Timepoint_Group == denom_tp], na.rm = TRUE),
      n_num      = sum(Timepoint_Group == num_tp),
      n_denom    = sum(Timepoint_Group == denom_tp),
      .groups    = "drop"
    ) %>%
    mutate(
      pooled_sd = sqrt(((n_num - 1) * sd_num^2 + (n_denom - 1) * sd_denom^2) / (n_num + n_denom - 2)),
      cohens_d  = (mean_num - mean_denom) / pooled_sd,
      direction = ifelse(cohens_d > 0,
                         paste0("\u2191 ", gsub("PostART_", "", num_tp)),
                         paste0("\u2191 ", gsub("PostART_", "", denom_tp)))
    )
  
  fc$Comparison <- comp_label
  
  stars <- stat_fc_df %>%
    filter(Comparison == comp_label) %>%
    select(Module, CD8_Annotation, star)
  fc <- fc %>% left_join(stars, by = c("Module", "CD8_Annotation"))
  fc$star[is.na(fc$star)] <- ""
  
  fc$Module_Label <- factor(fc$Module_Label, levels = rev(module_labels))
  fc$CD8_Annotation <- factor(fc$CD8_Annotation, levels = effector_clusters)
  fc
}

fc_sup_unsup <- build_fc(plot_data, "PostART_Suppressed", "PostART_Unsuppressed",
                         "Suppressed vs Unsuppressed", stat_fc, module_labels, effector_clusters)
fc_pre_unsup <- build_fc(plot_data, "PreART_Entry", "PostART_Unsuppressed",
                         "Pre-ART vs Unsuppressed", stat_fc, module_labels, effector_clusters)
fc_sup_pre   <- build_fc(plot_data, "PostART_Suppressed", "PreART_Entry",
                         "Suppressed vs Pre-ART", stat_fc, module_labels, effector_clusters)

# ── Helper: plot one FC bar chart ────────────────────────────────────────────
plot_fc_bars <- function(fc_df, title_text, fill_colors) {
  x_range <- max(abs(fc_df$cohens_d), na.rm = TRUE)
  offset <- x_range * 0.06
  
  fc_df$star_x <- ifelse(fc_df$cohens_d > 0,
                         fc_df$cohens_d + offset,
                         fc_df$cohens_d - offset)
  fc_df$star_hjust <- ifelse(fc_df$cohens_d > 0, 0, 1)
  
  ggplot(fc_df, aes(x = cohens_d, y = Module_Label, fill = direction)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(aes(x = star_x, label = star, hjust = star_hjust),
              size = 4.5, color = "black", vjust = 0.4) +
    facet_wrap(~ CD8_Annotation, nrow = 1) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    coord_cartesian(clip = "off") +
    labs(
      x     = "Cohen's d (effect size)",
      y     = NULL,
      title = title_text
    ) +
    theme_cowplot(font_size = 12) +
    theme(
      strip.text        = element_text(size = 11, face = "bold"),
      strip.background  = element_rect(fill = "grey95", color = NA),
      axis.text.y       = element_text(size = 10),
      legend.position   = "bottom",
      legend.text       = element_text(size = 10),
      plot.margin       = margin(10, 30, 10, 10),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
}

p_fc1 <- plot_fc_bars(fc_sup_unsup,
                      "Suppressed vs Unsuppressed — expanding clones",
                      c("\u2191 Suppressed" = "#52B788", "\u2191 Unsuppressed" = "#E76F51"))

ggsave(file.path(module_dir, "Fig4E_DiffMean_Sup_vs_Unsup.png"),
       plot = p_fc1, width = 14, height = 7, dpi = 300, bg = "white")

p_fc2 <- plot_fc_bars(fc_pre_unsup,
                      "Pre-ART vs Unsuppressed — expanding clones",
                      c("\u2191 PreART_Entry" = "#4A90D9", "\u2191 Unsuppressed" = "#E76F51"))

ggsave(file.path(module_dir, "Fig4E_DiffMean_PreART_vs_Unsup.png"),
       plot = p_fc2, width = 14, height = 7, dpi = 300, bg = "white")

p_fc3 <- plot_fc_bars(fc_sup_pre,
                      "Suppressed vs Pre-ART — expanding clones",
                      c("\u2191 Suppressed" = "#52B788", "\u2191 PreART_Entry" = "#4A90D9"))

ggsave(file.path(module_dir, "Fig4E_DiffMean_Sup_vs_PreART.png"),
       plot = p_fc3, width = 14, height = 7, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 1B: EXHAUSTION + STEMNESS — all 3 comparisons
# ══════════════════════════════════════════════════════════════════════════════
cat("  Drawing exhaustion + stemness focus plots...\n")

comp_colors <- c(
  "Suppressed vs Unsuppressed" = "#52B788",
  "Pre-ART vs Unsuppressed"    = "#4A90D9",
  "Suppressed vs Pre-ART"      = "#9C6FD6"
)

for (focus_mod in c("Exhaustion", "Stemness")) {
  focus_fc <- bind_rows(
    fc_sup_unsup %>% filter(Module == focus_mod),
    fc_pre_unsup %>% filter(Module == focus_mod),
    fc_sup_pre   %>% filter(Module == focus_mod)
  )
  focus_fc$Comparison <- factor(focus_fc$Comparison,
                                levels = c("Suppressed vs Unsuppressed", "Pre-ART vs Unsuppressed", "Suppressed vs Pre-ART"))
  focus_fc$star_x <- ifelse(focus_fc$cohens_d > 0,
                            focus_fc$cohens_d + max(abs(focus_fc$cohens_d)) * 0.08,
                            focus_fc$cohens_d - max(abs(focus_fc$cohens_d)) * 0.08)
  focus_fc$star_hjust <- ifelse(focus_fc$cohens_d > 0, 0, 1)
  
  p_focus <- ggplot(focus_fc, aes(x = cohens_d, y = Comparison, fill = Comparison)) +
    geom_col(width = 0.55, alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(aes(x = star_x, label = star, hjust = star_hjust),
              size = 5, color = "black", vjust = 0.4) +
    facet_wrap(~ CD8_Annotation, nrow = 1) +
    scale_fill_manual(values = comp_colors, name = NULL) +
    scale_y_discrete(limits = rev(levels(focus_fc$Comparison))) +
    coord_cartesian(clip = "off") +
    labs(
      x     = paste0("Cohen's d (", tolower(focus_mod), " module)"),
      y     = NULL,
      title = paste0(focus_mod, " module — all pairwise comparisons")
    ) +
    theme_cowplot(font_size = 12) +
    theme(
      strip.text        = element_text(size = 11, face = "bold"),
      strip.background  = element_rect(fill = "grey95", color = NA),
      axis.text.y       = element_text(size = 10),
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(module_dir, paste0("Fig4E_", focus_mod, "_AllComparisons.png")),
         plot = p_focus, width = 14, height = 5, dpi = 300, bg = "white")
}

cat("  All Cohen's d bar plots saved.\n")

# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION 2: VOLCANO PLOTS
# ══════════════════════════════════════════════════════════════════════════════
cat("  Drawing volcano plots...\n")

volcano_dir <- file.path(dir_10_modules, "volcanos")

exhaustion_genes   <- c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1")
stemness_genes     <- c("TCF7", "SELL", "CCR7", "LEF1", "BACH2", "BCL2", "IL7R", "S1PR1", "KLF2")
cytotoxic_genes    <- c("GZMB", "GNLY", "PRF1", "NKG7", "GZMA", "GZMH", "FGFBP2")
stress_genes       <- c("HSPA1A", "HSPA1B", "DNAJB1", "IFI27", "IFI44L")
ifn_memory_genes   <- c("IFIT1", "IFIT3", "ISG15", "MX1")
chemokine_genes    <- c("CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5")
terminal_genes     <- c("ZEB2", "PRDM1", "TBX21", "CX3CR1", "S1PR5", "ID2", "EOMES")
activation_genes   <- c("HLA-DRA", "CD38", "FAS", "CD69")

highlight_cols <- c(
  "Exhaustion"        = "#A32D2D",
  "Stemness"          = "#185FA5",
  "Cytotoxicity"      = "#3B6D11",
  "Stress response"   = "#BA7517",
  "Type I IFN memory" = "#0F6E56",
  "Chemokines"        = "#993556",
  "Terminal diff."    = "#534AB7",
  "Activation"        = "#D85A30",
  "Other"             = "grey80"
)

# UPDATED DGE file names
dge_files <- list(
  "TEM CD8"               = file.path(dge_dir, "DGE_MAST_Expanding_TEMCD8_Sup_vs_Unsup.csv"),
  "TEMRA CD8"             = file.path(dge_dir, "DGE_MAST_Expanding_TEMRACD8_Sup_vs_Unsup.csv"),
  "Transitional Tem CD8"  = file.path(dge_dir, "DGE_MAST_Expanding_TransitionalTem_Sup_vs_Unsup.csv")
)

for (cl_name in names(dge_files)) {
  fpath <- dge_files[[cl_name]]
  if (!file.exists(fpath)) {
    cat("    Skipping", cl_name, "— DGE file not found\n")
    next
  }
  
  dge_data <- read.csv(fpath, row.names = 1)
  
  dge_data$highlight <- "Other"
  dge_data$highlight[rownames(dge_data) %in% exhaustion_genes]  <- "Exhaustion"
  dge_data$highlight[rownames(dge_data) %in% stemness_genes]    <- "Stemness"
  dge_data$highlight[rownames(dge_data) %in% cytotoxic_genes]   <- "Cytotoxicity"
  dge_data$highlight[rownames(dge_data) %in% stress_genes]      <- "Stress response"
  dge_data$highlight[rownames(dge_data) %in% ifn_memory_genes]  <- "Type I IFN memory"
  dge_data$highlight[rownames(dge_data) %in% chemokine_genes]   <- "Chemokines"
  dge_data$highlight[rownames(dge_data) %in% terminal_genes]    <- "Terminal diff."
  dge_data$highlight[rownames(dge_data) %in% activation_genes]  <- "Activation"
  
  genes_to_label <- c(exhaustion_genes, stemness_genes, cytotoxic_genes,
                      stress_genes, ifn_memory_genes, chemokine_genes,
                      terminal_genes, activation_genes)
  genes_to_label <- genes_to_label[genes_to_label %in% rownames(dge_data)]
  
  dge_data$gene <- rownames(dge_data)
  dge_data$neg_log10_padj <- -log10(dge_data$p_val_adj + 1e-300)
  dge_data$sig_cat <- "ns"
  dge_data$sig_cat[dge_data$p_val_adj < 0.05 & abs(dge_data$avg_log2FC) > 0.25] <- "sig"
  
  dge_data$label <- ""
  sig_highlight <- dge_data$gene %in% genes_to_label & dge_data$p_val_adj < 0.05
  dge_data$label[sig_highlight] <- dge_data$gene[sig_highlight]
  nonsig_highlight <- dge_data$gene %in% genes_to_label & dge_data$p_val_adj >= 0.05
  dge_data$label[nonsig_highlight] <- dge_data$gene[nonsig_highlight]
  
  dge_data$highlight <- factor(dge_data$highlight,
                               levels = c("Exhaustion", "Stemness", "Cytotoxicity", "Stress response",
                                          "Type I IFN memory", "Chemokines", "Terminal diff.", "Activation", "Other"))
  
  p_volc <- ggplot(dge_data, aes(x = avg_log2FC, y = neg_log10_padj, color = highlight)) +
    geom_point(data = dge_data %>% filter(highlight == "Other"),
               size = 0.5, alpha = 0.3) +
    geom_point(data = dge_data %>% filter(highlight != "Other"),
               size = 2.5, alpha = 0.85) +
    ggrepel::geom_text_repel(
      data = dge_data %>% filter(label != "" & p_val_adj < 0.05),
      aes(label = label), size = 3.5, fontface = "italic",
      max.overlaps = 40, segment.size = 0.3, segment.color = "grey50",
      min.segment.length = 0.1, box.padding = 0.4, point.padding = 0.2, show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = dge_data %>% filter(label != "" & p_val_adj >= 0.05),
      aes(label = label), size = 2.8, fontface = "italic", color = "grey50",
      max.overlaps = 20, segment.size = 0.2, segment.color = "grey70",
      segment.linetype = "dashed", min.segment.length = 0.1, box.padding = 0.3, show.legend = FALSE
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = highlight_cols, name = "Gene category", drop = FALSE) +
    labs(
      x        = expression(log[2]~fold~change~(Suppressed / Unsuppressed)),
      y        = expression(-log[10]~adjusted~italic(p)~value),
      title    = paste0(cl_name, " — suppressed vs unsuppressed"),
      subtitle = "Expanding clones | positive log2FC = ↑suppressed"
    ) +
    theme_cowplot(font_size = 12) +
    theme(
      legend.position  = "right",
      legend.text      = element_text(size = 9),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  safe_cl <- gsub("[^A-Za-z0-9]", "", cl_name)
  ggsave(file.path(volcano_dir, paste0("Volcano_", safe_cl, "_Sup_vs_Unsup.png")),
         plot = p_volc, width = 12, height = 10, dpi = 300, bg = "white")
  
  cat("    Volcano saved for", cl_name, "\n")
}

# ── Volcano: all clusters combined ──────────────────────────────────────────
all_dge_path <- file.path(dge_dir, "DGE_MAST_Expanding_AllClusters_Sup_vs_Unsup.csv")
if (file.exists(all_dge_path)) {
  dge_all <- read.csv(all_dge_path, row.names = 1)
  
  dge_all$highlight <- "Other"
  dge_all$highlight[rownames(dge_all) %in% exhaustion_genes]  <- "Exhaustion"
  dge_all$highlight[rownames(dge_all) %in% stemness_genes]    <- "Stemness"
  dge_all$highlight[rownames(dge_all) %in% cytotoxic_genes]   <- "Cytotoxicity"
  dge_all$highlight[rownames(dge_all) %in% stress_genes]      <- "Stress response"
  dge_all$highlight[rownames(dge_all) %in% ifn_memory_genes]  <- "Type I IFN memory"
  dge_all$highlight[rownames(dge_all) %in% chemokine_genes]   <- "Chemokines"
  dge_all$highlight[rownames(dge_all) %in% terminal_genes]    <- "Terminal diff."
  dge_all$highlight[rownames(dge_all) %in% activation_genes]  <- "Activation"
  
  genes_to_label_all <- c(exhaustion_genes, stemness_genes, cytotoxic_genes,
                          stress_genes, ifn_memory_genes, chemokine_genes,
                          terminal_genes, activation_genes)
  genes_to_label_all <- genes_to_label_all[genes_to_label_all %in% rownames(dge_all)]
  
  dge_all$gene <- rownames(dge_all)
  dge_all$neg_log10_padj <- -log10(dge_all$p_val_adj + 1e-300)
  dge_all$sig_cat <- "ns"
  dge_all$sig_cat[dge_all$p_val_adj < 0.05 & abs(dge_all$avg_log2FC) > 0.25] <- "sig"
  
  dge_all$label <- ""
  sig_hl_all <- dge_all$gene %in% genes_to_label_all & dge_all$p_val_adj < 0.05
  dge_all$label[sig_hl_all] <- dge_all$gene[sig_hl_all]
  nonsig_hl_all <- dge_all$gene %in% genes_to_label_all & dge_all$p_val_adj >= 0.05
  dge_all$label[nonsig_hl_all] <- dge_all$gene[nonsig_hl_all]
  
  dge_all$highlight <- factor(dge_all$highlight,
                              levels = c("Exhaustion", "Stemness", "Cytotoxicity", "Stress response",
                                         "Type I IFN memory", "Chemokines", "Terminal diff.", "Activation", "Other"))
  
  p_volc_all <- ggplot(dge_all, aes(x = avg_log2FC, y = neg_log10_padj, color = highlight)) +
    geom_point(data = dge_all %>% filter(highlight == "Other"), size = 0.5, alpha = 0.3) +
    geom_point(data = dge_all %>% filter(highlight != "Other"), size = 2.5, alpha = 0.85) +
    ggrepel::geom_text_repel(
      data = dge_all %>% filter(label != "" & p_val_adj < 0.05),
      aes(label = label), size = 3.5, fontface = "italic",
      max.overlaps = 40, segment.size = 0.3, segment.color = "grey50",
      min.segment.length = 0.1, box.padding = 0.4, point.padding = 0.2, show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = dge_all %>% filter(label != "" & p_val_adj >= 0.05),
      aes(label = label), size = 2.8, fontface = "italic", color = "grey50",
      max.overlaps = 20, segment.size = 0.2, segment.color = "grey70",
      segment.linetype = "dashed", min.segment.length = 0.1, box.padding = 0.3, show.legend = FALSE
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = highlight_cols, name = "Gene category", drop = FALSE) +
    labs(
      x        = expression(log[2]~fold~change~(Suppressed / Unsuppressed)),
      y        = expression(-log[10]~adjusted~italic(p)~value),
      title    = "All effector clusters — suppressed vs unsuppressed",
      subtitle = "Expanding clones | positive log2FC = ↑suppressed"
    ) +
    theme_cowplot(font_size = 12) +
    theme(
      legend.position  = "right", legend.text = element_text(size = 9),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  ggsave(file.path(volcano_dir, "Volcano_AllClusters_Sup_vs_Unsup.png"),
         plot = p_volc_all, width = 12, height = 10, dpi = 300, bg = "white")
  cat("    Volcano saved for all clusters combined\n")
}

legend_df <- data.frame(Category = names(highlight_cols), Color = highlight_cols, stringsAsFactors = FALSE)
write.csv(legend_df, file.path(volcano_dir, "Volcano_color_legend.csv"), row.names = FALSE)

cat("  All volcano plots saved to:", volcano_dir, "\n")

# ── Export module score metadata ─────────────────────────────────────────────
write.csv(
  expand_eff@meta.data %>%
    select(CD8_Annotation, Timepoint_Group, cloneSize, all_of(score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)])),
  file.path(module_dir, "ModuleScores_per_cell.csv"),
  row.names = TRUE
)

gene_list_df <- do.call(rbind, lapply(names(module_gene_lists), function(m) {
  genes <- module_gene_lists[[m]]
  present <- genes %in% rownames(expand_eff[["RNA"]])
  data.frame(Module = m, Gene = genes, Present = present, stringsAsFactors = FALSE)
}))
write.csv(gene_list_df, file.path(module_dir, "ModuleScore_gene_lists.csv"), row.names = FALSE)

cat("\n  Module score plots and stats saved to:", module_dir, "\n")

################################################################################
# 8H (cont). MODULE-DEPENDENT PLOTS: G3, G3b, viral load, pseudotime × modules
################################################################################
cat("\n=== Running deferred module-dependent analyses ===\n")

# ── G3: Exhaustion by clone size × ART status ────────────────────────────────
exh_clone_df <- expand_eff@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters &
           !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") %>%
  mutate(
    Timepoint_Group = factor(Timepoint_Group,
                             levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed")),
    cloneSize = factor(cloneSize, levels = c(
      "Small (1 < X <= 5)", "Medium (5 < X <= 20)",
      "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 500)"))
  )

if ("Exhaustion" %in% colnames(expand_eff@meta.data)) {
  exh_clone_df$Exhaustion <- expand_eff@meta.data[rownames(exh_clone_df), "Exhaustion"]
  
  p_g3 <- ggplot(exh_clone_df,
                 aes(x = cloneSize, y = Exhaustion, fill = Timepoint_Group)) +
    geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.7,
                 position = position_dodge(width = 0.75)) +
    facet_wrap(~ CD8_Annotation, nrow = 1) +
    scale_fill_manual(values = art_colors, name = "ART Status") +
    scale_x_discrete(labels = c("Small", "Medium", "Large", "Hyper")) +
    labs(x = "Clone Size", y = "Exhaustion Module Score",
         title = "Exhaustion by clone size and ART status — controlling for expansion level") +
    theme_cowplot(font_size = 12) +
    theme(
      axis.text.x       = element_text(size = 9, angle = 30, hjust = 1),
      strip.text        = element_text(size = 11, face = "bold"),
      strip.background  = element_rect(fill = "grey95", color = NA),
      legend.position   = "bottom",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(dir_08_expansion, "G3_Exhaustion_by_CloneSize_and_ART.png"),
         plot = p_g3, width = 14, height = 6, dpi = 300, bg = "white")
  cat("  G3 saved.\n")
}

if ("Stemness" %in% colnames(expand_eff@meta.data)) {
  exh_clone_df$Stemness <- expand_eff@meta.data[rownames(exh_clone_df), "Stemness"]
  
  p_g3b <- ggplot(exh_clone_df,
                  aes(x = cloneSize, y = Stemness, fill = Timepoint_Group)) +
    geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.7,
                 position = position_dodge(width = 0.75)) +
    facet_wrap(~ CD8_Annotation, nrow = 1) +
    scale_fill_manual(values = art_colors, name = "ART Status") +
    scale_x_discrete(labels = c("Small", "Medium", "Large", "Hyper")) +
    labs(x = "Clone Size", y = "Stemness Module Score",
         title = "Stemness by clone size and ART status") +
    theme_cowplot(font_size = 12) +
    theme(
      axis.text.x       = element_text(size = 9, angle = 30, hjust = 1),
      strip.text        = element_text(size = 11, face = "bold"),
      strip.background  = element_rect(fill = "grey95", color = NA),
      legend.position   = "bottom",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(dir_08_expansion, "G3b_Stemness_by_CloneSize_and_ART.png"),
         plot = p_g3b, width = 14, height = 6, dpi = 300, bg = "white")
  cat("  G3b saved.\n")
}

# ── Viral load correlations ──────────────────────────────────────────────────
cat("  Computing viral load correlations...\n")
vl_dir <- dir_09_viralload

available_modules <- score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)]

sample_scores <- expand_eff@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters) %>%
  group_by(orig.ident, PID, Timepoint_Group, Viral_Load_num) %>%
  summarise(
    n_expanding = n(),
    across(all_of(available_modules), mean, na.rm = TRUE),
    .groups = "drop"
  )

sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)

if ("Exhaustion" %in% available_modules) {
  cor_exh <- cor.test(sample_scores$log10_VL, sample_scores$Exhaustion,
                      method = "spearman", exact = FALSE)
  
  p_h1 <- ggplot(sample_scores, aes(x = log10_VL, y = Exhaustion)) +
    geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = PID), size = 2.5, max.overlaps = 15, color = "grey40") +
    scale_color_manual(values = art_colors, name = "ART Status") +
    scale_size_continuous(name = "N expanding", range = c(2, 7)) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.2f, p = %.2g", cor_exh$estimate, cor_exh$p.value),
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(x = expression(log[10]~viral~load), y = "Mean exhaustion score",
         title = "Exhaustion vs viral load") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "right", plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(vl_dir, "H1_Exhaustion_vs_ViralLoad.png"),
         plot = p_h1, width = 10, height = 8, dpi = 300, bg = "white")
}

if ("Stemness" %in% available_modules) {
  cor_stem <- cor.test(sample_scores$log10_VL, sample_scores$Stemness,
                       method = "spearman", exact = FALSE)
  
  p_h2 <- ggplot(sample_scores, aes(x = log10_VL, y = Stemness)) +
    geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = PID), size = 2.5, max.overlaps = 15, color = "grey40") +
    scale_color_manual(values = art_colors, name = "ART Status") +
    scale_size_continuous(name = "N expanding", range = c(2, 7)) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.2f, p = %.2g", cor_stem$estimate, cor_stem$p.value),
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(x = expression(log[10]~viral~load), y = "Mean stemness score",
         title = "Stemness vs viral load") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "right", plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(vl_dir, "H2_Stemness_vs_ViralLoad.png"),
         plot = p_h2, width = 10, height = 8, dpi = 300, bg = "white")
}

# All modules vs VL
vl_long <- sample_scores %>%
  select(orig.ident, PID, Timepoint_Group, log10_VL, n_expanding, all_of(available_modules)) %>%
  tidyr::pivot_longer(cols = all_of(available_modules), names_to = "Module", values_to = "Score")

vl_long$Module_Label <- module_labels[vl_long$Module]
vl_long$Module_Label <- factor(vl_long$Module_Label, levels = module_labels)

p_h3 <- ggplot(vl_long, aes(x = log10_VL, y = Score)) +
  geom_point(aes(color = Timepoint_Group), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~ Module_Label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = art_colors, name = "ART Status") +
  labs(x = expression(log[10]~viral~load), y = "Mean module score (per sample)",
       title = "Module scores vs viral load — expanding effector clones") +
  theme_cowplot(font_size = 11) +
  theme(strip.text = element_text(size = 10), strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(vl_dir, "H3_AllModules_vs_ViralLoad.png"),
       plot = p_h3, width = 12, height = 14, dpi = 300, bg = "white")

cor_results <- lapply(available_modules, function(mod) {
  ct <- cor.test(sample_scores$log10_VL, sample_scores[[mod]], method = "spearman", exact = FALSE)
  data.frame(Module = mod, rho = ct$estimate, p_value = ct$p.value, stringsAsFactors = FALSE)
})
cor_df <- do.call(rbind, cor_results)
cor_df$p_adj <- p.adjust(cor_df$p_value, method = "BH")
write.csv(cor_df, file.path(vl_dir, "ViralLoad_Spearman_correlations.csv"), row.names = FALSE)

paired_patients <- sample_scores %>%
  group_by(PID) %>% filter(n_distinct(Timepoint_Group) > 1) %>% ungroup()

if (nrow(paired_patients) > 2 & "Exhaustion" %in% available_modules) {
  p_h4 <- ggplot(paired_patients, aes(x = log10_VL, y = Exhaustion)) +
    geom_point(aes(color = Timepoint_Group), size = 4, alpha = 0.8) +
    geom_line(aes(group = PID), color = "grey50", linewidth = 0.5, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = paste0(PID, "\n",
                                                c("PreART_Entry"="Pre","PostART_Suppressed"="Sup","PostART_Unsuppressed"="Unsup")[as.character(Timepoint_Group)])),
                             size = 3, max.overlaps = 20) +
    scale_color_manual(values = art_colors, name = "ART Status") +
    labs(x = expression(log[10]~viral~load), y = "Mean exhaustion score",
         title = "Paired patient analysis — exhaustion tracks viral load") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "right", plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(vl_dir, "H4_Paired_Exhaustion_vs_ViralLoad.png"),
         plot = p_h4, width = 10, height = 8, dpi = 300, bg = "white")
}

write.csv(sample_scores, file.path(vl_dir, "PerSample_ModuleScores_and_ViralLoad.csv"), row.names = FALSE)
cat("  Viral load analysis saved to:", vl_dir, "\n")

# ── Module scores along Monocle3 pseudotime (M6) ────────────────────────────
if (exists("traj_cells") & exists("cds")) {
  cat("  Plotting module scores along Monocle3 pseudotime...\n")
  monocle_dir <- file.path(dir_07_trajectory, "monocle3")
  
  common_cells_m3 <- intersect(colnames(traj_cells), colnames(expand_eff))
  
  if (length(common_cells_m3) > 100) {
    m3_module_df <- data.frame(
      pseudotime = traj_cells$monocle3_pseudotime[match(common_cells_m3, colnames(traj_cells))],
      Timepoint  = expand_eff@meta.data[common_cells_m3, "Timepoint_Group"],
      Cluster    = expand_eff@meta.data[common_cells_m3, "CD8_Annotation"]
    )
    for (mod in available_modules) {
      m3_module_df[[mod]] <- expand_eff@meta.data[common_cells_m3, mod]
    }
    m3_module_df <- m3_module_df %>% filter(!is.na(pseudotime) & is.finite(pseudotime))
    
    m3_exh_stem <- m3_module_df %>%
      select(pseudotime, Timepoint, Exhaustion, Stemness) %>%
      tidyr::pivot_longer(cols = c(Exhaustion, Stemness), names_to = "Module", values_to = "Score")
    
    p_m3_exh_stem <- ggplot(m3_exh_stem, aes(x = pseudotime, y = Score, color = Timepoint)) +
      geom_point(size = 0.2, alpha = 0.15) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.4, alpha = 0.2, span = 0.4) +
      facet_wrap(~ Module, nrow = 1, scales = "free_y") +
      scale_color_manual(values = art_colors, name = "ART Status") +
      labs(x = "Monocle3 Pseudotime", y = "Module Score",
           title = "Exhaustion and stemness along trajectory — expanding clones") +
      theme_cowplot(font_size = 12) +
      theme(strip.text = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "grey95", color = NA),
            legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA)) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(monocle_dir, "Monocle3_ExhaustionStemness_along_pseudotime.png"),
           plot = p_m3_exh_stem, width = 14, height = 6, dpi = 300, bg = "white")
    
    m3_all_mods <- m3_module_df %>%
      select(pseudotime, Timepoint, all_of(available_modules)) %>%
      tidyr::pivot_longer(cols = all_of(available_modules), names_to = "Module", values_to = "Score")
    m3_all_mods$Module_Label <- module_labels[m3_all_mods$Module]
    m3_all_mods$Module_Label <- factor(m3_all_mods$Module_Label, levels = module_labels)
    
    p_m3_all <- ggplot(m3_all_mods, aes(x = pseudotime, y = Score, color = Timepoint)) +
      geom_point(size = 0.1, alpha = 0.1) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2, span = 0.4) +
      facet_wrap(~ Module_Label, ncol = 2, scales = "free_y") +
      scale_color_manual(values = art_colors, name = "ART Status") +
      labs(x = "Monocle3 Pseudotime", y = "Module Score") +
      theme_cowplot(font_size = 11) +
      theme(strip.text = element_text(size = 10), legend.position = "bottom",
            plot.background = element_rect(fill = "white", color = NA)) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(monocle_dir, "Monocle3_AllModules_along_pseudotime.png"),
           plot = p_m3_all, width = 12, height = 16, dpi = 300, bg = "white")
    
    cat("  Module score × pseudotime plots saved.\n")
  }
}

################################################################################
# 9. SAVE FINAL OBJECT
################################################################################
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Final annotated HEI CD8 object saved.\n")

################################################################################
# SESSION INFO
################################################################################
cat("\n=== Done: All Figure 4-5 panels + analyses generated (HEI only) ===\n")
cat("Output directory:", out_dir, "\n")
cat("Analysis subfolders:\n")
cat("  DGE:              ", dge_dir, "\n")
cat("  DPE:              ", dpe_dir, "\n")
cat("  Proportions:      ", prop_dir, "\n")
cat("  Clonal expansion: ", clone_dir, "\n")
cat("  Module scores:    ", module_dir, "\n")
cat("  Trajectory:       ", traj_dir, "\n")
cat("  Expansion ctrl:   ", expansion_dir, "\n")
cat("  Viral load:       ", vl_dir, "\n")
sessionInfo()

