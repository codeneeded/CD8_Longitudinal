################################################################################
# FIGURE 4-5 PREP: CD8 T Cell Sub-clustering, Annotation & All Analyses
#   — TARA Cohort, HEI only (v5 — Updated Annotations)
#
# KEY CHANGES FROM v4:
#   - Data: TARA_ALL_annotated_final.qs2 (new annotation pipeline)
#   - Column: Manual_Annotation_refined → Annotation
#   - CD8 cluster names: no number prefixes
#   - Output: Analysis/CD8_subset (not Manuscript/Fig 4-5)
#   - ADT averaging: manual rowMeans (not AverageExpression — expm1 bug)
#   - Checkpoint: exports CSVs for annotation verification before proceeding
#
# Workflow:
#   1. Load TARA_ALL → filter to HEI → subset CD8 + rescue
#   2. Reprocess: RNA → ADT → WNN → UMAP → clustering  [SAVE]
#   3. Explore clusters (avg expression CSVs + heatmaps)
#   ═══ CHECKPOINT: Review CSVs, then continue ═══
#   4. Annotate
#   5. Save annotated object
#   8A-H. Downstream analyses (DGE, trajectory, modules, etc.)
#   9. Save final object
#
# Input:  TARA_ALL_annotated_final.qs2
# Output: Analysis/CD8_subset/ + saved_R_data/
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
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
# CHANGED: output to Analysis/CD8_subset
analysis_dir <- file.path(base_dir, "Analysis", "CD8_subset")

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

for (d in c(analysis_dir, dir_01_integration, dir_02_clustering, dir_03_dge,
            dir_04_dpe, dir_05_proportions, dir_06_clonal, dir_07_trajectory,
            dir_08_expansion, dir_09_viralload, dir_10_modules)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Load ──────────────────────────────────────────────────────────────────────
# CHANGED: new annotated object
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_annotated_final.qs2"))

# Filter to HEI only
TARA_ALL <- subset(TARA_ALL, Condition == "HEI")
cat("Filtered to HEI only:", ncol(TARA_ALL), "cells\n")

################################################################################
# 1. SUBSET: CD8 clusters + rescue CD8+ cells from other clusters
################################################################################

# CHANGED: new annotation names (no number prefixes)
cd8_cluster_labels <- c(
  "Naive 1 CD8",
  "Naive 2 CD8",
  "Naive Intermediate CD8",
  "Tscm CD8",
  "TEMRA/CTL CD8",
  "Tex CD8"
)

# CHANGED: uses Annotation column
in_cd8_cluster <- TARA_ALL$Annotation %in% cd8_cluster_labels

# ── Diagnostic: CD8 ADT distribution ─────────────────────────────────────────
DefaultAssay(TARA_ALL) <- "ADT"
cd8_adt_name <- "CD8A"
cd8_adt_vals <- FetchData(TARA_ALL, vars = cd8_adt_name)[, 1]

diag_df <- data.frame(
  CD8_ADT = cd8_adt_vals,
  is_cd8  = ifelse(in_cd8_cluster, "CD8 Clusters", "Other Clusters"),
  cluster = as.character(TARA_ALL$Annotation)
)

p_diag <- ggplot(diag_df, aes(x = CD8_ADT, fill = is_cd8)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "CD8 ADT Expression: CD8 clusters vs Others",
       subtitle = "Red dashed line = proposed rescue threshold",
       x = paste0(cd8_adt_name, " (DSB-normalized)"),
       y = "Density", fill = NULL) +
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(dir_01_integration, "Diagnostic_CD8_ADT_distribution.png"),
       plot = p_diag, width = 10, height = 6, dpi = 300, bg = "white")

# ── Rescue: tight multi-gate ─────────────────────────────────────────────────
cd8a_threshold <- 2.0
cd3d_threshold <- 1.0
cd4_max        <- 1.0
cd14_max       <- 0.5
cd19_max       <- 0.5

gate_data <- FetchData(TARA_ALL, vars = c("CD8A", "CD3D", "CD4", "CD14", "CD19"),
                       layer = "data")

rescue_mask <- (!in_cd8_cluster) &
  (gate_data$CD8A > cd8a_threshold) &
  (gate_data$CD3D > cd3d_threshold) &
  (gate_data$CD4  < cd4_max) &
  (gate_data$CD14 < cd14_max) &
  (gate_data$CD19 < cd19_max)

cat("Cells in CD8 clusters:", sum(in_cd8_cluster), "\n")
cat("Cells rescued from other clusters:", sum(rescue_mask), "\n")
cat("  Rescued from:\n")
# CHANGED: uses Annotation
print(table(TARA_ALL$Annotation[rescue_mask]))

# Rescue diagnostic plot
gate_data$group <- case_when(
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

# ── Combine and subset ───────────────────────────────────────────────────────
keep_cells <- in_cd8_cluster | rescue_mask
TARA_cd8 <- subset(TARA_ALL, cells = colnames(TARA_ALL)[keep_cells])
cat("\nTotal CD8 subset (before γδ/MAIT filter):", ncol(TARA_cd8), "cells\n")

# ── POST-SUBSET FILTER: Remove γδ and MAIT cells ────────────────────────────
# γδ T cells are a separate lineage (different TCR, different biology) and
# distort the WNN embedding. They are not analyzed in Figures 4-5 (trajectory,
# module scores, expanding clone DGE all use αβ CD8 only). Removing them
# before sub-clustering produces a cleaner embedding and avoids wasting
# cluster slots on populations excluded downstream.
#
# Identification:
#   γδ: TRDV1 > 0 OR TRDV2 > 0 OR TRGV9 > 0 OR TRDC > 0 (log-normalized RNA)
#       TRDC = delta chain constant region, shared by ALL γδ T cells regardless
#       of V-gene usage. αβ T cells use TRAC, not TRDC. Any detectable TRDC = γδ.
#   MAIT: TCR-vA7.2 ADT > 1.0 DSB (TRAV1-2, definitive MAIT marker)

DefaultAssay(TARA_cd8) <- "RNA"
rna_mat <- GetAssayData(TARA_cd8, assay = "RNA", layer = "data")

# γδ markers: V-genes + constant region (catches ALL γδ regardless of V-gene)
gd_genes <- c("TRDV1", "TRDV2", "TRGV9", "TRDC")
gd_expr <- list()
for (g in gd_genes) {
  if (g %in% rownames(rna_mat)) {
    gd_expr[[g]] <- as.numeric(rna_mat[g, ])
  } else {
    gd_expr[[g]] <- rep(0, ncol(TARA_cd8))
  }
}

# Any detectable expression of ANY γδ gene = γδ T cell
# TRDC = delta constant region (all γδ cells), TRDV1/TRDV2 = V-delta genes,
# TRGV9 = V-gamma gene (Vγ9Vδ2 cells)
is_gd <- (gd_expr[["TRDV1"]] > 0) | (gd_expr[["TRDV2"]] > 0) |
  (gd_expr[["TRGV9"]] > 0) | (gd_expr[["TRDC"]] > 0)

# MAIT: ADT-based (more reliable than RNA for surface markers)
DefaultAssay(TARA_cd8) <- "ADT"
mait_adt_threshold <- 1.0  # DSB
is_mait <- rep(FALSE, ncol(TARA_cd8))
if ("TCR-vA7.2" %in% rownames(TARA_cd8[["ADT"]])) {
  va72_vals <- as.numeric(GetAssayData(TARA_cd8, assay = "ADT", layer = "data")["TCR-vA7.2", ])
  is_mait <- va72_vals > mait_adt_threshold
}

exclude_mask <- is_gd | is_mait

cat("\n=== γδ/MAIT exclusion filter ===\n")
cat("  γδ (TRDV1/TRDV2/TRGV9/TRDC > 0):", sum(is_gd), "cells\n")
cat("    Per gene:\n")
for (g in gd_genes) {
  cat(sprintf("      %s > 0: %d cells\n", g, sum(gd_expr[[g]] > 0)))
}
cat("  MAIT (TCR-vA7.2 ADT >", mait_adt_threshold, "):", sum(is_mait), "cells\n")
cat("  Overlap (γδ AND MAIT):", sum(is_gd & is_mait), "\n")
cat("  Total excluded:", sum(exclude_mask), "\n")

# Print which PBMC annotations the excluded cells came from
cat("\n  Excluded cells by PBMC annotation:\n")
print(table(TARA_cd8$Annotation[exclude_mask]))

# Also check: are we losing any cells with αβ TCR clonotype data?
if ("CTstrict" %in% colnames(TARA_cd8@meta.data)) {
  has_clone <- !is.na(TARA_cd8$CTstrict) & TARA_cd8$CTstrict != ""
  cat("\n  Of excluded cells, have αβ clonotype:", sum(exclude_mask & has_clone), "\n")
  cat("  (These are likely dual-expressors or misassigned — acceptable loss)\n")
}

# Apply filter
TARA_cd8 <- subset(TARA_cd8, cells = colnames(TARA_cd8)[!exclude_mask])
cat("\nFinal αβ CD8 subset:", ncol(TARA_cd8), "cells\n")
cat("  Cells removed:", sum(exclude_mask), "(", 
    round(sum(exclude_mask) / sum(keep_cells) * 100, 1), "%)\n\n")

rm(TARA_ALL); gc()

################################################################################
# 2. REPROCESS + INTEGRATE: RNA → ADT → WNN → UMAP → Clustering
################################################################################

SEED <- 42
set.seed(SEED)

# ── RNA: Split → Normalize → VariableFeatures → Scale → PCA → FastMNN ───────
DefaultAssay(TARA_cd8) <- "RNA"
TARA_cd8[["RNA"]] <- split(TARA_cd8[["RNA"]], f = TARA_cd8$orig.ident)
TARA_cd8 <- NormalizeData(TARA_cd8, assay = "RNA")
TARA_cd8 <- FindVariableFeatures(TARA_cd8, assay = "RNA",
                                 selection.method = "vst", nfeatures = 3000)
TARA_cd8 <- ScaleData(TARA_cd8, assay = "RNA")
TARA_cd8 <- RunPCA(TARA_cd8, assay = "RNA", reduction.name = "pca")

TARA_cd8 <- IntegrateLayers(
  TARA_cd8, method = FastMNNIntegration, assay = "RNA",
  new.reduction = "integrated.mnn.rna"
)
TARA_cd8[["RNA"]] <- JoinLayers(TARA_cd8[["RNA"]])

# Diagnostic UMAPs
TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "pca", dims = 1:30, seed.use = SEED,
                    reduction.name = "umap.unintegrated", reduction.key = "unintUMAP_")
TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "integrated.mnn.rna", dims = 1:30, seed.use = SEED,
                    reduction.name = "umap.mnn.rna", reduction.key = "mnnUMAP_")

p1 <- DimPlot2(TARA_cd8, reduction = "umap.unintegrated", group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("RNA Unintegrated")
p2 <- DimPlot2(TARA_cd8, reduction = "umap.mnn.rna", group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("RNA FastMNN")
ggsave(file.path(dir_01_integration, "CD8_RNA_integration_by_sample.png"),
       plot = p1 | p2, width = 22, height = 8, dpi = 300, bg = "white")

p3 <- DimPlot2(TARA_cd8, reduction = "umap.unintegrated", group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("RNA Unintegrated — by Timepoint")
p4 <- DimPlot2(TARA_cd8, reduction = "umap.mnn.rna", group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("RNA FastMNN — by Timepoint")
ggsave(file.path(dir_01_integration, "CD8_RNA_integration_by_timepoint_group.png"),
       plot = p3 | p4, width = 22, height = 8, dpi = 300, bg = "white")

cat("RNA integration diagnostics saved.\n")

# ── ADT: Split → Scale → PCA → FastMNN ──────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"
TARA_cd8[["ADT"]] <- split(TARA_cd8[["ADT"]], f = TARA_cd8$orig.ident)
VariableFeatures(TARA_cd8) <- rownames(TARA_cd8[["ADT"]])
TARA_cd8 <- ScaleData(TARA_cd8, assay = "ADT",
                      features = VariableFeatures(TARA_cd8), verbose = FALSE)
TARA_cd8 <- RunPCA(TARA_cd8, assay = "ADT",
                   features = VariableFeatures(TARA_cd8),
                   reduction.name = "apca")

TARA_cd8 <- IntegrateLayers(
  TARA_cd8, method = FastMNNIntegration, assay = "ADT",
  new.reduction = "integrated.mnn.adt"
)
TARA_cd8[["ADT"]] <- JoinLayers(TARA_cd8[["ADT"]])

# ADT diagnostic UMAPs
TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "apca", dims = 1:20, seed.use = SEED,
                    reduction.name = "umap.adt.unint", reduction.key = "adtunintUMAP_")
TARA_cd8 <- RunUMAP(TARA_cd8, reduction = "integrated.mnn.adt", dims = 1:20, seed.use = SEED,
                    reduction.name = "umap.mnn.adt", reduction.key = "mnnADTUMAP_")

p5 <- DimPlot2(TARA_cd8, reduction = "umap.adt.unint", group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("ADT Unintegrated")
p6 <- DimPlot2(TARA_cd8, reduction = "umap.mnn.adt", group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("ADT FastMNN")
ggsave(file.path(dir_01_integration, "CD8_ADT_integration_by_sample.png"),
       plot = p5 | p6, width = 22, height = 8, dpi = 300, bg = "white")

p7 <- DimPlot2(TARA_cd8, reduction = "umap.adt.unint", group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("ADT Unintegrated — by Timepoint")
p8 <- DimPlot2(TARA_cd8, reduction = "umap.mnn.adt", group.by = "Timepoint_Group", pt.size = 0.3) +
  ggtitle("ADT FastMNN — by Timepoint")
ggsave(file.path(dir_01_integration, "CD8_ADT_integration_by_timepoint_group.png"),
       plot = p7 | p8, width = 22, height = 8, dpi = 300, bg = "white")

cat("ADT integration diagnostics saved.\n")

# ── WNN → UMAP → Clustering ─────────────────────────────────────────────────
TARA_cd8 <- FindMultiModalNeighbors(
  TARA_cd8,
  reduction.list = list("integrated.mnn.rna", "integrated.mnn.adt"),
  dims.list      = list(1:30, 1:20)
)

TARA_cd8 <- RunUMAP(TARA_cd8, nn.name = "weighted.nn",
                    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                    seed.use = SEED)

resolutions <- seq(0.2, 2.0, by = 0.2)
for (res in resolutions) {
  TARA_cd8 <- FindClusters(TARA_cd8, graph.name = "wsnn", resolution = res,
                           algorithm = 3, random.seed = SEED, verbose = FALSE)
}

for (col in grep("^wsnn_res\\.", colnames(TARA_cd8@meta.data), value = TRUE)) {
  lvls <- as.character(sort(as.numeric(levels(TARA_cd8[[col]][, 1]))))
  TARA_cd8[[col]] <- factor(TARA_cd8[[col]][, 1], levels = lvls)
}

# Clustree
p_ct <- clustree(TARA_cd8, prefix = "wsnn_res.")
ggsave(file.path(dir_02_clustering, "CD8_clustree.png"),
       plot = p_ct, width = 15, height = 9, dpi = 300, bg = "white")

p_ct2 <- clustree(TARA_cd8, prefix = "wsnn_res.", node_colour = "sc3_stability")
ggsave(file.path(dir_02_clustering, "CD8_clustree_stability.png"),
       plot = p_ct2, width = 15, height = 9, dpi = 300, bg = "white")

# Pick resolution
chosen_res <- 0.4
TARA_cd8$seurat_clusters <- TARA_cd8[[paste0("wsnn_res.", chosen_res)]][, 1]
Idents(TARA_cd8) <- "seurat_clusters"

cat("Chosen resolution:", chosen_res, "\n")
cat("Number of clusters:", length(levels(Idents(TARA_cd8))), "\n")
print(table(Idents(TARA_cd8)))

# Resolution comparison UMAPs
for (res in c(0.4, 0.6, 0.8, 1.0)) {
  p <- DimPlot2(TARA_cd8, reduction = "wnn.umap",
                group.by = paste0("wsnn_res.", res),
                cols = "light", label = TRUE, box = TRUE, repel = TRUE,
                label.color = "black", pt.size = 0.6) +
    ggtitle(paste0("Resolution: ", res))
  ggsave(file.path(dir_02_clustering, paste0("CD8_UMAP_res", res, ".png")),
         plot = p, width = 8, height = 7, dpi = 300, bg = "white")
}

# CHECKPOINT: Save after integration
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_integrated.qs2"))
cat("CHECKPOINT: Integrated object saved.\n")


################################################################################
# 3. ANNOTATION EXPLORATION: Average expression per cluster
#    FIXED: ADT uses manual rowMeans (not AverageExpression which applies expm1)
################################################################################

annotation_markers_rna <- c(
  "CCR7", "SELL", "TCF7", "LEF1", "IL7R",
  "CD27", "CD28", "BCL2", "BACH2",
  "GZMK", "EOMES", "CXCR3",
  "GZMB", "GNLY", "PRF1", "NKG7", "TBX21", "CX3CR1",
  "TOX", "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CXCL13",
  "MKI67", "TOP2A", "STMN1",
  "ITGAE", "CXCR6", "CD69",
  "TRDV1", "TRGV9", "TRDC",
  "TYROBP", "KLRD1", "FCGR3A", "NCAM1",
  "HLA-DRA", "CD38"
)

annotation_markers_adt <- c(
  "CD45RA", "SELL", "CD7",
  "CD45RO", "IL7R", "CD27", "CD28",
  "B3GAT1", "KLRG1", "KIR3DL1",
  "TIGIT", "PDCD1", "LAG3",
  "CD38", "ENTPD1", "NT5E",
  "NCAM1", "FCGR3A", "SIGLEC7",
  "TCR-AB", "TCR-vA7.2", "TCR-vD2"
)

# ── RNA: AverageExpression is fine for log-normalized RNA ────────────────────
DefaultAssay(TARA_cd8) <- "RNA"
rna_feats <- annotation_markers_rna[annotation_markers_rna %in% rownames(TARA_cd8[["RNA"]])]
avg_rna <- AverageExpression(
  TARA_cd8, assays = "RNA", features = rna_feats,
  group.by = "seurat_clusters", slot = "data"
)$RNA

# ── ADT: MANUAL rowMeans from data slot ──────────────────────────────────────
# NOTE: Seurat v5 AverageExpression() applies expm1() to DSB-normalized ADT,
# producing exponentially inflated values. Use manual rowMeans instead.
DefaultAssay(TARA_cd8) <- "ADT"
adt_feats <- annotation_markers_adt[annotation_markers_adt %in% rownames(TARA_cd8[["ADT"]])]
adt_mat_raw <- GetAssayData(TARA_cd8, slot = "data")[adt_feats, , drop = FALSE]
clusters_vec <- TARA_cd8$seurat_clusters

avg_adt <- sapply(levels(factor(clusters_vec)), function(cl) {
  cells <- which(clusters_vec == cl)
  if (length(cells) == 0) return(rep(0, length(adt_feats)))
  rowMeans(adt_mat_raw[, cells, drop = FALSE])
})
rownames(avg_adt) <- adt_feats

# ── Scale for heatmap ────────────────────────────────────────────────────────
scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)))
}

avg_rna_scaled <- scale_01(avg_rna)
# CHANGED: no log10 for DSB ADT — already on meaningful scale
avg_adt_scaled <- scale_01(as.matrix(avg_adt))

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
pheatmap(avg_rna_scaled, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heatmap_colors, border_color = "white",
         cellwidth = 40, cellheight = 14, fontsize = 10,
         main = "RNA Average Expression (scaled) per CD8 Sub-cluster")
dev.off()

png(file.path(dir_02_clustering, "CD8_AvgExpression_ADT_exploration.png"),
    width = 14, height = 12, units = "in", res = 300, bg = "white")
pheatmap(avg_adt_scaled, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heatmap_colors, border_color = "white",
         cellwidth = 40, cellheight = 14, fontsize = 10,
         main = "ADT Average Expression (DSB, scaled) per CD8 Sub-cluster")
dev.off()

# ── Export raw CSVs for verification ─────────────────────────────────────────
checkpoint_dir <- file.path(analysis_dir, "02B_annotation_checkpoint")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(round(as.data.frame(avg_rna), 4),
          file.path(checkpoint_dir, "AvgExpr_RNA_raw.csv"))
write.csv(round(as.data.frame(avg_adt), 4),
          file.path(checkpoint_dir, "AvgExpr_ADT_raw_DSB.csv"))
write.csv(round(avg_rna_scaled, 4),
          file.path(checkpoint_dir, "AvgExpr_RNA_scaled.csv"))
write.csv(round(avg_adt_scaled, 4),
          file.path(checkpoint_dir, "AvgExpr_ADT_scaled.csv"))

# Cluster composition summary
cluster_summary <- TARA_cd8@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    pct_preART = round(sum(Timepoint_Group == "PreART_Entry") / n() * 100, 1),
    pct_suppressed = round(sum(Timepoint_Group == "PostART_Suppressed") / n() * 100, 1),
    pct_unsuppressed = round(sum(Timepoint_Group == "PostART_Unsuppressed") / n() * 100, 1),
    pct_expanded = round(sum(!is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") / n() * 100, 1),
    pct_TCR = round(sum(has_TCR == TRUE, na.rm = TRUE) / n() * 100, 1),
    .groups = "drop"
  )
write.csv(cluster_summary,
          file.path(checkpoint_dir, "Cluster_Composition_Summary.csv"), row.names = FALSE)

cat("\n=== RNA Average Expression (scaled) ===\n")
print(round(avg_rna_scaled, 2))
cat("\n=== ADT Average Expression (DSB, scaled) ===\n")
print(round(avg_adt_scaled, 2))
cat("\n=== Cluster Composition ===\n")
print(as.data.frame(cluster_summary))

cat("\n")
cat("================================================================\n")
cat("  CHECKPOINT: Annotation verification files exported\n")
cat("  Location:", checkpoint_dir, "\n")
cat("\n")
cat("  Files:\n")
cat("    AvgExpr_RNA_raw.csv        — raw log-normalized RNA averages\n")
cat("    AvgExpr_ADT_raw_DSB.csv    — raw DSB-normalized ADT averages\n")
cat("    AvgExpr_RNA_scaled.csv     — [0,1] scaled RNA\n")
cat("    AvgExpr_ADT_scaled.csv     — [0,1] scaled ADT (no log10)\n")
cat("    Cluster_Composition_Summary.csv — n_cells, %preART, %expanded, etc.\n")
cat("\n")
cat("  Upload these CSVs for annotation verification.\n")
cat("  Then continue from Section 4 below.\n")
cat("================================================================\n")

# ▼▼▼ UNCOMMENT TO STOP HERE ▼▼▼
# stop("CHECKPOINT — review CSVs before continuing to annotation")

# To resume:
#TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_integrated.qs2"))

################################################################################
# 4. ANNOTATION: αβ CD8, Resolution 0.4 (10 clusters: 0–9)
#    WITH ADT SUB-GATING OF NAÏVE CLUSTERS
#
#    POST γδ/MAIT REMOVAL — all clusters are αβ CD8 lineage.
#    No TCR-based splits needed.
#
#    CLUSTER MAP:
#      0 → Naïve CD8          (classical naive, 64% preART)
#      1 → TEMRA CD8          (terminal effector, B3GAT1+CX3CR1+, 69% expanded)
#      2 → Naïve CD8 2        (preART-dominated, CD38↑CD7↑, 90% preART)
#      3 → TEM CD8            (GZMK+CXCR3+TIGIT+, cycling, 29% expanded)
#      4 → Naïve CD8 3        (BACH2+BCL2+, quiescent, 51% suppressed)
#      5 → CD16+ Effector CD8 (FCGR3A↑, ADCC-capable, 41% expanded)
#      6 → Naïve CD8 4        (naive + GZMK low, STMN1↑, 47% suppressed)
#      7 → KIR+ innate-like CD8 (KIR3DL1↑↑, 65% suppressed)
#      8 → REMOVE (APC/DC contaminant: HLA-DRA↑↑ TYROBP↑↑ SIGLEC7↑↑)
#      9 → REMOVE (59 cells, 6.8% TCR, unresolvable)
#
#    ADT SUB-GATING (naïve clusters 0, 2, 4, 6):
#      CD45RO+             → Naïve Intermediate CD8
#      CD45RO− AND FAS+    → Tscm CD8
#      CD45RO− AND FAS−    → retains original naïve identity
#
#    FINAL ANNOTATION SET (10 populations):
#      Naïve CD8, Naïve CD8 2, Naïve CD8 3, Naïve CD8 4,
#      Tscm CD8, Naïve Intermediate CD8,
#      TEM CD8, TEMRA CD8, CD16+ Effector CD8,
#      KIR+ innate-like CD8
################################################################################

# ── Step 1: Remove contaminants ──────────────────────────────────────────────
clusters_to_remove <- c("8", "9")
n_before <- ncol(TARA_cd8)
TARA_cd8 <- subset(TARA_cd8, seurat_clusters %in% clusters_to_remove, invert = TRUE)
TARA_cd8$seurat_clusters <- droplevels(TARA_cd8$seurat_clusters)

cat("Removed clusters:", paste(clusters_to_remove, collapse = ", "), "\n")
cat("  Cells before:", n_before, "→ after:", ncol(TARA_cd8), "\n")
cat("  Removed:", n_before - ncol(TARA_cd8), "cells\n")
cat("  Remaining clusters:", paste(levels(TARA_cd8$seurat_clusters), collapse = ", "), "\n\n")

# ── Step 2: Initial annotation ───────────────────────────────────────────────
cd8_annotations <- c(
  "0" = "Naïve CD8",
  "1" = "TEMRA CD8",
  "2" = "Naïve CD8 2",
  "3" = "TEM CD8",
  "4" = "Naïve CD8 3",
  "5" = "CD16+ Effector CD8",
  "6" = "Naïve CD8 4",
  "7" = "KIR+ innate-like CD8"
)

annot_vec <- cd8_annotations[as.character(TARA_cd8$seurat_clusters)]
names(annot_vec) <- colnames(TARA_cd8)
TARA_cd8 <- AddMetaData(TARA_cd8, metadata = annot_vec, col.name = "CD8_Annotation")

cat("Initial annotation applied:\n")
print(sort(table(TARA_cd8$CD8_Annotation), decreasing = TRUE))
cat("\n")

################################################################################
# Step 3: ADT SUB-GATING OF NAÏVE CLUSTERS
#
# Gate each naïve population independently. Tscm and Naïve Intermediate
# are pooled across source clusters; true naïve cells keep their identity.
#
# Thresholds (same as Fig 1):
#   CD45RO > 0  → Naïve Intermediate CD8
#   CD45RO ≤ 0 AND FAS > 0 → Tscm CD8
#   CD45RO ≤ 0 AND FAS ≤ 0 → keep original naïve label
################################################################################

cat(strrep("═", 70), "\n")
cat("  ADT SUB-GATING: Recovering Tscm from naïve clusters\n")
cat(strrep("═", 70), "\n\n")

# Store source cluster before gating
TARA_cd8$source_cluster <- as.character(TARA_cd8$seurat_clusters)

# All four naïve populations
naive_labels <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3", "Naïve CD8 4")
naive_mask <- TARA_cd8$CD8_Annotation %in% naive_labels
naive_cells <- colnames(TARA_cd8)[naive_mask]

cat("Total naïve-lineage cells to sub-gate:", length(naive_cells), "\n")
cat("  By current annotation:\n")
print(table(TARA_cd8$CD8_Annotation[naive_mask]))
cat("\n")

# Extract ADT expression
DefaultAssay(TARA_cd8) <- "ADT"
adt_data <- GetAssayData(TARA_cd8, slot = "data")

stopifnot("CD45RO not in ADT" = "CD45RO" %in% rownames(adt_data))
stopifnot("FAS not in ADT"    = "FAS" %in% rownames(adt_data))

cd45ro_vals <- as.numeric(adt_data["CD45RO", naive_cells])
fas_vals    <- as.numeric(adt_data["FAS", naive_cells])
cd45ra_vals <- as.numeric(adt_data["CD45RA", naive_cells])
current_labels <- TARA_cd8$CD8_Annotation[naive_cells]

# Thresholds
cd45ro_threshold <- 0
fas_threshold    <- 0

cat("Gating thresholds (same as Fig 1):\n")
cat("  CD45RO:", cd45ro_threshold, "\n")
cat("  FAS:   ", fas_threshold, "\n\n")

# Diagnostics per naïve population
for (lbl in naive_labels) {
  lbl_mask <- current_labels == lbl
  if (sum(lbl_mask) == 0) next
  lbl_ro  <- cd45ro_vals[lbl_mask]
  lbl_fas <- fas_vals[lbl_mask]
  cat(sprintf("  %s (%d cells):\n", lbl, sum(lbl_mask)))
  cat(sprintf("    CD45RO: median=%.3f, %%>0=%.1f%%\n",
              median(lbl_ro), mean(lbl_ro > 0) * 100))
  cat(sprintf("    FAS:    median=%.3f, %%>0=%.1f%%\n",
              median(lbl_fas), mean(lbl_fas > 0) * 100))
}
cat("\n")

# Apply gating
cd45ro_pos <- cd45ro_vals > cd45ro_threshold
fas_pos    <- fas_vals > fas_threshold

gate_result <- case_when(
  cd45ro_pos                  ~ "Naïve Intermediate CD8",
  !cd45ro_pos &  fas_pos      ~ "Tscm CD8",
  !cd45ro_pos & !fas_pos      ~ as.character(current_labels)
)
names(gate_result) <- naive_cells

# Print results
cat("=== Overall gating results ===\n")
gt <- table(gate_result)
gp <- round(prop.table(gt) * 100, 1)
for (g in sort(names(gt))) {
  cat(sprintf("  %-30s %5d (%5.1f%%)\n", g, gt[g], gp[g]))
}
cat("\n")

cat("=== Gating results per source population ===\n")
gate_df <- data.frame(
  cell         = naive_cells,
  source_label = as.character(current_labels),
  source_clust = as.character(TARA_cd8$source_cluster[naive_cells]),
  gate         = gate_result,
  stringsAsFactors = FALSE
)

for (lbl in naive_labels) {
  lbl_gates <- gate_df %>% filter(source_label == lbl)
  if (nrow(lbl_gates) == 0) next
  cat(sprintf("  %s (%d cells):\n", lbl, nrow(lbl_gates)))
  lbl_gt <- table(lbl_gates$gate)
  lbl_gp <- round(prop.table(lbl_gt) * 100, 1)
  for (g in sort(names(lbl_gt))) {
    cat(sprintf("    %-28s %5d (%5.1f%%)\n", g, lbl_gt[g], lbl_gp[g]))
  }
}
cat("\n")

# Apply to metadata
TARA_cd8$CD8_Annotation[match(names(gate_result), colnames(TARA_cd8))] <- gate_result

# Store gate info
adt_gate_vec <- rep(NA_character_, ncol(TARA_cd8))
names(adt_gate_vec) <- colnames(TARA_cd8)
adt_gate_vec[names(gate_result)] <- gate_result
TARA_cd8 <- AddMetaData(TARA_cd8, metadata = adt_gate_vec, col.name = "ADT_gate")

# Gate × Timepoint
cat("=== Gating results × Timepoint_Group ===\n")
print(table(gate_result, TARA_cd8$Timepoint_Group[match(names(gate_result), colnames(TARA_cd8))]))
cat("\n")

# ── Step 4: Finalize ─────────────────────────────────────────────────────────
TARA_cd8$CD8_Annotation <- factor(TARA_cd8$CD8_Annotation)
Idents(TARA_cd8) <- "CD8_Annotation"

cat("=== Final CD8 sub-cluster sizes ===\n")
print(sort(table(TARA_cd8$CD8_Annotation), decreasing = TRUE))

cat("\n=== CD8_Annotation × Timepoint_Group ===\n")
print(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group))

cat("\n=== Tscm CD8 breakdown ===\n")
tscm_mask <- TARA_cd8$CD8_Annotation == "Tscm CD8"
if (sum(tscm_mask) > 0) {
  cat("  By source cluster:\n")
  print(table(TARA_cd8$source_cluster[tscm_mask]))
  cat("  By timepoint:\n")
  print(table(TARA_cd8$Timepoint_Group[tscm_mask]))
} else {
  cat("  WARNING: No Tscm CD8 cells found.\n")
}

cat("\n=== Naïve Intermediate CD8 breakdown ===\n")
naint_mask <- TARA_cd8$CD8_Annotation == "Naïve Intermediate CD8"
if (sum(naint_mask) > 0) {
  cat("  By source cluster:\n")
  print(table(TARA_cd8$source_cluster[naint_mask]))
  cat("  By timepoint:\n")
  print(table(TARA_cd8$Timepoint_Group[naint_mask]))
}

################################################################################
# Step 5: GATING DIAGNOSTIC PLOTS
################################################################################

cat("\n=== Generating gating diagnostic plots ===\n")

gate_plot_dir <- file.path(analysis_dir, "03B_naive_subgating")
dir.create(gate_plot_dir, recursive = TRUE, showWarnings = FALSE)

gate_plot_df <- data.frame(
  cell         = naive_cells,
  source_label = as.character(current_labels),
  source_clust = as.character(TARA_cd8$source_cluster[naive_cells]),
  CD45RO       = cd45ro_vals,
  CD45RA       = cd45ra_vals,
  FAS          = fas_vals,
  gate         = gate_result[naive_cells],
  Timepoint    = TARA_cd8$Timepoint_Group[naive_cells],
  stringsAsFactors = FALSE
)

gate_cols <- c(
  "Naïve CD8"              = "#4E79A7",
  "Naïve CD8 2"            = "#76B7B2",
  "Naïve CD8 3"            = "#B07AA1",
  "Naïve CD8 4"            = "#EDC948",
  "Tscm CD8"               = "#59A14F",
  "Naïve Intermediate CD8" = "#E15759"
)

source_labels_map <- c(
  "Naïve CD8"   = "Naïve CD8\n(classical, Cl.0)",
  "Naïve CD8 2" = "Naïve CD8 2\n(preART, Cl.2)",
  "Naïve CD8 3" = "Naïve CD8 3\n(BACH2+, Cl.4)",
  "Naïve CD8 4" = "Naïve CD8 4\n(GZMK-lo, Cl.6)"
)

# Plot 1: CD45RO vs FAS
p1 <- ggplot(gate_plot_df, aes(x = FAS, y = CD45RO, color = gate)) +
  geom_point(size = 0.4, alpha = 0.4) +
  scale_color_manual(values = gate_cols, name = "ADT Gate") +
  geom_hline(yintercept = cd45ro_threshold, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = fas_threshold, linetype = "dashed", color = "grey40") +
  facet_wrap(~ source_label, nrow = 1,
             labeller = labeller(source_label = source_labels_map)) +
  labs(title = "CD45RO vs FAS/CD95 — Tscm = CD45RO- FAS+",
       x = "FAS/CD95 (DSB)", y = "CD45RO (DSB)") +
  theme_cowplot(font_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(gate_plot_dir, "01_Gate_CD45RO_vs_FAS.png"),
       p1, width = 18, height = 6, dpi = 300, bg = "white")

# Plot 2: Histograms
p_h1 <- ggplot(gate_plot_df, aes(x = CD45RO, fill = gate)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = gate_cols) +
  geom_vline(xintercept = cd45ro_threshold, linetype = "dashed") +
  facet_wrap(~ source_label, nrow = 1, scales = "free_y",
             labeller = labeller(source_label = source_labels_map)) +
  labs(title = "CD45RO distribution", x = "CD45RO (DSB)") +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

p_h2 <- ggplot(gate_plot_df, aes(x = FAS, fill = gate)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = gate_cols) +
  geom_vline(xintercept = fas_threshold, linetype = "dashed") +
  facet_wrap(~ source_label, nrow = 1, scales = "free_y",
             labeller = labeller(source_label = source_labels_map)) +
  labs(title = "FAS/CD95 distribution", x = "FAS (DSB)") +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA))

p_hists <- p_h1 / p_h2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(gate_plot_dir, "02_Gate_Histograms.png"),
       p_hists, width = 18, height = 10, dpi = 300, bg = "white")

# Export gating CSV
write.csv(gate_plot_df %>%
            select(cell, source_label, source_clust, CD45RO, CD45RA, FAS, gate, Timepoint),
          file.path(gate_plot_dir, "03_Gating_per_cell.csv"), row.names = FALSE)

cat("  Gating diagnostic plots saved to:", gate_plot_dir, "\n\n")

################################################################################
# Step 6: ANNOTATION VERIFICATION
################################################################################

cat("\n")
cat(strrep("█", 70), "\n")
cat("  ANNOTATION VERIFICATION\n")
cat(strrep("█", 70), "\n\n")

# V1: Cross-table
cat("=== V1: seurat_clusters × CD8_Annotation ===\n")
print(table(TARA_cd8$seurat_clusters, TARA_cd8$CD8_Annotation))

# V2: Per-annotation summary
cat("\n=== V2: Per-annotation biological verification ===\n")
DefaultAssay(TARA_cd8) <- "RNA"
verify_df <- TARA_cd8@meta.data %>%
  group_by(CD8_Annotation) %>%
  summarise(
    n_cells        = n(),
    pct_preART     = round(sum(Timepoint_Group == "PreART_Entry") / n() * 100, 1),
    pct_suppressed = round(sum(Timepoint_Group == "PostART_Suppressed") / n() * 100, 1),
    pct_unsuppressed = round(sum(Timepoint_Group == "PostART_Unsuppressed") / n() * 100, 1),
    pct_expanded   = round(sum(!is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") / n() * 100, 1),
    pct_TCR        = round(sum(has_TCR == TRUE, na.rm = TRUE) / n() * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

cat(sprintf("%-28s %7s %8s %8s %8s %8s %8s\n",
            "Annotation", "N", "%preART", "%Sup", "%Unsup", "%Expand", "%TCR"))
cat(strrep("-", 95), "\n")
for (i in seq_len(nrow(verify_df))) {
  r <- verify_df[i, ]
  cat(sprintf("%-28s %7d %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n",
              r$CD8_Annotation, r$n_cells, r$pct_preART, r$pct_suppressed,
              r$pct_unsuppressed, r$pct_expanded, r$pct_TCR))
}

# V3: Diagnostic UMAPs
cat("\n=== V3: Generating diagnostic UMAPs ===\n")
diag_umap_dir <- file.path(analysis_dir, "04_annotation_verification")
dir.create(diag_umap_dir, recursive = TRUE, showWarnings = FALSE)

p_v1 <- DimPlot(TARA_cd8, reduction = "wnn.umap", group.by = "CD8_Annotation",
                label = TRUE, label.size = 3.5, repel = TRUE, pt.size = 0.4) +
  ggtitle("CD8_Annotation (final)") + NoAxes() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 8))

ggsave(file.path(diag_umap_dir, "V1_UMAP_CD8_Annotation.png"),
       plot = p_v1, width = 14, height = 10, dpi = 200, bg = "white")

art_colors <- c("PreART_Entry" = "#4A90D9", "PostART_Suppressed" = "#52B788",
                "PostART_Unsuppressed" = "#E76F51")

p_v2 <- DimPlot(TARA_cd8, reduction = "wnn.umap", group.by = "Timepoint_Group",
                pt.size = 0.4) +
  scale_color_manual(values = art_colors) +
  ggtitle("ART status") + NoAxes() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(diag_umap_dir, "V2_UMAP_ART_status.png"),
       plot = p_v2, width = 12, height = 10, dpi = 200, bg = "white")

# Key marker feature plots
DefaultAssay(TARA_cd8) <- "RNA"
for (feat in c("TCF7", "GZMK", "GZMB", "BACH2", "TOX", "FCGR3A")) {
  if (!feat %in% rownames(TARA_cd8[["RNA"]])) next
  p <- FeaturePlot(TARA_cd8, features = feat, reduction = "wnn.umap",
                   pt.size = 0.3, order = TRUE) +
    scale_color_viridis_c(option = "magma") + NoAxes() +
    ggtitle(feat) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(diag_umap_dir, paste0("V3_", feat, ".png")),
         plot = p, width = 8, height = 7, dpi = 200, bg = "white")
}

DefaultAssay(TARA_cd8) <- "ADT"
for (feat in c("KIR3DL1", "FAS", "CD45RO", "B3GAT1")) {
  if (!feat %in% rownames(TARA_cd8[["ADT"]])) next
  p <- FeaturePlot(TARA_cd8, features = feat, reduction = "wnn.umap",
                   pt.size = 0.3, order = TRUE) +
    scale_color_viridis_c(option = "inferno") + NoAxes() +
    ggtitle(paste0(feat, " (ADT)")) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(diag_umap_dir, paste0("V3_ADT_", feat, ".png")),
         plot = p, width = 8, height = 7, dpi = 200, bg = "white")
}

cat("  Diagnostic UMAPs saved to:", diag_umap_dir, "\n")

cat("\n")
cat(strrep("█", 70), "\n")
cat("  VERIFICATION COMPLETE — review output before proceeding\n")
cat(strrep("█", 70), "\n\n")

################################################################################
# Canonical cluster order + analysis group definitions
################################################################################
col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3", "Naïve CD8 4",
  "Tscm CD8", "Naïve Intermediate CD8",
  "TEM CD8", "TEMRA CD8", "CD16+ Effector CD8",
  "KIR+ innate-like CD8"
)

# Three effector clusters for expanding clone DGE / module score analyses
effector_clusters <- c("TEM CD8", "TEMRA CD8", "CD16+ Effector CD8", "KIR+ innate-like CD8")

# Naïve clusters (for pairwise DGE in section 8A)
naive_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3", "Naïve CD8 4",
                    "Tscm CD8", "Naïve Intermediate CD8")

# All clusters for trajectory (no γδ/MAIT exclusion needed — already removed)
ab_clusters <- col_order_cd8

################################################################################
# 5. SAVE ANNOTATED OBJECT
################################################################################
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_annotated.qs2"))
cat("Annotated HEI CD8 object saved.\n")
cat("  Cells:", ncol(TARA_cd8), "\n")
cat("  Annotations:", length(levels(TARA_cd8$CD8_Annotation)), "\n")
cat("  Effector clusters:", paste(effector_clusters, collapse = ", "), "\n\n")
#TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated.qs2"))

################################################################################
# 8. DOWNSTREAM ANALYSES
#    Updated for v5 annotation: no γδ/MAIT, CD16+ Effector CD8 replaces
#    Transitional Tem CD8 as third effector cluster.
################################################################################

# ── Shared aesthetics ────────────────────────────────────────────────────────
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

# These were defined in Section 4 but restated here for clarity if resuming
effector_clusters <- c("TEM CD8", "TEMRA CD8", "CD16+ Effector CD8", "KIR+ innate-like CD8")

col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3", "Naïve CD8 4",
  "Tscm CD8", "Naïve Intermediate CD8",
  "TEM CD8", "TEMRA CD8", "CD16+ Effector CD8",
  "KIR+ innate-like CD8"
)

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

dge_dir   <- dir_03_dge
dpe_dir   <- dir_04_dpe
prop_dir  <- dir_05_proportions
clone_dir <- dir_06_clonal

################################################################################
# 8A. DGE/DPE: Pairwise naïve cluster comparisons
################################################################################
cat("\n=== Running DGE/DPE (MAST) for naïve cluster comparisons ===\n")

Idents(TARA_cd8) <- "CD8_Annotation"

naive_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3", "Naïve CD8 4",
                    "Tscm CD8", "Naïve Intermediate CD8")

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

# TEM vs TEMRA
cat("  DGE: TEM CD8 vs TEMRA CD8...\n")
tryCatch({
  DefaultAssay(TARA_cd8) <- "RNA"
  dge_eff <- FindMarkers(TARA_cd8, ident.1 = "TEM CD8", ident.2 = "TEMRA CD8",
                         test.use = "MAST", logfc.threshold = 0.25,
                         min.pct = 0.1, verbose = FALSE)
  write.csv(dge_eff, file.path(dge_dir, "DGE_MAST_RNA_TEMCD8_vs_TEMRACD8.csv"))
  
  DefaultAssay(TARA_cd8) <- "ADT"
  dpe_eff <- FindMarkers(TARA_cd8, ident.1 = "TEM CD8", ident.2 = "TEMRA CD8",
                         test.use = "MAST", logfc.threshold = 0.1,
                         min.pct = 0.05, verbose = FALSE)
  write.csv(dpe_eff, file.path(dpe_dir, "DPE_MAST_ADT_TEMCD8_vs_TEMRACD8.csv"))
}, error = function(e) cat("    ERROR:", e$message, "\n"))

cat("  Pairwise DGE/DPE saved.\n")

################################################################################
# 8B. DGE/DPE: Expanding clones — Suppressed vs Unsuppressed
################################################################################
cat("\n=== Running DGE/DPE: Expanding clones Suppressed vs Unsuppressed ===\n")

run_expand_dge <- function(obj, cluster_name, label, dge_dir, dpe_dir) {
  if (!is.null(cluster_name)) {
    cells <- subset(obj,
                    CD8_Annotation == cluster_name &
                      !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)" &
                      Timepoint_Group %in% c("PostART_Suppressed", "PostART_Unsuppressed"))
  } else {
    cells <- subset(obj,
                    !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)" &
                      Timepoint_Group %in% c("PostART_Suppressed", "PostART_Unsuppressed"))
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

cat("  All expanding clones...\n")
run_expand_dge(TARA_cd8, NULL, "AllClusters", dge_dir, dpe_dir)

cat("  Within TEM CD8...\n")
run_expand_dge(TARA_cd8, "TEM CD8", "TEMCD8", dge_dir, dpe_dir)

cat("  Within TEMRA CD8...\n")
run_expand_dge(TARA_cd8, "TEMRA CD8", "TEMRACD8", dge_dir, dpe_dir)

# CHANGED: CD16+ Effector CD8 replaces Transitional Tem CD8
cat("  Within CD16+ Effector CD8...\n")
run_expand_dge(TARA_cd8, "CD16+ Effector CD8", "CD16EffectorCD8", dge_dir, dpe_dir)

cat("  Within KIR+ innate-like CD8...\n")
run_expand_dge(TARA_cd8, "KIR+ innate-like CD8", "KIRinnateLikeCD8", dge_dir, dpe_dir)

cat("  Expanding clone DGE/DPE saved.\n")

################################################################################
# 8C. CLUSTER PROPORTIONS
################################################################################
cat("\n=== Computing cluster proportions ===\n")

prop_tp <- as.data.frame(prop.table(
  table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 2))
colnames(prop_tp) <- c("Cluster", "Timepoint_Group", "Proportion")
write.csv(prop_tp, file.path(prop_dir, "ClusterProportions_by_TimepointGroup.csv"),
          row.names = FALSE)

p_prop_tp <- ggplot(prop_tp, aes(x = Timepoint_Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion", fill = "CD8 Sub-cluster",
       title = "CD8 Sub-cluster Composition by ART Status") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 11, angle = 35, hjust = 1),
        legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"),
        plot.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(prop_dir, "ClusterProportions_by_TimepointGroup.png"),
       plot = p_prop_tp, width = 12, height = 8, dpi = 300, bg = "white")

# By patient
pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
names(pid_vec) <- colnames(TARA_cd8)
TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")

prop_pid <- as.data.frame(prop.table(
  table(TARA_cd8$CD8_Annotation, TARA_cd8$PID), margin = 2))
colnames(prop_pid) <- c("Cluster", "PID", "Proportion")
write.csv(prop_pid, file.path(prop_dir, "ClusterProportions_by_Patient.csv"),
          row.names = FALSE)

write.csv(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group),
          file.path(prop_dir, "ClusterCounts_by_TimepointGroup.csv"))
write.csv(table(TARA_cd8$CD8_Annotation, TARA_cd8$PID),
          file.path(prop_dir, "ClusterCounts_by_Patient.csv"))

cat("  Proportion CSVs saved to:", prop_dir, "\n")

################################################################################
# 8D. CLONAL EXPANSION
################################################################################
cat("\n=== Plotting clonal expansion ===\n")

clone_size_levels <- c("Single (0 < X <= 1)", "Small (1 < X <= 5)",
                       "Medium (5 < X <= 20)", "Large (20 < X <= 100)",
                       "Hyperexpanded (100 < X <= 500)")
clone_size_colors <- setNames(c(colorblind_vector[c(1, 3, 4, 5, 7)]), clone_size_levels)

p_clone_umap <- DimPlot_scCustom(
  TARA_cd8, group.by = "cloneSize", reduction = "wnn.umap", pt.size = 0.5) +
  scale_color_manual(values = clone_size_colors) +
  ggtitle("Clonal Expansion — αβ CD8 Sub-clusters") +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(clone_dir, "CD8_UMAP_cloneSize.png"),
       plot = p_clone_umap, width = 10, height = 8, dpi = 300, bg = "white")

# CloneSize per cluster
clone_meta <- TARA_cd8@meta.data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(CD8_Annotation, cloneSize) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CD8_Annotation) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

clone_meta$cloneSize <- factor(clone_meta$cloneSize, levels = clone_size_levels)
clone_meta$CD8_Annotation <- factor(clone_meta$CD8_Annotation, levels = col_order_cd8)

p_clone_bar <- ggplot(clone_meta, aes(x = CD8_Annotation, y = n, fill = cloneSize)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.75, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = clone_size_colors, name = "Clone Size", drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of Cells") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 10, angle = 40, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(clone_dir, "CD8_CloneSize_per_cluster.png"),
       plot = p_clone_bar, width = 14, height = 8, dpi = 300, bg = "white")

write.csv(clone_meta, file.path(clone_dir, "CD8_CloneSize_per_cluster.csv"), row.names = FALSE)
cat("  Clonal expansion plots saved to:", clone_dir, "\n")

################################################################################
# 8E. TRAJECTORY ANALYSIS: Slingshot + Monocle3
#     All clusters are αβ CD8 — no exclusion needed.
################################################################################
cat("\n=== Running trajectory analysis ===\n")

if (!requireNamespace("slingshot", quietly = TRUE)) BiocManager::install("slingshot")
library(slingshot)

traj_dir <- dir_07_trajectory

# All clusters are αβ CD8 now
traj_cells <- TARA_cd8
cat("  Trajectory cells:", ncol(traj_cells), "(all clusters)\n")

# Slingshot
umap_coords <- Embeddings(traj_cells, reduction = "wnn.umap")
cluster_labels <- traj_cells$CD8_Annotation

cat("  Running Slingshot (start = Naïve CD8)...\n")
sling <- slingshot(
  data = umap_coords,
  clusterLabels = cluster_labels,
  start.clus = "Naïve CD8",
  stretch = 0
)

pseudotime_mat <- slingPseudotime(sling)
cat("  Lineages found:", ncol(pseudotime_mat), "\n")
for (i in seq_len(ncol(pseudotime_mat))) {
  cat("    Lineage", i, ":", paste(slingLineages(sling)[[i]], collapse = " → "), "\n")
}

pt_min <- apply(pseudotime_mat, 1, min, na.rm = TRUE)
pt_min[is.infinite(pt_min)] <- NA
names(pt_min) <- colnames(traj_cells)
traj_cells <- AddMetaData(traj_cells, metadata = pt_min, col.name = "pseudotime")

# Slingshot UMAP
umap_df <- data.frame(
  UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2],
  pseudotime = traj_cells$pseudotime,
  Cluster = traj_cells$CD8_Annotation,
  Timepoint = traj_cells$Timepoint_Group
)

curves_list <- slingCurves(sling)
curves_df <- do.call(rbind, lapply(seq_along(curves_list), function(i) {
  crv <- curves_list[[i]]$s[curves_list[[i]]$ord, ]
  data.frame(UMAP_1 = crv[, 1], UMAP_2 = crv[, 2], Lineage = paste0("Lineage ", i))
}))

p_pt <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.3, alpha = 0.6) +
  geom_path(data = curves_df, aes(x = UMAP_1, y = UMAP_2, group = Lineage),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_viridis_c(option = "inferno", na.value = "grey80") +
  labs(title = "Slingshot pseudotime — αβ CD8") +
  theme_cowplot(font_size = 12) + NoAxes() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(traj_dir, "Trajectory_UMAP_pseudotime.png"),
       plot = p_pt, width = 10, height = 8, dpi = 300, bg = "white")

# Pseudotime density by ART
pt_density_df <- data.frame(
  pseudotime = traj_cells$pseudotime,
  Timepoint = traj_cells$Timepoint_Group
) %>% filter(!is.na(pseudotime))

p_pt_dens <- ggplot(pt_density_df, aes(x = pseudotime, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = art_colors) +
  labs(x = "Pseudotime", y = "Density", title = "Pseudotime by ART status") +
  theme_cowplot(font_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(traj_dir, "Trajectory_pseudotime_density_byART.png"),
       plot = p_pt_dens, width = 10, height = 5, dpi = 300, bg = "white")

# Pseudotime by ART within effector clusters
pt_eff_df <- data.frame(
  pseudotime = traj_cells$pseudotime,
  Cluster = traj_cells$CD8_Annotation,
  Timepoint = traj_cells$Timepoint_Group
) %>% filter(!is.na(pseudotime), Cluster %in% effector_clusters)

p_pt_eff <- ggplot(pt_eff_df, aes(x = Timepoint, y = pseudotime, fill = Timepoint)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, alpha = 0.7) +
  stat_compare_means(
    comparisons = list(c("PreART_Entry", "PostART_Suppressed"),
                       c("PostART_Suppressed", "PostART_Unsuppressed"),
                       c("PreART_Entry", "PostART_Unsuppressed")),
    method = "wilcox.test", label = "p.signif",
    size = 3.5, step.increase = 0.12, tip.length = 0.01, hide.ns = TRUE) +
  facet_wrap(~ Cluster, nrow = 1) +
  scale_fill_manual(values = art_colors) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "Pseudotime", title = "Pseudotime by ART — effector clusters") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(traj_dir, "Trajectory_pseudotime_ART_effectorClusters.png"),
       plot = p_pt_eff, width = 12, height = 6, dpi = 300, bg = "white")

write.csv(
  data.frame(cell = colnames(traj_cells), CD8_Annotation = traj_cells$CD8_Annotation,
             Timepoint = traj_cells$Timepoint_Group, pseudotime = traj_cells$pseudotime,
             pseudotime_mat),
  file.path(traj_dir, "Slingshot_pseudotime_per_cell.csv"), row.names = FALSE)

cat("  Slingshot saved to:", traj_dir, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 8E (cont). MONOCLE3
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Running Monocle3 ===\n")

if (!requireNamespace("monocle3", quietly = TRUE)) BiocManager::install("monocle3")
library(monocle3)
library(SeuratWrappers)

monocle_dir <- file.path(dir_07_trajectory, "monocle3")
dir.create(monocle_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(traj_cells) <- "RNA"
cds <- as.cell_data_set(traj_cells)
reducedDims(cds)[["UMAP"]] <- Embeddings(traj_cells, "wnn.umap")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

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

m3_pt <- pseudotime(cds)[colnames(traj_cells)]
traj_cells <- AddMetaData(traj_cells, metadata = m3_pt, col.name = "monocle3_pseudotime")

# Monocle3 plots
p_m3_traj <- plot_cells(cds, color_cells_by = "pseudotime",
                        label_groups_by_cluster = FALSE, label_leaves = FALSE,
                        label_branch_points = TRUE, cell_size = 0.5) +
  scale_color_viridis_c(option = "inferno") +
  ggtitle("Monocle3 pseudotime") +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(monocle_dir, "Monocle3_UMAP_pseudotime.png"),
       plot = p_m3_traj, width = 10, height = 8, dpi = 300, bg = "white")

p_m3_clust <- plot_cells(cds, color_cells_by = "CD8_Annotation",
                         label_groups_by_cluster = TRUE, label_leaves = FALSE,
                         label_branch_points = FALSE, group_label_size = 3.5, cell_size = 0.5) +
  ggtitle("Monocle3 — CD8 annotations") +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(monocle_dir, "Monocle3_UMAP_clusters.png"),
       plot = p_m3_clust, width = 12, height = 8, dpi = 300, bg = "white")

colData(cds)$Timepoint_Group <- traj_cells$Timepoint_Group
p_m3_art <- plot_cells(cds, color_cells_by = "Timepoint_Group",
                       label_groups_by_cluster = FALSE, label_leaves = FALSE,
                       label_branch_points = FALSE, cell_size = 0.5) +
  scale_color_manual(values = art_colors) +
  ggtitle("Monocle3 — ART status") +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(monocle_dir, "Monocle3_UMAP_ART_status.png"),
       plot = p_m3_art, width = 10, height = 8, dpi = 300, bg = "white")

# Monocle3 density
m3_pt_df <- data.frame(pseudotime = traj_cells$monocle3_pseudotime,
                       Timepoint = traj_cells$Timepoint_Group) %>%
  filter(!is.na(pseudotime) & is.finite(pseudotime))

p_m3_dens <- ggplot(m3_pt_df, aes(x = pseudotime, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = art_colors) +
  labs(x = "Monocle3 Pseudotime", y = "Density") +
  theme_cowplot(font_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(monocle_dir, "Monocle3_pseudotime_density_byART.png"),
       plot = p_m3_dens, width = 10, height = 5, dpi = 300, bg = "white")

# Monocle3 effector clusters
m3_eff_df <- data.frame(pseudotime = traj_cells$monocle3_pseudotime,
                        Cluster = traj_cells$CD8_Annotation,
                        Timepoint = traj_cells$Timepoint_Group) %>%
  filter(!is.na(pseudotime) & is.finite(pseudotime), Cluster %in% effector_clusters)

p_m3_eff <- ggplot(m3_eff_df, aes(x = Timepoint, y = pseudotime, fill = Timepoint)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, alpha = 0.7) +
  stat_compare_means(
    comparisons = list(c("PreART_Entry", "PostART_Suppressed"),
                       c("PostART_Suppressed", "PostART_Unsuppressed"),
                       c("PreART_Entry", "PostART_Unsuppressed")),
    method = "wilcox.test", label = "p.signif",
    size = 3.5, step.increase = 0.12, tip.length = 0.01, hide.ns = TRUE) +
  facet_wrap(~ Cluster, nrow = 1) +
  scale_fill_manual(values = art_colors) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "Monocle3 Pseudotime") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(monocle_dir, "Monocle3_pseudotime_ART_effectorClusters.png"),
       plot = p_m3_eff, width = 12, height = 6, dpi = 300, bg = "white")

# graph_test
cat("  Running graph_test...\n")
tryCatch({
  graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  graph_test_sig <- graph_test_res %>% filter(q_value < 0.05) %>% arrange(q_value)
  write.csv(graph_test_sig, file.path(monocle_dir, "Monocle3_graph_test_significant.csv"))
  cat("    Significant genes:", nrow(graph_test_sig), "\n")
}, error = function(e) cat("    graph_test error:", e$message, "\n"))

# Gene expression along pseudotime helper
plot_genes_along_pt <- function(obj, genes, ncol_plot = 4, title_text = "",
                                filename, width, height) {
  DefaultAssay(obj) <- "RNA"
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  df <- data.frame(
    pseudotime = obj$monocle3_pseudotime, Timepoint = obj$Timepoint_Group,
    t(as.matrix(expr_mat)), check.names = FALSE
  ) %>%
    filter(!is.na(pseudotime) & is.finite(pseudotime)) %>%
    tidyr::pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")
  df$Gene <- factor(df$Gene, levels = genes)
  df$Timepoint <- factor(df$Timepoint, levels = names(art_colors))
  
  p <- ggplot(df, aes(x = pseudotime, y = Expression, color = Timepoint)) +
    geom_point(size = 0.1, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2, span = 0.5) +
    facet_wrap(~ Gene, ncol = ncol_plot, scales = "free_y") +
    scale_color_manual(values = art_colors) +
    labs(x = "Monocle3 Pseudotime", y = "Expression", title = title_text) +
    theme_cowplot(font_size = 11) +
    theme(strip.text = element_text(face = "italic", size = 10),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggsave(filename, plot = p, width = width, height = height, dpi = 300, bg = "white")
}

all_genes_avail <- rownames(traj_cells[["RNA"]])
filter_present <- function(genes) genes[genes %in% all_genes_avail]

narrative_genes <- filter_present(c("TOX","PDCD1","HAVCR2","TIGIT",
                                    "TCF7","BACH2","SELL","IL7R",
                                    "GNLY","GZMB","PRF1",
                                    "HSPA1A","CCL3","IFIT1"))

if (length(narrative_genes) >= 4) {
  plot_genes_along_pt(traj_cells, narrative_genes, ncol_plot = 4,
                      title_text = "Key narrative genes along pseudotime",
                      filename = file.path(monocle_dir, "Monocle3_NarrativeGenes_in_pseudotime.png"),
                      width = 16, height = ceiling(length(narrative_genes) / 4) * 3.5)
}

# Export Monocle3 pseudotime
write.csv(
  data.frame(cell = colnames(traj_cells), CD8_Annotation = traj_cells$CD8_Annotation,
             Timepoint = traj_cells$Timepoint_Group,
             monocle3_pseudotime = traj_cells$monocle3_pseudotime),
  file.path(monocle_dir, "Monocle3_pseudotime_per_cell.csv"), row.names = FALSE)

cat("  Monocle3 saved to:", monocle_dir, "\n")

# ── CHECKPOINT: Save after trajectory ────────────────────────────────────────
# Transfer pseudotime back to main object
TARA_cd8$monocle3_pseudotime <- traj_cells$monocle3_pseudotime[colnames(TARA_cd8)]
TARA_cd8$pseudotime <- traj_cells$pseudotime[colnames(TARA_cd8)]

qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_with_trajectory.qs2"))
cat("CHECKPOINT: Object with trajectory saved.\n")

################################################################################
# 8F. CLONAL EXPANSION vs EXHAUSTION CONTROL
################################################################################
cat("\n=== Clonal expansion vs exhaustion control ===\n")

expansion_dir <- dir_08_expansion

if (!"PID" %in% colnames(TARA_cd8@meta.data)) {
  pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
  names(pid_vec) <- colnames(TARA_cd8)
  TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")
}

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
  mutate(n_expanding = replace_na(n_expanding, 0),
         pct_expanding = n_expanding / total_cells * 100)

sample_expansion$Timepoint_Group <- factor(sample_expansion$Timepoint_Group,
                                           levels = names(art_colors))
sample_expansion$CD8_Annotation <- factor(sample_expansion$CD8_Annotation,
                                          levels = effector_clusters)

p_g1b <- ggplot(sample_expansion, aes(x = Timepoint_Group, y = pct_expanding, fill = Timepoint_Group)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(c("PreART_Entry", "PostART_Suppressed"),
                       c("PostART_Suppressed", "PostART_Unsuppressed"),
                       c("PreART_Entry", "PostART_Unsuppressed")),
    method = "wilcox.test", label = "p.signif",
    size = 3.5, step.increase = 0.12, tip.length = 0.01, hide.ns = TRUE) +
  facet_wrap(~ CD8_Annotation, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = art_colors) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  labs(x = NULL, y = "% clonally expanded",
       title = "Expansion level per sample by ART status") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(expansion_dir, "G1b_PctExpanding_perSample_byART.png"),
       plot = p_g1b, width = 14, height = 6, dpi = 300, bg = "white")

write.csv(sample_expansion, file.path(expansion_dir, "PerSample_expansion_counts.csv"),
          row.names = FALSE)

cat("  Expansion plots saved to:", expansion_dir, "\n")

################################################################################
# 8H. MODULE SCORES
################################################################################
cat("\n=== Computing module scores ===\n")

module_dir <- dir_10_modules
dir.create(file.path(module_dir, "volcanos"), recursive = TRUE, showWarnings = FALSE)

module_gene_lists <- list(
  Exhaustion = c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1"),
  Stemness   = c("TCF7", "SELL", "CCR7", "LEF1", "BACH2", "BCL2", "IL7R", "S1PR1", "KLF2"),
  Cytotoxicity = c("GZMB", "GNLY", "PRF1", "NKG7", "GZMA", "GZMH", "FGFBP2"),
  Stress_IFN = c("HSPA1A", "HSPA1B", "DNAJB1", "IFI27", "IFI44L"),
  TypeI_IFN_Memory = c("IFIT1", "IFIT3", "ISG15", "MX1"),
  Inflammatory_Chemokines = c("CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5"),
  Terminal_Differentiation = c("ZEB2", "PRDM1", "TBX21", "CX3CR1", "S1PR5", "ID2")
)

expand_eff <- subset(TARA_cd8,
                     CD8_Annotation %in% effector_clusters &
                       !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)")

cat("  Expanding effector cells:", ncol(expand_eff), "\n")
cat("  Breakdown:\n")
print(table(expand_eff$CD8_Annotation, expand_eff$Timepoint_Group))

DefaultAssay(expand_eff) <- "RNA"
for (mod_name in names(module_gene_lists)) {
  genes <- module_gene_lists[[mod_name]]
  genes_present <- genes[genes %in% rownames(expand_eff[["RNA"]])]
  if (length(genes_present) < 2) { cat("  Skip", mod_name, "\n"); next }
  cat("  Computing", mod_name, "(", length(genes_present), "/", length(genes), ")...\n")
  expand_eff <- AddModuleScore(expand_eff, features = list(genes_present),
                               name = mod_name, ctrl = min(100, length(genes_present) * 5))
}

# Clean column names
score_cols <- paste0(names(module_gene_lists), "1")
score_cols_clean <- names(module_gene_lists)
for (i in seq_along(score_cols)) {
  if (score_cols[i] %in% colnames(expand_eff@meta.data)) {
    colnames(expand_eff@meta.data)[colnames(expand_eff@meta.data) == score_cols[i]] <- score_cols_clean[i]
  }
}

module_labels <- c(
  "Exhaustion" = "Exhaustion", "Stemness" = "Stemness / naïve",
  "Cytotoxicity" = "Cytotoxicity", "Stress_IFN" = "Stress / IFN (acute)",
  "TypeI_IFN_Memory" = "Type I IFN memory",
  "Inflammatory_Chemokines" = "Inflammatory chemokines",
  "Terminal_Differentiation" = "Terminal differentiation"
)

available_modules <- score_cols_clean[score_cols_clean %in% colnames(expand_eff@meta.data)]

# Build plot data
plot_data <- expand_eff@meta.data %>%
  select(CD8_Annotation, Timepoint_Group, all_of(available_modules)) %>%
  tidyr::pivot_longer(cols = all_of(available_modules), names_to = "Module", values_to = "Score")
plot_data$Timepoint_Group <- factor(plot_data$Timepoint_Group, levels = names(art_colors))
plot_data$CD8_Annotation <- factor(plot_data$CD8_Annotation, levels = effector_clusters)
plot_data$Module_Label <- module_labels[plot_data$Module]
plot_data$Module_Label <- factor(plot_data$Module_Label, levels = module_labels)

# ── Cohen's d + Wilcoxon ─────────────────────────────────────────────────────
cat("  Computing Cohen's d...\n")

comparisons_list <- list(
  c("PreART_Entry", "PostART_Suppressed"),
  c("PostART_Suppressed", "PostART_Unsuppressed"),
  c("PreART_Entry", "PostART_Unsuppressed")
)

stat_results_fc <- list()
for (mod in available_modules) {
  for (cl in effector_clusters) {
    cl_data <- expand_eff@meta.data %>% filter(CD8_Annotation == cl)
    sup   <- cl_data %>% filter(Timepoint_Group == "PostART_Suppressed") %>% pull(!!sym(mod))
    unsup <- cl_data %>% filter(Timepoint_Group == "PostART_Unsuppressed") %>% pull(!!sym(mod))
    pre   <- cl_data %>% filter(Timepoint_Group == "PreART_Entry") %>% pull(!!sym(mod))
    
    pairs <- list(
      list(sup, unsup, "Suppressed vs Unsuppressed"),
      list(pre, unsup, "Pre-ART vs Unsuppressed"),
      list(sup, pre, "Suppressed vs Pre-ART")
    )
    for (pr in pairs) {
      if (length(pr[[1]]) >= 5 & length(pr[[2]]) >= 5) {
        wt <- wilcox.test(pr[[1]], pr[[2]])
        stat_results_fc[[paste(mod, cl, pr[[3]])]] <- data.frame(
          Module = mod, CD8_Annotation = cl, Comparison = pr[[3]],
          p_value = wt$p.value, stringsAsFactors = FALSE)
      }
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

write.csv(stat_fc, file.path(module_dir, "ModuleScore_Wilcoxon_all_comparisons.csv"),
          row.names = FALSE)

cat("\n=== Module score Wilcoxon results ===\n")
print(stat_fc %>% arrange(Comparison, Module, CD8_Annotation), n = 100)

# ── Cohen's d bar plot helper ────────────────────────────────────────────────
build_fc <- function(plot_data_raw, num_tp, denom_tp, comp_label, stat_fc_df) {
  fc <- plot_data_raw %>%
    filter(Timepoint_Group %in% c(num_tp, denom_tp)) %>%
    group_by(Module, Module_Label, CD8_Annotation) %>%
    summarise(
      mean_num = mean(Score[Timepoint_Group == num_tp], na.rm = TRUE),
      mean_denom = mean(Score[Timepoint_Group == denom_tp], na.rm = TRUE),
      sd_num = sd(Score[Timepoint_Group == num_tp], na.rm = TRUE),
      sd_denom = sd(Score[Timepoint_Group == denom_tp], na.rm = TRUE),
      n_num = sum(Timepoint_Group == num_tp),
      n_denom = sum(Timepoint_Group == denom_tp), .groups = "drop"
    ) %>%
    mutate(pooled_sd = sqrt(((n_num - 1) * sd_num^2 + (n_denom - 1) * sd_denom^2) / (n_num + n_denom - 2)),
           cohens_d = (mean_num - mean_denom) / pooled_sd,
           direction = ifelse(cohens_d > 0,
                              paste0("\u2191 ", gsub("PostART_", "", num_tp)),
                              paste0("\u2191 ", gsub("PostART_", "", denom_tp))))
  fc$Comparison <- comp_label
  stars <- stat_fc_df %>% filter(Comparison == comp_label) %>% select(Module, CD8_Annotation, star)
  fc <- fc %>% left_join(stars, by = c("Module", "CD8_Annotation"))
  fc$star[is.na(fc$star)] <- ""
  fc$Module_Label <- factor(fc$Module_Label, levels = rev(module_labels))
  fc$CD8_Annotation <- factor(fc$CD8_Annotation, levels = effector_clusters)
  fc
}

plot_fc_bars <- function(fc_df, title_text, fill_colors) {
  x_range <- max(abs(fc_df$cohens_d), na.rm = TRUE)
  offset <- x_range * 0.06
  fc_df$star_x <- ifelse(fc_df$cohens_d > 0, fc_df$cohens_d + offset, fc_df$cohens_d - offset)
  fc_df$star_hjust <- ifelse(fc_df$cohens_d > 0, 0, 1)
  
  ggplot(fc_df, aes(x = cohens_d, y = Module_Label, fill = direction)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(aes(x = star_x, label = star, hjust = star_hjust),
              size = 4.5, color = "black", vjust = 0.4) +
    facet_wrap(~ CD8_Annotation, nrow = 1) +
    scale_fill_manual(values = fill_colors, name = NULL) +
    coord_cartesian(clip = "off") +
    labs(x = "Cohen's d", y = NULL, title = title_text) +
    theme_cowplot(font_size = 12) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position = "bottom",
          plot.margin = margin(10, 30, 10, 10),
          plot.background = element_rect(fill = "white", color = NA))
}

fc_sup_unsup <- build_fc(plot_data, "PostART_Suppressed", "PostART_Unsuppressed",
                         "Suppressed vs Unsuppressed", stat_fc)
fc_pre_unsup <- build_fc(plot_data, "PreART_Entry", "PostART_Unsuppressed",
                         "Pre-ART vs Unsuppressed", stat_fc)
fc_sup_pre   <- build_fc(plot_data, "PostART_Suppressed", "PreART_Entry",
                         "Suppressed vs Pre-ART", stat_fc)

p_fc1 <- plot_fc_bars(fc_sup_unsup, "Suppressed vs Unsuppressed — expanding clones",
                      c("\u2191 Suppressed" = "#52B788", "\u2191 Unsuppressed" = "#E76F51"))
ggsave(file.path(module_dir, "Fig4E_DiffMean_Sup_vs_Unsup.png"),
       plot = p_fc1, width = 14, height = 7, dpi = 300, bg = "white")

p_fc2 <- plot_fc_bars(fc_pre_unsup, "Pre-ART vs Unsuppressed — expanding clones",
                      c("\u2191 PreART_Entry" = "#4A90D9", "\u2191 Unsuppressed" = "#E76F51"))
ggsave(file.path(module_dir, "Fig4E_DiffMean_PreART_vs_Unsup.png"),
       plot = p_fc2, width = 14, height = 7, dpi = 300, bg = "white")

p_fc3 <- plot_fc_bars(fc_sup_pre, "Suppressed vs Pre-ART — expanding clones",
                      c("\u2191 Suppressed" = "#52B788", "\u2191 PreART_Entry" = "#4A90D9"))
ggsave(file.path(module_dir, "Fig4E_DiffMean_Sup_vs_PreART.png"),
       plot = p_fc3, width = 14, height = 7, dpi = 300, bg = "white")

cat("  Cohen's d bar plots saved.\n")

# ── Volcano plots ────────────────────────────────────────────────────────────
cat("  Drawing volcano plots...\n")

volcano_dir <- file.path(module_dir, "volcanos")

exhaustion_genes   <- c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1")
stemness_genes     <- c("TCF7", "SELL", "CCR7", "LEF1", "BACH2", "BCL2", "IL7R", "S1PR1", "KLF2")
cytotoxic_genes    <- c("GZMB", "GNLY", "PRF1", "NKG7", "GZMA", "GZMH", "FGFBP2")
stress_genes       <- c("HSPA1A", "HSPA1B", "DNAJB1", "IFI27", "IFI44L")
ifn_memory_genes   <- c("IFIT1", "IFIT3", "ISG15", "MX1")
chemokine_genes    <- c("CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5")
terminal_genes     <- c("ZEB2", "PRDM1", "TBX21", "CX3CR1", "S1PR5", "ID2", "EOMES")
activation_genes   <- c("HLA-DRA", "CD38", "FAS", "CD69")

highlight_cols <- c(
  "Exhaustion" = "#A32D2D", "Stemness" = "#185FA5", "Cytotoxicity" = "#3B6D11",
  "Stress response" = "#BA7517", "Type I IFN memory" = "#0F6E56",
  "Chemokines" = "#993556", "Terminal diff." = "#534AB7",
  "Activation" = "#D85A30", "Other" = "grey80"
)

# UPDATED: DGE file names for new cluster names
dge_files <- list(
  "TEM CD8"              = file.path(dge_dir, "DGE_MAST_Expanding_TEMCD8_Sup_vs_Unsup.csv"),
  "TEMRA CD8"            = file.path(dge_dir, "DGE_MAST_Expanding_TEMRACD8_Sup_vs_Unsup.csv"),
  "CD16+ Effector CD8"   = file.path(dge_dir, "DGE_MAST_Expanding_CD16EffectorCD8_Sup_vs_Unsup.csv"),
  "KIR+ innate-like CD8" = file.path(dge_dir, "DGE_MAST_Expanding_KIRinnateLikeCD8_Sup_vs_Unsup.csv")
)

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

for (cl_name in names(dge_files)) {
  fpath <- dge_files[[cl_name]]
  if (!file.exists(fpath)) { cat("    Skipping", cl_name, "— DGE not found\n"); next }
  
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
  
  all_highlight <- c(exhaustion_genes, stemness_genes, cytotoxic_genes, stress_genes,
                     ifn_memory_genes, chemokine_genes, terminal_genes, activation_genes)
  
  dge_data$gene <- rownames(dge_data)
  dge_data$neg_log10_padj <- -log10(dge_data$p_val_adj + 1e-300)
  dge_data$label <- ""
  dge_data$label[dge_data$gene %in% all_highlight & dge_data$p_val_adj < 0.05] <-
    dge_data$gene[dge_data$gene %in% all_highlight & dge_data$p_val_adj < 0.05]
  dge_data$highlight <- factor(dge_data$highlight,
                               levels = names(highlight_cols))
  
  p_v <- ggplot(dge_data, aes(x = avg_log2FC, y = neg_log10_padj, color = highlight)) +
    geom_point(data = dge_data %>% filter(highlight == "Other"), size = 0.5, alpha = 0.3) +
    geom_point(data = dge_data %>% filter(highlight != "Other"), size = 2.5, alpha = 0.85) +
    ggrepel::geom_text_repel(
      data = dge_data %>% filter(label != ""), aes(label = label),
      size = 3.5, fontface = "italic", max.overlaps = 40,
      segment.size = 0.3, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = highlight_cols, drop = FALSE) +
    labs(x = expression(log[2]~FC~(Sup/Unsup)), y = expression(-log[10]~p[adj]),
         title = paste0(cl_name, " — Sup vs Unsup")) +
    theme_cowplot(font_size = 12) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  safe_cl <- gsub("[^A-Za-z0-9]", "", cl_name)
  ggsave(file.path(volcano_dir, paste0("Volcano_", safe_cl, "_Sup_vs_Unsup.png")),
         plot = p_v, width = 12, height = 10, dpi = 300, bg = "white")
  cat("    Volcano saved:", cl_name, "\n")
}

# Export module scores
write.csv(
  expand_eff@meta.data %>%
    select(CD8_Annotation, Timepoint_Group, cloneSize, all_of(available_modules)),
  file.path(module_dir, "ModuleScores_per_cell.csv"), row.names = TRUE)

cat("  Module scores saved to:", module_dir, "\n")

################################################################################
# 8H (cont). DEFERRED MODULE-DEPENDENT ANALYSES
################################################################################
cat("\n=== Module-dependent analyses ===\n")

# ── Exhaustion by clone size × ART ──────────────────────────────────────────
exh_clone_df <- expand_eff@meta.data %>%
  filter(!is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") %>%
  mutate(Timepoint_Group = factor(Timepoint_Group, levels = names(art_colors)),
         cloneSize = factor(cloneSize, levels = c(
           "Small (1 < X <= 5)", "Medium (5 < X <= 20)",
           "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 500)")))

if ("Exhaustion" %in% available_modules) {
  p_g3 <- ggplot(exh_clone_df, aes(x = cloneSize, y = Exhaustion, fill = Timepoint_Group)) +
    geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.7, position = position_dodge(0.75)) +
    facet_wrap(~ CD8_Annotation, nrow = 1) +
    scale_fill_manual(values = art_colors) +
    scale_x_discrete(labels = c("Small", "Medium", "Large", "Hyper")) +
    labs(x = "Clone Size", y = "Exhaustion Score",
         title = "Exhaustion by clone size and ART status") +
    theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
          strip.background = element_rect(fill = "grey95", color = NA),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(dir_08_expansion, "G3_Exhaustion_by_CloneSize_and_ART.png"),
         plot = p_g3, width = 14, height = 6, dpi = 300, bg = "white")
}

# ── Viral load correlations ──────────────────────────────────────────────────
cat("  Viral load correlations...\n")
vl_dir <- dir_09_viralload

sample_scores <- expand_eff@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters) %>%
  group_by(orig.ident, PID, Timepoint_Group, Viral_Load_num) %>%
  summarise(n_expanding = n(),
            across(all_of(available_modules), mean, na.rm = TRUE),
            .groups = "drop")
sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)

for (mod in c("Exhaustion", "Stemness")) {
  if (!mod %in% available_modules) next
  ct <- cor.test(sample_scores$log10_VL, sample_scores[[mod]], method = "spearman", exact = FALSE)
  
  p <- ggplot(sample_scores, aes(x = log10_VL, y = .data[[mod]])) +
    geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = PID), size = 2.5, max.overlaps = 15, color = "grey40") +
    scale_color_manual(values = art_colors) +
    scale_size_continuous(range = c(2, 7)) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.2f, p = %.2g", ct$estimate, ct$p.value),
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(x = expression(log[10]~VL), y = paste0("Mean ", tolower(mod), " score"),
         title = paste0(mod, " vs viral load")) +
    theme_cowplot(font_size = 12) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(vl_dir, paste0("H_", mod, "_vs_ViralLoad.png")),
         plot = p, width = 10, height = 8, dpi = 300, bg = "white")
}

cor_results <- lapply(available_modules, function(mod) {
  ct <- cor.test(sample_scores$log10_VL, sample_scores[[mod]], method = "spearman", exact = FALSE)
  data.frame(Module = mod, rho = ct$estimate, p_value = ct$p.value, stringsAsFactors = FALSE)
})
cor_df <- do.call(rbind, cor_results)
cor_df$p_adj <- p.adjust(cor_df$p_value, method = "BH")
write.csv(cor_df, file.path(vl_dir, "ViralLoad_Spearman_correlations.csv"), row.names = FALSE)
write.csv(sample_scores, file.path(vl_dir, "PerSample_ModuleScores_and_ViralLoad.csv"),
          row.names = FALSE)

cat("  Viral load saved to:", vl_dir, "\n")

# ── Module scores along Monocle3 pseudotime ──────────────────────────────────
if (exists("traj_cells")) {
  cat("  Module scores along pseudotime...\n")
  common_cells <- intersect(colnames(traj_cells), colnames(expand_eff))
  
  if (length(common_cells) > 100) {
    m3_df <- data.frame(
      pseudotime = traj_cells$monocle3_pseudotime[match(common_cells, colnames(traj_cells))],
      Timepoint  = expand_eff@meta.data[common_cells, "Timepoint_Group"],
      Cluster    = expand_eff@meta.data[common_cells, "CD8_Annotation"]
    )
    for (mod in available_modules) m3_df[[mod]] <- expand_eff@meta.data[common_cells, mod]
    m3_df <- m3_df %>% filter(!is.na(pseudotime) & is.finite(pseudotime))
    
    m3_long <- m3_df %>%
      select(pseudotime, Timepoint, Exhaustion, Stemness) %>%
      tidyr::pivot_longer(cols = c(Exhaustion, Stemness), names_to = "Module", values_to = "Score")
    
    p_m3_es <- ggplot(m3_long, aes(x = pseudotime, y = Score, color = Timepoint)) +
      geom_point(size = 0.2, alpha = 0.15) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.4, alpha = 0.2, span = 0.4) +
      facet_wrap(~ Module, nrow = 1, scales = "free_y") +
      scale_color_manual(values = art_colors) +
      labs(x = "Monocle3 Pseudotime", y = "Module Score",
           title = "Exhaustion + stemness along trajectory") +
      theme_cowplot(font_size = 12) +
      theme(strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = "grey95", color = NA),
            legend.position = "bottom",
            plot.background = element_rect(fill = "white", color = NA)) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(monocle_dir, "Monocle3_ExhaustionStemness_along_pseudotime.png"),
           plot = p_m3_es, width = 14, height = 6, dpi = 300, bg = "white")
    
    cat("  Pseudotime × module plots saved.\n")
  }
}

################################################################################
# 9. SAVE FINAL OBJECT
################################################################################
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("\n=== PIPELINE COMPLETE ===\n")
cat("  Final object:", ncol(TARA_cd8), "cells,",
    length(levels(TARA_cd8$CD8_Annotation)), "populations\n")
cat("  Saved to:", file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"), "\n")
cat("  Analysis:", analysis_dir, "\n")

# Print key tables for manuscript
cat("\n=== KEY TABLES ===\n")
cat("\nCluster × Timepoint:\n")
print(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group))
cat("\nWilcoxon results:\n")
print(stat_fc %>% arrange(Comparison, Module), n = 100)
cat("\nViral load correlations:\n")
print(cor_df)

sessionInfo()

