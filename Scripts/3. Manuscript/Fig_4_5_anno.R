################################################################################
# FIGURES 4–5 ANNOTATION VALIDATION
# CD8 Sub-cluster Annotation Support for Reviewer / Collaborator
#
# Loads the already-annotated HEI CD8 object (TARA_cd8_HEI_annotated.qs2)
# and generates:
#   1. VlnPlot2 (SeuratExtend) — per marker, grouped by CD8_Annotation
#      - ALL 108 ADT markers
#      - Expanded RNA markers (same as Fig 1 list: CD8 + non-CD8 lineage)
#   2. DimPlot2 feature plots — per marker on WNN UMAP
#      - ALL 108 ADT
#      - Expanded RNA
#   3. Proportion analyses:
#      - Cluster composition (% of total CD8)
#      - Per-cluster breakdown by Timepoint_Group (Pre-ART / Sup / Unsup)
#      - Per-cluster breakdown by donor (PID)
#      - Per-donor stacked bars showing cluster composition
#      - Per-timepoint stacked bars
#   4. Average expression CSVs (RNA + ADT, raw + scaled + %expressing)
#
# Input:  TARA_cd8_HEI_annotated.qs2
# Output: Fig 4-5/analysis/
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(SeuratExtend)
library(scCustomize)
library(cowplot)
library(patchwork)
library(pheatmap)
library(grid)
library(qs2)
library(scales)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir    <- "~/Documents/CD8_Longitudinal"
saved_dir   <- file.path(base_dir, "saved_R_data")
out_dir     <- file.path(base_dir, "Manuscript", "Fig 4-5")
analysis_dir <- file.path(out_dir, "analysis")

# Subfolders
out_vln_rna   <- file.path(analysis_dir, "Violins_RNA")
out_vln_adt   <- file.path(analysis_dir, "Violins_ADT")
out_feat_rna  <- file.path(analysis_dir, "FeaturePlots_RNA")
out_feat_adt  <- file.path(analysis_dir, "FeaturePlots_ADT")
out_props     <- file.path(analysis_dir, "Proportions")
out_avgexpr   <- file.path(analysis_dir, "AvgExpression")
out_dotplots  <- file.path(analysis_dir, "DotPlots")

for (d in c(out_vln_rna, out_vln_adt, out_feat_rna, out_feat_adt,
            out_props, out_avgexpr, out_dotplots)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Load annotated object ────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated.qs2"))

cat("Loaded annotated HEI CD8 object:", ncol(TARA_cd8), "cells\n")
cat("Annotations:\n")
print(table(TARA_cd8$CD8_Annotation))
cat("\n")

# ── Ensure PID exists ─────────────────────────────────────────────────────────
if (!"PID" %in% colnames(TARA_cd8@meta.data)) {
  TARA_cd8$PID <- sub("_.*$", "", TARA_cd8$orig.ident)
}

# ── Canonical cluster order (biological progression) ─────────────────────────
col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 (innate-like)", "Naïve CD8 (post-ART)",
  "Transitional Tem CD8", "Effector CD8", "Terminal TEMRA CD8",
  "KIR+ innate-like CD8", "MAIT-like Trm",
  "TRDV1+ γδ T cell", "Naïve-like TRDV1+ γδ", "Vγ9Vδ2 γδ T cell"
)

# Set factor levels for consistent ordering
TARA_cd8$CD8_Annotation <- factor(TARA_cd8$CD8_Annotation, levels = col_order_cd8)
Idents(TARA_cd8) <- "CD8_Annotation"

# ART status colors
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

cat("Setup complete. Starting analyses...\n\n")


################################################################################
# 1. MARKER LISTS
################################################################################

# ── ADT: ALL markers in the assay ────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"
all_adt_markers <- sort(rownames(TARA_cd8))
cat("Total ADT markers:", length(all_adt_markers), "\n")

# ── RNA: Expanded list (CD8 + non-CD8 lineage markers) ──────────────────────
rna_markers <- list(
  `CD8 Identity` = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  `Naive - Quiescence` = c("CCR7", "SELL", "TCF7", "LEF1", "BACH2", "IL7R", "BCL2",
                           "S1PR1", "KLF2", "FOXP1", "SATB1"),
  `Effector - Cytotoxicity` = c("GZMB", "GZMA", "GZMH", "GZMK", "GZMM",
                                "PRF1", "GNLY", "NKG7", "KLRG1", "CX3CR1", "FGFBP2", "FCGR3A"),
  `Effector TFs` = c("TBX21", "EOMES", "RUNX3", "PRDM1", "ZEB2", "ID2"),
  `Exhaustion` = c("TOX", "TOX2", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
                   "CTLA4", "ENTPD1", "CD160", "CD244"),
  `Memory - Stemness` = c("IL7R", "CD27", "CD28", "BACH2", "TCF7", "BCL2",
                          "BCL6", "CXCR5", "CCR7"),
  `Activation - Proliferation` = c("MKI67", "HLA-DRA", "HLA-DRB1", "CD38", "FAS",
                                   "ICOS", "TNFRSF9", "CXCR6", "CD69"),
  `Tissue Residency` = c("ITGAE", "ITGA1", "CD69", "CXCR6", "ZNF683", "PRDM1"),
  `GammaDelta TCR` = c("TRDV1", "TRDV2", "TRGV9", "TRDC", "TRGC1", "TRGC2"),
  `MAIT` = c("TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "RORC"),
  `NK-like - Innate` = c("TYROBP", "KLRD1", "KIR2DL3", "KIR3DL1", "NCAM1",
                         "KLRF1", "KLRC1", "KLRC2"),
  `Type I IFN` = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44L", "MX1",
                   "OAS1", "STAT1", "IRF7"),
  `Stress - Heat Shock` = c("HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "IFI27"),
  `Chemokines - Cytokines` = c("CCL3", "CCL4", "CCL4L2", "CCL5", "XCL1", "XCL2",
                               "IFNG", "TNF", "IL2", "CSF2"),
  `Terminal Differentiation` = c("B3GAT1", "KLRG1", "CX3CR1", "ZEB2", "GZMB"),
  # Non-CD8 lineage (contamination check)
  `CD4 T cell` = c("CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "GATA3", "RORC", "MAF"),
  `B cell` = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5"),
  `NK cell` = c("NCAM1", "KLRF1", "NCR1", "NCR3", "FCGR3A", "TYROBP", "SPON2"),
  `Monocyte` = c("CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1"),
  `pDC` = c("LILRA4", "CLEC4C", "IRF7", "TCF4", "IL3RA"),
  `Platelet` = c("PPBP", "PF4", "GP9"),
  `Erythrocyte` = c("HBB", "HBA1", "HBA2")
)

rna_markers_flat <- unique(unlist(rna_markers))

DefaultAssay(TARA_cd8) <- "RNA"
rna_available <- rna_markers_flat[rna_markers_flat %in% rownames(TARA_cd8)]
rna_missing   <- setdiff(rna_markers_flat, rna_available)

cat("RNA markers:", length(rna_available), "available,", length(rna_missing), "missing\n")
if (length(rna_missing) > 0) cat("Missing:", paste(rna_missing, collapse = ", "), "\n")
cat("\n")


################################################################################
# 2. AVERAGE EXPRESSION TABLES
################################################################################

message("Exporting average expression tables...")

# ── RNA ──────────────────────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "RNA"

avg_rna_raw <- AverageExpression(TARA_cd8, assays = "RNA", features = rna_available,
                                 group.by = "CD8_Annotation", slot = "data")$RNA

avg_rna_full <- AverageExpression(TARA_cd8, assays = "RNA",
                                  group.by = "CD8_Annotation", slot = "data")$RNA

scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x)) / rng
  }))
}

avg_rna_scaled <- scale_01(as.matrix(avg_rna_raw))
colnames(avg_rna_raw)    <- gsub("^g\\s*", "", colnames(avg_rna_raw))
colnames(avg_rna_scaled) <- gsub("^g\\s*", "", colnames(avg_rna_scaled))
colnames(avg_rna_full)   <- gsub("^g\\s*", "", colnames(avg_rna_full))

write.csv(avg_rna_raw, file.path(out_avgexpr, "AvgExpr_RNA_per_cluster.csv"))
write.csv(round(avg_rna_scaled, 4), file.path(out_avgexpr, "AvgExpr_RNA_per_cluster_scaled.csv"))
write.csv(avg_rna_full, file.path(out_avgexpr, "AvgExpr_RNA_AllGenes_per_cluster.csv"))

# ── ADT ──────────────────────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"

avg_adt_raw <- AverageExpression(TARA_cd8, assays = "ADT", features = all_adt_markers,
                                 group.by = "CD8_Annotation", slot = "data")$ADT
avg_adt_scaled <- scale_01(as.matrix(log10(avg_adt_raw + 1)))
colnames(avg_adt_raw)    <- gsub("^g\\s*", "", colnames(avg_adt_raw))
colnames(avg_adt_scaled) <- gsub("^g\\s*", "", colnames(avg_adt_scaled))

write.csv(avg_adt_raw, file.path(out_avgexpr, "AvgExpr_ADT_per_cluster.csv"))
write.csv(round(avg_adt_scaled, 4), file.path(out_avgexpr, "AvgExpr_ADT_per_cluster_scaled.csv"))

# ── % expressing ─────────────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "RNA"
clusters <- levels(droplevels(TARA_cd8$CD8_Annotation))

pct_rna <- as.data.frame(t(sapply(clusters, function(cl) {
  cells <- WhichCells(TARA_cd8, expression = CD8_Annotation == cl)
  mat   <- GetAssayData(TARA_cd8, slot = "data")[rna_available, cells, drop = FALSE]
  rowMeans(mat > 0) * 100
})))
write.csv(round(pct_rna, 2), file.path(out_avgexpr, "PctExpr_RNA_per_cluster.csv"))

DefaultAssay(TARA_cd8) <- "ADT"
pct_adt <- as.data.frame(t(sapply(clusters, function(cl) {
  cells <- WhichCells(TARA_cd8, expression = CD8_Annotation == cl)
  mat   <- GetAssayData(TARA_cd8, slot = "data")[all_adt_markers, cells, drop = FALSE]
  rowMeans(mat > 0) * 100
})))
write.csv(round(pct_adt, 2), file.path(out_avgexpr, "PctExpr_ADT_per_cluster.csv"))

message("✓ Average expression tables saved\n")


################################################################################
# 3. VIOLIN PLOTS — ADT (all 108 markers)
################################################################################

message("Generating ADT violin plots (", length(all_adt_markers), " markers)...")

for (feat in all_adt_markers) {
  
  p <- VlnPlot2(
    TARA_cd8,
    features     = feat,
    group.by     = "CD8_Annotation",
    assay        = "ADT",
    cols         = "light",
    show.mean    = TRUE,
    mean_colors  = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(
    file.path(out_vln_adt, paste0("VlnADT_", safe_name, ".png")),
    p, width = 14, height = 5, dpi = 300, bg = "white"
  )
}

message("✓ All ADT violin plots saved\n")


################################################################################
# 4. VIOLIN PLOTS — RNA (expanded markers)
################################################################################

message("Generating RNA violin plots (", length(rna_available), " markers)...")

for (gene in rna_available) {
  
  if (!gene %in% rownames(TARA_cd8[["RNA"]])) next
  
  p <- VlnPlot2(
    TARA_cd8,
    features     = gene,
    group.by     = "CD8_Annotation",
    assay        = "RNA",
    cols         = "light",
    show.mean    = TRUE,
    mean_colors  = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(
    file.path(out_vln_rna, paste0("VlnRNA_", safe_name, ".png")),
    p, width = 14, height = 5, dpi = 300, bg = "white"
  )
}

message("✓ All RNA violin plots saved\n")


################################################################################
# 5. FEATURE PLOTS — ADT (all 108, DimPlot2 on WNN UMAP)
################################################################################

message("Generating ADT feature plots (", length(all_adt_markers), " markers)...")

DefaultAssay(TARA_cd8) <- "ADT"

for (feat in all_adt_markers) {
  
  if (!feat %in% rownames(TARA_cd8[["ADT"]])) next
  
  fp <- DimPlot2(
    TARA_cd8,
    features  = feat,
    reduction = "wnn.umap"
  ) + ggtitle(paste("ADT |", feat))
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(
    file.path(out_feat_adt, paste0("FeatADT_", safe_name, ".png")),
    fp, dpi = 500, width = 8, height = 6, bg = "white"
  )
}

message("✓ All ADT feature plots saved\n")


################################################################################
# 6. FEATURE PLOTS — RNA (expanded markers, DimPlot2 on WNN UMAP)
################################################################################

message("Generating RNA feature plots (", length(rna_available), " markers)...")

DefaultAssay(TARA_cd8) <- "RNA"

for (gene in rna_available) {
  
  if (!gene %in% rownames(TARA_cd8[["RNA"]])) next
  
  fp <- DimPlot2(
    TARA_cd8,
    features  = gene,
    reduction = "wnn.umap"
  ) + ggtitle(paste("RNA |", gene))
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(
    file.path(out_feat_rna, paste0("FeatRNA_", safe_name, ".png")),
    fp, dpi = 500, width = 8, height = 6, bg = "white"
  )
}

message("✓ All RNA feature plots saved\n")


################################################################################
# 7. DOT PLOTS (overview)
################################################################################

message("Generating DotPlots...")

# RNA DotPlot — curated key markers
rna_dotplot_genes <- intersect(
  c("CD8A", "CD8B", "CD3D", "CD3E",
    "CCR7", "SELL", "TCF7", "LEF1", "BACH2", "S1PR1", "KLF2",
    "IL7R", "BCL2",
    "GZMB", "GZMA", "GZMK", "GZMH", "PRF1", "GNLY", "NKG7", "FGFBP2",
    "TBX21", "EOMES", "RUNX3", "ZEB2",
    "TOX", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
    "B3GAT1", "CX3CR1", "KLRG1",
    "TRDV1", "TRGV9", "TRDC",
    "SLC4A10", "KLRB1", "ZBTB16",
    "TYROBP", "KLRD1", "KIR2DL3", "KIR3DL1",
    "MKI67", "HLA-DRA", "CD69",
    "ISG15", "IFIT1", "IFI44L",
    "CCL3", "CCL4", "CCL5", "IFNG",
    "HSPA1A", "HSPA1B",
    "CD4", "FOXP3", "CD19", "MS4A1", "CD14", "LYZ", "PPBP", "HBB"),
  rna_available
)

DefaultAssay(TARA_cd8) <- "RNA"
p_dot_rna <- DotPlot(TARA_cd8, features = rna_dotplot_genes, group.by = "CD8_Annotation",
                     cols = c("lightgrey", "#D73027"), dot.scale = 6) +
  RotatedAxis() +
  ggtitle("RNA: CD8 Sub-cluster Markers") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 8))

ggsave(file.path(out_dotplots, "DotPlot_RNA.png"),
       p_dot_rna, width = 26, height = 8, dpi = 300, bg = "white")

# ADT DotPlot — all 108
DefaultAssay(TARA_cd8) <- "ADT"
p_dot_adt <- DotPlot(TARA_cd8, features = all_adt_markers, group.by = "CD8_Annotation",
                     cols = c("lightgrey", "#08519C"), dot.scale = 5) +
  RotatedAxis() +
  ggtitle("ADT: All Surface Protein Markers") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 8))

ggsave(file.path(out_dotplots, "DotPlot_ADT_All.png"),
       p_dot_adt, width = 38, height = 8, dpi = 300, bg = "white")

message("✓ DotPlots saved\n")


################################################################################
# 8. PROPORTIONS — comprehensive breakdown
################################################################################

message("Generating proportion analyses...")

meta <- TARA_cd8@meta.data %>%
  mutate(
    CD8_Annotation = factor(CD8_Annotation, levels = col_order_cd8),
    Timepoint_Group = factor(Timepoint_Group,
                             levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
  )

# ── 8a. Overall cluster composition (% of total CD8) ─────────────────────────
cluster_counts <- meta %>%
  group_by(CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = n / sum(n) * 100)

write.csv(cluster_counts, file.path(out_props, "Cluster_overall_composition.csv"), row.names = FALSE)

p_overall <- ggplot(cluster_counts, aes(x = CD8_Annotation, y = pct, fill = CD8_Annotation)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")),
            vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% of total CD8 T cells",
       title = "CD8 Sub-cluster Composition (HEI)") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 40, hjust = 1),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out_props, "Prop_overall_cluster_composition.png"),
       p_overall, width = 14, height = 7, dpi = 300, bg = "white")

# ── 8b. Per-cluster: Timepoint_Group breakdown ──────────────────────────────
# What fraction of each cluster comes from each timepoint?
tp_per_cluster <- meta %>%
  group_by(CD8_Annotation, Timepoint_Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CD8_Annotation) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

write.csv(tp_per_cluster, file.path(out_props, "TimepointGroup_per_cluster.csv"), row.names = FALSE)

p_tp_per_cl <- ggplot(tp_per_cluster,
                      aes(x = CD8_Annotation, y = pct, fill = Timepoint_Group)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = art_colors, name = "Timepoint Group",
                    labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of cluster",
       title = "Timepoint Group Composition Within Each CD8 Sub-cluster") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 40, hjust = 1),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out_props, "Prop_TimepointGroup_per_cluster.png"),
       p_tp_per_cl, width = 16, height = 7, dpi = 300, bg = "white")

# ── 8c. Per-cluster: Donor (PID) breakdown ──────────────────────────────────
# What fraction of each cluster comes from each donor?
pid_per_cluster <- meta %>%
  group_by(CD8_Annotation, PID) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CD8_Annotation) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

write.csv(pid_per_cluster, file.path(out_props, "Donor_per_cluster.csv"), row.names = FALSE)

p_pid_per_cl <- ggplot(pid_per_cluster,
                       aes(x = CD8_Annotation, y = pct, fill = PID)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of cluster",
       title = "Donor Composition Within Each CD8 Sub-cluster") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 9, angle = 40, hjust = 1),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out_props, "Prop_Donor_per_cluster.png"),
       p_pid_per_cl, width = 16, height = 7, dpi = 300, bg = "white")

# ── 8d. Per-timepoint: cluster composition ───────────────────────────────────
# What does the CD8 compartment look like at each timepoint?
cl_per_tp <- meta %>%
  group_by(Timepoint_Group, CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Timepoint_Group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

write.csv(cl_per_tp, file.path(out_props, "Cluster_per_TimepointGroup.csv"), row.names = FALSE)

p_cl_per_tp <- ggplot(cl_per_tp,
                      aes(x = Timepoint_Group, y = pct, fill = CD8_Annotation)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of CD8 T cells", fill = "CD8 Sub-cluster",
       title = "CD8 Sub-cluster Composition by ART Status") +
  theme_cowplot(font_size = 12) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"),
        plot.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(out_props, "Prop_Cluster_per_TimepointGroup.png"),
       p_cl_per_tp, width = 12, height = 8, dpi = 300, bg = "white")

# ── 8e. Per-donor: cluster composition ───────────────────────────────────────
cl_per_pid <- meta %>%
  group_by(PID, CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PID) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

write.csv(cl_per_pid, file.path(out_props, "Cluster_per_Donor.csv"), row.names = FALSE)

p_cl_per_pid <- ggplot(cl_per_pid,
                       aes(x = PID, y = pct, fill = CD8_Annotation)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of CD8 T cells", fill = "CD8 Sub-cluster",
       title = "CD8 Sub-cluster Composition by Donor") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"),
        plot.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(out_props, "Prop_Cluster_per_Donor.png"),
       p_cl_per_pid, width = 14, height = 8, dpi = 300, bg = "white")

# ── 8f. Per-donor × per-timepoint: cluster composition ──────────────────────
cl_per_pid_tp <- meta %>%
  group_by(PID, Timepoint_Group, CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PID, Timepoint_Group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

write.csv(cl_per_pid_tp, file.path(out_props, "Cluster_per_Donor_per_Timepoint.csv"), row.names = FALSE)

# Create a combined donor_timepoint label
cl_per_pid_tp$Donor_TP <- paste0(cl_per_pid_tp$PID, "\n",
                                 c("PreART_Entry" = "Pre-ART",
                                   "PostART_Suppressed" = "Sup",
                                   "PostART_Unsuppressed" = "Unsup")[as.character(cl_per_pid_tp$Timepoint_Group)])

p_cl_pid_tp <- ggplot(cl_per_pid_tp,
                      aes(x = Donor_TP, y = pct, fill = CD8_Annotation)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of CD8 T cells", fill = "CD8 Sub-cluster",
       title = "CD8 Sub-cluster Composition by Donor × Timepoint") +
  theme_cowplot(font_size = 11) +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"),
        plot.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(out_props, "Prop_Cluster_per_Donor_per_Timepoint.png"),
       p_cl_pid_tp, width = 20, height = 8, dpi = 300, bg = "white")

# ── 8g. Cell counts table ────────────────────────────────────────────────────
count_table_tp <- as.data.frame.matrix(table(meta$CD8_Annotation, meta$Timepoint_Group))
count_table_pid <- as.data.frame.matrix(table(meta$CD8_Annotation, meta$PID))

write.csv(count_table_tp, file.path(out_props, "CellCounts_Cluster_x_Timepoint.csv"))
write.csv(count_table_pid, file.path(out_props, "CellCounts_Cluster_x_Donor.csv"))

# ── 8h. Per-cluster: % from each timepoint (heatmap view) ───────────────────
tp_pct_mat <- tp_per_cluster %>%
  select(CD8_Annotation, Timepoint_Group, pct) %>%
  pivot_wider(names_from = Timepoint_Group, values_from = pct, values_fill = 0) %>%
  tibble::column_to_rownames("CD8_Annotation")

png(file.path(out_props, "Heatmap_TimepointPct_per_cluster.png"),
    width = 8, height = 8, units = "in", res = 300, bg = "white")
pheatmap(
  as.matrix(tp_pct_mat),
  cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = TRUE, number_format = "%.1f",
  color = colorRampPalette(c("white", "#52B788", "#1B4332"))(50),
  border_color = "white",
  cellwidth = 60, cellheight = 22, fontsize = 10,
  main = "% of Each Cluster from Each Timepoint Group",
  labels_col = c("Pre-ART", "Suppressed", "Unsuppressed")
)
dev.off()

message("✓ All proportion analyses saved\n")


################################################################################
# 9. UMAP REFERENCE — cluster labels + split by timepoint
################################################################################

message("Generating reference UMAPs...")

# Cluster-labeled UMAP
p_umap_cl <- DimPlot2(
  TARA_cd8, reduction = "wnn.umap", group.by = "CD8_Annotation",
  label = TRUE, repel = TRUE, label.size = 3.5, pt.size = 0.4
) + ggtitle("CD8 Sub-clusters (HEI)")

ggsave(file.path(analysis_dir, "UMAP_CD8_Annotations.png"),
       p_umap_cl, width = 10, height = 8, dpi = 300, bg = "white")

# Split by Timepoint_Group
p_umap_tp <- DimPlot2(
  TARA_cd8, reduction = "wnn.umap", group.by = "CD8_Annotation",
  split.by = "Timepoint_Group",
  label = TRUE, repel = TRUE, label.size = 2.5, pt.size = 0.3, ncol = 3
) + ggtitle("CD8 Sub-clusters Split by ART Status")

ggsave(file.path(analysis_dir, "UMAP_CD8_Annotations_split_by_Timepoint.png"),
       p_umap_tp, width = 22, height = 8, dpi = 300, bg = "white")

# Colored by Timepoint_Group
p_umap_art <- DimPlot2(
  TARA_cd8, reduction = "wnn.umap", group.by = "Timepoint_Group",
  cols = art_colors, pt.size = 0.4
) + ggtitle("ART Status on CD8 UMAP")

ggsave(file.path(analysis_dir, "UMAP_CD8_by_ART_status.png"),
       p_umap_art, width = 10, height = 8, dpi = 300, bg = "white")

# Colored by donor
p_umap_pid <- DimPlot2(
  TARA_cd8, reduction = "wnn.umap", group.by = "PID", pt.size = 0.4
) + ggtitle("Donor on CD8 UMAP")

ggsave(file.path(analysis_dir, "UMAP_CD8_by_Donor.png"),
       p_umap_pid, width = 10, height = 8, dpi = 300, bg = "white")

message("✓ Reference UMAPs saved\n")


################################################################################
# SUMMARY
################################################################################

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " Fig 4-5 annotation validation saved to:\n",
        " ", out_dir, "/analysis/\n",
        "══════════════════════════════════════════════════════════════\n\n",
        " Contents:\n",
        "   Violins_RNA/         — VlnPlot2 per gene (", length(rna_available), " RNA markers)\n",
        "   Violins_ADT/         — VlnPlot2 per marker (", length(all_adt_markers), " ADT markers)\n",
        "   FeaturePlots_RNA/    — DimPlot2 per gene on WNN UMAP\n",
        "   FeaturePlots_ADT/    — DimPlot2 per marker on WNN UMAP\n",
        "   DotPlots/            — RNA + ADT overview dotplots\n",
        "   Proportions/         — All proportion breakdowns:\n",
        "     • Overall cluster composition\n",
        "     • Timepoint breakdown per cluster (← shows why 'post-ART naïve')\n",
        "     • Donor breakdown per cluster\n",
        "     • Cluster composition per timepoint\n",
        "     • Cluster composition per donor\n",
        "     • Cluster composition per donor × timepoint\n",
        "     • Heatmap of timepoint % per cluster\n",
        "     • Raw cell count tables\n",
        "   AvgExpression/       — CSVs (raw, scaled, % expressing)\n",
        "   UMAP_*.png           — Reference UMAPs (clusters, ART status, donor)\n"
)

