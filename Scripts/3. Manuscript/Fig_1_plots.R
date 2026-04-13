################################################################################
# FIGURE 1 — MANUSCRIPT PANELS (v4)
#
# LAYOUT:
#   Fig 1A — Total PBMC UMAP (broad cell type labels)
#   Fig 1B — CD3/CD8 feature plots (mRNA + protein, inferno)
#   Fig 1C — ADT + RNA heatmap (horizontal, shared group legend)
#   Fig 1D — CD8 cluster UMAP (all PBMC, non-CD8 grayed out)
#   Fig 1E — TEMRA/CTL cluster proportion (significant only)
#   Fig 1F — Tex VL correlation (significant only)
#
# SUPPLEMENTARY:
#   All 8 cluster proportions + all 8 VL correlations (PBMC + CD8 denom)
#
# REQUIRES: TARA_ALL_sorted_v4.qs2
################################################################################


# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(SeuratExtend)
library(scCustomize)
library(qs2)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(rstatix)
library(grid)
library(scales)
library(viridis)


# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
fig1_dir  <- file.path(base_dir, "Manuscript", "Fig 1")
supp1_dir <- file.path(base_dir, "Manuscript", "Supplementary 1")

for (d in c(fig1_dir, supp1_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ── Load data ─────────────────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted_v4.qs2"))


################################################################################
# STEP 0: SUBSET TO HEI / HEU / HUU
################################################################################

keep_groups <- c("PreART_Entry", "HEU", "HUU")

TARA_sub <- subset(TARA_ALL, subset = Timepoint_Group %in% keep_groups)
TARA_sub$Condition <- factor(
  case_when(
    TARA_sub$Timepoint_Group == "PreART_Entry" ~ "HEI",
    TARA_sub$Timepoint_Group == "HEU" ~ "HEU",
    TARA_sub$Timepoint_Group == "HUU" ~ "HUU"
  ),
  levels = c("HUU", "HEU", "HEI")
)

cat("=== Working subset ===\n")
print(table(TARA_sub$Condition))
cat("\n")


################################################################################
# STEP 1: DEFINE CD8 CLUSTERS & PALETTES
################################################################################

# ORDER: Naive 1, Naive 2, Naive Intermediate, Tscm, Transitional, TEMRA/CTL, Tex, gd
cd8_cluster_names <- c(
  "1a: Naive 1 CD8",
  "6: Naive 2 CD8",
  "1c: Naive Intermediate CD8",
  "1b: Tscm CD8",
  "12a: Transitional CD8",
  "8: TEMRA/CTL CD8",
  "27: Tex CD8",
  "9: γδ T cell"
)

cd8_short_labels <- setNames(
  c("Naive 1", "Naive 2", "Naive Intermediate",
    "Tscm", "Transitional", "TEMRA/CTL",
    "Tex", "γδ T cell"),
  cd8_cluster_names
)

# IMPROVED: Publication-quality distinct colors
cluster_cols <- c(
  "1a: Naive 1 CD8"            = "#2166AC",
  "6: Naive 2 CD8"             = "#67A9CF",
  "1c: Naive Intermediate CD8" = "#D6604D",
  "1b: Tscm CD8"               = "#1A9850",
  "12a: Transitional CD8"      = "#FDAE61",
  "8: TEMRA/CTL CD8"           = "#B2182B",
  "27: Tex CD8"                = "#762A83",
  "9: γδ T cell"               = "#CC79A7"
)

cond_cols <- c("HUU" = "#4E79A7", "HEU" = "#F28E2B", "HEI" = "#E15759")

# ── Subset CD8 ────────────────────────────────────────────────────────────────
TARA_cd8 <- subset(TARA_sub,
                   subset = Manual_Annotation_refined %in% cd8_cluster_names)
TARA_cd8$Manual_Annotation_refined <- droplevels(
  factor(TARA_cd8$Manual_Annotation_refined, levels = cd8_cluster_names)
)

cat("=== CD8 subset ===\n")
print(table(TARA_cd8$Manual_Annotation_refined))
cat("\n")

# ── Helpers ───────────────────────────────────────────────────────────────────
scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x)) / rng
  }))
}

rename_seurat_cols <- function(mat) {
  number_to_cluster <- setNames(
    cd8_cluster_names,
    gsub("^([0-9a-z]+):.*", "\\1", cd8_cluster_names)
  )
  col_nums <- sub(".*_([0-9a-z]+)$", "\\1", colnames(mat))
  colnames(mat) <- number_to_cluster[col_nums]
  col_order <- cd8_cluster_names[cd8_cluster_names %in% colnames(mat)]
  mat[, col_order, drop = FALSE]
}


################################################################################
# FIG 1A — TOTAL PBMC UMAP (IMPROVED FOR PUBLICATION)
################################################################################

message("Generating Fig 1A...")

broad_labels <- case_when(
  grepl("^0:|^2:|^5:|^7:|^10:|^12b:|^20:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "CD4 T cells",
  TARA_sub$Manual_Annotation_refined %in% cd8_cluster_names
  ~ "CD8 T cells",
  grepl("^16:|^21:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "γδ T cells",
  grepl("^18:|^23:|^26:|^28:|^29:|^30:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "DN T cells",
  grepl("^3:|^15:|^17:|^25:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "NK cells",
  grepl("^4:|^11:|^13:|^22:|^31:|^35:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "B cells",
  grepl("^32:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "Plasmablasts",
  grepl("^14:|^19:|^24:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "Monocytes",
  grepl("^33:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "pDC",
  grepl("^34:", as.character(TARA_sub$Manual_Annotation_refined))
  ~ "APC",
  TRUE ~ "Other"
)
TARA_sub$Broad_CellType <- factor(broad_labels)

# IMPROVED: Publication-quality distinct colors
broad_cols <- c(
  "CD4 T cells"  = "#1B9E77",
  "CD8 T cells"  = "#D95F02",
  "γδ T cells"   = "#7570B3",
  "DN T cells"   = "#E7298A",
  "NK cells"     = "#66A61E",
  "B cells"      = "#E6AB02",
  "Plasmablasts" = "#A6761D",
  "Monocytes"    = "#666666",
  "pDC"          = "#1F78B4",
  "APC"          = "#FB9A99",
  "Other"        = "#CCCCCC"
)

# IMPROVED: Build UMAP with ggplot2 for full control
umap_coords <- as.data.frame(Embeddings(TARA_sub, reduction = "wnn.umap"))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$CellType <- TARA_sub$Broad_CellType

# Calculate cluster centroids for labels
centroids <- umap_coords %>%
  group_by(CellType) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  ) %>%
  filter(CellType != "Other")

# UMAP arrow positions
x_arrow <- min(umap_coords$UMAP1, na.rm = TRUE) + 1
y_arrow <- min(umap_coords$UMAP2, na.rm = TRUE) + 1

p_total_umap_colored <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = broad_cols) +
  # Labels: white fill, colored text and border matching cluster
  geom_label_repel(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, label = CellType, color = CellType),
    fill = "white",
    size = 12,
    fontface = "bold",
    label.size = 0.8,
    label.padding = unit(0.4, "lines"),
    box.padding = 1.2,
    point.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.8,
    min.segment.length = 0,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  # UMAP arrows - bottom left, large
  annotate("segment", 
           x = x_arrow, xend = x_arrow + 3.5,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow + 1.75, y = y_arrow - 1.3,
           label = "UMAP1", size = 7, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + 3.5,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow - 1.3, y = y_arrow + 1.75,
           label = "UMAP2", size = 7, fontface = "bold", angle = 90) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

p_total_umap_black <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = broad_cols) +
  # Labels: white fill, black text and border
  geom_label_repel(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, label = CellType),
    fill = "white",
    color = "black",
    size = 12,
    fontface = "bold",
    label.size = 0.8,
    label.padding = unit(0.4, "lines"),
    box.padding = 1.2,
    point.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.8,
    min.segment.length = 0,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  # UMAP arrows - bottom left, large
  annotate("segment", 
           x = x_arrow, xend = x_arrow + 3.5,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow + 1.75, y = y_arrow - 1.3,
           label = "UMAP1", size = 7, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + 3.5,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow - 1.3, y = y_arrow + 1.75,
           label = "UMAP2", size = 7, fontface = "bold", angle = 90) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

ggsave(file.path(fig1_dir, "Fig1A_Total_UMAP_colored.png"),
       p_total_umap_colored, width = 12, height = 11, dpi = 300, bg = "white")

ggsave(file.path(fig1_dir, "Fig1A_Total_UMAP_black.png"),
       p_total_umap_black, width = 12, height = 11, dpi = 300, bg = "white")

message("✓ Fig 1A saved (both colored and black versions)\n")


################################################################################
# FIG 1B — CD3/CD8 FEATURE PLOTS (IMPROVED - LARGER TITLES + UMAP ARROWS)
################################################################################

message("Generating Fig 1B...")

# IMPROVED: Much larger title and legend text
feat_theme <- theme(
  plot.title = element_text(size = 48, face = "bold", hjust = 0.5),
  legend.key.size = unit(1.5, "cm"),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 20, face = "bold"),
  plot.margin = margin(10, 10, 30, 30)
)

# Get UMAP range for arrow positioning
umap_range <- as.data.frame(Embeddings(TARA_sub, reduction = "wnn.umap"))
x_arrow_feat <- min(umap_range[,1], na.rm = TRUE) + 2
y_arrow_feat <- min(umap_range[,2], na.rm = TRUE) + 2

# UMAP arrow annotation layer (to add to all plots)
umap_arrows <- list(
  annotate("segment", 
           x = x_arrow_feat, xend = x_arrow_feat + 3,
           y = y_arrow_feat, yend = y_arrow_feat,
           arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
           linewidth = 1.2, color = "black"),
  annotate("text", x = x_arrow_feat + 1.5, y = y_arrow_feat - 1.5,
           label = "UMAP1", size = 6, fontface = "bold"),
  annotate("segment",
           x = x_arrow_feat, xend = x_arrow_feat,
           y = y_arrow_feat, yend = y_arrow_feat + 3,
           arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
           linewidth = 1.2, color = "black"),
  annotate("text", x = x_arrow_feat - 1.5, y = y_arrow_feat + 1.5,
           label = "UMAP2", size = 6, fontface = "bold", angle = 90)
)

DefaultAssay(TARA_sub) <- "RNA"
p_cd3_rna <- FeaturePlot(TARA_sub, features = "CD3E", reduction = "wnn.umap",
                         pt.size = 0.3, order = TRUE, raster = FALSE) +
  scale_color_viridis(option = "inferno") +
  NoAxes() + ggtitle("CD3E (mRNA)") + feat_theme + umap_arrows

p_cd8_rna <- FeaturePlot(TARA_sub, features = "CD8A", reduction = "wnn.umap",
                         pt.size = 0.3, order = TRUE, raster = FALSE) +
  scale_color_viridis(option = "inferno") +
  NoAxes() + ggtitle("CD8A (mRNA)") + feat_theme + umap_arrows

DefaultAssay(TARA_sub) <- "ADT"
p_cd3_adt <- FeaturePlot(TARA_sub, features = "CD3D", reduction = "wnn.umap",
                         pt.size = 0.3, order = TRUE, raster = FALSE) +
  scale_color_viridis(option = "inferno") +
  NoAxes() + ggtitle("CD3 (protein)") + feat_theme + umap_arrows

p_cd8_adt <- FeaturePlot(TARA_sub, features = "CD8A", reduction = "wnn.umap",
                         pt.size = 0.3, order = TRUE, raster = FALSE) +
  scale_color_viridis(option = "inferno") +
  NoAxes() + ggtitle("CD8a (protein)") + feat_theme + umap_arrows

p_features <- (p_cd3_rna | p_cd8_rna) / (p_cd3_adt | p_cd8_adt)

ggsave(file.path(fig1_dir, "Fig1B_CD3_CD8_FeaturePlots.png"),
       p_features, width = 16, height = 14, dpi = 300, bg = "white")

message("✓ Fig 1B saved\n")


################################################################################
# FIG 1C — ADT + RNA HEATMAP (horizontal/stacked, shared group legend)
################################################################################

message("Generating Fig 1C heatmaps...")

# ── ADT markers (33): assay name → group ─────────────────────────────────────
adt_heatmap_markers <- c(
  "CD7"       = "Naive/Quiescence",     "SELL"      = "Naive/Quiescence",
  "CD45RA"    = "Naive/Quiescence",     "IL7R"      = "Naive/Quiescence",
  "NT5E"      = "Memory/Stemness",      "CD45RO"    = "Memory/Stemness",
  "FAS"       = "Memory/Stemness",      "ENTPD1"    = "Memory/Stemness",
  "SLAMF6"    = "Memory/Stemness",
  "CD27"      = "Co-stimulation/Homing", "CD28"     = "Co-stimulation/Homing",
  "ICOS"      = "Co-stimulation/Homing", "ITGB7"    = "Co-stimulation/Homing",
  "B3GAT1"    = "Effector/TEMRA",        "KLRG1"    = "Effector/TEMRA",
  "KIR3DL1"   = "Effector/TEMRA",        "CX3CR1"   = "Effector/TEMRA",
  "NCAM1"     = "NK-like/Innate",        "SIGLEC7"  = "NK-like/Innate",
  "KLRD1"     = "NK-like/Innate",        "KIR2DL3"  = "NK-like/Innate",
  "PDCD1"     = "Exhaustion",            "TIGIT"    = "Exhaustion",
  "LAG3"      = "Exhaustion",            "CTLA4"    = "Exhaustion",
  "CD38"      = "Activation",            "HLA-DRA"  = "Activation",
  "CD69"      = "Activation",            "ITGA4"    = "Activation",
  "TCR-AB"    = "TCR Identity",          "TCR-vA7.2" = "TCR Identity",
  "TCR-vD2"   = "TCR Identity",          "CD8A"     = "TCR Identity"
)

adt_gene_to_protein <- c(
  "CD7" = "CD7", "SELL" = "CD62L", "CD45RA" = "CD45RA", "IL7R" = "CD127",
  "NT5E" = "CD73", "CD45RO" = "CD45RO", "FAS" = "CD95", "ENTPD1" = "CD39",
  "SLAMF6" = "SLAMF6",
  "CD27" = "CD27", "CD28" = "CD28", "ICOS" = "ICOS", "ITGB7" = "Integrin b7",
  "B3GAT1" = "CD57", "KLRG1" = "KLRG1", "KIR3DL1" = "KIR3DL1", "CX3CR1" = "CX3CR1",
  "NCAM1" = "CD56", "SIGLEC7" = "Siglec-7", "KLRD1" = "CD94", "KIR2DL3" = "KIR2DL3",
  "PDCD1" = "PD-1", "TIGIT" = "TIGIT", "LAG3" = "LAG-3", "CTLA4" = "CTLA-4",
  "CD38" = "CD38", "HLA-DRA" = "HLA-DR", "CD69" = "CD69", "ITGA4" = "CD49d",
  "TCR-AB" = "TCRab", "TCR-vA7.2" = "TCR Va7.2", "TCR-vD2" = "TCR Vd2",
  "CD8A" = "CD8a"
)

# ── RNA markers (33): gene → group ───────────────────────────────────────────
rna_heatmap_markers <- c(
  "CCR7"    = "Naive/Quiescence",      "SELL"    = "Naive/Quiescence",
  "TCF7"    = "Naive/Quiescence",      "LEF1"    = "Naive/Quiescence",
  "IL7R"    = "Memory/Stemness",       "BCL2"    = "Memory/Stemness",
  "BACH2"   = "Memory/Stemness",       "S1PR1"   = "Memory/Stemness",
  "FAS"     = "Memory/Stemness",
  "CD27"    = "Co-stimulation/Homing", "CD28"    = "Co-stimulation/Homing",
  "ICOS"    = "Co-stimulation/Homing", "ITGB7"   = "Co-stimulation/Homing",
  "GZMB"    = "Effector/TEMRA",        "GNLY"    = "Effector/TEMRA",
  "PRF1"    = "Effector/TEMRA",        "NKG7"    = "Effector/TEMRA",
  "TYROBP"  = "NK-like/Innate",        "KLRD1"   = "NK-like/Innate",
  "KIR2DL3" = "NK-like/Innate",        "FCGR3A"  = "NK-like/Innate",
  "TOX"     = "Exhaustion",            "PDCD1"   = "Exhaustion",
  "TIGIT"   = "Exhaustion",            "HAVCR2"  = "Exhaustion",
  "MKI67"   = "Activation",            "HLA-DRA" = "Activation",
  "CD69"    = "Activation",            "CD38"    = "Activation",
  "TRDV1"   = "TCR Identity",          "TRGV9"   = "TCR Identity",
  "TRDC"    = "TCR Identity",          "CD8A"    = "TCR Identity"
)

# ── Shared group colors ──────────────────────────────────────────────────────
group_colors <- c(
  "Naive/Quiescence"      = "#74C2E1",
  "Memory/Stemness"       = "#66BB6A",
  "Co-stimulation/Homing" = "#AED581",
  "Effector/TEMRA"        = "#EF5350",
  "NK-like/Innate"        = "#F06292",
  "Exhaustion"            = "#A1887F",
  "Activation"            = "#FFD54F",
  "TCR Identity"          = "#9C6FD6"
)

# ── Filter (preserve defined order) ──────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"
adt_heatmap_feats <- names(adt_heatmap_markers)[names(adt_heatmap_markers) %in% rownames(TARA_cd8)]

DefaultAssay(TARA_cd8) <- "RNA"
rna_heatmap_feats <- names(rna_heatmap_markers)[names(rna_heatmap_markers) %in% rownames(TARA_cd8)]

cat("ADT markers:", length(adt_heatmap_feats),
    "/ RNA markers:", length(rna_heatmap_feats), "\n")

# ── Average expression ───────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"
avg_adt_hm <- AverageExpression(
  TARA_cd8, assays = "ADT", features = adt_heatmap_feats,
  group.by = "Manual_Annotation_refined", slot = "data"
)$ADT

DefaultAssay(TARA_cd8) <- "RNA"
avg_rna_hm <- AverageExpression(
  TARA_cd8, assays = "RNA", features = rna_heatmap_feats,
  group.by = "Manual_Annotation_refined", slot = "data"
)$RNA

avg_adt_hm_sc <- scale_01(as.matrix(log10(avg_adt_hm + 1)))
avg_rna_hm_sc <- scale_01(as.matrix(avg_rna_hm))

# ── Fix columns ──────────────────────────────────────────────────────────────
avg_adt_hm_sc <- rename_seurat_cols(avg_adt_hm_sc)
avg_rna_hm_sc <- rename_seurat_cols(avg_rna_hm_sc)

# Short labels with line breaks for long names
cd8_heatmap_labels <- setNames(
  c("Naive 1", "Naive 2", "Naive\nIntermediate",
    "Tscm", "Transitional", "TEMRA/\nCTL",
    "Tex", "γδ T cell"),
  cd8_cluster_names
)

colnames(avg_adt_hm_sc) <- cd8_heatmap_labels[colnames(avg_adt_hm_sc)]
colnames(avg_rna_hm_sc) <- cd8_heatmap_labels[colnames(avg_rna_hm_sc)]

# ── Rename ADT rows to protein names + enforce order ─────────────────────────
rownames(avg_adt_hm_sc) <- ifelse(
  rownames(avg_adt_hm_sc) %in% names(adt_gene_to_protein),
  adt_gene_to_protein[rownames(avg_adt_hm_sc)],
  rownames(avg_adt_hm_sc)
)

adt_row_order <- ifelse(
  adt_heatmap_feats %in% names(adt_gene_to_protein),
  adt_gene_to_protein[adt_heatmap_feats],
  adt_heatmap_feats
)
adt_row_order <- adt_row_order[adt_row_order %in% rownames(avg_adt_hm_sc)]
avg_adt_hm_sc <- avg_adt_hm_sc[adt_row_order, , drop = FALSE]

rna_row_order <- rna_heatmap_feats[rna_heatmap_feats %in% rownames(avg_rna_hm_sc)]
avg_rna_hm_sc <- avg_rna_hm_sc[rna_row_order, , drop = FALSE]

# ── Row annotations ───────────────────────────────────────────────────────────
protein_to_group <- setNames(
  adt_heatmap_markers[adt_heatmap_feats],
  ifelse(adt_heatmap_feats %in% names(adt_gene_to_protein),
         adt_gene_to_protein[adt_heatmap_feats],
         adt_heatmap_feats)
)

adt_annot <- data.frame(
  Group = factor(protein_to_group[rownames(avg_adt_hm_sc)],
                 levels = names(group_colors)),
  row.names = rownames(avg_adt_hm_sc)
)

rna_annot <- data.frame(
  Group = factor(rna_heatmap_markers[rownames(avg_rna_hm_sc)],
                 levels = names(group_colors)),
  row.names = rownames(avg_rna_hm_sc)
)

get_gaps <- function(annot_df) {
  grps <- as.character(annot_df$Group)
  which(grps[-length(grps)] != grps[-1])
}

annot_colors <- list(Group = group_colors)

heatmap_colors <- colorRampPalette(
  c("#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C")
)(100)

# ── Clean matrices ────────────────────────────────────────────────────────────
adt_mat <- matrix(as.numeric(avg_adt_hm_sc), nrow = nrow(avg_adt_hm_sc),
                  dimnames = dimnames(avg_adt_hm_sc))
rna_mat <- matrix(as.numeric(avg_rna_hm_sc), nrow = nrow(avg_rna_hm_sc),
                  dimnames = dimnames(avg_rna_hm_sc))

cat("ADT matrix:", nrow(adt_mat), "x", ncol(adt_mat),
    "range:", range(adt_mat), "\n")
cat("RNA matrix:", nrow(rna_mat), "x", ncol(rna_mat),
    "range:", range(rna_mat), "\n")

# ── Heatmaps (large text, horizontal column labels) ───────────────────────────
p_adt <- pheatmap(
  adt_mat, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
  color = heatmap_colors, border_color = "white",
  annotation_row = adt_annot, annotation_colors = annot_colors,
  annotation_names_row = FALSE, gaps_row = get_gaps(adt_annot),
  cellwidth = 80, cellheight = 24, fontsize = 18, fontsize_row = 16,
  fontsize_col = 16, angle_col = 0, main = "Protein (ADT)", silent = TRUE
)

p_rna <- pheatmap(
  rna_mat, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
  color = heatmap_colors, border_color = "white",
  annotation_row = rna_annot, annotation_colors = annot_colors,
  annotation_names_row = FALSE, gaps_row = get_gaps(rna_annot),
  cellwidth = 80, cellheight = 24, fontsize = 18, fontsize_row = 16,
  fontsize_col = 16, angle_col = 0, main = "mRNA (RNA)", silent = TRUE
)

# ── Shared group legend (standalone) ──────────────────────────────────────────
legend_df <- data.frame(
  Group = factor(names(group_colors), levels = names(group_colors)),
  x = seq_along(group_colors),
  y = 1
)

p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = Group)) +
  geom_tile(width = 0.8, height = 0.5) +
  scale_fill_manual(values = group_colors, name = "Functional Group") +
  theme_void(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.6, "cm")) +
  guides(fill = guide_legend(nrow = 1))

# Extract just the legend
legend_grob <- cowplot::get_legend(p_legend)

# Horizontal (stacked: ADT on top, RNA below, shared legend at bottom)
png(file.path(fig1_dir, "Fig1C_Heatmap.png"),
    width = 20, height = 44, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(19, 19, 1.5), "in"))))
pushViewport(viewport(layout.pos.row = 1)); grid.draw(p_adt$gtable); popViewport()
pushViewport(viewport(layout.pos.row = 2)); grid.draw(p_rna$gtable); popViewport()
pushViewport(viewport(layout.pos.row = 3)); grid.draw(legend_grob); popViewport()
dev.off()

# Individual panels
png(file.path(fig1_dir, "Fig1C_ADT_Heatmap.png"),
    width = 18, height = 22, units = "in", res = 300, bg = "white")
grid.draw(p_adt$gtable)
dev.off()

png(file.path(fig1_dir, "Fig1C_RNA_Heatmap.png"),
    width = 18, height = 22, units = "in", res = 300, bg = "white")
grid.draw(p_rna$gtable)
dev.off()

# ── Standalone group legend (horizontal, for bottom of figure) ────────────────
# Add line breaks to long names
legend_labels <- c(
  "Naive/\nQuiescence", "Memory/\nStemness", "Co-stimulation/\nHoming",
  "Effector/\nTEMRA", "NK-like/\nInnate", "Exhaustion",
  "Activation", "TCR\nIdentity"
)

legend_df <- data.frame(
  Group = factor(names(group_colors), levels = names(group_colors)),
  x = 1, y = seq_along(group_colors)
)

legend_df <- data.frame(
  Group = factor(names(group_colors), levels = names(group_colors)),
  y = seq(1, by = 0.7, length.out = length(group_colors)),
  label = legend_labels
)

p_legend_only <- ggplot(legend_df, aes(x = 1, y = y)) +
  geom_point(aes(fill = Group), shape = 22, size = 10, stroke = 0) +
  geom_text(aes(x = 0.82, label = label), size = 5.5, fontface = "bold",
            lineheight = 0.8, vjust = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(limits = c(0.6, 1.1)) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(0, 2, 0, 2))

ggsave(file.path(fig1_dir, "Fig1C_Group_Legend.png"),
       p_legend_only, width = 14, height = 1.4, dpi = 300, bg = "white")

message("✓ Fig 1C heatmaps + legend saved\n")


################################################################################
# FIG 1D — CD8 CLUSTER UMAP (IMPROVED FOR PUBLICATION)
################################################################################

message("Generating Fig 1D...")

# Use short labels for display
TARA_sub$CD8_highlight <- ifelse(
  TARA_sub$Manual_Annotation_refined %in% cd8_cluster_names,
  cd8_short_labels[as.character(TARA_sub$Manual_Annotation_refined)],
  "Other"
)
TARA_sub$CD8_highlight <- factor(
  TARA_sub$CD8_highlight,
  levels = c(cd8_short_labels, "Other")
)

# IMPROVED: Colors keyed by short labels - distinct and visible
highlight_cols <- setNames(
  c("#2166AC", "#67A9CF", "#D6604D", "#1A9850", 
    "#FDAE61", "#B2182B", "#762A83", "#CC79A7", "#E0E0E0"),
  c("Naive 1", "Naive 2", "Naive Intermediate", "Tscm", 
    "Transitional", "TEMRA/CTL", "Tex", "γδ T cell", "Other")
)

# IMPROVED: Build UMAP with ggplot2 for full control
umap_cd8_coords <- as.data.frame(Embeddings(TARA_sub, reduction = "wnn.umap"))
colnames(umap_cd8_coords) <- c("UMAP1", "UMAP2")
umap_cd8_coords$Cluster <- TARA_sub$CD8_highlight

# Split into CD8 and Other for layered plotting
umap_other <- umap_cd8_coords %>% filter(Cluster == "Other")
umap_cd8_only <- umap_cd8_coords %>% filter(Cluster != "Other")

# Calculate centroids for CD8 clusters only
centroids_cd8 <- umap_cd8_only %>%
  group_by(Cluster) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .groups = "drop"
  )

# UMAP arrow positions
x_arrow <- min(umap_cd8_coords$UMAP1, na.rm = TRUE) + 1
y_arrow <- min(umap_cd8_coords$UMAP2, na.rm = TRUE) + 1

p_cd8_umap_colored <- ggplot() +
  # Background: Other cells (gray)
  geom_point(data = umap_other, 
             aes(x = UMAP1, y = UMAP2),
             color = "#E0E0E0", size = 0.5, alpha = 0.4) +
  # Foreground: CD8 clusters - larger points
  geom_point(data = umap_cd8_only,
             aes(x = UMAP1, y = UMAP2, color = Cluster),
             size = 1.2, alpha = 0.8) +
  scale_color_manual(values = highlight_cols) +
  # Labels: white fill, colored text and border matching cluster
  geom_label_repel(
    data = centroids_cd8,
    aes(x = UMAP1, y = UMAP2, label = Cluster, color = Cluster),
    fill = "white",
    size = 12,
    fontface = "bold",
    label.size = 0.8,
    label.padding = unit(0.4, "lines"),
    box.padding = 1.2,
    point.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.8,
    min.segment.length = 0,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  # UMAP arrows - bottom left, large
  annotate("segment", 
           x = x_arrow, xend = x_arrow + 3.5,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow + 1.75, y = y_arrow - 1.3,
           label = "UMAP1", size = 7, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + 3.5,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow - 1.3, y = y_arrow + 1.75,
           label = "UMAP2", size = 7, fontface = "bold", angle = 90) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

p_cd8_umap_black <- ggplot() +
  # Background: Other cells (gray)
  geom_point(data = umap_other, 
             aes(x = UMAP1, y = UMAP2),
             color = "#E0E0E0", size = 0.5, alpha = 0.4) +
  # Foreground: CD8 clusters - larger points
  geom_point(data = umap_cd8_only,
             aes(x = UMAP1, y = UMAP2, color = Cluster),
             size = 1.2, alpha = 0.8) +
  scale_color_manual(values = highlight_cols) +
  # Labels: white fill, black text and border
  geom_label_repel(
    data = centroids_cd8,
    aes(x = UMAP1, y = UMAP2, label = Cluster),
    fill = "white",
    color = "black",
    size = 12,
    fontface = "bold",
    label.size = 0.8,
    label.padding = unit(0.4, "lines"),
    box.padding = 1.2,
    point.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.8,
    min.segment.length = 0,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  # UMAP arrows - bottom left, large
  annotate("segment", 
           x = x_arrow, xend = x_arrow + 3.5,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow + 1.75, y = y_arrow - 1.3,
           label = "UMAP1", size = 7, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + 3.5,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
           linewidth = 1.2, color = "black") +
  annotate("text", x = x_arrow - 1.3, y = y_arrow + 1.75,
           label = "UMAP2", size = 7, fontface = "bold", angle = 90) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

ggsave(file.path(fig1_dir, "Fig1D_CD8_Cluster_UMAP_colored.png"),
       p_cd8_umap_colored, width = 12, height = 11, dpi = 300, bg = "white")

ggsave(file.path(fig1_dir, "Fig1D_CD8_Cluster_UMAP_black.png"),
       p_cd8_umap_black, width = 12, height = 11, dpi = 300, bg = "white")

message("✓ Fig 1D saved (both colored and black versions)\n")


################################################################################
# FIG 1E — TEMRA/CTL PROPORTION (significant cluster only)
# FIG 1F — Tex VL CORRELATION (significant cluster only)
#
# Full versions → Supplementary 1
################################################################################

message("Generating Fig 1D-F + Supplementary...")

# ── Build frequency tables ────────────────────────────────────────────────────
meta_sub <- as.data.frame(TARA_sub@meta.data)

sample_totals <- meta_sub %>%
  group_by(orig.ident, Condition) %>%
  summarise(total_pbmc = n(), .groups = "drop")

cluster_freq <- meta_sub %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  group_by(orig.ident, Condition, Manual_Annotation_refined) %>%
  summarise(cluster_cells = n(), .groups = "drop") %>%
  group_by(orig.ident, Condition) %>%
  complete(
    Manual_Annotation_refined = cd8_cluster_names,
    fill = list(cluster_cells = 0)
  ) %>%
  ungroup() %>%
  left_join(sample_totals, by = c("orig.ident", "Condition")) %>%
  mutate(
    percent = 100 * cluster_cells / total_pbmc,
    Facet_Label = factor(
      cd8_short_labels[as.character(Manual_Annotation_refined)],
      levels = cd8_short_labels
    )
  )

# ── Statistics (all clusters) ─────────────────────────────────────────────────
facet_max <- cluster_freq %>%
  group_by(Facet_Label) %>%
  summarise(max_percent = max(percent, na.rm = TRUE), .groups = "drop")

stat_df_all <- cluster_freq %>%
  group_by(Facet_Label) %>%
  wilcox_test(percent ~ Condition) %>%
  left_join(facet_max, by = "Facet_Label") %>%
  mutate(
    xmin = case_when(
      group1 == "HUU" & group2 == "HEU" ~ 1,
      group1 == "HUU" & group2 == "HEI" ~ 1,
      group1 == "HEU" & group2 == "HEI" ~ 2
    ),
    xmax = case_when(
      group1 == "HUU" & group2 == "HEU" ~ 2,
      group1 == "HUU" & group2 == "HEI" ~ 3,
      group1 == "HEU" & group2 == "HEI" ~ 3
    ),
    y.position = case_when(
      group1 == "HUU" & group2 == "HEU" ~ max_percent * 1.10,
      group1 == "HUU" & group2 == "HEI" ~ max_percent * 1.22,
      group1 == "HEU" & group2 == "HEI" ~ max_percent * 1.16
    ),
    p.signif = case_when(
      p <= 0.0001 ~ "****", p <= 0.001 ~ "***",
      p <= 0.01   ~ "**",   p <= 0.05  ~ "*",
      TRUE        ~ "ns"
    )
  ) %>%
  ungroup()

stat_df_sig <- stat_df_all %>% filter(p.signif != "ns")

# ── Fig 1E: TEMRA/CTL only ───────────────────────────────────────────────────
temra_freq <- cluster_freq %>% filter(Facet_Label == "TEMRA/CTL")
temra_stat <- stat_df_sig %>% filter(Facet_Label == "TEMRA/CTL")

p_temra <- ggplot(temra_freq, aes(x = Condition, y = percent,
                                  fill = Condition, color = Condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.8) +
  stat_pvalue_manual(temra_stat, label = "p.signif", tip.length = 0.01,
                     bracket.size = 0.6, size = 6) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  labs(x = NULL, y = "% of total PBMC", title = "TEMRA/CTL CD8") +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x  = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(face = "bold", size = 16)
  )

ggsave(file.path(fig1_dir, "Fig1E_TEMRA_Proportion.png"),
       p_temra, width = 5, height = 6, dpi = 600, bg = "white")

# ── Supplementary: All cluster proportions ────────────────────────────────────
p_prop_all <- ggplot(cluster_freq, aes(x = Condition, y = percent,
                                       fill = Condition, color = Condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", nrow = 2) +
  stat_pvalue_manual(stat_df_sig, label = "p.signif", tip.length = 0.01,
                     bracket.size = 0.6, size = 5) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  labs(x = NULL, y = "% of total PBMC") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 13),
    axis.text.x     = element_text(face = "bold", size = 12),
    axis.text.y     = element_text(size = 12),
    axis.title.y    = element_text(face = "bold", size = 14)
  )

ggsave(file.path(supp1_dir, "SuppFig1_All_Cluster_Proportions.png"),
       p_prop_all, width = 24, height = 12, dpi = 600, bg = "white")

message("✓ Fig 1E + Supplementary proportions saved\n")


################################################################################
# FIG 1F — Tex VL CORRELATION + Supplementary (all clusters)
################################################################################

message("Generating Fig 1F + Supplementary VL correlations...")

meta_hei <- meta_sub %>%
  filter(Condition == "HEI", !is.na(Viral_Load_num), Viral_Load_num > 0)

sample_totals_hei <- meta_hei %>%
  group_by(orig.ident) %>%
  summarise(total_pbmc = n(), .groups = "drop")

sample_cd8_totals <- meta_hei %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  group_by(orig.ident) %>%
  summarise(total_cd8 = n(), .groups = "drop")

sample_vl <- meta_hei %>%
  group_by(orig.ident) %>%
  summarise(Viral_Load_num = first(Viral_Load_num), .groups = "drop") %>%
  mutate(log10_VL = log10(Viral_Load_num))

cluster_vl_df <- meta_hei %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  group_by(orig.ident, Manual_Annotation_refined) %>%
  summarise(cluster_cells = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  complete(
    Manual_Annotation_refined = cd8_cluster_names,
    fill = list(cluster_cells = 0)
  ) %>%
  ungroup() %>%
  left_join(sample_totals_hei, by = "orig.ident") %>%
  left_join(sample_cd8_totals, by = "orig.ident") %>%
  left_join(sample_vl, by = "orig.ident") %>%
  mutate(
    pct_pbmc = 100 * cluster_cells / total_pbmc,
    pct_cd8  = 100 * cluster_cells / total_cd8,
    Facet_Label = factor(
      cd8_short_labels[as.character(Manual_Annotation_refined)],
      levels = cd8_short_labels
    )
  )

# ── Fig 1F: Tex only ─────────────────────────────────────────────────────────
tex_vl <- cluster_vl_df %>% filter(Facet_Label == "Tex")

tex_cor <- tex_vl %>%
  cor_test(log10_VL, pct_pbmc, method = "spearman")

p_tex_vl <- ggplot(tex_vl, aes(x = log10_VL, y = pct_pbmc)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              color = "grey40", fill = "grey85", alpha = 0.3) +
  geom_point(shape = 21, size = 4, stroke = 0.5,
             fill = "#C44E52", color = "white") +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("rho = ", sprintf("%.2f", tex_cor$cor),
                          "\np = ", signif(tex_cor$p, 2)),
           hjust = -0.1, vjust = 1.3, size = 6, fontface = "bold",
           color = "grey15", lineheight = 1.2) +
  labs(x = expression(log[10] ~ "viral load"),
       y = "% of total PBMC",
       title = "Tex CD8") +
  theme_classic(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text    = element_text(size = 14, color = "black"),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 6)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 6)),
    legend.position = "none"
  )

ggsave(file.path(fig1_dir, "Fig1F_Tex_VL_Correlation.png"),
       p_tex_vl, width = 6, height = 6, dpi = 600, bg = "white")

cat("=== Fig 1F: Tex VL correlation ===\n")
print(tex_cor[, c("cor", "statistic", "p")])
cat("\n")

# ── Supplementary: All clusters VL correlations ──────────────────────────────
make_vl_plot <- function(df, y_var, y_label, filename) {
  
  cor_stats <- df %>%
    group_by(Facet_Label) %>%
    cor_test(log10_VL, !!sym(y_var), method = "spearman") %>%
    mutate(
      stat_label = paste0("rho = ", sprintf("%.2f", cor),
                          "\np = ", signif(p, 2)),
      x = -Inf, y = Inf
    ) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = log10_VL, y = .data[[y_var]])) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
                color = "grey40", fill = "grey85", alpha = 0.3) +
    geom_point(shape = 21, size = 3.5, stroke = 0.5,
               fill = "#C44E52", color = "white") +
    geom_text(data = cor_stats,
              aes(x = x, y = y, label = stat_label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.2,
              size = 5, lineheight = 1.3, color = "grey15",
              fontface = "bold") +
    facet_wrap(~ Facet_Label, scales = "free", nrow = 2) +
    scale_x_continuous(breaks = pretty_breaks(n = 4)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = expression(log[10] ~ "viral load"), y = y_label) +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 13),
      axis.text        = element_text(size = 12, color = "black"),
      axis.title.x     = element_text(face = "bold", size = 14,
                                      margin = margin(t = 6)),
      axis.title.y     = element_text(face = "bold", size = 14,
                                      margin = margin(r = 6)),
      legend.position  = "none",
      panel.spacing    = unit(0.8, "lines")
    )
  
  ggsave(file.path(supp1_dir, filename),
         p, width = 24, height = 12, dpi = 600, bg = "white")
  
  cat("=== Correlations:", y_label, "===\n")
  print(cor_stats[, c("Facet_Label", "cor", "p")])
  cat("\n")
  
  return(p)
}

# ── Supplementary S1G: VL correlations % PBMC (7 clusters, exclude Tex) ──────
cluster_vl_noTex <- cluster_vl_df %>% filter(Facet_Label != "Tex")

p_vl_pbmc <- make_vl_plot(
  cluster_vl_noTex, "pct_pbmc", "% of total PBMC",
  "S1G_VL_vs_PctPBMC.png"
)

# ── Supplementary S1H: VL correlations % CD8 (all 8 — none in main fig) ─────
p_vl_cd8 <- make_vl_plot(
  cluster_vl_df, "pct_cd8", "% of CD8 T cells",
  "S1H_VL_vs_PctCD8.png"
)

message("✓ Fig 1F + S1G-H VL correlations saved\n")


################################################################################
# SUPPLEMENTARY FIGURE S1A — BROAD PBMC ANNOTATION VALIDATION
#
# Dot plot of canonical lineage markers grouped by broad cell type (matching 1A)
################################################################################

message("Generating Supplementary S1A...")

DefaultAssay(TARA_sub) <- "ADT"

lineage_markers_adt <- c(
  # T cell
  "CD3D", "CD4", "CD8A",
  # B cell
  "CD19", "MS4A1", "CD79B", "IGHM", "IGHD",
  # NK
  "NCAM1", "KLRF1", "KLRD1",
  # Monocyte
  "CD14", "FCGR3A", "ITGAM",
  # DC
  "FCER1A", "CD1C", "CLEC4C",
  # Activation / identity
  "HLA-DRA", "CD38", "IL7R",
  # GammaDelta
  "TCR-AB", "TCR-vD2",
  # Other
  "CD56", "PDCD1", "TIGIT",
  "CD45RA", "CD45RO",
  "CD27", "CD28",
  "ITGB7", "FAS"
)

lineage_adt_available <- lineage_markers_adt[lineage_markers_adt %in% rownames(TARA_sub)]

# Map to protein display names
adt_display <- c(
  "CD3D" = "CD3", "CD4" = "CD4", "CD8A" = "CD8a",
  "CD19" = "CD19", "MS4A1" = "CD20", "CD79B" = "CD79b", "IGHM" = "IgM", "IGHD" = "IgD",
  "NCAM1" = "CD56", "KLRF1" = "NKp80", "KLRD1" = "CD94",
  "CD14" = "CD14", "FCGR3A" = "CD16", "ITGAM" = "CD11b",
  "FCER1A" = "FcER1a", "CD1C" = "CD1c", "CLEC4C" = "BDCA-2",
  "HLA-DRA" = "HLA-DR", "CD38" = "CD38", "IL7R" = "CD127",
  "TCR-AB" = "TCRab", "TCR-vD2" = "TCR Vd2",
  "CD56" = "CD56", "PDCD1" = "PD-1", "TIGIT" = "TIGIT",
  "CD45RA" = "CD45RA", "CD45RO" = "CD45RO",
  "CD27" = "CD27", "CD28" = "CD28",
  "ITGB7" = "Integrin b7", "FAS" = "CD95"
)

# Order broad cell types to match Fig 1A
broad_order <- c("CD4 T cells", "CD8 T cells", "γδ T cells", "DN T cells",
                 "NK cells", "B cells", "Plasmablasts", "Monocytes",
                 "pDC", "APC", "Other")
TARA_sub$Broad_CellType <- factor(TARA_sub$Broad_CellType, levels = broad_order)

p_s1a <- DotPlot(
  TARA_sub, features = lineage_adt_available,
  group.by = "Broad_CellType",
  cols = c("lightgrey", "#08519C"), dot.scale = 8
) + RotatedAxis() +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(adt_display),
                                               adt_display[x], x)) +
  labs(title = "S1A: Broad PBMC annotation — surface protein markers (ADT)",
       y = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 13, face = "bold")
  )

ggsave(file.path(supp1_dir, "S1A_PBMC_Lineage_DotPlot.png"),
       p_s1a, width = 18, height = 8, dpi = 300, bg = "white")

message("✓ S1A saved\n")


################################################################################
# SUPPLEMENTARY FIGURE S1B — CD8 SUB-GATING: CD45RO vs FAS
################################################################################

message("Generating Supplementary S1B...")

cl1_clusters <- c("1a: Naive 1 CD8", "1b: Tscm CD8", "1c: Naive Intermediate CD8")
TARA_cl1 <- subset(TARA_cd8,
                   subset = Manual_Annotation_refined %in% cl1_clusters)

DefaultAssay(TARA_cl1) <- "ADT"
adt_cl1 <- GetAssayData(TARA_cl1, slot = "data")

gate_df <- data.frame(
  CD45RO = as.numeric(adt_cl1["CD45RO", ]),
  CD45RA = as.numeric(adt_cl1["CD45RA", ]),
  FAS    = as.numeric(adt_cl1["FAS",    ]),
  Gate   = as.character(TARA_cl1$Manual_Annotation_refined),
  stringsAsFactors = FALSE
)
gate_df$Gate <- cd8_short_labels[gate_df$Gate]
gate_df$Gate <- factor(gate_df$Gate, levels = c("Naive 1", "Tscm", "Naive Intermediate"))

gate_cols_short <- c(
  "Naive 1"            = cluster_cols[["1a: Naive 1 CD8"]],
  "Tscm"               = cluster_cols[["1b: Tscm CD8"]],
  "Naive Intermediate" = cluster_cols[["1c: Naive Intermediate CD8"]]
)

p_s1b_1 <- ggplot(gate_df, aes(x = FAS, y = CD45RO, color = Gate)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_manual(values = gate_cols_short) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "CD45RO vs CD95",
       x = "CD95/FAS (DSB)", y = "CD45RO (DSB)") +
  theme_classic(base_size = 28) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 22),
        plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24, face = "bold"))

p_s1b_2 <- ggplot(gate_df, aes(x = CD45RA, y = CD45RO, color = Gate)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_manual(values = gate_cols_short) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "CD45RO vs CD45RA",
       x = "CD45RA (DSB)", y = "CD45RO (DSB)") +
  theme_classic(base_size = 28) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 22),
        plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24, face = "bold"))

p_s1b <- p_s1b_1 | p_s1b_2
p_s1b <- p_s1b + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(supp1_dir, "S1B_CD8_SubGating.png"),
       p_s1b, width = 18, height = 9, dpi = 300, bg = "white")

message("✓ S1B saved\n")


################################################################################
# SUPPLEMENTARY FIGURE S1C — NAIVE 1 vs NAIVE 2: PROTEIN DIFFERENCES
################################################################################

message("Generating Supplementary S1C...")

DefaultAssay(TARA_cd8) <- "ADT"
adt_cd8_data <- GetAssayData(TARA_cd8, slot = "data")

n1_cells <- WhichCells(TARA_cd8,
                       expression = Manual_Annotation_refined == "1a: Naive 1 CD8")
n2_cells <- WhichCells(TARA_cd8,
                       expression = Manual_Annotation_refined == "6: Naive 2 CD8")

pct_n1 <- rowMeans(adt_cd8_data[, n1_cells] > 0) * 100
pct_n2 <- rowMeans(adt_cd8_data[, n2_cells] > 0) * 100

naive_diff <- data.frame(
  Marker = names(pct_n1),
  Naive_1 = round(pct_n1, 1),
  Naive_2 = round(pct_n2, 1),
  Diff = round(pct_n1 - pct_n2, 1),
  stringsAsFactors = FALSE
)
naive_diff <- naive_diff[order(-abs(naive_diff$Diff)), ]

top20 <- head(naive_diff, 20)

top20$Protein <- ifelse(
  top20$Marker %in% names(adt_gene_to_protein),
  adt_gene_to_protein[top20$Marker],
  top20$Marker
)

top20$Protein <- factor(top20$Protein, levels = rev(top20$Protein))
top20$Direction <- ifelse(top20$Diff > 0, "Higher in Naive 1", "Higher in Naive 2")

p_s1c <- ggplot(top20, aes(x = Diff, y = Protein, fill = Direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Higher in Naive 1" = cluster_cols[["1a: Naive 1 CD8"]],
                               "Higher in Naive 2" = cluster_cols[["6: Naive 2 CD8"]])) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  labs(x = "Difference in % expressing (Naive 1 - Naive 2)",
       y = NULL) +
  theme_classic(base_size = 28) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 24, face = "bold")
  )

ggsave(file.path(supp1_dir, "S1C_Naive1_vs_Naive2_ADT.png"),
       p_s1c, width = 14, height = 14, dpi = 300, bg = "white")

message("✓ S1C saved\n")


################################################################################
# SUPPLEMENTARY FIGURE S1D — CLUSTER PROPORTIONS (7, excl TEMRA/CTL, 1 row)
################################################################################

message("Generating Supplementary S1D...")

cluster_freq_noTEMRA <- cluster_freq %>%
  filter(Facet_Label != "TEMRA/CTL")

stat_df_noTEMRA <- stat_df_all %>%
  filter(Facet_Label != "TEMRA/CTL", p.signif != "ns")

p_s1d <- ggplot(cluster_freq_noTEMRA, aes(x = Condition, y = percent,
                                          fill = Condition, color = Condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 3.5, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", nrow = 1) +
  {if (nrow(stat_df_noTEMRA) > 0)
    stat_pvalue_manual(stat_df_noTEMRA, label = "p.signif", tip.length = 0.01,
                       bracket.size = 0.8, size = 8)} +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  labs(x = NULL, y = "% of total PBMC") +
  theme_classic(base_size = 28) +
  theme(
    legend.position = "none",
    strip.text   = element_text(face = "bold", size = 22),
    axis.text.x  = element_text(face = "bold", size = 20),
    axis.text.y  = element_text(size = 20),
    axis.title.y = element_text(face = "bold", size = 24)
  )

ggsave(file.path(supp1_dir, "S1D_Cluster_Proportions.png"),
       p_s1d, width = 36, height = 9, dpi = 600, bg = "white")

message("✓ S1D saved\n")


################################################################################
# SUPPLEMENTARY FIGURE S1E — VL CORRELATIONS % PBMC (7, excl Tex, 1 row)
################################################################################

message("Generating Supplementary S1E...")

cluster_vl_noTex <- cluster_vl_df %>% filter(Facet_Label != "Tex")

cor_stats_s1e <- cluster_vl_noTex %>%
  group_by(Facet_Label) %>%
  cor_test(log10_VL, pct_pbmc, method = "spearman") %>%
  mutate(
    stat_label = paste0("rho = ", sprintf("%.2f", cor),
                        "\np = ", signif(p, 2)),
    x = -Inf, y = Inf
  ) %>%
  ungroup()

p_s1e <- ggplot(cluster_vl_noTex, aes(x = log10_VL, y = pct_pbmc)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              color = "grey40", fill = "grey85", alpha = 0.3) +
  geom_point(shape = 21, size = 4.5, stroke = 0.5,
             fill = "#C44E52", color = "white") +
  geom_text(data = cor_stats_s1e,
            aes(x = x, y = y, label = stat_label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.2,
            size = 7, lineheight = 1.2, color = "grey15",
            fontface = "bold") +
  facet_wrap(~ Facet_Label, scales = "free", nrow = 1) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = expression(log[10] ~ "viral load"), y = "% of total PBMC") +
  theme_classic(base_size = 28) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 22),
    axis.text        = element_text(size = 20, color = "black"),
    axis.title.x     = element_text(face = "bold", size = 24,
                                    margin = margin(t = 6)),
    axis.title.y     = element_text(face = "bold", size = 24,
                                    margin = margin(r = 6)),
    legend.position  = "none",
    panel.spacing    = unit(0.8, "lines")
  )

ggsave(file.path(supp1_dir, "S1E_VL_vs_PctPBMC.png"),
       p_s1e, width = 36, height = 9, dpi = 600, bg = "white")

cat("=== S1E VL correlations ===\n")
print(cor_stats_s1e[, c("Facet_Label", "cor", "p")])
cat("\n")

message("✓ S1E saved\n")


################################################################################
# SUMMARY
################################################################################

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " All panels saved.\n",
        "══════════════════════════════════════════════════════════════\n",
        "\n",
        " MAIN FIGURE 1:  ", fig1_dir, "\n",
        "   1A: Total PBMC UMAP (broad cell types)\n",
        "   1B: CD3/CD8 feature plots (mRNA + protein)\n",
        "   1C: ADT + RNA heatmap (32 markers each)\n",
        "   1D: CD8 cluster UMAP (non-CD8 grayed)\n",
        "   1E: TEMRA/CTL proportion (significant)\n",
        "   1F: Tex VL correlation (significant)\n",
        "\n",
        " SUPPLEMENTARY S1:  ", supp1_dir, "\n",
        "   S1A: Broad PBMC lineage protein dot plot\n",
        "   S1B: CD8 sub-gating (CD45RO vs FAS)\n",
        "   S1C: Naive 1 vs Naive 2 ADT protein differences\n",
        "   S1D: Cluster proportions (7, excl TEMRA/CTL)\n",
        "   S1E: VL correlations % PBMC (7, excl Tex)\n",
        "══════════════════════════════════════════════════════════════\n"
)