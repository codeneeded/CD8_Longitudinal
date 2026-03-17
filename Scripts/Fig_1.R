################################################################################
# FIGURE 1 — CD8 T Cell Landscape in HIV-Exposed Infants
# Manuscript: CD8 Longitudinal (TARA Cohort)
#
# Figure order:
#   Fig 1A — Methods schematic (created manually, not coded here)
#   Fig 1B — UMAP: all clusters with numeric labels
#   Fig 1C — ADT + RNA heatmap: CD8 cluster identity validation
#   Fig 1D — CD8 cluster proportions (% of PBMC): HUU vs HEU vs HEI
#   Fig 1E — CD8 cluster % vs viral load correlations (HEI only)
#
# Cluster ordering throughout (by cluster number):
#   Cluster 1  → Tcm/Tscm CD8 (Memory CD8 T cell)
#   Cluster 6  → Naïve CD8 T cell
#   Cluster 8  → TEMRA/CTL
#   Cluster 9  → TRDV1+ γδ T cell
#   Cluster 27 → Tpex CD8 (GZMK+ CD8 T cell)
################################################################################

# ── Libraries ──────────────────────────────────────────────────────────────────
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
library(scRepertoire)
library(monocle3)
library(SeuratWrappers)
library(qs2)
library(Matrix)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(EnhancedVolcano)
library(rstatix)
library(grid)
library(scales)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")

fig1_dir  <- file.path(base_dir, "Manuscript", "Fig 1")

out_umap  <- file.path(fig1_dir, "01_UMAP")
out_heat  <- file.path(fig1_dir, "02_ADT_RNA_Heatmap")
out_prop  <- file.path(fig1_dir, "03_Cluster_Proportions")
out_vl    <- file.path(fig1_dir, "04_Viral_Load_Correlations")

for (d in c(out_umap, out_heat, out_prop, out_vl)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Read data ──────────────────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted.qs2"))


################################################################################
# CLUSTER ANNOTATION: Merge αβ TCR+ contaminants from Cluster 9 into Cluster 8
#
# Rationale:
#   Cluster 9 is predominantly TRDV1+ γδ T cells. However, a subset of cells
#   in this cluster carry αβ TCRs (TRBV+, zero TRDV/TRDC) and are
#   misassigned. These are reassigned to Cluster 8 (CTL-like / TEMRA).
#   Remaining Cluster 9 cells (no αβ TCR) are confirmed genuine γδ T cells.
################################################################################

TARA_ALL$Manual_Annotation_refined <- case_when(
  # Cluster 9, no αβ TCR → confirmed TRDV1+ γδ T cell
  TARA_ALL$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_ALL$has_TCR == FALSE
  ~ "9: TRDV1+ γδ T cell",
  # Cluster 9, αβ TCR present → contaminant; reassign to Cluster 8
  TARA_ALL$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_ALL$has_TCR == TRUE
  ~ "8: CTL-like",
  # All other clusters unchanged
  TRUE
  ~ as.character(TARA_ALL$Manual_Annotation)
)

TARA_ALL$Manual_Annotation_refined <- factor(TARA_ALL$Manual_Annotation_refined)

# Verify merge
cat("=== Post-merge cluster sizes (Clusters 8 & 9) ===\n")
print(table(TARA_ALL$Manual_Annotation_refined)[
  grepl("^8|^9", names(table(TARA_ALL$Manual_Annotation_refined)))
])

# Numeric cluster label (strip everything after the colon) — used for UMAP
TARA_ALL$Cluster_Number_refined <- gsub("^([0-9]+):.*", "\\1",
                                        TARA_ALL$Manual_Annotation_refined)


################################################################################
# SUBSETS
#
# TARA_fig1   — Fig 1B only: PreART_Entry, HEU, HUU (relabelled as HEI/HEU/HUU)
# TARA_cd8    — All CD8 clusters (used for heatmap, proportions, VL correlations)
################################################################################

# ── Fig 1B subset ─────────────────────────────────────────────────────────────
keep_groups <- c("PreART_Entry", "HEU", "HUU")

TARA_fig1 <- subset(TARA_ALL, subset = Timepoint_Group %in% keep_groups)
TARA_fig1$Condition <- factor(
  case_when(
    TARA_fig1$Timepoint_Group == "PreART_Entry" ~ "HEI",
    TRUE ~ as.character(TARA_fig1$Timepoint_Group)
  ),
  levels = c("HEI", "HEU", "HUU")
)

# ── CD8 cluster subset ────────────────────────────────────────────────────────
# Ordered by cluster number: 1, 6, 8, 9, 27
cd8_cluster_names <- c(
  "1: Memory CD8 T cell",
  "6: Naïve CD8 T cell",
  "8: CTL-like",
  "9: TRDV1+ γδ T cell",
  "27: GZMK+ CD8 T cell"
)

TARA_cd8 <- subset(TARA_ALL,
                   subset = Manual_Annotation_refined %in% cd8_cluster_names)

# Display labels used in all figures (ordered by cluster number)
cd8_display_labels <- c(
  "1: Memory CD8 T cell"   = "Tcm/Tscm CD8\n(Cluster 1)",
  "6: Naïve CD8 T cell"    = "Naïve CD8\n(Cluster 6)",
  "8: CTL-like"            = "TEMRA/CTL\n(Cluster 8)",
  "9: TRDV1+ γδ T cell"    = "TRDV1+ γδ\n(Cluster 9)",
  "27: GZMK+ CD8 T cell"   = "Tpex CD8\n(Cluster 27)"
)

# Short labels for strip text in faceted plots (same order)
cd8_short_labels <- c(
  "1: Memory CD8 T cell"   = "Tcm/Tscm CD8",
  "6: Naïve CD8 T cell"    = "Naïve CD8",
  "8: CTL-like"            = "TEMRA/CTL",
  "9: TRDV1+ γδ T cell"    = "TRDV1+ γδ",
  "27: GZMK+ CD8 T cell"   = "Tpex CD8"
)

# ── Shared color palettes ──────────────────────────────────────────────────────
cond_cols <- c("HUU" = "#4E79A7", "HEU" = "#F28E2B", "HEI" = "#E15759")


################################################################################
# FIG 1B — UMAP: All clusters, numeric labels
################################################################################

p_umap_all <- DimPlot2(
  TARA_fig1,
  reduction  = "wnn.umap",
  group.by   = "Cluster_Number_refined",
  cols       = "default",
  label      = TRUE,
  box        = TRUE,
  label.size = 6,
  repel      = TRUE,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
)

ggsave(file.path(out_umap, "Fig1B_UMAP_AllClusters.png"),
       p_umap_all, width = 6, height = 6, dpi = 300, bg = "white")

message("✓ Fig 1B saved")


################################################################################
# FIG 1C — ADT + RNA heatmap: CD8 cluster identity validation
#
# Columns ordered by cluster number: 1, 6, 8, 9, 27
# Rows grouped by biological function; within-group clustering by Ward D2
################################################################################

# ── Marker definitions (order = display order in heatmap rows) ─────────────────
# Named vector: name = marker feature, value = group label
adt_marker_groups <- c(
  # Naïve
  "CD7"       = "Naïve",
  "SELL"      = "Naïve",
  "CD45RA"    = "Naïve",
  # TCR Identity
  "TCR-vA7.2" = "TCR Identity",
  "TCR-AB"    = "TCR Identity",
  "TCR-vD2"   = "TCR Identity",
  # Memory / Stemness
  "IL7R"      = "Memory/Stemness",
  "CD38"      = "Memory/Stemness",
  "NT5E"      = "Memory/Stemness",
  "ENTPD1"    = "Memory/Stemness",
  "CD45RO"    = "Memory/Stemness",
  # Co-stimulation
  "CD27"      = "Co-stimulation",
  "CD28"      = "Co-stimulation",
  # Effector / TEMRA
  "B3GAT1"    = "Effector/TEMRA",
  "KLRG1"     = "Effector/TEMRA",
  "KIR3DL1"   = "Effector/TEMRA",
  # NK-like / Innate
  "NCAM1"     = "NK-like/Innate",
  "FCGR3A"    = "NK-like/Innate",
  "SIGLEC7"   = "NK-like/Innate",
  # Exhaustion
  "TIGIT"     = "Exhaustion",
  "PDCD1"     = "Exhaustion",
  "LAG3"      = "Exhaustion"
)

rna_marker_groups <- c(
  # Naïve / Stemness
  "CCR7"    = "Naïve/Stemness",
  "SELL"    = "Naïve/Stemness",
  "TCF7"    = "Naïve/Stemness",
  "LEF1"    = "Naïve/Stemness",
  # γδ TCR
  "TRDV1"   = "γδ TCR",
  "TRGV9"   = "γδ TCR",
  "TRDC"    = "γδ TCR",
  # Memory / Survival
  "IL7R"    = "Memory/Stemness",
  "BCL2"    = "Memory/Stemness",
  "BACH2"   = "Memory/Stemness",
  # Cytotoxic / Effector
  "GZMK"    = "Effector/TEMRA",
  "GZMB"    = "Effector/TEMRA",
  "GNLY"    = "Effector/TEMRA",
  "PRF1"    = "Effector/TEMRA",
  "NKG7"    = "Effector/TEMRA",
  # Effector Transcription Factors
  "TBX21"   = "Effector TFs",
  "EOMES"   = "Effector TFs",
  "RUNX3"   = "Effector TFs",
  # Exhaustion
  "TOX"     = "Exhaustion",
  "PDCD1"   = "Exhaustion",
  "TIGIT"   = "Exhaustion",
  "HAVCR2"  = "Exhaustion",
  # NK / Innate
  "TYROBP"  = "NK-like/Innate",
  "KLRD1"   = "NK-like/Innate",
  "KIR2DL3" = "NK-like/Innate",
  # Activation / Proliferation
  "MKI67"   = "Activation",
  "CXCR6"   = "Activation",
  "HLA-DRA" = "Activation"
)

# ── Shared group colors (consistent across ADT and RNA panels) ─────────────────
group_colors <- c(
  "Naïve"           = "#74C2E1",
  "Naïve/Stemness"  = "#74C2E1",
  "TCR Identity"    = "#9C6FD6",
  "γδ TCR"          = "#9C6FD6",
  "Memory/Stemness" = "#66BB6A",
  "Co-stimulation"  = "#AED581",
  "Effector/TEMRA"  = "#EF5350",
  "Effector TFs"    = "#FFA726",
  "NK-like/Innate"  = "#F06292",
  "Exhaustion"      = "#A1887F",
  "Activation"      = "#FFD54F"
)

# ── Compute per-cluster average expression ────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"
avg_adt <- AverageExpression(
  TARA_cd8,
  assays   = "ADT",
  features = names(adt_marker_groups),
  group.by = "Manual_Annotation_refined",
  slot     = "data"
)$ADT

DefaultAssay(TARA_cd8) <- "RNA"
avg_rna <- AverageExpression(
  TARA_cd8,
  assays   = "RNA",
  features = names(rna_marker_groups),
  group.by = "Manual_Annotation_refined",
  slot     = "data"
)$RNA

# ── Min-max scale each marker (row) to [0, 1] ─────────────────────────────────
scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)))
  colnames(scaled) <- colnames(mat)
  scaled
}

avg_adt_scaled <- scale_01(log10(avg_adt + 1))
avg_rna_scaled <- scale_01(avg_rna)

# Strip any leading "g " prefix that Seurat may prepend to colnames
colnames(avg_adt_scaled) <- gsub("^g ", "", colnames(avg_adt_scaled))
colnames(avg_rna_scaled) <- gsub("^g ", "", colnames(avg_rna_scaled))

# ── Reorder columns by cluster number (1, 6, 8, 9, 27) ───────────────────────
# Raw Seurat colnames match the Manual_Annotation_refined levels without prefix
col_order_raw <- c(
  "Memory CD8 T cell_1",
  "Naïve CD8 T cell_6",
  "CTL-like_8",
  "TRDV1+ γδ T cell_9",
  "GZMK+ CD8 T cell_27"
)

# Column labels — two-line: name (size 11) over cluster tag (size 9).
# pheatmap cannot mix font sizes within one label, so we:
#   1. Pass a plain \n-joined string so pheatmap sizes the column space correctly.
#   2. After pheatmap renders, replace its column label grobs with custom
#      richtext-style grobs that render the two lines at different sizes.
col_name_top <- c("Tcm/Tscm CD8", "Naïve CD8",   "TEMRA/CTL",
                  "TRDV1+ γδ",    "Tpex CD8")
col_name_sub <- c("(Cluster 1)",  "(Cluster 6)",  "(Cluster 8)",
                  "(Cluster 9)",  "(Cluster 27)")

col_display_labels <- paste0(col_name_top, "\n", col_name_sub)

# Helper: replace pheatmap column label row with two-font-size text grobs.
# `ph`        — pheatmap object (contains $gtable)
# `top_strs`  — character vector, top line per column
# `sub_strs`  — character vector, bottom line per column
# `top_size`  — font size for top line (pt)
# `sub_size`  — font size for bottom line (pt)
replace_col_labels <- function(ph, top_strs, sub_strs,
                               top_size = 11, sub_size = 9) {
  gt <- ph$gtable
  
  # Identify the column-label row: the grob named "col_names"
  col_label_idx <- which(gt$layout$name == "col_names")
  if (length(col_label_idx) == 0) return(ph)   # safety: nothing to replace
  
  # Build a replacement grob for each column label
  n <- length(top_strs)
  label_grobs <- lapply(seq_len(n), function(i) {
    grid::grobTree(
      grid::textGrob(
        top_strs[i],
        x = 0.5, y = 0.65,
        just = "centre",
        gp = grid::gpar(fontsize = top_size, fontface = "plain")
      ),
      grid::textGrob(
        sub_strs[i],
        x = 0.5, y = 0.25,
        just = "centre",
        gp = grid::gpar(fontsize = sub_size, fontface = "plain",
                        col = "grey40")
      )
    )
  })
  
  # The col_names grob is typically a gTree with one text grob per column.
  # Replace it with our custom gTree.
  gt$grobs[[col_label_idx]] <- grid::grobTree(
    children = do.call(grid::gList, label_grobs)
  )
  
  ph$gtable <- gt
  ph
}


avg_adt_scaled <- avg_adt_scaled[, col_order_raw]
colnames(avg_adt_scaled) <- col_display_labels

avg_rna_scaled <- avg_rna_scaled[, col_order_raw]
colnames(avg_rna_scaled) <- col_display_labels

# ── Row annotations ────────────────────────────────────────────────────────────
adt_annot <- data.frame(
  Group = factor(adt_marker_groups[rownames(avg_adt_scaled)],
                 levels = names(group_colors)),
  row.names = rownames(avg_adt_scaled)
)

rna_annot <- data.frame(
  Group = factor(rna_marker_groups[rownames(avg_rna_scaled)],
                 levels = names(group_colors)),
  row.names = rownames(avg_rna_scaled)
)

# ── Cluster rows within each group (Ward D2), keep group order fixed ───────────
cluster_within_groups <- function(mat, annot_df) {
  groups    <- levels(droplevels(annot_df$Group))
  row_order <- character(0)
  for (grp in groups) {
    rows <- rownames(annot_df)[annot_df$Group == grp]
    if (length(rows) == 1) {
      row_order <- c(row_order, rows)
    } else {
      hc        <- hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2")
      row_order <- c(row_order, rows[hc$order])
    }
  }
  row_order
}

adt_row_order <- cluster_within_groups(avg_adt_scaled, adt_annot)
rna_row_order <- cluster_within_groups(avg_rna_scaled, rna_annot)

avg_adt_ordered   <- avg_adt_scaled[adt_row_order, ]
avg_rna_ordered   <- avg_rna_scaled[rna_row_order, ]
adt_annot_ordered <- adt_annot[adt_row_order, , drop = FALSE]
rna_annot_ordered <- rna_annot[rna_row_order, , drop = FALSE]

# ── Gap positions between groups ──────────────────────────────────────────────
get_gaps <- function(annot_df) {
  grps <- as.character(annot_df$Group)
  which(grps[-length(grps)] != grps[-1])
}

adt_gaps <- get_gaps(adt_annot_ordered)
rna_gaps <- get_gaps(rna_annot_ordered)

# ── Annotation color lists ────────────────────────────────────────────────────
annot_colors_adt <- list(Group = group_colors[levels(droplevels(adt_annot_ordered$Group))])
annot_colors_rna <- list(Group = group_colors[levels(droplevels(rna_annot_ordered$Group))])

heatmap_colors <- colorRampPalette(c(
  "#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"
))(100)

# ── Draw heatmaps ──────────────────────────────────────────────────────────────
# Column labels are kept horizontal (angle_col = 0).
# The "Group" annotation name is suppressed (annotation_names_row = FALSE) since
# the colored sidebar + legend already conveys this — this was the source of the
# bottom-left overlap with the horizontal column labels.
p_adt <- pheatmap(
  avg_adt_ordered,
  cluster_rows         = FALSE,
  cluster_cols         = FALSE,
  scale                = "none",
  color                = heatmap_colors,
  border_color         = "white",
  annotation_row       = adt_annot_ordered,
  annotation_colors    = annot_colors_adt,
  annotation_names_row = FALSE,  # removes "Group" label → fixes bottom-left overlap
  gaps_row             = adt_gaps,
  cellwidth            = 70,
  cellheight           = 16,
  fontsize             = 11,
  fontsize_row         = 10,
  fontsize_col         = 11,
  angle_col            = 0,      # horizontal labels
  main                 = "Protein (ADT)",
  silent               = TRUE
)
p_rna <- pheatmap(
  avg_rna_ordered,
  cluster_rows         = FALSE,
  cluster_cols         = FALSE,
  scale                = "none",
  color                = heatmap_colors,
  border_color         = "white",
  annotation_row       = rna_annot_ordered,
  annotation_colors    = annot_colors_rna,
  annotation_names_row = FALSE,  # removes "Group" label → fixes bottom-left overlap
  gaps_row             = rna_gaps,
  cellwidth            = 70,
  cellheight           = 16,
  fontsize             = 11,
  fontsize_row         = 10,
  fontsize_col         = 11,
  angle_col            = 0,
  main                 = "mRNA (RNA)",
  silent               = TRUE
)



# ── Save: combined panel + individual panels ──────────────────────────────────
png(file.path(out_heat, "Fig1C_ADT_RNA_Heatmap_CD8_Combined.png"),
    width = 22, height = 15, units = "in", res = 300, bg = "white")
grid.newpage()
grid.text("CD8 T Cell Cluster Validation: Protein (ADT) & mRNA",
          x = 0.5, y = 0.975,
          gp = gpar(fontsize = 18, fontface = "bold"))
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_adt$gtable); popViewport()
pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_rna$gtable); popViewport()
dev.off()

png(file.path(out_heat, "Fig1C_ADT_Heatmap_CD8.png"),
    width = 11, height = 10, units = "in", res = 300, bg = "white")
grid.draw(p_adt$gtable)
dev.off()

png(file.path(out_heat, "Fig1C_RNA_Heatmap_CD8.png"),
    width = 11, height = 13, units = "in", res = 300, bg = "white")
grid.draw(p_rna$gtable)
dev.off()

message("✓ Fig 1C heatmaps saved")


################################################################################
# FIG 1D — CD8 Cluster Proportions (% of total PBMC): HUU vs HEU vs HEI
#
# All 5 CD8 clusters shown, ordered by cluster number (1, 6, 8, 9, 27)
# Wilcoxon test brackets shown only for significant pairs
################################################################################

meta_sub <- as.data.frame(TARA_ALL@meta.data) %>%
  filter(as.character(Timepoint_Group) %in% keep_groups) %>%
  mutate(
    Condition = factor(case_when(
      Timepoint_Group == "PreART_Entry" ~ "HEI",
      TRUE ~ as.character(Timepoint_Group)
    ), levels = c("HUU", "HEU", "HEI"))
  )

# ── Per-sample % of total PBMC for each CD8 cluster ──────────────────────────
cluster_freq <- meta_sub %>%
  group_by(orig.ident, Condition, Manual_Annotation_refined) %>%
  summarise(cluster_cells = dplyr::n(), .groups = "drop") %>%
  # Fill 0 for samples with no cells in a given cluster
  group_by(orig.ident, Condition) %>%
  tidyr::complete(
    Manual_Annotation_refined = cd8_cluster_names,
    fill = list(cluster_cells = 0)
  ) %>%
  ungroup() %>%
  # Denominator = all PBMC per sample (before CD8 subsetting)
  group_by(orig.ident, Condition) %>%
  mutate(total_cells = sum(cluster_cells)) %>%
  ungroup() %>%
  mutate(percent = 100 * cluster_cells / total_cells) %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  mutate(
    Facet_Label = factor(
      cd8_short_labels[as.character(Manual_Annotation_refined)],
      levels = cd8_short_labels   # preserves cluster-number order
    )
  )

# ── Wilcoxon pairwise stats per facet ─────────────────────────────────────────
facet_max <- cluster_freq %>%
  group_by(Facet_Label) %>%
  summarise(max_percent = max(percent, na.rm = TRUE), .groups = "drop")

stat_df <- cluster_freq %>%
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
      p <= 0.0001 ~ "****",
      p <= 0.001  ~ "***",
      p <= 0.01   ~ "**",
      p <= 0.05   ~ "*",
      TRUE        ~ "ns"
    )
  ) %>%
  ungroup() %>%
  filter(p.signif != "ns")

# ── Plot ───────────────────────────────────────────────────────────────────────
p_prop <- ggplot(
  cluster_freq,
  aes(x = Condition, y = percent, fill = Condition, color = Condition)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 2.3, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", nrow = 1) +
  stat_pvalue_manual(
    stat_df,
    label        = "p.signif",
    tip.length   = 0.01,
    bracket.size = 0.6,
    size         = 4.5
  ) +
  scale_fill_manual(values  = cond_cols) +
  scale_color_manual(values = cond_cols) +
  labs(x = NULL, y = "% of total PBMC") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(face = "bold", size = 11),
    axis.text.y     = element_text(size = 11),
    axis.title.y    = element_text(face = "bold", size = 13)
  )

ggsave(
  file.path(out_prop, "Fig1D_CD8_Cluster_Proportions.png"),
  p_prop,
  width  = 18,   # wide to accommodate 5 facets
  height = 6,
  dpi    = 600,
  bg     = "white"
)

message("✓ Fig 1D proportions saved")


################################################################################
# FIG 1E — CD8 Cluster % vs HIV Viral Load (HEI only, Spearman correlation)
#
# All 5 CD8 clusters, ordered by cluster number (1, 6, 8, 9, 27)
# Points colored by viral load (low→blue, high→red) via scale_color_gradientn
# Stat labels pinned to top-left of each facet using Inf/-Inf coordinates —
# this is immune to axis scale differences across facets (no more random placement)
################################################################################

# ── PBMC metadata: HEI patients with known viral load ─────────────────────────
meta_pbmc <- as.data.frame(TARA_ALL@meta.data) %>%
  filter(
    Timepoint_Group == "PreART_Entry",
    !is.na(Viral_Load_num),
    Viral_Load_num > 0
  )

# Denominator: total PBMC per sample
sample_totals <- meta_pbmc %>%
  group_by(orig.ident) %>%
  summarise(total_pbmc = dplyr::n(), .groups = "drop")

# One viral load value per sample (all cells share the same metadata value)
sample_vl <- meta_pbmc %>%
  group_by(orig.ident) %>%
  summarise(Viral_Load_num = dplyr::first(Viral_Load_num), .groups = "drop") %>%
  mutate(log10_VL = log10(Viral_Load_num))

# ── Numerator: CD8 cluster cell counts per sample ─────────────────────────────
meta_cd8_hei <- as.data.frame(TARA_cd8@meta.data) %>%
  filter(orig.ident %in% sample_totals$orig.ident)

cluster_vl_df <- meta_cd8_hei %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  group_by(orig.ident, Manual_Annotation_refined) %>%
  summarise(cluster_cells = dplyr::n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  tidyr::complete(
    Manual_Annotation_refined = cd8_cluster_names,
    fill = list(cluster_cells = 0)
  ) %>%
  ungroup() %>%
  left_join(sample_totals, by = "orig.ident") %>%
  left_join(sample_vl,     by = "orig.ident") %>%
  mutate(
    percent_pbmc = 100 * cluster_cells / total_pbmc,
    Facet_Label  = factor(
      cd8_short_labels[as.character(Manual_Annotation_refined)],
      levels = cd8_short_labels
    )
  )

# ── Spearman correlation per cluster ──────────────────────────────────────────
cor_stats <- cluster_vl_df %>%
  group_by(Facet_Label) %>%
  cor_test(log10_VL, percent_pbmc, method = "spearman") %>%
  mutate(
    stat_label = paste0(
      "Spearman's \u03c1 = ", sprintf("%.2f", cor),
      "\np-value = ",         signif(p, 2)
    )
  ) %>%
  ungroup()

# Stat labels anchored to top-left of every facet via Inf/-Inf coordinates.
# ggplot maps -Inf → axis minimum and Inf → axis maximum after scale expansion,
# so the label always lands in the same corner regardless of free scales.
cor_stats_plot <- cor_stats %>%
  mutate(
    x = -Inf,
    y =  Inf
  )

# ── Plot ───────────────────────────────────────────────────────────────────────
p_vl <- ggplot(cluster_vl_df, aes(x = log10_VL, y = percent_pbmc)) +
  
  # Trend line drawn first so points sit on top
  geom_smooth(
    method    = "lm",
    se        = TRUE,
    linewidth = 0.8,
    color     = "grey40",
    fill      = "grey85",
    alpha     = 0.3
  ) +
  
  geom_point(
    shape  = 21,
    size   = 3.2,
    stroke = 0.4,
    fill   = "#C44E52",  # single clean color; VL already encoded on x-axis
    color  = "white"     # thin white halo separates overlapping points
  ) +
  
  # Stat label: hjust=0 + vjust=1 + x=-Inf + y=Inf → always top-left corner
  geom_text(
    data        = cor_stats_plot,
    aes(x = x, y = y, label = stat_label),
    inherit.aes = FALSE,
    hjust       = -0.05,
    vjust       = 1.2,
    size        = 5,     # clearly readable at publication size
    lineheight  = 1.4,
    color       = "grey15",
    fontface='bold'
  ) +
  
  facet_wrap(~ Facet_Label, scales = "free", nrow = 1) +
  
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  
  labs(
    x = expression(log[10]~"viral load"),
    y = "% of total PBMC"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 11,
                                    margin = margin(4, 0, 4, 0)),
    axis.text        = element_text(size = 10, color = "black"),
    axis.title.x     = element_text(face = "bold", size = 12, margin = margin(t = 6)),
    axis.title.y     = element_text(face = "bold", size = 12, margin = margin(r = 6)),
    legend.position  = "none",
    panel.spacing    = unit(0.8, "lines"),
    plot.margin      = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_vl, "Fig1E_CD8_Clusters_ViralLoad_Correlation.png"),
  p_vl,
  width  = 18,
  height = 5.5,
  dpi    = 600,
  bg     = "white"
)

message("✓ Fig 1E viral load correlations saved")

# Print correlation summary to console
print(cor_stats[, c("Facet_Label", "cor", "statistic", "p")])
