################################################################################
# FIGURE 4 PLOTS: CD8 sub-cluster landscape & ART-dependent functional divergence
#
# Main Panels:
#   4A — UMAP annotated CD8 sub-clusters (13 populations)
#   4B — UMAP clonal expansion overlay
#   4C — UMAP colored by ART status
#   4D — % clonally expanded per sample
#   4E — Volcano (tighter axes, directional arrows, no title)
#   4F — Cohen's d: all 3 pairwise comparisons, all effector clusters
#
# Supplementary Panels (saved to Supplementary/ subfolder):
#   S_ADT_gating    — FAS vs CD45RO scatter plots showing Tscm/Intermediate gating
#   S_annotation    — ADT + RNA heatmap validating 13-population annotation
#   S_cluster_props — Cluster proportions by ART status and by patient
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(patchwork)
library(qs2)
library(scCustomize)
library(SeuratExtend)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(grid)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
manuscript   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 4-5"
fig4_dir     <- file.path(manuscript, "Figure4_panels")
supp_dir     <- file.path(manuscript, "Supplementary_panels")
analysis_dir <- file.path(manuscript, "analysis")

dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load annotated object ────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells,", length(unique(TARA_cd8$CD8_Annotation)), "clusters\n")
cat("Annotations present:\n")
print(sort(table(TARA_cd8$CD8_Annotation), decreasing = TRUE))

# ── Shared aesthetics (CONSISTENT across Fig 4, Fig 5, Supplementary) ────────

# ART status colors — these 3 colors are RESERVED for conditions only
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

art_labels <- c(
  "PreART_Entry"           = "Pre-ART",
  "PostART_Suppressed"     = "Suppressed",
  "PostART_Unsuppressed"   = "Unsuppressed"
)

# Comparison colors — DISTINCT from art_colors to avoid confusion
comp_colors <- c(
  "Suppressed vs Unsuppressed" = "#D4A03C",
  "Pre-ART vs Unsuppressed"    = "#8B5CF6",
  "Suppressed vs Pre-ART"      = "#EC4899"
)

# Cluster colors — 13 populations, grouped by lineage
cluster_colors <- c(
  # Naïve lineage (blues/teals)
  "Naïve CD8"              = "#4E79A7",
  "Naïve CD8 2"            = "#76B7B2",
  "Naïve CD8 3"            = "#A0CBE8",
  "Tscm CD8"               = "#59A14F",
  "Naïve Intermediate CD8" = "#8CD17D",
  # Effector lineage (warm)
  "Transitional Tem CD8"   = "#FF9DA7",
  "TEM CD8"                = "#F28E2B",
  "TEMRA CD8"              = "#EDC948",
  # Innate-like
  "KIR+ innate-like CD8"   = "#BAB0AC",
  "MAIT-like Trm"          = "#9C755F",
  # γδ
  "γδ1 T cell"             = "#B07AA1",
  "Naïve γδ1 T cell"       = "#D4A6C8",
  "γδ2 T cell"             = "#E15759"
)

# Canonical cluster order
col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
  "Tscm CD8", "Naïve Intermediate CD8",
  "Transitional Tem CD8", "TEM CD8", "TEMRA CD8",
  "KIR+ innate-like CD8", "MAIT-like Trm",
  "γδ1 T cell", "Naïve γδ1 T cell", "γδ2 T cell"
)

# Short labels for Figure 4A (remove "CD8" and "T cell")
short_labels <- c(
  "Naïve CD8"              = "Naïve",
  "Naïve CD8 2"            = "Naïve 2",
  "Naïve CD8 3"            = "Naïve 3",
  "Tscm CD8"               = "Tscm",
  "Naïve Intermediate CD8" = "Naïve Intermediate",
  "Transitional Tem CD8"   = "Transitional Tem",
  "TEM CD8"                = "TEM",
  "TEMRA CD8"              = "TEMRA",
  "KIR+ innate-like CD8"   = "KIR+ innate-like",
  "MAIT-like Trm"          = "MAIT-like Trm",
  "γδ1 T cell"             = "γδ1",
  "Naïve γδ1 T cell"       = "Naïve γδ1",
  "γδ2 T cell"             = "γδ2"
)

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)
clone_size_levels <- c("Single (0 < X <= 1)", "Small (1 < X <= 5)",
                       "Medium (5 < X <= 20)", "Large (20 < X <= 100)",
                       "Hyperexpanded (100 < X <= 500)")
clone_size_colors <- setNames(colorblind_vector[c(1, 3, 4, 5, 7)], clone_size_levels)

effector_clusters <- c("TEM CD8", "TEMRA CD8", "Transitional Tem CD8")

module_labels <- c(
  "Exhaustion" = "Exhaustion", "Stemness" = "Stemness / naïve",
  "Cytotoxicity" = "Cytotoxicity", "Stress_IFN" = "Stress / IFN (acute)",
  "TypeI_IFN_Memory" = "Type I IFN memory",
  "Inflammatory_Chemokines" = "Inflammatory chemokines",
  "Terminal_Differentiation" = "Terminal differentiation"
)

# Set cluster order as factor for consistent plotting
TARA_cd8$CD8_Annotation <- factor(TARA_cd8$CD8_Annotation,
                                  levels = col_order_cd8[col_order_cd8 %in% unique(TARA_cd8$CD8_Annotation)])

if (!"PID" %in% colnames(TARA_cd8@meta.data)) {
  pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
  names(pid_vec) <- colnames(TARA_cd8)
  TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")
}

################################################################################
# Panel A: UMAP sub-clusters — Pure ggplot2 with labels and arrows
################################################################################
cat("  Panel A: UMAP sub-clusters (13 populations)...\n")

# Extract UMAP coordinates
umap_coords_A <- as.data.frame(Embeddings(TARA_cd8, reduction = "wnn.umap"))
colnames(umap_coords_A) <- c("UMAP1", "UMAP2")
umap_coords_A$Cluster <- TARA_cd8$CD8_Annotation

# Calculate cluster centroids for labels
centroids_A <- umap_coords_A %>%
  group_by(Cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop") %>%
  mutate(Short_Label = short_labels[as.character(Cluster)]) %>%
  # Nudge TEMRA down slightly to separate from Transitional Tem
  
  mutate(UMAP2 = ifelse(Cluster == "TEMRA CD8", UMAP2 - 0.8, UMAP2))

# Calculate data ranges for proportional arrows
x_range <- diff(range(umap_coords_A$UMAP1, na.rm = TRUE))
y_range <- diff(range(umap_coords_A$UMAP2, na.rm = TRUE))

# Arrow lengths proportional to ranges (same visual length on 10x10 output)
arrow_frac <- 0.12  # 12% of each axis range
x_arrow_len <- x_range * arrow_frac
y_arrow_len <- y_range * arrow_frac

# Arrow position - outside data
x_arrow <- min(umap_coords_A$UMAP1, na.rm = TRUE) - x_range * 0.05
y_arrow <- min(umap_coords_A$UMAP2, na.rm = TRUE) - y_range * 0.08

p_A <- ggplot(umap_coords_A, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = cluster_colors) +
  # Labels: WHITE fill, COLORED text, using SHORT labels
  geom_label_repel(
    data = centroids_A,
    aes(x = UMAP1, y = UMAP2, label = Short_Label, color = Cluster),
    fill = "white",
    size = 10,
    fontface = "bold",
    label.size = 0.6,
    label.padding = unit(0.35, "lines"),
    box.padding = 1.5,
    point.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.6,
    min.segment.length = 0,
    max.overlaps = 30,
    force = 5,
    force_pull = 0.2,
    seed = 42,
    show.legend = FALSE
  ) +
  # UMAP arrows - proportional lengths for equal visual size
  annotate("segment", 
           x = x_arrow, xend = x_arrow + x_arrow_len,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, color = "black") +
  annotate("text", x = x_arrow + x_arrow_len/2, y = y_arrow - y_range * 0.04,
           label = "UMAP1", size = 5, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + y_arrow_len,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, color = "black") +
  annotate("text", x = x_arrow - x_range * 0.04, y = y_arrow + y_arrow_len/2,
           label = "UMAP2", size = 5, fontface = "bold", angle = 90) +
  ggtitle("CD8 T cell sub-clusters") +
  theme_void() +
  theme(
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 25, 25),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig4_dir, "Fig4A_UMAP_CD8_subclusters.png"),
       plot = p_A, width = 10, height = 10, dpi = 300, bg = "white")

################################################################################
# Panel B: UMAP cloneSize — Pure ggplot2 with arrows, NO legend (separate saved)
################################################################################
cat("  Panel B: UMAP cloneSize...\n")

umap_coords_B <- umap_coords_A
umap_coords_B$cloneSize <- TARA_cd8$cloneSize

# Order by clone size so larger clones plot on top
umap_coords_B$cloneSize <- factor(umap_coords_B$cloneSize, levels = clone_size_levels)
umap_ordered_B <- umap_coords_B %>%
  arrange(cloneSize)

p_B <- ggplot(umap_ordered_B, aes(x = UMAP1, y = UMAP2, color = cloneSize)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = clone_size_colors, na.value = "grey80") +
  # UMAP arrows - proportional lengths (same as A)
  annotate("segment", 
           x = x_arrow, xend = x_arrow + x_arrow_len,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, color = "black") +
  annotate("text", x = x_arrow + x_arrow_len/2, y = y_arrow - y_range * 0.04,
           label = "UMAP1", size = 5, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + y_arrow_len,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, color = "black") +
  annotate("text", x = x_arrow - x_range * 0.04, y = y_arrow + y_arrow_len/2,
           label = "UMAP2", size = 5, fontface = "bold", angle = 90) +
  ggtitle("Clonal expansion") +
  theme_void() +
  theme(
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 25, 25),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig4_dir, "Fig4B_UMAP_cloneSize.png"),
       plot = p_B, width = 10, height = 10, dpi = 300, bg = "white")

# Save separate legend for B (horizontal, for placing below plot)
p_B_legend <- ggplot(umap_ordered_B, aes(x = UMAP1, y = UMAP2, color = cloneSize)) +
  geom_point(size = 5) +
  scale_color_manual(values = clone_size_colors, name = "Clone Size", na.value = "grey80",
                     labels = c("Single", "Small", "Medium", "Large", "Hyperexpanded")) +
  theme_void() +
  theme(
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 20),
    legend.title      = element_text(size = 22, face = "bold"),
    legend.key.size   = unit(1.2, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 7), nrow = 1))

legend_B <- cowplot::get_legend(p_B_legend)
ggsave(file.path(fig4_dir, "Fig4B_legend_cloneSize.png"),
       plot = as_ggplot(legend_B), width = 14, height = 2.5, dpi = 300, bg = "white")

################################################################################
# Panel C: UMAP by ART status — Pure ggplot2 with arrows, NO legend (separate saved)
################################################################################
cat("  Panel C: UMAP ART status...\n")

umap_coords_C <- umap_coords_A
umap_coords_C$ART_Status <- TARA_cd8$Timepoint_Group

p_C <- ggplot(umap_coords_C, aes(x = UMAP1, y = UMAP2, color = ART_Status)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = art_colors) +
  # UMAP arrows - proportional lengths (same as A)
  annotate("segment", 
           x = x_arrow, xend = x_arrow + x_arrow_len,
           y = y_arrow, yend = y_arrow,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, color = "black") +
  annotate("text", x = x_arrow + x_arrow_len/2, y = y_arrow - y_range * 0.04,
           label = "UMAP1", size = 5, fontface = "bold") +
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow, yend = y_arrow + y_arrow_len,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, color = "black") +
  annotate("text", x = x_arrow - x_range * 0.04, y = y_arrow + y_arrow_len/2,
           label = "UMAP2", size = 5, fontface = "bold", angle = 90) +
  ggtitle("ART status") +
  theme_void() +
  theme(
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 25, 25),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig4_dir, "Fig4C_UMAP_ART_status.png"),
       plot = p_C, width = 10, height = 10, dpi = 300, bg = "white")

# Save separate legend for C (horizontal, for placing below plot)
p_C_legend <- ggplot(umap_coords_C, aes(x = UMAP1, y = UMAP2, color = ART_Status)) +
  geom_point(size = 5) +
  scale_color_manual(values = art_colors, name = "ART Status", labels = art_labels) +
  theme_void() +
  theme(
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 20),
    legend.title      = element_text(size = 22, face = "bold"),
    legend.key.size   = unit(1.2, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 7, alpha = 1), nrow = 1))

legend_C <- cowplot::get_legend(p_C_legend)
ggsave(file.path(fig4_dir, "Fig4C_legend_ART_status.png"),
       plot = as_ggplot(legend_C), width = 12, height = 2, dpi = 300, bg = "white")

################################################################################
# Panel D: % expanding per sample
################################################################################
cat("  Panel D: % expanding per sample...\n")

total_per_sample <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters) %>%
  group_by(orig.ident, PID, Timepoint_Group) %>%
  summarise(total_cells = n(), .groups = "drop")

expand_per_sample <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% effector_clusters &
           !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)") %>%
  group_by(orig.ident, PID, Timepoint_Group) %>%
  summarise(n_expanding = n(), .groups = "drop")

sample_expansion <- total_per_sample %>%
  left_join(expand_per_sample, by = c("orig.ident", "PID", "Timepoint_Group")) %>%
  mutate(
    n_expanding = replace_na(n_expanding, 0),
    pct_expanding = n_expanding / total_cells * 100
  )

sample_expansion$Timepoint_Group <- factor(sample_expansion$Timepoint_Group,
                                           levels = names(art_colors))

p_D <- ggplot(sample_expansion, aes(x = Timepoint_Group, y = pct_expanding, fill = Timepoint_Group)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(c("PostART_Suppressed", "PostART_Unsuppressed")),
    method = "wilcox.test", label = "p.signif",
    size = 6, step.increase = 0.12, tip.length = 0.01, hide.ns = FALSE
  ) +
  scale_fill_manual(values = art_colors, name = NULL) +
  scale_x_discrete(labels = art_labels) +
  labs(x = NULL, y = "% cells clonally expanded",
       title = "Expansion level per sample") +
  theme_cowplot(font_size = 20) +
  theme(
    axis.text.x      = element_text(size = 14),
    axis.text.y      = element_text(size = 13),
    axis.title.y     = element_text(size = 15),
    plot.title       = element_text(size = 16),
    legend.position  = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(fig4_dir, "Fig4D_PctExpanding_perSample.png"),
       plot = p_D, width = 6, height = 7, dpi = 300, bg = "white")

################################################################################
# Panel E: Volcano — NO TITLE, BIGGER AXIS TEXT, DIRECTIONAL ARROWS (TIGHTENED)
#   EDIT: Arrows repositioned closer to axis; reduced bottom margin
#   EDIT: H3F3A and H3F3B force-labeled regardless of significance
################################################################################
cat("  Panel E: Volcano...\n")

dge_all_path <- file.path(analysis_dir, "03_DGE", "DGE_MAST_Expanding_AllClusters_Sup_vs_Unsup.csv")

if (file.exists(dge_all_path)) {
  dge_all <- read.csv(dge_all_path, row.names = 1)
  
  exhaustion_genes   <- c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1")
  stemness_genes     <- c("TCF7", "SELL", "CCR7", "LEF1", "BACH2", "BCL2", "IL7R", "S1PR1", "KLF2")
  cytotoxic_genes    <- c("GZMB", "GNLY", "PRF1", "NKG7", "GZMA", "GZMH", "FGFBP2")
  stress_genes       <- c("HSPA1A", "HSPA1B", "DNAJB1", "IFI27", "IFI44L")
  ifn_memory_genes   <- c("IFIT1", "IFIT3", "ISG15", "MX1")
  chemokine_genes    <- c("CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5")
  terminal_genes     <- c("ZEB2", "PRDM1", "TBX21", "CX3CR1", "S1PR5", "ID2", "EOMES")
  activation_genes   <- c("HLA-DRA", "CD38", "FAS", "CD69")
  
  highlight_cols <- c(
    "Exhaustion"        = "#A32D2D",  "Stemness"          = "#185FA5",
    "Cytotoxicity"      = "#3B6D11",  "Stress response"   = "#BA7517",
    "Type I IFN memory" = "#0F6E56",  "Chemokines"        = "#993556",
    "Terminal diff."    = "#534AB7",  "Activation"        = "#D85A30",
    "Other"             = "grey80"
  )
  
  dge_all$gene <- rownames(dge_all)
  dge_all$neg_log10_padj <- -log10(dge_all$p_val_adj + 1e-300)
  
  dge_all$highlight <- "Other"
  dge_all$highlight[dge_all$gene %in% exhaustion_genes]  <- "Exhaustion"
  dge_all$highlight[dge_all$gene %in% stemness_genes]    <- "Stemness"
  dge_all$highlight[dge_all$gene %in% cytotoxic_genes]   <- "Cytotoxicity"
  dge_all$highlight[dge_all$gene %in% stress_genes]      <- "Stress response"
  dge_all$highlight[dge_all$gene %in% ifn_memory_genes]  <- "Type I IFN memory"
  dge_all$highlight[dge_all$gene %in% chemokine_genes]   <- "Chemokines"
  dge_all$highlight[dge_all$gene %in% terminal_genes]    <- "Terminal diff."
  dge_all$highlight[dge_all$gene %in% activation_genes]  <- "Activation"
  dge_all$highlight <- factor(dge_all$highlight, levels = names(highlight_cols))
  
  # ── Genes to label: significant category genes + force-labeled histones ────
  genes_to_label <- c(exhaustion_genes, stemness_genes, cytotoxic_genes,
                      stress_genes, ifn_memory_genes, chemokine_genes,
                      terminal_genes, activation_genes)
  genes_to_label <- genes_to_label[genes_to_label %in% dge_all$gene]
  
  dge_all$label <- ""
  dge_all$label[dge_all$gene %in% genes_to_label & dge_all$p_val_adj < 0.05] <-
    dge_all$gene[dge_all$gene %in% genes_to_label & dge_all$p_val_adj < 0.05]
  
  x_vals <- abs(dge_all$avg_log2FC[is.finite(dge_all$avg_log2FC)])
  x_lim  <- quantile(x_vals, 0.995) * 1.15
  dge_all$avg_log2FC_plot <- pmax(pmin(dge_all$avg_log2FC, x_lim * 0.98), -x_lim * 0.98)
  
  y_vals <- dge_all$neg_log10_padj[is.finite(dge_all$neg_log10_padj)]
  y_cap  <- quantile(y_vals, 0.995) * 1.10
  dge_all$neg_log10_padj_plot <- pmin(dge_all$neg_log10_padj, y_cap)
  
  p_E <- ggplot(dge_all, aes(x = avg_log2FC_plot, y = neg_log10_padj_plot, color = highlight)) +
    geom_point(data = dge_all %>% filter(highlight == "Other"), size = 2, alpha = 0.3) +
    geom_point(data = dge_all %>% filter(highlight != "Other"), size = 7, alpha = 0.85) +
    geom_label_repel(
      data = dge_all %>% filter(label != ""),
      aes(label = label, fill = highlight), color = "white",
      size = 8, fontface = "bold.italic",
      max.overlaps = 50, segment.size = 0.4, segment.color = "grey40",
      min.segment.length = 0.1, box.padding = 0.5,
      label.size = 0.2, show.legend = FALSE,
      force = 2, force_pull = 0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = highlight_cols, name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = highlight_cols, guide = "none") +
    scale_x_continuous(limits = c(-x_lim, x_lim), oob = scales::squish) +
    scale_y_continuous(limits = c(0, y_cap), oob = scales::squish) +
    labs(x = expression(log[2]~fold~change),
         y = expression(-log[10]~adjusted~italic(p))) +
    theme_cowplot(font_size = 20) +
    theme(
      plot.title       = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size = 22),
      legend.title     = element_text(size = 23, face = "bold"),
      legend.key.size  = unit(1.5, "cm"),
      legend.spacing.y = unit(0.3, "cm"),
      axis.text        = element_text(size = 18),
      axis.title       = element_text(size = 20),
      axis.title.x     = element_text(size = 18),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 9, alpha = 1)))
  
  # ── Directional arrows — repositioned tighter beneath the x-axis ──────────
  arrow_right <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(1.0, "cm"), type = "closed"),
              gp = gpar(col = "#52B788", lwd = 10)),
    textGrob("Higher in suppressed", x = unit(0.74, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#52B788", fontsize = 28, fontface = "bold"))
  )
  arrow_left <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(1.0, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 10)),
    textGrob("Higher in unsuppressed", x = unit(0.26, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 28, fontface = "bold"))
  )
  
  p_E_final <- p_E +
    theme(plot.margin = margin(10, 10, 80, 10)) +
    annotation_custom(grob = arrow_right, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.24, ymax = -y_cap * 0.10) +
    annotation_custom(grob = arrow_left, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.24, ymax = -y_cap * 0.10) +
    coord_cartesian(ylim = c(0, y_cap), xlim = c(-x_lim, x_lim), clip = "off")
  
  ggsave(file.path(fig4_dir, "Fig4E_Volcano_AllClusters.png"),
         plot = p_E_final, width = 14, height = 11, dpi = 300, bg = "white")
} else {
  cat("    WARNING: DGE file not found at", dge_all_path, "\n")
}

################################################################################
# Panel F: Cohen's d — ALL 3 PAIRWISE COMPARISONS, ALL EFFECTOR CLUSTERS
################################################################################
cat("  Panel F: Cohen's d (all 3 comparisons)...\n")

module_path <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")

if (file.exists(module_path)) {
  mod_scores <- read.csv(module_path, row.names = 1)
  available_mods <- names(module_labels)[names(module_labels) %in% colnames(mod_scores)]
  
  cat("    Clusters in ModuleScores CSV:\n")
  print(table(mod_scores$CD8_Annotation))
  
  comparisons <- list(
    list(num = "PostART_Suppressed",   denom = "PostART_Unsuppressed", label = "Suppressed vs Unsuppressed"),
    list(num = "PreART_Entry",         denom = "PostART_Unsuppressed", label = "Pre-ART vs Unsuppressed"),
    list(num = "PostART_Suppressed",   denom = "PreART_Entry",         label = "Suppressed vs Pre-ART")
  )
  
  all_cd_results <- list()
  min_cells <- 3
  
  for (comp in comparisons) {
    for (mod in available_mods) {
      for (cl in effector_clusters) {
        cl_data <- mod_scores %>% filter(CD8_Annotation == cl)
        vals_num   <- cl_data %>% filter(Timepoint_Group == comp$num) %>% pull(!!sym(mod))
        vals_denom <- cl_data %>% filter(Timepoint_Group == comp$denom) %>% pull(!!sym(mod))
        
        if (length(vals_num) >= min_cells & length(vals_denom) >= min_cells) {
          pooled_sd <- sqrt(((length(vals_num)-1)*sd(vals_num)^2 +
                               (length(vals_denom)-1)*sd(vals_denom)^2) /
                              (length(vals_num)+length(vals_denom)-2))
          cd <- (mean(vals_num) - mean(vals_denom)) / pooled_sd
          wt <- wilcox.test(vals_num, vals_denom)
          all_cd_results[[paste(comp$label, mod, cl)]] <- data.frame(
            Comparison = comp$label, Module = mod, CD8_Annotation = cl,
            cohens_d = cd, p_value = wt$p.value, stringsAsFactors = FALSE)
        } else {
          cat(sprintf("    SKIPPED: %s | %s | %s (n_num=%d, n_denom=%d)\n",
                      comp$label, cl, mod, length(vals_num), length(vals_denom)))
        }
      }
    }
  }
  
  cd_all <- do.call(rbind, all_cd_results)
  
  cd_all <- cd_all %>%
    group_by(Comparison) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  cd_all$star <- ifelse(cd_all$p_adj < 0.001, "***",
                        ifelse(cd_all$p_adj < 0.01, "**",
                               ifelse(cd_all$p_adj < 0.05, "*", "")))
  
  cd_all$Module_Label <- module_labels[cd_all$Module]
  module_order <- c(
    "Exhaustion", "Stemness / naïve", "Cytotoxicity",
    "Terminal differentiation", "Stress / IFN (acute)",
    "Inflammatory chemokines", "Type I IFN memory"
  )
  cd_all$Module_Label <- factor(cd_all$Module_Label, levels = rev(module_order))
  cd_all$CD8_Annotation <- factor(cd_all$CD8_Annotation,
                                  levels = c("Transitional Tem CD8", "TEM CD8", "TEMRA CD8"))
  cd_all$Comparison <- factor(cd_all$Comparison,
                              levels = c("Suppressed vs Unsuppressed",
                                         "Pre-ART vs Unsuppressed",
                                         "Suppressed vs Pre-ART"))
  
  cd_all <- cd_all %>% filter(!is.na(CD8_Annotation) & !is.na(Comparison) & !is.na(Module_Label))
  
  x_range_all <- max(abs(cd_all$cohens_d), na.rm = TRUE)
  cd_all$star_x <- ifelse(cd_all$cohens_d > 0,
                          cd_all$cohens_d + x_range_all * 0.04,
                          cd_all$cohens_d - x_range_all * 0.04)
  cd_all$star_hjust <- ifelse(cd_all$cohens_d > 0, 0, 1)
  
  p_F <- ggplot(cd_all, aes(x = cohens_d, y = Module_Label, fill = Comparison)) +
    geom_col(width = 0.75, alpha = 0.85,
             position = position_dodge(width = 0.8)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(aes(x = star_x, label = star, hjust = star_hjust, group = Comparison),
              position = position_dodge(width = 0.8),
              size = 4.5, color = "black", vjust = 0.35) +
    facet_wrap(~ CD8_Annotation, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = comp_colors, name = "Comparison", drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = "Cohen's d (effect size)", y = NULL,
         title = "Module score differences across ART conditions") +
    theme_cowplot(font_size = 18) +
    theme(
      strip.text        = element_text(size = 14, face = "bold"),
      strip.background  = element_rect(fill = "grey95", color = NA),
      axis.text.y       = element_text(size = 12),
      axis.text.x       = element_text(size = 12),
      axis.title.x      = element_text(size = 14),
      legend.position   = "bottom",
      legend.text       = element_text(size = 13),
      legend.title      = element_text(size = 14, face = "bold"),
      plot.title        = element_text(size = 17, face = "bold"),
      plot.margin       = margin(10, 30, 10, 10),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(fig4_dir, "Fig4F_CohenD_AllComparisons.png"),
         plot = p_F, width = 18, height = 9, dpi = 300, bg = "white")
  
  write.csv(cd_all %>% select(Comparison, Module, CD8_Annotation, cohens_d, p_value, p_adj, star),
            file.path(fig4_dir, "Fig4F_CohenD_stats.csv"), row.names = FALSE)
} else {
  cat("    WARNING: Module scores file not found at", module_path, "\n")
}

################################################################################
# SUPPLEMENTARY: ADT Gating Validation
################################################################################
cat("\n=== Generating supplementary panels ===\n")

# ── S1: ADT gating scatter plots (FAS vs CD45RO) ────────────────────────────
cat("  Supp: ADT gating validation...\n")

if ("ADT_gate" %in% colnames(TARA_cd8@meta.data) & "source_cluster" %in% colnames(TARA_cd8@meta.data)) {
  
  # Get naïve-lineage cells (those that were gated)
  gated_mask <- !is.na(TARA_cd8$ADT_gate)
  if (sum(gated_mask) > 0) {
    gated_cells <- colnames(TARA_cd8)[gated_mask]
    
    DefaultAssay(TARA_cd8) <- "ADT"
    adt_data <- GetAssayData(TARA_cd8, slot = "data")
    
    gate_plot_df <- data.frame(
      CD45RO    = as.numeric(adt_data["CD45RO", gated_cells]),
      CD45RA    = as.numeric(adt_data["CD45RA", gated_cells]),
      FAS       = as.numeric(adt_data["FAS", gated_cells]),
      gate      = TARA_cd8$CD8_Annotation[gated_cells],
      source    = TARA_cd8$source_cluster[gated_cells],
      stringsAsFactors = FALSE
    )
    
    # Map source clusters to meaningful labels
    source_map <- c("0" = "Cluster 0\n(classical Naïve)",
                    "2" = "Cluster 2\n(KIR+/pre-ART)",
                    "5" = "Cluster 5\n(BACH2+/post-ART)",
                    "7" = "Cluster 7\n(rescued αβ)")
    gate_plot_df$source_label <- source_map[gate_plot_df$source]
    gate_plot_df <- gate_plot_df %>% filter(!is.na(source_label))
    
    gate_cols_supp <- c(
      "Naïve CD8"              = "#4E79A7",
      "Naïve CD8 2"            = "#76B7B2",
      "Naïve CD8 3"            = "#A0CBE8",
      "Tscm CD8"               = "#59A14F",
      "Naïve Intermediate CD8" = "#8CD17D"
    )
    
    p_supp_gate <- ggplot(gate_plot_df, aes(x = FAS, y = CD45RO, color = gate)) +
      geom_point(size = 0.5, alpha = 0.4) +
      scale_color_manual(values = gate_cols_supp, name = "ADT Gate") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      facet_wrap(~ source_label, nrow = 1) +
      labs(title = "ADT sub-gating of naïve clusters",
           subtitle = "Tscm = CD45RO- FAS+  |  Naïve Intermediate = CD45RO+  |  True naïve = CD45RO- FAS-",
           x = "FAS/CD95 (DSB)", y = "CD45RO (DSB)") +
      theme_cowplot(font_size = 12) +
      theme(
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text       = element_text(size = 11),
        legend.position  = "bottom",
        legend.text      = element_text(size = 11),
        plot.background  = element_rect(fill = "white", color = NA)
      ) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
    
    ggsave(file.path(supp_dir, "SuppFig_ADT_gating_validation.png"),
           plot = p_supp_gate, width = 18, height = 6, dpi = 300, bg = "white")
    cat("    ADT gating validation saved.\n")
  }
} else {
  cat("    ADT_gate or source_cluster columns not found — skipping gating plot.\n")
  cat("    (This is expected if running from a pre-gating saved object.)\n")
}

# ── S2: Cluster proportions by ART status ────────────────────────────────────
cat("  Supp: Cluster proportions...\n")

prop_tp <- as.data.frame(prop.table(
  table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 2
))
colnames(prop_tp) <- c("Cluster", "Timepoint_Group", "Proportion")
prop_tp$Cluster <- factor(prop_tp$Cluster, levels = col_order_cd8)

p_supp_prop <- ggplot(prop_tp, aes(x = Timepoint_Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_fill_manual(values = cluster_colors, name = "CD8 Sub-cluster") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_discrete(labels = art_labels) +
  labs(x = NULL, y = "Proportion",
       title = "CD8 sub-cluster composition by ART status") +
  theme_cowplot(font_size = 14) +
  theme(
    axis.text.x      = element_text(size = 12),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.4, "cm"),
    plot.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(ncol = 2))

ggsave(file.path(supp_dir, "SuppFig_ClusterProportions_byART.png"),
       plot = p_supp_prop, width = 10, height = 8, dpi = 300, bg = "white")

# ── S3: Cluster proportions by patient ───────────────────────────────────────
prop_pid <- as.data.frame(prop.table(
  table(TARA_cd8$CD8_Annotation, TARA_cd8$PID), margin = 2
))
colnames(prop_pid) <- c("Cluster", "PID", "Proportion")
prop_pid$Cluster <- factor(prop_pid$Cluster, levels = col_order_cd8)

p_supp_pid <- ggplot(prop_pid, aes(x = PID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_fill_manual(values = cluster_colors, name = "CD8 Sub-cluster") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion",
       title = "CD8 sub-cluster composition by patient") +
  theme_cowplot(font_size = 14) +
  theme(
    axis.text.x      = element_text(size = 11, angle = 35, hjust = 1),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.4, "cm"),
    plot.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(ncol = 2))

ggsave(file.path(supp_dir, "SuppFig_ClusterProportions_byPatient.png"),
       plot = p_supp_pid, width = 12, height = 8, dpi = 300, bg = "white")

cat("\n=== Figure 4 complete ===\n")
cat("  Main panels:", length(list.files(fig4_dir, "\\.png$")), "files in", fig4_dir, "\n")
cat("  Supplementary:", length(list.files(supp_dir, "\\.png$")), "files in", supp_dir, "\n")