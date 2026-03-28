################################################################################
# FIGURE 4 PLOTS: CD8 sub-cluster landscape & ART-dependent functional divergence
#
# Panels:
#   4A — UMAP annotated CD8 sub-clusters (large fonts)
#   4B — UMAP clonal expansion overlay (large legend)
#   4C — UMAP colored by ART status
#   4D — % clonally expanded per sample
#   4E — Volcano (tighter axes, directional arrows, no title)
#   4F — Cohen's d: all 3 pairwise comparisons, all effector clusters
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

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
manuscript   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 4-5"
fig4_dir     <- file.path(manuscript, "Figure4_panels")
analysis_dir <- file.path(manuscript, "analysis")

dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load annotated object ────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells,", length(unique(TARA_cd8$CD8_Annotation)), "clusters\n")

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
# These are used only for pairwise comparison labels in Cohen's d plots
comp_colors <- c(
  "Suppressed vs Unsuppressed" = "#D4A03C",   # gold/amber
  "Pre-ART vs Unsuppressed"    = "#8B5CF6",   # violet
  "Suppressed vs Pre-ART"      = "#EC4899"    # pink
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

if (!"PID" %in% colnames(TARA_cd8@meta.data)) {
  pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
  names(pid_vec) <- colnames(TARA_cd8)
  TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")
}

################################################################################
# Panel A: UMAP sub-clusters — LARGER FONTS
################################################################################
cat("  Panel A: UMAP sub-clusters...\n")

p_A <- DimPlot2(
  TARA_cd8,
  reduction  = "wnn.umap",
  group.by   = "CD8_Annotation",
  cols       = "default",
  label      = TRUE,
  box        = TRUE,
  label.size = 9,
  repel      = TRUE,
  pt.size    = 0.6,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
) +
  ggtitle("CD8 T cell sub-clusters") +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5))

ggsave(file.path(fig4_dir, "Fig4A_UMAP_CD8_subclusters.png"),
       plot = p_A, width = 12, height = 10, dpi = 300, bg = "white")

################################################################################
# Panel B: UMAP cloneSize — LARGER LEGEND
################################################################################
cat("  Panel B: UMAP cloneSize...\n")

p_B <- DimPlot_scCustom(
  TARA_cd8,
  group.by  = "cloneSize",
  reduction = "wnn.umap",
  pt.size   = 0.6
) +
  scale_color_manual(values = clone_size_colors, name = "Clone Size") +
  NoAxes() +
  theme_umap_arrows() +
  ggtitle("Clonal expansion") +
  theme(
    plot.title        = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.text       = element_text(size = 24),
    legend.title      = element_text(size = 26, face = "bold"),
    legend.key.size   = unit(2, "cm"),
    legend.key.height = unit(1.8, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 7)))

ggsave(file.path(fig4_dir, "Fig4B_UMAP_cloneSize.png"),
       plot = p_B, width = 14, height = 10, dpi = 300, bg = "white")

################################################################################
# Panel C: UMAP by ART status — NEW
################################################################################
cat("  Panel C: UMAP ART status...\n")

p_C <- DimPlot_scCustom(
  TARA_cd8,
  group.by   = "Timepoint_Group",
  reduction  = "wnn.umap",
  pt.size    = 0.6,
  colors_use = art_colors
) +
  NoAxes() +
  theme_umap_arrows() +
  ggtitle("ART status") +
  scale_color_manual(values = art_colors, name = "ART Status", labels = art_labels) +
  theme(
    plot.title        = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.text       = element_text(size = 24),
    legend.title      = element_text(size = 26, face = "bold"),
    legend.key.size   = unit(2, "cm"),
    legend.key.height = unit(1.8, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 7, alpha = 1)))

ggsave(file.path(fig4_dir, "Fig4C_UMAP_ART_status.png"),
       plot = p_C, width = 12, height = 10, dpi = 300, bg = "white")

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
# Panel E: Volcano — NO TITLE, BIGGER AXIS TEXT, DIRECTIONAL ARROWS
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
  
  genes_to_label <- c(exhaustion_genes, stemness_genes, cytotoxic_genes,
                      stress_genes, ifn_memory_genes, chemokine_genes,
                      terminal_genes, activation_genes)
  genes_to_label <- genes_to_label[genes_to_label %in% dge_all$gene]
  dge_all$label <- ""
  dge_all$label[dge_all$gene %in% genes_to_label & dge_all$p_val_adj < 0.05] <-
    dge_all$gene[dge_all$gene %in% genes_to_label & dge_all$p_val_adj < 0.05]
  
  # Tight x-axis: cap at 99.5th percentile of |log2FC| (most genes 0–6, outliers at 10–20)
  x_vals <- abs(dge_all$avg_log2FC[is.finite(dge_all$avg_log2FC)])
  x_lim  <- quantile(x_vals, 0.995) * 1.15
  # Clamp outlier genes to the edge so they're visible but don't stretch the axis
  dge_all$avg_log2FC_plot <- pmax(pmin(dge_all$avg_log2FC, x_lim * 0.98), -x_lim * 0.98)
  
  # Tight y-axis: cap at 99.5th percentile
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
      label.size = 0.2, show.legend = FALSE
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = highlight_cols, name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = highlight_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_cap), clip = "off") +
    # ── Directional arrows BELOW the plot area (in margin space) ─────────────
    # Negative y values place these below y=0, visible because clip = "off"
    annotate("segment", x = 0.15, xend = x_lim * 0.85,
             y = -y_cap * 0.12, yend = -y_cap * 0.12,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "#52B788", linewidth = 1.3) +
    annotate("text", x = x_lim * 0.52, y = -y_cap * 0.17,
             label = "Higher in suppressed", color = "#52B788",
             size = 6.5, fontface = "bold") +
    annotate("segment", x = -0.15, xend = -x_lim * 0.85,
             y = -y_cap * 0.12, yend = -y_cap * 0.12,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             color = "#E76F51", linewidth = 1.3) +
    annotate("text", x = -x_lim * 0.52, y = -y_cap * 0.17,
             label = "Higher in unsuppressed", color = "#E76F51",
             size = 6.5, fontface = "bold") +
    labs(
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~italic(p))
    ) +
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
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(10, 10, 75, 10)  # extra bottom margin for arrows below axis title
    ) +
    guides(color = guide_legend(override.aes = list(size = 9, alpha = 1)))
  
  ggsave(file.path(fig4_dir, "Fig4E_Volcano_AllClusters.png"),
         plot = p_E, width = 14, height = 11, dpi = 300, bg = "white")
} else {
  cat("    WARNING: DGE file not found at", dge_all_path, "\n")
}

################################################################################
# Panel F: Cohen's d — ALL 3 PAIRWISE COMPARISONS, ALL EFFECTOR CLUSTERS
#   Comparison colors are DISTINCT from art_colors (gold/violet/pink)
#   so they don't get confused with the condition colors (blue/green/orange)
################################################################################
cat("  Panel F: Cohen's d (all 3 comparisons)...\n")

module_path <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")

if (file.exists(module_path)) {
  mod_scores <- read.csv(module_path, row.names = 1)
  available_mods <- names(module_labels)[names(module_labels) %in% colnames(mod_scores)]
  
  # ── DIAGNOSTIC: what clusters and conditions are in the data? ──────────────
  cat("    Clusters in ModuleScores CSV:\n")
  print(table(mod_scores$CD8_Annotation))
  cat("    Clusters × Timepoint_Group:\n")
  print(table(mod_scores$CD8_Annotation, mod_scores$Timepoint_Group))
  cat("    Looking for effector_clusters:", paste(effector_clusters, collapse = ", "), "\n")
  
  # Check if old names are still in the CSV (prep not re-run with updated annotations)
  old_names_present <- intersect(
    c("Effector CD8", "Terminal TEMRA CD8"),
    unique(mod_scores$CD8_Annotation)
  )
  if (length(old_names_present) > 0) {
    cat("\n    *** WARNING: Old cluster names found in CSV:", paste(old_names_present, collapse = ", "), "\n")
    cat("    *** Remapping old names to new names for this plot...\n\n")
    mod_scores$CD8_Annotation <- dplyr::recode(mod_scores$CD8_Annotation,
                                               "Effector CD8"        = "TEM CD8",
                                               "Terminal TEMRA CD8"  = "TEMRA CD8"
    )
  }
  
  comparisons <- list(
    list(num = "PostART_Suppressed",   denom = "PostART_Unsuppressed", label = "Suppressed vs Unsuppressed"),
    list(num = "PreART_Entry",         denom = "PostART_Unsuppressed", label = "Pre-ART vs Unsuppressed"),
    list(num = "PostART_Suppressed",   denom = "PreART_Entry",         label = "Suppressed vs Pre-ART")
  )
  
  all_cd_results <- list()
  min_cells <- 3  # minimum cells per group (lowered from 5 to capture small clusters)
  
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
  
  # BH-correct within each comparison
  cd_all <- cd_all %>%
    group_by(Comparison) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  cd_all$star <- ifelse(cd_all$p_adj < 0.001, "***",
                        ifelse(cd_all$p_adj < 0.01, "**",
                               ifelse(cd_all$p_adj < 0.05, "*", "")))
  
  cd_all$Module_Label <- module_labels[cd_all$Module]
  # Order modules: key narrative pair first, then effector, then secondary
  module_order <- c(
    "Exhaustion",
    "Stemness / naïve",
    "Cytotoxicity",
    "Terminal differentiation",
    "Stress / IFN (acute)",
    "Inflammatory chemokines",
    "Type I IFN memory"
  )
  cd_all$Module_Label <- factor(cd_all$Module_Label, levels = rev(module_order))
  
  # Order clusters by differentiation: Transitional → TEM → TEMRA
  cd_all$CD8_Annotation <- factor(cd_all$CD8_Annotation,
                                  levels = c("Transitional Tem CD8", "TEM CD8", "TEMRA CD8"))
  
  # Order comparisons: main comparison first, then the other two
  cd_all$Comparison <- factor(cd_all$Comparison,
                              levels = c("Suppressed vs Unsuppressed",
                                         "Pre-ART vs Unsuppressed",
                                         "Suppressed vs Pre-ART"))
  
  # ── Diagnostic: verify data before plotting ─────────────────────────────────
  cat("    Cohen's d data summary:\n")
  cat("    Rows:", nrow(cd_all), "\n")
  cat("    Comparisons:", paste(levels(cd_all$Comparison), collapse = ", "), "\n")
  cat("    Clusters:", paste(levels(cd_all$CD8_Annotation), collapse = ", "), "\n")
  cat("    Non-NA clusters:", sum(!is.na(cd_all$CD8_Annotation)), "/", nrow(cd_all), "\n")
  cat("    Non-NA comparisons:", sum(!is.na(cd_all$Comparison)), "/", nrow(cd_all), "\n")
  
  # Drop any rows where factoring created NAs
  cd_all <- cd_all %>% filter(!is.na(CD8_Annotation) & !is.na(Comparison) & !is.na(Module_Label))
  
  # Star positions
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
  
  # Export stats
  write.csv(cd_all %>% select(Comparison, Module, CD8_Annotation, cohens_d, p_value, p_adj, star),
            file.path(fig4_dir, "Fig4F_CohenD_stats.csv"), row.names = FALSE)
  
  cat("    Clusters with data in F:\n")
  print(table(cd_all$CD8_Annotation, cd_all$Comparison))
} else {
  cat("    WARNING: Module scores file not found at", module_path, "\n")
}

cat("\n=== Figure 4 complete:", length(list.files(fig4_dir, "\\.png$")), "panels ===\n")
cat("Output:", fig4_dir, "\n")