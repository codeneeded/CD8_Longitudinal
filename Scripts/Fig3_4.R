################################################################################
# FIGURE 3 PLOTS: Main Figure + Supplementary Figure
#   — TARA Cohort, HEI-only CD8 sub-clusters
#
# Prerequisite: Run Figure3_Prep.R first to generate:
#   - TARA_cd8_HEI_annotated_final.qs2
#   - All analysis CSVs in analysis/ subfolders
#
# This script ONLY generates publication-ready figure panels.
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(patchwork)
library(pheatmap)
library(grid)
library(qs2)
library(scCustomize)
library(SeuratExtend)
library(ggrepel)
library(ggpubr)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir    <- "~/Documents/CD8_Longitudinal"
saved_dir   <- file.path(base_dir, "saved_R_data")
out_dir     <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 3/"
fig3_dir    <- file.path(out_dir, "Figure3_panels")
fig4_dir    <- file.path(out_dir, "Figure4_panels")
fig_supp    <- file.path(out_dir, "supplementary")
analysis_dir <- file.path(out_dir, "analysis")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp, recursive = TRUE, showWarnings = FALSE)

# ── Load annotated object ────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells,", length(unique(TARA_cd8$CD8_Annotation)), "clusters\n")

# ── Shared aesthetics ────────────────────────────────────────────────────────
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)
clone_size_levels <- c("Single (0 < X <= 1)", "Small (1 < X <= 5)",
                       "Medium (5 < X <= 20)", "Large (20 < X <= 100)",
                       "Hyperexpanded (100 < X <= 500)")
clone_size_colors <- setNames(c(colorblind_vector[c(1, 3, 4, 5, 7)]), clone_size_levels)

################################################################################
# ══════════════════════════════════════════════════════════════════════════════
# MAIN FIGURE — 9 panels (A–I)
# ══════════════════════════════════════════════════════════════════════════════
################################################################################

# ── Panel A: UMAP with annotated CD8 sub-clusters ───────────────────────────
p_A <- DimPlot2(
  TARA_cd8,
  reduction  = "wnn.umap",
  group.by   = "CD8_Annotation",
  cols       = "default",
  label      = TRUE,
  box        = TRUE,
  label.size = 7,
  repel      = TRUE,
  pt.size    = 0.6,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
)

ggsave(file.path(fig3_dir, "Fig3A_UMAP_CD8_subclusters.png"),
       plot = p_A, width = 10, height = 9, dpi = 300, bg = "white")

# ── Panel B: UMAP with cloneSize overlay ─────────────────────────────────────
p_B <- DimPlot_scCustom(
  TARA_cd8,
  group.by  = "cloneSize",
  reduction = "wnn.umap",
  pt.size   = 0.6
) +
  scale_color_manual(values = clone_size_colors, name = "Clone Size") +
  NoAxes() +
  theme(
    legend.text      = element_text(size = 22),
    legend.title     = element_text(size = 23, face = "bold"),
    legend.key.size  = unit(1.5, "cm"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(file.path(fig3_dir, "Fig3B_UMAP_cloneSize.png"),
       plot = p_B, width = 12, height = 9, dpi = 300, bg = "white")

# ── Panel C: % expanding per sample by ART status ───────────────────────────
effector_clusters <- c("Effector CD8", "Terminal TEMRA CD8", "Transitional Tem CD8")

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
  mutate(
    n_expanding = replace_na(n_expanding, 0),
    pct_expanding = n_expanding / total_cells * 100
  )

sample_expansion$Timepoint_Group <- factor(sample_expansion$Timepoint_Group,
                                           levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))

p_C <- ggplot(sample_expansion, aes(x = Timepoint_Group, y = pct_expanding, fill = Timepoint_Group)) +
  geom_boxplot(width = 0.5, outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(c("PostART_Suppressed", "PostART_Unsuppressed")),
    method = "wilcox.test", label = "p.signif",
    size = 6, step.increase = 0.12, tip.length = 0.01, hide.ns = FALSE
  ) +
  scale_fill_manual(values = art_colors, name = NULL) +
  scale_x_discrete(labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
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

ggsave(file.path(fig3_dir, "Fig3C_PctExpanding_perSample.png"),
       plot = p_C, width = 6, height = 7, dpi = 300, bg = "white")

# ── Panel D: Volcano plot (all clusters combined) ───────────────────────────
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
  
  dge_all$highlight <- factor(dge_all$highlight,
                              levels = c("Exhaustion", "Stemness", "Cytotoxicity", "Stress response",
                                         "Type I IFN memory", "Chemokines", "Terminal diff.", "Activation", "Other"))
  
  genes_to_label <- c(exhaustion_genes, stemness_genes, cytotoxic_genes,
                      stress_genes, ifn_memory_genes, chemokine_genes,
                      terminal_genes, activation_genes)
  genes_to_label <- genes_to_label[genes_to_label %in% dge_all$gene]
  
  dge_all$label <- ""
  dge_all$label[dge_all$gene %in% genes_to_label & dge_all$p_val_adj < 0.05] <- 
    dge_all$gene[dge_all$gene %in% genes_to_label & dge_all$p_val_adj < 0.05]
  
  # Trim x-axis to data range + small buffer (removes empty space)
  x_lim <- max(abs(dge_all$avg_log2FC[is.finite(dge_all$avg_log2FC)]), na.rm = TRUE) * 1.15
  
  p_D <- ggplot(dge_all, aes(x = avg_log2FC, y = neg_log10_padj, color = highlight)) +
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
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    labs(
      x = expression(log[2]~FC~(suppressed / unsuppressed)),
      y = expression(-log[10]~adjusted~italic(p)),
      title = "Expanding clones: suppressed vs unsuppressed"
    ) +
    theme_cowplot(font_size = 20) +
    theme(
      legend.position  = "right",
      legend.text      = element_text(size = 22),
      legend.title     = element_text(size = 23, face = "bold"),
      legend.key.size  = unit(1.5, "cm"),
      legend.spacing.y = unit(0.3, "cm"),
      axis.text        = element_text(size = 13),
      axis.title       = element_text(size = 14),
      plot.title       = element_text(size = 16),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 9, alpha = 1)))
  
  ggsave(file.path(fig3_dir, "Fig3D_Volcano_AllClusters.png"),
         plot = p_D, width = 14, height = 10, dpi = 300, bg = "white")
}

# ── Fig3 Panel E: Cohen's d bar plot (all modules, Sup vs Unsup) ────────────
cohens_path <- file.path(analysis_dir, "10_module_scores", "ModuleScore_Wilcoxon_all_comparisons.csv")

if (file.exists(cohens_path)) {
  stat_fc <- read.csv(cohens_path)
  
  # Recompute Cohen's d from the module score data
  module_path_e <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")
  if (file.exists(module_path_e)) {
    mod_scores_e <- read.csv(module_path_e, row.names = 1)
    
    module_labels <- c(
      "Exhaustion" = "Exhaustion", "Stemness" = "Stemness / naïve",
      "Cytotoxicity" = "Cytotoxicity", "Stress_IFN" = "Stress / IFN (acute)",
      "TypeI_IFN_Memory" = "Type I IFN memory",
      "Inflammatory_Chemokines" = "Inflammatory chemokines",
      "Terminal_Differentiation" = "Terminal differentiation"
    )
    
    available_mods <- names(module_labels)[names(module_labels) %in% colnames(mod_scores_e)]
    eff_clusters <- c("Effector CD8", "Terminal TEMRA CD8", "Transitional Tem CD8")
    
    cd_results <- list()
    for (mod in available_mods) {
      for (cl in eff_clusters) {
        cl_data <- mod_scores_e %>% filter(CD8_Annotation == cl)
        sup <- cl_data %>% filter(Timepoint_Group == "PostART_Suppressed") %>% pull(!!sym(mod))
        unsup <- cl_data %>% filter(Timepoint_Group == "PostART_Unsuppressed") %>% pull(!!sym(mod))
        if (length(sup) >= 5 & length(unsup) >= 5) {
          pooled_sd <- sqrt(((length(sup)-1)*sd(sup)^2 + (length(unsup)-1)*sd(unsup)^2) / (length(sup)+length(unsup)-2))
          cd <- (mean(sup) - mean(unsup)) / pooled_sd
          wt <- wilcox.test(sup, unsup)
          p_adj <- wt$p.value  # will BH-correct below
          cd_results[[paste(mod, cl)]] <- data.frame(
            Module = mod, CD8_Annotation = cl, cohens_d = cd, p_value = p_adj)
        }
      }
    }
    cd_df <- do.call(rbind, cd_results)
    cd_df$p_adj <- p.adjust(cd_df$p_value, method = "BH")
    cd_df$star <- ifelse(cd_df$p_adj < 0.001, "***",
                         ifelse(cd_df$p_adj < 0.01, "**",
                                ifelse(cd_df$p_adj < 0.05, "*", "")))
    cd_df$direction <- ifelse(cd_df$cohens_d > 0, "Suppressed", "Unsuppressed")
    cd_df$Module_Label <- module_labels[cd_df$Module]
    cd_df$Module_Label <- factor(cd_df$Module_Label, levels = rev(module_labels))
    cd_df$CD8_Annotation <- factor(cd_df$CD8_Annotation, levels = eff_clusters)
    
    x_range <- max(abs(cd_df$cohens_d), na.rm = TRUE)
    cd_df$star_x <- ifelse(cd_df$cohens_d > 0,
                           cd_df$cohens_d + x_range * 0.06,
                           cd_df$cohens_d - x_range * 0.06)
    cd_df$star_hjust <- ifelse(cd_df$cohens_d > 0, 0, 1)
    
    bar_colors_cd <- c("Suppressed" = "#52B788", "Unsuppressed" = "#E76F51")
    
    p_E <- ggplot(cd_df, aes(x = cohens_d, y = Module_Label, fill = direction)) +
      geom_col(width = 0.6, alpha = 0.85) +
      geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
      geom_text(aes(x = star_x, label = star, hjust = star_hjust),
                size = 5, color = "black", vjust = 0.4) +
      facet_wrap(~ CD8_Annotation, nrow = 1) +
      scale_fill_manual(values = bar_colors_cd, name = NULL) +
      coord_cartesian(clip = "off") +
      labs(x = "Cohen's d (← higher in unsuppressed | higher in suppressed →)", y = NULL,
           title = "Module score differences: suppressed vs unsuppressed") +
      theme_cowplot(font_size = 20) +
      theme(
        strip.text        = element_text(size = 13, face = "bold"),
        strip.background  = element_rect(fill = "grey95", color = NA),
        axis.text.y       = element_text(size = 12),
        axis.text.x       = element_text(size = 12),
        axis.title.x      = element_text(size = 13),
        legend.position   = "bottom",
        legend.text       = element_text(size = 12),
        plot.margin       = margin(10, 30, 10, 10),
        plot.background   = element_rect(fill = "white", color = NA),
        panel.background  = element_rect(fill = "white", color = NA)
      )
    
    ggsave(file.path(fig3_dir, "Fig3E_CohenD_AllModules.png"),
           plot = p_E, width = 16, height = 8, dpi = 300, bg = "white")
  }
}

cat("\n=== Figure 3 panels saved to:", fig3_dir, "===\n")

################################################################################
# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Differentiation trajectory & ART-dependent functional divergence
#   A: Monocle3 pseudotime UMAP with trajectory arrows
#   B: Paired patient exhaustion vs viral load
#   C: Exhaustion along pseudotime
#   D: Stemness along pseudotime
#   E: Exhaustion vs viral load (all samples)
#   F: Stemness vs viral load (all samples)
# ══════════════════════════════════════════════════════════════════════════════
################################################################################

# ── Fig4A: Monocle3 trajectory UMAP (pseudotime + arrows) ───────────────────
m3_pt_path <- file.path(analysis_dir, "07_trajectory", "monocle3", "Monocle3_pseudotime_per_cell.csv")

if (file.exists(m3_pt_path)) {
  m3_pt <- read.csv(m3_pt_path)
  
  umap_emb <- Embeddings(TARA_cd8, "wnn.umap")
  pt_umap_df <- data.frame(
    UMAP_1 = umap_emb[, 1], UMAP_2 = umap_emb[, 2],
    cell = colnames(TARA_cd8)
  ) %>% left_join(m3_pt %>% select(cell, monocle3_pseudotime), by = "cell")
  
  # Cluster centroids for trajectory paths
  cluster_centroids <- pt_umap_df %>%
    left_join(data.frame(cell = colnames(TARA_cd8), Cluster = TARA_cd8$CD8_Annotation), by = "cell") %>%
    filter(!is.na(monocle3_pseudotime)) %>%
    group_by(Cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),
              median_pt = median(monocle3_pseudotime, na.rm = TRUE), .groups = "drop")
  
  make_path_df <- function(centroids, path_clusters, lineage_name) {
    path <- centroids %>% filter(Cluster %in% path_clusters)
    path$Cluster <- factor(path$Cluster, levels = path_clusters)
    path <- path %>% arrange(Cluster)
    if (nrow(path) >= 2) { path$Lineage <- lineage_name; return(path) }
    NULL
  }
  
  path1 <- make_path_df(cluster_centroids,
                        c("Naïve CD8", "Naïve CD8 (innate-like)", "Naïve CD8 (post-ART)"), "Naïve trajectory")
  path2 <- make_path_df(cluster_centroids,
                        c("Naïve CD8", "Transitional Tem CD8", "Effector CD8", "Terminal TEMRA CD8"), "Effector trajectory")
  traj_paths <- bind_rows(path1, path2)
  
  p_4A <- ggplot(pt_umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = monocle3_pseudotime), size = 0.5, alpha = 0.6) +
    geom_path(data = traj_paths, aes(group = Lineage),
              color = "black", linewidth = 1.8, arrow = arrow(length = unit(0.35, "cm"), type = "closed")) +
    geom_point(data = traj_paths, color = "black", size = 4) +
    scale_color_viridis_c(option = "inferno", name = "Pseudotime", na.value = "grey85") +
    labs(title = "Monocle3 pseudotime with trajectories") +
    NoAxes() +
    theme_cowplot(font_size = 22) +
    theme(
      legend.text       = element_text(size = 16),
      legend.title      = element_text(size = 18, face = "bold"),
      legend.key.height = unit(1.5, "cm"),
      plot.title        = element_text(size = 22, face = "bold"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(fig4_dir, "Fig4A_Monocle3_pseudotime_UMAP.png"),
         plot = p_4A, width = 11, height = 9, dpi = 300, bg = "white")
}

# ── Fig4B: Paired patient exhaustion vs viral load ───────────────────────────
vl_path <- file.path(analysis_dir, "09_viral_load", "PerSample_ModuleScores_and_ViralLoad.csv")
paired_vl_path <- file.path(analysis_dir, "09_viral_load", "H4_Paired_Exhaustion_vs_ViralLoad.png")

if (file.exists(vl_path)) {
  sample_scores <- read.csv(vl_path)
  sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)
  
  # Paired patients (>1 timepoint)
  paired_patients <- sample_scores %>%
    group_by(PID) %>% filter(n_distinct(Timepoint_Group) > 1) %>% ungroup()
  
  if (nrow(paired_patients) > 2 & "Exhaustion" %in% colnames(paired_patients)) {
    p_4B <- ggplot(paired_patients, aes(x = log10_VL, y = Exhaustion)) +
      geom_line(aes(group = PID), color = "grey40", linewidth = 0.8, linetype = "dashed") +
      geom_point(aes(color = Timepoint_Group), size = 6, alpha = 0.9) +
      geom_text_repel(aes(label = paste0(PID, "\n",
                                         c("PreART_Entry"="Pre","PostART_Suppressed"="Sup","PostART_Unsuppressed"="Unsup")[as.character(Timepoint_Group)])),
                      size = 5, max.overlaps = 20, lineheight = 0.8) +
      scale_color_manual(values = art_colors, name = "ART Status",
                         labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
      labs(x = expression(log[10]~viral~load), y = "Mean exhaustion score",
           title = "Paired patients: exhaustion tracks viral load") +
      theme_cowplot(font_size = 22) +
      theme(
        axis.text        = element_text(size = 16),
        axis.title       = element_text(size = 18),
        plot.title       = element_text(size = 20, face = "bold"),
        legend.text      = element_text(size = 16),
        legend.title     = element_text(size = 17, face = "bold"),
        legend.position  = "right",
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      guides(color = guide_legend(override.aes = list(size = 5)))
    
    ggsave(file.path(fig4_dir, "Fig4B_Paired_Exhaustion_vs_ViralLoad.png"),
           plot = p_4B, width = 10, height = 8, dpi = 300, bg = "white")
  }
}

# ── Fig4C: Exhaustion along pseudotime ───────────────────────────────────────
module_path <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")

if (file.exists(m3_pt_path) & file.exists(module_path)) {
  mod_scores <- read.csv(module_path, row.names = 1)
  
  common <- intersect(m3_pt$cell, rownames(mod_scores))
  pt_mod <- data.frame(
    pseudotime = m3_pt$monocle3_pseudotime[match(common, m3_pt$cell)],
    Timepoint  = mod_scores[common, "Timepoint_Group"],
    Exhaustion = mod_scores[common, "Exhaustion"],
    Stemness   = mod_scores[common, "Stemness"]
  ) %>% filter(!is.na(pseudotime) & is.finite(pseudotime))
  
  p_4C <- ggplot(pt_mod, aes(x = pseudotime, y = Exhaustion, color = Timepoint)) +
    geom_point(size = 0.15, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    labs(x = "Monocle3 pseudotime", y = "Exhaustion score",
         title = "Exhaustion along differentiation") +
    theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "bottom",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  ggsave(file.path(fig4_dir, "Fig4C_Exhaustion_pseudotime.png"),
         plot = p_4C, width = 9, height = 7, dpi = 300, bg = "white")
  
  # ── Fig4D: Stemness along pseudotime ─────────────────────────────────────
  p_4D <- ggplot(pt_mod, aes(x = pseudotime, y = Stemness, color = Timepoint)) +
    geom_point(size = 0.15, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    labs(x = "Monocle3 pseudotime", y = "Stemness score",
         title = "Stemness along differentiation") +
    theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "bottom",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  ggsave(file.path(fig4_dir, "Fig4D_Stemness_pseudotime.png"),
         plot = p_4D, width = 9, height = 7, dpi = 300, bg = "white")
}

# ── Fig4E: Exhaustion vs viral load (all samples) ───────────────────────────
if (file.exists(vl_path)) {
  if (!exists("sample_scores")) {
    sample_scores <- read.csv(vl_path)
    sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)
  }
  
  cor_exh <- cor.test(sample_scores$log10_VL, sample_scores$Exhaustion,
                      method = "spearman", exact = FALSE)
  
  p_4E <- ggplot(sample_scores, aes(x = log10_VL, y = Exhaustion)) +
    geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
    geom_text_repel(aes(label = PID), size = 5, max.overlaps = 15, color = "grey30") +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    scale_size_continuous(name = "N expanding", range = c(3, 9)) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("rho = %.2f\np = %.2g", cor_exh$estimate, cor_exh$p.value),
             hjust = -0.1, vjust = 1.3, size = 8, fontface = "bold") +
    labs(x = expression(log[10]~viral~load), y = "Mean exhaustion score",
         title = "Exhaustion vs viral load") +
    theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  ggsave(file.path(fig4_dir, "Fig4E_Exhaustion_vs_ViralLoad.png"),
         plot = p_4E, width = 11, height = 8, dpi = 300, bg = "white")
  
  # ── Fig4F: Stemness vs viral load ────────────────────────────────────────
  cor_stem <- cor.test(sample_scores$log10_VL, sample_scores$Stemness,
                       method = "spearman", exact = FALSE)
  
  p_4F <- ggplot(sample_scores, aes(x = log10_VL, y = Stemness)) +
    geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
    geom_text_repel(aes(label = PID), size = 5, max.overlaps = 15, color = "grey30") +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    scale_size_continuous(name = "N expanding", range = c(3, 9)) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("rho = %.2f\np = %.2g", cor_stem$estimate, cor_stem$p.value),
             hjust = -0.1, vjust = 1.3, size = 8, fontface = "bold") +
    labs(x = expression(log[10]~viral~load), y = "Mean stemness score",
         title = "Stemness vs viral load") +
    theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  ggsave(file.path(fig4_dir, "Fig4F_Stemness_vs_ViralLoad.png"),
         plot = p_4F, width = 11, height = 8, dpi = 300, bg = "white")
}

cat("\n=== Figure 4 panels saved to:", fig4_dir, "===\n")

################################################################################
# ══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY FIGURE
# ══════════════════════════════════════════════════════════════════════════════
################################################################################

# ── S1: ADT + RNA annotation heatmap (generated fresh) ───────────────────────
cat("  Generating annotation heatmap...\n")

adt_marker_groups_cd8 <- c(
  "CD45RA"="Naïve", "SELL"="Naïve", "CD7"="Naïve",
  "FAS"="Memory/Tscm", "IL2RB"="Memory/Tscm", "CD45RO"="Memory/Tscm", "IL7R"="Memory/Tscm", "CD44"="Memory/Tscm",
  "CD27"="Co-stimulation", "CD28"="Co-stimulation", "ICOS"="Co-stimulation",
  "B3GAT1"="Effector/TEMRA", "KLRG1"="Effector/TEMRA", "KIR3DL1"="Effector/TEMRA", "CX3CR1"="Effector/TEMRA",
  "TIGIT"="Exhaustion", "PDCD1"="Exhaustion", "LAG3"="Exhaustion",
  "CD38"="Activation", "CD69"="Activation", "ENTPD1"="Activation", "NT5E"="Activation",
  "CXCR3"="Homing", "ITGB7"="Homing", "ITGA1"="Homing",
  "NCAM1"="NK-like/Innate", "FCGR3A"="NK-like/Innate", "SIGLEC7"="NK-like/Innate", "KLRD1"="NK-like/Innate", "KLRB1"="NK-like/Innate",
  "TCR-AB"="TCR Identity", "TCR-vA7.2"="TCR Identity", "TCR-vD2"="TCR Identity"
)

rna_marker_groups_cd8 <- c(
  "CCR7"="Naïve/Stemness", "SELL"="Naïve/Stemness", "TCF7"="Naïve/Stemness", "LEF1"="Naïve/Stemness",
  "IL7R"="Naïve/Stemness", "KLF2"="Naïve/Stemness", "S1PR1"="Naïve/Stemness", "FOXP1"="Naïve/Stemness",
  "BACH2"="Memory/Survival", "BCL2"="Memory/Survival", "CD27"="Memory/Survival", "CD28"="Memory/Survival", "ID3"="Memory/Survival",
  "FAS"="Transitional", "IL2RB"="Transitional", "CD44"="Transitional", "CXCR3"="Transitional",
  "GZMK"="Effector Memory", "EOMES"="Effector Memory",
  "GZMB"="TEMRA/Cytotoxic", "GNLY"="TEMRA/Cytotoxic", "PRF1"="TEMRA/Cytotoxic", "NKG7"="TEMRA/Cytotoxic",
  "TBX21"="TEMRA/Cytotoxic", "CX3CR1"="TEMRA/Cytotoxic", "FGFBP2"="TEMRA/Cytotoxic",
  "GZMA"="TEMRA/Cytotoxic", "GZMH"="TEMRA/Cytotoxic", "GZMM"="TEMRA/Cytotoxic",
  "RUNX3"="Effector TFs", "ZEB2"="Effector TFs", "PRDM1"="Effector TFs", "ID2"="Effector TFs",
  "TOX"="Exhaustion", "PDCD1"="Exhaustion", "TIGIT"="Exhaustion", "HAVCR2"="Exhaustion", "LAG3"="Exhaustion",
  "MKI67"="Proliferation", "TOP2A"="Proliferation",
  "ITGAE"="Tissue Residency", "CXCR6"="Tissue Residency", "CD69"="Tissue Residency",
  "IFNG"="Cytokines", "TNF"="Cytokines",
  "TRDV1"="γδ TCR", "TRDV2"="γδ TCR", "TRGV9"="γδ TCR", "TRDC"="γδ TCR",
  "TYROBP"="NK-like/Innate", "KLRD1"="NK-like/Innate", "KLRB1"="NK-like/Innate", "FCGR3A"="NK-like/Innate",
  "SLC4A10"="MAIT", "ZBTB16"="MAIT", "NCR3"="MAIT",
  "HLA-DRA"="Activation", "CD38"="Activation", "S1PR5"="Activation"
)

group_colors_cd8 <- c(
  "Naïve"="#74C2E1", "Naïve/Stemness"="#74C2E1", "Memory/Tscm"="#66BB6A", "Memory/Survival"="#66BB6A",
  "Co-stimulation"="#AED581", "Transitional"="#FFCC80", "Effector Memory"="#FFA726",
  "Effector/TEMRA"="#EF5350", "TEMRA/Cytotoxic"="#EF5350", "Effector TFs"="#FF8A65",
  "Exhaustion"="#A1887F", "Proliferation"="#CE93D8", "Tissue Residency"="#80DEEA",
  "Cytokines"="#FFE082", "Activation"="#FFD54F", "Homing"="#B2DFDB",
  "γδ TCR"="#9C6FD6", "NK-like/Innate"="#F06292", "MAIT"="#DCE775", "TCR Identity"="#B0BEC5"
)

col_order_cd8 <- c("Naïve CD8", "Naïve CD8 (innate-like)", "Naïve CD8 (post-ART)",
                   "Transitional Tem CD8", "Effector CD8", "Terminal TEMRA CD8",
                   "KIR+ innate-like CD8", "MAIT-like Trm", "TRDV1+ γδ T cell",
                   "Naïve-like TRDV1+ γδ", "Vγ9Vδ2 γδ T cell")

# Compute averages
DefaultAssay(TARA_cd8) <- "ADT"
avg_adt <- AverageExpression(TARA_cd8, assays = "ADT", features = names(adt_marker_groups_cd8),
                             group.by = "CD8_Annotation", slot = "data")$ADT
DefaultAssay(TARA_cd8) <- "RNA"
avg_rna <- AverageExpression(TARA_cd8, assays = "RNA", features = names(rna_marker_groups_cd8),
                             group.by = "CD8_Annotation", slot = "data")$RNA

scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) { r <- max(x) - min(x); if (r == 0) rep(0.5, length(x)) else (x - min(x)) / r }))
}

avg_adt_s <- scale_01(log10(avg_adt + 1))
avg_rna_s <- scale_01(avg_rna)
colnames(avg_adt_s) <- gsub("^g ", "", colnames(avg_adt_s))
colnames(avg_rna_s) <- gsub("^g ", "", colnames(avg_rna_s))

col_ord <- col_order_cd8[col_order_cd8 %in% colnames(avg_adt_s)]
avg_adt_s <- avg_adt_s[, col_ord]; avg_rna_s <- avg_rna_s[, col_ord]

# Row annotations + within-group clustering + gaps
make_annot_gaps <- function(mat, marker_groups, grp_colors) {
  annot <- data.frame(Group = factor(marker_groups[rownames(mat)], levels = names(grp_colors)), row.names = rownames(mat))
  
  # Cluster within groups (Ward D2)
  groups <- levels(droplevels(annot$Group))
  row_order <- character(0)
  for (grp in groups) {
    rows <- rownames(annot)[annot$Group == grp]
    if (length(rows) == 1) { row_order <- c(row_order, rows)
    } else {
      hc <- hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2")
      row_order <- c(row_order, rows[hc$order])
    }
  }
  
  annot_ordered <- annot[row_order, , drop = FALSE]
  grps <- as.character(annot_ordered$Group)
  gaps <- which(grps[-length(grps)] != grps[-1])
  
  list(mat = mat[row_order, ], annot = annot_ordered, gaps = gaps,
       colors = list(Group = grp_colors[levels(droplevels(annot_ordered$Group))]))
}

ag_adt <- make_annot_gaps(avg_adt_s, adt_marker_groups_cd8, group_colors_cd8)
ag_rna <- make_annot_gaps(avg_rna_s, rna_marker_groups_cd8, group_colors_cd8)

hm_cols <- colorRampPalette(c("#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"))(100)

# MUCH larger fonts for readability
p_adt_hm <- pheatmap(ag_adt$mat, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
                     color = hm_cols, border_color = "white", annotation_row = ag_adt$annot,
                     annotation_colors = ag_adt$colors, annotation_names_row = FALSE,
                     gaps_row = ag_adt$gaps, cellwidth = 60, cellheight = 22,
                     fontsize = 16, fontsize_row = 14, fontsize_col = 14, angle_col = 45,
                     main = "Protein (ADT)", silent = TRUE)

p_rna_hm <- pheatmap(ag_rna$mat, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
                     color = hm_cols, border_color = "white", annotation_row = ag_rna$annot,
                     annotation_colors = ag_rna$colors, annotation_names_row = FALSE,
                     gaps_row = ag_rna$gaps, cellwidth = 60, cellheight = 22,
                     fontsize = 16, fontsize_row = 14, fontsize_col = 14, angle_col = 45,
                     main = "mRNA (RNA)", silent = TRUE)

# Save combined
png(file.path(fig_supp, "S1_ADT_RNA_Annotation_Heatmap.png"),
    width = 34, height = 22, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_adt_hm$gtable); popViewport()
pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_rna_hm$gtable); popViewport()
dev.off()

# Also save individually
png(file.path(fig_supp, "S1a_ADT_Heatmap.png"), width = 17, height = 16, units = "in", res = 300, bg = "white")
grid.draw(p_adt_hm$gtable); dev.off()

png(file.path(fig_supp, "S1b_RNA_Heatmap.png"), width = 17, height = 24, units = "in", res = 300, bg = "white")
grid.draw(p_rna_hm$gtable); dev.off()

cat("  Annotation heatmaps saved.\n")

# ── S2: Monocle3 individual module score along pseudotime plots ──────────────
m3_dir <- file.path(analysis_dir, "07_trajectory", "monocle3")

# Copy individual module plots (one per module, NOT the combined AllModules faceted plot)
m3_individual_files <- list.files(m3_dir, pattern = "^Monocle3_[A-Z].*_along_pseudotime\\.png$", full.names = TRUE)
# Exclude the combined AllModules and ExhaustionStemness (those are in main figure)
m3_individual_files <- m3_individual_files[!grepl("AllModules|ExhaustionStemness", m3_individual_files)]

for (f in m3_individual_files) {
  file.copy(f, file.path(fig_supp, paste0("S2_", basename(f))), overwrite = TRUE)
}

# Also copy gene expression along pseudotime plots and trajectory UMAPs
m3_gene_files <- list.files(m3_dir, pattern = "Genes_in_pseudotime|NarrativeGenes|TopVariable|TopGenes_UMAP|UMAP_clusters|UMAP_ART", full.names = TRUE)
for (f in m3_gene_files) {
  file.copy(f, file.path(fig_supp, paste0("S2_", basename(f))), overwrite = TRUE)
}

cat("  Individual module + gene pseudotime plots copied to supplementary.\n")

# ── S3: Naïve trajectory volcano (Naïve CD8 vs innate-like) ────────────────
cat("  Generating naïve trajectory volcano...\n")
naive_dge_path <- file.path(analysis_dir, "03_DGE", "DGE_MAST_RNA_NaveCD8_vs_NaveCD8innatelike.csv")
if (!file.exists(naive_dge_path)) {
  naive_dge_path <- "/mnt/user-data/uploads/DGE_MAST_RNA_NaveCD8_vs_NaveCD8innatelike.csv"
}

if (file.exists(naive_dge_path)) {
  dge_naive <- read.csv(naive_dge_path, row.names = 1)
  dge_naive$gene <- rownames(dge_naive)
  dge_naive$neg_log10_padj <- -log10(dge_naive$p_val_adj + 1e-300)
  
  quiescence_genes  <- c("KLF2", "S1PR1", "SELL", "BACH2", "BCL2", "IL7R", "FOXP1", "CCR7", "LEF1", "TCF7")
  hypoxia_genes     <- c("IGFBP2", "HILPDA", "PFKFB4", "EGLN3", "GBP5", "AK4", "HIF1A")
  innate_genes      <- c("TYROBP", "KLRD1", "KLRB1", "SIGLEC7", "KIR3DL1", "FCGR3A")
  ifn_genes_naive   <- c("MX1", "IFIT1", "ISG15", "IFI44L")
  
  naive_highlight_cols <- c(
    "Quiescence / homing" = "#185FA5", "Hypoxia / metabolism" = "#A32D2D",
    "Innate-like" = "#BA7517", "Type I IFN" = "#0F6E56", "Other" = "grey80"
  )
  
  dge_naive$highlight <- "Other"
  dge_naive$highlight[dge_naive$gene %in% quiescence_genes]  <- "Quiescence / homing"
  dge_naive$highlight[dge_naive$gene %in% hypoxia_genes]     <- "Hypoxia / metabolism"
  dge_naive$highlight[dge_naive$gene %in% innate_genes]      <- "Innate-like"
  dge_naive$highlight[dge_naive$gene %in% ifn_genes_naive]   <- "Type I IFN"
  dge_naive$highlight <- factor(dge_naive$highlight,
                                levels = c("Quiescence / homing", "Hypoxia / metabolism", "Innate-like", "Type I IFN", "Other"))
  
  genes_label_naive <- c(quiescence_genes, hypoxia_genes, innate_genes, ifn_genes_naive)
  genes_label_naive <- genes_label_naive[genes_label_naive %in% dge_naive$gene]
  dge_naive$label <- ""
  dge_naive$label[dge_naive$gene %in% genes_label_naive & dge_naive$p_val_adj < 0.05] <-
    dge_naive$gene[dge_naive$gene %in% genes_label_naive & dge_naive$p_val_adj < 0.05]
  
  x_lim_n <- max(abs(dge_naive$avg_log2FC[is.finite(dge_naive$avg_log2FC)]), na.rm = TRUE) * 1.15
  
  p_naive_volc <- ggplot(dge_naive, aes(x = avg_log2FC, y = neg_log10_padj, color = highlight)) +
    geom_point(data = dge_naive %>% filter(highlight == "Other"), size = 2, alpha = 0.3) +
    geom_point(data = dge_naive %>% filter(highlight != "Other"), size = 7, alpha = 0.85) +
    geom_label_repel(
      data = dge_naive %>% filter(label != ""),
      aes(label = label, fill = highlight), color = "white",
      size = 8, fontface = "bold.italic",
      max.overlaps = 50, segment.size = 0.4, segment.color = "grey40",
      min.segment.length = 0.1, box.padding = 0.5,
      label.size = 0.2, show.legend = FALSE
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = naive_highlight_cols, name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = naive_highlight_cols, guide = "none") +
    coord_cartesian(xlim = c(-x_lim_n, x_lim_n)) +
    labs(
      x = expression(log[2]~FC~(classical~naïve / innate-like~naïve)),
      y = expression(-log[10]~adjusted~italic(p)),
      title = "Naïve trajectory: classical naïve vs innate-like"
    ) +
    theme_cowplot(font_size = 20) +
    theme(
      legend.position  = "right",
      legend.text      = element_text(size = 22),
      legend.title     = element_text(size = 23, face = "bold"),
      legend.key.size  = unit(1.5, "cm"),
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 9, alpha = 1)))
  
  ggsave(file.path(fig_supp, "S3_Volcano_Naive_vs_InnateLike.png"),
         plot = p_naive_volc, width = 14, height = 10, dpi = 300, bg = "white")
  cat("  Naïve volcano saved.\n")
}

cat("\n=== Supplementary panels saved to:", fig_supp, "===\n")
cat("Figure 3 panels:", length(list.files(fig3_dir, pattern = "\\.png$")), "\n")
cat("Figure 4 panels:", length(list.files(fig4_dir, pattern = "\\.png$")), "\n")
cat("Supplementary panels:", length(list.files(fig_supp, pattern = "\\.png$|\\.csv$")), "\n")

sessionInfo()