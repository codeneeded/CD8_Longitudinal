################################################################################
# FIGURE 5 PLOTS: Differentiation trajectory & viral load correlations
#
# Panels (per manuscript):
#   5A — Monocle3 pseudotime UMAP with trajectory arrows
#   5B — Pseudotime density by ART status
#   5C — Exhaustion along pseudotime (loess by ART status)
#   5D — Stemness along pseudotime (loess by ART status)
#   5E — Exhaustion vs viral load (per-sample scatter)
#   5F — Stemness vs viral load (per-sample scatter)
#
# Prerequisite: Run Fig4_5_Prep_updated.R first to generate:
#   - TARA_cd8_HEI_annotated_final.qs2
#   - Monocle3_pseudotime_per_cell.csv
#   - ModuleScores_per_cell.csv
#   - PerSample_ModuleScores_and_ViralLoad.csv
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(qs2)
library(SeuratExtend)
library(ggrepel)
library(ggpubr)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
manuscript   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 4-5"
fig5_dir     <- file.path(manuscript, "Figure5_panels")
analysis_dir <- file.path(manuscript, "analysis")

dir.create(fig5_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load annotated object ────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells\n")

# ── Shared aesthetics ────────────────────────────────────────────────────────
art_colors <- c(
  "PreART_Entry"           = "#4A90D9",
  "PostART_Suppressed"     = "#52B788",
  "PostART_Unsuppressed"   = "#E76F51"
)

# UPDATED cluster names
col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
  "Transitional Tem CD8", "TEM CD8", "TEMRA CD8",
  "KIR+ innate-like CD8", "MAIT-like Trm",
  "γδ1 T cell", "Naïve γδ1 T cell", "γδ2 T cell"
)

# ── Data paths ───────────────────────────────────────────────────────────────
m3_pt_path   <- file.path(analysis_dir, "07_trajectory", "monocle3", "Monocle3_pseudotime_per_cell.csv")
module_path  <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")
vl_path      <- file.path(analysis_dir, "09_viral_load", "PerSample_ModuleScores_and_ViralLoad.csv")

################################################################################
# Panel A: Monocle3 trajectory UMAP (pseudotime + arrows)
################################################################################
cat("  Panel A: Monocle3 trajectory UMAP...\n")

if (file.exists(m3_pt_path)) {
  m3_pt <- read.csv(m3_pt_path)
  
  umap_emb <- Embeddings(TARA_cd8, "wnn.umap")
  pt_umap_df <- data.frame(
    UMAP_1 = umap_emb[, 1], UMAP_2 = umap_emb[, 2],
    cell = colnames(TARA_cd8)
  ) %>% left_join(m3_pt %>% select(cell, monocle3_pseudotime), by = "cell")
  
  # Cluster centroids for trajectory path arrows
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
  
  # UPDATED trajectory paths with new cluster names
  path1 <- make_path_df(cluster_centroids,
                        c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3"), "Naïve trajectory")
  path2 <- make_path_df(cluster_centroids,
                        c("Naïve CD8", "Transitional Tem CD8", "TEM CD8", "TEMRA CD8"), "Effector trajectory")
  traj_paths <- bind_rows(path1, path2)
  
  # Compute arrow position — use each axis range so arrows look equal on screen
  x_range <- range(pt_umap_df$UMAP_1, na.rm = TRUE)
  y_range <- range(pt_umap_df$UMAP_2, na.rm = TRUE)
  arrow_x0 <- x_range[1] + diff(x_range) * 0.02
  arrow_y0 <- y_range[1] + diff(y_range) * 0.02
  arrow_len_x <- diff(x_range) * 0.15   # horizontal arrow uses x-axis scale
  arrow_len_y <- diff(y_range) * 0.15   # vertical arrow uses y-axis scale
  
  p_5A <- ggplot(pt_umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = monocle3_pseudotime), size = 0.5, alpha = 0.6) +
    geom_path(data = traj_paths, aes(group = Lineage),
              color = "black", linewidth = 1.8, arrow = arrow(length = unit(0.35, "cm"), type = "closed")) +
    geom_point(data = traj_paths, color = "black", size = 4) +
    scale_color_viridis_c(option = "inferno", name = "Pseudotime", na.value = "grey85") +
    # wnnUMAP arrows in bottom-left corner
    annotate("segment", x = arrow_x0, xend = arrow_x0 + arrow_len_x, y = arrow_y0, yend = arrow_y0,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.8) +
    annotate("segment", x = arrow_x0, xend = arrow_x0, y = arrow_y0, yend = arrow_y0 + arrow_len_y,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.8) +
    annotate("text", x = arrow_x0 + arrow_len_x / 2, y = arrow_y0 - diff(y_range) * 0.03,
             label = "wnnUMAP 1", size = 5, fontface = "plain") +
    annotate("text", x = arrow_x0 - diff(x_range) * 0.03, y = arrow_y0 + arrow_len_y / 2,
             label = "wnnUMAP 2", size = 5, fontface = "plain", angle = 90) +
    labs(title = "Monocle3 pseudotime with trajectories", x = NULL, y = NULL) +
    theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      legend.text       = element_text(size = 16),
      legend.title      = element_text(size = 18, face = "bold"),
      legend.key.height = unit(1.5, "cm"),
      plot.title        = element_text(size = 22, face = "bold"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(fig5_dir, "Fig5A_Monocle3_pseudotime_UMAP.png"),
         plot = p_5A, width = 11, height = 9, dpi = 300, bg = "white")
} else {
  cat("    WARNING: Monocle3 pseudotime file not found\n")
}

################################################################################
# Panel B: Pseudotime density by ART status
################################################################################
cat("  Panel B: Pseudotime density...\n")

if (exists("m3_pt")) {
  pt_density_df <- m3_pt %>%
    filter(!is.na(monocle3_pseudotime) & is.finite(monocle3_pseudotime))
  
  # Handle column name variation
  if ("Timepoint" %in% colnames(pt_density_df) & !"Timepoint_Group" %in% colnames(pt_density_df)) {
    pt_density_df$Timepoint_Group <- pt_density_df$Timepoint
  }
  
  pt_density_df$Timepoint_Group <- factor(pt_density_df$Timepoint_Group,
                                          levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
  pt_density_df <- pt_density_df %>% filter(!is.na(Timepoint_Group))
  
  # Median pseudotime per group for vertical lines
  medians_pt <- pt_density_df %>%
    group_by(Timepoint_Group) %>%
    summarise(median_pt = median(monocle3_pseudotime, na.rm = TRUE), .groups = "drop")
  
  p_5B <- ggplot(pt_density_df, aes(x = monocle3_pseudotime, fill = Timepoint_Group,
                                    color = Timepoint_Group)) +
    geom_density(alpha = 0.35, linewidth = 1.2) +
    geom_vline(data = medians_pt, aes(xintercept = median_pt, color = Timepoint_Group),
               linewidth = 1, linetype = "dashed", show.legend = FALSE) +
    scale_fill_manual(values = art_colors, name = "ART Status",
                      labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    labs(x = "Monocle3 pseudotime", y = "Density",
         title = "Distribution along differentiation trajectory",
         subtitle = "Dashed lines = median pseudotime per group") +
    theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      plot.subtitle    = element_text(size = 14, color = "grey40"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "bottom",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.6)), color = "none")
  
  ggsave(file.path(fig5_dir, "Fig5B_Pseudotime_density_by_ART.png"),
         plot = p_5B, width = 10, height = 7, dpi = 300, bg = "white")
}

################################################################################
# Panels C, D & E: Exhaustion, Stemness, and Type I IFN memory along pseudotime
################################################################################
cat("  Panels C, D & E: Module scores along pseudotime...\n")

if (exists("m3_pt") & file.exists(module_path)) {
  mod_scores <- read.csv(module_path, row.names = 1)
  
  common <- intersect(m3_pt$cell, rownames(mod_scores))
  pt_mod <- data.frame(
    pseudotime = m3_pt$monocle3_pseudotime[match(common, m3_pt$cell)],
    Timepoint  = mod_scores[common, "Timepoint_Group"],
    Exhaustion = mod_scores[common, "Exhaustion"],
    Stemness   = mod_scores[common, "Stemness"],
    TypeI_IFN_Memory = mod_scores[common, "TypeI_IFN_Memory"]
  ) %>% filter(!is.na(pseudotime) & is.finite(pseudotime))
  
  pt_mod$Timepoint <- factor(pt_mod$Timepoint,
                             levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
  
  # Shared theme for pseudotime panels
  pt_theme <- theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "bottom",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # ── Panel C: Exhaustion ────────────────────────────────────────────────────
  p_5C <- ggplot(pt_mod, aes(x = pseudotime, y = Exhaustion, color = Timepoint)) +
    geom_point(size = 0.15, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    labs(x = "Monocle3 pseudotime", y = "Exhaustion score",
         title = "Exhaustion along differentiation") +
    pt_theme +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  ggsave(file.path(fig5_dir, "Fig5C_Exhaustion_pseudotime.png"),
         plot = p_5C, width = 9, height = 7, dpi = 300, bg = "white")
  
  # ── Panel D: Stemness ──────────────────────────────────────────────────────
  p_5D <- ggplot(pt_mod, aes(x = pseudotime, y = Stemness, color = Timepoint)) +
    geom_point(size = 0.15, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
    labs(x = "Monocle3 pseudotime", y = "Stemness score",
         title = "Stemness along differentiation") +
    pt_theme +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
  
  ggsave(file.path(fig5_dir, "Fig5D_Stemness_pseudotime.png"),
         plot = p_5D, width = 9, height = 7, dpi = 300, bg = "white")
  
  # ── Panel E: Type I IFN memory ─────────────────────────────────────────────
  if ("TypeI_IFN_Memory" %in% colnames(pt_mod) && any(!is.na(pt_mod$TypeI_IFN_Memory))) {
    p_5E <- ggplot(pt_mod, aes(x = pseudotime, y = TypeI_IFN_Memory, color = Timepoint)) +
      geom_point(size = 0.15, alpha = 0.1) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
      scale_color_manual(values = art_colors, name = "ART Status",
                         labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
      labs(x = "Monocle3 pseudotime", y = "Type I IFN memory score",
           title = "Type I IFN memory along differentiation") +
      pt_theme +
      guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
    
    ggsave(file.path(fig5_dir, "Fig5E_TypeI_IFN_Memory_pseudotime.png"),
           plot = p_5E, width = 9, height = 7, dpi = 300, bg = "white")
  } else {
    cat("    WARNING: TypeI_IFN_Memory not found in module scores — skipping 5E\n")
  }
} else {
  cat("    WARNING: Cannot generate C/D/E — missing pseudotime or module score files\n")
}

################################################################################
# Panels F, G & H: Exhaustion, Stemness, and IFN memory vs viral load
################################################################################
cat("  Panels F, G & H: Module scores vs viral load...\n")

if (file.exists(vl_path)) {
  sample_scores <- read.csv(vl_path)
  sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)
  
  sample_scores$Timepoint_Group <- factor(sample_scores$Timepoint_Group,
                                          levels = c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
  
  # Shared theme for VL panels
  vl_theme <- theme_cowplot(font_size = 22) +
    theme(
      axis.text        = element_text(size = 16),
      axis.title       = element_text(size = 18),
      plot.title       = element_text(size = 20, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 17, face = "bold"),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # ── Panel F: Exhaustion vs viral load ──────────────────────────────────────
  if ("Exhaustion" %in% colnames(sample_scores)) {
    cor_exh <- cor.test(sample_scores$log10_VL, sample_scores$Exhaustion,
                        method = "spearman", exact = FALSE)
    
    p_5F <- ggplot(sample_scores, aes(x = log10_VL, y = Exhaustion)) +
      geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
      geom_text_repel(aes(label = PID), size = 5, max.overlaps = 15, color = "grey30") +
      scale_color_manual(values = art_colors, name = "ART Status",
                         labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
      scale_size_continuous(name = "N expanding", range = c(3, 9)) +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("\u03c1 = %.2f\np = %.2g", cor_exh$estimate, cor_exh$p.value),
               hjust = -0.1, vjust = 1.3, size = 8, fontface = "bold") +
      labs(x = expression(log[10]~viral~load), y = "Mean exhaustion score",
           title = "Exhaustion vs viral load") +
      vl_theme +
      guides(color = guide_legend(override.aes = list(size = 5)))
    
    ggsave(file.path(fig5_dir, "Fig5F_Exhaustion_vs_ViralLoad.png"),
           plot = p_5F, width = 11, height = 8, dpi = 300, bg = "white")
  }
  
  # ── Panel G: Stemness vs viral load ────────────────────────────────────────
  if ("Stemness" %in% colnames(sample_scores)) {
    cor_stem <- cor.test(sample_scores$log10_VL, sample_scores$Stemness,
                         method = "spearman", exact = FALSE)
    
    p_5G <- ggplot(sample_scores, aes(x = log10_VL, y = Stemness)) +
      geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
      geom_text_repel(aes(label = PID), size = 5, max.overlaps = 15, color = "grey30") +
      scale_color_manual(values = art_colors, name = "ART Status",
                         labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
      scale_size_continuous(name = "N expanding", range = c(3, 9)) +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("\u03c1 = %.2f\np = %.2g", cor_stem$estimate, cor_stem$p.value),
               hjust = -0.1, vjust = 1.3, size = 8, fontface = "bold") +
      labs(x = expression(log[10]~viral~load), y = "Mean stemness score",
           title = "Stemness vs viral load") +
      vl_theme +
      guides(color = guide_legend(override.aes = list(size = 5)))
    
    ggsave(file.path(fig5_dir, "Fig5G_Stemness_vs_ViralLoad.png"),
           plot = p_5G, width = 11, height = 8, dpi = 300, bg = "white")
  }
  
  # ── Panel H: Type I IFN memory vs viral load ──────────────────────────────
  if ("TypeI_IFN_Memory" %in% colnames(sample_scores)) {
    cor_ifn <- cor.test(sample_scores$log10_VL, sample_scores$TypeI_IFN_Memory,
                        method = "spearman", exact = FALSE)
    
    p_5H <- ggplot(sample_scores, aes(x = log10_VL, y = TypeI_IFN_Memory)) +
      geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
      geom_text_repel(aes(label = PID), size = 5, max.overlaps = 15, color = "grey30") +
      scale_color_manual(values = art_colors, name = "ART Status",
                         labels = c("Pre-ART", "Suppressed", "Unsuppressed")) +
      scale_size_continuous(name = "N expanding", range = c(3, 9)) +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("\u03c1 = %.2f\np = %.2g", cor_ifn$estimate, cor_ifn$p.value),
               hjust = -0.1, vjust = 1.3, size = 8, fontface = "bold") +
      labs(x = expression(log[10]~viral~load), y = "Mean Type I IFN memory score",
           title = "Type I IFN memory vs viral load") +
      vl_theme +
      guides(color = guide_legend(override.aes = list(size = 5)))
    
    ggsave(file.path(fig5_dir, "Fig5H_IFNmemory_vs_ViralLoad.png"),
           plot = p_5H, width = 11, height = 8, dpi = 300, bg = "white")
  }
} else {
  cat("    WARNING: Viral load file not found at", vl_path, "\n")
}

cat("\n=== Figure 5 complete:", length(list.files(fig5_dir, "\\.png$")), "panels saved to:", fig5_dir, "===\n")

