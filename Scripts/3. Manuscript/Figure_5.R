################################################################################
# FIGURE 5 PLOTS: Differentiation trajectory & ART-dependent functional divergence
#   v5b — Narrative reframe:
#     5D: Type I IFN memory along PT (was stemness)
#     5F: Naive branch composition by ART (was stemness vs VL)
#     Stemness panels moved to S8
#
# 6 panels, 3 rows:
#   Row 1: 5A (trajectory UMAP) + 5B (pseudotime density)
#   Row 2: 5C (exhaustion along PT) + 5D (IFN memory along PT)
#   Row 3: 5E (VL vs exhaustion) + 5F (naive branch composition)
################################################################################

library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
library(cowplot); library(patchwork); library(qs2)
library(SeuratExtend); library(ggrepel); library(ggpubr)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
fig5_dir     <- file.path(base_dir, "Manuscript", "Figure 5")
analysis_dir <- file.path(base_dir, "Analysis", "CD8_subset")

dir.create(fig5_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells\n")

# ── Shared aesthetics ────────────────────────────────────────────────────────
art_colors <- c(
  "PreART_Entry" = "#4A90D9", "PostART_Suppressed" = "#52B788",
  "PostART_Unsuppressed" = "#E76F51"
)
art_labels_vec <- c(
  "PreART_Entry" = "Pre-ART", "PostART_Suppressed" = "Suppressed",
  "PostART_Unsuppressed" = "Unsuppressed"
)

effector_clusters <- c("TEM CD8", "TEMRA CD8", "CD16+ Effector CD8")

if (!"PID" %in% colnames(TARA_cd8@meta.data)) {
  pid_vec <- sub("_.*$", "", TARA_cd8$orig.ident)
  names(pid_vec) <- colnames(TARA_cd8)
  TARA_cd8 <- AddMetaData(TARA_cd8, metadata = pid_vec, col.name = "PID")
}

# ── Data paths ───────────────────────────────────────────────────────────────
m3_pt_path  <- file.path(analysis_dir, "07_trajectory", "monocle3", "Monocle3_pseudotime_per_cell.csv")
module_path <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")
vl_path     <- file.path(analysis_dir, "09_viral_load", "PerSample_ModuleScores_and_ViralLoad.csv")

# ── Shared themes ────────────────────────────────────────────────────────────
pt_theme <- theme_cowplot(font_size = 18) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

vl_theme <- theme_cowplot(font_size = 18) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA))

################################################################################
# Panel A: Monocle3 trajectory UMAP (3 branches annotated)
################################################################################
cat("  Panel A: Monocle3 trajectory UMAP...\n")

if (file.exists(m3_pt_path)) {
  m3_pt <- read.csv(m3_pt_path)
  if (!"cell" %in% colnames(m3_pt)) m3_pt$cell <- m3_pt[, 1]
  
  umap_emb <- Embeddings(TARA_cd8, "wnn.umap")
  pt_umap_df <- data.frame(
    UMAP_1 = umap_emb[, 1], UMAP_2 = umap_emb[, 2],
    cell = colnames(TARA_cd8)
  ) %>% left_join(m3_pt %>% select(cell, monocle3_pseudotime), by = "cell")
  
  cluster_centroids <- pt_umap_df %>%
    left_join(data.frame(cell = colnames(TARA_cd8),
                         Cluster = TARA_cd8$CD8_Annotation), by = "cell") %>%
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
                        c("Naïve CD8", "Naïve CD8 2"), "Naive branch 1 (viremia)")
  path2 <- make_path_df(cluster_centroids,
                        c("Naïve CD8", "Naïve CD8 3"), "Naive branch 2 (ART)")
  path3 <- make_path_df(cluster_centroids,
                        c("Naïve CD8", "Naïve CD8 4", "TEM CD8", "TEMRA CD8",
                          "KIR+ innate-like CD8", "CD16+ Effector CD8"), "Effector trajectory")
  traj_paths <- bind_rows(path1, path2, path3)
  
  x_range_val <- diff(range(pt_umap_df$UMAP_1, na.rm = TRUE))
  y_range_val <- diff(range(pt_umap_df$UMAP_2, na.rm = TRUE))
  arrow_frac <- 0.12
  x_arrow_len <- x_range_val * arrow_frac
  y_arrow_len <- y_range_val * arrow_frac
  x_arrow <- min(pt_umap_df$UMAP_1, na.rm = TRUE) - x_range_val * 0.05
  y_arrow <- min(pt_umap_df$UMAP_2, na.rm = TRUE) - y_range_val * 0.08
  
  p_5A <- ggplot(pt_umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = monocle3_pseudotime), size = 0.5, alpha = 0.6) +
    geom_path(data = traj_paths, aes(x = UMAP_1, y = UMAP_2, group = Lineage),
              color = "black", linewidth = 1.8,
              arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
              inherit.aes = FALSE) +
    geom_point(data = traj_paths, aes(x = UMAP_1, y = UMAP_2),
               color = "black", size = 4, inherit.aes = FALSE) +
    scale_color_viridis_c(option = "inferno", name = "Pseudotime", na.value = "grey85") +
    annotate("segment", x = x_arrow, xend = x_arrow + x_arrow_len,
             y = y_arrow, yend = y_arrow,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"), linewidth = 0.8) +
    annotate("text", x = x_arrow + x_arrow_len/2, y = y_arrow - y_range_val * 0.04,
             label = "UMAP 1", size = 5, fontface = "bold") +
    annotate("segment", x = x_arrow, xend = x_arrow,
             y = y_arrow, yend = y_arrow + y_arrow_len,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"), linewidth = 0.8) +
    annotate("text", x = x_arrow - x_range_val * 0.04, y = y_arrow + y_arrow_len/2,
             label = "UMAP 2", size = 5, fontface = "bold", angle = 90) +
    labs(title = "Monocle3 pseudotime", x = NULL, y = NULL) +
    theme_void() +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16, face = "bold"),
          legend.key.height = unit(1.2, "cm"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(10, 10, 25, 25)) +
    coord_cartesian(clip = "off")
  
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
  if ("Timepoint" %in% colnames(pt_density_df) & !"Timepoint_Group" %in% colnames(pt_density_df)) {
    pt_density_df$Timepoint_Group <- pt_density_df$Timepoint
  }
  pt_density_df$Timepoint_Group <- factor(pt_density_df$Timepoint_Group, levels = names(art_colors))
  pt_density_df <- pt_density_df %>% filter(!is.na(Timepoint_Group))
  
  medians_pt <- pt_density_df %>%
    group_by(Timepoint_Group) %>%
    summarise(median_pt = median(monocle3_pseudotime, na.rm = TRUE), .groups = "drop")
  
  p_5B <- ggplot(pt_density_df, aes(x = monocle3_pseudotime, fill = Timepoint_Group,
                                    color = Timepoint_Group)) +
    geom_density(alpha = 0.35, linewidth = 1.2) +
    geom_vline(data = medians_pt, aes(xintercept = median_pt, color = Timepoint_Group),
               linewidth = 1, linetype = "dashed", show.legend = FALSE) +
    scale_fill_manual(values = art_colors, name = "ART Status", labels = art_labels_vec) +
    scale_color_manual(values = art_colors, name = "ART Status", labels = art_labels_vec) +
    labs(x = "Monocle3 pseudotime", y = "Density", title = "Pseudotime distribution") +
    theme_cowplot(font_size = 18) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA)) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.6)), color = "none")
  
  ggsave(file.path(fig5_dir, "Fig5B_Pseudotime_density.png"),
         plot = p_5B, width = 10, height = 7, dpi = 300, bg = "white")
}

################################################################################
# Panels C, D: Exhaustion + Type I IFN memory along pseudotime — SHARED LEGEND
################################################################################
cat("  Panels C, D: Exhaustion + IFN memory along pseudotime...\n")

if (exists("m3_pt") & file.exists(module_path)) {
  mod_scores <- read.csv(module_path, row.names = 1)
  common <- intersect(m3_pt$cell, rownames(mod_scores))
  
  pt_mod <- data.frame(
    pseudotime      = m3_pt$monocle3_pseudotime[match(common, m3_pt$cell)],
    Timepoint       = mod_scores[common, "Timepoint_Group"],
    Exhaustion      = mod_scores[common, "Exhaustion"]
  )
  if ("TypeI_IFN_Memory" %in% colnames(mod_scores)) {
    pt_mod$TypeI_IFN_Memory <- mod_scores[common, "TypeI_IFN_Memory"]
  }
  pt_mod <- pt_mod %>% filter(!is.na(pseudotime) & is.finite(pseudotime))
  pt_mod$Timepoint <- factor(pt_mod$Timepoint, levels = names(art_colors))
  
  p_5C <- ggplot(pt_mod, aes(x = pseudotime, y = Exhaustion, color = Timepoint)) +
    geom_point(size = 0.15, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
    scale_color_manual(values = art_colors) +
    labs(x = "Pseudotime", y = "Exhaustion score", title = "Exhaustion") +
    pt_theme
  
  p_5D <- ggplot(pt_mod, aes(x = pseudotime, y = TypeI_IFN_Memory, color = Timepoint)) +
    geom_point(size = 0.15, alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
    scale_color_manual(values = art_colors) +
    labs(x = "Pseudotime", y = "Type I IFN memory score", title = "Type I IFN memory") +
    pt_theme
  
  ggsave(file.path(fig5_dir, "Fig5C_Exhaustion_pseudotime.png"),
         plot = p_5C, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(fig5_dir, "Fig5D_IFN_memory_pseudotime.png"),
         plot = p_5D, width = 7, height = 6, dpi = 300, bg = "white")
  
  p_legend_source <- ggplot(pt_mod, aes(x = pseudotime, y = Exhaustion, color = Timepoint)) +
    geom_smooth(method = "loess", linewidth = 1.8) +
    scale_color_manual(values = art_colors, name = "ART Status", labels = art_labels_vec, drop = FALSE) +
    theme(legend.direction = "horizontal",
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20, face = "bold"),
          legend.key.size = unit(1.2, "cm")) +
    guides(color = guide_legend(override.aes = list(linewidth = 3), nrow = 1))
  
  legend_row2 <- cowplot::get_legend(p_legend_source)
  ggsave(file.path(fig5_dir, "Fig5_Row2_shared_legend.png"),
         plot = as_ggplot(legend_row2), width = 12, height = 2, dpi = 300, bg = "white")
  
  row2 <- (p_5C | p_5D) / as_ggplot(legend_row2) + plot_layout(heights = c(1, 0.08))
  ggsave(file.path(fig5_dir, "Fig5_Row2_CD_combined.png"),
         plot = row2, width = 14, height = 7.5, dpi = 300, bg = "white")
}

################################################################################
# Panel E: Exhaustion vs viral load
################################################################################
cat("  Panel E: Exhaustion vs viral load...\n")

if (file.exists(vl_path)) {
  sample_scores <- read.csv(vl_path)
  sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)
  sample_scores$Timepoint_Group <- factor(sample_scores$Timepoint_Group, levels = names(art_colors))
  
  cor_test <- cor.test(sample_scores$log10_VL, sample_scores$Exhaustion,
                       method = "spearman", exact = FALSE)
  
  p_5E <- ggplot(sample_scores, aes(x = log10_VL, y = Exhaustion)) +
    geom_point(aes(color = Timepoint_Group, size = n_expanding), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8, linetype = "dashed") +
    geom_text_repel(aes(label = PID), size = 4, max.overlaps = 15, color = "grey30") +
    scale_color_manual(values = art_colors) +
    scale_size_continuous(name = "N expanding", range = c(3, 9)) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("\u03c1 = %.2f\np = %.2g", cor_test$estimate, cor_test$p.value),
             hjust = -0.1, vjust = 1.3, size = 7, fontface = "bold") +
    labs(x = expression(log[10]~viral~load), y = "Mean exhaustion score",
         title = "Exhaustion vs viral load") +
    vl_theme + guides(size = "none")
  
  ggsave(file.path(fig5_dir, "Fig5E_Exhaustion_vs_VL.png"),
         plot = p_5E, width = 8, height = 7, dpi = 300, bg = "white")
}

################################################################################
# Panel F: Naive branch composition by ART status (NEW)
#   Fisher's exact enrichment tests saved to CSV for text/supplementary
################################################################################
cat("  Panel F: Naive branch composition by ART...\n")

branch_pops <- c("Naïve Intermediate CD8", "Tscm CD8", "Naïve CD8 2", "Naïve CD8 3")

# Per-cluster composition
branch_df <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% branch_pops) %>%
  group_by(CD8_Annotation, Timepoint_Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CD8_Annotation) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

branch_df$CD8_Annotation <- factor(branch_df$CD8_Annotation,
                                   levels = c("Naïve Intermediate CD8", "Tscm CD8",
                                              "Naïve CD8 2", "Naïve CD8 3"))
branch_df$Timepoint_Group <- factor(branch_df$Timepoint_Group, levels = names(art_colors))

branch_df$Label <- dplyr::recode(as.character(branch_df$CD8_Annotation),
                                 "Naïve Intermediate CD8" = "Intermediate",
                                 "Tscm CD8" = "Tscm",
                                 "Naïve CD8 2" = "Naïve CD8 2",
                                 "Naïve CD8 3" = "Naïve CD8 3")
branch_df$Label <- factor(branch_df$Label,
                          levels = c("Intermediate", "Tscm", "Naïve CD8 2", "Naïve CD8 3"))

branch_df$pct_label <- ifelse(branch_df$pct >= 5, paste0(round(branch_df$pct, 0), "%"), "")

# Fisher's exact enrichment tests (saved for text/supplementary)
fisher_results <- list()
expected_pct <- prop.table(table(TARA_cd8$Timepoint_Group))
for (pop in branch_pops) {
  pop_cells <- TARA_cd8@meta.data %>% filter(CD8_Annotation == pop)
  n_pop <- nrow(pop_cells)
  for (cond in names(art_colors)) {
    obs_in <- sum(pop_cells$Timepoint_Group == cond)
    ft <- fisher.test(matrix(c(obs_in, n_pop - obs_in,
                               table(TARA_cd8$Timepoint_Group)[cond] - obs_in,
                               ncol(TARA_cd8) - n_pop - (table(TARA_cd8$Timepoint_Group)[cond] - obs_in)),
                             nrow = 2))
    obs_pct <- obs_in / n_pop * 100
    exp_pct <- as.numeric(expected_pct[cond]) * 100
    fisher_results[[length(fisher_results) + 1]] <- data.frame(
      Cluster = pop, Condition = cond,
      obs_pct = round(obs_pct, 1), exp_pct = round(exp_pct, 1),
      fold_enrichment = round(obs_pct / exp_pct, 2),
      fisher_p = ft$p.value,
      star = ifelse(ft$p.value < 0.001, "***",
                    ifelse(ft$p.value < 0.01, "**",
                           ifelse(ft$p.value < 0.05, "*", "ns")))
    )
  }
}
fisher_df <- bind_rows(fisher_results)
write.csv(fisher_df, file.path(fig5_dir, "Fig5F_enrichment_fisher_tests.csv"), row.names = FALSE)
cat("  Fisher's exact enrichment tests saved.\n")

p_5F <- ggplot(branch_df, aes(x = Label, y = pct, fill = Timepoint_Group)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = pct_label),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 5.5) +
  # Dashed separator between branch points and endpoints
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values = art_colors, name = "ART Status", labels = art_labels_vec) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  # Bracket labels above bars
  annotate("text", x = 1.5, y = 102, label = "Branch point",
           size = 4.5, fontface = "italic", color = "grey30") +
  annotate("text", x = 3.5, y = 102, label = "Branch endpoints",
           size = 4.5, fontface = "italic", color = "grey30") +
  labs(x = "", y = "% of cluster",
       title = "Naïve branch composition by ART status") +
  theme_cowplot(font_size = 18) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA)) +
  coord_cartesian(clip = "off", ylim = c(0, 100))

ggsave(file.path(fig5_dir, "Fig5F_Naive_branch_composition.png"),
       plot = p_5F, width = 10, height = 7, dpi = 300, bg = "white")

# Combined Row 3
row3 <- (p_5E | p_5F)
ggsave(file.path(fig5_dir, "Fig5_Row3_EF_combined.png"),
       plot = row3, width = 18, height = 8, dpi = 300, bg = "white")

cat("\n=== Figure 5 complete ===\n")
cat("Output:", fig5_dir, "\n")