############################################################
# Figure 3 — Panel Generation Script (revised)
# Functional validation of HIV-specific CD8+ T cell clonotypes
#
# PANELS:
#   A — Experimental schematic (Illustrator — not generated here)
#   B — Identification of HIV-specific clones in CP003 stimulation
#       B1: UMAP colored by cell cycle phase
#       B2: UMAP colored by HIV-Specific TCR status
#       B3: Cell cycle phase proportion by cluster (stacked bar)
#   C — Alluvial: HIV-specific clones tracked across all timepoints
#   D — Summary: persistence histogram + cluster composition in TARA
#
# OUTPUT: Fig3/panels/ (PNG only)
############################################################

library(Seurat)
library(SeuratExtend)
library(scRepertoire)
library(scCustomize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(RColorBrewer)
library(qs2)

# =============================================================
# 0) PATHS
# =============================================================
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
fig3_dir  <- file.path(base_dir, "Manuscript/Fig3")
ana_dir   <- file.path(fig3_dir, "analysis")
panel_dir <- file.path(fig3_dir, "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1) LOAD DATA
# =============================================================
cat("Loading objects...\n")
cp003      <- qs_read(file.path(saved_dir, "cp003_with_cc_epitope.qs2"))
tara_cp003 <- qs_read(file.path(saved_dir, "tara_cp003_fig3_annotated.qs2"))

# Analysis CSVs
alluvial_df    <- read.csv(file.path(ana_dir, "Alluvial_HIV_clones_combined_timepoints.csv"))
persistence_df <- read.csv(file.path(ana_dir, "HIV_clone_total_persistence_summary.csv"))
cluster_status <- read.csv(file.path(ana_dir, "TARA_CP003_cluster_by_match_status.csv"))

# =============================================================
# PANEL B: Identification of HIV-specific clones in CP003
# Shows HOW we identified them: cell cycle + clonal expansion
# =============================================================
cat("Generating Panel B...\n")

red_use <- "umap.mnn.rna"

# --- B1: UMAP colored by cluster number (reference for B3 cell cycle plot) ---
panelB1 <- DimPlot2(
  cp003,
  reduction = red_use,
  group.by = "mnn_clusters_rna",
  cols = "default",
  label = TRUE,
  box = TRUE,
  label.size = 10,
  repel = TRUE,
  pt.size = 0.4,
  raster = FALSE,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("CP003 clusters") +
  theme(plot.title = element_text(size = 18, face = "bold"))

ggsave(file.path(panel_dir, "Fig3B_UMAP_Clusters.png"),
       panelB1, width = 6, height = 6, dpi = 400, bg = "white")

# --- B2: UMAP colored by HIV-Specific TCR status ---
hiv_stim_cols <- c("Other" = "grey85", "HIV-Specific TCR" = "#E63946")

panelB2 <- DimPlot2(
  cp003,
  features = "HIV_Specific_TCR",
  reduction = red_use,
  cols = hiv_stim_cols,
  pt.size = 0.4,
  raster = FALSE,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("HIV-specific clonal expansion") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.6, "cm")
  )

ggsave(file.path(panel_dir, "Fig3D_UMAP_HIV_Specific.png"),
       panelB2, width = 6, height = 6, dpi = 400, bg = "white")

# --- B3: Cell cycle phase proportion by cluster (stacked bar) ---
phase_cols <- c("G1" = "#AEC6CF", "S" = "#FFD700", "G2M" = "#E63946")

cc_cluster <- as.data.frame(cp003@meta.data) %>%
  filter(!is.na(Phase), !is.na(mnn_snn_res.0.6)) %>%
  group_by(mnn_snn_res.0.6, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mnn_snn_res.0.6) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

panelB3 <- ggplot(cc_cluster,
                  aes(x = as.factor(mnn_snn_res.0.6),
                      y = proportion,
                      fill = Phase)) +
  geom_col() +
  scale_fill_manual(values = phase_cols) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Cluster",
    y = "Proportion of cells",
    fill = "Cell cycle phase",
    title = "Cell cycle by cluster"
  )

ggsave(file.path(panel_dir, "Fig3C_CellCycle_by_Cluster.png"),
       panelB3, width = 10, height = 6, dpi = 400, bg = "white")

# =============================================================
# PANEL E: Alluvial — HIV-specific clones across all timepoints
# =============================================================
cat("Generating Panel C...\n")

# Relabel timepoints: both entry and 2m stim are "2 months" but stim gets subtitle
alluvial_df$timepoint_label <- dplyr::recode(
  alluvial_df$timepoint_label,
  "Entry (~1-2m)" = "2 months",
  "2m (stim)"     = "2 months\n(HIV-stim)",
  "12m"           = "12 months",
  "24m"           = "24 months",
  "101m (stim)"   = "101 months\n(HIV-stim)"
)

alluvial_df$timepoint_label <- factor(
  alluvial_df$timepoint_label,
  levels = c("2 months\n(HIV-stim)", "2 months", "12 months", "24 months", "101 months\n(HIV-stim)")
)

# Keep clones at 2+ timepoints for cleaner alluvial
multi_tp_clones <- alluvial_df %>%
  group_by(CTstrict) %>%
  filter(n_distinct(timepoint_label) >= 2) %>%
  ungroup()

if (nrow(multi_tp_clones) > 0) {
  
  # Assign colors by total abundance
  clone_order <- multi_tp_clones %>%
    group_by(CTstrict) %>%
    summarise(total = sum(n_cells), .groups = "drop") %>%
    arrange(desc(total))
  
  n_clones <- nrow(clone_order)
  clone_pal <- colorRampPalette(brewer.pal(min(n_clones, 12), "Set3"))(n_clones)
  names(clone_pal) <- clone_order$CTstrict
  
  panelC <- ggplot(multi_tp_clones,
                   aes(x = timepoint_label, y = n_cells,
                       stratum = CTstrict, alluvium = CTstrict,
                       fill = CTstrict)) +
    geom_stratum(width = 0.3, color = "grey30", linewidth = 0.3) +
    geom_flow(alpha = 0.5, width = 0.3) +
    scale_fill_manual(values = clone_pal, guide = "none") +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(size = 17),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 14)
    ) +
    labs(
      x = NULL,
      y = "Number of cells",
      title = "HIV-specific clonotype persistence across timepoints",
      subtitle = paste0(n_clones, " clonotypes detected at \u22652 timepoints")
    )
  
  ggsave(file.path(panel_dir, "Fig3E_alluvial_HIV_clones.png"),
         panelC, width = 11, height = 7, dpi = 400, bg = "white")
} else {
  cat("WARNING: No multi-timepoint clones found for alluvial.\n")
}

# =============================================================
# PANEL D: Clone persistence histogram
# =============================================================
cat("Generating Panel D...\n")

persist_summary <- persistence_df %>%
  dplyr::count(n_timepoints) %>%
  mutate(n_timepoints = factor(n_timepoints))

panelD <- ggplot(persist_summary,
                 aes(x = n_timepoints, y = n)) +
  geom_col(fill = "#E63946", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Number of timepoints detected",
    y = "Number of clonotypes",
    title = "HIV-specific clonotype persistence"
  )

ggsave(file.path(panel_dir, "Fig3F_persistence_histogram.png"),
       panelD, width = 7, height = 6, dpi = 400, bg = "white")

# =============================================================
# DONE
# =============================================================
cat("\n=== All panels saved to:", panel_dir, "===\n")
cat("Panel A: Experimental schematic (Illustrator)\n")
cat("Panel B: CP003 cluster UMAP (Fig3B_UMAP_Clusters.png)\n")
cat("Panel C: Cell cycle by cluster (Fig3C_CellCycle_by_Cluster.png)\n")
cat("Panel D: HIV-specific UMAP (Fig3D_UMAP_HIV_Specific.png)\n")
cat("Panel E: Alluvial tracking (Fig3E_alluvial_HIV_clones.png)\n")
cat("Panel F: Persistence histogram (Fig3F_persistence_histogram.png)\n")

sessionInfo()