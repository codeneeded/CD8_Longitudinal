################################################################################
# FIGURE 2: TCR Clonal Expansion Analysis — TARA Cohort
#
# Panel A: DimPlot of clonally expanding cells split by HEI / HEU / HUU
# Panel B: Bar plot — cluster composition of expanded clones by cloneSize
# Panel C: [Placeholder — reserved for future panel]
# Panel D: Alluvial plots (clonalCompare) for 5 patients across timepoints
#
# Input:  TARA_ALL_sorted.qs2  +  combined.TCR.TARA (built from VDJ contigs)
# Output: /home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 2/
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(scCustomize)
library(scRepertoire)
library(RColorBrewer)
library(Polychrome)
library(qs2)

# ── Paths ─────────────────────────────────────────────────────────────────────
saved_dir  <- "~/Documents/CD8_Longitudinal/saved_R_data/"
out_dir    <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 2/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Data ─────────────────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted.qs2"))

################################################################################
# 1. CLEANING: Merge αβ TCR+ contaminants from Cluster 9 → Cluster 8
################################################################################
TARA_ALL$Manual_Annotation_refined <- case_when(
  TARA_ALL$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_ALL$has_TCR == FALSE
  ~ "9: TRDV1+ γδ T cell",
  TARA_ALL$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_ALL$has_TCR == TRUE
  ~ "8: CTL-like",
  TRUE
  ~ as.character(TARA_ALL$Manual_Annotation)
)
TARA_ALL$Manual_Annotation_refined <- factor(TARA_ALL$Manual_Annotation_refined)

# Numeric cluster label (for UMAP)
TARA_ALL$Cluster_Number_refined <- gsub("^([0-9]+):.*", "\\1",
                                        TARA_ALL$Manual_Annotation_refined)

# Verify
cat("=== Post-merge cluster sizes (Clusters 8 & 9) ===\n")
print(table(TARA_ALL$Manual_Annotation_refined)[
  grepl("^8|^9", names(table(TARA_ALL$Manual_Annotation_refined)))
])

################################################################################
# 2. UPDATE ANNOTATION DISPLAY NAMES
################################################################################

# Specific CD8 cluster renames
cd8_display_labels <- c(
  "1: Memory CD8 T cell"   = "Tcm/Tscm CD8",
  "6: Naïve CD8 T cell"    = "Naïve CD8",
  "8: CTL-like"            = "TEMRA/CTL",
  "9: TRDV1+ γδ T cell"    = "TRDV1+ γδ",
  "27: GZMK+ CD8 T cell"   = "Tpex CD8"
)

# For all OTHER clusters: strip "N: " prefix, keep annotation name only
current_levels <- levels(TARA_ALL$Manual_Annotation_refined)
new_levels     <- current_levels  # start with a copy

for (i in seq_along(current_levels)) {
  lbl <- current_levels[i]
  if (lbl %in% names(cd8_display_labels)) {
    # Apply specific CD8 display name
    new_levels[i] <- cd8_display_labels[lbl]
  } else {
    # Strip cluster number prefix ("N: ") for all other clusters
    new_levels[i] <- gsub("^[0-9]+:\\s*", "", lbl)
  }
}

levels(TARA_ALL$Manual_Annotation_refined) <- new_levels

cat("\n=== Updated annotation levels ===\n")
print(levels(TARA_ALL$Manual_Annotation_refined))

################################################################################
# 3. SAVE CLEANED OBJECT
################################################################################
qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_sorted_refined.qs2"))
cat("Saved refined object to:", file.path(saved_dir, "TARA_ALL_sorted_refined.qs2"), "\n")

################################################################################
# FIGURE 2A — UMAP: Clonal expansion split by Condition (HEI / HEU / HUU)
#              Uses DimPlot_scCustom; grey cells = no TCR detected
################################################################################

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

fig2a <- DimPlot_scCustom(
  TARA_ALL,
  group.by     = "cloneSize",
  reduction    = "wnn.umap",
  split.by     = "Condition",
  split_seurat = TRUE,
  pt.size      = 0.4
) +
  scale_color_manual(
    values  = c(colorblind_vector[c(1, 3, 4, 5, 7)])
  ) +
  theme(
    plot.title       = element_blank(),
    strip.text       = element_text(size = 14, face = "bold"),
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.position  = "right",
    axis.title       = element_text(size = 12),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "Fig2A_ClonalExpansion_UMAP_byCondition.png"),
       plot = fig2a, width = 18, height = 6, dpi = 300, bg = "white")

cat("Fig 2A saved.\n")

################################################################################
# FIGURE 2B — Stacked bar plot: Clonal expansion categories per CD8 cluster
#              Same plasma palette as Fig 2A; includes γδ to show 0 clones
#              Cell counts shown inside bars (omitted if segment too small)
################################################################################

# ── Define cloneSize levels & matching plasma colors (same as 2A) ─────────────
clone_order <- c("Single (0 < X <= 1)",
                 "Small (1 < X <= 5)",
                 "Medium (5 < X <= 20)",
                 "Large (20 < X <= 100)",
                 "Hyperexpanded (100 < X <= 500)")

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)
clone_colors <- setNames(
  c(colorblind_vector[c(1, 3, 4, 5, 7)]),
  clone_order
)

# ── Subset to 5 CD8 clusters ─────────────────────────────────────────────────
cd8_labels <- c("Tcm/Tscm CD8", "Naïve CD8", "TEMRA/CTL",
                "TRDV1+ γδ", "Tpex CD8")

TARA_CD8 <- subset(TARA_ALL, Manual_Annotation_refined %in% cd8_labels)
TARA_CD8$Manual_Annotation_refined <- droplevels(TARA_CD8$Manual_Annotation_refined)

# ── Build bar data: count cells per cluster × cloneSize ──────────────────────
meta_cd8 <- TARA_CD8@meta.data %>%
  mutate(
    Cluster   = Manual_Annotation_refined,
    cloneSize = factor(cloneSize, levels = clone_order)
  ) %>%
  filter(!is.na(cloneSize))

bar_data <- meta_cd8 %>%
  group_by(Cluster, cloneSize, .drop = FALSE) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ── Cluster order (biological) ────────────────────────────────────────────────
bar_data$Cluster <- factor(bar_data$Cluster,
                           levels = c("Tcm/Tscm CD8", "Naïve CD8", "TEMRA/CTL",
                                      "TRDV1+ γδ", "Tpex CD8"))

# ── Compute label positions using raw counts (Single at bottom visually) ──────
# position_stack(reverse = TRUE) stacks factor levels in order: Single at bottom,
# Hyperexpanded at top. Label cumsum must match that same bottom-to-top order.
bar_data <- bar_data %>%
  arrange(Cluster, cloneSize) %>%
  group_by(Cluster) %>%
  mutate(
    total     = sum(n),
    cum_n     = cumsum(n),
    y_mid     = cum_n - n / 2
  ) %>%
  ungroup()

# Raw cell count label — only if the segment is tall enough in absolute terms
# Adjust min_height to control: this is in the same units as the y-axis (cell count)
# e.g. min_height = 50 means segments shorter than 50 cells get no label
min_height <- 50   # ← EDIT THIS NUMBER to tune which segments get labels
bar_data <- bar_data %>%
  mutate(label = ifelse(n >= min_height, as.character(n), ""))

# ── Plot ─────────────────────────────────────────────────────────────────────
fig2b <- ggplot(bar_data, aes(x = Cluster, y = n, fill = cloneSize)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(aes(y = y_mid, label = label),
            size = 4.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = clone_colors, name = "Clone Size", drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of Cells") +
  theme_cowplot(font_size = 14) +
  theme(
    axis.text.x       = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y       = element_text(size = 13),
    axis.title.y      = element_text(size = 15, face = "bold"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 13),
    legend.key.size   = unit(0.6, "cm"),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(out_dir, "Fig2B_ClonalOccupy_CD8clusters.png"),
       plot = fig2b, width = 10, height = 7, dpi = 300, bg = "white")

cat("Fig 2B saved.\n")

################################################################################
# FIGURE 2C — [PLACEHOLDER]
################################################################################
cat("Fig 2C: Placeholder — panel reserved for future content.\n")

################################################################################
# FIGURE 2D — Alluvial Plots: Clonal tracking across timepoints (per patient)
#              Requires combined.TCR.TARA from scRepertoire
################################################################################

# ── Rebuild combined.TCR.TARA ────────────────────────────────────────────────
in.path <- "~/Documents/10x_Genomics/"

get_folder_names <- function(in.path) {
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  folder_names <- folder_names[folder_names != ""]
  return(folder_names)
}

f_names    <- get_folder_names(in.path)
TARA_names <- f_names[grepl("^C", f_names)]

for (name in TARA_names) {
  t_file <- paste0(in.path, name, "/per_sample_outs/TCR/filtered_contig_annotations.csv")
  assign(paste0(name, ".TCR"), read.csv(t_file))
}

TARA.TCR            <- paste0(TARA_names, ".TCR")
TARA.contig_list.TCR <- as.list(mget(TARA.TCR))
combined.TCR.TARA   <- combineTCR(TARA.contig_list.TCR, samples = TARA_names)

# ── Patient-specific alluvial plots ──────────────────────────────────────────

# Define patient info: sample names per patient (in temporal order)
patient_info <- list(
  CP003 = c("CP003_entry", "CP003_V12", "CP003_V24"),
  CP006 = c("CP006_entry", "CP006_12m", "CP006_V24"),
  CP013 = c("CP013_1m",    "CP013_12m", "CP013_24m"),
  CP018 = c("CP018_entry", "CP018_V24", "CP018_42m"),
  CP020 = c("CP020_V1",    "CP020_V12", "CP020_V44")
)

# Unified display labels for each sample's timepoint
timepoint_labels <- c(
  "CP003_entry" = "1-2m pre-ART",  "CP003_V12" = "12m",   "CP003_V24" = "24m",
  "CP006_entry" = "1-2m pre-ART",  "CP006_12m" = "12m",   "CP006_V24" = "24m",
  "CP013_1m"    = "1-2m pre-ART",  "CP013_12m" = "12m",   "CP013_24m" = "24m",
  "CP018_entry" = "1-2m pre-ART",  "CP018_V24" = "24m",   "CP018_42m" = "42m",
  "CP020_V1"    = "1-2m pre-ART",  "CP020_V12" = "12m",   "CP020_V44" = "44m"
)

alluvial_plots <- list()

for (pt in names(patient_info)) {
  samps <- patient_info[[pt]]
  display_labels <- timepoint_labels[samps]
  
  # Plot
  p <- clonalCompare(
    combined.TCR.TARA,
    top.clones     = 20,
    samples        = samps,
    order.by       = samps,
    cloneCall      = "strict",
    relabel.clones = TRUE,
    proportion     = FALSE,
    graph          = "alluvial"
  ) +
    scale_x_discrete(labels = display_labels) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      axis.text.x      = element_text(size = 48, face = "bold"),
      axis.text.y      = element_text(size = 44, face = "bold"),
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  alluvial_plots[[pt]] <- p
  
  # Save individual PNG
  ggsave(file.path(out_dir, paste0("Fig2D_Alluvial_", pt, ".png")),
         plot = p, width = 15, height = 11, dpi = 300, bg = "white")
}

# Combined multi-panel alluvial figure
fig2d_combined <- wrap_plots(alluvial_plots, ncol = 3) +
  plot_annotation(
    title = "Fig 2D: Longitudinal Clonal Tracking by Patient",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(out_dir, "Fig2D_Alluvial_AllPatients.png"),
       plot = fig2d_combined, width = 36, height = 22, dpi = 300, bg = "white")

cat("Fig 2D saved (individual + combined).\n")

################################################################################
# SESSION INFO
################################################################################
cat("\n=== Done: All Figure 2 panels generated ===\n")
cat("Output directory:", out_dir, "\n")
sessionInfo()