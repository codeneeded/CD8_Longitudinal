################################################################################
# FIGURE 2: TCR Clonal Expansion Analysis — TARA Cohort
#
# Panel A: (i)  DimPlot — clonal expansion split by HEI / HEU / HUU
#          (ii) Stacked bar — cloneSize per CD8 cluster
# Panel B: Alluvial plots — longitudinal clonal tracking per patient
# Panel C: Waffle charts — Trex epitope specificity per HEI patient
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
library(Trex)

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

TARA_ALL$Cluster_Number_refined <- gsub("^([0-9]+):.*", "\\1",
                                        TARA_ALL$Manual_Annotation_refined)

cat("=== Post-merge cluster sizes (Clusters 8 & 9) ===\n")
print(table(TARA_ALL$Manual_Annotation_refined)[
  grepl("^8|^9", names(table(TARA_ALL$Manual_Annotation_refined)))
])

################################################################################
# 2. UPDATE ANNOTATION DISPLAY NAMES
################################################################################
cd8_display_labels <- c(
  "1: Memory CD8 T cell"   = "Tcm/Tscm CD8",
  "6: Naïve CD8 T cell"    = "Naïve CD8",
  "8: CTL-like"            = "TEMRA/CTL",
  "9: TRDV1+ γδ T cell"    = "TRDV1+ γδ",
  "27: GZMK+ CD8 T cell"   = "Tpex CD8"
)

current_levels <- levels(TARA_ALL$Manual_Annotation_refined)
new_levels     <- current_levels

for (i in seq_along(current_levels)) {
  lbl <- current_levels[i]
  if (lbl %in% names(cd8_display_labels)) {
    new_levels[i] <- cd8_display_labels[lbl]
  } else {
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
cat("Saved refined object.\n")

################################################################################
# FIGURE 2A(i) — UMAP: Clonal expansion split by HEI / HEU / HUU
################################################################################
colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

fig2a_umap <- DimPlot_scCustom(
  TARA_ALL,
  group.by     = "cloneSize",
  reduction    = "wnn.umap",
  split.by     = "Condition",
  split_seurat = TRUE,
  pt.size      = 0.4
) +
  scale_color_manual(values = c(colorblind_vector[c(1, 3, 4, 5, 7)])) +
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

ggsave(file.path(out_dir, "Fig2A_i_ClonalExpansion_UMAP_byCondition.png"),
       plot = fig2a_umap, width = 18, height = 6, dpi = 300, bg = "white")
cat("Fig 2A(i) saved.\n")

################################################################################
# FIGURE 2A(ii) — Stacked bar: cloneSize per CD8 cluster
################################################################################
clone_order <- c("Single (0 < X <= 1)",
                 "Small (1 < X <= 5)",
                 "Medium (5 < X <= 20)",
                 "Large (20 < X <= 100)",
                 "Hyperexpanded (100 < X <= 500)")

clone_colors <- setNames(c(colorblind_vector[c(1, 3, 4, 5, 7)]), clone_order)

cd8_labels <- c("Tcm/Tscm CD8", "Naïve CD8", "TEMRA/CTL",
                "TRDV1+ γδ", "Tpex CD8")

TARA_CD8 <- subset(TARA_ALL, Manual_Annotation_refined %in% cd8_labels)
TARA_CD8$Manual_Annotation_refined <- droplevels(TARA_CD8$Manual_Annotation_refined)

meta_cd8 <- TARA_CD8@meta.data %>%
  mutate(Cluster = Manual_Annotation_refined,
         cloneSize = factor(cloneSize, levels = clone_order)) %>%
  filter(!is.na(cloneSize))

bar_data <- meta_cd8 %>%
  group_by(Cluster, cloneSize, .drop = FALSE) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

bar_data$Cluster <- factor(bar_data$Cluster,
                           levels = c("Tcm/Tscm CD8", "Naïve CD8", "TEMRA/CTL",
                                      "TRDV1+ γδ", "Tpex CD8"))

bar_data <- bar_data %>%
  arrange(Cluster, cloneSize) %>%
  group_by(Cluster) %>%
  mutate(total = sum(n), cum_n = cumsum(n), y_mid = cum_n - n / 2) %>%
  ungroup()

min_height <- 50   # ← EDIT THIS NUMBER to tune label visibility
bar_data <- bar_data %>%
  mutate(label = ifelse(n >= min_height, as.character(n), ""))

fig2a_bar <- ggplot(bar_data, aes(x = Cluster, y = n, fill = cloneSize)) +
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

ggsave(file.path(out_dir, "Fig2A_ii_ClonalOccupy_CD8clusters.png"),
       plot = fig2a_bar, width = 10, height = 7, dpi = 300, bg = "white")
cat("Fig 2A(ii) saved.\n")

################################################################################
# FIGURE 2B — Alluvial Plots: Longitudinal clonal tracking per patient
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

TARA.TCR             <- paste0(TARA_names, ".TCR")
TARA.contig_list.TCR <- as.list(mget(TARA.TCR))
combined.TCR.TARA    <- combineTCR(TARA.contig_list.TCR, samples = TARA_names)

# ── Patient alluvial plots ───────────────────────────────────────────────────
patient_info <- list(
  CP003 = c("CP003_entry", "CP003_V12", "CP003_V24"),
  CP006 = c("CP006_entry", "CP006_12m", "CP006_V24"),
  CP013 = c("CP013_1m",    "CP013_12m", "CP013_24m"),
  CP018 = c("CP018_entry", "CP018_V24", "CP018_42m"),
  CP020 = c("CP020_V1",    "CP020_V12", "CP020_V44")
)

timepoint_labels <- c(
  "CP003_entry" = "1-2m pre-ART", "CP003_V12" = "12m",  "CP003_V24" = "24m",
  "CP006_entry" = "1-2m pre-ART", "CP006_12m" = "12m",  "CP006_V24" = "24m",
  "CP013_1m"    = "1-2m pre-ART", "CP013_12m" = "12m",  "CP013_24m" = "24m",
  "CP018_entry" = "1-2m pre-ART", "CP018_V24" = "24m",  "CP018_42m" = "42m",
  "CP020_V1"    = "1-2m pre-ART", "CP020_V12" = "12m",  "CP020_V44" = "44m"
)

alluvial_plots <- list()

for (pt in names(patient_info)) {
  samps <- patient_info[[pt]]
  display_labels <- timepoint_labels[samps]
  
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
  
  ggsave(file.path(out_dir, paste0("Fig2B_Alluvial_", pt, ".png")),
         plot = p, width = 15, height = 11, dpi = 300, bg = "white")
}

# Combined with axis labels
y_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Count", size = 16,
           angle = 90, fontface = "bold") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

x_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Timepoint", size = 16,
           fontface = "bold") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

alluvial_grid <- wrap_plots(alluvial_plots, ncol = 3)

fig2b_combined <- (
  (y_label | alluvial_grid) + plot_layout(widths = c(0.03, 1))
) / (
  (plot_spacer() | x_label) + plot_layout(widths = c(0.03, 1))
) + plot_layout(heights = c(1, 0.04))

ggsave(file.path(out_dir, "Fig2B_Alluvial_AllPatients.png"),
       plot = fig2b_combined, width = 36, height = 24, dpi = 300, bg = "white")
cat("Fig 2B saved (individual + combined).\n")

################################################################################
# FIGURE 2C — Waffle Charts: Epitope specificity (HEI only, ED2)
#              1 square = 1 unique expanded clone. All grids same size.
################################################################################

# ── Run Trex (ED2) ───────────────────────────────────────────────────────────
TARA_ALL_TRB_2 <- annotateDB(TARA_ALL, chains = "TRB", edit.distance = 2)

trex_meta <- TARA_ALL_TRB_2@meta.data %>%
  filter(!is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(
    PID = sub("_.*$", "", orig.ident),
    Epitope_Species = ifelse(is.na(TRB_Epitope.species), "Unknown", TRB_Epitope.species)
  )

# ── Simplify species ─────────────────────────────────────────────────────────
waffle_data <- trex_meta %>%
  distinct(PID, CTstrict, Epitope_Species) %>%
  mutate(Epitope_Simple = sapply(strsplit(Epitope_Species, ";"), `[`, 1)) %>%
  mutate(Epitope_Clean = case_when(
    Epitope_Simple == "Unknown"                                                    ~ "Unknown",
    grepl("HIV", Epitope_Simple, ignore.case = TRUE)                               ~ "HIV",
    grepl("CMV|Cytomegalovirus", Epitope_Simple, ignore.case = TRUE)               ~ "CMV",
    grepl("EBV|Epstein", Epitope_Simple, ignore.case = TRUE)                       ~ "EBV",
    grepl("Influenza", Epitope_Simple, ignore.case = TRUE)                         ~ "Influenza",
    grepl("tuberculosis|M\\.tuberculosis|Mtb", Epitope_Simple, ignore.case = TRUE) ~ "TB",
    TRUE ~ "Other"
  ))

# ── Palette ──────────────────────────────────────────────────────────────────
epitope_colors <- c(
  "HIV"       = "#E41A1C",
  "CMV"       = "#377EB8",
  "EBV"       = "#4DAF4A",
  "Influenza" = "#FF7F00",
  "TB"        = "#984EA3",
  "Other"     = "#A65628",
  "Unknown"   = "#D9D9D9"
)
epitope_order <- c("HIV", "CMV", "EBV", "Influenza", "TB", "Other", "Unknown")

# ── Counts ───────────────────────────────────────────────────────────────────
waffle_counts <- waffle_data %>%
  group_by(PID, Epitope_Clean) %>%
  summarise(n_clones = n(), .groups = "drop") %>%
  mutate(Epitope_Clean = factor(Epitope_Clean, levels = epitope_order))

# ── Uniform grid: find max clones, pad all patients to same size ─────────────
hei_patients <- c("CP003", "CP006", "CP013", "CP018", "CP020")
ncol_waffle  <- 10

max_clones <- waffle_counts %>%
  filter(PID %in% hei_patients) %>%
  group_by(PID) %>%
  summarise(total = sum(n_clones)) %>%
  pull(total) %>%
  max()

max_squares <- ceiling(max_clones / ncol_waffle) * ncol_waffle
nrow_waffle <- max_squares / ncol_waffle
cat("Waffle grid:", ncol_waffle, "cols x", nrow_waffle, "rows =", max_squares, "squares\n")

# ── Helper: make tile df padded to uniform grid ──────────────────────────────
make_waffle_df <- function(counts, ncol, total_squares) {
  categories <- rep(names(counts), times = counts)
  n_pad      <- total_squares - length(categories)
  if (n_pad > 0) categories <- c(categories, rep("_empty", n_pad))
  n <- length(categories)
  data.frame(
    category = factor(categories, levels = c(epitope_order, "_empty")),
    x = ((seq_len(n) - 1) %% ncol) + 1,
    y = ((seq_len(n) - 1) %/% ncol) + 1
  )
}

epitope_colors_ext <- c(epitope_colors, "_empty" = NA)

# ── Build per-patient waffles ────────────────────────────────────────────────
waffle_plots <- list()

for (pid in hei_patients) {
  
  pid_data <- waffle_counts %>%
    filter(PID == pid) %>%
    arrange(Epitope_Clean)
  
  waf_vec      <- setNames(pid_data$n_clones, as.character(pid_data$Epitope_Clean))
  total_clones <- sum(waf_vec)
  tile_df      <- make_waffle_df(waf_vec, ncol = ncol_waffle,
                                 total_squares = max_squares)
  
  p <- ggplot(tile_df, aes(x = x, y = y, fill = category)) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_fill_manual(values = epitope_colors_ext, drop = FALSE,
                      na.value = "white", guide = "none") +
    scale_y_reverse() +
    coord_equal(xlim = c(0.5, ncol_waffle + 0.5),
                ylim = c(nrow_waffle + 0.5, 0.5)) +
    labs(title = pid, subtitle = paste0("n = ", total_clones, " clones")) +
    theme_void() +
    theme(
      plot.title       = element_text(size = 32, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 24, hjust = 0.5, color = "grey40",
                                      margin = margin(t = 4)),
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin      = margin(8, 8, 8, 8)
    )
  
  waffle_plots[[pid]] <- p
  
  # Individual PNG
  ggsave(file.path(out_dir, paste0("Fig2C_Waffle_", pid, ".png")),
         plot = p, width = 7, height = 7, dpi = 300, bg = "white")
}

# ── Shared legend panel ──────────────────────────────────────────────────────
# ── Shared legend panel — matched to waffle grid width ───────────────────────
# Spread 7 legend items evenly across positions 1–5 (matching 5 waffle columns)
n_items <- length(epitope_order)
legend_df <- data.frame(
  Epitope = factor(epitope_order, levels = epitope_order),
  x = seq(from = 1, to = 5, length.out = n_items),
  y = 1
)

legend_panel <- ggplot(legend_df, aes(x = x, y = y, fill = Epitope)) +
  geom_point(shape = 22, size = 14, color = "grey30", stroke = 0.5) +
  geom_text(aes(label = Epitope), vjust = 2.8, size = 8, fontface = "bold") +
  scale_fill_manual(values = epitope_colors) +
  scale_x_continuous(limits = c(-1, 7), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.2, 1.8), clip = "off") +
  theme_void() +
  theme(
    legend.position  = "none",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(5, 8, 20, 8)
  )

# ── Combined figure: all 5 in one row ────────────────────────────────────────
waffle_grid <- wrap_plots(waffle_plots, ncol = 5)

fig2c <- waffle_grid / legend_panel +
  plot_layout(heights = c(1, 0.15))

ggsave(file.path(out_dir, "Fig2C_Waffle_EpitopeSpecificity_HEI.png"),
       plot = fig2c, width = 40, height = 10, dpi = 300, bg = "white")

cat("Fig 2C saved (individual + combined).\n")

################################################################################
# SESSION INFO
################################################################################
cat("\n=== Done: All Figure 2 panels generated ===\n")
cat("Output directory:", out_dir, "\n")
sessionInfo()