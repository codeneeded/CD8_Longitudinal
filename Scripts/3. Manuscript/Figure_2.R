################################################################################
# FIGURE 2 — TCR Clonal Expansion Analysis (v5 — Updated Annotations)
#
# Updated to match TARA Unified WNN Annotation Pipeline (Final)
#
# KEY CHANGES FROM v4:
#   - Data file: TARA_ALL_sorted_v4.qs2 → TARA_ALL_annotated_final.qs2
#   - Annotation column: Manual_Annotation_refined → Annotation
#   - CD8 clusters: 6 (no Transitional CD8, no γδ in CD8 set)
#   - Output dirs: "Figure 2" / "Supplementary 2" (matching Fig 1 style)
#   - Added manuscript stats extraction section
#
# LAYOUT:
#   Fig 2A(i):  UMAP — clonal expansion split by HEI / HEU / HUU
#   Fig 2A(ii): Stacked bar — cloneSize per CD8 cluster
#   Fig 2B:     Alluvial plots + VL/CD4 clinical data per patient
#   Fig 2C:     Waffle charts — Trex epitope specificity
#
# SUPPLEMENTARY:
#   S2A: Clone size distribution by condition
#   S2B: Clonal persistence per patient
#   S2C: Cross-cluster clonotype sharing
#   S2D: Temporal epitope specificity
#
# REQUIRES: TARA_ALL_annotated_final.qs2, TARA_VL_CD4.xlsx
################################################################################


# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(scCustomize)
library(scRepertoire)
library(SeuratExtend)
library(qs2)
library(Trex)
library(readxl)
library(RColorBrewer)
library(ggrepel)


# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
fig2_dir   <- file.path(base_dir, "Manuscript", "Figure 2")
supp2_dir  <- file.path(base_dir, "Manuscript", "Supplementary 2")
vl_path    <- file.path(base_dir, "TARA_VL_CD4.xlsx")

for (d in c(fig2_dir, supp2_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ── Load Data ─────────────────────────────────────────────────────────────────
# CHANGED: now using final annotated object
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_annotated_final.qs2"))


# ── CD8 cluster definitions (match Fig 1 v5) ─────────────────────────────────
# CHANGED: 6 clusters, no number prefixes, no Transitional, no γδ
cd8_cluster_names <- c(
  "Naive 1 CD8", "Naive 2 CD8", "Naive Intermediate CD8",
  "Tscm CD8", "TEMRA/CTL CD8", "Tex CD8"
)

cd8_short_labels <- setNames(
  c("Naive 1", "Naive 2", "Naive Intermediate",
    "Tscm", "TEMRA/CTL", "Tex"),
  cd8_cluster_names
)

cluster_cols <- c(
  "Naive 1 CD8"            = "#2166AC",
  "Naive 2 CD8"            = "#67A9CF",
  "Naive Intermediate CD8" = "#D6604D",
  "Tscm CD8"               = "#1A9850",
  "TEMRA/CTL CD8"          = "#B2182B",
  "Tex CD8"                = "#762A83"
)

cond_cols <- c("HUU" = "#4E79A7", "HEU" = "#F28E2B", "HEI" = "#E15759")


# ── Clone size palette ────────────────────────────────────────────────────────
colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)
clone_order <- c("Single (0 < X <= 1)",
                 "Small (1 < X <= 5)",
                 "Medium (5 < X <= 20)",
                 "Large (20 < X <= 100)",
                 "Hyperexpanded (100 < X <= 500)")
clone_colors <- setNames(colorblind_vector[c(1, 3, 4, 5, 7)], clone_order)


################################################################################
# FIG 2A(i) — UMAP: Clonal expansion split by condition
################################################################################

message("Generating Fig 2A(i)...")

conditions <- c("HEI", "HEU", "HUU")
umap_list <- list()

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  sub_obj <- subset(TARA_ALL, subset = Condition == cond)
  
  umap_coords <- as.data.frame(Embeddings(sub_obj, reduction = "wnn.umap"))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$cloneSize <- sub_obj$cloneSize
  
  umap_coords <- umap_coords %>%
    mutate(
      cloneSize_plot = ifelse(is.na(cloneSize), "No TCR", as.character(cloneSize)),
      cloneSize_plot = factor(cloneSize_plot, levels = c(clone_order, "No TCR"))
    )
  
  clone_colors_ext <- c(clone_colors, "No TCR" = "#E0E0E0")
  
  x_arrow <- min(umap_coords$UMAP1, na.rm = TRUE) + 1
  y_arrow <- min(umap_coords$UMAP2, na.rm = TRUE) + 1
  
  umap_gray  <- umap_coords %>% filter(cloneSize_plot == "No TCR")
  umap_color <- umap_coords %>% filter(cloneSize_plot != "No TCR")
  
  p <- ggplot() +
    geom_point(data = umap_gray, aes(x = UMAP1, y = UMAP2),
               color = "#E0E0E0", size = 0.6, alpha = 0.5) +
    geom_point(data = umap_color, aes(x = UMAP1, y = UMAP2, color = cloneSize_plot),
               size = 1.2, alpha = 0.8) +
    scale_color_manual(values = clone_colors_ext, drop = FALSE) +
    ggtitle(cond) +
    theme_void() +
    theme(
      plot.title = element_text(size = 48, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_fixed()
  
  if (i == 1) {
    p <- p +
      annotate("segment",
               x = x_arrow, xend = x_arrow + 4,
               y = y_arrow, yend = y_arrow,
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               linewidth = 1.5, color = "black") +
      annotate("text", x = x_arrow + 2, y = y_arrow - 1.8,
               label = "UMAP1", size = 10, fontface = "bold") +
      annotate("segment",
               x = x_arrow, xend = x_arrow,
               y = y_arrow, yend = y_arrow + 4,
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               linewidth = 1.5, color = "black") +
      annotate("text", x = x_arrow - 1.8, y = y_arrow + 2,
               label = "UMAP2", size = 10, fontface = "bold", angle = 90)
  }
  
  umap_list[[cond]] <- p
}

# Shared legend
legend_plot <- ggplot(data.frame(x = 1:5, y = 1, cs = factor(clone_order, levels = clone_order)),
                      aes(x = x, y = y, fill = cs)) +
  geom_tile() +
  scale_fill_manual(values = clone_colors, name = "Clone Size") +
  theme_void(base_size = 32) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 32),
        legend.text = element_text(size = 26),
        legend.key.size = unit(1.5, "cm")) +
  guides(fill = guide_legend(nrow = 1))

legend_grob <- cowplot::get_legend(legend_plot)

fig2a_umap <- (umap_list[["HEI"]] | umap_list[["HEU"]] | umap_list[["HUU"]]) /
  wrap_elements(legend_grob) +
  plot_layout(heights = c(1, 0.1))

ggsave(file.path(fig2_dir, "Fig2A_i_ClonalExpansion_UMAP.png"),
       fig2a_umap, width = 30, height = 14, dpi = 300, bg = "white")

message("\u2713 Fig 2A(i) saved\n")


################################################################################
# FIG 2A(ii) — Stacked bar: cloneSize per CD8 cluster
################################################################################

message("Generating Fig 2A(ii)...")

# CHANGED: uses Annotation column, 6 clusters
TARA_CD8 <- subset(TARA_ALL, subset = Annotation %in% cd8_cluster_names)
TARA_CD8$Annotation <- droplevels(
  factor(TARA_CD8$Annotation, levels = cd8_cluster_names)
)

meta_cd8 <- TARA_CD8@meta.data %>%
  mutate(
    Cluster = cd8_short_labels[as.character(Annotation)],
    Cluster = factor(Cluster, levels = cd8_short_labels),
    cloneSize = factor(cloneSize, levels = clone_order)
  ) %>%
  filter(!is.na(cloneSize))

bar_data <- meta_cd8 %>%
  group_by(Cluster, cloneSize, .drop = FALSE) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

fig2a_bar <- ggplot(bar_data, aes(x = Cluster, y = n, fill = cloneSize)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.75, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = clone_colors, name = "Clone Size", drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of Cells") +
  theme_classic(base_size = 32) +
  theme(
    axis.text.x  = element_text(size = 26, angle = 30, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 26),
    axis.title.y = element_text(size = 30, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(fig2_dir, "Fig2A_ii_ClonalOccupy_CD8clusters.png"),
       fig2a_bar, width = 14, height = 10, dpi = 300, bg = "white")

# Separate bar legend
bar_legend_df <- data.frame(
  cs = factor(clone_order, levels = clone_order),
  x = seq_along(clone_order),
  y = 1
)

p_bar_legend <- ggplot(bar_legend_df, aes(x = x, y = y)) +
  geom_point(aes(fill = cs), shape = 22, size = 12, stroke = 0.5) +
  geom_text(aes(label = cs), vjust = -1.5, size = 7, fontface = "bold") +
  scale_fill_manual(values = clone_colors) +
  scale_x_continuous(limits = c(0, length(clone_order) + 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 2.5), clip = "off") +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(5, 20, 20, 20))

ggsave(file.path(fig2_dir, "Fig2A_ii_Legend.png"),
       p_bar_legend, width = 24, height = 2.5, dpi = 300, bg = "white")

message("\u2713 Fig 2A(ii) saved\n")


################################################################################
# FIG 2B — Alluvial Plots + VL/CD4 clinical data (SEPARATED)
################################################################################

message("Generating Fig 2B...")

# ── Load clinical data ────────────────────────────────────────────────────────
vl_cd4 <- read_excel(vl_path) %>%
  rename(
    PID = `STUDY ID`,
    Age_months = `AGE IN MONTHS`,
    VL = `HIV cp/mL`,
    VL_log = `VL Log`,
    CD4_mm3 = mm3,
    CD4_pct = `CD4%`
  ) %>%
  filter(PID %in% c("CP003", "CP006", "CP013", "CP018", "CP020")) %>%
  mutate(
    Age_months = as.numeric(Age_months),
    VL = as.numeric(VL),
    CD4_mm3 = as.numeric(CD4_mm3)
  ) %>%
  filter(!is.na(Age_months))

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
  if (file.exists(t_file)) assign(paste0(name, ".TCR"), read.csv(t_file))
}

TARA.TCR             <- paste0(TARA_names, ".TCR")
TARA.contig_list.TCR <- as.list(mget(TARA.TCR))
combined.TCR.TARA    <- combineTCR(TARA.contig_list.TCR, samples = TARA_names)

# ── Patient info (TCR folder names → new orig.ident mapping) ─────────────────
patient_info <- list(
  CP003 = list(samples = c("CP003_entry", "CP003_V12", "CP003_V24"), max_age = 24),
  CP006 = list(samples = c("CP006_entry", "CP006_12m", "CP006_V24"), max_age = 24),
  CP013 = list(samples = c("CP013_1m", "CP013_12m", "CP013_24m"),    max_age = 25),
  CP018 = list(samples = c("CP018_entry", "CP018_V24", "CP018_42m"), max_age = 42),
  CP020 = list(samples = c("CP020_V1", "CP020_V12", "CP020_V44"),    max_age = 44)
)

# CHANGED: maps TCR folder names → new age-based orig.ident
tcr_to_orig <- c(
  "CP003_entry" = "CP003_1m",  "CP003_V12" = "CP003_12m", "CP003_V24" = "CP003_24m",
  "CP006_entry" = "CP006_1m",  "CP006_12m" = "CP006_12m", "CP006_V24" = "CP006_24m",
  "CP013_1m"    = "CP013_1m",  "CP013_12m" = "CP013_12m", "CP013_24m" = "CP013_25m",
  "CP018_entry" = "CP018_1m",  "CP018_V24" = "CP018_25m", "CP018_42m" = "CP018_42m",
  "CP020_V1"    = "CP020_1m",  "CP020_V12" = "CP020_12m", "CP020_V44" = "CP020_44m"
)

# CHANGED: uses new orig.ident (age-based) for age lookup
age_lookup <- TARA_ALL@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Age = median(Age, na.rm = TRUE), .groups = "drop") %>%
  deframe()

for (pt in names(patient_info)) {
  samps <- patient_info[[pt]]$samples
  orig_names <- tcr_to_orig[samps]
  ages <- age_lookup[orig_names]
  patient_info[[pt]]$labels <- paste0(round(ages), "m")
  cat(pt, ":", paste(samps, "\u2192", orig_names, "\u2192", patient_info[[pt]]$labels), "\n")
}

# ── Alluvial palette ─────────────────────────────────────────────────────────
alluvial_pal <- c(
  "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#264653",
  "#F4A261", "#606C38", "#A8DADC", "#D62828", "#023E8A",
  "#6D6875", "#BC6C25", "#118AB2", "#EF476F", "#06D6A0",
  "#8338EC", "#3A86FF", "#FB5607", "#FF006E", "#FFBE0B",
  hcl.colors(30, palette = "Dark 3")
)

# ── Generate plots per patient ────────────────────────────────────────────────
alluvial_plots <- list()
vl_plots <- list()
cd4_plots <- list()
clinical_combined_plots <- list()
clinical_combined_grobs <- list()

for (pt in names(patient_info)) {
  info <- patient_info[[pt]]
  
  # ── ALLUVIAL PLOT ──────────────────────────────────────────────────────────
  p_alluv <- clonalCompare(
    combined.TCR.TARA,
    top.clones     = 20,
    samples        = info$samples,
    order.by       = info$samples,
    cloneCall      = "strict",
    relabel.clones = TRUE,
    proportion     = FALSE,
    graph          = "alluvial"
  )
  
  build <- ggplot_build(p_alluv)
  raw_dat <- p_alluv$data
  fill_col_name <- tryCatch(
    rlang::as_name(p_alluv$mapping$fill),
    error = function(e) {
      char_cols <- names(raw_dat)[sapply(raw_dat, is.character) | sapply(raw_dat, is.factor)]
      char_cols[which.max(sapply(char_cols, function(x) length(unique(raw_dat[[x]]))))]
    }
  )
  fill_lvls <- unique(as.character(raw_dat[[fill_col_name]]))
  n_fills <- length(fill_lvls)
  pal_named <- setNames(rep_len(alluvial_pal, n_fills), fill_lvls)
  
  p_alluv <- p_alluv +
    scale_x_discrete(labels = info$labels) +
    scale_fill_manual(values = pal_named) +
    labs(x = NULL, y = "Count", title = pt) +
    theme_classic(base_size = 48) +
    theme(
      plot.title   = element_text(size = 72, face = "bold", hjust = 0.5),
      axis.text.x  = element_text(size = 64, face = "bold"),
      axis.text.y  = element_text(size = 56),
      axis.title.y = element_text(size = 64, face = "bold"),
      legend.position = "none"
    )
  
  alluvial_plots[[pt]] <- p_alluv
  ggsave(file.path(fig2_dir, paste0("Fig2B_Alluvial_", pt, ".png")),
         p_alluv, width = 14, height = 12, dpi = 300, bg = "white")
  
  # ── CLINICAL DATA ──────────────────────────────────────────────────────────
  pt_clin <- vl_cd4 %>%
    filter(PID == pt, Age_months > 0, Age_months <= info$max_age)
  pt_vl  <- pt_clin %>% filter(!is.na(VL))
  pt_cd4 <- pt_clin %>% filter(!is.na(CD4_mm3))
  
  message(sprintf("  %s: VL points = %d, CD4 points = %d", pt, nrow(pt_vl), nrow(pt_cd4)))
  
  # ── VL PLOT ────────────────────────────────────────────────────────────────
  p_vl <- ggplot(pt_vl, aes(x = Age_months, y = log10(pmax(VL, 1)))) +
    geom_hline(yintercept = log10(200), linetype = "dashed",
               color = "#666666", linewidth = 2.5) +
    geom_line(color = "#B2182B", linewidth = 4) +
    geom_point(color = "#B2182B", size = 10, shape = 16) +
    scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), expand = c(0.02, 0)) +
    scale_x_continuous(
      limits = c(0, info$max_age),
      breaks = seq(0, info$max_age, by = ifelse(info$max_age > 30, 12, 6)),
      expand = c(0.02, 0)
    ) +
    labs(x = "Age (months)", y = expression("Viral Load (log"[10]*" cp/mL)"),
         title = paste0(pt, " - Viral Load")) +
    theme_bw(base_size = 40) +
    theme(
      plot.title       = element_text(size = 48, face = "bold", hjust = 0.5),
      axis.text        = element_text(size = 36, color = "black", face = "bold"),
      axis.title       = element_text(size = 40, face = "bold"),
      axis.title.y     = element_text(color = "#B2182B"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", linewidth = 2, fill = NA),
      plot.margin      = margin(15, 15, 15, 15)
    )
  
  vl_plots[[pt]] <- p_vl
  ggsave(file.path(fig2_dir, paste0("Fig2B_VL_", pt, ".png")),
         p_vl, width = 10, height = 8, dpi = 300, bg = "white")
  
  # ── CD4 PLOT ───────────────────────────────────────────────────────────────
  p_cd4 <- ggplot(pt_cd4, aes(x = Age_months, y = CD4_mm3)) +
    geom_line(color = "#2166AC", linewidth = 4) +
    geom_point(color = "#2166AC", size = 10, shape = 16) +
    scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 2000), expand = c(0.02, 0)) +
    scale_x_continuous(
      limits = c(0, info$max_age),
      breaks = seq(0, info$max_age, by = ifelse(info$max_age > 30, 12, 6)),
      expand = c(0.02, 0)
    ) +
    labs(x = "Age (months)", y = expression("CD4 Count (cells/"*mu*"L)"),
         title = paste0(pt, " - CD4 Count")) +
    theme_bw(base_size = 40) +
    theme(
      plot.title       = element_text(size = 48, face = "bold", hjust = 0.5),
      axis.text        = element_text(size = 36, color = "black", face = "bold"),
      axis.title       = element_text(size = 40, face = "bold"),
      axis.title.y     = element_text(color = "#2166AC"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", linewidth = 2, fill = NA),
      plot.margin      = margin(15, 15, 15, 15)
    )
  
  cd4_plots[[pt]] <- p_cd4
  ggsave(file.path(fig2_dir, paste0("Fig2B_CD4_", pt, ".png")),
         p_cd4, width = 10, height = 8, dpi = 300, bg = "white")
  
  # ── COMBINED VL + CD4 (dual y-axis) ────────────────────────────────────────
  cd4_max <- 6000
  vl_max_log <- 8
  scale_factor_val <- cd4_max / vl_max_log
  
  pt_vl_plot  <- pt_vl %>% mutate(VL_log = log10(pmax(VL, 1)))
  pt_cd4_plot <- pt_cd4 %>% mutate(CD4_scaled = CD4_mm3 / scale_factor_val)
  
  p_clinical <- ggplot() +
    geom_hline(yintercept = log10(200), linetype = "dashed",
               color = "#666666", linewidth = 2.5) +
    geom_line(data = pt_vl_plot, aes(x = Age_months, y = VL_log),
              color = "#B2182B", linewidth = 4) +
    geom_point(data = pt_vl_plot, aes(x = Age_months, y = VL_log),
               color = "#B2182B", size = 10, shape = 16) +
    geom_line(data = pt_cd4_plot, aes(x = Age_months, y = CD4_scaled),
              color = "#2166AC", linewidth = 4) +
    geom_point(data = pt_cd4_plot, aes(x = Age_months, y = CD4_scaled),
               color = "#2166AC", size = 10, shape = 16) +
    scale_y_continuous(
      name = expression("VL (log"[10]*" cp/mL)"),
      limits = c(0, vl_max_log),
      breaks = seq(0, 8, by = 2),
      sec.axis = sec_axis(~ . * scale_factor_val,
                          name = expression("CD4 (cells/"*mu*"L)"),
                          breaks = seq(0, 6000, by = 2000))
    ) +
    scale_x_continuous(
      limits = c(0, info$max_age),
      breaks = seq(0, info$max_age, by = ifelse(info$max_age > 30, 12, 6)),
      expand = c(0.02, 0)
    ) +
    labs(x = "Age (months)", title = pt) +
    theme_bw(base_size = 52) +
    theme(
      plot.title         = element_text(size = 80, face = "bold", hjust = 0.5),
      axis.text          = element_text(size = 48, color = "black", face = "bold"),
      axis.text.x        = element_text(size = 48, face = "bold"),
      axis.title.x       = element_text(size = 56, face = "bold"),
      axis.title.y.left  = element_text(size = 56, face = "bold", color = "#B2182B"),
      axis.title.y.right = element_text(size = 56, face = "bold", color = "#2166AC"),
      axis.text.y.left   = element_text(size = 48, color = "#B2182B", face = "bold"),
      axis.text.y.right  = element_text(size = 48, color = "#2166AC", face = "bold"),
      panel.grid.major   = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", linewidth = 2, fill = NA),
      plot.margin        = margin(20, 25, 20, 25)
    )
  
  clinical_combined_plots[[pt]] <- p_clinical
  clinical_combined_grobs[[pt]] <- ggplotGrob(p_clinical)
  
  ggsave(file.path(fig2_dir, paste0("Fig2B_Clinical_", pt, ".png")),
         p_clinical, width = 12, height = 8, dpi = 300, bg = "white")
}

# ── Save combined grids ──────────────────────────────────────────────────────
alluvial_grid <- wrap_plots(alluvial_plots, ncol = 5)
ggsave(file.path(fig2_dir, "Fig2B_Alluvial_AllPatients.png"),
       alluvial_grid, width = 60, height = 12, dpi = 200, bg = "white",
       limitsize = FALSE)

vl_grid <- wrap_plots(vl_plots, ncol = 5)
ggsave(file.path(fig2_dir, "Fig2B_VL_AllPatients.png"),
       vl_grid, width = 50, height = 8, dpi = 300, bg = "white",
       limitsize = FALSE)

cd4_grid <- wrap_plots(cd4_plots, ncol = 5)
ggsave(file.path(fig2_dir, "Fig2B_CD4_AllPatients.png"),
       cd4_grid, width = 50, height = 8, dpi = 300, bg = "white",
       limitsize = FALSE)

clinical_grid <- cowplot::plot_grid(
  clinical_combined_grobs[["CP003"]],
  clinical_combined_grobs[["CP006"]],
  clinical_combined_grobs[["CP013"]],
  clinical_combined_grobs[["CP018"]],
  clinical_combined_grobs[["CP020"]],
  ncol = 5, align = "h", axis = "tb"
)
ggsave(file.path(fig2_dir, "Fig2B_Clinical_AllPatients.png"),
       clinical_grid, width = 80, height = 14, dpi = 300, bg = "white",
       limitsize = FALSE)

message("\u2713 Fig 2B saved (all individual + combined grids)\n")


################################################################################
# FIG 2C — Waffle Charts: Epitope specificity (HEI only, ED2)
################################################################################

message("Generating Fig 2C...")

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

epitope_colors <- c(
  "HIV" = "#E41A1C", "CMV" = "#377EB8", "EBV" = "#4DAF4A",
  "Influenza" = "#FF7F00", "TB" = "#984EA3",
  "Other" = "#A65628", "Unknown" = "#D9D9D9"
)
epitope_order <- c("HIV", "CMV", "EBV", "Influenza", "TB", "Other", "Unknown")

waffle_counts <- waffle_data %>%
  group_by(PID, Epitope_Clean) %>%
  summarise(n_clones = n(), .groups = "drop") %>%
  mutate(Epitope_Clean = factor(Epitope_Clean, levels = epitope_order))

# ── Uniform grid ─────────────────────────────────────────────────────────────
hei_patients <- c("CP003", "CP006", "CP013", "CP018", "CP020")
ncol_waffle  <- 10

max_clones <- waffle_counts %>%
  filter(PID %in% hei_patients) %>%
  group_by(PID) %>%
  summarise(total = sum(n_clones)) %>%
  pull(total) %>%
  max()

max_squares  <- ceiling(max_clones / ncol_waffle) * ncol_waffle
nrow_waffle  <- max_squares / ncol_waffle

# ── Helper ────────────────────────────────────────────────────────────────────
make_waffle_df <- function(counts, ncol, total_squares) {
  categories <- rep(names(counts), times = counts)
  n_pad <- total_squares - length(categories)
  if (n_pad > 0) categories <- c(categories, rep("_empty", n_pad))
  n <- length(categories)
  data.frame(
    category = factor(categories, levels = c(epitope_order, "_empty")),
    x = ((seq_len(n) - 1) %% ncol) + 1,
    y = ((seq_len(n) - 1) %/% ncol) + 1
  )
}

epitope_colors_ext <- c(epitope_colors, "_empty" = NA)

# ── Per-patient waffles ──────────────────────────────────────────────────────
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
      plot.title    = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5, color = "grey40",
                                   margin = margin(t = 4)),
      legend.position = "none",
      plot.margin     = margin(8, 8, 8, 8)
    )
  
  waffle_plots[[pid]] <- p
  ggsave(file.path(fig2_dir, paste0("Fig2C_Waffle_", pid, ".png")),
         p, width = 7, height = 7, dpi = 300, bg = "white")
}

# Separate legend
legend_df <- data.frame(
  Epitope = factor(epitope_order, levels = epitope_order),
  y = seq_along(epitope_order)
)

p_waffle_legend <- ggplot(legend_df, aes(x = 1, y = y)) +
  geom_point(aes(fill = Epitope), shape = 22, size = 10, stroke = 0) +
  geom_text(aes(x = 0.82, label = as.character(Epitope)), size = 6,
            fontface = "bold", lineheight = 0.8, vjust = 0.5) +
  scale_fill_manual(values = epitope_colors) +
  scale_x_continuous(limits = c(0.6, 1.1)) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(0, 2, 0, 2))

ggsave(file.path(fig2_dir, "Fig2C_Waffle_Legend.png"),
       p_waffle_legend, width = 14, height = 1.4, dpi = 300, bg = "white")

waffle_grid <- wrap_plots(waffle_plots, ncol = 5)
ggsave(file.path(fig2_dir, "Fig2C_Waffle_Combined.png"),
       waffle_grid, width = 40, height = 9, dpi = 300, bg = "white")

message("\u2713 Fig 2C saved (individual + combined + legend)\n")


################################################################################
# SUPPLEMENTARY FIGURE S2
################################################################################

message("Generating Supplementary Figure S2 panels...")


# ══════════════════════════════════════════════════════════════════════════════
# S2A: Clone size distribution by condition
# ══════════════════════════════════════════════════════════════════════════════

message("  S2A: Clone size distribution by condition...")

# CHANGED: uses Annotation column
clonesize_df <- TARA_ALL@meta.data %>%
  filter(Annotation %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency)) %>%
  mutate(Condition = factor(Condition, levels = c("HEI", "HEU", "HUU")))

p_s2a_box <- ggplot(clonesize_df, aes(x = Condition, y = clonalFrequency, fill = Condition)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = cond_cols) +
  labs(x = NULL, y = "Clonal Frequency",
       title = "Clonal frequency by condition") +
  theme_classic(base_size = 18) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

expand_summary <- clonesize_df %>%
  group_by(Condition) %>%
  summarise(
    n_total    = n(),
    n_expanded = sum(clonalFrequency > 1),
    pct_expanded = n_expanded / n_total * 100,
    .groups = "drop"
  )

p_s2a_bar <- ggplot(expand_summary, aes(x = Condition, y = pct_expanded, fill = Condition)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct_expanded)), vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = cond_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(expand_summary$pct_expanded) * 1.2)) +
  labs(x = NULL, y = "% Expanded Cells",
       title = "% clonally expanded") +
  theme_classic(base_size = 18) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

p_s2a <- p_s2a_box | p_s2a_bar
ggsave(file.path(supp2_dir, "S2A_CloneSizeDistribution.png"),
       p_s2a, width = 10, height = 5, dpi = 300, bg = "white")

message("  \u2713 S2A saved")


# ══════════════════════════════════════════════════════════════════════════════
# S2B: Clonal persistence per patient
# ══════════════════════════════════════════════════════════════════════════════

message("  S2B: Clonal persistence summary...")

hei_5 <- c("CP003", "CP006", "CP013", "CP018", "CP020")

persistence_df <- TARA_ALL@meta.data %>%
  filter(!is.na(CTstrict) & !is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(PID = sub("_.*$", "", orig.ident)) %>%
  filter(PID %in% hei_5) %>%
  group_by(PID, CTstrict) %>%
  summarise(
    n_timepoints = n_distinct(orig.ident),
    n_cells      = n(),
    .groups      = "drop"
  )

persistence_summary <- persistence_df %>%
  group_by(PID) %>%
  summarise(
    total_expanded_clones  = n(),
    clones_1_timepoint     = sum(n_timepoints == 1),
    clones_2_timepoints    = sum(n_timepoints == 2),
    clones_3_timepoints    = sum(n_timepoints == 3),
    pct_persistent         = sum(n_timepoints >= 2) / n() * 100,
    .groups = "drop"
  )

cat("\n=== Clonal persistence summary (expanded clones, HEI) ===\n")
print(as.data.frame(persistence_summary))

persist_long <- persistence_df %>%
  mutate(Persistence = ifelse(n_timepoints >= 2, "Persistent (2+ timepoints)", "Single timepoint")) %>%
  group_by(PID, Persistence) %>%
  summarise(n_clones = n(), .groups = "drop") %>%
  mutate(PID = factor(PID, levels = hei_5))

p_s2b <- ggplot(persist_long, aes(x = PID, y = n_clones, fill = Persistence)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Persistent (2+ timepoints)" = "#E63946",
                               "Single timepoint" = "#A8DADC"),
                    name = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "N Expanded Clonotypes",
       title = "Clonal persistence") +
  theme_classic(base_size = 18) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.text    = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top"
  )

ggsave(file.path(supp2_dir, "S2B_ClonalPersistence.png"),
       p_s2b, width = 6, height = 5, dpi = 300, bg = "white")

message("  \u2713 S2B saved")


# ══════════════════════════════════════════════════════════════════════════════
# S2C: Cross-cluster clonotype sharing
# ══════════════════════════════════════════════════════════════════════════════

message("  S2C: Cross-cluster clonotype sharing...")

# CHANGED: uses Annotation column, 6 clusters
expanded_meta <- TARA_ALL@meta.data %>%
  filter(Annotation %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(
    Cluster = cd8_short_labels[as.character(Annotation)],
    Cluster = factor(Cluster, levels = cd8_short_labels)
  )

clone_cluster_span <- expanded_meta %>%
  group_by(CTstrict) %>%
  summarise(
    n_clusters = n_distinct(Cluster),
    clusters   = paste(sort(unique(as.character(Cluster))), collapse = " + "),
    n_cells    = n(),
    .groups    = "drop"
  )

cat("\n=== Expanded clone cluster spanning (HEI) ===\n")
cat("Clones in 1 cluster: ", sum(clone_cluster_span$n_clusters == 1), "\n")
cat("Clones in 2 clusters:", sum(clone_cluster_span$n_clusters == 2), "\n")
cat("Clones in 3+ clusters:", sum(clone_cluster_span$n_clusters >= 3), "\n")

if (sum(clone_cluster_span$n_clusters >= 2) > 0) {
  cat("\nMulti-cluster clone combinations:\n")
  print(table(clone_cluster_span$clusters[clone_cluster_span$n_clusters >= 2]))
}

p_s2c_span <- ggplot(clone_cluster_span, aes(x = factor(n_clusters))) +
  geom_bar(fill = "#4E79A7", width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 6, fontface = "bold") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(x = "Clusters per Clone", y = "N Clonotypes",
       title = "Cluster span") +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold")
  )

if (sum(clone_cluster_span$n_clusters >= 2) > 0) {
  multi_clones <- clone_cluster_span %>% filter(n_clusters >= 2)
  
  cluster_pairs <- expanded_meta %>%
    filter(CTstrict %in% multi_clones$CTstrict) %>%
    group_by(CTstrict) %>%
    summarise(clusters = list(sort(unique(as.character(Cluster)))), .groups = "drop")
  
  all_clusters <- levels(expanded_meta$Cluster)
  cooccur_mat <- matrix(0, length(all_clusters), length(all_clusters),
                        dimnames = list(all_clusters, all_clusters))
  
  for (cl_list in cluster_pairs$clusters) {
    if (length(cl_list) >= 2) {
      pairs <- combn(cl_list, 2)
      for (j in seq_len(ncol(pairs))) {
        cooccur_mat[pairs[1, j], pairs[2, j]] <- cooccur_mat[pairs[1, j], pairs[2, j]] + 1
        cooccur_mat[pairs[2, j], pairs[1, j]] <- cooccur_mat[pairs[2, j], pairs[1, j]] + 1
      }
    }
  }
  
  cooccur_df <- as.data.frame(as.table(cooccur_mat)) %>%
    rename(Cluster1 = Var1, Cluster2 = Var2, n_shared = Freq) %>%
    filter(as.numeric(Cluster1) <= as.numeric(Cluster2))
  
  p_s2c_cooccur <- ggplot(cooccur_df, aes(x = Cluster1, y = Cluster2, fill = n_shared)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(n_shared > 0, n_shared, "")), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#E15759", name = "Shared\nClones") +
    labs(title = "Cluster co-occurrence", x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title  = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12)
    )
  
  p_s2c <- p_s2c_span | p_s2c_cooccur
  ggsave(file.path(supp2_dir, "S2C_CrossCluster_ClonalSharing.png"),
         p_s2c, width = 12, height = 5, dpi = 300, bg = "white")
} else {
  ggsave(file.path(supp2_dir, "S2C_CrossCluster_ClonalSharing.png"),
         p_s2c_span, width = 6, height = 5, dpi = 300, bg = "white")
}

message("  \u2713 S2C saved")


# ══════════════════════════════════════════════════════════════════════════════
# S2D: Temporal epitope specificity
# ══════════════════════════════════════════════════════════════════════════════

message("  S2D: Temporal epitope specificity...")

trex_temporal <- TARA_ALL_TRB_2@meta.data %>%
  filter(!is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(
    PID = sub("_.*$", "", orig.ident),
    Epitope_Species = ifelse(is.na(TRB_Epitope.species), "Unknown", TRB_Epitope.species)
  ) %>%
  distinct(PID, orig.ident, CTstrict, Epitope_Species) %>%
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

# CHANGED: uses new age-based orig.ident for lookup
trex_temporal$Age <- age_lookup[trex_temporal$orig.ident]

temporal_counts <- trex_temporal %>%
  filter(PID %in% hei_5) %>%
  group_by(PID, orig.ident, Age, Epitope_Clean) %>%
  summarise(n_clones = n(), .groups = "drop") %>%
  group_by(PID, orig.ident, Age) %>%
  mutate(pct = n_clones / sum(n_clones) * 100,
         total = sum(n_clones)) %>%
  ungroup() %>%
  mutate(
    Epitope_Clean = factor(Epitope_Clean, levels = epitope_order),
    Age_label = paste0(round(Age), "m"),
    PID = factor(PID, levels = hei_5)
  )

p_s2d <- ggplot(temporal_counts,
                aes(x = reorder(Age_label, Age), y = pct, fill = Epitope_Clean)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(data = temporal_counts %>% distinct(PID, orig.ident, Age, Age_label, total),
            aes(x = reorder(Age_label, Age), y = 108, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 5, fontface = "italic") +
  facet_wrap(~PID, ncol = 5, scales = "free_x") +
  scale_fill_manual(values = epitope_colors, name = "Epitope") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) +
  labs(x = "Age (months)", y = "% of Expanded Clones") +
  theme_classic(base_size = 18) +
  theme(
    strip.text      = element_text(size = 20, face = "bold"),
    axis.text.x     = element_text(size = 16, face = "bold"),
    axis.text.y     = element_text(size = 16),
    axis.title      = element_text(size = 18, face = "bold"),
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "right"
  )

ggsave(file.path(supp2_dir, "S2D_Temporal_EpitopeSpecificity.png"),
       p_s2d, width = 14, height = 4.5, dpi = 300, bg = "white")

message("  \u2713 S2D saved")


################################################################################
# MANUSCRIPT STATS EXTRACTION — Figure 2
#
# Prints all values needed to update the manuscript text
################################################################################

cat("\n")
cat("================================================================\n")
cat("  FIGURE 2 — MANUSCRIPT STATISTICS\n")
cat("================================================================\n\n")

# ── All CD8 cells with TCR, by condition ──────────────────────────────────────
meta_all <- TARA_ALL@meta.data %>%
  filter(Annotation %in% cd8_cluster_names)

cat("── CD8 CLONAL EXPANSION BY CONDITION ──\n\n")

for (cond in c("HEI", "HEU", "HUU")) {
  cond_cells <- meta_all %>% filter(Condition == cond)
  tcr_cells  <- cond_cells %>% filter(!is.na(clonalFrequency))
  expanded   <- tcr_cells %>% filter(clonalFrequency > 1)
  n_clonotypes <- length(unique(expanded$CTstrict))
  
  cat(sprintf("%s:\n", cond))
  cat(sprintf("  Total CD8 cells: %d\n", nrow(cond_cells)))
  cat(sprintf("  With productive TCR: %d\n", nrow(tcr_cells)))
  cat(sprintf("  Expanded cells: %d (%.1f%% of CD8)\n",
              nrow(expanded), 100 * nrow(expanded) / nrow(cond_cells)))
  cat(sprintf("  Unique expanded clonotypes: %d\n\n", n_clonotypes))
}

cat("── PER-CLUSTER EXPANSION (HEI, all timepoints) ──\n\n")

cluster_expansion_all <- meta_all %>%
  filter(Condition == "HEI") %>%
  mutate(
    is_expanded = !is.na(clonalFrequency) & clonalFrequency > 1,
    Cluster = cd8_short_labels[as.character(Annotation)]
  ) %>%
  group_by(Cluster) %>%
  summarise(
    total    = n(),
    expanded = sum(is_expanded),
    pct      = round(100 * mean(is_expanded), 1),
    .groups  = "drop"
  ) %>%
  arrange(desc(pct))

print(as.data.frame(cluster_expansion_all))

cat("\n── CLONAL PERSISTENCE (per patient) ──\n\n")
print(as.data.frame(persistence_summary))

cat("\n── CROSS-CLUSTER SHARING ──\n\n")
cat("Clones in 1 cluster: ", sum(clone_cluster_span$n_clusters == 1), "\n")
cat("Clones in 2 clusters:", sum(clone_cluster_span$n_clusters == 2), "\n")
cat("Clones in 3+ clusters:", sum(clone_cluster_span$n_clusters >= 3), "\n")
cat("Total expanded clonotypes:", nrow(clone_cluster_span), "\n")

if (sum(clone_cluster_span$n_clusters >= 2) > 0) {
  cat("\nTop multi-cluster combinations:\n")
  combo_table <- sort(table(clone_cluster_span$clusters[clone_cluster_span$n_clusters >= 2]),
                      decreasing = TRUE)
  print(head(combo_table, 10))
}

cat("\n── EPITOPE SPECIFICITY (per patient, waffle chart values) ──\n\n")

epitope_by_patient <- waffle_data %>%
  filter(PID %in% hei_5) %>%
  group_by(PID, Epitope_Clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Epitope_Clean, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(-PID)))

print(as.data.frame(epitope_by_patient))
################################################################################
# RESERVOIR EXTENSION — TABLE S1 + CSVs
#
# Run AFTER Figure_2.R in the same R session. Requires:
#   TARA_ALL, persistence_df, waffle_data, age_lookup, hei_5
#
# OUTPUTS (to Manuscript/Supplementary 2/):
#   Supplementary_Table_S1.csv   — per-PID reservoir + Jaccard + persistence
#   Jaccard_Consecutive_Timepoints.csv
#   ReservoirExtension_SpearmanCor.csv
#
# No supplementary figure renumbering needed — this is just a table.
################################################################################

library(tidyverse)

stopifnot(exists("TARA_ALL"), exists("persistence_df"), exists("hei_5"),
          exists("waffle_data"), exists("age_lookup"))

supp2_dir <- file.path("~/Documents/CD8_Longitudinal", "Manuscript", "Supplementary 2")
dir.create(supp2_dir, recursive = TRUE, showWarnings = FALSE)


################################################################################
# 1. RESERVOIR DATA
################################################################################

reservoir_summary <- tribble(
  ~PID,    ~ART_status,            ~Reservoir_timepoint_m, ~Intact_per_mil, ~Pct_intact, ~LTR_per_mil, ~Total_HIV_1m,
  "CP003", "Intermittent / rebound",          12,             282.73,          41.19,      690.75,       3988.42,
  "CP006", "Rapidly suppressed",              12,               4.39,          14.27,       65.51,      31491.73,
  "CP013", "Suppressed (low seeding)",        87,               0.00,           0.00,        0.00,         41.47,
  "CP018", "Suppressed",                      80,                 NA,             NA,       18.32,       3470.07,
  "CP020", "Never suppressed",                84,               0.00,           0.00,      312.30,       5578.01
)


################################################################################
# 2. PERSISTENCE METRICS
################################################################################

persistence_metrics <- persistence_df %>%
  group_by(PID) %>%
  summarise(
    n_expanded_clones = n(),
    n_persist_2p      = sum(n_timepoints >= 2),
    pct_persist_2p    = round(100 * n_persist_2p / n_expanded_clones, 1),
    max_timepoints    = max(n_timepoints),
    .groups = "drop"
  )


################################################################################
# 3. JACCARD OVERLAP
################################################################################

expanded_by_sample <- TARA_ALL@meta.data %>%
  filter(Condition == "HEI",
         !is.na(CTstrict), !is.na(clonalFrequency), clonalFrequency > 1) %>%
  mutate(PID = sub("_.*$", "", orig.ident)) %>%
  filter(PID %in% hei_5) %>%
  mutate(Age = age_lookup[orig.ident]) %>%
  distinct(PID, orig.ident, Age, CTstrict)

jaccard_transitions <- expanded_by_sample %>%
  group_by(PID) %>%
  arrange(Age, .by_group = TRUE) %>%
  summarise(samples = list(unique(orig.ident)),
            ages    = list(sort(unique(Age))),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(res = list({
    samps <- samples
    if (length(samps) < 2) {
      tibble(Timepoint_1 = character(), Timepoint_2 = character(),
             Clones_tp1 = integer(), Clones_tp2 = integer(),
             Shared = integer(), Jaccard = numeric())
    } else {
      cur_pid <- PID
      df_order <- expanded_by_sample %>% filter(PID == cur_pid) %>%
        distinct(orig.ident, Age) %>% arrange(Age)
      samps_ord <- df_order$orig.ident
      ages_ord  <- df_order$Age
      map_dfr(seq_len(length(samps_ord) - 1), function(i) {
        a <- samps_ord[i]; b <- samps_ord[i + 1]
        A <- expanded_by_sample %>% filter(orig.ident == a) %>% pull(CTstrict) %>% unique()
        B <- expanded_by_sample %>% filter(orig.ident == b) %>% pull(CTstrict) %>% unique()
        inter <- length(intersect(A, B))
        uni   <- length(union(A, B))
        tibble(Timepoint_1 = paste0(round(ages_ord[i]), "m"),
               Timepoint_2 = paste0(round(ages_ord[i + 1]), "m"),
               Clones_tp1 = length(A), Clones_tp2 = length(B),
               Shared = inter,
               Jaccard = if (uni == 0) NA_real_ else round(inter / uni, 4))
      })
    }
  })) %>%
  ungroup() %>%
  select(PID, res) %>%
  unnest(res)


################################################################################
# 4. HIV CLONE COUNTS
################################################################################

hiv_clones_per_pid <- waffle_data %>%
  filter(PID %in% hei_5) %>%
  group_by(PID) %>%
  summarise(
    n_HIV   = sum(Epitope_Clean == "HIV"),
    pct_HIV = round(100 * n_HIV / n(), 1),
    .groups = "drop"
  )


################################################################################
# 5. SPEARMAN CORRELATIONS (reservoir vs TCR)
################################################################################

combined_for_cor <- reservoir_summary %>%
  left_join(persistence_metrics, by = "PID") %>%
  left_join(hiv_clones_per_pid, by = "PID")

resv_cols <- c("Intact_per_mil", "Pct_intact", "LTR_per_mil")
tcr_cols  <- c("n_persist_2p", "pct_persist_2p", "n_HIV")
numeric_cols <- combined_for_cor %>% select(where(is.numeric)) %>% colnames()

cor_rows <- list()
for (rm in intersect(resv_cols, numeric_cols)) {
  for (tm in intersect(tcr_cols, numeric_cols)) {
    x <- combined_for_cor[[rm]]; y <- combined_for_cor[[tm]]
    ok <- complete.cases(x, y)
    if (sum(ok) >= 3) {
      test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      cor_rows[[length(cor_rows) + 1]] <- tibble(
        reservoir_metric = rm, tcr_metric = tm,
        n = sum(ok), rho = round(unname(test$estimate), 3), p = round(test$p.value, 4))
    }
  }
}
cor_table <- bind_rows(cor_rows) %>% arrange(p)


################################################################################
# 6. TABLE S1: COMBINED
################################################################################

# Pivot Jaccard into wide columns
jaccard_wide <- jaccard_transitions %>%
  mutate(transition = paste0(Timepoint_1, "\u2192", Timepoint_2)) %>%
  select(PID, transition, Jaccard) %>%
  pivot_wider(names_from = transition, values_from = Jaccard,
              names_prefix = "Jaccard_")

table_s1 <- reservoir_summary %>%
  left_join(persistence_metrics, by = "PID") %>%
  left_join(jaccard_wide, by = "PID") %>%
  left_join(hiv_clones_per_pid, by = "PID") %>%
  select(PID, ART_status, Reservoir_timepoint_m,
         Intact_per_mil, Pct_intact, LTR_per_mil, Total_HIV_1m,
         starts_with("Jaccard_"),
         n_expanded_clones, n_persist_2p, pct_persist_2p, max_timepoints,
         n_HIV, pct_HIV)


################################################################################
# 7. SAVE
################################################################################

write.csv(table_s1, file.path(supp2_dir, "Supplementary_Table_S1.csv"), row.names = FALSE)
write.csv(jaccard_transitions, file.path(supp2_dir, "Jaccard_Consecutive_Timepoints.csv"), row.names = FALSE)
write.csv(cor_table, file.path(supp2_dir, "ReservoirExtension_SpearmanCor.csv"), row.names = FALSE)

cat("\n=== SUPPLEMENTARY TABLE S1 ===\n")
print(as.data.frame(table_s1))

cat("\n=== JACCARD TRANSITIONS ===\n")
print(as.data.frame(jaccard_transitions))

cat("\n=== SPEARMAN CORRELATIONS (p < 0.1) ===\n")
print(as.data.frame(cor_table %>% filter(p < 0.1)))

cat("\n=== MANUSCRIPT TEXT (copy-paste) ===\n\n")
cat("Jaccard overlap between the 1-month and 12-month expanded clonotype\n")
cat("repertoires was non-zero only for CP003 (J = 0.106; J = 0.000 for all\n")
cat("other participants), while all participants showed moderate overlap at\n")
cat("their second transition (J = 0.13-0.28; Supplementary Table S1).\n")

cat("\n  DONE — 3 CSVs saved to:", supp2_dir, "\n")

