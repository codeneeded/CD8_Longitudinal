################################################################################
# FIGURE 2: TCR Clonal Expansion Analysis — TARA Cohort (v4 annotations)
#
# Panel A: (i)  UMAP — clonal expansion split by HEI / HEU / HUU
#          (ii) Stacked bar — cloneSize per CD8 cluster
# Panel B: Alluvial plots + VL/CD4 clinical data per patient
# Panel C: Waffle charts — Trex epitope specificity + separate legend
#
# REQUIRES: TARA_ALL_sorted_v4.qs2, TARA_VL_CD4.xlsx
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

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
fig2_dir   <- file.path(base_dir, "Manuscript", "Fig 2")
vl_path    <- file.path(base_dir, "TARA_VL_CD4.xlsx")
dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load Data ─────────────────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted_v4.qs2"))

# ── CD8 cluster definitions (match Fig 1 v4) ─────────────────────────────────
cd8_cluster_names <- c(
  "1a: Naive 1 CD8", "6: Naive 2 CD8", "1c: Naive Intermediate CD8",
  "1b: Tscm CD8", "12a: Transitional CD8", "8: TEMRA/CTL CD8",
  "27: Tex CD8", "9: γδ T cell"
)

cd8_short_labels <- setNames(
  c("Naive 1", "Naive 2", "Naive Intermediate",
    "Tscm", "Transitional", "TEMRA/CTL", "Tex", "γδ T cell"),
  cd8_cluster_names
)

cluster_cols <- c(
  "1a: Naive 1 CD8"            = "#4E79A7",
  "6: Naive 2 CD8"             = "#76B7B2",
  "1c: Naive Intermediate CD8" = "#E15759",
  "1b: Tscm CD8"               = "#59A14F",
  "12a: Transitional CD8"      = "#FF9DA7",
  "8: TEMRA/CTL CD8"           = "#EDC948",
  "27: Tex CD8"                = "#9C755F",
  "9: γδ T cell"               = "#B07AA1"
)

# ── Use full object (includes longitudinal data) ─────────────────────────────
# TARA_ALL already has a Condition column from the annotation script.
# For the UMAP we show HEI, HEU, HUU (HEI = all HEI timepoints combined)

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
#              Order: HEI, HEU, HUU. UMAP arrows on last panel only.
################################################################################

message("Generating Fig 2A(i)...")

conditions <- c("HEI", "HEU", "HUU")
umap_list <- list()

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  sub_obj <- subset(TARA_ALL, subset = Condition == cond)
  
  # UMAP arrows only on HEI (first panel)
  if (i == 1) {
    theme_list <- list(NoAxes(), theme_umap_arrows(),
                       ggtitle(cond),
                       theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
                             legend.position = "none"))
  } else {
    theme_list <- list(NoAxes(),
                       ggtitle(cond),
                       theme(plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
                             legend.position = "none"))
  }
  
  p <- DimPlot2(
    sub_obj,
    reduction  = "wnn.umap",
    group.by   = "cloneSize",
    cols       = clone_colors,
    pt.size    = 0.4,
    raster     = FALSE,
    theme      = theme_list
  )
  umap_list[[cond]] <- p
}

# Shared legend
legend_plot <- ggplot(data.frame(x = 1:5, y = 1, cs = factor(clone_order, levels = clone_order)),
                      aes(x = x, y = y, fill = cs)) +
  geom_tile() +
  scale_fill_manual(values = clone_colors, name = "Clone Size") +
  theme_void(base_size = 26) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 26),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1, "cm")) +
  guides(fill = guide_legend(nrow = 1))

legend_grob <- cowplot::get_legend(legend_plot)

fig2a_umap <- (umap_list[["HEI"]] | umap_list[["HEU"]] | umap_list[["HUU"]]) /
  wrap_elements(legend_grob) +
  plot_layout(heights = c(1, 0.08))

ggsave(file.path(fig2_dir, "Fig2A_i_ClonalExpansion_UMAP.png"),
       fig2a_umap, width = 28, height = 12, dpi = 300, bg = "white")

message("✓ Fig 2A(i) saved\n")


################################################################################
# FIG 2A(ii) — Stacked bar: cloneSize per CD8 cluster (bigger fonts, no labels)
################################################################################

message("Generating Fig 2A(ii)...")

TARA_CD8 <- subset(TARA_ALL, subset = Manual_Annotation_refined %in% cd8_cluster_names)
TARA_CD8$Manual_Annotation_refined <- droplevels(
  factor(TARA_CD8$Manual_Annotation_refined, levels = cd8_cluster_names)
)

meta_cd8 <- TARA_CD8@meta.data %>%
  mutate(
    Cluster = cd8_short_labels[as.character(Manual_Annotation_refined)],
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
    axis.text.x       = element_text(size = 26, angle = 30, hjust = 1, face = "bold"),
    axis.text.y       = element_text(size = 26),
    axis.title.y      = element_text(size = 30, face = "bold"),
    legend.position   = "none"
  )

ggsave(file.path(fig2_dir, "Fig2A_ii_ClonalOccupy_CD8clusters.png"),
       fig2a_bar, width = 16, height = 10, dpi = 300, bg = "white")

# ── Separate bar legend ──────────────────────────────────────────────────────
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

message("✓ Fig 2A(ii) saved\n")


################################################################################
# FIG 2B — Alluvial Plots + VL/CD4 clinical data
#
# Plot only up to the latest timepoint in the alluvial.
# Below each alluvial: VL and CD4 trajectories from TARA_VL_CD4.xlsx
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

# ── Patient info: alluvial samples, max age for clinical ─────────────────────
# Timepoint labels will be pulled from TARA_ALL$Age metadata
patient_info <- list(
  CP003 = list(samples = c("CP003_entry", "CP003_V12", "CP003_V24"), max_age = 24),
  CP006 = list(samples = c("CP006_entry", "CP006_12m", "CP006_V24"), max_age = 24),
  CP013 = list(samples = c("CP013_1m", "CP013_12m", "CP013_24m"),    max_age = 25),
  CP018 = list(samples = c("CP018_entry", "CP018_V24", "CP018_42m"), max_age = 42),
  CP020 = list(samples = c("CP020_V1", "CP020_V12", "CP020_V44"),    max_age = 44)
)

# ── Build labels from TARA_ALL$Age metadata ──────────────────────────────────
# orig.ident uses names like CP003_1m, CP003_12m, CP003_24m
# TCR/alluvial uses names like CP003_entry, CP003_V12, CP003_V24
# Map TCR names → orig.ident names to look up ages

tcr_to_orig <- c(
  "CP003_entry" = "CP003_1m",  "CP003_V12" = "CP003_12m", "CP003_V24" = "CP003_24m",
  "CP006_entry" = "CP006_1m",  "CP006_12m" = "CP006_12m", "CP006_V24" = "CP006_24m",
  "CP013_1m"    = "CP013_1m",  "CP013_12m" = "CP013_12m", "CP013_24m" = "CP013_25m",
  "CP018_entry" = "CP018_1m",  "CP018_V24" = "CP018_25m", "CP018_42m" = "CP018_42m",
  "CP020_V1"    = "CP020_1m",  "CP020_V12" = "CP020_12m", "CP020_V44" = "CP020_44m"
)

age_lookup <- TARA_ALL@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Age = median(Age, na.rm = TRUE), .groups = "drop") %>%
  deframe()

for (pt in names(patient_info)) {
  samps <- patient_info[[pt]]$samples
  orig_names <- tcr_to_orig[samps]
  ages <- age_lookup[orig_names]
  patient_info[[pt]]$labels <- paste0(round(ages), "m")
  cat(pt, ":", paste(samps, "→", orig_names, "→", patient_info[[pt]]$labels), "\n")
}

# ── Custom alluvial palette — 50 distinct colors ─────────────────────────────
alluvial_pal <- c(
  "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#264653",
  "#F4A261", "#606C38", "#A8DADC", "#D62828", "#023E8A",
  "#6D6875", "#BC6C25", "#118AB2", "#EF476F", "#06D6A0",
  "#8338EC", "#3A86FF", "#FB5607", "#FF006E", "#FFBE0B",
  hcl.colors(30, palette = "Dark 3")
)

# ── Generate per-patient panels: alluvial + combined VL/CD4 ──────────────────
patient_panels <- list()

for (pt in names(patient_info)) {
  info <- patient_info[[pt]]
  
  # ── Alluvial plot with custom colors ────────────────────────────────────
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
  
  # Extract fill levels from the ggplot object and build named palette
  build <- ggplot_build(p_alluv)
  fill_scale <- build$plot$scales$get_scales("fill")
  # Get the levels from the data — try the fill column in the raw data
  raw_dat <- p_alluv$data
  # Find which column is mapped to fill
  fill_col_name <- tryCatch(
    rlang::as_name(p_alluv$mapping$fill),
    error = function(e) {
      # Fallback: find column with most unique character values
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
  
  # ── Combined VL + CD4 on same plot (dual y-axis) ───────────────────────
  pt_clin <- vl_cd4 %>%
    filter(PID == pt, Age_months <= info$max_age)
  
  pt_vl  <- pt_clin %>% filter(!is.na(VL) & VL > 0)
  pt_cd4 <- pt_clin %>% filter(!is.na(CD4_mm3) & CD4_mm3 > 0)
  
  cd4_max  <- max(pt_cd4$CD4_mm3, 5000, na.rm = TRUE)
  vl_max_log <- 8
  scale_factor <- cd4_max / vl_max_log
  
  p_clinical <- ggplot() +
    geom_hline(yintercept = log10(200), linetype = "dashed",
               color = "grey50", linewidth = 1.2) +
    geom_line(data = pt_vl, aes(x = Age_months, y = log10(pmax(VL, 1))),
              color = "#E15759", linewidth = 2.5) +
    geom_point(data = pt_vl, aes(x = Age_months, y = log10(pmax(VL, 1))),
               color = "#E15759", size = 7) +
    geom_line(data = pt_cd4, aes(x = Age_months, y = CD4_mm3 / scale_factor),
              color = "#4E79A7", linewidth = 2.5) +
    geom_point(data = pt_cd4, aes(x = Age_months, y = CD4_mm3 / scale_factor),
               color = "#4E79A7", size = 7) +
    scale_y_continuous(
      name = expression("VL (log"[10]*" cp/mL)"),
      limits = c(0, vl_max_log),
      sec.axis = sec_axis(~ . * scale_factor,
                          name = expression("CD4 (cells/"*mu*"L)"))
    ) +
    scale_x_continuous(limits = c(0, info$max_age),
                       breaks = seq(0, info$max_age, by = ifelse(info$max_age > 30, 12, 6))) +
    labs(x = "Age (months)") +
    theme_classic(base_size = 48) +
    theme(
      axis.text          = element_text(size = 48),
      axis.title.x       = element_text(size = 52, face = "bold"),
      axis.title.y.left  = element_text(size = 52, face = "bold", color = "#E15759"),
      axis.title.y.right = element_text(size = 52, face = "bold", color = "#4E79A7"),
      axis.text.y.left   = element_text(size = 48, color = "#E15759"),
      axis.text.y.right  = element_text(size = 48, color = "#4E79A7"),
      axis.text.x        = element_text(size = 48)
    )
  
  # Stack: alluvial on top, clinical below
  p_combined <- p_alluv / p_clinical +
    plot_layout(heights = c(2.5, 1))
  
  patient_panels[[pt]] <- p_combined
  
  ggsave(file.path(fig2_dir, paste0("Fig2B_", pt, ".png")),
         p_combined, width = 20, height = 16, dpi = 300, bg = "white")
}

# ── Combined grid (all 5 patients side by side) ──────────────────────────────
fig2b <- wrap_plots(patient_panels, ncol = 5)

ggsave(file.path(fig2_dir, "Fig2B_Alluvial_AllPatients.png"),
       fig2b, width = 50, height = 16, dpi = 200, bg = "white",
       limitsize = FALSE)

message("✓ Fig 2B saved (individual + combined)\n")


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

# ── Uniform grid ──────────────────────────────────────────────────────────────
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

# ── Per-patient waffles ───────────────────────────────────────────────────────
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

# ── Separate legend (same approach as Fig 1C legend) ──────────────────────────
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

# ── Combined waffles (no legend — legend is separate file) ───────────────────
waffle_grid <- wrap_plots(waffle_plots, ncol = 5)

ggsave(file.path(fig2_dir, "Fig2C_Waffle_Combined.png"),
       waffle_grid, width = 40, height = 9, dpi = 300, bg = "white")

message("✓ Fig 2C saved (individual + combined + legend)\n")


################################################################################
# SUPPLEMENTARY FIGURE S2 — Supporting analyses for Figure 2
#
# Panel order matches manuscript narrative:
#   S2A: Clone size distribution by condition (HEI vs HEU vs HUU)
#   S2B: Clonal persistence per patient
#   S2C: Cross-cluster clonotype sharing (cluster span + co-occurrence)
#   S2D: Temporal epitope specificity (stacked bars per timepoint)
#
# Output: ~/Documents/CD8_Longitudinal/Manuscript/Supplementary 2/
################################################################################

s2_dir <- file.path(base_dir, "Manuscript", "Supplementary 2")
dir.create(s2_dir, recursive = TRUE, showWarnings = FALSE)

message("Generating Supplementary Figure S2 panels...")


# ══════════════════════════════════════════════════════════════════════════════
# S2A: Clone size distribution by condition (HEI vs HEU vs HUU)
# ══════════════════════════════════════════════════════════════════════════════

message("  S2A: Clone size distribution by condition...")

clonesize_df <- TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency)) %>%
  mutate(Condition = factor(Condition, levels = c("HEI", "HEU", "HUU")))

p_s2a_box <- ggplot(clonesize_df, aes(x = Condition, y = clonalFrequency, fill = Condition)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = c("HEI" = "#E15759", "HEU" = "#4E79A7", "HUU" = "#59A14F")) +
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
  scale_fill_manual(values = c("HEI" = "#E15759", "HEU" = "#4E79A7", "HUU" = "#59A14F")) +
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

ggsave(file.path(s2_dir, "S2A_CloneSizeDistribution.png"),
       p_s2a, width = 10, height = 5, dpi = 300, bg = "white")

message("  ✓ S2A saved")


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
  mutate(Persistence = ifelse(n_timepoints >= 2, "Persistent (2+ visits)", "Single visit")) %>%
  group_by(PID, Persistence) %>%
  summarise(n_clones = n(), .groups = "drop") %>%
  mutate(PID = factor(PID, levels = hei_5))

p_s2b <- ggplot(persist_long, aes(x = PID, y = n_clones, fill = Persistence)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Persistent (2+ visits)" = "#E63946",
                               "Single visit" = "#A8DADC"),
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

ggsave(file.path(s2_dir, "S2B_ClonalPersistence.png"),
       p_s2b, width = 6, height = 5, dpi = 300, bg = "white")

message("  ✓ S2B saved")


# ══════════════════════════════════════════════════════════════════════════════
# S2C: Cross-cluster clonotype sharing
# ══════════════════════════════════════════════════════════════════════════════

message("  S2C: Cross-cluster clonotype sharing...")

expanded_meta <- TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(
    Cluster = cd8_short_labels[as.character(Manual_Annotation_refined)],
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
    plot.title   = element_text(size = 20, face = "bold"),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 18, face = "bold")
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
      plot.title   = element_text(size = 20, face = "bold"),
      axis.text.x  = element_text(size = 13, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12)
    )
  
  p_s2c <- p_s2c_span | p_s2c_cooccur
  ggsave(file.path(s2_dir, "S2C_CrossCluster_ClonalSharing.png"),
         p_s2c, width = 12, height = 5, dpi = 300, bg = "white")
} else {
  p_s2c <- p_s2c_span
  ggsave(file.path(s2_dir, "S2C_CrossCluster_ClonalSharing.png"),
         p_s2c, width = 6, height = 5, dpi = 300, bg = "white")
}

message("  ✓ S2C saved")


# ══════════════════════════════════════════════════════════════════════════════
# S2D: Temporal epitope specificity — stacked bar per timepoint per patient
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
    strip.text       = element_text(size = 20, face = "bold"),
    axis.text.x      = element_text(size = 16, face = "bold"),
    axis.text.y      = element_text(size = 16),
    axis.title       = element_text(size = 18, face = "bold"),
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(0.5, "cm"),
    legend.position  = "right"
  )

ggsave(file.path(s2_dir, "S2D_Temporal_EpitopeSpecificity.png"),
       p_s2d, width = 14, height = 4.5, dpi = 300, bg = "white")

message("  ✓ S2D saved")


################################################################################
# SUMMARY
################################################################################

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " Supplementary S2 panels saved to: ", s2_dir, "\n",
        "══════════════════════════════════════════════════════════════\n",
        "\n",
        " SUPPLEMENTARY S2 (matches manuscript narrative order):\n",
        " S2A: Clone size distribution by condition (box + bar)\n",
        " S2B: Clonal persistence (persistent vs single-visit)\n",
        " S2C: Cross-cluster clonotype sharing (cluster co-occurrence)\n",
        " S2D: Temporal epitope specificity (stacked bars per timepoint)\n",
        "══════════════════════════════════════════════════════════════\n"
)

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " All Figure 2 panels saved to: ", fig2_dir, "\n",
        " Supplementary S2 panels saved to: ", s2_dir, "\n",
        "══════════════════════════════════════════════════════════════\n",
        "\n",
        " MAIN FIGURE 2:\n",
        " 2A(i):  Clonal expansion UMAP (SeuratExtend, 3 panels)\n",
        " 2A(ii): Stacked bar — cloneSize per CD8 cluster\n",
        " 2B:     Alluvial + VL/CD4 per patient (5 panels)\n",
        " 2C:     Waffle charts + separate legend\n",
        "\n",
        " SUPPLEMENTARY S2:\n",
        " S2A: Temporal epitope specificity (stacked bars per timepoint)\n",
        " S2B: Clone size distribution by condition (box + bar)\n",
        " S2C: Clonal persistence summary (persistent vs single-visit)\n",
        " S2D: Cross-cluster clonotype sharing (cluster co-occurrence)\n",
        "══════════════════════════════════════════════════════════════\n"
)

################ Diagonostic Outputs #############3
################################################################################
# CONSOLE OUTPUT BLOCK 
################################################################################

library(tidyverse)

# ── 1. Clonal expansion summary by condition ─────────────────────────────────
cat("\n\n========== 1. CLONAL EXPANSION BY CONDITION ==========\n")
TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency)) %>%
  group_by(Condition) %>%
  summarise(
    n_total = n(),
    n_expanded = sum(clonalFrequency > 1),
    pct_expanded = round(n_expanded / n_total * 100, 1),
    n_unique_clonotypes_expanded = n_distinct(CTstrict[clonalFrequency > 1])
  ) %>%
  as.data.frame() %>%
  print()

# ── 2. Expanded clonotypes per CD8 cluster ───────────────────────────────────
cat("\n\n========== 2. EXPANSION PER CD8 CLUSTER ==========\n")
TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency)) %>%
  mutate(Cluster = cd8_short_labels[as.character(Manual_Annotation_refined)]) %>%
  group_by(Cluster) %>%
  summarise(
    n_total = n(),
    n_expanded = sum(clonalFrequency > 1),
    pct_expanded = round(n_expanded / n_total * 100, 1)
  ) %>%
  arrange(desc(pct_expanded)) %>%
  as.data.frame() %>%
  print()

# ── 3. Clonal persistence summary ────────────────────────────────────────────
cat("\n\n========== 3. CLONAL PERSISTENCE SUMMARY ==========\n")
persistence_df <- TARA_ALL@meta.data %>%
  filter(!is.na(CTstrict) & !is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(PID = sub("_.*$", "", orig.ident)) %>%
  filter(PID %in% c("CP003", "CP006", "CP013", "CP018", "CP020")) %>%
  group_by(PID, CTstrict) %>%
  summarise(
    n_timepoints = n_distinct(orig.ident),
    n_cells = n(),
    .groups = "drop"
  )

persistence_df %>%
  group_by(PID) %>%
  summarise(
    total_expanded_clones = n(),
    clones_1_timepoint = sum(n_timepoints == 1),
    clones_2_timepoints = sum(n_timepoints == 2),
    clones_3_timepoints = sum(n_timepoints == 3),
    pct_persistent = round(sum(n_timepoints >= 2) / n() * 100, 1),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print()

# ── 4. Cross-cluster spanning ────────────────────────────────────────────────
cat("\n\n========== 4. CROSS-CLUSTER SPANNING (HEI) ==========\n")
expanded_meta <- TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(
    Cluster = cd8_short_labels[as.character(Manual_Annotation_refined)],
    Cluster = factor(Cluster, levels = cd8_short_labels)
  )

clone_cluster_span <- expanded_meta %>%
  group_by(CTstrict) %>%
  summarise(
    n_clusters = n_distinct(Cluster),
    clusters = paste(sort(unique(as.character(Cluster))), collapse = " + "),
    n_cells = n(),
    .groups = "drop"
  )

cat("Clones in 1 cluster: ", sum(clone_cluster_span$n_clusters == 1), "\n")
cat("Clones in 2 clusters:", sum(clone_cluster_span$n_clusters == 2), "\n")
cat("Clones in 3+ clusters:", sum(clone_cluster_span$n_clusters >= 3), "\n")

if (sum(clone_cluster_span$n_clusters >= 2) > 0) {
  cat("\nMulti-cluster clone combinations:\n")
  print(table(clone_cluster_span$clusters[clone_cluster_span$n_clusters >= 2]))
}

# ── 5. Trex epitope counts per patient (waffle numbers) ──────────────────────
cat("\n\n========== 5. TREX EPITOPE COUNTS PER PATIENT ==========\n")
waffle_counts %>%
  filter(PID %in% c("CP003", "CP006", "CP013", "CP018", "CP020")) %>%
  pivot_wider(names_from = Epitope_Clean, values_from = n_clones, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  as.data.frame() %>%
  print()

# ── 6. Expand summary by condition (from S2B) ────────────────────────────────
cat("\n\n========== 6. EXPANSION SUMMARY BY CONDITION (S2B) ==========\n")
TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency)) %>%
  mutate(Condition = factor(Condition, levels = c("HEI", "HEU", "HUU"))) %>%
  group_by(Condition) %>%
  summarise(
    n_total = n(),
    n_expanded = sum(clonalFrequency > 1),
    pct_expanded = round(n_expanded / n_total * 100, 1),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print()

# ── 7. Patient timepoint labels ──────────────────────────────────────────────
cat("\n\n========== 7. PATIENT TIMEPOINT LABELS ==========\n")
tcr_to_orig <- c(
  "CP003_entry" = "CP003_1m",  "CP003_V12" = "CP003_12m", "CP003_V24" = "CP003_24m",
  "CP006_entry" = "CP006_1m",  "CP006_12m" = "CP006_12m", "CP006_V24" = "CP006_24m",
  "CP013_1m"    = "CP013_1m",  "CP013_12m" = "CP013_12m", "CP013_24m" = "CP013_25m",
  "CP018_entry" = "CP018_1m",  "CP018_V24" = "CP018_25m", "CP018_42m" = "CP018_42m",
  "CP020_V1"    = "CP020_1m",  "CP020_V12" = "CP020_12m", "CP020_V44" = "CP020_44m"
)

age_lookup <- TARA_ALL@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Age = median(Age, na.rm = TRUE), .groups = "drop") %>%
  deframe()

patient_info <- list(
  CP003 = list(samples = c("CP003_entry", "CP003_V12", "CP003_V24")),
  CP006 = list(samples = c("CP006_entry", "CP006_12m", "CP006_V24")),
  CP013 = list(samples = c("CP013_1m", "CP013_12m", "CP013_24m")),
  CP018 = list(samples = c("CP018_entry", "CP018_V24", "CP018_42m")),
  CP020 = list(samples = c("CP020_V1", "CP020_V12", "CP020_V44"))
)

for (pt in names(patient_info)) {
  samps <- patient_info[[pt]]$samples
  orig_names <- tcr_to_orig[samps]
  ages <- age_lookup[orig_names]
  cat(pt, ":", paste0(round(ages), "m", collapse = " → "), "\n")
}

# ── 8. CD8 cluster definitions (for reference) ──────────────────────────────
cat("\n\n========== 8. CD8 CLUSTER DEFINITIONS (v4) ==========\n")
data.frame(
  Full_Name = cd8_cluster_names,
  Short_Label = cd8_short_labels
) %>% print()

# ── 9. Total cells in dataset ────────────────────────────────────────────────
cat("\n\n========== 9. DATASET OVERVIEW ==========\n")
cat("Total cells in TARA_ALL:", ncol(TARA_ALL), "\n")
cat("Cells per condition:\n")
print(table(TARA_ALL$Condition))
cat("\nCD8 cells per condition:\n")
TARA_ALL@meta.data %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  pull(Condition) %>%
  table() %>%
  print()
cat("\nSamples (orig.ident):\n")
print(sort(unique(TARA_ALL$orig.ident)))

