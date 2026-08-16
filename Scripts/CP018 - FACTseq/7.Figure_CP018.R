################################################################################
# CP018 FACT-Seq — Figure 3 equivalent
#
# Mirrors "Scripts/3. Manuscript/Figure_3.R" panel for panel, using the same
# annotation vocabulary and palettes so a CP018 figure can sit beside the
# published CP003 one (or replace it as a two-participant figure).
#
# PANELS
#   B1 — UMAP: annotated clusters
#   B2 — UMAP: HIV-specific TCR overlay
#   C  — Alluvial: clonotype persistence across stim + unstimulated timepoints
#   D  — Functional gene expression of HIV-validated cells in unstimulated data
#   E  — Persistence histogram (timepoints per clonotype)
#   F  — Unstimulated CD8 subset composition of HIV-validated clones
#
# CP018-SPECIFIC PANELS (no CP003 counterpart)
#   G  — Replicate concordance: clonotype detection A vs B
#   H  — Type I IFN memory module in HIV-specific cells (the suppression test)
#
# REQUIRES (in order):
#   1-5  the CP018 pipeline scripts
#   6    6.HIV_Specific_TCR_TARA_bulk.R  -> tara_cp018_crossref.qs2
#   config: CP018_CLUSTER_ANNOTATION must be set
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(Seurat)
library(SeuratExtend)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
library(qs2)

FIG_DIR  <- file.path(CP018_ANALYSIS, "Figures", "Figure_CP018")
SUPP_DIR <- file.path(CP018_ANALYSIS, "Figures", "Supplementary_CP018")
mkdirs(FIG_DIR, SUPP_DIR)

# -----------------------------------------------------------------------------
# Load + gate on the manual annotation step
# -----------------------------------------------------------------------------
require_file(QS_ANNOTATED, "CP018 annotated object (script 5)")
cp018 <- qs_read(QS_ANNOTATED)

if (is.null(CP018_CLUSTER_ANNOTATION)) {
  stop("CP018_CLUSTER_ANNOTATION is NULL in config.R.\n",
       "Map each MNN cluster to the Figure 3 vocabulary first — see the block\n",
       "at the bottom of config.R and the cluster profile CSV from script 5.",
       call. = FALSE)
}

clusters_present <- levels(cp018$mnn_clusters_rna)
unmapped <- setdiff(clusters_present, names(CP018_CLUSTER_ANNOTATION))
if (length(unmapped) > 0) {
  stop("Clusters with no annotation in CP018_CLUSTER_ANNOTATION: ",
       paste(unmapped, collapse = ", "), call. = FALSE)
}
bad_labels <- setdiff(unique(CP018_CLUSTER_ANNOTATION), ANNO_ORDER)
if (length(bad_labels) > 0) {
  stop("Annotation labels not in ANNO_ORDER: ", paste(bad_labels, collapse = ", "),
       "\nUse only: ", paste(ANNO_ORDER, collapse = ", "), call. = FALSE)
}

cp018$Fig_Annotation <- factor(
  CP018_CLUSTER_ANNOTATION[as.character(cp018$mnn_clusters_rna)],
  levels = ANNO_ORDER
)

cat("Annotation assigned:\n"); print(table(cp018$Fig_Annotation))

# Unstimulated cross-reference object (script 6). Panels D/F/H need it.
crossref_path <- file.path(SAVED_DIR, "tara_cp018_crossref.qs2")
have_crossref <- file.exists(crossref_path)
if (have_crossref) {
  tara_cp018 <- qs_read(crossref_path)
  if (!"clone_match" %in% colnames(tara_cp018@meta.data)) {
    have_crossref <- FALSE
    message("crossref object lacks clone_match; skipping panels D/F/H.")
  }
} else {
  message("tara_cp018_crossref.qs2 not found — run script 6. Panels D/F/H skipped.")
}

red_use <- "umap.mnn.rna"

# -----------------------------------------------------------------------------
# Shared UMAP axis-arrow helper (matches the CP003 figure style)
# -----------------------------------------------------------------------------
umap_df <- as.data.frame(Embeddings(cp018, reduction = red_use))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

x_range <- diff(range(umap_df$UMAP1, na.rm = TRUE))
y_range <- diff(range(umap_df$UMAP2, na.rm = TRUE))
x_arrow <- min(umap_df$UMAP1, na.rm = TRUE) - x_range * 0.05
y_arrow <- min(umap_df$UMAP2, na.rm = TRUE) - y_range * 0.08
x_len   <- x_range * 0.12
y_len   <- y_range * 0.12

add_umap_arrows <- function(p) {
  p +
    annotate("segment", x = x_arrow, xend = x_arrow + x_len,
             y = y_arrow, yend = y_arrow,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             linewidth = 0.8, color = "black") +
    annotate("text", x = x_arrow + x_len / 2, y = y_arrow - y_range * 0.04,
             label = "UMAP1", size = 5, fontface = "bold") +
    annotate("segment", x = x_arrow, xend = x_arrow,
             y = y_arrow, yend = y_arrow + y_len,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             linewidth = 0.8, color = "black") +
    annotate("text", x = x_arrow - x_range * 0.04, y = y_arrow + y_len / 2,
             label = "UMAP2", size = 5, fontface = "bold", angle = 90) +
    coord_cartesian(clip = "off")
}

# ── Panel B1: annotated UMAP ─────────────────────────────────────────────────
umap_df$Annotation <- cp018$Fig_Annotation

centers <- umap_df %>%
  group_by(Annotation) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

panelB1 <- add_umap_arrows(
  ggplot(umap_df, aes(UMAP1, UMAP2, color = Annotation)) +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_color_manual(values = ANNO_COLS, drop = FALSE) +
    geom_label_repel(data = centers,
                     aes(label = Annotation, color = Annotation),
                     fill = "white", size = 7, fontface = "bold",
                     label.size = 0.4, box.padding = 0.6,
                     segment.color = "grey50", max.overlaps = 20,
                     show.legend = FALSE) +
    ggtitle("CP018 Stimulation") +
    theme_void() +
    theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
          legend.position = "none", plot.margin = margin(10, 10, 20, 20))
)

ggsave(file.path(FIG_DIR, "FigCP018_B1_UMAP_Annotated.png"), panelB1,
       width = 10, height = 9, dpi = 400, bg = "white")

# ── Panel B2: HIV-specific overlay ───────────────────────────────────────────
umap_df$HIV <- factor(as.character(cp018$HIV_Specific_TCR),
                      levels = c("Other", "HIV-Specific TCR"))

panelB2 <- add_umap_arrows(
  ggplot(umap_df[order(umap_df$HIV), ], aes(UMAP1, UMAP2, color = HIV)) +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_color_manual(values = HIV_COLS, name = NULL) +
    ggtitle("HIV-Specific TCR") +
    theme_void() +
    theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 16, face = "bold"),
          legend.position = "bottom", plot.margin = margin(10, 10, 20, 20))
)

ggsave(file.path(FIG_DIR, "FigCP018_B2_UMAP_HIV_Specific.png"), panelB2,
       width = 10, height = 9, dpi = 400, bg = "white")

# ── Panel C: alluvial clonotype persistence ──────────────────────────────────
stim_all <- cp018@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict), CTstrict != "") %>%
  mutate(timepoint_label = paste0(Timepoint, "\n(HIV-stim)")) %>%
  count(CTstrict, timepoint_label, name = "n_cells")

if (have_crossref) {
  tara_all <- tara_cp018@meta.data %>%
    filter(clone_match == "HIV-Specific Match", !is.na(CTstrict), CTstrict != "") %>%
    mutate(timepoint_label = sub("^CP018_", "", orig.ident)) %>%
    count(CTstrict, timepoint_label, name = "n_cells")
} else {
  tara_all <- stim_all[0, ]
}

alluvial_df <- bind_rows(stim_all, tara_all)

# Order timepoints chronologically: stim 2m, then unstimulated, then stim 90m
tp_levels <- c("2m\n(HIV-stim)",
               intersect(c("1m", "25m", "42m"), unique(alluvial_df$timepoint_label)),
               "90m\n(HIV-stim)")
alluvial_df$timepoint_label <- factor(alluvial_df$timepoint_label,
                                     levels = intersect(tp_levels,
                                                        unique(alluvial_df$timepoint_label)))

multi_tp <- alluvial_df %>%
  group_by(CTstrict) %>% filter(n_distinct(timepoint_label) >= 2) %>% ungroup()

if (nrow(multi_tp) > 0) {
  clone_order <- multi_tp %>% count(CTstrict, wt = n_cells, name = "total") %>%
    arrange(desc(total))
  n_clones  <- nrow(clone_order)
  clone_pal <- colorRampPalette(brewer.pal(min(max(n_clones, 3), 12), "Set3"))(n_clones)
  names(clone_pal) <- clone_order$CTstrict

  panelC <- ggplot(multi_tp,
                   aes(x = timepoint_label, y = n_cells,
                       stratum = CTstrict, alluvium = CTstrict, fill = CTstrict)) +
    geom_stratum(width = 0.3, color = "grey30", linewidth = 0.3) +
    geom_flow(alpha = 0.5, width = 0.3) +
    scale_fill_manual(values = clone_pal, guide = "none") +
    scale_x_discrete(expand = c(0.08, 0.08), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_classic(base_size = 18) +
    theme(axis.text.x = element_text(size = 15, hjust = 0.5, lineheight = 0.9,
                                    face = "bold"),
          plot.title    = element_text(size = 24, face = "bold"),
          plot.subtitle = element_text(size = 16, color = "#E63946", face = "bold")) +
    labs(x = NULL, y = "Number of cells",
         title = "CP018 HIV-Specific Clonotype Persistence",
         subtitle = paste0(n_clones, " clonotypes detected at \u22652 timepoints"))

  ggsave(file.path(FIG_DIR, "FigCP018_C_Alluvial.png"), panelC,
         width = 13, height = 9, dpi = 400, bg = "white")
} else {
  message("Panel C skipped — no clonotype detected at >=2 timepoints.")
}

# ── Panel E: persistence histogram ───────────────────────────────────────────
persist <- alluvial_df %>%
  group_by(CTstrict) %>%
  summarise(n_timepoints = n_distinct(timepoint_label), .groups = "drop") %>%
  count(n_timepoints, name = "n_clonotypes")

if (nrow(persist) > 0) {
  panelE <- ggplot(persist, aes(factor(n_timepoints), n_clonotypes)) +
    geom_col(fill = "#E63946", color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = n_clonotypes), vjust = -0.4, size = 6, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Timepoints detected", y = "HIV-specific clonotypes",
         title = "Clonotype persistence") +
    theme_classic(base_size = 18) +
    theme(plot.title = element_text(size = 24, face = "bold"))

  ggsave(file.path(FIG_DIR, "FigCP018_E_Persistence.png"), panelE,
         width = 7, height = 6, dpi = 400, bg = "white")
}

# ── Panels D / F / H: require the unstimulated cross-reference ───────────────
if (have_crossref) {

  DefaultAssay(tara_cp018) <- "RNA"
  matched_bc <- colnames(tara_cp018)[tara_cp018$clone_match == "HIV-Specific Match"]

  # Panel D — functional gene expression in matched cells
  func_groups <- list(
    Cytotoxicity        = c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7"),
    Chemokines          = c("CCL3", "CCL4", "CCL5", "IFNG"),
    `Effector TF`       = c("TBX21", "ZEB2", "PRDM1", "ID2"),
    Exhaustion          = c("PDCD1", "LAG3", "TIGIT", "TOX"),
    `Stemness`          = c("TCF7", "SELL", "BACH2", "IL7R"),
    `Type I IFN memory` = c("IFIT1", "IFIT3", "ISG15", "MX1")
  )

  expr <- LayerData(tara_cp018, assay = "RNA", layer = "data")

  func_tab <- bind_rows(lapply(names(func_groups), function(g) {
    genes <- intersect(func_groups[[g]], rownames(expr))
    if (length(genes) == 0) return(NULL)
    data.frame(
      Group = g, Gene = genes,
      pct = round(100 * Matrix::rowMeans(expr[genes, matched_bc, drop = FALSE] > 0), 1)
    )
  }))

  if (nrow(func_tab) > 0) {
    func_tab$Group <- factor(func_tab$Group, levels = names(func_groups))
    func_tab <- func_tab %>% arrange(Group, desc(pct))
    func_tab$Gene <- factor(func_tab$Gene, levels = func_tab$Gene)

    panelD <- ggplot(func_tab, aes(Gene, pct, fill = Group)) +
      geom_col(color = "black", linewidth = 0.25) +
      facet_grid(~ Group, scales = "free_x", space = "free_x") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = "% of cells expressing",
           title = "CP018 HIV-specific clonotypes at rest",
           subtitle = paste0(length(matched_bc), " matched cells (unstimulated)")) +
      theme_bw(base_size = 15) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
            legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 20, face = "bold"))

    ggsave(file.path(FIG_DIR, "FigCP018_D_Functional_Profile.png"), panelD,
           width = 13, height = 6, dpi = 400, bg = "white")
    write.csv(func_tab, file.path(FIG_DIR, "FigCP018_D_source_data.csv"),
              row.names = FALSE)
  }

  # Panel F — which unstimulated CD8 subset the matched clones occupy
  anno_col <- intersect(c("Annotation", "CD8_Annotation", "Manual_Annotation_refined"),
                        colnames(tara_cp018@meta.data))[1]
  if (!is.na(anno_col)) {
    comp <- tara_cp018@meta.data[matched_bc, , drop = FALSE] %>%
      count(.data[[anno_col]], name = "n_cells") %>%
      rename(Subset = 1) %>%
      mutate(pct = round(100 * n_cells / sum(n_cells), 1)) %>%
      arrange(desc(n_cells))

    panelF <- ggplot(comp, aes(reorder(Subset, n_cells), n_cells)) +
      geom_col(fill = "#E63946", color = "black", linewidth = 0.25) +
      geom_text(aes(label = paste0(n_cells, " (", pct, "%)")),
                hjust = -0.1, size = 5) +
      coord_flip(clip = "off") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
      labs(x = NULL, y = "Cells",
           title = "Resting subset of HIV-specific clones") +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(size = 20, face = "bold"),
            plot.margin = margin(5, 40, 5, 5))

    ggsave(file.path(FIG_DIR, "FigCP018_F_Subset_Composition.png"), panelF,
           width = 9, height = 6, dpi = 400, bg = "white")
    write.csv(comp, file.path(FIG_DIR, "FigCP018_F_source_data.csv"), row.names = FALSE)
  }

  # Panel H — the suppression test: IFN memory in matched vs unmatched cells.
  # CP018 suppressed (VL 20 at 25m and 42m); the Figure 4/5 claim is that
  # suppression installs a type I IFN memory signature. If HIV-specific clones
  # carry it here, that is the CP003-unavailable comparison.
  ifn_genes <- intersect(c("IFIT1", "IFIT3", "ISG15", "MX1"), rownames(tara_cp018))
  stem_genes <- intersect(c("TCF7", "SELL", "BACH2", "LEF1", "IL7R"), rownames(tara_cp018))

  if (length(ifn_genes) >= 2) {
    tara_cp018 <- AddModuleScore(tara_cp018, features = list(ifn_genes),
                                 name = "IFN_memory")
    tara_cp018 <- AddModuleScore(tara_cp018, features = list(stem_genes),
                                 name = "Stemness")

    mod_df <- tara_cp018@meta.data %>%
      transmute(clone_match, Timepoint_Group = if ("Timepoint_Group" %in% names(.))
                  Timepoint_Group else NA_character_,
                IFN = IFN_memory1, Stem = Stemness1) %>%
      pivot_longer(c(IFN, Stem), names_to = "module", values_to = "score") %>%
      mutate(module = recode(module, IFN = "Type I IFN memory", Stem = "Stemness"))

    stats_tab <- mod_df %>%
      group_by(module) %>%
      summarise(
        mean_matched   = mean(score[clone_match == "HIV-Specific Match"], na.rm = TRUE),
        mean_other     = mean(score[clone_match == "Other"], na.rm = TRUE),
        p_wilcox = tryCatch(
          wilcox.test(score[clone_match == "HIV-Specific Match"],
                      score[clone_match == "Other"])$p.value,
          error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(delta = round(mean_matched - mean_other, 4),
             p_adj = p.adjust(p_wilcox, method = "BH"))

    write.csv(stats_tab, file.path(FIG_DIR, "FigCP018_H_module_stats.csv"),
              row.names = FALSE)
    cat("\n=== Panel H: module scores, matched vs other (unstimulated) ===\n")
    print(as.data.frame(stats_tab))

    panelH <- ggplot(mod_df, aes(clone_match, score, fill = clone_match)) +
      geom_violin(scale = "width", alpha = 0.8, colour = "grey25") +
      geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
      facet_wrap(~ module, scales = "free_y") +
      scale_fill_manual(values = c("Other" = "#BDBDBD",
                                   "HIV-Specific Match" = "#E63946")) +
      labs(x = NULL, y = "Module score",
           title = "CP018 (ART-suppressed): resting state of HIV-specific clones") +
      theme_bw(base_size = 15) +
      theme(legend.position = "none",
            axis.text.x = element_text(face = "bold"),
            plot.title  = element_text(size = 18, face = "bold"))

    ggsave(file.path(FIG_DIR, "FigCP018_H_IFN_Stemness.png"), panelH,
           width = 10, height = 6, dpi = 400, bg = "white")
  }
}

# ── Panel G: replicate concordance (CP018-only) ──────────────────────────────
# CP003 had no replicate at its late timepoint, so this panel has no counterpart.
rep_df <- cp018@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict), CTstrict != "") %>%
  count(CTstrict, Timepoint, Replicate, name = "n_cells") %>%
  pivot_wider(names_from = Replicate, values_from = n_cells, values_fill = 0)

if (all(c("A", "B") %in% names(rep_df)) && nrow(rep_df) > 0) {
  rep_df$detected_in <- with(rep_df, ifelse(A > 0 & B > 0, "Both replicates",
                                     ifelse(A > 0, "A only", "B only")))

  lab <- rep_df %>% group_by(Timepoint) %>%
    summarise(pct_both = round(100 * mean(detected_in == "Both replicates"), 1),
              .groups = "drop")

  panelG <- ggplot(rep_df, aes(A + 1, B + 1, colour = detected_in)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.6, alpha = 0.8) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("Both replicates" = "#1A9850",
                                   "A only" = "#F4A261", "B only" = "#9D4EDD"),
                        name = NULL) +
    facet_wrap(~ Timepoint) +
    labs(title = "HIV-specific clonotype detection across technical replicates",
         subtitle = paste(sprintf("%s: %.1f%% in both", lab$Timepoint, lab$pct_both),
                          collapse = "   |   "),
         x = "Cells in replicate A (+1, log)", y = "Cells in replicate B (+1, log)") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(size = 17, face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(FIG_DIR, "FigCP018_G_Replicate_Concordance.png"), panelG,
         width = 10, height = 6, dpi = 400, bg = "white")
  write.csv(rep_df, file.path(FIG_DIR, "FigCP018_G_source_data.csv"), row.names = FALSE)
  print(as.data.frame(lab))
}

################################################################################
#                            SUPPLEMENTARY PANELS
################################################################################

# S-A: cell cycle by annotation
if ("Phase" %in% colnames(cp018@meta.data)) {
  phase_df <- cp018@meta.data %>%
    count(Fig_Annotation, Phase, name = "n") %>%
    group_by(Fig_Annotation) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()

  ggsave(file.path(SUPP_DIR, "SuppCP018_A_CellCycle_by_Annotation.png"),
         ggplot(phase_df, aes(Fig_Annotation, pct, fill = Phase)) +
           geom_col(color = "black", linewidth = 0.25) +
           scale_fill_manual(values = PHASE_COLS) +
           labs(x = NULL, y = "% of cells", title = "Cell cycle phase by annotation") +
           theme_classic(base_size = 14) +
           theme(axis.text.x = element_text(angle = 40, hjust = 1),
                 plot.title = element_text(face = "bold")),
         width = 9, height = 6, dpi = 400, bg = "white")
}

# S-B: HIV-specific fraction per annotation
hiv_by_anno <- cp018@meta.data %>%
  group_by(Fig_Annotation) %>%
  summarise(n_cells = n(),
            pct_hiv = 100 * mean(HIV_Specific_TCR == "HIV-Specific TCR"),
            .groups = "drop")

ggsave(file.path(SUPP_DIR, "SuppCP018_B_HIV_pct_by_Annotation.png"),
       ggplot(hiv_by_anno, aes(Fig_Annotation, pct_hiv, fill = Fig_Annotation)) +
         geom_col(color = "black", linewidth = 0.25) +
         geom_text(aes(label = paste0(round(pct_hiv, 1), "%")), vjust = -0.4, size = 4.5) +
         scale_fill_manual(values = ANNO_COLS, guide = "none") +
         scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
         labs(x = NULL, y = "% HIV-specific", title = "HIV-specific fraction by annotation") +
         theme_classic(base_size = 14) +
         theme(axis.text.x = element_text(angle = 40, hjust = 1),
               plot.title = element_text(face = "bold")),
       width = 9, height = 6, dpi = 400, bg = "white")
write.csv(hiv_by_anno, file.path(SUPP_DIR, "SuppCP018_B_source_data.csv"),
          row.names = FALSE)

# S-C: timepoint composition per annotation
comp_tp <- cp018@meta.data %>%
  count(Fig_Annotation, Timepoint, name = "n") %>%
  group_by(Fig_Annotation) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()

ggsave(file.path(SUPP_DIR, "SuppCP018_C_Timepoint_Composition.png"),
       ggplot(comp_tp, aes(Fig_Annotation, pct, fill = Timepoint)) +
         geom_col(color = "black", linewidth = 0.25) +
         scale_fill_manual(values = TP_COLORS) +
         labs(x = NULL, y = "% of cells", title = "Timepoint composition by annotation") +
         theme_classic(base_size = 14) +
         theme(axis.text.x = element_text(angle = 40, hjust = 1),
               plot.title = element_text(face = "bold")),
       width = 9, height = 6, dpi = 400, bg = "white")

# S-D: dot plot of annotation markers
marker_panel <- intersect(
  c("TCF7","LEF1","SELL","CCR7","IL7R","GZMK","TOX","PDCD1","HAVCR2","LAG3",
    "GZMB","GNLY","PRF1","FGFBP2","CX3CR1","MKI67","TOP2A","IFIT1","ISG15"),
  rownames(cp018))

if (length(marker_panel) > 0) {
  Idents(cp018) <- "Fig_Annotation"
  ggsave(file.path(SUPP_DIR, "SuppCP018_D_Marker_DotPlot.png"),
         DotPlot(cp018, features = marker_panel) +
           scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#E63946") +
           labs(title = "Annotation markers", x = NULL, y = NULL) +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                 plot.title = element_text(face = "bold")),
         width = 13, height = 6, dpi = 400, bg = "white")
}

# S-E: module scores within HIV-specific cells, 2m vs 90m
mod_cols <- grep("_score1$", colnames(cp018@meta.data), value = TRUE)
if (length(mod_cols) > 0) {
  hiv_mod <- cp018@meta.data %>%
    filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
    select(Timepoint, all_of(mod_cols)) %>%
    pivot_longer(all_of(mod_cols), names_to = "module", values_to = "score") %>%
    mutate(module = sub("_score1$", "", module))

  ggsave(file.path(SUPP_DIR, "SuppCP018_E_Modules_2m_vs_90m.png"),
         ggplot(hiv_mod, aes(Timepoint, score, fill = Timepoint)) +
           geom_violin(scale = "width", alpha = 0.8, colour = "grey25") +
           geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
           facet_wrap(~ module, scales = "free_y") +
           scale_fill_manual(values = TP_COLORS, guide = "none") +
           labs(x = NULL, y = "Module score",
                title = "HIV-specific cells: 2m vs 90m") +
           theme_bw(base_size = 14) +
           theme(plot.title = element_text(face = "bold")),
         width = 10, height = 7, dpi = 400, bg = "white")

  hiv_mod_stats <- hiv_mod %>%
    group_by(module) %>%
    summarise(mean_2m  = mean(score[Timepoint == "2m"], na.rm = TRUE),
              mean_90m = mean(score[Timepoint == "90m"], na.rm = TRUE),
              p_wilcox = tryCatch(wilcox.test(score[Timepoint == "2m"],
                                              score[Timepoint == "90m"])$p.value,
                                  error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(p_adj = p.adjust(p_wilcox, method = "BH"))
  write.csv(hiv_mod_stats, file.path(SUPP_DIR, "SuppCP018_E_module_stats.csv"),
            row.names = FALSE)
  print(as.data.frame(hiv_mod_stats))
}

# -----------------------------------------------------------------------------
# Manuscript numbers
# -----------------------------------------------------------------------------
stats_out <- list(
  total_cells            = ncol(cp018),
  cells_2m               = sum(cp018$Timepoint == "2m"),
  cells_90m              = sum(cp018$Timepoint == "90m"),
  cells_with_TCR         = sum(!is.na(cp018$CTstrict)),
  hiv_specific_cells     = sum(cp018$HIV_Specific_TCR == "HIV-Specific TCR"),
  hiv_specific_clonotypes= n_distinct(cp018$CTstrict[
                             cp018$HIV_Specific_TCR == "HIV-Specific TCR" &
                             !is.na(cp018$CTstrict)]),
  matched_in_unstim      = if (have_crossref)
                             sum(tara_cp018$clone_match == "HIV-Specific Match") else NA
)
stats_df <- data.frame(metric = names(stats_out),
                       value  = unlist(stats_out), row.names = NULL)
write.csv(stats_df, file.path(FIG_DIR, "CP018_manuscript_numbers.csv"), row.names = FALSE)

cat("\n=== CP018 figure numbers ===\n"); print(stats_df)
cat("\nFigures in:\n  ", FIG_DIR, "\n  ", SUPP_DIR, "\n")
