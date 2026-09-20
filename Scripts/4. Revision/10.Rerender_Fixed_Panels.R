################################################################################
# Re-render ONLY the three panels fixed for the reviewers
#
#   Fig 1E   Rev 1 minor 3 / Rev 2 minor 4  - cropped y-axis title
#   S2C      Rev 2 minor 8                  - first bar's count label clipped
#   S4C      Rev 2 minor 16                 - facets named, not cluster numbers
#
# WHY THIS EXISTS RATHER THAN RE-RUNNING THE THREE FIGURE SCRIPTS
#   1. All three hardcode  base_dir <- "~/Documents/CD8_Longitudinal", which
#      does not exist on this machine (the repo lives under
#      ~/Desktop/Projects/Repositories/). They would fail or write nowhere.
#   2. Figure_2.R also reads 10x contigs from "~/Documents/10x_Genomics/" and
#      runs combineTCR + Trex before it reaches S2C - minutes of work that has
#      nothing to do with the panel being fixed.
#   3. Supliment_4_9.R makes 12 marker/DGE calls for other panels.
#   The blocks below are copied verbatim from those scripts (same data, same
#   aesthetics) with only the reviewer fixes applied, so the panels are
#   reproducible without re-deriving anything.
#
# RUN
#   cd ~/Desktop/Projects/Repositories/CD8_Longitudinal
#   Rscript "Scripts/4. Revision/10.Rerender_Fixed_Panels.R"
#
# Needs Seurat, qs2, ggplot2, dplyr, tidyr, rstatix, ggpubr, patchwork,
# cowplot. It does NOT need scRepertoire, Trex, SeuratExtend or monocle3.
################################################################################

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(ggplot2); library(dplyr)
  library(tidyr);  library(rstatix); library(ggpubr); library(patchwork)
  library(cowplot)
})

# ── Paths: resolved from the repo, not hardcoded ──────────────────────────────
REPO <- normalizePath(file.path(dirname(dirname(getwd())), basename(getwd())),
                      mustWork = FALSE)
if (!dir.exists(file.path(REPO, "saved_R_data"))) REPO <- normalizePath(getwd())
if (!dir.exists(file.path(REPO, "saved_R_data")))
  REPO <- path.expand("~/Desktop/Projects/Repositories/CD8_Longitudinal")
stopifnot(dir.exists(file.path(REPO, "saved_R_data")))

SAVED <- file.path(REPO, "saved_R_data")
FIG1  <- file.path(REPO, "Manuscript", "Figure 1")
SUPP2 <- file.path(REPO, "Manuscript", "Supplementary 2")
SUPP4 <- file.path(REPO, "Manuscript", "Supplementary 4")
CSV   <- file.path(REPO, "Analysis", "CSV_outputs")
invisible(lapply(c(FIG1, SUPP2, SUPP4, CSV), dir.create,
                 recursive = TRUE, showWarnings = FALSE))
cat("repo:  ", REPO, "\n")

# Keep the originals recoverable: the reviewer response cites these panels, so
# a bad re-render must not destroy the version the reviewer saw.
BACKUP <- file.path(REPO, "Manuscript", "_pre_revision_panels")
dir.create(BACKUP, recursive = TRUE, showWarnings = FALSE)
for (f in c(file.path(FIG1,  "Fig1E_TEMRA_Proportion.png"),
            file.path(SUPP2, "S2C_CrossCluster_ClonalSharing.png"),
            file.path(SUPP4, "S4C_ADT_gating.png"))) {
  if (file.exists(f) && !file.exists(file.path(BACKUP, basename(f))))
    file.copy(f, file.path(BACKUP, basename(f)))
}
cat("originals backed up to Manuscript/_pre_revision_panels/\n\n")

################################################################################
# PANEL 1 + 2  -  need TARA_ALL_annotated_final.qs2
################################################################################
cat("loading TARA_ALL_annotated_final.qs2 (1.7 GB, ~1 min)...\n")
TARA_ALL <- qs_read(file.path(SAVED, "TARA_ALL_annotated_final.qs2"))

# --- shared definitions, verbatim from Figure_1.R -----------------------------
cd8_cluster_names <- c("Naive 1 CD8", "Naive 2 CD8", "Naive Intermediate CD8",
                       "Tscm CD8", "TEMRA/CTL CD8", "Tex CD8")
cd8_short_labels <- setNames(
  c("Naive 1", "Naive 2", "Naive Intermediate", "Tscm", "TEMRA/CTL", "Tex"),
  cd8_cluster_names)
cond_cols <- c("HUU" = "#4E79A7", "HEU" = "#F28E2B", "HEI" = "#E15759")

stopifnot(all(c("Annotation", "Condition", "Timepoint_Group",
                "orig.ident", "clonalFrequency") %in%
              colnames(TARA_ALL@meta.data)))

################################################################################
# FIG 1E  -  Rev 1 minor 3 / Rev 2 minor 4
################################################################################
cat("\n=== Fig 1E: TEMRA/CTL proportion ===\n")

TARA_sub <- subset(TARA_ALL, subset =
                     (Condition == "HEU") |
                     (Condition == "HUU") |
                     (Condition == "HEI" & Timepoint_Group == "PreART_Entry"))
TARA_sub$Condition <- factor(TARA_sub$Condition, levels = c("HUU", "HEU", "HEI"))
meta_sub <- as.data.frame(TARA_sub@meta.data)

sample_totals <- meta_sub %>%
  group_by(orig.ident, Condition) %>%
  summarise(total_pbmc = n(), .groups = "drop")

cluster_freq <- meta_sub %>%
  filter(Annotation %in% cd8_cluster_names) %>%
  group_by(orig.ident, Condition, Annotation) %>%
  summarise(cluster_cells = n(), .groups = "drop") %>%
  group_by(orig.ident, Condition) %>%
  complete(Annotation = cd8_cluster_names, fill = list(cluster_cells = 0)) %>%
  ungroup() %>%
  left_join(sample_totals, by = c("orig.ident", "Condition")) %>%
  mutate(percent = 100 * cluster_cells / total_pbmc,
         Facet_Label = factor(cd8_short_labels[as.character(Annotation)],
                              levels = cd8_short_labels))

facet_max <- cluster_freq %>%
  group_by(Facet_Label) %>%
  summarise(max_percent = max(percent, na.rm = TRUE), .groups = "drop")

stat_df_all <- cluster_freq %>%
  group_by(Facet_Label) %>%
  wilcox_test(percent ~ Condition) %>%
  left_join(facet_max, by = "Facet_Label") %>%
  mutate(
    xmin = case_when(group1 == "HUU" & group2 == "HEU" ~ 1,
                     group1 == "HUU" & group2 == "HEI" ~ 1,
                     group1 == "HEU" & group2 == "HEI" ~ 2),
    xmax = case_when(group1 == "HUU" & group2 == "HEU" ~ 2,
                     group1 == "HUU" & group2 == "HEI" ~ 3,
                     group1 == "HEU" & group2 == "HEI" ~ 3),
    y.position = case_when(
      group1 == "HUU" & group2 == "HEU" ~ max_percent * 1.10,
      group1 == "HUU" & group2 == "HEI" ~ max_percent * 1.22,
      group1 == "HEU" & group2 == "HEI" ~ max_percent * 1.16),
    p.signif = case_when(p <= 0.0001 ~ "****", p <= 0.001 ~ "***",
                         p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "ns")) %>%
  ungroup()

temra_freq <- cluster_freq %>% filter(Facet_Label == "TEMRA/CTL")
temra_stat <- stat_df_all %>% filter(p.signif != "ns", Facet_Label == "TEMRA/CTL")

p_temra <- ggplot(temra_freq, aes(x = Condition, y = percent,
                                  fill = Condition, color = Condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.8) +
  stat_pvalue_manual(temra_stat, label = "p.signif", tip.length = 0.01,
                     bracket.size = 0.6, size = 6) +
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  labs(x = NULL, y = "% of total PBMC", title = "TEMRA/CTL CD8") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.text.x  = element_text(face = "bold", size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        # THE FIX: at 5 in wide a 16 pt bold rotated y-title runs past the left
        # edge and is cropped. A left margin plus a wider canvas gives it room.
        plot.margin  = margin(6, 8, 6, 16))

ggsave(file.path(FIG1, "Fig1E_TEMRA_Proportion.png"), p_temra,
       width = 5.9, height = 6, dpi = 600, bg = "white")
write.csv(temra_freq, file.path(CSV, "Fig1E_TEMRA_proportion_data.csv"),
          row.names = FALSE)
cat("  wrote Fig1E_TEMRA_Proportion.png (", nrow(temra_freq), "samples )\n")

################################################################################
# S2C  -  Rev 2 minor 8
################################################################################
cat("\n=== S2C: cross-cluster clonotype sharing ===\n")

expanded_meta <- TARA_ALL@meta.data %>%
  filter(Annotation %in% cd8_cluster_names) %>%
  filter(!is.na(clonalFrequency) & clonalFrequency > 1) %>%
  filter(Condition == "HEI") %>%
  mutate(Cluster = cd8_short_labels[as.character(Annotation)],
         Cluster = factor(Cluster, levels = cd8_short_labels))

clone_cluster_span <- expanded_meta %>%
  group_by(CTstrict) %>%
  summarise(n_clusters = n_distinct(Cluster),
            clusters   = paste(sort(unique(as.character(Cluster))), collapse = " + "),
            n_cells    = n(), .groups = "drop")

cat("  clones in 1 cluster: ", sum(clone_cluster_span$n_clusters == 1), "\n")
cat("  clones in 2 clusters:", sum(clone_cluster_span$n_clusters == 2), "\n")
cat("  clones in 3+ clusters:", sum(clone_cluster_span$n_clusters >= 3), "\n")

p_s2c_span <- ggplot(clone_cluster_span, aes(x = factor(n_clusters))) +
  geom_bar(fill = "#4E79A7", width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5, size = 6, fontface = "bold") +
  # THE FIX: expand = c(0, 0) with limits = c(0, NA) put the upper limit at
  # exactly the tallest bar, so ONLY that bar's label fell outside the panel.
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, NA)) +
  labs(x = "Clusters per Clone", y = "N Clonotypes", title = "Cluster span") +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text  = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"))

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
      prs <- combn(cl_list, 2)
      for (j in seq_len(ncol(prs))) {
        cooccur_mat[prs[1, j], prs[2, j]] <- cooccur_mat[prs[1, j], prs[2, j]] + 1
        cooccur_mat[prs[2, j], prs[1, j]] <- cooccur_mat[prs[2, j], prs[1, j]] + 1
      }
    }
  }
  cooccur_df <- as.data.frame(as.table(cooccur_mat)) %>%
    rename(Cluster1 = Var1, Cluster2 = Var2, n_shared = Freq) %>%
    filter(as.numeric(Cluster1) <= as.numeric(Cluster2))

  p_s2c_cooccur <- ggplot(cooccur_df, aes(x = Cluster1, y = Cluster2, fill = n_shared)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(n_shared > 0, n_shared, "")),
              size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#E15759", name = "Shared\nClones") +
    labs(title = "Cluster co-occurrence", x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(plot.title  = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12))

  ggsave(file.path(SUPP2, "S2C_CrossCluster_ClonalSharing.png"),
         p_s2c_span | p_s2c_cooccur, width = 12, height = 5, dpi = 300, bg = "white")
  write.csv(cooccur_df, file.path(CSV, "S2C_cluster_cooccurrence.csv"), row.names = FALSE)
} else {
  ggsave(file.path(SUPP2, "S2C_CrossCluster_ClonalSharing.png"),
         p_s2c_span, width = 6, height = 5, dpi = 300, bg = "white")
}
write.csv(clone_cluster_span, file.path(CSV, "S2C_clone_cluster_span.csv"),
          row.names = FALSE)
cat("  wrote S2C_CrossCluster_ClonalSharing.png\n")

rm(TARA_ALL, TARA_sub, meta_sub); invisible(gc(verbose = FALSE))

################################################################################
# S4C  -  Rev 2 minor 16
################################################################################
cat("\n=== S4C: ADT sub-gating of naive clusters ===\n")
cat("loading TARA_cd8_HEI_annotated_final.qs2 (332 MB)...\n")
TARA_cd8 <- qs_read(file.path(SAVED, "TARA_cd8_HEI_annotated_final.qs2"))

naive_clusters <- c("Na\u00efve CD8", "Na\u00efve CD8 2", "Na\u00efve CD8 3",
                    "Na\u00efve CD8 4", "Tscm CD8", "Na\u00efve Intermediate CD8")
naive_mask <- TARA_cd8$CD8_Annotation %in% naive_clusters
stopifnot(sum(naive_mask) > 0)

naive_cells <- colnames(TARA_cd8)[naive_mask]
DefaultAssay(TARA_cd8) <- "ADT"
# `slot=` is defunct in SeuratObject >= 5.0.0 (the original script predates
# that), but `layer=` does not exist on older versions. Try the modern form and
# fall back, so this runs on either.
adt_data <- tryCatch(
  GetAssayData(TARA_cd8, layer = "data"),
  error = function(e) GetAssayData(TARA_cd8, slot = "data"))
stopifnot(all(c("CD45RO", "FAS") %in% rownames(adt_data)))

gate_plot_df <- data.frame(
  CD45RO = as.numeric(adt_data["CD45RO", naive_cells]),
  FAS    = as.numeric(adt_data["FAS", naive_cells]),
  annotation = TARA_cd8$CD8_Annotation[naive_cells],
  stringsAsFactors = FALSE)

gate_plot_df$plot_group <- dplyr::case_when(
  gate_plot_df$annotation == "Tscm CD8" ~ "Tscm CD8",
  gate_plot_df$annotation == "Na\u00efve Intermediate CD8" ~ "Na\u00efve Intermediate CD8",
  TRUE ~ "Na\u00efve (ungated)")
gate_plot_df$plot_group <- factor(
  gate_plot_df$plot_group,
  levels = c("Na\u00efve (ungated)", "Tscm CD8", "Na\u00efve Intermediate CD8"))

if ("source_cluster" %in% colnames(TARA_cd8@meta.data)) {
  gate_plot_df$source <- TARA_cd8$source_cluster[naive_cells]
  # THE FIX: facets read "Source cluster 0 / 2 / 4 / 6". A cluster NUMBER means
  # nothing to a reader and the manuscript uses numbers nowhere else. Label with
  # the subset name(s) the cluster resolves to, number kept in parentheses.
  .src_names <- tapply(as.character(gate_plot_df$annotation),
                       as.character(gate_plot_df$source),
                       function(v) paste(sort(unique(v)), collapse = " / "))
  gate_plot_df$facet <- paste0(.src_names[as.character(gate_plot_df$source)],
                               "\n(source cluster ", gate_plot_df$source, ")")
  cat("  facet labels now:\n")
  for (f in unique(gate_plot_df$facet)) cat("    ", gsub("\n", "  ", f), "\n")
} else {
  gate_plot_df$facet <- "All na\u00efve clusters"
  cat("  NOTE: no source_cluster column - single facet\n")
}

# VERBATIM from Supliment_4_9.R: the original filters the ungated cells OUT
# and plots only the two gated populations, at point size 2.5 and font 22. The
# reviewer asked ONLY about the facet labels, so nothing else changes - a
# different point size or a different cell set would make this panel
# inconsistent with the published figure set.
gate_cols <- c("Tscm CD8" = "#E31A1C", "Na\u00efve Intermediate CD8" = "#1F78B4")
fg_df <- gate_plot_df %>% dplyr::filter(plot_group != "Na\u00efve (ungated)")

p_s4c <- ggplot(fg_df, aes(x = FAS, y = CD45RO, color = plot_group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = gate_cols, name = "Population",
                     labels = c("Tscm CD8" = "Tscm (CD45RO\u2212 FAS+)",
                                "Na\u00efve Intermediate CD8" = "Intermediate (CD45RO+)")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  facet_wrap(~ facet, nrow = 1) +
  labs(title = "ADT sub-gating of na\u00efve clusters",
       x = "FAS/CD95 (DSB)", y = "CD45RO (DSB)") +
  theme_cowplot(font_size = 22) +
  theme(strip.background = element_rect(fill = "grey95", color = NA),
        # only change beyond the facet labels: the strips now carry NAMES, not
        # single digits, so they need a smaller face and room for two lines.
        strip.text = element_text(size = 13, face = "bold", lineheight = 1.05),
        legend.position = "bottom", legend.text = element_text(size = 20),
        plot.background = element_rect(fill = "white", color = NA)) +
  guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))

ggsave(file.path(SUPP4, "S4C_ADT_gating.png"), plot = p_s4c,
       width = 20, height = 7, dpi = 300, bg = "white")
cat("  gated cells plotted:", nrow(fg_df), "of", nrow(gate_plot_df), "naive cells\n")
write.csv(fg_df, file.path(CSV, "S4C_ADT_gating_data.csv"), row.names = FALSE)
cat("  wrote S4C_ADT_gating.png (", nrow(fg_df), "naive cells )\n")

################################################################################
# Verify: every panel written, and no ink on the canvas edge
################################################################################
cat("\n=== written files ===\n")
for (f in c(file.path(FIG1,  "Fig1E_TEMRA_Proportion.png"),
            file.path(SUPP2, "S2C_CrossCluster_ClonalSharing.png"),
            file.path(SUPP4, "S4C_ADT_gating.png")))
  cat(sprintf("  %-42s %s  %.0f KB\n", basename(f),
              if (file.exists(f)) "OK " else "MISSING",
              if (file.exists(f)) file.size(f) / 1024 else 0))

cat("\nOriginals kept in Manuscript/_pre_revision_panels/ - delete once happy.\n")
cat("Then tell the agent and it will fold these into TARA_Revised_Figures.pdf\n")
cat("and switch the four responses (R1-m3, R2-m4, R2-m8, R2-m16) to past tense.\n")
