############################################################
# CP018 FACT-Seq — technical replicate concordance (A vs B)
#
# NEW for CP018: both timepoints were run in duplicate (2m = OB1/OB2,
# 90m = OB3/OB4), so replicate reproducibility can be measured directly.
# CP003 had replicates only at 2m and none at 101m, so this was not possible.
#
# Establishes, BEFORE merging replicates, that A and B agree on:
#   1. cell recovery and QC distributions
#   2. pseudobulk expression (per-gene correlation)
#   3. cell-type composition (Azimuth)
#   4. cell cycle phase composition
#
# Input : QS_SINGLETS  (from 1.QC_Full.R)
# Output: CONC_ROOT plots + CSVs
############################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(qs2)
library(Matrix)

mkdirs(CONC_ROOT)

require_file(QS_SINGLETS, "post-QC singlets object")
seu <- qs_read(QS_SINGLETS)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# -----------------------------
# 1) Cell recovery per replicate
# -----------------------------
recovery <- seu@meta.data %>%
  count(Timepoint, Replicate, Sample_Subfolder, name = "n_cells") %>%
  group_by(Timepoint) %>%
  mutate(
    pct_of_timepoint = round(100 * n_cells / sum(n_cells), 1),
    # ratio of the smaller to larger replicate: 1.0 = perfect balance
    balance_ratio    = round(min(n_cells) / max(n_cells), 3)
  ) %>%
  ungroup()

write.csv(recovery, file.path(CONC_ROOT, "CP018_replicate_cell_recovery.csv"),
          row.names = FALSE)
print(recovery)

p_recovery <- ggplot(recovery, aes(x = Timepoint, y = n_cells, fill = Replicate)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_text(aes(label = n_cells), position = position_dodge(0.8),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = REP_COLORS) +
  labs(title = "CP018 cell recovery by technical replicate",
       y = "Cells passing QC", x = NULL) +
  theme_classic(base_size = 12)

ggsave(file.path(CONC_ROOT, "Cell_recovery_by_replicate.png"), p_recovery,
       width = 6, height = 4.5, dpi = 400, bg = "white")

# -----------------------------
# 2) QC metric distributions, A vs B
# -----------------------------
qc_feats <- intersect(c("nCount_RNA", "nFeature_RNA", "percent_mito",
                        "percent_ribo", "log10GenesPerUMI"), colnames(seu@meta.data))

qc_long <- seu@meta.data %>%
  select(Timepoint, Replicate, all_of(qc_feats)) %>%
  pivot_longer(all_of(qc_feats), names_to = "metric", values_to = "value")

p_qc <- ggplot(qc_long, aes(x = Replicate, y = value, fill = Replicate)) +
  geom_violin(scale = "width", alpha = 0.75, colour = "grey20") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  facet_grid(metric ~ Timepoint, scales = "free_y") +
  scale_fill_manual(values = REP_COLORS) +
  labs(title = "QC metrics by technical replicate", x = NULL, y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "none")

ggsave(file.path(CONC_ROOT, "QC_metrics_by_replicate.png"), p_qc,
       width = 7, height = 10, dpi = 400, bg = "white")

# Wilcoxon A vs B within each timepoint (descriptive: with thousands of cells,
# tiny differences reach significance — report effect size alongside)
qc_tests <- qc_long %>%
  group_by(Timepoint, metric) %>%
  summarise(
    median_A = median(value[Replicate == "A"], na.rm = TRUE),
    median_B = median(value[Replicate == "B"], na.rm = TRUE),
    pct_diff = round(100 * (median_B - median_A) / median_A, 2),
    p_wilcox = tryCatch(
      wilcox.test(value[Replicate == "A"], value[Replicate == "B"])$p.value,
      error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH"))

write.csv(qc_tests, file.path(CONC_ROOT, "CP018_replicate_QC_tests.csv"),
          row.names = FALSE)
print(qc_tests)

# -----------------------------
# 3) Pseudobulk expression correlation
# -----------------------------
# Sum raw counts per library, CPM-normalise, log1p, then correlate A vs B
counts <- LayerData(seu, assay = "RNA", layer = "counts")
libs   <- unique(seu$Sample_Subfolder)

pb <- sapply(libs, function(l) {
  Matrix::rowSums(counts[, which(seu$Sample_Subfolder == l), drop = FALSE])
})
pb_cpm <- log1p(sweep(pb, 2, colSums(pb), "/") * 1e6)

# keep genes with signal in at least one library
pb_cpm <- pb_cpm[rowSums(pb_cpm) > 0, , drop = FALSE]

cor_pearson  <- cor(pb_cpm, method = "pearson")
cor_spearman <- cor(pb_cpm, method = "spearman")

write.csv(cor_pearson,  file.path(CONC_ROOT, "CP018_pseudobulk_cor_pearson.csv"))
write.csv(cor_spearman, file.path(CONC_ROOT, "CP018_pseudobulk_cor_spearman.csv"))

cat("\nPseudobulk Pearson correlation (log1p CPM):\n"); print(round(cor_pearson, 4))

# Scatter for each replicate pair
for (tp in unique(CP018_SAMPLES$Timepoint)) {
  a <- paste0("CP018_", tp, "_A"); b <- paste0("CP018_", tp, "_B")
  if (!all(c(a, b) %in% colnames(pb_cpm))) next

  df <- data.frame(A = pb_cpm[, a], B = pb_cpm[, b])
  r_p <- cor(df$A, df$B, method = "pearson")
  r_s <- cor(df$A, df$B, method = "spearman")

  p <- ggplot(df, aes(A, B)) +
    geom_point(alpha = 0.15, size = 0.5, colour = TP_COLORS[[tp]]) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title = paste0("CP018 ", tp, " — pseudobulk replicate concordance"),
         subtitle = sprintf("Pearson r = %.4f   Spearman rho = %.4f   (%d genes)",
                            r_p, r_s, nrow(df)),
         x = paste0(tp, " replicate A  (log1p CPM)"),
         y = paste0(tp, " replicate B  (log1p CPM)")) +
    theme_classic(base_size = 12)

  p <- p + theme(plot.margin = margin(10, 18, 10, 10),
                 plot.title.position = "plot")
  ggsave(file.path(CONC_ROOT, paste0("Pseudobulk_scatter_", tp, ".png")), p,
         width = 5.5, height = 5.5, dpi = 400, bg = "white")
}

# -----------------------------
# 4) Composition concordance (Azimuth + Phase)
# -----------------------------
comp_concordance <- function(col, label) {
  if (!col %in% colnames(seu@meta.data)) {
    message("Skipping ", label, " — column '", col, "' not present"); return(invisible(NULL))
  }

  comp <- seu@meta.data %>%
    count(Timepoint, Replicate, !!rlang::sym(col), name = "n") %>%
    rename(Category = !!rlang::sym(col)) %>%
    group_by(Timepoint, Replicate) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  write.csv(comp, file.path(CONC_ROOT, paste0("CP018_composition_", label, ".csv")),
            row.names = FALSE)

  p_bar <- ggplot(comp, aes(x = Replicate, y = pct, fill = Category)) +
    geom_col(colour = "black", linewidth = 0.2) +
    facet_wrap(~ Timepoint) +
    labs(title = paste0("CP018 composition by replicate — ", label),
         y = "% of cells", x = NULL) +
    theme_classic(base_size = 11)
  ggsave(file.path(CONC_ROOT, paste0("Composition_stacked_", label, ".png")), p_bar,
         width = 8, height = 5.5, dpi = 400, bg = "white")

  # A vs B percentage scatter, per timepoint
  wide <- comp %>%
    select(Timepoint, Replicate, Category, pct) %>%
    pivot_wider(names_from = Replicate, values_from = pct, values_fill = 0)

  if (all(c("A", "B") %in% names(wide))) {
    lab_df <- wide %>% group_by(Timepoint) %>%
      summarise(r = cor(A, B, method = "pearson"), .groups = "drop")

    p_sc <- ggplot(wide, aes(A, B, colour = Timepoint)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(size = 2.5, alpha = 0.85) +
      scale_colour_manual(values = TP_COLORS) +
      labs(title = paste0("Replicate composition agreement — ", label),
           subtitle = paste(sprintf("%s: r = %.3f", lab_df$Timepoint, lab_df$r),
                            collapse = "   |   "),
           x = "% in replicate A", y = "% in replicate B") +
      theme_classic(base_size = 12)

    ggsave(file.path(CONC_ROOT, paste0("Composition_scatter_", label, ".png")), p_sc,
           width = 6, height = 5, dpi = 400, bg = "white")

    write.csv(lab_df, file.path(CONC_ROOT,
              paste0("CP018_composition_correlation_", label, ".csv")), row.names = FALSE)
    print(lab_df)
  }
}

comp_concordance("predicted.celltype.l2", "Azimuth_l2")
comp_concordance("Phase", "CellCyclePhase")

# -----------------------------
# 5) Verdict
# -----------------------------
# Pooling replicates is justified when pseudobulk r is high (>0.95 is typical
# for technical replicates) and composition tracks the diagonal. Record the
# numbers so the decision is documented rather than assumed.
verdict <- data.frame(
  Timepoint = character(), pearson_r = numeric(),
  spearman_rho = numeric(), balance_ratio = numeric(),
  stringsAsFactors = FALSE
)
for (tp in unique(CP018_SAMPLES$Timepoint)) {
  a <- paste0("CP018_", tp, "_A"); b <- paste0("CP018_", tp, "_B")
  if (!all(c(a, b) %in% colnames(pb_cpm))) next
  verdict <- rbind(verdict, data.frame(
    Timepoint     = tp,
    pearson_r     = round(cor_pearson[a, b], 4),
    spearman_rho  = round(cor_spearman[a, b], 4),
    balance_ratio = recovery$balance_ratio[recovery$Timepoint == tp][1]
  ))
}
write.csv(verdict, file.path(CONC_ROOT, "CP018_replicate_concordance_verdict.csv"),
          row.names = FALSE)

# full per-cell metadata for this stage
export_metadata(seu, file.path(CONC_ROOT, "CP018_cell_metadata_postQC.csv"))

cat("\n=== Replicate concordance summary ===\n")
print(verdict)
cat("\nPooling replicates is reasonable when pearson_r is high and composition\n",
    "scatters sit on the diagonal. Inspect the plots in:\n  ", CONC_ROOT, "\n")
