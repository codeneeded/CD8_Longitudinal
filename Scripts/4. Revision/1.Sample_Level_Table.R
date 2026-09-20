################################################################################
# Revision script 1 — sample-level cell counts and clinical table
#
# ANSWERS
#   Rev 2 major 1a  "How many CD8 T cells in total? Provide total PBMC and
#                    sub-setted CD8 counts per participant"
#   Rev 2 minor 1   "Supplementary table with participant-level clinical info
#                    and classification of each timepoint sampled"
#   Rev 2 minor 10  "Fig 3B: number of cells analysed, per-sample breakdown"
#   Editor          "number of data points stated clearly in all figure legends"
#
# OUTPUT
#   Tables/TableS1_sample_level_counts.csv   <- one row PER SAMPLE (not per PID)
#   Tables/TableS1_participant_summary.csv
#   Tables/Cell_counts_by_lineage_per_sample.csv
#   Figures/R2m2_01_cells_per_sample.png
#
# The reviewer explicitly suggests converting Table 1 into a sample-level table,
# one line per sample. That is what TableS1_sample_level_counts.csv is.
################################################################################

# Works whether run interactively (RStudio) or via Rscript. `sys.frame()$ofile`
# is NULL in RStudio's console, so resolve the repo instead of the script dir.
.cfg <- Sys.getenv("CD8_REV_CONFIG", "")
if (!nzchar(.cfg)) {
  for (b in c("~/Documents/CD8_Longitudinal",
              "~/Desktop/Projects/Repositories/CD8_Longitudinal")) {
    p <- file.path(path.expand(b), "Scripts", "4. Revision", "config_revision.R")
    if (file.exists(p)) { .cfg <- p; break }
  }
}
source(.cfg)

################################################################################
# 1. Whole-PBMC object: per-sample totals and lineage breakdown
################################################################################
require_obj(QS_ALL)
all_obj <- qs_read(QS_ALL)
md <- all_obj@meta.data

cat("\nWhole-PBMC object:", ncol(all_obj), "cells,", length(unique(md$orig.ident)),
    "samples,", length(unique(md$PID)), "participants\n")

# Azimuth level-1 lineage is the defensible way to answer "how many CD8 T cells"
# without depending on the manual annotation the reviewer is questioning.
lineage_tab <- md %>%
  count(orig.ident, PID, Timepoint_Group, Condition, Cohort, Age,
        lineage = predicted.celltype.l1, name = "n_cells") %>%
  tidyr::pivot_wider(names_from = lineage, values_from = n_cells, values_fill = 0)

per_sample <- md %>%
  group_by(Sample = orig.ident) %>%
  summarise(
    PID              = first(PID),
    Age_months       = first(Age),
    Timepoint_Group  = first(Timepoint_Group),
    Condition        = first(Condition),
    Cohort           = first(Cohort),
    Viral_Load       = first(Viral_Load),
    n_PBMC_QC_pass   = n(),
    n_CD8_T_azimuth  = sum(predicted.celltype.l1 == "CD8 T"),
    n_CD4_T_azimuth  = sum(predicted.celltype.l1 == "CD4 T"),
    n_B_azimuth      = sum(predicted.celltype.l1 == "B"),
    n_NK_azimuth     = sum(predicted.celltype.l1 == "NK"),
    n_Mono_azimuth   = sum(predicted.celltype.l1 == "Mono"),
    n_DC_azimuth     = sum(predicted.celltype.l1 == "DC"),
    n_with_TCR       = sum(has_TCR %in% c(TRUE, "TRUE"), na.rm = TRUE),
    median_nFeature  = round(median(nFeature_RNA)),
    median_nCount    = round(median(nCount_RNA)),
    median_pct_mito  = round(median(percent.mt), 2),
    .groups = "drop"
  ) %>%
  mutate(pct_CD8_of_PBMC = round(100 * n_CD8_T_azimuth / n_PBMC_QC_pass, 1))

################################################################################
# 2. Dedicated CD8 object: per-sample counts in the object Figures 4/5 use
################################################################################
require_obj(QS_CD8)
cd8 <- qs_read(QS_CD8)
cmd <- cd8@meta.data

cat("Dedicated CD8 object:", ncol(cd8), "cells (HEI only),",
    length(unique(cmd$PID)), "participants\n")

cd8_per_sample <- cmd %>%
  group_by(Sample = orig.ident) %>%
  summarise(
    n_CD8_subclustered = n(),
    n_CD8_with_TCR     = sum(has_TCR %in% c(TRUE, "TRUE"), na.rm = TRUE),
    n_CD8_expanded     = sum(clonalFrequency > 1, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(pct_expanded = round(100 * n_CD8_expanded / n_CD8_with_TCR, 1))

sample_table <- per_sample %>%
  left_join(cd8_per_sample, by = "Sample") %>%
  arrange(PID, Age_months)

save_tab("TableS1_sample_level_counts", sample_table)

# Participant-level roll-up, for the Table 1 revision
participant_table <- sample_table %>%
  group_by(PID) %>%
  summarise(
    Condition        = first(Condition),
    Cohort           = first(Cohort),
    n_timepoints     = n(),
    timepoints       = paste(sort(Age_months), collapse = ", "),
    groups_sampled   = paste(sort(unique(Timepoint_Group)), collapse = "; "),
    total_PBMC       = sum(n_PBMC_QC_pass),
    total_CD8_T      = sum(n_CD8_T_azimuth),
    total_CD8_subclust = sum(n_CD8_subclustered, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(total_PBMC))
save_tab("TableS1_participant_summary", participant_table)

save_tab("Cell_counts_by_lineage_per_sample", lineage_tab)

################################################################################
# 3. Figure: cells recovered per sample, so N is visible not tabulated
################################################################################
plot_df <- sample_table %>%
  select(Sample, PID, Timepoint_Group, n_PBMC_QC_pass, n_CD8_T_azimuth) %>%
  tidyr::pivot_longer(c(n_PBMC_QC_pass, n_CD8_T_azimuth),
                      names_to = "compartment", values_to = "n") %>%
  mutate(compartment = recode(compartment,
           n_PBMC_QC_pass  = "All PBMC (QC-pass)",
           n_CD8_T_azimuth = "CD8 T cells"))

p <- ggplot(plot_df, aes(x = reorder(Sample, -n), y = n, fill = compartment)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("All PBMC (QC-pass)" = "#999999",
                               "CD8 T cells" = "#D55E00"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = scales::comma) +
  labs(
    title = "Cells recovered per sample after quality control",
    subtitle = paste0(nrow(sample_table), " samples from ",
                      length(unique(sample_table$PID)), " participants; ",
                      format(sum(sample_table$n_PBMC_QC_pass), big.mark = ","),
                      " PBMC total, ",
                      format(sum(sample_table$n_CD8_T_azimuth), big.mark = ","),
                      " CD8 T cells"),
    x = NULL, y = "Cells",
    caption = "CD8 T cells assigned by Azimuth reference mapping (level 1).") +
  theme_rev(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

save_fig("R2m2_01_cells_per_sample", p, width = 11, height = 5.5)
save_src("R2m2_01_cells_per_sample", plot_df)

################################################################################
# 4. Report the headline numbers the revision text needs
################################################################################
cat("\n=== NUMBERS FOR THE RESPONSE LETTER ===\n")
cat("  samples:                    ", nrow(sample_table), "\n")
cat("  participants:               ", length(unique(sample_table$PID)), "\n")
cat("  total PBMC (QC-pass):       ", format(sum(sample_table$n_PBMC_QC_pass), big.mark=","), "\n")
cat("  total CD8 T (Azimuth l1):   ", format(sum(sample_table$n_CD8_T_azimuth), big.mark=","), "\n")
cat("  CD8 in dedicated object:    ", format(sum(sample_table$n_CD8_subclustered, na.rm=TRUE), big.mark=","), "\n")
cat("  median CD8 per sample:      ", median(sample_table$n_CD8_T_azimuth), "\n")
cat("  range CD8 per sample:       ", min(sample_table$n_CD8_T_azimuth), "-",
    max(sample_table$n_CD8_T_azimuth), "\n")

cat("\n  Discrepancy to explain in Methods: the dedicated CD8 object holds\n")
cat("  ", ncol(cd8), " cells but Azimuth calls ", sum(md$predicted.celltype.l1=="CD8 T"),
    " CD8 T in the full object.\n", sep = "")
cat("  Difference is HEI-only restriction plus the manual-annotation gates in\n")
cat("  4.CD8_Subseting.R. State both.\n")
