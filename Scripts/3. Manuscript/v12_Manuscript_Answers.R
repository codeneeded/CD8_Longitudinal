################################################################################
# v12 MANUSCRIPT ANSWERS — compute the data points needed to resolve Lesley's
# open comments and Akshay's outstanding text placeholders.
#
# Run this and paste the printed values + the contents of v12_manuscript_answers.csv
# back to me so I can apply them to v12.
#
# This script answers five questions for v12:
#
#   Q1. Pre-/post-ART contraction % for Tscm CD8 and Naïve Intermediate CD8
#       (Lesley Comment 466 — "both of which resolved after ART initiation")
#
#   Q2. Exact Spearman ρ and p-value for the Stemness ~ viral load correlation
#       (Lesley Comment 570 — "not statistically significant")
#       Reports per-cluster results for Tscm, Naïve Intermediate, and the full
#       naïve compartment; the manuscript text refers to "stemness erosion with
#       increasing viral load", so we report whichever cluster the original
#       claim was anchored on plus the aggregate.
#
#   Q3 & Q4. Empirical VL cutoffs used to classify samples as PostART_Suppressed
#       vs PostART_Unsuppressed (the "VL<xxx" placeholders in Results §3.2).
#       Derived from the existing Timepoint_Group labels in the metadata.
#
#   Q5. Confirm the population that was sub-gated to produce Tscm and
#       Naïve Intermediate ("Sub-gating of XXX") — pulled directly from the
#       Supliment_4_9.R workflow (the six naïve-lineage clusters).
################################################################################

library(Seurat); library(dplyr); library(tidyr); library(qs2)

# ── Paths (match existing scripts) ─────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
analysis_dir <- file.path(base_dir, "Analysis", "CD8_subset")
out_dir      <- file.path(base_dir, "Manuscript", "v12_answers")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_path <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")
vl_path     <- file.path(analysis_dir, "09_viral_load", "PerSample_ModuleScores_and_ViralLoad.csv")

# Container for all the results we'll save to a single CSV at the end
answers <- list()

# ── Load CD8 Seurat object ────────────────────────────────────────────────────
cat("Loading TARA_cd8 Seurat object...\n")
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("  Loaded:", ncol(TARA_cd8), "cells\n")

art_levels <- c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed")
TARA_cd8$Timepoint_Group <- factor(TARA_cd8$Timepoint_Group, levels = art_levels)

################################################################################
# Q1 — Tscm CD8 and Naïve Intermediate CD8 proportions pre- vs post-ART
################################################################################
cat("\n=========================================================\n")
cat("Q1. Tscm and Naïve Intermediate proportions across ART status\n")
cat("=========================================================\n")

cluster_pop <- c("Tscm CD8", "Naïve Intermediate CD8")

# Two complementary views — both reported because the manuscript phrasing
# could mean either:
#   (a) "X% pre-ART → Y% post-ART of THIS cluster's cells"   (cluster composition)
#   (b) "the cluster represents X% of all CD8 pre-ART → Y% post-ART" (abundance)

# (a) Cluster composition: of all cells in a given cluster, what fraction came from each ART status?
comp_df <- TARA_cd8@meta.data %>%
  filter(CD8_Annotation %in% cluster_pop) %>%
  group_by(CD8_Annotation, Timepoint_Group) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(CD8_Annotation) %>%
  mutate(pct_of_cluster = round(n_cells / sum(n_cells) * 100, 1)) %>%
  ungroup() %>%
  select(CD8_Annotation, Timepoint_Group, n_cells, pct_of_cluster) %>%
  arrange(CD8_Annotation, Timepoint_Group)

cat("\n(a) Cluster composition (what fraction of each cluster's cells came from each ART status)\n")
print(comp_df)

# (b) Cluster abundance: of all CD8 cells in a given ART status, what fraction is this cluster?
abund_df <- TARA_cd8@meta.data %>%
  group_by(Timepoint_Group, CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Timepoint_Group) %>%
  mutate(pct_of_total_cd8 = round(n / sum(n) * 100, 2)) %>%
  ungroup() %>%
  filter(CD8_Annotation %in% cluster_pop) %>%
  select(CD8_Annotation, Timepoint_Group, n, pct_of_total_cd8) %>%
  arrange(CD8_Annotation, Timepoint_Group)

cat("\n(b) Cluster abundance (fraction of all CD8 cells at each ART status that belong to this cluster)\n")
print(abund_df)

# Fisher's exact: is the contraction post-ART significant for each cluster?
cat("\nFisher's exact: cluster x (Pre-ART vs Post-ART) — pooling Suppressed + Unsuppressed as Post\n")
fisher_q1 <- list()
for (cl in cluster_pop) {
  is_cl  <- TARA_cd8$CD8_Annotation == cl
  is_pre <- TARA_cd8$Timepoint_Group == "PreART_Entry"
  mat <- matrix(c(
    sum(is_cl &  is_pre),  sum(is_cl & !is_pre),
    sum(!is_cl &  is_pre), sum(!is_cl & !is_pre)
  ), nrow = 2, byrow = TRUE,
     dimnames = list(c("Cluster", "Other"), c("PreART", "PostART")))
  ft <- fisher.test(mat)
  cat(sprintf("  %s: OR = %.2f (95%% CI %.2f–%.2f), p = %.3g\n",
              cl, ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
  fisher_q1[[cl]] <- data.frame(
    question = "Q1_fisher",
    cluster = cl,
    odds_ratio_pre_vs_post = round(ft$estimate, 3),
    ci_low = round(ft$conf.int[1], 3),
    ci_high = round(ft$conf.int[2], 3),
    p_value = ft$p.value
  )
}

answers$Q1_composition  <- comp_df  %>% mutate(question = "Q1_composition")
answers$Q1_abundance    <- abund_df %>% mutate(question = "Q1_abundance")
answers$Q1_fisher       <- bind_rows(fisher_q1)

################################################################################
# Q2 — Stemness ~ viral load correlation (Spearman, exact p-value)
################################################################################
cat("\n=========================================================\n")
cat("Q2. Stemness module score ~ plasma viral load (Spearman)\n")
cat("=========================================================\n")

stemness_results <- data.frame()

# Path A: per-sample summary file (preferred — matches Figure 5E pipeline)
if (file.exists(vl_path)) {
  sample_scores <- read.csv(vl_path)
  if (!"log10_VL" %in% colnames(sample_scores)) {
    sample_scores$log10_VL <- log10(sample_scores$Viral_Load_num + 1)
  }

  if ("Stemness" %in% colnames(sample_scores)) {
    ct_overall <- cor.test(sample_scores$log10_VL, sample_scores$Stemness,
                           method = "spearman", exact = FALSE)
    cat(sprintf("\nOverall (per-sample mean Stemness vs log10 VL):\n"))
    cat(sprintf("  n = %d samples,  rho = %.3f,  p = %.4g\n",
                nrow(sample_scores), ct_overall$estimate, ct_overall$p.value))

    stemness_results <- rbind(stemness_results, data.frame(
      question = "Q2_stemness_VL",
      level = "per_sample_overall",
      cluster = "all_cells",
      n = nrow(sample_scores),
      spearman_rho = round(ct_overall$estimate, 4),
      p_value = ct_overall$p.value,
      n_below_p_0.05 = as.integer(ct_overall$p.value < 0.05),
      n_below_p_0.10 = as.integer(ct_overall$p.value < 0.10)
    ))
  } else {
    cat("\nNOTE: 'Stemness' column not found in PerSample file. Available:\n")
    cat("  ", paste(colnames(sample_scores), collapse = ", "), "\n")
  }
}

# Path B: per-cell module scores aggregated per sample for specific clusters
if (file.exists(module_path)) {
  mod_scores <- read.csv(module_path, row.names = 1)
  cat(sprintf("\nPer-cell module scores: %d cells, columns: %s\n",
              nrow(mod_scores), paste(colnames(mod_scores), collapse = ", ")))

  if ("Stemness" %in% colnames(mod_scores) &&
      "Viral_Load_num" %in% colnames(mod_scores)) {

    # Per-sample mean stemness within each naïve / stem-like cluster
    naive_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
                        "Naïve CD8 4", "Tscm CD8", "Naïve Intermediate CD8")

    for (cl in c(naive_clusters, "ALL_NAIVE_COMPARTMENT")) {
      cl_data <- if (cl == "ALL_NAIVE_COMPARTMENT") {
        mod_scores %>% filter(CD8_Annotation %in% naive_clusters)
      } else {
        mod_scores %>% filter(CD8_Annotation == cl)
      }
      per_sample <- cl_data %>%
        group_by(orig.ident, Timepoint_Group) %>%
        summarise(mean_stemness = mean(Stemness, na.rm = TRUE),
                  mean_VL = mean(Viral_Load_num, na.rm = TRUE),
                  n_cells = n(), .groups = "drop") %>%
        filter(!is.na(mean_VL) & n_cells >= 5)
      if (nrow(per_sample) >= 4) {
        per_sample$log10_VL <- log10(per_sample$mean_VL + 1)
        ct <- cor.test(per_sample$log10_VL, per_sample$mean_stemness,
                       method = "spearman", exact = FALSE)
        cat(sprintf("  %-30s  n=%2d  rho=%6.3f  p=%.4g\n",
                    cl, nrow(per_sample), ct$estimate, ct$p.value))
        stemness_results <- rbind(stemness_results, data.frame(
          question = "Q2_stemness_VL",
          level = "per_sample_per_cluster",
          cluster = cl,
          n = nrow(per_sample),
          spearman_rho = round(ct$estimate, 4),
          p_value = ct$p.value,
          n_below_p_0.05 = as.integer(ct$p.value < 0.05),
          n_below_p_0.10 = as.integer(ct$p.value < 0.10)
        ))
      }
    }
  } else {
    cat("\nNOTE: Stemness or Viral_Load_num column missing in per-cell module scores file.\n")
    cat("  Available:", paste(colnames(mod_scores), collapse = ", "), "\n")
  }
}

if (nrow(stemness_results) > 0) answers$Q2_stemness_VL <- stemness_results

################################################################################
# Q3 / Q4 — VL cutoffs for suppressed and unsuppressed classification
################################################################################
cat("\n=========================================================\n")
cat("Q3/Q4. Empirical VL cutoffs from Timepoint_Group labels\n")
cat("=========================================================\n")

# Try to locate a viral-load-per-cell or per-sample column on the Seurat object
md <- TARA_cd8@meta.data
vl_col_candidates <- c("Viral_Load_num", "ViralLoad_num", "VL_num",
                       "Viral_Load", "ViralLoad", "VL", "log10_VL")
vl_col <- vl_col_candidates[vl_col_candidates %in% colnames(md)][1]

vl_summary <- NULL
if (!is.na(vl_col)) {
  cat(sprintf("\nUsing VL column from metadata: '%s'\n", vl_col))
  # Per-sample (orig.ident) VL — one value per sample
  vl_per_sample <- md %>%
    filter(!is.na(.data[[vl_col]])) %>%
    group_by(orig.ident, PID, Timepoint_Group) %>%
    summarise(VL = mean(.data[[vl_col]], na.rm = TRUE), .groups = "drop")
  vl_summary <- vl_per_sample %>%
    group_by(Timepoint_Group) %>%
    summarise(
      n_samples = n(),
      VL_min  = min(VL, na.rm = TRUE),
      VL_med  = median(VL, na.rm = TRUE),
      VL_max  = max(VL, na.rm = TRUE),
      .groups = "drop"
    )
  cat("\nPer-sample VL (one per orig.ident) by Timepoint_Group:\n")
  print(vl_summary)
  cat("\nFull per-sample VL listing:\n")
  print(vl_per_sample %>% arrange(Timepoint_Group, VL))
} else {
  # Fall back to the per-sample VL CSV used by Figure 5E
  if (file.exists(vl_path)) {
    sample_scores <- read.csv(vl_path)
    cat(sprintf("\nNo VL column in Seurat metadata. Reading from per-sample CSV: %s\n", vl_path))
    if ("Viral_Load_num" %in% colnames(sample_scores) &&
        "Timepoint_Group" %in% colnames(sample_scores)) {
      vl_summary <- sample_scores %>%
        group_by(Timepoint_Group) %>%
        summarise(
          n_samples = n(),
          VL_min  = min(Viral_Load_num, na.rm = TRUE),
          VL_med  = median(Viral_Load_num, na.rm = TRUE),
          VL_max  = max(Viral_Load_num, na.rm = TRUE),
          .groups = "drop"
        )
      cat("\nPer-sample VL by Timepoint_Group:\n")
      print(vl_summary)
      cat("\nFull per-sample VL listing:\n")
      print(sample_scores %>% select(any_of(c("PID","orig.ident","Timepoint_Group","Viral_Load_num"))) %>%
              arrange(Timepoint_Group, Viral_Load_num))
    }
  }
}

if (!is.null(vl_summary)) {
  cat("\nInferred cutoff for Suppressed = MAX VL among PostART_Suppressed samples.\n")
  cat("Inferred cutoff for Unsuppressed = MIN VL among PostART_Unsuppressed samples.\n")
  cat("If they straddle a clean threshold (e.g., 50, 200, 1000), that is the cutoff used.\n")
  answers$Q3_Q4_VL_cutoffs <- vl_summary %>% mutate(question = "Q3_Q4_VL_cutoffs")
}

################################################################################
# Q5 — What population was sub-gated? (Comment 1000 "XXX" placeholder)
################################################################################
cat("\n=========================================================\n")
cat("Q5. Sub-gating parent population (XXX placeholder)\n")
cat("=========================================================\n")
cat("From Supliment_4_9.R, the FAS x CD45RO sub-gating in Supplementary Fig S4C\n")
cat("is applied to the union of these six clusters:\n")
naive_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3", "Naïve CD8 4",
                    "Tscm CD8", "Naïve Intermediate CD8")
print(naive_clusters)
cat("\nSuggested replacement for 'Sub-gating of XXX' in the manuscript:\n")
cat("  'Sub-gating of the naïve CD8 compartment by surface CD45RA, CD45RO and FAS/CD95'\n")
cat("OR more specifically:\n")
cat("  'Sub-gating of Naïve CD8 1-4 by surface CD45RA, CD45RO and FAS/CD95'\n")

# Per-cluster sub-gating counts within naïve compartment, for reference
adt_assay <- "ADT"
if (adt_assay %in% names(TARA_cd8@assays)) {
  DefaultAssay(TARA_cd8) <- adt_assay
  if (all(c("CD45RO", "FAS") %in% rownames(TARA_cd8))) {
    naive_mask <- TARA_cd8$CD8_Annotation %in% naive_clusters
    gate_df <- data.frame(
      cell = colnames(TARA_cd8)[naive_mask],
      cluster = TARA_cd8$CD8_Annotation[naive_mask],
      CD45RO = GetAssayData(TARA_cd8, slot = "data")["CD45RO", naive_mask],
      FAS = GetAssayData(TARA_cd8, slot = "data")["FAS", naive_mask]
    )
    gate_df$gate <- with(gate_df, ifelse(CD45RO > 0, "Naïve Intermediate-like (CD45RO+)",
                                  ifelse(FAS > 0, "Tscm-like (CD45RO- FAS+)",
                                                  "Other naïve (CD45RO- FAS-)")))
    cat("\nADT gate counts within the naïve compartment:\n")
    print(table(gate_df$cluster, gate_df$gate))
  }
}

DefaultAssay(TARA_cd8) <- "RNA"

################################################################################
# Save all answers to a single CSV
################################################################################
answers_long <- bind_rows(lapply(answers, function(df) {
  df$value_pairs <- apply(df, 1, function(row) {
    paste(names(row), "=", row, collapse = " | ")
  })
  df %>% select(question, value_pairs)
}))

write.csv(answers_long, file.path(out_dir, "v12_manuscript_answers.csv"), row.names = FALSE)

# Save individual data frames for clarity too
for (nm in names(answers)) {
  write.csv(answers[[nm]], file.path(out_dir, paste0("v12_", nm, ".csv")), row.names = FALSE)
}

cat(sprintf("\n=== Done. Results written to: %s ===\n", out_dir))
cat("Files:\n")
print(list.files(out_dir))

cat("\n=========================================================\n")
cat("SUMMARY FOR AKSHAY — paste these values back to me:\n")
cat("=========================================================\n")
cat("  Q1 — Tscm CD8 pre-ART %    :   <fill from Q1_composition>\n")
cat("  Q1 — Tscm CD8 post-ART %   :   <fill from Q1_composition>\n")
cat("  Q1 — Naïve Int. pre-ART %  :   <fill from Q1_composition>\n")
cat("  Q1 — Naïve Int. post-ART % :   <fill from Q1_composition>\n")
cat("  Q2 — Stemness ~ VL rho/p   :   <fill from Q2_stemness_VL>\n")
cat("  Q3 — VL_suppressed cutoff  :   <max VL in PostART_Suppressed>\n")
cat("  Q4 — VL_unsuppressed cutoff:   <min VL in PostART_Unsuppressed>\n")
cat("  Q5 — XXX = the naïve CD8 compartment (Naïve CD8 1-4)\n")
cat("=========================================================\n")

