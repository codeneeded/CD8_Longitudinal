################################################################################
# CP003 vs CP018 — the suppressed / unsuppressed FACT-Seq comparison
#
# THE QUESTION
#   Both participants had HIV-peptide-stimulated FACT-Seq at two timepoints.
#   They differ in ART response within the sampled window:
#
#     CP003  NEVER suppressed   VL >500,000 at both 12m and 24m
#     CP018  SUPPRESSED         pre-ART VL 176,970 -> VL 20 from day 65
#
#   The manuscript's Figure 4/5 claim is that suppression resolves exhaustion
#   and preserves stemness. That was shown across the unstimulated cohort.
#   Here it is tested in ANTIGEN-VALIDATED cells from two individuals.
#
# WHAT THIS SCRIPT DOES *NOT* DO
#   It does not merge or co-cluster the two datasets. n = 1 per group, so
#   cross-participant differences are DESCRIPTIVE. Any difference also carries
#   a Cell Ranger version component (CP003 = 9.0.1, CP018 = 10.1.0) and
#   different stimulation timepoints (2m/101m vs 2m/90m). Treat every number
#   here as hypothesis-generating, not as a statistical test between groups.
#   The per-cell Wilcoxon tests below describe the two cell populations that
#   were observed; they are NOT tests of a participant-level effect.
#
# Input : QS_ANNOTATED (CP018, from script 5), QS_CP003_ANNOT
# Output: CMP_ROOT plots + tables
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggrepel)

CMP_ROOT <- file.path(CP018_ANALYSIS, "CP003_vs_CP018")
csv_dir  <- file.path(CMP_ROOT, "tables")
mkdirs(CMP_ROOT, csv_dir)

# ── Load both ────────────────────────────────────────────────────────────────
require_file(QS_ANNOTATED,   "CP018 annotated object (script 5)")
# Prefer the object that carries Fig3_Annotation - without it the MAIT
# exclusion cannot be applied to CP003 and the "matched" recall is not matched.
cp003_path <- if (file.exists(QS_CP003_ANNOT_FIG3)) QS_CP003_ANNOT_FIG3 else QS_CP003_ANNOT
require_file(cp003_path, "CP003 annotated object")
message("CP003 object: ", basename(cp003_path))

cp18 <- qs_read(QS_ANNOTATED)
cp3  <- qs_read(cp003_path)
if (!"Fig3_Annotation" %in% colnames(cp3@meta.data))
  warning("CP003 object lacks Fig3_Annotation - MAIT cells will NOT be excluded ",
          "from its HIV-specific recall; treat that count as an upper bound.")

DefaultAssay(cp18) <- "RNA"; cp18 <- JoinLayers(cp18)
DefaultAssay(cp3)  <- "RNA"; cp3  <- JoinLayers(cp3)

message("CP018: ", ncol(cp18), " cells | CP003: ", ncol(cp3), " cells")

# ── Participant metadata (from the manuscript Table 1) ───────────────────────
participants <- data.frame(
  participant       = c("CP003", "CP018"),
  ART_status        = c("Never suppressed", "Suppressed"),
  preART_VL         = c(NA, 176970),
  VL_at_suppression = c(NA, 20),
  days_to_suppress  = c(NA, 65),
  stim_timepoints   = c("2m, 101m", "2m, 90m"),
  cellranger        = c("9.0.1", "10.1.0"),
  n_cells           = c(ncol(cp3), ncol(cp18))
)
save_csv(participants, file.path(csv_dir, "participant_metadata.csv"))

# ── Signature scores, computed IDENTICALLY on both objects ───────────────────
# Recomputed here rather than reusing each object's stored scores, because
# AddModuleScore control-gene sampling differs between runs. Same genes, same
# method, both datasets -> the only honest way to compare.
SIGS <- list(
  Stemness     = c("TCF7","LEF1","SELL","CCR7","IL7R","ID3","BACH2"),
  Exhaustion   = c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1"),
  Cytotoxicity = c("GZMB","GNLY","PRF1","NKG7","GZMH","FGFBP2","KLRG1"),
  IFN_I        = c("IFIT1","IFIT3","ISG15","MX1","OAS1","STAT1","IFI6"),
  Activation   = c("CD38","HLA-DRA","MKI67","TNFRSF9","IL2RA")
)

sig_scores <- function(obj, label) {
  X <- GetAssayData(obj, layer = "data")
  out <- lapply(names(SIGS), function(s) {
    g <- intersect(SIGS[[s]], rownames(X))
    Matrix::colMeans(X[g, , drop = FALSE])
  })
  names(out) <- names(SIGS)
  df <- as.data.frame(out)
  df$cell_barcode <- colnames(obj)
  df$participant  <- label
  md <- obj@meta.data
  df$Timepoint      <- as.character(md$Timepoint)
  df$HIV_Specific   <- if ("HIV_Specific_TCR" %in% names(md))
                         as.character(md$HIV_Specific_TCR) else NA_character_
  df$clonalFrequency <- if ("clonalFrequency" %in% names(md)) md$clonalFrequency else NA
  df$CTstrict       <- if ("CTstrict" %in% names(md)) as.character(md$CTstrict) else NA_character_
  df$Phase          <- as.character(md$Phase)
  df
}

sc <- bind_rows(sig_scores(cp3, "CP003"), sig_scores(cp18, "CP018"))
sc$participant <- factor(sc$participant, levels = c("CP003", "CP018"))
save_csv(sc, file.path(csv_dir, "signature_scores_per_cell.csv"))

# Genes actually found in both objects — report so a missing gene is visible.
gene_report <- do.call(rbind, lapply(names(SIGS), function(s) data.frame(
  signature = s, gene = SIGS[[s]],
  in_CP003  = SIGS[[s]] %in% rownames(cp3),
  in_CP018  = SIGS[[s]] %in% rownames(cp18))))
save_csv(gene_report, file.path(csv_dir, "signature_gene_coverage.csv"))

# ── 1) Whole-sample signature comparison ─────────────────────────────────────
long <- sc %>% pivot_longer(all_of(names(SIGS)),
                            names_to = "signature", values_to = "score")

sig_summary <- long %>%
  group_by(signature, participant) %>%
  summarise(n = n(), mean = mean(score), median = median(score),
            sd = sd(score), .groups = "drop") %>%
  mutate(across(c(mean, median, sd), ~round(.x, 4)))
save_csv(sig_summary, file.path(csv_dir, "signature_summary_by_participant.csv"))

# Descriptive comparison of the observed cell populations (see header caveat).
sig_tests <- long %>%
  group_by(signature) %>%
  summarise(
    mean_CP003 = round(mean(score[participant == "CP003"]), 4),
    mean_CP018 = round(mean(score[participant == "CP018"]), 4),
    difference = round(mean(score[participant == "CP018"]) -
                       mean(score[participant == "CP003"]), 4),
    cohens_d   = round((mean(score[participant == "CP018"]) -
                        mean(score[participant == "CP003"])) /
                       sqrt((var(score[participant == "CP018"]) +
                             var(score[participant == "CP003"])) / 2), 3),
    wilcox_p   = wilcox.test(score ~ participant)$p.value,
    .groups = "drop") %>%
  mutate(wilcox_p_BH = p.adjust(wilcox_p, method = "BH"),
         note = "cell-level descriptive comparison; n=1 participant per group")
save_csv(sig_tests, file.path(csv_dir, "signature_comparison.csv"))
cat("\n=== Signature comparison (CP018 suppressed - CP003 unsuppressed) ===\n")
print(as.data.frame(sig_tests), row.names = FALSE)

p_sig <- ggplot(long, aes(x = participant, y = score, fill = participant)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(CP003 = "#B23A48", CP018 = "#1F78B4")) +
  labs(x = NULL, y = "Mean normalised expression",
       title = "Signature scores: unsuppressed (CP003) vs suppressed (CP018)",
       subtitle = "All stimulated cells. n = 1 participant per group - descriptive only.") +
  theme_cp018(base_size = 15) + theme(legend.position = "none")
save_png(file.path(CMP_ROOT, "Signature_scores_by_participant.png"),
         p_sig, width = 16, height = 5.5)

# ── 2) The key contrast: HIV-specific cells only ─────────────────────────────
hiv <- sc %>% filter(HIV_Specific == "HIV-Specific TCR")

if (nrow(hiv) > 0 && dplyr::n_distinct(hiv$participant) == 2) {
  hiv_long <- hiv %>% pivot_longer(all_of(names(SIGS)),
                                   names_to = "signature", values_to = "score")
  hiv_tests <- hiv_long %>%
    group_by(signature) %>%
    summarise(n_CP003 = sum(participant == "CP003"),
              n_CP018 = sum(participant == "CP018"),
              mean_CP003 = round(mean(score[participant == "CP003"]), 4),
              mean_CP018 = round(mean(score[participant == "CP018"]), 4),
              difference = round(mean(score[participant == "CP018"]) -
                                 mean(score[participant == "CP003"]), 4),
              wilcox_p = wilcox.test(score ~ participant)$p.value,
              .groups = "drop") %>%
    mutate(wilcox_p_BH = p.adjust(wilcox_p, method = "BH"))
  save_csv(hiv_tests, file.path(csv_dir, "signature_comparison_HIVspecific.csv"))
  cat("\n=== HIV-specific cells only ===\n")
  print(as.data.frame(hiv_tests), row.names = FALSE)

  p_hiv <- ggplot(hiv_long, aes(x = participant, y = score, fill = participant)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.4) +
    facet_wrap(~ signature, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(CP003 = "#B23A48", CP018 = "#1F78B4")) +
    labs(x = NULL, y = "Mean normalised expression",
         title = "HIV-specific cells: unsuppressed vs suppressed",
         subtitle = "Antigen-validated clonotypes only") +
    theme_cp018(base_size = 15) + theme(legend.position = "none")
  save_png(file.path(CMP_ROOT, "Signature_scores_HIVspecific.png"),
           p_hiv, width = 16, height = 5.5)
} else {
  message("HIV_Specific_TCR missing in one object - skipping the HIV-specific contrast.")
}

# ── 3) Exhaustion / stemness balance per cluster ─────────────────────────────
# The single clearest read: does any cluster in either participant sit in the
# high-exhaustion / low-stemness corner (= terminal exhaustion)?
cp3_cl  <- if ("Fig3_Annotation" %in% colnames(cp3@meta.data))
             as.character(cp3$Fig3_Annotation) else as.character(cp3$mnn_clusters_rna)
cp18_cl <- if ("Fig_Annotation" %in% colnames(cp18@meta.data))
             as.character(cp18$Fig_Annotation) else as.character(cp18$mnn_clusters_rna)
clus_bal <- bind_rows(
  cbind(data.frame(participant = "CP003", cluster = cp3_cl, stringsAsFactors = FALSE),
        sig_scores(cp3, "CP003")[, c("Stemness","Exhaustion")]),
  cbind(data.frame(participant = "CP018", cluster = cp18_cl, stringsAsFactors = FALSE),
        sig_scores(cp18, "CP018")[, c("Stemness","Exhaustion")])
) %>%
  group_by(participant, cluster) %>%
  summarise(n_cells = n(),
            stemness   = round(mean(Stemness), 4),
            exhaustion = round(mean(Exhaustion), 4),
            ratio_ex_stem = round(mean(Exhaustion) / mean(Stemness), 4),
            .groups = "drop") %>%
  arrange(participant, desc(exhaustion))
save_csv(clus_bal, file.path(csv_dir, "cluster_exhaustion_stemness.csv"))
cat("\n=== Cluster exhaustion/stemness balance ===\n")
print(as.data.frame(clus_bal), row.names = FALSE)

p_bal <- ggplot(clus_bal, aes(x = stemness, y = exhaustion,
                              colour = participant, size = n_cells)) +
  geom_point(alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = cluster), size = 4.2, show.legend = FALSE,
                           max.overlaps = 20) +
  scale_colour_manual(values = c(CP003 = "#B23A48", CP018 = "#1F78B4")) +
  scale_size_continuous(range = c(3, 11)) +
  labs(x = "Stemness score", y = "Exhaustion score",
       title = "Cluster position in exhaustion / stemness space",
       subtitle = "Top-left = terminal exhaustion. Point size = cluster size.") +
  theme_cp018(base_size = 15)
save_png(file.path(CMP_ROOT, "Cluster_exhaustion_vs_stemness.png"),
         p_bal, width = 10, height = 7.5)

# ── 4) Repertoire comparison ─────────────────────────────────────────────────
rep_stats <- function(obj, label) {
  md <- obj@meta.data
  ct <- md$CTstrict[!is.na(md$CTstrict) & md$CTstrict != ""]
  tb <- table(ct)
  hiv_n <- if ("HIV_Specific_TCR" %in% names(md))
    dplyr::n_distinct(md$CTstrict[md$HIV_Specific_TCR == "HIV-Specific TCR" &
                                  !is.na(md$CTstrict)]) else NA
  data.frame(
    participant        = label,
    n_cells            = ncol(obj),
    n_cells_with_TCR   = length(ct),
    pct_with_TCR       = round(100 * length(ct) / ncol(obj), 1),
    n_clonotypes       = length(tb),
    n_expanded         = sum(tb > 1),
    pct_expanded_cells = round(100 * sum(tb[tb > 1]) / length(ct), 1),
    largest_clone      = if (length(tb)) max(tb) else NA,
    # Chao1 richness and Shannon diversity, comparable across unequal depth
    shannon            = round(-sum((tb/sum(tb)) * log(tb/sum(tb))), 3),
    n_HIV_specific_clonotypes = hiv_n   # stored call; see matched recall in section 4b
  )
}
rep_cmp <- bind_rows(rep_stats(cp3, "CP003"), rep_stats(cp18, "CP018"))
save_csv(rep_cmp, file.path(csv_dir, "repertoire_comparison.csv"))
cat("\n=== Repertoire comparison ===\n")
print(as.data.frame(rep_cmp), row.names = FALSE)

################################################################################
# 4b. Matched HIV-specific recall — same rule applied to BOTH participants
################################################################################
# The published CP003 call used "clonalFrequency > 1 AND not in bystander
# clusters" (hand-picked cluster numbers). CP018 now uses expansion only
# (config: HIV_MIN_CLONE_SIZE, HIV_PER_TIMEPOINT). Comparing 179 clonotypes
# called one way against ~73 called another is not a comparison, so both are
# recalled here under the CP018 rule. The published CP003 numbers are retained
# alongside so the effect of the rule change is visible rather than hidden.
recall_hiv <- function(obj, label) {
  m <- obj@meta.data
  ct <- as.character(m$CTstrict)
  tp <- as.character(m$Timepoint)
  keep <- !is.na(ct) & ct != ""

  if (isTRUE(HIV_PER_TIMEPOINT)) {
    key <- paste(ct, tp, sep = "||")
    sz  <- table(key[keep]); csize <- as.integer(sz[key])
  } else {
    sz  <- table(ct[keep]);  csize <- as.integer(sz[ct])
  }
  csize[!keep] <- NA_integer_

  # MAIT exclusion by LABEL - never by cluster number (that was the CP018 bug).
  anno_col <- intersect(c("Fig_Annotation", "Fig3_Annotation"), colnames(m))
  is_mait <- rep(FALSE, nrow(m))
  if (isTRUE(HIV_EXCLUDE_MAIT) && length(anno_col)) {
    is_mait <- grepl("MAIT", as.character(m[[anno_col[1]]]), ignore.case = TRUE)
  }
  called <- !is.na(csize) & csize > HIV_MIN_CLONE_SIZE & !is_mait

  data.frame(
    participant          = label,
    rule                 = paste0("clone size > ", HIV_MIN_CLONE_SIZE,
                                  if (isTRUE(HIV_PER_TIMEPOINT)) " within timepoint" else " pooled"),
    n_cells              = nrow(m),
    n_cells_with_TCR     = sum(keep),
    n_HIV_specific_cells = sum(called),
    n_HIV_specific_clonotypes = dplyr::n_distinct(ct[called]),
    n_MAIT_excluded      = sum(is_mait),
    mait_exclusion_possible = length(anno_col) > 0,
    stringsAsFactors = FALSE)
}

matched <- bind_rows(recall_hiv(cp3, "CP003"), recall_hiv(cp18, "CP018"))
save_csv(matched, file.path(csv_dir, "HIV_specific_matched_recall.csv"))
cat("\n=== HIV-specific, SAME rule applied to both participants ===\n")
print(as.data.frame(matched), row.names = FALSE)

if (any(!matched$mait_exclusion_possible)) {
  cat("\nNOTE: no annotation column found for ",
      paste(matched$participant[!matched$mait_exclusion_possible], collapse = ", "),
      " - MAIT cells were NOT excluded there. The recall is therefore not fully\n",
      "matched; treat that participant's count as an upper bound.\n", sep = "")
}

# Published-definition counts, for comparison with the matched recall above.
published <- data.frame(
  participant = c("CP003", "CP018"),
  published_definition = c("clonalFrequency > 1 AND not in bystander clusters",
                           "clonalFrequency > 1 AND not in bystander clusters (INVALID - cluster numbers came from a different resolution)"),
  n_clonotypes_published = c(
    if ("HIV_Specific_TCR" %in% colnames(cp3@meta.data))
      dplyr::n_distinct(cp3$CTstrict[cp3$HIV_Specific_TCR == "HIV-Specific TCR" &
                                     !is.na(cp3$CTstrict)]) else NA,
    if ("HIV_Specific_TCR" %in% colnames(cp18@meta.data))
      dplyr::n_distinct(cp18$CTstrict[cp18$HIV_Specific_TCR == "HIV-Specific TCR" &
                                      !is.na(cp18$CTstrict)]) else NA),
  stringsAsFactors = FALSE)
published$n_clonotypes_matched_rule <- matched$n_HIV_specific_clonotypes[
  match(published$participant, matched$participant)]
save_csv(published, file.path(csv_dir, "HIV_specific_definition_comparison.csv"))
cat("\n=== Published vs matched definition ===\n")
print(published[, c("participant","n_clonotypes_published","n_clonotypes_matched_rule")],
      row.names = FALSE)

# ── 5) Shared clonotypes between participants ────────────────────────────────
# Any overlap is either a public TCR specificity or a barcode artefact — both
# worth knowing about before anything is claimed.
cl3  <- unique(cp3$CTstrict[!is.na(cp3$CTstrict)])
cl18 <- unique(cp18$CTstrict[!is.na(cp18$CTstrict)])
shared <- intersect(cl3, cl18)
cat("\nClonotypes shared between CP003 and CP018:", length(shared), "\n")
# data.frame() recycles a length-1 column against a length-0 one and errors, so
# the empty case needs its own branch. Zero overlap is the EXPECTED result for
# two unrelated individuals - it is evidence against a barcode/index artefact.
shared_df <- if (length(shared) > 0) {
  data.frame(CTstrict = shared,
             note = "public specificity OR barcode artefact - inspect before use",
             stringsAsFactors = FALSE)
} else {
  data.frame(CTstrict = character(0), note = character(0), stringsAsFactors = FALSE)
}
save_csv(shared_df, file.path(csv_dir, "shared_clonotypes_CP003_CP018.csv"))
if (length(shared) == 0) {
  cat("  -> no overlap, as expected between unrelated participants;\n",
      "     also rules out barcode bleed-through between the two datasets.\n", sep = "")
}

# ── 6) Cell-cycle / activation composition ───────────────────────────────────
# table() will not recycle a length-1 label against a length-n vector, so build
# the participant column at full length first.
phase_cmp <- bind_rows(
  data.frame(participant = "CP003", Phase = as.character(cp3$Phase),
             stringsAsFactors = FALSE),
  data.frame(participant = "CP018", Phase = as.character(cp18$Phase),
             stringsAsFactors = FALSE)
) %>%
  count(participant, Phase, name = "Freq") %>%
  group_by(participant) %>%
  mutate(pct = round(100 * Freq / sum(Freq), 1)) %>%
  ungroup()
save_csv(phase_cmp, file.path(csv_dir, "cellcycle_composition.csv"))

p_phase <- ggplot(phase_cmp, aes(x = participant, y = pct, fill = Phase)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = paste0(pct, "%")), position = position_stack(vjust = 0.5),
            size = 4.5, fontface = "bold") +
  labs(x = NULL, y = "% of cells", title = "Cell cycle composition") +
  theme_cp018(base_size = 15)
save_png(file.path(CMP_ROOT, "CellCycle_composition.png"), p_phase,
         width = 7, height = 5.5)

# ── 7) One-page summary ──────────────────────────────────────────────────────
summary_tbl <- data.frame(
  metric = c("ART status", "Stimulation timepoints", "Cell Ranger version",
             "Cells", "Cells with TCR", "Clonotypes", "Expanded clonotypes",
             "HIV-specific clonotypes",
             "Mean stemness", "Mean exhaustion", "Mean cytotoxicity",
             "Mean type-I IFN", "Max cluster exhaustion"),
  CP003 = c("Never suppressed", "2m, 101m", "9.0.1",
            rep_cmp$n_cells[1], rep_cmp$n_cells_with_TCR[1],
            rep_cmp$n_clonotypes[1], rep_cmp$n_expanded[1],
            rep_cmp$n_HIV_specific_clonotypes[1],
            sig_summary$mean[sig_summary$signature=="Stemness"     & sig_summary$participant=="CP003"],
            sig_summary$mean[sig_summary$signature=="Exhaustion"   & sig_summary$participant=="CP003"],
            sig_summary$mean[sig_summary$signature=="Cytotoxicity" & sig_summary$participant=="CP003"],
            sig_summary$mean[sig_summary$signature=="IFN_I"        & sig_summary$participant=="CP003"],
            max(clus_bal$exhaustion[clus_bal$participant=="CP003"])),
  CP018 = c("Suppressed (day 65)", "2m, 90m", "10.1.0",
            rep_cmp$n_cells[2], rep_cmp$n_cells_with_TCR[2],
            rep_cmp$n_clonotypes[2], rep_cmp$n_expanded[2],
            rep_cmp$n_HIV_specific_clonotypes[2],
            sig_summary$mean[sig_summary$signature=="Stemness"     & sig_summary$participant=="CP018"],
            sig_summary$mean[sig_summary$signature=="Exhaustion"   & sig_summary$participant=="CP018"],
            sig_summary$mean[sig_summary$signature=="Cytotoxicity" & sig_summary$participant=="CP018"],
            sig_summary$mean[sig_summary$signature=="IFN_I"        & sig_summary$participant=="CP018"],
            max(clus_bal$exhaustion[clus_bal$participant=="CP018"]))
)
save_csv(summary_tbl, file.path(csv_dir, "CP003_vs_CP018_summary.csv"))
cat("\n=== SUMMARY ===\n"); print(summary_tbl, row.names = FALSE)

cat("\nOutputs:\n  ", CMP_ROOT, "\n  ", csv_dir, "\n")
cat("\nREMINDER: n = 1 participant per group. Differences are descriptive and\n",
    "carry a Cell Ranger version confound (9.0.1 vs 10.1.0). State both in Methods.\n", sep = "")
