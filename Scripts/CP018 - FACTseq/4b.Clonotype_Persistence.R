################################################################################
# CP018 — clonotype persistence 2m -> 90m (alluvial)
#
# THE QUESTION
#   Are clonotypes detected at 2m still present at 90m (7.3 years later)?
#
# THE DETECTION PROBLEM — READ THIS BEFORE INTERPRETING ANY OUTPUT
#   Technical replicates of the SAME sample share only 3.1% (2m) and 15.2%
#   (90m) of clonotypes, and 95.8% of clonotypes are observed as single cells.
#   Single-cell TCR sequencing samples a few thousand cells from a repertoire of
#   ~10^6-10^7 clonotypes, so ABSENCE AT 90m IS NOT EVIDENCE OF LOSS. The
#   within-sample replicate overlap is the detection ceiling: cross-timepoint
#   persistence can only be interpreted against it, never in absolute terms.
#
#   This script therefore reports:
#     (a) observed persistence,
#     (b) the replicate-overlap ceiling as the positive control, and
#     (c) persistence stratified by clone size, which is the only comparison
#         with enough detection power to mean anything.
#
# Input : QS_WITH_TCR (script 4) — or QS_ANNOTATED if script 5 has run
# Output: VDJ_ROOT/Persistence — alluvial plots + CSVs
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(scales)

PERS_ROOT <- file.path(VDJ_ROOT, "Persistence")
csv_dir   <- file.path(PERS_ROOT, "tables")
mkdirs(PERS_ROOT, csv_dir)

# Prefer the annotated object (has HIV_Specific_TCR); fall back to the TCR object.
obj_path <- if (file.exists(QS_ANNOTATED)) QS_ANNOTATED else QS_WITH_TCR
require_file(obj_path, "TCR-attached object")
message("Using: ", basename(obj_path))
seu <- qs_read(obj_path)

md <- seu@meta.data
md$cell_barcode <- rownames(md)
has_hiv  <- "HIV_Specific_TCR" %in% names(md)
has_anno <- "Fig_Annotation"   %in% names(md)
rm(seu); gc(verbose = FALSE)          # keep memory free - object is ~180 MB

tcr <- md %>% filter(!is.na(CTstrict), CTstrict != "")
message(nrow(tcr), " cells with a clonotype")

################################################################################
# 1. Per-clonotype summary
################################################################################
clono <- tcr %>%
  group_by(CTstrict) %>%
  summarise(
    CTaa        = dplyr::first(CTaa),
    n_total     = n(),
    n_2m        = sum(Timepoint == "2m"),
    n_90m       = sum(Timepoint == "90m"),
    n_2m_A      = sum(Sample_Subfolder == "CP018_2m_A"),
    n_2m_B      = sum(Sample_Subfolder == "CP018_2m_B"),
    n_90m_A     = sum(Sample_Subfolder == "CP018_90m_A"),
    n_90m_B     = sum(Sample_Subfolder == "CP018_90m_B"),
    clusters    = paste(sort(unique(as.character(mnn_clusters_rna))), collapse = ";"),
    annotation  = if (has_anno)
                    paste(sort(unique(as.character(Fig_Annotation))), collapse = ";") else NA,
    hiv_specific= if (has_hiv)
                    any(HIV_Specific_TCR == "HIV-Specific TCR", na.rm = TRUE) else NA,
    .groups = "drop") %>%
  mutate(
    fate = case_when(
      n_2m > 0 & n_90m > 0 ~ "Persistent (2m + 90m)",
      n_2m > 0             ~ "2m only",
      TRUE                 ~ "90m only"),
    expanded_2m  = n_2m  > 1,
    expanded_90m = n_90m > 1,
    # a clonotype seen in BOTH replicates of a timepoint is not a library artefact
    rep_confirmed_2m  = n_2m_A  > 0 & n_2m_B  > 0,
    rep_confirmed_90m = n_90m_A > 0 & n_90m_B > 0)

save_csv(clono, file.path(csv_dir, "CP018_clonotype_fates.csv"))

################################################################################
# 2. Detection-power controls — the numbers persistence must be read against
################################################################################
rep_overlap <- function(a, b, label) {
  ca <- clono$CTstrict[clono[[a]] > 0]
  cb <- clono$CTstrict[clono[[b]] > 0]
  sh <- length(intersect(ca, cb))
  data.frame(comparison = label, n_A = length(ca), n_B = length(cb),
             n_shared = sh,
             pct_of_smaller = round(100 * sh / min(length(ca), length(cb)), 2),
             jaccard = round(sh / length(union(ca, cb)), 4))
}
ctrl <- bind_rows(
  rep_overlap("n_2m_A",  "n_2m_B",  "2m A vs B (same sample, technical)"),
  rep_overlap("n_90m_A", "n_90m_B", "90m A vs B (same sample, technical)"),
  rep_overlap("n_2m",    "n_90m",   "2m vs 90m (BIOLOGICAL, 88 months apart)")
)
save_csv(ctrl, file.path(csv_dir, "CP018_detection_power_controls.csv"))
cat("\n=== Detection power ===\n"); print(ctrl, row.names = FALSE)
cat("\nThe two technical rows are the CEILING. Cross-timepoint persistence\n",
    "cannot exceed what replicates of the same sample recover.\n", sep = "")

# clone-size distribution — why the ceiling is so low
size_dist <- clono %>% count(n_total, name = "n_clonotypes") %>%
  mutate(pct = round(100 * n_clonotypes / sum(n_clonotypes), 2))
save_csv(size_dist, file.path(csv_dir, "CP018_clone_size_distribution.csv"))
cat(sprintf("\nSingletons: %d / %d clonotypes (%.1f%%)\n",
            size_dist$n_clonotypes[size_dist$n_total == 1], nrow(clono),
            100 * size_dist$n_clonotypes[size_dist$n_total == 1] / nrow(clono)))

################################################################################
# 3. Persistence stratified by clone size — the interpretable comparison
################################################################################
# Detection probability scales with clone size, so persistence is only
# comparable WITHIN a size stratum.
strat <- clono %>%
  filter(n_2m > 0) %>%
  mutate(size_bin = cut(n_2m, breaks = c(0, 1, 2, 4, 9, Inf),
                        labels = c("1", "2", "3-4", "5-9", "10+"))) %>%
  group_by(size_bin) %>%
  summarise(n_clonotypes = n(),
            n_redetected_90m = sum(n_90m > 0),
            pct_redetected = round(100 * mean(n_90m > 0), 2),
            .groups = "drop")
save_csv(strat, file.path(csv_dir, "CP018_persistence_by_clone_size.csv"))
cat("\n=== Persistence by clone size at 2m ===\n"); print(as.data.frame(strat), row.names = FALSE)

p_strat <- ggplot(strat, aes(x = size_bin, y = pct_redetected)) +
  geom_col(fill = "#1F78B4", colour = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = paste0(n_redetected_90m, "/", n_clonotypes)),
            vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Clone size at 2m (cells)", y = "% re-detected at 90m",
       title = "Persistence scales with clonal expansion",
       subtitle = "Larger clones are more likely to be re-sampled 88 months later") +
  theme_cp018()
save_png(file.path(PERS_ROOT, "Persistence_by_clone_size.png"), p_strat,
         width = 8, height = 5.5)

################################################################################
# 4. Alluvial — clonotype flow 2m -> 90m
################################################################################
# Every clonotype is one stratum, categorised at each timepoint by expansion
# state. "Not detected" is an explicit category, NOT an absence, because the
# controls above show non-detection is the expected outcome for a singleton.
size_cat <- function(n) {
  dplyr::case_when(n == 0 ~ "Not detected",
                   n == 1 ~ "Single cell",
                   n <= 4 ~ "Small (2-4)",
                   n <= 9 ~ "Medium (5-9)",
                   TRUE   ~ "Large (10+)")
}
LEVELS <- c("Large (10+)", "Medium (5-9)", "Small (2-4)", "Single cell", "Not detected")

allu <- clono %>%
  transmute(CTstrict,
            `2m`  = factor(size_cat(n_2m),  levels = LEVELS),
            `90m` = factor(size_cat(n_90m), levels = LEVELS),
            fate, hiv_specific) %>%
  count(`2m`, `90m`, fate, name = "n_clonotypes")
save_csv(allu, file.path(csv_dir, "CP018_alluvial_source_data.csv"))

allu_long <- allu %>%
  mutate(alluvium = row_number()) %>%
  pivot_longer(c(`2m`, `90m`), names_to = "Timepoint", values_to = "state") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("2m", "90m")),
         state = factor(state, levels = LEVELS))

STATE_COLS <- c("Large (10+)"  = "#B2182B", "Medium (5-9)" = "#EF8A62",
                "Small (2-4)"  = "#F4A261", "Single cell"  = "#92C5DE",
                "Not detected" = "grey88")

p_allu <- ggplot(allu_long,
                 aes(x = Timepoint, stratum = state, alluvium = alluvium,
                     y = n_clonotypes, fill = state)) +
  geom_flow(alpha = 0.55, width = 0.35, colour = NA) +
  geom_stratum(width = 0.35, colour = "black", linewidth = 0.35) +
  scale_fill_manual(values = STATE_COLS, name = "Clone size") +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(x = NULL, y = "Clonotypes",
       title = "CP018 clonotype flow, 2m to 90m",
       subtitle = paste0("All ", nrow(clono), " clonotypes. 'Not detected' reflects ",
                         "sampling depth, not clonal loss -\nsee detection-power controls.")) +
  theme_cp018(base_size = 15)
save_png(file.path(PERS_ROOT, "Alluvial_all_clonotypes.png"), p_allu,
         width = 9, height = 7)

# The informative view: EXPANDED clonotypes only. Singletons dominate the plot
# above and carry almost no detection power.
exp_clono <- clono %>% filter(n_total > 1)
if (nrow(exp_clono) > 0) {
  allu_e <- exp_clono %>%
    transmute(`2m`  = factor(size_cat(n_2m),  levels = LEVELS),
              `90m` = factor(size_cat(n_90m), levels = LEVELS)) %>%
    count(`2m`, `90m`, name = "n_clonotypes") %>%
    mutate(alluvium = row_number()) %>%
    pivot_longer(c(`2m`, `90m`), names_to = "Timepoint", values_to = "state") %>%
    mutate(Timepoint = factor(Timepoint, levels = c("2m", "90m")),
           state = factor(state, levels = LEVELS))

  p_e <- ggplot(allu_e, aes(x = Timepoint, stratum = state, alluvium = alluvium,
                            y = n_clonotypes, fill = state)) +
    geom_flow(alpha = 0.6, width = 0.35, colour = NA) +
    geom_stratum(width = 0.35, colour = "black", linewidth = 0.35) +
    geom_text(stat = "stratum", aes(label = after_stat(count)), size = 4) +
    scale_fill_manual(values = STATE_COLS, name = "Clone size") +
    labs(x = NULL, y = "Clonotypes",
         title = paste0("Expanded clonotypes only (n = ", nrow(exp_clono), ")"),
         subtitle = "Clonotypes observed in >1 cell - the population with detection power") +
    theme_cp018(base_size = 15)
  save_png(file.path(PERS_ROOT, "Alluvial_expanded_clonotypes.png"), p_e,
           width = 9, height = 7)
}

# HIV-specific clonotypes, if script 5 has run
if (has_hiv) {
  hiv_clono <- clono %>% filter(hiv_specific)
  if (nrow(hiv_clono) > 0) {
    allu_h <- hiv_clono %>%
      transmute(`2m`  = factor(size_cat(n_2m),  levels = LEVELS),
                `90m` = factor(size_cat(n_90m), levels = LEVELS)) %>%
      count(`2m`, `90m`, name = "n_clonotypes") %>%
      mutate(alluvium = row_number()) %>%
      pivot_longer(c(`2m`, `90m`), names_to = "Timepoint", values_to = "state") %>%
      mutate(Timepoint = factor(Timepoint, levels = c("2m", "90m")),
             state = factor(state, levels = LEVELS))

    p_h <- ggplot(allu_h, aes(x = Timepoint, stratum = state, alluvium = alluvium,
                              y = n_clonotypes, fill = state)) +
      geom_flow(alpha = 0.6, width = 0.35, colour = NA) +
      geom_stratum(width = 0.35, colour = "black", linewidth = 0.35) +
      geom_text(stat = "stratum", aes(label = after_stat(count)), size = 4) +
      scale_fill_manual(values = STATE_COLS, name = "Clone size") +
      labs(x = NULL, y = "Clonotypes",
           title = paste0("HIV-specific clonotypes (n = ", nrow(hiv_clono), ")"),
           subtitle = "Antigen-validated clonotypes across 88 months") +
      theme_cp018(base_size = 15)
    save_png(file.path(PERS_ROOT, "Alluvial_HIV_specific.png"), p_h,
             width = 9, height = 7)

    save_csv(hiv_clono %>% arrange(desc(n_total)),
             file.path(csv_dir, "CP018_HIVspecific_clonotype_fates.csv"))
    cat(sprintf("\nHIV-specific clonotypes: %d | persisting 2m->90m: %d\n",
                nrow(hiv_clono), sum(hiv_clono$fate == "Persistent (2m + 90m)")))
  }
}

################################################################################
# 4c. Persistence-focused views
################################################################################
# The whole-repertoire alluvial above is dominated by ~4,300 clonotypes flowing
# into "Not detected", so the handful that ARE maintained are sub-pixel ribbons.
# These panels make maintenance the subject:
#   (i)  every persistent clonotype as its own labelled ribbon
#   (ii) the full repertoire recoloured by FATE, so persistent flows are visible
#   (iii) cell-weighted flow, where a 15-cell clone is not equal to a singleton
#
# NOTE ON WHAT "MAINTAINED" MEANS HERE: 3 clonotypes is what exact-sequence
# matching gives. That count is robust - CTnt, CTaa and CTstrict all return the
# same 3 - so it is not an artefact of over-strict matching.

pers_ct <- clono %>% filter(fate == "Persistent (2m + 90m)") %>% arrange(desc(n_total))

if (nrow(pers_ct) > 0) {

  # ---- (i) each persistent clonotype as its own ribbon -----------------------
  # y = CELLS, so a 15-cell clone reads as 15x a singleton. Short CDR3b label
  # keeps the axis readable; full sequences are in the CSV.
  # Some cells carry two beta chains (";"-separated) - a real biological
  # occurrence (allelic inclusion) as well as a possible doublet signature.
  # Keep the first chain for the label and mark dual-chain clonotypes with
  # "(2\u03b2)" rather than truncating mid-sequence, which reads as data loss.
  short_lab <- function(aa) {
    b     <- sub("^[^_]*_", "", aa)
    dual  <- grepl(";", b)
    first <- sub(";.*$", "", b)
    ifelse(dual, paste0(first, " (2\u03b2)"), first)
  }
  pers_long <- pers_ct %>%
    mutate(clone_id = paste0(short_lab(CTaa), "  (n=", n_total, ")"),
           clone_id = factor(clone_id, levels = rev(clone_id))) %>%
    select(clone_id, `2m` = n_2m, `90m` = n_90m, hiv_specific) %>%
    pivot_longer(c(`2m`, `90m`), names_to = "Timepoint", values_to = "n_cells") %>%
    mutate(Timepoint = factor(Timepoint, levels = c("2m", "90m")))

  p_pers <- ggplot(pers_long,
                   aes(x = Timepoint, y = n_cells, stratum = clone_id,
                       alluvium = clone_id, fill = clone_id)) +
    geom_flow(alpha = 0.7, width = 0.3, colour = "black", linewidth = 0.25) +
    geom_stratum(width = 0.3, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = ifelse(n_cells > 0, n_cells, "")),
              stat = "stratum", size = 4.2, fontface = "bold") +
    scale_fill_brewer(palette = "Set2", name = "Clonotype (CDR3\u03b2)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Cells",
         title = paste0("Maintained clonotypes (n = ", nrow(pers_ct), ")"),
         subtitle = paste0("Detected at both 2m and 90m, ", 88, " months apart.\n",
                           "Ribbon width = cells. See detection-power caveat.")) +
    theme_cp018(base_size = 14) +
    theme(plot.title.position = "plot",
          plot.title    = element_text(hjust = 0, face = "bold"),
          plot.subtitle = element_text(hjust = 0))
  save_png(file.path(PERS_ROOT, "Alluvial_persistent_only.png"), p_pers,
           width = 10.5, height = 6.5)

  # ---- companion: paired size plot, easier to read than 3 ribbons ------------
  p_pair <- ggplot(pers_long, aes(x = Timepoint, y = n_cells, group = clone_id,
                                  colour = clone_id)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    geom_point(size = 3.4) +
    geom_text(data = subset(pers_long, Timepoint == "2m"),
              aes(label = n_cells), hjust = 1.6, size = 4, show.legend = FALSE) +
    geom_text(data = subset(pers_long, Timepoint == "90m"),
              aes(label = n_cells), hjust = -0.8, size = 4, show.legend = FALSE) +
    scale_colour_brewer(palette = "Set2", name = "Clonotype (CDR3\u03b2)") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Cells detected",
         title = "Size of each maintained clonotype at both timepoints") +
    theme_cp018(base_size = 14) +
  theme(plot.title.position = "plot")
  save_png(file.path(PERS_ROOT, "Persistent_clone_sizes.png"), p_pair,
           width = 8, height = 5.5)

  save_csv(pers_long, file.path(csv_dir, "CP018_persistent_alluvial_source.csv"))
}

# ---- (ii) full repertoire, coloured by FATE --------------------------------
# Same data as the main alluvial, but fill = fate. Persistent flows get a
# saturated colour against grey, so "how much is maintained" is readable even
# when the answer is "very little" - which is itself the finding.
FATE_COLS <- c("Persistent (2m + 90m)" = "#B2182B",
               "2m only"               = "#92C5DE",
               "90m only"              = "#F4A261")
allu_fate <- clono %>%
  transmute(fate,
            `2m`  = factor(size_cat(n_2m),  levels = LEVELS),
            `90m` = factor(size_cat(n_90m), levels = LEVELS)) %>%
  count(fate, `2m`, `90m`, name = "n_clonotypes") %>%
  mutate(alluvium = row_number()) %>%
  pivot_longer(c(`2m`, `90m`), names_to = "Timepoint", values_to = "state") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("2m", "90m")),
         state = factor(state, levels = LEVELS))

p_fate <- ggplot(allu_fate,
                 aes(x = Timepoint, stratum = state, alluvium = alluvium,
                     y = n_clonotypes, fill = fate)) +
  geom_flow(alpha = 0.75, width = 0.35, colour = NA) +
  geom_stratum(width = 0.35, colour = "black", linewidth = 0.3, fill = "grey95") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.4) +
  scale_fill_manual(values = FATE_COLS, name = "Fate") +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(x = NULL, y = "Clonotypes",
       title = "Clonotype fate across 88 months",  # persistent flows in red
       subtitle = paste0(sum(clono$fate == "Persistent (2m + 90m)"),
                         " of ", sum(clono$n_2m > 0),
                         " clonotypes seen at 2m were re-detected at 90m")) +
  theme_cp018(base_size = 14) + theme(plot.title.position = "plot")
save_png(file.path(PERS_ROOT, "Alluvial_by_fate.png"), p_fate, width = 9.5, height = 7)

# ---- (iii) cell-weighted: how much of the REPERTOIRE is maintained? ---------
# Weighting by cells rather than clonotypes asks a different question: what
# fraction of the sampled T-cell pool sits in a maintained clone?
cellw <- clono %>%
  group_by(fate) %>%
  summarise(cells_2m = sum(n_2m), cells_90m = sum(n_90m), .groups = "drop") %>%
  pivot_longer(c(cells_2m, cells_90m), names_to = "Timepoint", values_to = "n_cells") %>%
  mutate(Timepoint = factor(sub("cells_", "", Timepoint), levels = c("2m", "90m")))
save_csv(cellw, file.path(csv_dir, "CP018_cells_by_fate.csv"))

p_cellw <- ggplot(cellw, aes(x = Timepoint, y = n_cells, fill = fate)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.65, position = "fill") +
  scale_fill_manual(values = FATE_COLS, name = "Fate") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of cells with a TCR",
       title = "Share of sampled cells in a maintained clonotype") +
  theme_cp018(base_size = 14) +
  theme(plot.title.position = "plot")
save_png(file.path(PERS_ROOT, "Cells_by_fate.png"), p_cellw, width = 7, height = 5.5)

cat("\n=== Cells by fate ===\n")
print(as.data.frame(cellw), row.names = FALSE)

################################################################################
# 5. The persistent clonotypes, in full
################################################################################
pers <- clono %>% filter(fate == "Persistent (2m + 90m)") %>% arrange(desc(n_total))

# Report WHICH annotation each persistent clonotype occupies, and note when a
# clonotype spans more than one - a clonotype is a set of cells, and those cells
# need not share a phenotype.
if (has_anno && nrow(pers) > 0) {
  pers <- pers %>%
    mutate(spans_multiple_states = grepl(";", annotation),
           n_states = lengths(strsplit(annotation, ";")))
}
save_csv(pers, file.path(csv_dir, "CP018_persistent_clonotypes.csv"))

# Per-cell breakdown: which state is each persistent clone's cells in, and at
# which timepoint? The 2m and 90m cells of one clonotype can differ.
if (has_anno && nrow(pers) > 0) {
  pers_cells <- tcr %>%
    filter(CTstrict %in% pers$CTstrict) %>%
    count(CTstrict, CTaa, Timepoint, Fig_Annotation, name = "n_cells") %>%
    arrange(CTstrict, Timepoint)
  save_csv(pers_cells, file.path(csv_dir, "CP018_persistent_clonotype_states.csv"))
  cat("\n=== State of each persistent clonotype's cells ===\n")
  print(as.data.frame(pers_cells), row.names = FALSE)

  multi <- sum(pers$spans_multiple_states, na.rm = TRUE)
  cat("\n", nrow(pers), " persistent clonotype(s); ", multi,
      " span more than one annotated state.\n",
      "Do not describe them as sharing a single phenotype unless n_states == 1 ",
      "for all.\n", sep = "")
}
cat("\n=== Persistent clonotypes (detected at BOTH timepoints) ===\n")
if (nrow(pers)) {
  print(as.data.frame(pers[, c("CTaa","n_total","n_2m","n_90m",
                               "rep_confirmed_2m","rep_confirmed_90m",
                               "clusters","hiv_specific")]), row.names = FALSE)
} else cat("  none\n")

summary_tbl <- data.frame(
  metric = c("Total clonotypes", "Detected at 2m", "Detected at 90m",
             "Persistent (both)", "Persistence rate (% of 2m)",
             "Singleton clonotypes (%)",
             "Technical replicate overlap 2m (%)",
             "Technical replicate overlap 90m (%)",
             "Interval between timepoints (months)"),
  value = c(nrow(clono), sum(clono$n_2m > 0), sum(clono$n_90m > 0),
            nrow(pers),
            round(100 * nrow(pers) / sum(clono$n_2m > 0), 2),
            round(100 * mean(clono$n_total == 1), 1),
            ctrl$pct_of_smaller[1], ctrl$pct_of_smaller[2], 88))
save_csv(summary_tbl, file.path(csv_dir, "CP018_persistence_summary.csv"))
cat("\n=== Summary ===\n"); print(summary_tbl, row.names = FALSE)

cat("\nINTERPRETATION GUARD:\n",
    "Persistence (", round(100*nrow(pers)/sum(clono$n_2m>0), 2), "% of 2m clonotypes) must be read\n",
    "against the technical replicate ceiling (", ctrl$pct_of_smaller[1], "% / ",
    ctrl$pct_of_smaller[2], "%). Splitting ONE sample in two\n",
    "recovers only a small fraction of the same clonotypes, so non-detection at\n",
    "90m is the expected outcome for any clone that is not large. Report\n",
    "persistence of EXPANDED clones; do not claim clonal loss from these data.\n", sep = "")

cat("\nOutputs:\n  ", PERS_ROOT, "\n")
