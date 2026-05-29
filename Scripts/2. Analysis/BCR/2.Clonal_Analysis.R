# ═══════════════════════════════════════════════════════════════════════════
# 02_bcr_clonal_analysis.R
#
# BCR-only clonal analysis (no Seurat). Follows the Immcantation 10x tutorial.
# Reads Stage 1 IgBLAST/MakeDb output, applies QC, calls clones per participant
# at the threshold determined in 01_threshold_inspection.R (0.10), reconstructs
# germlines, then computes the analyses this (minimally-expanded, polyclonal)
# repertoire actually supports: clonal abundance/diversity, SHM frequency, and
# isotype usage.
#
# CONTEXT — what this data is:
#   01_threshold_inspection.R found a UNIMODAL nearest-neighbour distribution
#   (median 0.31; only 1.5% of cells below 0.15) => minimal clonal expansion.
#   Clone-calling at 0.10 therefore yields mostly singletons, which is the
#   correct representation of a highly polyclonal infant blood B-cell
#   repertoire. Consequently:
#     - diversity / SHM / isotype = the real, supported readouts.
#     - lineage trees & expansion-based clonal tracking = not supported; the
#       script runs the checks and DOCUMENTS their absence rather than forcing
#       an empty/ misleading analysis.
#
# Seurat-dependent questions (SHM vs viral load, BCR phenotype on the UMAP,
# isotype/SHM by cell-type annotation) are handled in the NEXT script, which
# merges this output with the Seurat object by cell barcode.
#
# Clones are defined PER PARTICIPANT (across that participant's timepoints).
#
# PREREQUISITES:
#   install.packages(c("dplyr","tidyr","ggplot2","patchwork","alakazam","shazam","scoper","qs2"))
#   IMGT germlines copied out of Docker (once):
#     docker create --name tmp immcantation/suite:4.5.0
#     docker cp tmp:/usr/local/share/germlines/imgt ~/Documents/CD8_Longitudinal/VDJ/germlines
#     docker rm tmp
# ═══════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(alakazam)
library(shazam)
library(scoper)
library(qs2)

# ── Paths ──────────────────────────────────────────────────────────────────
# Read Stage 1 output from the OLD location; write all analysis to the NEW
# Analysis tree, organized into per-analysis subfolders.
data_path    <- "~/Documents/CD8_Longitudinal/VDJ/BCR"                   # Stage 1 (read)
germline_dir <- "~/Documents/CD8_Longitudinal/VDJ/germlines/imgt/human/vdj"
base_path    <- "~/Documents/CD8_Longitudinal/Analysis/VDJ/BCR"          # analysis home (write)
meta_file    <- file.path(base_path, "sample_metadata.csv")
save_path    <- file.path(base_path, "saved_R_data")

# Per-analysis figure subfolders
fig_qc        <- file.path(base_path, "qc")
fig_diversity <- file.path(base_path, "diversity")
fig_shm       <- file.path(base_path, "shm")
fig_isotype   <- file.path(base_path, "isotype")
for (d in c(save_path, fig_qc, fig_diversity, fig_shm, fig_isotype)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Parameters ───────────────────────────────────────────────────────────────
# Threshold determined in 01_threshold_inspection.R (unimodal -> fixed 0.10).
CLONE_THRESHOLD <- 0.10

# Helper: save a ggplot as both PNG (viewing) and PDF (record)
save_fig <- function(plot, dir, name, width, height, dpi = 150) {
  ggsave(file.path(dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
}


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 — LOAD STAGE 1 OUTPUT
# ───────────────────────────────────────────────────────────────────────────
# NOTE: read the full _db-pass.tsv (all loci) and do all QC here in R, so the
# data is identical to what the threshold step saw. cell_id_unique prefixes the
# barcode with sample_id (10x barcodes are reused across libraries).
# ═══════════════════════════════════════════════════════════════════════════

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
cat("Metadata:", nrow(meta), "samples\n")

load_dbpass <- function(sample_id) {
  f <- file.path(data_path, sample_id, paste0(sample_id, "_db-pass.tsv"))
  if (!file.exists(f)) { cat("  MISSING:", f, "\n"); return(NULL) }
  db <- read.delim(f, stringsAsFactors = FALSE)
  if (nrow(db) == 0) return(NULL)
  db$sample_id <- sample_id
  db
}

cat("Loading per-sample db-pass tables...\n")
bcr <- bind_rows(lapply(meta$sample, load_dbpass)) %>%
  left_join(meta, by = c("sample_id" = "sample"))
# 10x sequence_id and cell barcodes are only unique WITHIN a sample, so prefix
# both with sample_id before pooling. createGermlines requires unique
# sequence_id; clone-calling requires unique cell_id.
bcr$sequence_id    <- paste(bcr$sample_id, bcr$sequence_id, sep = "_")
bcr$cell_id_unique <- paste(bcr$sample_id, bcr$cell_id, sep = "_")
cat("Total sequences:", nrow(bcr),
    "|", length(unique(bcr$sample_id)), "samples,",
    length(unique(bcr$participant)), "participants\n")
cat("Unique sequence_id:", length(unique(bcr$sequence_id)) == nrow(bcr), "\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2 — QC FILTERING (tutorial-matched; identical to threshold step)
# ═══════════════════════════════════════════════════════════════════════════

qc_log <- data.frame(step = character(), before = integer(),
                     after = integer(), removed = integer())

logstep <- function(label, before, after) {
  qc_log <<- rbind(qc_log, data.frame(step = label, before = before,
                                      after = after, removed = before - after))
  cat(sprintf("  %-26s %6d -> %6d  (removed %d)\n", label, before, after, before - after))
}

# 2a V/J consistency (C-gene checked only when present)
before <- nrow(bcr)
bcr <- bcr %>% filter(
  (grepl("^IGHV", v_call) & grepl("^IGHJ", j_call) & (is.na(c_call) | grepl("^IGH", c_call))) |
    (grepl("^IGKV", v_call) & grepl("^IGKJ", j_call) & (is.na(c_call) | grepl("^IGK", c_call))) |
    (grepl("^IGLV", v_call) & grepl("^IGLJ", j_call) & (is.na(c_call) | grepl("^IGL", c_call)))
)
logstep("V/J consistency", before, nrow(bcr))

# 2b productive
before <- nrow(bcr); bcr <- bcr %>% filter(productive)
logstep("Productive only", before, nrow(bcr))

# 2c remove multi-heavy cells
multi_heavy <- table(filter(bcr, locus == "IGH")$cell_id_unique)
before <- nrow(bcr)
bcr <- bcr %>% filter(!cell_id_unique %in% names(multi_heavy)[multi_heavy > 1])
logstep("Remove multi-heavy cells", before, nrow(bcr))

# 2d remove heavy-less cells
heavy_cells <- filter(bcr, locus == "IGH")$cell_id_unique
light_cells <- filter(bcr, locus %in% c("IGK", "IGL"))$cell_id_unique
before <- nrow(bcr)
bcr <- bcr %>% filter(!cell_id_unique %in% setdiff(light_cells, heavy_cells))
logstep("Remove heavy-less cells", before, nrow(bcr))

write.csv(qc_log, file.path(fig_qc, "qc_filtering_log.csv"), row.names = FALSE)
cat("\nQC log saved to", file.path(fig_qc, "qc_filtering_log.csv"), "\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 — DEFINE CLONES (per participant, heavy chain, threshold = 0.10)
# ═══════════════════════════════════════════════════════════════════════════

cat("Calling clones per participant (threshold =", CLONE_THRESHOLD, ")...\n")
results <- hierarchicalClones(
  bcr,
  cell_id          = "cell_id_unique",
  threshold        = CLONE_THRESHOLD,
  only_heavy       = TRUE,
  split_light      = FALSE,
  summarize_clones = FALSE,
  fields           = "participant"
)
n_clones <- length(unique(na.omit(results$clone_id)))
cat("Sequences with clone_id:", sum(!is.na(results$clone_id)),
    "| total clones:", n_clones, "\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 — RECONSTRUCT GERMLINES (per participant; needed for SHM)
# ═══════════════════════════════════════════════════════════════════════════

cat("Reading IMGT references...\n")
references <- readIMGT(dir = germline_dir)
cat("Reconstructing germlines...\n")
results <- createGermlines(results, references, fields = "participant", nproc = 1)

results_heavy <- results %>% filter(locus == "IGH")
cat("Heavy chains carried forward:", nrow(results_heavy), "\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 — CLONE-SIZE / EXPANSION SUMMARY  (documents the polyclonality)
# ═══════════════════════════════════════════════════════════════════════════

clone_tab <- results_heavy %>%
  filter(!is.na(clone_id)) %>%
  count(participant, clone_id, name = "clone_size")

size_spectrum <- as.data.frame(table(clone_size = clone_tab$clone_size))
pct_singleton <- round(100 * mean(clone_tab$clone_size == 1), 2)
n_expanded5   <- sum(clone_tab$clone_size >= 5)

cat("Clone-size spectrum:\n"); print(size_spectrum)
cat(sprintf("\nSingletons: %.2f%% of clones | clones with >=5 seqs: %d\n\n",
            pct_singleton, n_expanded5))

write.csv(clone_tab,     file.path(fig_qc, "clone_sizes_per_clone.csv"), row.names = FALSE)
write.csv(size_spectrum, file.path(fig_qc, "clone_size_spectrum.csv"),  row.names = FALSE)

# Clone-size bar plot per sample (overview)
clone_sizes_sample <- countClones(results_heavy, groups = c("participant", "sample_id"))
p_sizes <- ggplot(clone_sizes_sample, aes(x = seq_count)) +
  geom_bar(fill = "#4472C4") +
  facet_wrap(~ sample_id, scales = "free_y") +
  labs(x = "Sequences per clone", y = "Number of clones",
       title = "Clone-size distribution (per sample)",
       subtitle = sprintf("%.1f%% of clones are singletons — minimal expansion", pct_singleton)) +
  theme_bw(base_size = 9)
save_fig(p_sizes, fig_qc, "clone_size_distribution", width = 12, height = 9)
cat("Saved clone-size distribution plot\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6 — CLONAL ABUNDANCE & DIVERSITY
# ───────────────────────────────────────────────────────────────────────────
# NOTE: For a polyclonal repertoire these quantify HOW diverse (near-maximal,
# singleton-dominated). That is itself the finding. Grouped by sample; Hill
# numbers via alphaDiversity with bootstrap CIs.
# ═══════════════════════════════════════════════════════════════════════════

cat("Estimating abundance and diversity...\n")

abund <- estimateAbundance(results_heavy, group = "sample_id", nboot = 100)
p_abund <- plot(abund, silent = TRUE) +
  ggtitle("BCR clonal abundance (rank-abundance)") + theme_bw(base_size = 11)
save_fig(p_abund, fig_diversity, "clonal_abundance", width = 9, height = 6)

div <- alphaDiversity(results_heavy, group = "sample_id", nboot = 100)
p_div <- plot(div, silent = TRUE) +
  ggtitle("BCR clonal diversity (Hill numbers)") + theme_bw(base_size = 11)
save_fig(p_div, fig_diversity, "clonal_diversity", width = 9, height = 6)

# Diversity by condition too (HEI vs HUU)
div_cond <- tryCatch(alphaDiversity(results_heavy, group = "condition", nboot = 100),
                     error = function(e) NULL)
if (!is.null(div_cond)) {
  p_div_cond <- plot(div_cond, silent = TRUE) +
    ggtitle("BCR clonal diversity by condition") + theme_bw(base_size = 11)
  save_fig(p_div_cond, fig_diversity, "clonal_diversity_by_condition", width = 7, height = 6)
}
cat("Saved abundance + diversity plots\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 — SOMATIC HYPERMUTATION (SHM) in the V gene
# ───────────────────────────────────────────────────────────────────────────
# NOTE: per-sequence mutation frequency vs reconstructed germline over the IMGT
# V region. Independent of clonality, so fully valid here. This mu_freq column
# is the key quantity the NEXT (Seurat-merge) script uses for SHM-vs-viral-load
# and SHM-by-cell-type.
# ═══════════════════════════════════════════════════════════════════════════

cat("Calculating SHM frequency (IMGT V region)...\n")
results_heavy <- observedMutations(
  results_heavy,
  sequenceColumn   = "sequence_alignment",
  germlineColumn   = "germline_alignment_d_mask",
  regionDefinition = IMGT_V,
  frequency        = TRUE,
  combine          = TRUE,
  nproc            = 1
)

# SHM by participant/timepoint
p_shm <- ggplot(filter(results_heavy, !is.na(mu_freq)),
                aes(x = factor(timepoint_months), y = mu_freq, fill = condition)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  facet_wrap(~ participant, scales = "free_x") +
  labs(x = "Timepoint (months)", y = "V-gene mutation frequency",
       title = "Somatic hypermutation by participant / timepoint") +
  theme_bw(base_size = 10)
save_fig(p_shm, fig_shm, "shm_by_participant_timepoint", width = 11, height = 8)

# SHM by condition (HEI vs HUU)
p_shm_cond <- ggplot(filter(results_heavy, !is.na(mu_freq)),
                     aes(x = condition, y = mu_freq, fill = condition)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.25) +
  labs(x = NULL, y = "V-gene mutation frequency",
       title = "Somatic hypermutation by condition") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
save_fig(p_shm_cond, fig_shm, "shm_by_condition", width = 6, height = 6)
cat("Saved SHM plots\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 8 — ISOTYPE USAGE (class switching)
# ═══════════════════════════════════════════════════════════════════════════

cat("Summarizing isotype usage...\n")
iso <- results_heavy %>%
  filter(!is.na(c_call), c_call != "") %>%
  mutate(isotype = sub("\\*.*", "", c_call)) %>%
  count(participant, sample_id, timepoint_months, condition, isotype) %>%
  group_by(sample_id) %>% mutate(prop = n / sum(n)) %>% ungroup()

p_iso <- ggplot(iso, aes(x = factor(timepoint_months), y = prop, fill = isotype)) +
  geom_col() +
  facet_wrap(~ participant, scales = "free_x") +
  labs(x = "Timepoint (months)", y = "Proportion", fill = "Isotype",
       title = "Isotype distribution over time") +
  theme_bw(base_size = 10)
save_fig(p_iso, fig_isotype, "isotype_over_time", width = 11, height = 8)

write.csv(iso, file.path(fig_isotype, "isotype_proportions.csv"), row.names = FALSE)
cat("Saved isotype plot + table\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 9 — LINEAGE TREES / EXPANSION TRACKING  (documented as not supported)
# ───────────────────────────────────────────────────────────────────────────
# NOTE: trees need clones with >=5 distinct sequences; expansion tracking needs
# clones shared across timepoints. Given the unimodal distance distribution and
# >99% singletons, neither is expected. We run the check and record the result
# rather than forcing an empty plot.
# ═══════════════════════════════════════════════════════════════════════════

cat("Checking feasibility of lineage trees / expansion tracking...\n")
expanded <- clone_tab %>% filter(clone_size >= 5)
multi_tp <- results_heavy %>%
  filter(!is.na(clone_id)) %>%
  distinct(participant, clone_id, timepoint_months) %>%
  count(participant, clone_id, name = "n_timepoints") %>%
  filter(n_timepoints >= 2)

tree_note <- sprintf(
  paste0("Lineage trees and expansion-based clonal tracking were not performed: ",
         "%d clones reached >=5 sequences and %d clones spanned >=2 timepoints, ",
         "consistent with the unimodal nearest-neighbour distribution and the ",
         "highly polyclonal, minimally-expanded character of the repertoire."),
  nrow(expanded), nrow(multi_tp))
cat(tree_note, "\n\n")
writeLines(tree_note, file.path(fig_qc, "lineage_feasibility_note.txt"))


# ═══════════════════════════════════════════════════════════════════════════
# STEP 10 — SAVE
# ═══════════════════════════════════════════════════════════════════════════

qs_save(results,       file.path(save_path, "bcr_results_full.qs2"))
qs_save(results_heavy, file.path(save_path, "bcr_results_heavy.qs2"))

cat("════════════════════════════════════════════════════════\n")
cat("DONE.\n")
cat("  Objects:", save_path, "\n")
cat("  QC:      ", fig_qc, "\n")
cat("  Diversity:", fig_diversity, "\n")
cat("  SHM:     ", fig_shm, "\n")
cat("  Isotype: ", fig_isotype, "\n")
cat("\nbcr_results_heavy.qs2 carries clone_id + mu_freq + isotype, ready for the\n")
cat("next script (Seurat merge: SHM-vs-viral-load, phenotype, UMAP).\n")
cat("════════════════════════════════════════════════════════\n")