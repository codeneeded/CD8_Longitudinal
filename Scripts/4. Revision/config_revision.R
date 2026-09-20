################################################################################
# Revision — shared configuration
#
# source() this at the top of every script in Scripts/4. Revision/.
#
# WHY THIS EXISTS
#   The Figure_1..5 scripts hardcode base_dir <- "~/Documents/CD8_Longitudinal".
#   That path does not exist on every machine the repo has been cloned to, so
#   this resolves the repo root instead of assuming it. Set CD8_BASE to override.
#
# OUTPUT POLICY (house rules)
#   - PNG only. No PDF, no SVG.
#   - Every figure writes a *_source_data.csv beside it.
#   - Nothing is written to Manuscript/. Revision output goes to
#     Reviewer_Comments/Figures and Reviewer_Comments/Tables, and the user
#     decides later what is folded into the manuscript figures.
################################################################################

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(qs2)
})

# ── Repo root ────────────────────────────────────────────────────────────────
BASE <- Sys.getenv("CD8_BASE", "")
if (!nzchar(BASE)) {
  cands <- c("~/Documents/CD8_Longitudinal",
             "~/Desktop/Projects/Repositories/CD8_Longitudinal")
  for (c in cands) if (dir.exists(file.path(path.expand(c), "saved_R_data"))) {
    BASE <- path.expand(c); break
  }
}
if (!nzchar(BASE)) stop(
  "Could not find the repo. Set CD8_BASE, e.g.\n",
  '  Sys.setenv(CD8_BASE = "~/Documents/CD8_Longitudinal")', call. = FALSE)
BASE <- path.expand(BASE)

SAVED    <- file.path(BASE, "saved_R_data")
REV      <- file.path(BASE, "Reviewer_Comments")
REV_FIG  <- file.path(REV, "Figures")
REV_TAB  <- file.path(REV, "Tables")
MANUSCRIPT <- file.path(BASE, "Manuscript")   # READ ONLY - never write here

for (d in c(REV_FIG, REV_TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Objects ──────────────────────────────────────────────────────────────────
# Whole-PBMC WNN object: 126,435 cells. Annotation column = `Annotation`.
# Backs Figures 1 and 2. This is the clustering the reviewers commented on.
QS_ALL <- file.path(SAVED, "TARA_ALL_annotated_final.qs2")

# Dedicated CD8-only re-analysis: 18,843 cells, HEI only, 10 populations in
# `CD8_Annotation`, WNN clustering swept over wsnn_res.0.2 .. wsnn_res.2.
# Backs Figures 4, 5 and Supplementary 4-9.
# NOTE for Reviewer 2 major 1d: this object IS the dedicated sub-clustering the
# reviewer asks for. It exists and is already used - just not for Figure 1.
QS_CD8 <- file.path(SAVED, "TARA_cd8_HEI_annotated_final.qs2")

QS_MODULE <- file.path(SAVED, "TARA_ALL_ModuleScored.qs2")
QS_CMV    <- file.path(SAVED, "TARA_TRB_ED1_CMVgroup.qs2")

require_obj <- function(p) {
  if (!file.exists(p)) stop("Object not found:\n  ", p,
    "\nCheck SAVED in config_revision.R", call. = FALSE)
  invisible(p)
}

# ── Aesthetics: inherited verbatim from Figure_4.R so revision panels sit
#    beside the published ones without a legend mismatch ─────────────────────
ART_COLORS <- c(PreART_Entry = "#4A90D9", PostART_Suppressed = "#52B788",
                PostART_Unsuppressed = "#E76F51")
ART_LABELS <- c(PreART_Entry = "Pre-ART", PostART_Suppressed = "Suppressed",
                PostART_Unsuppressed = "Unsuppressed")

CD8_COLORS <- c(
  "Naïve CD8" = "#2D5F8A", "Naïve CD8 2" = "#1B9E77", "Naïve CD8 3" = "#7570B3",
  "Naïve CD8 4" = "#6BAED6", "Tscm CD8" = "#2CA02C",
  "Naïve Intermediate CD8" = "#98D8C8", "TEM CD8" = "#D95F02",
  "TEMRA CD8" = "#E6AB02", "CD16+ Effector CD8" = "#D62728",
  "KIR+ innate-like CD8" = "#8C564B")

# Editor asked for colour-vision-deficiency-safe palettes. Okabe-Ito.
LINEAGE_COLORS <- c("CD4 T" = "#0072B2", "CD8 T" = "#D55E00", "B" = "#009E73",
                    "NK" = "#CC79A7", "Mono" = "#E69F00", "DC" = "#56B4E9",
                    "other T" = "#999999", "other" = "#666666")

GROUP_COLORS <- c(HEI = "#D55E00", HEU = "#0072B2", HUU = "#009E73")

theme_rev <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(
      plot.title    = element_text(face = "bold", size = base + 2, hjust = 0),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.caption  = element_text(size = base - 2.5, colour = "grey45", hjust = 0),
      axis.title    = element_text(face = "bold", size = base),
      legend.title  = element_text(face = "bold", size = base - 1),
      strip.background = element_blank(),
      strip.text    = element_text(face = "bold", size = base))
}

# ── Save helpers ─────────────────────────────────────────────────────────────
# PNG ONLY, per instruction. Sizes are in INCHES so type scales with the canvas
# rather than being sized for one figure and reused on another.
save_fig <- function(name, plot, width = 8, height = 6, dpi = 400) {
  p <- file.path(REV_FIG, paste0(name, ".png"))
  ggsave(p, plot, width = width, height = height, dpi = dpi,
         bg = "white", device = "png", limitsize = FALSE)
  message("  [png] ", basename(p), "  (", width, "x", height, " in)")
  invisible(p)
}

# Every figure ships its numbers.
save_src <- function(name, df, dir = REV_FIG) {
  p <- file.path(dir, paste0(name, "_source_data.csv"))
  write.csv(df, p, row.names = FALSE)
  message("  [csv] ", basename(p), "  (", nrow(df), " rows)")
  invisible(p)
}

save_tab <- function(name, df) {
  p <- file.path(REV_TAB, paste0(name, ".csv"))
  write.csv(df, p, row.names = FALSE)
  message("  [csv] ", basename(p), "  (", nrow(df), " rows)")
  invisible(p)
}

# Reviewer 2, minor 2: label UMAPs with the number of cells.
n_cells_label <- function(obj) paste0("n = ", format(ncol(obj), big.mark = ","), " cells")

# ── Directory helper ─────────────────────────────────────────────────────────
mkdirs <- function(...) invisible(lapply(c(...), dir.create,
                                        recursive = TRUE, showWarnings = FALSE))

# ── Intermediate object checkpoints ──────────────────────────────────────────
# HOUSE RULE: every new UMAP, re-clustering or annotation must be saved as a
# .qs2 so a figure can be regenerated without recomputing. Objects go OUTSIDE
# the repo (they are gitignored data, and some are GB-scale) into saved_R_data/,
# prefixed rev_ so revision products are distinguishable from the published set.
#
#   save_qs(obj, "cd8_reclustered")  -> saved_R_data/rev_cd8_reclustered.qs2
#   obj <- load_qs("cd8_reclustered")
#
# UMAP coordinates live INSIDE the object, so re-plotting from a checkpoint is
# deterministic - the stochastic step is never repeated.
REV_QS <- SAVED

rev_qs_path <- function(name) file.path(REV_QS, paste0("rev_", name, ".qs2"))

save_qs <- function(obj, name, note = NULL) {
  p <- rev_qs_path(name)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(obj, file = p)
  cat(sprintf("  [qs2] %s  (%s cells, %.0f MB)%s\n", basename(p),
              format(ncol(obj), big.mark = ","),
              file.size(p) / 1048576,
              if (is.null(note)) "" else paste0("  - ", note)))
  # Manifest so the PDF and the reader can see what exists and when it was made.
  man <- file.path(REV_TAB, "rev_object_manifest.csv")
  row <- data.frame(object = basename(p), name = name,
                    cells = ncol(obj), size_MB = round(file.size(p) / 1048576, 1),
                    reductions = paste(Reductions(obj), collapse = ";"),
                    note = if (is.null(note)) "" else note,
                    written = format(Sys.time(), "%Y-%m-%d %H:%M"))
  if (file.exists(man)) {
    old <- read.csv(man, check.names = FALSE)
    old <- old[old$object != basename(p), , drop = FALSE]
    row <- rbind(old, row)
  }
  write.csv(row, man, row.names = FALSE)
  invisible(p)
}

load_qs <- function(name) {
  p <- rev_qs_path(name)
  require_obj(p)
  qs2::qs_read(p)
}

# TRUE when a checkpoint already exists - lets a script skip an expensive
# recompute on re-run while still regenerating every figure.
has_qs <- function(name) file.exists(rev_qs_path(name))

# Guard: the house rule is that nothing in the revision writes to Manuscript/.
assert_not_manuscript <- function(path) {
  if (grepl(normalizePath(MANUSCRIPT, mustWork = FALSE), normalizePath(path, mustWork = FALSE),
            fixed = TRUE))
    stop("Refusing to write into Manuscript/: ", path, call. = FALSE)
  invisible(TRUE)
}

message("Revision config loaded | BASE: ", BASE)
message("  figures -> ", REV_FIG)
message("  tables  -> ", REV_TAB)
