################################################################################
# Stamp the canonical clonotype column (CT_CANON) onto every object
#
# WHY THIS EXISTS
#   Cell Ranger 10.1.0 writes an EMPTY D-gene field where 9.0.1 and the older
#   eCITE-seq alignment write a literal "NA":
#
#     CP003 stim (CR 9.0.1)      TRBV6-1.NA.TRBJ2-5.TRBC2;TGTGCC...
#     unstim eCITE-seq (older)   TRBV15.NA.TRBJ2-3.TRBC2;TGTGCC...
#     CP018 stim (CR 10.1.0)     TRBV29-1..TRBJ1-2.TRBC1;TGCAGC...   <- ".."
#
#   Any cross-experiment comparison on raw CTstrict therefore returns ZERO
#   matches, which is indistinguishable from real clonal absence. Rather than
#   each script remembering to normalise, this stamps a stored column once so
#   the objects themselves are consistent.
#
# IS CP003 SAFE? YES - measured, not assumed.
#   Normalisation only ever rewrites ".." -> ".NA.", so objects already using
#   ".NA." are untouched. This script PROVES that per object and refuses to
#   save if a legacy object would change (section 3), so published CP003
#   numbers cannot drift.
#
# IDEMPOTENT: safe to re-run. Running it twice changes nothing.
#
# Usage:  Rscript "Scripts/CP018 - FACTseq/0.Normalise_Clonotypes.R"
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr)})

OUT <- file.path(CP018_ANALYSIS, "Clonotype_Normalisation"); mkdirs(OUT)

# label            path                          expected convention
TARGETS <- list(
  list(lab = "CP018 stim (FACT-Seq)",   path = QS_WITH_TCR,             legacy = FALSE),
  list(lab = "CP018 stim annotated",    path = QS_ANNOTATED,            legacy = FALSE),
  list(lab = "CP003 stim (published)",  path = QS_CP003_ANNOT_FIG3,     legacy = TRUE),
  list(lab = "unstim CP018 crossref",
       path = file.path(SAVED_DIR, "tara_cp018_crossref.qs2"),          legacy = TRUE),
  list(lab = "CP003+CP018 integrated",
       path = file.path(SAVED_DIR, "CP003_CP018_integrated.qs2"),       legacy = FALSE)
)

report <- list()

for (tg in TARGETS) {
  if (!file.exists(tg$path)) {
    cat("SKIP (not built yet):", tg$lab, "\n")
    next
  }
  cat("\n", tg$lab, "\n", sep = "")
  obj <- qs_read(tg$path)

  if (!"CTstrict" %in% colnames(obj@meta.data)) {
    cat("  no CTstrict column - skipping (TCR not attached)\n")
    rm(obj); gc(verbose = FALSE); next
  }

  raw <- as.character(obj@meta.data$CTstrict)
  nm  <- normalise_ct(raw)
  has <- !is.na(raw) & raw != "NA" & raw != ""
  n_changed <- sum(has & raw != nm, na.rm = TRUE)

  ##############################################################################
  # 3. GUARD: a legacy object must NOT change. If it does, the assumption
  #    behind "CP003 is unaffected" is wrong and nothing should be written.
  ##############################################################################
  if (tg$legacy && n_changed > 0) {
    stop("REFUSING TO WRITE: ", tg$lab, " is expected to use the legacy\n",
         "'.NA.' convention, but normalisation would change ", n_changed,
         " of its ", sum(has), " clonotype strings.\n",
         "That would alter published numbers. Investigate before proceeding.",
         call. = FALSE)
  }

  obj <- add_ct_canon(obj)

  # Distinct-clonotype count must be preserved or REDUCED (two raw strings can
  # collapse to one after normalisation, never the reverse). An increase means
  # the transform is corrupting strings.
  n_raw  <- n_distinct(raw[has])
  n_norm <- n_distinct(nm[has])
  if (n_norm > n_raw) {
    stop("Normalisation INCREASED distinct clonotypes (", n_raw, " -> ", n_norm,
         ") for ", tg$lab, ". That is impossible for a collapsing transform; ",
         "aborting.", call. = FALSE)
  }

  qs_save(obj, file = tg$path)

  report[[length(report) + 1]] <- data.frame(
    object            = tg$lab,
    file              = basename(tg$path),
    cells             = ncol(obj),
    cells_with_TCR    = sum(has),
    strings_normalised = n_changed,
    pct_normalised    = round(if (sum(has)) 100 * n_changed / sum(has) else 0, 1),
    clonotypes_raw    = n_raw,
    clonotypes_canon  = n_norm,
    collapsed         = n_raw - n_norm,
    convention        = obj$ct_convention[1],
    stringsAsFactors  = FALSE)

  cat("  distinct clonotypes:", n_raw, "->", n_norm,
      if (n_raw > n_norm) paste0(" (", n_raw - n_norm, " collapsed)") else " (unchanged)", "\n")
  cat("  saved:", basename(tg$path), "\n")
  rm(obj); gc(verbose = FALSE)
}

REP <- do.call(rbind, report)
save_csv(REP, file.path(OUT, "clonotype_normalisation_report.csv"))
cat("\n================ SUMMARY ================\n")
print(REP[, c("object", "cells_with_TCR", "strings_normalised", "pct_normalised",
              "clonotypes_raw", "clonotypes_canon")], row.names = FALSE)

# Check the LEGACY-flagged objects specifically. A substring match on "CP003"
# would also catch "CP003+CP018 integrated", which legitimately contains CP018
# cells and is expected to change - hence matching on the flag, not the name.
legacy_labs <- vapply(TARGETS[vapply(TARGETS, `[[`, logical(1), "legacy")],
                      `[[`, character(1), "lab")
legacy_rows <- REP[REP$object %in% legacy_labs, , drop = FALSE]
cat("\nLegacy objects (must be unchanged):\n")
if (nrow(legacy_rows)) {
  print(legacy_rows[, c("object", "strings_normalised")], row.names = FALSE)
  if (sum(legacy_rows$strings_normalised) == 0) {
    cat("  -> 0 changed, as expected. Published CP003 numbers unaffected.\n")
  } else {
    cat("  -> UNEXPECTED: a legacy object changed. Investigate before using.\n")
  }
} else cat("  (none present)\n")

cat("\nCR10-convention objects (expected to change):\n")
cr10 <- REP[!REP$object %in% legacy_labs, , drop = FALSE]
if (nrow(cr10)) print(cr10[, c("object", "strings_normalised", "pct_normalised")],
                      row.names = FALSE)
cat("\nEvery script can now use CT_CANON ('", CT_CANON,
    "') for clonotype identity.\n", sep = "")
