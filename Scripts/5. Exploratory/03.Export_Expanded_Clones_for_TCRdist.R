################################################################################
# EXPLORATORY - export expanded clonotypes for tcrdist specificity scoring
#
# NOT FOR THE PAPER. The published specificity call in this project is Trex
# (Borcherding's R package, annotateDB()), which matches TRB CDR3 only at
# Levenshtein edit distance <= 2 against VDJdb, McPAS-TCR and IEDB, and is
# stored in the metadata column TRB_Epitope.species. This script exports the
# input for an INDEPENDENT scoring pass with tcrdist, so the authors can see
# how the two methods compare. It replaces nothing: TRB_Epitope.species is
# read here and carried through unmodified.
#
# NOTE: Trex is NOT TCRex. They are different tools and are not conflated
# anywhere in this analysis.
#
# WHERE THE TREX CALL COMES FROM
#   TARA_ALL_annotated_final.qs2 does NOT carry an epitope column - checked,
#   none of its 67 metadata columns matches /pitope|pecies|Trex/. The calls
#   live on the TCR objects in TARA_TCR_Combined.RData, one per Levenshtein
#   edit distance (TARA_ALL_TRB_0/1/2). The manuscript reports a threshold of
#   2, and Scripts/4. Revision/12.CMV_enrichment_in_CD16.R - the script behind
#   the CMV enrichment odds ratio in the revision - reads TARA_ALL_TRB_2 for
#   exactly that reason. This script therefore also reads TARA_ALL_TRB_2, so
#   the Trex side of the comparison is the one the paper actually reports.
#
#   CD8_Annotation is likewise absent from the ALL object; it is joined from
#   rev_cd8_annotated_scored.qs2 (the CD8 subset, 18,843 cells), so it is
#   populated only for CD8 cells and NA elsewhere. That is expected.
#
# SCOPE. 72,705 cells carry a TCR over ~68,000 clonotypes, but the specificity
# question is about clones that actually expanded. The unit of analysis is
# therefore the SAMPLE x CLONOTYPE pair (orig.ident x CTstrict) with >= 2
# cells, which is a few hundred rows rather than tens of thousands.
#
# PARSING. scRepertoire's columns have three traps, all handled below:
#
#   CTgene  TRAV38-2/DV8.TRAJ47.TRAC_TRBV7-9.NA.TRBJ2-5.TRBC2
#   CTaa    CAVSVLSGGGADGLTF_CASSLSMNTEAFF
#
#   1. Alpha and beta are separated by "_", segments within a chain by ".".
#      The BETA slot is V . D . J . C with D left as NA, so beta J is the
#      THIRD field while alpha J is the SECOND.
#   2. The object carries BOTH spellings of the dual TRAV/TRDV genes,
#      "TRAV14/DV4" and "TRAV14DV4". tcrdist's germline CDR reference keys on
#      the slashed form only - resolve_v_loops("TRAV14DV4") returns NULL while
#      "TRAV14/DV4" resolves - so the slashless form is normalised, not
#      discarded. It cost ~840 cells in the CoNGA export.
#   3. Multi-chain cells join a second CDR3 with ";" in CTaa, but CTgene
#      reports only the FIRST V/J pair for that chain, so the fields are NOT
#      parallel. The first element of every ";" list is taken so that V, J and
#      CDR3 stay internally consistent.
#
# Unlike the CoNGA export, cells are NOT required to have a complete pair.
# tcrdist scores paired where both chains are present and single-chain beta
# otherwise, so beta-only clonotypes are kept and flagged via `chain_mode`.
#
# Out: Projects/Analysis_Exploratory/TCRdist/expanded_clones.csv
################################################################################

suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr)})

REPO <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
OUT  <- "/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- Trex calls from the ED-2 object (cached; the RData is 10.3 GB) --------
trex_csv <- file.path(OUT, "trex_epitope_calls_TRB2.csv")
if (!file.exists(trex_csv)) {
  tcr_f <- file.path(REPO, "saved_R_data", "TARA_TCR_Combined.RData")
  stopifnot(file.exists(tcr_f))
  cat("loading", basename(tcr_f), sprintf("(%.1f GB) ...\n", file.size(tcr_f) / 1073741824))
  e <- new.env(); load(tcr_f, envir = e)
  tcr <- get("TARA_ALL_TRB_2", envir = e)
  stopifnot("TRB_Epitope.species" %in% names(tcr@meta.data))
  trex <- data.frame(barcode = colnames(tcr),
                     TRB_Epitope.species = as.character(tcr$TRB_Epitope.species),
                     TRB_Epitope.target  = as.character(tcr$TRB_Epitope.target),
                     stringsAsFactors = FALSE)
  rm(e, tcr); invisible(gc())
  write.csv(trex, trex_csv, row.names = FALSE)
} else {
  trex <- read.csv(trex_csv, stringsAsFactors = FALSE)
}
cat("Trex (TRB_2) rows:", nrow(trex), "\n")

seu <- qs_read(file.path(REPO, "saved_R_data", "TARA_ALL_annotated_final.qs2"))
m   <- seu@meta.data
m$barcode <- rownames(m)
rm(seu); invisible(gc())
cat("cells in object:", nrow(m), "\n")

n_before <- nrow(m)
m <- dplyr::left_join(m, trex, by = "barcode")
stopifnot(nrow(m) == n_before)
cat("cells matched to a Trex row:", sum(m$barcode %in% trex$barcode),
    "of", nrow(m), "\n")
if (sum(m$barcode %in% trex$barcode) < 0.5 * nrow(m))
  stop("Fewer than half the barcodes matched the TCR object; barcode ",
       "conventions differ between objects. Resolve before reporting anything.")

## ---- CD8_Annotation from the CD8 subset ------------------------------------
cd8 <- qs_read(file.path(REPO, "saved_R_data", "rev_cd8_annotated_scored.qs2"))
cd8_ann <- data.frame(barcode = colnames(cd8),
                      CD8_Annotation = as.character(cd8$CD8_Annotation),
                      stringsAsFactors = FALSE)
rm(cd8); invisible(gc())
m <- dplyr::left_join(m, cd8_ann, by = "barcode")
stopifnot(nrow(m) == n_before)
cat("cells with a CD8_Annotation:", sum(!is.na(m$CD8_Annotation)), "\n")

## ---- helpers ---------------------------------------------------------------
blank_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x %in% c("", "NA", "None", "NA.NA", "NA.NA.NA.NA")] <- NA_character_
  x
}
split2 <- function(x, i) vapply(strsplit(as.character(x), "_", fixed = TRUE),
                               function(p) if (length(p) >= i) p[[i]] else NA_character_, "")
seg    <- function(x, i) vapply(strsplit(as.character(x), ".", fixed = TRUE),
                               function(p) if (length(p) >= i) p[[i]] else NA_character_, "")
first_of <- function(x) vapply(strsplit(as.character(x), ";", fixed = TRUE),
                              function(p) if (length(p) >= 1) p[[1]] else NA_character_, "")
# Trap 2: normalise the slashless dual TRAV/TRDV spelling to the IMGT form
# tcrdist's reference table uses.
fix_dv <- function(g) sub("^(TRAV[0-9]+(?:-[0-9]+)?)(DV[0-9]+)$", "\\1/\\2", g, perl = TRUE)

## ---- keep cells carrying a receptor ----------------------------------------
m <- m[!is.na(m$CTstrict) & !is.na(m$CTgene) & !is.na(m$CTaa), ]
cat("cells with a TCR:", nrow(m), "\n")

## ---- parse chains (trap 1: beta J is the THIRD field) ----------------------
a_gene <- blank_to_na(split2(m$CTgene, 1))
b_gene <- blank_to_na(split2(m$CTgene, 2))
m$va_gene <- blank_to_na(seg(a_gene, 1)); m$ja_gene <- blank_to_na(seg(a_gene, 2))
m$vb_gene <- blank_to_na(seg(b_gene, 1)); m$jb_gene <- blank_to_na(seg(b_gene, 3))
m$cdr3a   <- blank_to_na(split2(m$CTaa, 1))
m$cdr3b   <- blank_to_na(split2(m$CTaa, 2))

## ---- trap 3: multi-chain cells, keep the primary chain ---------------------
n_dual_a <- sum(grepl(";", m$cdr3a, fixed = TRUE), na.rm = TRUE)
n_dual_b <- sum(grepl(";", m$cdr3b, fixed = TRUE), na.rm = TRUE)
for (cc in c("cdr3a", "cdr3b", "va_gene", "ja_gene", "vb_gene", "jb_gene"))
  m[[cc]] <- blank_to_na(first_of(m[[cc]]))
cat("cells with a second alpha:", n_dual_a, " second beta:", n_dual_b,
    "- primary chain kept\n")

## ---- trap 2: dual TRAV/TRDV spelling ---------------------------------------
n_slashless <- sum(grepl("^TRAV[0-9]+(-[0-9]+)?DV[0-9]+$", m$va_gene), na.rm = TRUE)
m$va_gene <- fix_dv(m$va_gene)
cat("alpha V genes normalised from the slashless DV spelling:", n_slashless, "\n")

## ---- collapse to sample x clonotype, keep the modal receptor ---------------
mode1 <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else
                       names(sort(table(x), decreasing = TRUE))[1] }

cl <- m %>%
  group_by(sample = orig.ident, CTstrict) %>%
  summarise(
    n_cells            = n(),
    cloneSize          = mode1(as.character(cloneSize)),
    va_gene            = mode1(va_gene), ja_gene = mode1(ja_gene), cdr3a = mode1(cdr3a),
    vb_gene            = mode1(vb_gene), jb_gene = mode1(jb_gene), cdr3b = mode1(cdr3b),
    TRB_Epitope.species = mode1(TRB_Epitope.species),
    TRB_Epitope.target  = mode1(TRB_Epitope.target),
    Annotation          = mode1(Annotation),
    CD8_Annotation      = mode1(CD8_Annotation),
    .groups = "drop")

cat("sample x clonotype pairs:", nrow(cl), "\n")

## ---- EXPANDED = >= 2 cells in a sample -------------------------------------
exp_cl <- cl %>% filter(n_cells >= 2) %>% arrange(desc(n_cells))
cat("EXPANDED (>= 2 cells in a sample):", nrow(exp_cl), "\n")

## ---- what tcrdist can score ------------------------------------------------
exp_cl <- exp_cl %>%
  mutate(chain_mode = case_when(
    !is.na(cdr3a) & !is.na(va_gene) & !is.na(cdr3b) & !is.na(vb_gene) ~ "paired",
    !is.na(cdr3b) & !is.na(vb_gene)                                   ~ "beta_only",
    TRUE                                                              ~ "unscorable"))
print(table(exp_cl$chain_mode, useNA = "ifany"))

## Trex call present? (same CMV pattern as Scripts/2. Analysis/5.CMV_vs_Expanding_Clones.R)
cmv_pattern <- "CMV|cytomegalo|herpesvirus 5|HHV-?5"
exp_cl <- exp_cl %>%
  mutate(trex_annotated = !is.na(TRB_Epitope.species) & TRB_Epitope.species != "",
         trex_cmv       = trex_annotated &
                          grepl(cmv_pattern, TRB_Epitope.species, ignore.case = TRUE))
cat("Trex-annotated:", sum(exp_cl$trex_annotated),
    "(", round(100 * mean(exp_cl$trex_annotated), 1), "% ) | unannotated:",
    sum(!exp_cl$trex_annotated),
    "(", round(100 * mean(!exp_cl$trex_annotated), 1), "% )\n")
cat("Trex CMV calls:", sum(exp_cl$trex_cmv), "\n")

write.csv(exp_cl, file.path(OUT, "expanded_clones.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "expanded_clones.csv"), "\n")
