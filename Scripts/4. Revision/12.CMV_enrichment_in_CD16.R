################################################################################
# Revision — are database-predicted CMV-specific clonotypes enriched in the
#            CD16+ effector population?        (Reviewer 2, major 3d and 5c)
#
# THE TWO COMMENTS THIS SERVES
#   3d: "The potential confounding effect of CMV co-infection is not even
#        mentioned."
#   5c: "A CD16+ effector CD8 population has been described in congenital CMV
#        infection (PMID 39531313)."
#
#   Script 11 showed our CD16+ population matches the congenital-CMV phenotype
#   on every testable feature. That raises the obvious question the reviewer is
#   really asking: is the population CMV-driven in OUR data? This tests the one
#   piece of CMV evidence the cohort contains.
#
# THE HONEST LIMIT, STATED BEFORE THE RESULT
#   CMV serostatus was NEVER MEASURED in this cohort. The only CMV signal
#   available is database annotation of TCR CDR3 sequences - which is exactly
#   the annotation Reviewer 2 (points 3a-3c) correctly calls unreliable, partly
#   because HIV-reactive CDR3s can be mis-assigned to CMV. So this test uses a
#   weak instrument, and it cannot establish CMV causation whichever way it
#   comes out. It is reported as supporting evidence for a Discussion
#   limitation, never as a result. An enrichment would be consistent with a CMV
#   contribution; an absence would not exclude one, because the annotation
#   misses private clonotypes entirely.
#
# OUT
#   Tables/R2M3d_02_CMV_enrichment_by_population.csv
#   Figures/R2M3d_01_CMV_enrichment.png
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(scales)
})
source("config_revision.R")

o <- load_qs("cd8_annotated_scored")
ANN <- if ("CD8_Annotation_Rev" %in% names(o@meta.data)) "CD8_Annotation_Rev" else "CD8_Annotation"
md <- o@meta.data
md$.pop <- as.character(md[[ANN]])

# ---- bring in the CMV annotation from the TCR object -----------------------
# The CD8 object does not carry epitope calls; they live on the TARA TCR
# objects as TRB_Epitope.species. There are three, one per Levenshtein edit
# distance (TRB_0/1/2). The manuscript reports a threshold of 2, so TRB_2 is
# the object that matches the published analysis - using a different threshold
# here would answer a question the paper did not ask.
tcr_f <- file.path(BASE, "saved_R_data", "TARA_TCR_Combined.RData")
stopifnot(file.exists(tcr_f))
cat("loading", basename(tcr_f), sprintf("(%.1f GB) ...\n", file.size(tcr_f)/1073741824))
e <- new.env(); load(tcr_f, envir = e)
tcr <- get("TARA_ALL_TRB_2", envir = e)
stopifnot("TRB_Epitope.species" %in% names(tcr@meta.data))

ann <- data.frame(.cell = colnames(tcr),
                  species = as.character(tcr$TRB_Epitope.species),
                  target  = as.character(tcr$TRB_Epitope.target),
                  stringsAsFactors = FALSE)
rm(e, tcr); invisible(gc())

md$.cell <- rownames(md)
n_before <- nrow(md)
md <- dplyr::left_join(md, ann, by = ".cell")
stopifnot(nrow(md) == n_before)
matched <- sum(!is.na(md$species) | md$.cell %in% ann$.cell)
cat("CD8 cells matched to a TCR-object row:", sum(md$.cell %in% ann$.cell),
    "of", nrow(md), "\n")
if (sum(md$.cell %in% ann$.cell) < 0.5 * nrow(md))
  stop("Fewer than half the CD8 barcodes matched the TCR object. Barcode ",
       "conventions differ between objects; resolve before reporting any ",
       "enrichment, since an unmatched cell would be counted as CMV-negative.")

# species is a semicolon-separated multi-assignment, e.g. "CMV;EBV"
md$.cmv <- !is.na(md$species) & grepl("CMV", md$species, ignore.case = TRUE)
cat("cells with a CMV-annotated clonotype:", sum(md$.cmv), "of", nrow(md),
    sprintf("(%.2f%%)\n", 100 * mean(md$.cmv)))

# Denominator must be cells that COULD be annotated, i.e. those with a
# recovered TCR. Including TCR-negative cells would dilute every population by
# its TCR recovery rate rather than its CMV content.
# A cell can only be called CMV if it has a recovered TRB that the database
# could match. Cells absent from the TCR object have no chance of a call, so
# including them would scale every population by its TCR recovery rate instead
# of its CMV content.
md$.hastcr <- md$.cell %in% ann$.cell
cat("denominator = CD8 cells present in the TCR object:", sum(md$.hastcr), "cells\n")
use <- md %>% filter(.hastcr)

tb <- use %>% group_by(.pop) %>%
  summarise(n_tcr_cells = n(), n_cmv = sum(.cmv),
            pct_cmv = round(100 * mean(.cmv), 2), .groups = "drop") %>%
  arrange(desc(pct_cmv))

# ---- test: CD16+ vs the rest ------------------------------------------------
TARGET <- grep("CD16", tb$.pop, value = TRUE)
stopifnot(length(TARGET) == 1)
a <- use %>% filter(.pop == TARGET); b <- use %>% filter(.pop != TARGET)
m <- matrix(c(sum(a$.cmv), nrow(a) - sum(a$.cmv),
              sum(b$.cmv), nrow(b) - sum(b$.cmv)), nrow = 2,
            dimnames = list(c("CMV", "notCMV"), c(TARGET, "other")))
cat("\n=== contingency ===\n"); print(m)
ft <- fisher.test(m)
cat(sprintf("\nFisher exact: OR = %.2f  (95%% CI %.2f-%.2f)  p = %.3g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

cat("\n=== CMV-annotated fraction by population (TCR-recovered cells) ===\n")
print(as.data.frame(tb), row.names = FALSE)

tb$is_target <- tb$.pop == TARGET
tb$fisher_OR_vs_rest <- ifelse(tb$is_target, round(ft$estimate, 3), NA)
tb$fisher_p_vs_rest  <- ifelse(tb$is_target, signif(ft$p.value, 3), NA)
tb$interpretation <- paste0(
  "Database-predicted CMV specificity only. CMV serostatus was NOT measured in ",
  "this cohort, and database annotation is itself unreliable (Reviewer 2, 3a-3c): ",
  "HIV-reactive CDR3s can be mis-assigned to CMV, and private clonotypes are ",
  "unannotatable. Supporting evidence for a Discussion limitation, not a result.")

mkdirs(REV_TAB); mkdirs(REV_FIG)
write.csv(tb, file.path(REV_TAB, "R2M3d_02_CMV_enrichment_by_population.csv"), row.names = FALSE)

pl <- tb %>% mutate(.pop = factor(.pop, levels = rev(.pop)))
p <- ggplot(pl, aes(x = pct_cmv, y = .pop, fill = is_target)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%  (%d/%d)", pct_cmv, n_cmv, n_tcr_cells)),
            hjust = -0.06, size = 4.4) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey65"), guide = "none") +
  scale_x_continuous(limits = c(0, max(pl$pct_cmv) * 1.45),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Database-predicted CMV-specific clonotypes by CD8 population",
       subtitle = paste0(TARGET, " highlighted; Fisher exact vs all other populations: OR = ",
                         sprintf("%.2f", ft$estimate), ", p = ", signif(ft$p.value, 3),
                         "\nDenominator is TCR-recovered cells. CMV serostatus was not ",
                         "measured; database annotation is a weak instrument."),
       x = "cells whose clonotype is CMV-annotated", y = NULL) +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 13.5, colour = "grey25"))

ggsave(file.path(REV_FIG, "R2M3d_01_CMV_enrichment.png"), plot = p,
       width = 14.5, height = 9, dpi = 300, bg = "white")
cat("\nwrote R2M3d_02_CMV_enrichment_by_population.csv and R2M3d_01_CMV_enrichment.png\n")

################################################################################
# A CONFOUND IN THIS ANALYSIS THAT MUST BE TESTED, NOT ASSUMED AWAY
#
# The CD16+ population is 37% gamma-delta TCR-positive (script 11), and the
# epitope databases are built almost entirely from alpha-beta TCRs. A
# gamma-delta cell therefore CANNOT receive a CMV call regardless of its true
# specificity, which would depress the CMV fraction of any gamma-delta-rich
# population for a purely technical reason. The observed DEPLETION is exactly
# what that artefact would produce, so the depletion cannot be interpreted
# until gamma-delta cells are excluded and the test repeated.
#
# Note the same applies to KIR+ innate-like CD8, the other low population.
################################################################################
cnt2 <- GetAssayData(o, assay = "RNA", layer = "counts")
gd_genes <- intersect(c("TRDC", "TRGC1", "TRGC2"), rownames(cnt2))
gd <- Matrix::colSums(cnt2[gd_genes, , drop = FALSE]) > 0
md$.gd <- gd[md$.cell]

cat("\n=== gamma-delta content by population ===\n")
gdtb <- md %>% filter(.hastcr) %>% group_by(.pop) %>%
  summarise(n = n(), pct_gd = round(100 * mean(.gd), 1), .groups = "drop") %>%
  arrange(desc(pct_gd))
print(as.data.frame(gdtb), row.names = FALSE)

cat("\n=== the same test, gamma-delta cells EXCLUDED ===\n")
ab <- md %>% filter(.hastcr, !.gd)
cat("alpha-beta-only cells:", nrow(ab), "of", sum(md$.hastcr), "\n")
a2 <- ab %>% filter(.pop == TARGET); b2 <- ab %>% filter(.pop != TARGET)
m2 <- matrix(c(sum(a2$.cmv), nrow(a2) - sum(a2$.cmv),
               sum(b2$.cmv), nrow(b2) - sum(b2$.cmv)), nrow = 2,
             dimnames = list(c("CMV","notCMV"), c(TARGET, "other")))
print(m2)
ft2 <- fisher.test(m2)
cat(sprintf("Fisher exact (alpha-beta only): OR = %.2f  (95%% CI %.2f-%.2f)  p = %.3g\n",
            ft2$estimate, ft2$conf.int[1], ft2$conf.int[2], ft2$p.value))
cat(sprintf("  CD16+ CMV fraction: %.2f%%   other populations: %.2f%%\n",
            100*mean(a2$.cmv), 100*mean(b2$.cmv)))

cat("\n=== verdict ===\n")
if (ft2$p.value > 0.05) {
  cat("The depletion does NOT survive excluding gamma-delta cells. It was a\n",
      "technical artefact of database coverage, not a biological result.\n",
      "Report: no CMV enrichment or depletion can be established.\n", sep="")
} else if (ft2$estimate < 1) {
  cat("The depletion PERSISTS among alpha-beta cells, so it is not explained by\n",
      "gamma-delta database blindness alone. Still not evidence against CMV:\n",
      "database annotation covers only public clonotypes.\n", sep="")
} else {
  cat("Direction REVERSES among alpha-beta cells - report the alpha-beta result\n",
      "as primary and the unrestricted one as confounded.\n", sep="")
}

gdtb$test <- "gamma-delta content"
write.csv(gdtb, file.path(REV_TAB, "R2M3d_03_gammadelta_by_population.csv"), row.names = FALSE)
write.csv(data.frame(
  analysis = c("all TCR-recovered cells", "alpha-beta cells only"),
  OR = c(round(ft$estimate,3), round(ft2$estimate,3)),
  ci_low = c(round(ft$conf.int[1],3), round(ft2$conf.int[1],3)),
  ci_high = c(round(ft$conf.int[2],3), round(ft2$conf.int[2],3)),
  p = c(signif(ft$p.value,3), signif(ft2$p.value,3)),
  pct_cmv_CD16pos = c(round(100*mean(a$.cmv),2), round(100*mean(a2$.cmv),2)),
  pct_cmv_other   = c(round(100*mean(b$.cmv),2), round(100*mean(b2$.cmv),2)),
  n_CD16pos = c(nrow(a), nrow(a2))),
  file.path(REV_TAB, "R2M3d_04_CMV_enrichment_gd_controlled.csv"), row.names = FALSE)
cat("\nwrote R2M3d_03_gammadelta_by_population.csv and R2M3d_04_CMV_enrichment_gd_controlled.csv\n")

