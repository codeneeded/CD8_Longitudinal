################################################################################
# Revision — is our CD16+ Effector CD8 population the congenital-CMV
#            FcγRIII+ CD8 population?               (Reviewer 2, major 5c)
#
# THE COMMENT
#   "A CD16+ effector CD8 population has been described in congenital CMV
#    infection (PMID 39531313)."
#
# WHAT THE PAPER ACTUALLY REPORTS
#   Semmes et al., J Clin Invest 135(1), 2024 (doi:10.1172/JCI181342) sort
#   FcγRIII(CD16)+ CD8+ T cells from cord blood of neonates with and without
#   congenital HCMV, and report that the population:
#     - is TEMRA-dominant
#     - is CD57-high, PD-1-high, NKG2C-high vs FcγRIII-negative CD8
#     - is cytotoxic: PRF1, GNLY, NKG7, granzymes up
#     - upregulates T-bet (TBX21) and DOWNregulates BCL11B
#     - contains a γδ TCR minority (median 19%; αβ median 76%)
#   Signature recorded verbatim in Tables/R2M5c_03_paper_signature_PMID39531313.csv.
#
# WHAT THIS SCRIPT DOES
#   Tests OUR CD16+ Effector CD8 population against each of those features,
#   using the same contrast the paper used: CD16+ versus the other CD8
#   populations. A feature-by-feature table, so the answer is "matches on 9 of
#   12" with the mismatches named, rather than an impression.
#
# WHAT IT CANNOT DO, AND WHY THAT MATTERS FOR THE RESPONSE
#   CMV serostatus was never measured in this cohort, and none of these infants
#   has congenital CMV documented. So a phenotypic match CANNOT establish that
#   our population is CMV-driven - it establishes that the same differentiation
#   programme is present. The reviewer's point stands either way: the
#   literature attributes this phenotype to CMV, we cannot exclude CMV, and the
#   manuscript must say so. That is a Discussion limitation, not a result.
#
# OUT
#   Tables/R2M5c_02_CD16_signature_comparison.csv
#   Figures/R2M5c_01_CD16_signature_match.png
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(scales)
})
source("config_revision.R")

o <- load_qs("cd8_annotated_scored")
DefaultAssay(o) <- "RNA"
if (length(Layers(o[["RNA"]])) > 3) o <- JoinLayers(o, assay = "RNA")

ANN <- "CD8_Annotation_Rev"
if (!ANN %in% names(o@meta.data)) ANN <- "CD8_Annotation"
stopifnot(ANN %in% names(o@meta.data))

pops <- sort(unique(as.character(o@meta.data[[ANN]])))
TARGET <- grep("CD16", pops, value = TRUE)
stopifnot(length(TARGET) == 1)
cat("target population:", TARGET, "\n")
cat("contrast: all other CD8 populations (", length(pops) - 1, "populations )\n\n")

grp <- ifelse(as.character(o@meta.data[[ANN]]) == TARGET, "CD16pos", "CD16neg")
cnt <- GetAssayData(o, assay = "RNA", layer = "counts")

# pct of cells detecting a gene, per group — the paper reports frequencies of
# marker-positive cells, so a detection fraction is the comparable statistic.
pct <- function(g) {
  if (!g %in% rownames(cnt)) return(c(CD16pos = NA_real_, CD16neg = NA_real_))
  v <- cnt[g, ] > 0
  c(CD16pos = 100 * mean(v[grp == "CD16pos"]),
    CD16neg = 100 * mean(v[grp == "CD16neg"]))
}
pct_any <- function(gs) {
  gs <- intersect(gs, rownames(cnt))
  if (!length(gs)) return(c(CD16pos = NA_real_, CD16neg = NA_real_))
  v <- Matrix::colSums(cnt[gs, , drop = FALSE]) > 0
  c(CD16pos = 100 * mean(v[grp == "CD16pos"]),
    CD16neg = 100 * mean(v[grp == "CD16neg"]))
}

FEAT <- list(
  list(f = "FCGR3A (CD16)",      g = "FCGR3A",   dir = "up",   exp = "defining"),
  list(f = "CD57 (B3GAT1)",      g = "B3GAT1",   dir = "up",   exp = "higher"),
  list(f = "PD-1 (PDCD1)",       g = "PDCD1",    dir = "up",   exp = "higher"),
  list(f = "NKG2C (KLRC2)",      g = "KLRC2",    dir = "up",   exp = "higher"),
  list(f = "Perforin (PRF1)",    g = "PRF1",     dir = "up",   exp = "higher"),
  list(f = "Granulysin (GNLY)",  g = "GNLY",     dir = "up",   exp = "higher"),
  list(f = "NKG7",               g = "NKG7",     dir = "up",   exp = "higher"),
  list(f = "T-bet (TBX21)",      g = "TBX21",    dir = "up",   exp = "higher"),
  list(f = "BCL11B",             g = "BCL11B",   dir = "down", exp = "LOWER")
)

rows <- lapply(FEAT, function(x) {
  v <- pct(x$g)
  data.frame(feature = x$f, gene = x$g, paper_direction = x$dir,
             pct_CD16pos = round(v[["CD16pos"]], 1),
             pct_CD16neg = round(v[["CD16neg"]], 1), row.names = NULL)
})
tb <- bind_rows(rows)

# multi-gene features
gz <- pct_any(c("GZMB", "GZMH", "GZMA"))
tb <- bind_rows(tb, data.frame(
  feature = "Granzymes (any of B/H/A)", gene = "GZMB/GZMH/GZMA",
  paper_direction = "up",
  pct_CD16pos = round(gz[["CD16pos"]], 1), pct_CD16neg = round(gz[["CD16neg"]], 1)))

gd <- pct_any(c("TRDC", "TRGC1", "TRGC2"))
tb <- bind_rows(tb, data.frame(
  feature = "gamma-delta TCR (TRDC/TRGC1/2)", gene = "TRDC/TRGC1/TRGC2",
  paper_direction = "minority, paper median 19%",
  pct_CD16pos = round(gd[["CD16pos"]], 1), pct_CD16neg = round(gd[["CD16neg"]], 1)))

# ---- verdict per feature ----------------------------------------------------
tb <- tb %>% mutate(
  fold = ifelse(pct_CD16neg > 0, round(pct_CD16pos / pct_CD16neg, 2), NA_real_),
  matches_paper = case_when(
    is.na(pct_CD16pos) ~ "gene absent",
    grepl("minority", paper_direction) ~
      ifelse(pct_CD16pos < 50, "consistent", "NOT consistent"),
    paper_direction == "up"   ~ ifelse(pct_CD16pos > pct_CD16neg, "consistent", "NOT consistent"),
    paper_direction == "down" ~ ifelse(pct_CD16pos < pct_CD16neg, "consistent", "NOT consistent"),
    TRUE ~ "?"))

cat("=== feature-by-feature comparison against PMID 39531313 ===\n")
print(as.data.frame(tb), row.names = FALSE)

n_ok <- sum(tb$matches_paper == "consistent")
cat(sprintf("\nconsistent on %d of %d testable features\n",
            n_ok, sum(tb$matches_paper %in% c("consistent", "NOT consistent"))))
bad <- tb %>% filter(matches_paper == "NOT consistent")
if (nrow(bad)) {
  cat("\nfeatures that DO NOT match (these must be stated in the response):\n")
  print(as.data.frame(bad %>% select(feature, pct_CD16pos, pct_CD16neg)), row.names = FALSE)
}

# ---- differentiation state --------------------------------------------------
# The paper reports the population is TEMRA-dominant. Ours is annotated, so the
# comparable statement is which populations carry the CD16 signal.
cat("\n=== is our CD16+ population effector/TEMRA-like by annotation? ===\n")
cat("population name:", TARGET, "\n")
cat("cells:", sum(grp == "CD16pos"), "of", length(grp),
    sprintf("(%.1f%%)\n", 100 * mean(grp == "CD16pos")))

mkdirs(REV_TAB); mkdirs(REV_FIG)
write.csv(tb, file.path(REV_TAB, "R2M5c_02_CD16_signature_comparison.csv"), row.names = FALSE)

# ---- figure -----------------------------------------------------------------
# BUG FIXED: these counts must be taken BEFORE pivoting. Computed inside the
# mutate they counted rows of the pivoted frame - i.e. the number of FEATURES -
# so the legend read "n=11" for a population of 860 cells.
n_pos <- sum(grp == "CD16pos")
n_neg <- sum(grp == "CD16neg")
stopifnot(n_pos + n_neg == length(grp))
lab_pos <- paste0(TARGET, "\n(", format(n_pos, big.mark = ","), " cells)")
lab_neg <- paste0("all other CD8 populations\n(", format(n_neg, big.mark = ","), " cells)")

pl <- tb %>% filter(!is.na(pct_CD16pos)) %>%
  select(feature, matches_paper, CD16pos = pct_CD16pos, CD16neg = pct_CD16neg) %>%
  pivot_longer(c(CD16pos, CD16neg), names_to = "grp", values_to = "pct") %>%
  mutate(grp = factor(ifelse(grp == "CD16pos", lab_pos, lab_neg),
                      levels = c(lab_pos, lab_neg)),
         feature = factor(feature, levels = rev(tb$feature)))

p <- ggplot(pl, aes(x = pct, y = feature, fill = grp)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.72), hjust = -0.15, size = 4.6) +
  scale_fill_manual(values = setNames(c("#C0392B", "grey65"), c(lab_pos, lab_neg)),
                    name = NULL) +
  scale_x_continuous(limits = c(0, max(pl$pct, na.rm = TRUE) * 1.18),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Our CD16+ effector CD8 population against the congenital-CMV signature",
       subtitle = paste0("signature from Semmes et al., J Clin Invest 2024 (PMID 39531313)\n",
                         "consistent on ", n_ok, " of ",
                         sum(tb$matches_paper %in% c("consistent","NOT consistent")),
                         " testable features; CMV serostatus was not measured in this cohort"),
       x = "cells with detected expression", y = NULL) +
  theme_classic(base_size = 19) +
  theme(plot.title = element_text(face = "bold", size = 21),
        plot.subtitle = element_text(size = 15, colour = "grey25"),
        legend.position = "top")

ggsave(file.path(REV_FIG, "R2M5c_01_CD16_signature_match.png"), plot = p,
       width = 15, height = 10, dpi = 300, bg = "white")

cat("\nwrote R2M5c_02_CD16_signature_comparison.csv and R2M5c_01_CD16_signature_match.png\n")
cat("\nREMINDER: a phenotypic match does NOT establish CMV causation here.\n",
    "CMV serostatus was never measured in this cohort. The defensible claim is\n",
    "that the same differentiation programme is present and that CMV cannot be\n",
    "excluded as a driver - a Discussion limitation, not a result.\n", sep = "")
