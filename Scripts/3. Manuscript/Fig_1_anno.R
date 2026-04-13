################################################################################
# FIGURE 1 — PART 1: Annotation Validation (v2 — REVISED)
# Manuscript: CD8 Longitudinal (TARA Cohort)
#
# CHANGES FROM v1:
#   1. METADATA CORRECTION: Applies per-sample ground-truth Timepoint_Group,
#      Viral_Load_num, and Age from clinical Excel, plus orig.ident renaming
#      to PID_AgeM format. Applied to full TARA_ALL.
#
#   2. NO SUBSETTING: All annotation work is done on the full TARA_ALL object.
#      Subsetting to HEI/HEU/HUU is deferred to Fig1_Plots.R.
#
#   3. CLUSTER 1 RE-ANNOTATION: Cluster 1 ("Tcm CD8") is sub-gated by
#      ADT CD45RO and CD95 (FAS) protein expression into:
#        - Naïve CD8 (CD45RA+ CD45RO- CD95-)
#        - Tscm CD8  (CD45RA+ CD45RO- CD95+)
#        - Tcm CD8   (CD45RO+ CCR7+)
#
#   4. CLUSTER 12a RE-ANNOTATION: Renamed from "Tem CD8" to
#      "Transitional CD8" based on high CD62L/CCR7/CD27 retention
#      with activation markers (CD69, CD38, MKI67).
#
#   5. All downstream outputs regenerated with corrected annotations.
#
# OUTPUT: Fig 1/00_Annotation_Validation/
#
# USAGE: Run this script FIRST. Then run Fig1_Plots.R for manuscript panels
#        (which will subset to HEI/HEU/HUU as needed).
################################################################################

# ── Libraries ──────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(SeuratExtend)
library(scCustomize)
library(tidyr)
library(scRepertoire)
library(qs2)
library(Matrix)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(rstatix)
library(grid)
library(scales)

# ── Paths ──────────────────────────────────────────────────────────────────────
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
fig1_dir  <- file.path(base_dir, "Manuscript", "Fig 1")

# Non-manuscript supporting figures (annotation validation)
out_annot      <- file.path(fig1_dir, "00_Annotation_Validation")
out_vln_rna    <- file.path(out_annot, "Violins_RNA")
out_vln_adt    <- file.path(out_annot, "Violins_ADT")
out_feat_rna   <- file.path(out_annot, "FeaturePlots_RNA")
out_feat_adt   <- file.path(out_annot, "FeaturePlots_ADT")
out_avgexpr    <- file.path(out_annot, "AvgExpression")
out_gating     <- file.path(out_annot, "Cluster1_Gating")

# Manuscript figure panel directories (directly under Fig 1/)
out_umap  <- file.path(fig1_dir, "01_UMAP")
out_heat  <- file.path(fig1_dir, "02_ADT_RNA_Heatmap")
out_prop  <- file.path(fig1_dir, "03_Cluster_Proportions")
out_vl    <- file.path(fig1_dir, "04_Viral_Load_Correlations")

for (d in c(out_vln_rna, out_vln_adt, out_feat_rna, out_feat_adt, out_avgexpr,
            out_umap, out_heat, out_prop, out_vl, out_gating)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Read data ──────────────────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted.qs2"))


################################################################################
# STEP 1: METADATA CORRECTION
#
# Apply per-sample ground-truth metadata from clinical Excel
# (TARA_VL_CD4.xlsx). This corrects:
#   - CP020_V44: VL was 158 in Seurat, actual is 534,772
#   - CP013_24m: actual age is 25m
#   - CP018_V24: actual age is 25m
#   - CP003 V12/V24: VL values aligned to Excel
#
# Applied to the FULL TARA_ALL object.
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 1: METADATA CORRECTION (full TARA_ALL)\n",
    "══════════════════════════════════════════════════════════════\n\n")

# Define per-sample ground truth: Timepoint_Group, Viral_Load_num, Age
sample_metadata <- list(
  # CP002
  "CP002_entry"  = list(group = "PreART_Entry",         vl = 4284389,  age = 1),
  # CP003
  "CP003_entry"  = list(group = "PreART_Entry",         vl = 656769,   age = 1),
  "CP003_V12"    = list(group = "PostART_Unsuppressed",  vl = 503497,   age = 12),
  "CP003_V24"    = list(group = "PostART_Unsuppressed",  vl = 489676,   age = 24),
  # CP006
  "CP006_entry"  = list(group = "PreART_Entry",         vl = 10000000, age = 1),
  "CP006_12m"    = list(group = "PostART_Suppressed",    vl = 73,       age = 12),
  "CP006_V24"    = list(group = "PostART_Suppressed",    vl = 20,       age = 24),
  # CP011
  "CP011_entry"  = list(group = "PreART_Entry",         vl = 36965,    age = 2),
  # CP013
  "CP013_entry"  = list(group = "PreART_Entry",         vl = 3434,     age = 1),
  "CP013_12m"    = list(group = "PostART_Suppressed",    vl = 20,       age = 12),
  "CP013_24m"    = list(group = "PostART_Suppressed",    vl = 20,       age = 25),
  # CP016
  "CP016_entry"  = list(group = "PreART_Entry",         vl = 1978332,  age = 2),
  # CP017
  "CP017_entry"  = list(group = "PreART_Entry",         vl = 3167384,  age = 2),
  # CP018
  "CP018_entry"  = list(group = "PreART_Entry",         vl = 176970,   age = 1),
  "CP018_V24"    = list(group = "PostART_Suppressed",    vl = 20,       age = 25),
  "CP018_42m"    = list(group = "PostART_Suppressed",    vl = 20,       age = 42),
  # CP020
  "CP020_V1"     = list(group = "PreART_Entry",         vl = 414409,   age = 1),
  "CP020_V12"    = list(group = "PostART_Unsuppressed",  vl = 6293122,  age = 12),
  "CP020_V44"    = list(group = "PostART_Unsuppressed",  vl = 534772,   age = 44),
  # CP022
  "CP022_entry"  = list(group = "PreART_Entry",         vl = 5075764,  age = 1),
  # CP042
  "CP042_entry"  = list(group = "PreART_Entry",         vl = 6113,     age = 1)
)

# ── Print BEFORE state ──
cat("=== METADATA BEFORE CORRECTION ===\n")
cat("Timepoint_Group:\n")
print(table(TARA_ALL$Timepoint_Group))

samples_in_data <- unique(TARA_ALL$orig.ident)
cat("\nSamples in data:", paste(sort(samples_in_data), collapse = ", "), "\n")
cat("Samples in lookup:", paste(sort(names(sample_metadata)), collapse = ", "), "\n")

missing_in_lookup <- setdiff(samples_in_data, names(sample_metadata))
missing_in_data   <- setdiff(names(sample_metadata), samples_in_data)
if (length(missing_in_lookup) > 0) {
  cat("NOTE: Samples in data but NOT in lookup (HEU/HUU — no correction needed):\n")
  cat("  ", paste(missing_in_lookup, collapse = ", "), "\n")
}
if (length(missing_in_data) > 0) {
  cat("NOTE: Samples in lookup but NOT in data:",
      paste(missing_in_data, collapse = ", "), "\n")
}

# ── Apply corrections ──
n_group_changed <- 0
n_vl_changed    <- 0
n_age_changed   <- 0

for (sample_name in names(sample_metadata)) {
  mask <- TARA_ALL$orig.ident == sample_name
  n_cells <- sum(mask)
  if (n_cells == 0) next
  
  meta <- sample_metadata[[sample_name]]
  
  # Timepoint_Group
  old_group <- unique(TARA_ALL$Timepoint_Group[mask])
  if (length(old_group) != 1 || old_group != meta$group) {
    n_group_changed <- n_group_changed + n_cells
    TARA_ALL$Timepoint_Group[mask] <- meta$group
    if (length(old_group) == 1 && old_group != meta$group) {
      cat(sprintf("  CHANGED group %s: %s -> %s (%d cells)\n",
                  sample_name, old_group, meta$group, n_cells))
    }
  }
  
  # Viral_Load_num
  old_vl <- unique(TARA_ALL$Viral_Load_num[mask])
  if (length(old_vl) != 1 || abs(old_vl - meta$vl) > 1) {
    n_vl_changed <- n_vl_changed + n_cells
    TARA_ALL$Viral_Load_num[mask] <- meta$vl
    if (length(old_vl) == 1 && abs(old_vl - meta$vl) > 1) {
      cat(sprintf("  CHANGED VL   %s: %s -> %s (%d cells)\n",
                  sample_name, old_vl, meta$vl, n_cells))
    }
  }
  
  # Age
  old_age <- unique(TARA_ALL$Age[mask])
  if (length(old_age) != 1 || old_age != meta$age) {
    n_age_changed <- n_age_changed + n_cells
    TARA_ALL$Age[mask] <- meta$age
    if (length(old_age) == 1 && old_age != meta$age) {
      cat(sprintf("  CHANGED Age  %s: %s -> %s (%d cells)\n",
                  sample_name, old_age, meta$age, n_cells))
    }
  }
}

cat(sprintf("\nCells with Timepoint_Group changed: %d\n", n_group_changed))
cat(sprintf("Cells with Viral_Load_num changed:  %d\n", n_vl_changed))
cat(sprintf("Cells with Age changed:             %d\n", n_age_changed))

# ── Rename orig.ident to PID_AgeM format ──────────────────────────────────────
cat("\n=== Renaming orig.ident to PID_AgeM format ===\n")

ident_rename <- sapply(names(sample_metadata), function(s) {
  pid <- sub("_.*$", "", s)
  age <- sample_metadata[[s]]$age
  paste0(pid, "_", age, "m")
})
names(ident_rename) <- names(sample_metadata)

cat("Rename map:\n")
for (old_name in sort(names(ident_rename))) {
  cat(sprintf("  %s -> %s\n", old_name, ident_rename[old_name]))
}

# Store original ident
TARA_ALL$orig.ident_raw <- TARA_ALL$orig.ident

# Apply rename (only to HEI samples that are in the lookup)
new_idents <- unname(ident_rename[TARA_ALL$orig.ident])
n_renamed  <- sum(!is.na(new_idents))
n_kept     <- sum(is.na(new_idents))

# Keep original ident for HEU/HUU samples not in the rename map
new_idents[is.na(new_idents)] <- TARA_ALL$orig.ident[is.na(new_idents)]
TARA_ALL$orig.ident <- new_idents

cat(sprintf("Renamed %d HEI cells to PID_AgeM format\n", n_renamed))
cat(sprintf("Kept original ident for %d HEU/HUU cells\n", n_kept))

# ── Print AFTER state ──
cat("\n=== METADATA AFTER CORRECTION ===\n")
cat("Timepoint_Group:\n")
print(table(TARA_ALL$Timepoint_Group))

# ── Verification table ──
cat("\n=== PER-SAMPLE VERIFICATION (HEI only) ===\n")
verify <- TARA_ALL@meta.data %>%
  filter(Timepoint_Group %in% c("PreART_Entry", "PostART_Suppressed",
                                "PostART_Unsuppressed")) %>%
  group_by(orig.ident, Timepoint_Group) %>%
  summarise(
    n_cells = n(),
    VL      = first(Viral_Load_num),
    Age     = first(Age),
    .groups = "drop"
  ) %>%
  arrange(orig.ident)
print(as.data.frame(verify), row.names = FALSE)


################################################################################
# STEP 2: CLUSTER ANNOTATION REFINEMENT — INITIAL PASS
#
# Operates on full TARA_ALL (all conditions, all timepoints).
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 2: INITIAL CLUSTER ANNOTATION REFINEMENT (full object)\n",
    "══════════════════════════════════════════════════════════════\n\n")

# ── Diagnostic: check Cluster 12 ─────────────────────────────────────────────
cat("=== All unique annotations containing '12' ===\n")
all_annot <- as.character(TARA_ALL$Manual_Annotation)
cl12_raw <- unique(all_annot[grepl("12", all_annot)])
cat("  Manual_Annotation matches:",
    if (length(cl12_raw) > 0) paste(cl12_raw, collapse = "; ") else "NONE", "\n")
cat("  Cell count:", sum(grepl("12", all_annot)), "\n\n")

TARA_ALL$Manual_Annotation_refined <- case_when(
  # Cluster 9: αβ TCR+ contaminants → merge into Cluster 8
  TARA_ALL$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_ALL$has_TCR == TRUE
  ~ "8: TEMRA/CTL CD8",
  # Cluster 9: confirmed γδ
  TARA_ALL$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_ALL$has_TCR == FALSE
  ~ "9: γδ T cell",
  # Cluster 1 — temporary label pending sub-gating
  grepl("^1: ", as.character(TARA_ALL$Manual_Annotation))
  ~ "1: Naïve/Tcm CD8",
  # Cluster 6
  grepl("^6: ", as.character(TARA_ALL$Manual_Annotation))
  ~ "6: Naïve CD8",
  # Cluster 8
  grepl("^8: ", as.character(TARA_ALL$Manual_Annotation))
  ~ "8: TEMRA/CTL CD8",
  # Cluster 27
  grepl("^27: ", as.character(TARA_ALL$Manual_Annotation))
  ~ "27: Tex CD8",
  # Cluster 12: leave as-is for now (will split below)
  grepl("^12: ", as.character(TARA_ALL$Manual_Annotation))
  ~ "12: Mixed CD4/CD8",
  # All other clusters unchanged
  TRUE ~ as.character(TARA_ALL$Manual_Annotation)
)

# ── Split Cluster 12 by CD8A expression ───────────────────────────────────────
DefaultAssay(TARA_ALL) <- "RNA"
cd8a_expr <- GetAssayData(TARA_ALL, slot = "data")["CD8A", ]

is_cl12 <- TARA_ALL$Manual_Annotation_refined == "12: Mixed CD4/CD8"
cat("=== CD8A expression in Cluster 12 ===\n")
cat("  Cells in '12: Mixed CD4/CD8':", sum(is_cl12), "\n")
cat("  CD8A > 0:", sum(cd8a_expr[is_cl12] > 0, na.rm = TRUE), "\n")
cat("  CD8A == 0:", sum(cd8a_expr[is_cl12] == 0, na.rm = TRUE), "\n\n")

# 12a renamed to Transitional CD8 (not Tem)
TARA_ALL$Manual_Annotation_refined[is_cl12 & cd8a_expr > 0] <- "12a: Transitional CD8"
TARA_ALL$Manual_Annotation_refined[is_cl12 & (cd8a_expr <= 0 | is.na(cd8a_expr))] <- "12b: CD4"

cat("=== Cluster 12 split results ===\n")
n_12a <- sum(TARA_ALL$Manual_Annotation_refined == "12a: Transitional CD8")
n_12b <- sum(TARA_ALL$Manual_Annotation_refined == "12b: CD4")
cat("  12a (Transitional CD8):", n_12a, "cells\n")
cat("  12b (CD4):             ", n_12b, "cells\n\n")


################################################################################
# STEP 3: SUB-GATE CLUSTER 1 BY ADT PROTEIN
#
# Flow-based gating strategy on Cluster 1:
#   1. CD45RO+ → Tcm CD8
#   2. CD45RO- AND CD95/FAS+ → Tscm CD8
#   3. CD45RO- AND CD95- → Naïve CD8
#
# ADJUST THRESHOLDS after inspecting the diagnostic plots.
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 3: SUB-GATING CLUSTER 1 (Naïve/Tcm → Naïve + Tscm + Tcm)\n",
    "══════════════════════════════════════════════════════════════\n\n")

# ── Extract Cluster 1 cells ──────────────────────────────────────────────────
is_cl1 <- TARA_ALL$Manual_Annotation_refined == "1: Naïve/Tcm CD8"
cl1_cells <- colnames(TARA_ALL)[is_cl1]
cat("Cluster 1 cell count:", length(cl1_cells), "\n\n")

# ── Get ADT expression for gating markers ────────────────────────────────────
DefaultAssay(TARA_ALL) <- "ADT"
adt_data <- GetAssayData(TARA_ALL, slot = "data")

cd45ro_cl1 <- adt_data["CD45RO", cl1_cells]
cd45ra_cl1 <- adt_data["CD45RA", cl1_cells]
fas_cl1    <- adt_data["FAS",    cl1_cells]
cd62l_cl1  <- adt_data["SELL",   cl1_cells]

# ── Diagnostic: distribution of gating markers ───────────────────────────────
cat("=== CD45RO distribution in Cluster 1 ===\n")
cat("  Summary:\n"); print(summary(as.numeric(cd45ro_cl1)))
cat("  Quantiles:\n"); print(quantile(as.numeric(cd45ro_cl1),
                                      probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))

cat("\n=== FAS (CD95) distribution in Cluster 1 ===\n")
cat("  Summary:\n"); print(summary(as.numeric(fas_cl1)))
cat("  Quantiles:\n"); print(quantile(as.numeric(fas_cl1),
                                      probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))

cat("\n=== CD45RA distribution in Cluster 1 ===\n")
cat("  Summary:\n"); print(summary(as.numeric(cd45ra_cl1)))
cat("  Quantiles:\n"); print(quantile(as.numeric(cd45ra_cl1),
                                      probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
cat("\n")

# ── Thresholds ────────────────────────────────────────────────────────────────
# Default for DSB-normalized data. REVIEW THE DIAGNOSTIC PLOTS AND ADJUST.
cd45ro_threshold <- 0   # CD45RO > 0 = positive (memory)
fas_threshold    <- 0   # FAS > 0 = positive (antigen-experienced)

cat("=== Gating thresholds ===\n")
cat("  CD45RO threshold:", cd45ro_threshold, "\n")
cat("  FAS threshold:   ", fas_threshold, "\n\n")

# ── Apply gating logic ────────────────────────────────────────────────────────
cd45ro_pos <- as.numeric(cd45ro_cl1) > cd45ro_threshold
fas_pos    <- as.numeric(fas_cl1)    > fas_threshold

gate_assignment <- case_when(
  cd45ro_pos                  ~ "1c: Tcm CD8",
  !cd45ro_pos &  fas_pos      ~ "1b: Tscm CD8",
  !cd45ro_pos & !fas_pos      ~ "1a: Naïve CD8",
  TRUE                        ~ "1a: Naïve CD8"  # fallback
)
names(gate_assignment) <- cl1_cells

cat("=== Cluster 1 sub-gating results ===\n")
print(table(gate_assignment))
cat("\nProportions:\n")
print(round(prop.table(table(gate_assignment)) * 100, 1))
cat("\n")

# ── Apply to metadata ────────────────────────────────────────────────────────
TARA_ALL$Manual_Annotation_refined[match(names(gate_assignment),
                                         colnames(TARA_ALL))] <- gate_assignment

# ── Verify no cells were lost ─────────────────────────────────────────────────
n_cl1_after <- sum(grepl("^1[abc]:", TARA_ALL$Manual_Annotation_refined))
cat("Cluster 1 cells before gating:", length(cl1_cells), "\n")
cat("Cluster 1 cells after gating: ", n_cl1_after, "\n")
stopifnot("Cell count mismatch after Cluster 1 gating!" =
            n_cl1_after == length(cl1_cells))

# ── Generate gating diagnostic plots ──────────────────────────────────────────
message("Generating Cluster 1 gating diagnostic plots...")

gate_df <- data.frame(
  CD45RO = as.numeric(cd45ro_cl1),
  CD45RA = as.numeric(cd45ra_cl1),
  FAS    = as.numeric(fas_cl1),
  CD62L  = as.numeric(cd62l_cl1),
  Gate   = gate_assignment[cl1_cells],
  stringsAsFactors = FALSE
)

gate_cols <- c("1a: Naïve CD8" = "#4E79A7",
               "1b: Tscm CD8"  = "#59A14F",
               "1c: Tcm CD8"   = "#E15759")

# Plot 1: CD45RO vs CD45RA
p_gate1 <- ggplot(gate_df, aes(x = CD45RA, y = CD45RO, color = Gate)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = gate_cols) +
  geom_hline(yintercept = cd45ro_threshold, linetype = "dashed", color = "grey40") +
  labs(title = "Cluster 1: CD45RO vs CD45RA (ADT)",
       subtitle = paste("Dashed line = CD45RO threshold:", cd45ro_threshold),
       x = "CD45RA (DSB)", y = "CD45RO (DSB)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_gating, "Gate_CD45RO_vs_CD45RA.png"),
       p_gate1, width = 8, height = 7, dpi = 300, bg = "white")

# Plot 2: CD45RO vs FAS
p_gate2 <- ggplot(gate_df, aes(x = FAS, y = CD45RO, color = Gate)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = gate_cols) +
  geom_hline(yintercept = cd45ro_threshold, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = fas_threshold, linetype = "dashed", color = "grey40") +
  labs(title = "Cluster 1: CD45RO vs FAS/CD95 (ADT)",
       subtitle = paste("Thresholds — CD45RO:", cd45ro_threshold,
                        "| FAS:", fas_threshold),
       x = "FAS/CD95 (DSB)", y = "CD45RO (DSB)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_gating, "Gate_CD45RO_vs_FAS.png"),
       p_gate2, width = 8, height = 7, dpi = 300, bg = "white")

# Plot 3: FAS vs CD45RA (Tscm separation)
p_gate3 <- ggplot(gate_df, aes(x = CD45RA, y = FAS, color = Gate)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = gate_cols) +
  geom_hline(yintercept = fas_threshold, linetype = "dashed", color = "grey40") +
  labs(title = "Cluster 1: FAS/CD95 vs CD45RA (ADT)",
       subtitle = "Tscm = CD45RA+ CD95+ CD45RO-",
       x = "CD45RA (DSB)", y = "FAS/CD95 (DSB)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_gating, "Gate_FAS_vs_CD45RA.png"),
       p_gate3, width = 8, height = 7, dpi = 300, bg = "white")

# Plot 4: Histograms
p_hist1 <- ggplot(gate_df, aes(x = CD45RO, fill = Gate)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = gate_cols) +
  geom_vline(xintercept = cd45ro_threshold, linetype = "dashed") +
  labs(title = "CD45RO distribution in Cluster 1", x = "CD45RO (DSB)") +
  theme_minimal(base_size = 12)

p_hist2 <- ggplot(gate_df, aes(x = FAS, fill = Gate)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = gate_cols) +
  geom_vline(xintercept = fas_threshold, linetype = "dashed") +
  labs(title = "FAS/CD95 distribution in Cluster 1", x = "FAS (DSB)") +
  theme_minimal(base_size = 12)

p_hists <- p_hist1 / p_hist2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_gating, "Gate_Histograms.png"),
       p_hists, width = 10, height = 10, dpi = 300, bg = "white")

message("✓ Gating diagnostic plots saved to: ", out_gating)


################################################################################
# STEP 4: FINALIZE CLUSTER DEFINITIONS & EXTRACT CD8 SUBSET
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 4: FINALIZE CD8 CLUSTER DEFINITIONS\n",
    "══════════════════════════════════════════════════════════════\n\n")

TARA_ALL$Manual_Annotation_refined <- factor(TARA_ALL$Manual_Annotation_refined)

cd8_cluster_names <- c(
  "1a: Naïve CD8",
  "1b: Tscm CD8",
  "1c: Tcm CD8",
  "6: Naïve CD8",
  "8: TEMRA/CTL CD8",
  "9: γδ T cell",
  "12a: Transitional CD8",
  "27: Tex CD8"
)

cd8_short_labels <- setNames(
  c("Naïve CD8 (1a)", "Tscm CD8", "Tcm CD8", "Naïve CD8 (6)",
    "TEMRA/CTL", "γδ T cell", "Transitional CD8", "Tex CD8"),
  cd8_cluster_names
)

# ── Subset CD8 clusters from the FULL object ─────────────────────────────────
TARA_cd8 <- subset(TARA_ALL,
                   subset = Manual_Annotation_refined %in% cd8_cluster_names)

TARA_cd8$Manual_Annotation_refined <- droplevels(TARA_cd8$Manual_Annotation_refined)
TARA_cd8$Manual_Annotation_refined <- factor(
  TARA_cd8$Manual_Annotation_refined,
  levels = cd8_cluster_names[cd8_cluster_names %in%
                               levels(TARA_cd8$Manual_Annotation_refined)]
)

cat("=== CD8 subset cell counts (REVISED, all conditions) ===\n")
print(table(TARA_cd8$Manual_Annotation_refined))
cat("\n")

cat("=== CD8 clusters × Timepoint_Group ===\n")
print(table(TARA_cd8$Manual_Annotation_refined, TARA_cd8$Timepoint_Group))
cat("\n")

missing_clusters <- setdiff(cd8_cluster_names,
                            levels(TARA_cd8$Manual_Annotation_refined))
if (length(missing_clusters) > 0) {
  cat("WARNING: Missing clusters:", paste(missing_clusters, collapse = ", "), "\n")
  cat("  Check gating thresholds in Step 3.\n\n")
}

# ── Shared color palettes ────────────────────────────────────────────────────
cond_cols <- c("HUU" = "#4E79A7", "HEU" = "#F28E2B", "HEI" = "#E15759")

cluster_cols <- c(
  "1a: Naïve CD8"        = "#4E79A7",
  "1b: Tscm CD8"         = "#59A14F",
  "1c: Tcm CD8"          = "#E15759",
  "6: Naïve CD8"         = "#76B7B2",
  "8: TEMRA/CTL CD8"     = "#EDC948",
  "9: γδ T cell"         = "#B07AA1",
  "12a: Transitional CD8"= "#FF9DA7",
  "27: Tex CD8"          = "#9C755F"
)


################################################################################
# STEP 5: GENE / PROTEIN MARKER LISTS
################################################################################

DefaultAssay(TARA_cd8) <- "ADT"
all_adt_markers <- sort(rownames(TARA_cd8))
cat("=== Total ADT markers in assay:", length(all_adt_markers), "===\n\n")

rna_markers <- list(
  `CD8 Identity` = c(
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"
  ),
  `Naive - Quiescence` = c(
    "CCR7", "SELL", "TCF7", "LEF1", "BACH2", "IL7R", "BCL2",
    "S1PR1", "KLF2", "FOXP1", "SATB1"
  ),
  `Effector - Cytotoxicity` = c(
    "GZMB", "GZMA", "GZMH", "GZMK", "GZMM",
    "PRF1", "GNLY", "NKG7", "KLRG1", "CX3CR1",
    "FGFBP2", "FCGR3A"
  ),
  `Effector TFs` = c(
    "TBX21", "EOMES", "RUNX3", "PRDM1", "ZEB2", "ID2"
  ),
  `Exhaustion` = c(
    "TOX", "TOX2", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
    "CTLA4", "ENTPD1", "CD160", "CD244"
  ),
  `Memory - Stemness` = c(
    "IL7R", "CD27", "CD28", "BACH2", "TCF7", "BCL2",
    "BCL6", "CXCR5", "CCR7"
  ),
  `Activation - Proliferation` = c(
    "MKI67", "HLA-DRA", "HLA-DRB1", "CD38", "FAS",
    "ICOS", "TNFRSF9", "CXCR6", "CD69"
  ),
  `Tissue Residency` = c(
    "ITGAE", "ITGA1", "CD69", "CXCR6", "ZNF683", "PRDM1"
  ),
  `GammaDelta TCR` = c(
    "TRDV1", "TRDV2", "TRGV9", "TRDC", "TRGC1", "TRGC2"
  ),
  `MAIT` = c(
    "TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "RORC"
  ),
  `NK-like - Innate` = c(
    "TYROBP", "KLRD1", "KIR2DL3", "KIR3DL1", "NCAM1",
    "KLRF1", "KLRC1", "KLRC2"
  ),
  `Type I IFN` = c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44L", "MX1",
    "OAS1", "STAT1", "IRF7"
  ),
  `Stress - Heat Shock` = c(
    "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "IFI27"
  ),
  `Chemokines - Cytokines` = c(
    "CCL3", "CCL4", "CCL4L2", "CCL5", "XCL1", "XCL2",
    "IFNG", "TNF", "IL2", "CSF2"
  ),
  `Terminal Differentiation` = c(
    "B3GAT1", "KLRG1", "CX3CR1", "ZEB2", "GZMB"
  ),
  `CD4 T cell` = c(
    "CD4", "IL7R", "FOXP3", "IL2RA", "CTLA4", "IKZF2",
    "CCR4", "RORC", "GATA3", "TBX21", "BCL6",
    "CXCR5", "MAF", "ICOS", "BATF"
  ),
  `Treg` = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT"),
  `B cell` = c(
    "CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1",
    "BLK", "BLNK", "TCL1A", "IGHM", "IGHD", "IGKC", "IGLC2"
  ),
  `Plasma cell` = c("SDC1", "XBP1", "MZB1", "JCHAIN", "IGHG1", "IGHA1", "PRDM1", "IRF4"),
  `NK cell` = c(
    "NCAM1", "KLRF1", "KLRC1", "KLRB1", "NCR1", "NCR3",
    "FCGR3A", "TYROBP", "FCER1G", "SH2D1B", "SPON2", "CLIC3", "MYOM2"
  ),
  `Monocyte - Macrophage` = c(
    "CD14", "LYZ", "S100A8", "S100A9", "S100A12",
    "FCGR3A", "CSF1R", "CD68", "MARCO", "VCAN", "FCN1", "MNDA"
  ),
  `cDC` = c("FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1", "BATF3", "IRF8", "IRF4", "LAMP3"),
  `pDC` = c("LILRA4", "CLEC4C", "IRF7", "TCF4", "GZMB", "IL3RA", "JCHAIN"),
  `Platelet` = c("PPBP", "PF4", "GP9", "ITGA2B", "TUBB1", "TREML1"),
  `Erythrocyte` = c("HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA"),
  `HSC - Progenitor` = c("CD34", "SPINK2", "AVP", "CRHBP", "CYTL1", "SOX4", "HOPX"),
  `ILC` = c("IL7R", "KIT", "GATA3", "RORC", "IL1RL1", "AREG", "IL22", "NCR2")
)

rna_markers_flat <- unique(unlist(rna_markers))
cat("Total unique RNA markers:", length(rna_markers_flat), "\n")
cat("Total ADT markers:", length(all_adt_markers), "\n\n")


################################################################################
# STEP 6: CHECK WHICH MARKERS EXIST
################################################################################

DefaultAssay(TARA_cd8) <- "RNA"
rna_available <- rna_markers_flat[rna_markers_flat %in% rownames(TARA_cd8)]
rna_missing   <- setdiff(rna_markers_flat, rna_available)

cat("=== RNA: ", length(rna_available), " available, ",
    length(rna_missing), " missing ===\n")
if (length(rna_missing) > 0) cat("Missing:", paste(rna_missing, collapse = ", "), "\n\n")

DefaultAssay(TARA_cd8) <- "ADT"
adt_available <- all_adt_markers
cat("=== ADT: ", length(adt_available), " (all in assay) ===\n\n")


################################################################################
# STEP 7: EXPORT AVERAGE EXPRESSION TABLES
################################################################################

message("Exporting average expression tables...")

DefaultAssay(TARA_cd8) <- "RNA"

avg_rna_raw <- AverageExpression(
  TARA_cd8, assays = "RNA", features = rna_available,
  group.by = "Manual_Annotation_refined", slot = "data"
)$RNA

avg_rna_full <- AverageExpression(
  TARA_cd8, assays = "RNA",
  group.by = "Manual_Annotation_refined", slot = "data"
)$RNA

scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x)) / rng
  }))
  colnames(scaled) <- colnames(mat)
  scaled
}

avg_rna_scaled <- scale_01(as.matrix(avg_rna_raw))

colnames(avg_rna_raw)    <- gsub("^g ", "", colnames(avg_rna_raw))
colnames(avg_rna_scaled) <- gsub("^g ", "", colnames(avg_rna_scaled))
colnames(avg_rna_full)   <- gsub("^g ", "", colnames(avg_rna_full))

write.csv(avg_rna_raw,    file.path(out_avgexpr, "AvgExpr_RNA_per_cluster.csv"))
write.csv(round(avg_rna_scaled, 4), file.path(out_avgexpr, "AvgExpr_RNA_per_cluster_scaled.csv"))
write.csv(avg_rna_full,   file.path(out_avgexpr, "AvgExpr_RNA_AllGenes_per_cluster.csv"))

DefaultAssay(TARA_cd8) <- "ADT"

avg_adt_raw <- AverageExpression(
  TARA_cd8, assays = "ADT", features = adt_available,
  group.by = "Manual_Annotation_refined", slot = "data"
)$ADT

avg_adt_scaled <- scale_01(as.matrix(log10(avg_adt_raw + 1)))

colnames(avg_adt_raw)    <- gsub("^g ", "", colnames(avg_adt_raw))
colnames(avg_adt_scaled) <- gsub("^g ", "", colnames(avg_adt_scaled))

write.csv(avg_adt_raw,    file.path(out_avgexpr, "AvgExpr_ADT_per_cluster.csv"))
write.csv(round(avg_adt_scaled, 4), file.path(out_avgexpr, "AvgExpr_ADT_per_cluster_scaled.csv"))

# Percent expression
DefaultAssay(TARA_cd8) <- "RNA"
pct_rna <- as.data.frame(
  t(sapply(levels(droplevels(TARA_cd8$Manual_Annotation_refined)), function(cl) {
    cells <- WhichCells(TARA_cd8, expression = Manual_Annotation_refined == cl)
    mat   <- GetAssayData(TARA_cd8, slot = "data")[rna_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_rna, 2), file.path(out_avgexpr, "PctExpr_RNA_per_cluster.csv"))

DefaultAssay(TARA_cd8) <- "ADT"
pct_adt <- as.data.frame(
  t(sapply(levels(droplevels(TARA_cd8$Manual_Annotation_refined)), function(cl) {
    cells <- WhichCells(TARA_cd8, expression = Manual_Annotation_refined == cl)
    mat   <- GetAssayData(TARA_cd8, slot = "data")[adt_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_adt, 2), file.path(out_avgexpr, "PctExpr_ADT_per_cluster.csv"))

message("✓ Average expression tables saved")


################################################################################
# STEP 8: VIOLIN PLOTS — ADT
################################################################################

message("Generating ADT violin plots (", length(adt_available), " markers)...")

for (feat in adt_available) {
  p <- VlnPlot2(
    TARA_cd8, features = feat, group.by = "Manual_Annotation_refined",
    assay = "ADT", cols = "light", show.mean = TRUE, mean_colors = c("red", "blue")
  )
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(file.path(out_vln_adt, paste0("VlnADT_", safe_name, ".png")),
         p, width = 14, height = 5, dpi = 300, bg = "white")
}
message("✓ ADT violin plots saved\n")


################################################################################
# STEP 9: VIOLIN PLOTS — RNA
################################################################################

message("Generating RNA violin plots (", length(rna_available), " markers)...")

for (gene in rna_available) {
  if (!gene %in% rownames(TARA_cd8[["RNA"]])) next
  p <- VlnPlot2(
    TARA_cd8, features = gene, group.by = "Manual_Annotation_refined",
    assay = "RNA", cols = "light", show.mean = TRUE, mean_colors = c("red", "blue")
  )
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(file.path(out_vln_rna, paste0("VlnRNA_", safe_name, ".png")),
         p, width = 14, height = 5, dpi = 300, bg = "white")
}
message("✓ RNA violin plots saved\n")


################################################################################
# STEP 10: FEATURE PLOTS — ADT
################################################################################

message("Generating ADT feature plots...")
DefaultAssay(TARA_cd8) <- "ADT"
for (feat in adt_available) {
  if (!feat %in% rownames(TARA_cd8[["ADT"]])) next
  fp <- DimPlot2(TARA_cd8, features = feat, reduction = "wnn.umap") +
    ggtitle(paste("ADT |", feat))
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(file.path(out_feat_adt, paste0("FeatADT_", safe_name, ".png")),
         fp, dpi = 500, width = 8, height = 6, bg = "white")
}
message("✓ ADT feature plots saved\n")


################################################################################
# STEP 11: FEATURE PLOTS — RNA
################################################################################

message("Generating RNA feature plots...")
DefaultAssay(TARA_cd8) <- "RNA"
for (gene in rna_available) {
  if (!gene %in% rownames(TARA_cd8[["RNA"]])) next
  fp <- DimPlot2(TARA_cd8, features = gene, reduction = "wnn.umap") +
    ggtitle(paste("RNA |", gene))
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(file.path(out_feat_rna, paste0("FeatRNA_", safe_name, ".png")),
         fp, dpi = 500, width = 8, height = 6, bg = "white")
}
message("✓ RNA feature plots saved\n")


################################################################################
# STEP 12: COMBINED OVERVIEW — DotPlots
################################################################################

message("Generating DotPlots...")

rna_dotplot_genes <- intersect(
  c("CD8A", "CD8B", "CD3D", "CD3E",
    "CCR7", "SELL", "TCF7", "LEF1", "BACH2", "S1PR1", "KLF2",
    "IL7R", "BCL2",
    "GZMB", "GZMA", "GZMK", "GZMH", "PRF1", "GNLY", "NKG7", "FGFBP2",
    "TBX21", "EOMES", "RUNX3", "ZEB2",
    "TOX", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
    "B3GAT1", "CX3CR1", "KLRG1",
    "TRDV1", "TRGV9", "TRDC",
    "SLC4A10", "KLRB1", "ZBTB16",
    "TYROBP", "KLRD1", "KIR2DL3", "KIR3DL1",
    "MKI67", "HLA-DRA", "CD69", "FAS",
    "ISG15", "IFIT1", "IFI44L",
    "CCL3", "CCL4", "CCL5", "IFNG",
    "CD4", "FOXP3", "CD19", "MS4A1", "CD79A",
    "CD14", "LYZ", "S100A8", "NCAM1", "KLRF1",
    "PPBP", "PF4", "HBB"),
  rna_available
)

DefaultAssay(TARA_cd8) <- "RNA"
p_dot_rna <- DotPlot(
  TARA_cd8, features = rna_dotplot_genes,
  group.by = "Manual_Annotation_refined",
  cols = c("lightgrey", "#D73027"), dot.scale = 6
) + RotatedAxis() +
  ggtitle("RNA: CD8 Cluster Identity (REVISED)") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 9))

ggsave(file.path(out_annot, "DotPlot_RNA_CD8_Clusters_REVISED.png"),
       p_dot_rna, width = 26, height = 9, dpi = 300, bg = "white")

DefaultAssay(TARA_cd8) <- "ADT"
p_dot_adt <- DotPlot(
  TARA_cd8, features = adt_available,
  group.by = "Manual_Annotation_refined",
  cols = c("lightgrey", "#08519C"), dot.scale = 5
) + RotatedAxis() +
  ggtitle("ADT: All Surface Protein Markers (REVISED)") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 9))

ggsave(file.path(out_annot, "DotPlot_ADT_All_CD8_Clusters_REVISED.png"),
       p_dot_adt, width = 38, height = 9, dpi = 300, bg = "white")

message("✓ DotPlots saved\n")


################################################################################
# STEP 13: CLUSTER 12a (Transitional CD8) — ANNOTATION DIAGNOSTIC
################################################################################

message("Generating Cluster 12a diagnostics...")
DefaultAssay(TARA_cd8) <- "RNA"

TARA_cd8$is_cl12a <- ifelse(
  TARA_cd8$Manual_Annotation_refined == "12a: Transitional CD8",
  "Cluster 12a", "Other CD8"
)

cl12a_markers <- FindMarkers(
  TARA_cd8, ident.1 = "Cluster 12a", group.by = "is_cl12a",
  test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25
)

write.csv(cl12a_markers, file.path(out_annot, "Cluster12a_vs_OtherCD8_DEgenes_REVISED.csv"))

cat("=== Cluster 12a (Transitional CD8) — Top 30 UP ===\n")
print(head(rownames(cl12a_markers[cl12a_markers$avg_log2FC > 0, ]), 30))
cat("\n=== Top 30 DOWN ===\n")
print(head(rownames(cl12a_markers[cl12a_markers$avg_log2FC < 0, ]), 30))
cat("\n")

top6 <- head(rownames(cl12a_markers[cl12a_markers$avg_log2FC > 0, ]), 6)
if (length(top6) > 0) {
  cl12a_feat_plots <- lapply(top6, function(g) {
    DimPlot2(TARA_cd8, features = g, reduction = "wnn.umap") + ggtitle(paste("RNA |", g))
  })
  ggsave(file.path(out_annot, "Cluster12a_Top6_FeaturePlot_REVISED.png"),
         wrap_plots(cl12a_feat_plots, ncol = 3),
         width = 20, height = 12, dpi = 300, bg = "white")
}

TARA_cd8$is_cl12a <- NULL
message("✓ Cluster 12a diagnostics saved\n")


################################################################################
# STEP 14: VALIDATE NEW SUB-GATES
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 14: VALIDATE CLUSTER 1 SUB-GATES\n",
    "══════════════════════════════════════════════════════════════\n\n")

DefaultAssay(TARA_cd8) <- "RNA"

run_pairwise_de <- function(obj, ident1, ident2, outdir, label) {
  sub_obj <- subset(obj, subset = Manual_Annotation_refined %in% c(ident1, ident2))
  n1 <- sum(sub_obj$Manual_Annotation_refined == ident1)
  n2 <- sum(sub_obj$Manual_Annotation_refined == ident2)
  cat(sprintf("%s: %d vs %d cells\n", label, n1, n2))
  if (n1 < 10 | n2 < 10) { cat("  ⚠ Insufficient cells — skipping\n\n"); return(NULL) }
  de <- tryCatch({
    FindMarkers(sub_obj, ident.1 = ident1, group.by = "Manual_Annotation_refined",
                test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.15)
  }, error = function(e) { message("  ⚠ DE failed: ", e$message); NULL })
  if (!is.null(de)) {
    safe_label <- gsub("[^A-Za-z0-9._-]", "_", label)
    write.csv(de, file.path(outdir, paste0("DE_", safe_label, ".csv")))
    cat(sprintf("  Sig genes (p_adj < 0.05): %d\n", sum(de$p_val_adj < 0.05, na.rm = TRUE)))
    cat("  Top 10 UP:\n"); print(head(rownames(de[de$avg_log2FC > 0, ]), 10))
    cat("  Top 10 DOWN:\n"); print(head(rownames(de[de$avg_log2FC < 0, ]), 10))
    cat("\n")
  }
  return(de)
}

de_1b_vs_1a <- run_pairwise_de(TARA_cd8, "1b: Tscm CD8", "1a: Naïve CD8",
                               out_gating, "Tscm_vs_Naive_1b_vs_1a")

de_1c_vs_1a <- run_pairwise_de(TARA_cd8, "1c: Tcm CD8", "1a: Naïve CD8",
                               out_gating, "Tcm_vs_Naive_1c_vs_1a")

de_1b_vs_1c <- run_pairwise_de(TARA_cd8, "1b: Tscm CD8", "1c: Tcm CD8",
                               out_gating, "Tscm_vs_Tcm_1b_vs_1c")

cat("=== Comparing sub-gated Naïve (1a) vs original Naïve (6) ===\n")
de_1a_vs_6 <- run_pairwise_de(TARA_cd8, "1a: Naïve CD8", "6: Naïve CD8",
                              out_gating, "Naive1a_vs_Naive6")

if (!is.null(de_1a_vs_6)) {
  n_sig <- sum(de_1a_vs_6$p_val_adj < 0.05, na.rm = TRUE)
  if (n_sig < 20) {
    cat("  → Few DE genes: 1a and 6 are transcriptionally similar.\n")
    cat("  → Consider merging into a single Naïve CD8 cluster.\n\n")
  } else {
    cat("  → Substantial DE: 1a and 6 are distinct naïve sub-populations.\n")
    cat("  → Keep separate or investigate what distinguishes them.\n\n")
  }
}


################################################################################
# STEP 15: SAVE REVISED OBJECTS
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 15: SAVING REVISED OBJECTS\n",
    "══════════════════════════════════════════════════════════════\n\n")

qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_sorted_v2.qs2"))
cat("✓ Saved: TARA_ALL_sorted_v2.qs2\n")

qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_revised_v2.qs2"))
cat("✓ Saved: TARA_cd8_revised_v2.qs2\n")



################################################################################
# STEPS 16–20 REPLACEMENT (v4 — FINAL NAMING)
#
# Run this AFTER Step 15 has completed (v2 objects saved).
# Assumes TARA_ALL and TARA_cd8 are in memory with v2 annotations:
#   1a: Naïve CD8, 1b: Tscm CD8, 1c: Tcm CD8, 6: Naïve CD8,
#   8: TEMRA/CTL CD8, 9: γδ T cell, 12a: Transitional CD8, 27: Tex CD8
#
# CHANGES FROM v2/v3:
#   • CD45RO threshold stays at 0 (consistent with all other DSB gates)
#   • No reclassification of 1c cells
#   • RELABELLING ONLY:
#       1a: Naïve CD8       → 1a: Naïve 1 CD8
#       1c: Tcm CD8         → 1c: Naïve Intermediate CD8
#       6:  Naïve CD8       → 6:  Naïve 2 CD8
#   • 1b: Tscm CD8, 8, 9, 12a, 27 — unchanged
#   • Expression tables re-exported with final labels
################################################################################

TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted_v2.qs2"))

scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x)) / rng
  }))
  colnames(scaled) <- colnames(mat)
  scaled
}

################################################################################
# STEP 16: RELABEL CLUSTERS TO FINAL NOMENCLATURE
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 16: RELABEL CLUSTERS (v4 — final naming)\n",
    "══════════════════════════════════════════════════════════════\n\n")

# ── Extract the actual strings from the object (avoids encoding issues) ───────
TARA_ALL$Manual_Annotation_refined <- as.character(TARA_ALL$Manual_Annotation_refined)

old_1a <- unique(TARA_ALL$Manual_Annotation_refined)[
  grepl("^1a:", unique(TARA_ALL$Manual_Annotation_refined))]
old_1c <- unique(TARA_ALL$Manual_Annotation_refined)[
  grepl("^1c:", unique(TARA_ALL$Manual_Annotation_refined))]
old_6  <- unique(TARA_ALL$Manual_Annotation_refined)[
  grepl("^6:", unique(TARA_ALL$Manual_Annotation_refined))]

cat("Renaming:\n")
cat("  ", old_1a, " → 1a: Naive 1 CD8\n")
cat("  ", old_1c, " → 1c: Naive Intermediate CD8\n")
cat("  ", old_6,  " → 6: Naive 2 CD8\n\n")

# ── Apply renaming ────────────────────────────────────────────────────────────
TARA_ALL$Manual_Annotation_refined[TARA_ALL$Manual_Annotation_refined == old_1a] <- "1a: Naive 1 CD8"
TARA_ALL$Manual_Annotation_refined[TARA_ALL$Manual_Annotation_refined == old_1c] <- "1c: Naive Intermediate CD8"
TARA_ALL$Manual_Annotation_refined[TARA_ALL$Manual_Annotation_refined == old_6]  <- "6: Naive 2 CD8"

TARA_ALL$Manual_Annotation_refined <- factor(TARA_ALL$Manual_Annotation_refined)

# ── Updated cluster definitions ───────────────────────────────────────────────
cd8_cluster_names <- c(
  "1a: Naive 1 CD8",
  "1b: Tscm CD8",
  "1c: Naive Intermediate CD8",
  "6: Naive 2 CD8",
  "8: TEMRA/CTL CD8",
  "9: γδ T cell",
  "12a: Transitional CD8",
  "27: Tex CD8"
)

cluster_cols <- c(
  "1a: Naive 1 CD8"             = "#4E79A7",
  "1b: Tscm CD8"                = "#59A14F",
  "1c: Naive Intermediate CD8"  = "#E15759",
  "6: Naive 2 CD8"              = "#76B7B2",
  "8: TEMRA/CTL CD8"            = "#EDC948",
  "9: γδ T cell"                = "#B07AA1",
  "12a: Transitional CD8"       = "#FF9DA7",
  "27: Tex CD8"                 = "#9C755F"
)

# ── Verify ────────────────────────────────────────────────────────────────────
cat("=== TARA_ALL CD8 cluster counts (final naming) ===\n")
cd8_counts <- table(TARA_ALL$Manual_Annotation_refined)[
  names(table(TARA_ALL$Manual_Annotation_refined)) %in% cd8_cluster_names
]
print(cd8_counts)
cat("\n")


################################################################################
# STEP 17: REBUILD TARA_cd8 WITH FINAL LABELS
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 17: REBUILD TARA_cd8 WITH FINAL LABELS\n",
    "══════════════════════════════════════════════════════════════\n\n")

TARA_cd8 <- subset(TARA_ALL,
                   subset = Manual_Annotation_refined %in% cd8_cluster_names)

TARA_cd8$Manual_Annotation_refined <- droplevels(
  factor(TARA_cd8$Manual_Annotation_refined, levels = cd8_cluster_names)
)
DefaultAssay(TARA_cd8) <- "RNA"
rna_available <- rownames(TARA_cd8)  # or reload your curated list
DefaultAssay(TARA_cd8) <- "ADT"
adt_available <- rownames(TARA_cd8)
cat("=== TARA_cd8 cluster counts ===\n")
print(table(TARA_cd8$Manual_Annotation_refined))
cat("\n")

cat("=== CD8 clusters × Timepoint_Group ===\n")
print(table(TARA_cd8$Manual_Annotation_refined, TARA_cd8$Timepoint_Group))
cat("\n")

# ── Naive Intermediate (1c) breakdown by condition ────────────────────────────
cat("=== Naive Intermediate (1c) × Timepoint_Group ===\n")
is_1c <- TARA_cd8$Manual_Annotation_refined == "1c: Naive Intermediate CD8"
print(table(TARA_cd8$Timepoint_Group[is_1c]))
cat("\n")
print(round(prop.table(table(TARA_cd8$Timepoint_Group[is_1c])) * 100, 1))
cat("\n")


################################################################################
# STEP 18: RE-EXPORT EXPRESSION TABLES (v4 — final labels)
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 18: RE-EXPORT EXPRESSION TABLES (v4)\n",
    "══════════════════════════════════════════════════════════════\n\n")

out_avgexpr_v4 <- file.path(out_annot, "AvgExpression_v4")
dir.create(out_avgexpr_v4, recursive = TRUE, showWarnings = FALSE)

# ── RNA ───────────────────────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "RNA"

avg_rna_raw_v4 <- AverageExpression(
  TARA_cd8, assays = "RNA", features = rna_available,
  group.by = "Manual_Annotation_refined", slot = "data"
)$RNA

avg_rna_scaled_v4 <- scale_01(as.matrix(avg_rna_raw_v4))

colnames(avg_rna_raw_v4)    <- gsub("^g ", "", colnames(avg_rna_raw_v4))
colnames(avg_rna_scaled_v4) <- gsub("^g ", "", colnames(avg_rna_scaled_v4))

write.csv(avg_rna_raw_v4,
          file.path(out_avgexpr_v4, "AvgExpr_RNA_per_cluster.csv"))
write.csv(round(avg_rna_scaled_v4, 4),
          file.path(out_avgexpr_v4, "AvgExpr_RNA_per_cluster_scaled.csv"))

# Full transcriptome average
avg_rna_full_v4 <- AverageExpression(
  TARA_cd8, assays = "RNA",
  group.by = "Manual_Annotation_refined", slot = "data"
)$RNA
colnames(avg_rna_full_v4) <- gsub("^g ", "", colnames(avg_rna_full_v4))
write.csv(avg_rna_full_v4,
          file.path(out_avgexpr_v4, "AvgExpr_RNA_AllGenes_per_cluster.csv"))

# ── ADT ───────────────────────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"

avg_adt_raw_v4 <- AverageExpression(
  TARA_cd8, assays = "ADT", features = adt_available,
  group.by = "Manual_Annotation_refined", slot = "data"
)$ADT

avg_adt_scaled_v4 <- scale_01(as.matrix(log10(avg_adt_raw_v4 + 1)))

colnames(avg_adt_raw_v4)    <- gsub("^g ", "", colnames(avg_adt_raw_v4))
colnames(avg_adt_scaled_v4) <- gsub("^g ", "", colnames(avg_adt_scaled_v4))

write.csv(avg_adt_raw_v4,
          file.path(out_avgexpr_v4, "AvgExpr_ADT_per_cluster.csv"))
write.csv(round(avg_adt_scaled_v4, 4),
          file.path(out_avgexpr_v4, "AvgExpr_ADT_per_cluster_scaled.csv"))

# ── Percent expression ────────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "RNA"
pct_rna_v4 <- as.data.frame(
  t(sapply(levels(droplevels(TARA_cd8$Manual_Annotation_refined)), function(cl) {
    cells <- WhichCells(TARA_cd8, expression = Manual_Annotation_refined == cl)
    mat   <- GetAssayData(TARA_cd8, slot = "data")[rna_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_rna_v4, 2),
          file.path(out_avgexpr_v4, "PctExpr_RNA_per_cluster.csv"))

DefaultAssay(TARA_cd8) <- "ADT"
pct_adt_v4 <- as.data.frame(
  t(sapply(levels(droplevels(TARA_cd8$Manual_Annotation_refined)), function(cl) {
    cells <- WhichCells(TARA_cd8, expression = Manual_Annotation_refined == cl)
    mat   <- GetAssayData(TARA_cd8, slot = "data")[adt_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_adt_v4, 2),
          file.path(out_avgexpr_v4, "PctExpr_ADT_per_cluster.csv"))

message("✓ v4 expression tables saved to: ", out_avgexpr_v4)


################################################################################
# STEP 19: DIAGNOSTIC — NAIVE 1 vs NAIVE 2 COMPARISON
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 19: NAIVE 1 vs NAIVE 2 — KEY DIFFERENCES\n",
    "══════════════════════════════════════════════════════════════\n\n")

DefaultAssay(TARA_cd8) <- "ADT"
adt_cd8 <- GetAssayData(TARA_cd8, slot = "data")

naive1_cells <- WhichCells(TARA_cd8,
                           expression = Manual_Annotation_refined == "1a: Naive 1 CD8")
naive2_cells <- WhichCells(TARA_cd8,
                           expression = Manual_Annotation_refined == "6: Naive 2 CD8")

pct_n1 <- rowMeans(adt_cd8[, naive1_cells] > 0) * 100
pct_n2 <- rowMeans(adt_cd8[, naive2_cells] > 0) * 100

adt_diff <- data.frame(
  Marker  = names(pct_n1),
  Naive1  = round(pct_n1, 1),
  Naive2  = round(pct_n2, 1),
  Diff    = round(pct_n1 - pct_n2, 1),
  stringsAsFactors = FALSE
)
adt_diff <- adt_diff[order(-abs(adt_diff$Diff)), ]

cat("=== Top 20 ADT protein differences (Naive 1 vs Naive 2) ===\n")
cat(sprintf("%-16s %8s %8s %8s\n", "Marker", "Naive 1", "Naive 2", "Diff"))
cat(strrep("-", 44), "\n")
for (i in 1:min(20, nrow(adt_diff))) {
  r <- adt_diff[i, ]
  cat(sprintf("%-16s %7.1f%% %7.1f%% %+7.1f\n",
              r$Marker, r$Naive1, r$Naive2, r$Diff))
}
cat("\n")


################################################################################
# STEP 20: SAVE FINAL OBJECTS (v4)
################################################################################

cat("\n",
    "══════════════════════════════════════════════════════════════\n",
    " STEP 20: SAVING FINAL OBJECTS (v4)\n",
    "══════════════════════════════════════════════════════════════\n\n")

qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_sorted_v4.qs2"))
cat("✓ Saved: TARA_ALL_sorted_v4.qs2\n")

qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_revised_v4.qs2"))
cat("✓ Saved: TARA_cd8_revised_v4.qs2\n")


################################################################################
# DONE
################################################################################

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " Annotation v4 complete (final naming).\n",
        "══════════════════════════════════════════════════════════════\n",
        "\n",
        " FINAL CD8 CLUSTER DEFINITIONS:\n",
        "   1a: Naive 1 CD8            — CD45RA+ CD45RO- CD95- CCR7+\n",
        "                                 (CD73hi CD127hi ITGB7hi ICOShi)\n",
        "   1b: Tscm CD8               — CD45RA+ CD45RO- CD95+ CCR7+\n",
        "   1c: Naive Intermediate CD8  — CD45RA+ CD45RO+ (DP) CCR7+\n",
        "   6:  Naive 2 CD8            — CD45RA+ CD45RO- CCR7+ (quiescent)\n",
        "   8:  TEMRA/CTL CD8          — CD45RA+ CCR7- CD27- CD57+\n",
        "   9:  gd T cell              — TRDV1+/TRGV9+ ab TCR-\n",
        "  12a: Transitional CD8       — Activated, intermediate effector\n",
        "  27:  Tex CD8                — PD-1+ CD45RO+ TOX+ exhausted\n",
        "\n",
        " SAVED:\n",
        "   TARA_ALL_sorted_v4.qs2  — final full object\n",
        "   TARA_cd8_revised_v4.qs2 — final CD8 subset\n",
        "══════════════════════════════════════════════════════════════\n"
)

