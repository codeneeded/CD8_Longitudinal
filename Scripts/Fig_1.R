################################################################################
# FIGURE 1 — CD8 T Cell Annotation Validation & Landscape
# Manuscript: CD8 Longitudinal (TARA Cohort)
#
# PURPOSE:
#   1. Pull out Cluster 12 as additional CD8 cluster for evaluation
#   2. Generate violin + feature plots for ALL 108 ADT markers
#   3. Generate violin + feature plots for expanded RNA (CD8 + non-CD8 markers)
#   4. Export average expression tables (RNA + ADT) per cluster for annotation
#   5. Ensure ALL panels use the HEI/HEU/HUU subset (PreART_Entry, HEU, HUU)
#
# OUTPUT STRUCTURE:
#   Fig 1/
#     ├── 00_Annotation_Validation/        ← Supporting figures
#     │   ├── Violins_RNA/
#     │   ├── Violins_ADT/
#     │   ├── FeaturePlots_RNA/
#     │   ├── FeaturePlots_ADT/
#     │   ├── AvgExpression/
#     │   ├── DotPlot_*.png
#     │   └── Cluster12a_*.png / .csv
#     ├── 01_UMAP/                         ← Manuscript panels
#     ├── 02_ADT_RNA_Heatmap/
#     ├── 03_Cluster_Proportions/
#     └── 04_Viral_Load_Correlations/
#
# CLUSTER ORDERING (by cluster number):
#   Cluster 1  → Tcm CD8
#   Cluster 6  → Naïve CD8
#   Cluster 8  → TEMRA/CTL CD8
#   Cluster 9  → γδ T cell
#   Cluster 12 → Mixed → split into 12a (CD8) / 12b (CD4)
#   Cluster 27 → Tex CD8
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

# Manuscript figure panel directories (directly under Fig 1/)
out_umap  <- file.path(fig1_dir, "01_UMAP")
out_heat  <- file.path(fig1_dir, "02_ADT_RNA_Heatmap")
out_prop  <- file.path(fig1_dir, "03_Cluster_Proportions")
out_vl    <- file.path(fig1_dir, "04_Viral_Load_Correlations")

for (d in c(out_vln_rna, out_vln_adt, out_feat_rna, out_feat_adt, out_avgexpr,
            out_umap, out_heat, out_prop, out_vl)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Read data ──────────────────────────────────────────────────────────────────
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_sorted.qs2"))


################################################################################
# STEP 0: SUBSET TO HEI / HEU / HUU ONLY
#
# CRITICAL: All downstream analyses use this subset, not the full dataset.
################################################################################

keep_groups <- c("PreART_Entry", "HEU", "HUU")

TARA_sub <- subset(TARA_ALL, subset = Timepoint_Group %in% keep_groups)
TARA_sub$Condition <- factor(
  case_when(
    TARA_sub$Timepoint_Group == "PreART_Entry" ~ "HEI",
    TRUE ~ as.character(TARA_sub$Timepoint_Group)
  ),
  levels = c("HUU", "HEU", "HEI")
)

cat("=== Working subset (HEI/HEU/HUU only) ===\n")
print(table(TARA_sub$Condition))
cat("\n")


################################################################################
# STEP 1: CLUSTER ANNOTATION REFINEMENT
#
# Final annotations:
#   Cluster 1  → Tcm CD8
#   Cluster 6  → Naïve CD8
#   Cluster 8  → TEMRA/CTL CD8
#   Cluster 9  → γδ T cell
#   Cluster 12 → Mixed CD4/CD8 → split by CD8A expression into:
#                  12a: CD8 (CD8A+) — to be further annotated
#                  12b: CD4 (CD8A−) — excluded from CD8 analyses
#   Cluster 27 → Tex CD8
################################################################################

# ── Diagnostic: check what cluster 12 looks like in raw annotations ───────────
cat("=== All unique annotations containing '12' ===\n")
all_annot <- as.character(TARA_sub$Manual_Annotation)
cl12_raw <- unique(all_annot[grepl("12", all_annot)])
cat("  Manual_Annotation matches:", if (length(cl12_raw) > 0) paste(cl12_raw, collapse = "; ") else "NONE", "\n")
cat("  Cell count:", sum(grepl("12", all_annot)), "\n\n")

TARA_sub$Manual_Annotation_refined <- case_when(
  # Cluster 9: αβ TCR+ contaminants → merge into Cluster 8
  TARA_sub$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_sub$has_TCR == TRUE
  ~ "8: TEMRA/CTL CD8",
  # Cluster 9: confirmed γδ
  TARA_sub$Manual_Annotation == "9: TRDV1+ CTL-like" & TARA_sub$has_TCR == FALSE
  ~ "9: γδ T cell",
  # Cluster 1 (exact: "1: ..." not "10: ..." or "12: ...")
  grepl("^1: ", as.character(TARA_sub$Manual_Annotation))
  ~ "1: Tcm CD8",
  # Cluster 6
  grepl("^6: ", as.character(TARA_sub$Manual_Annotation))
  ~ "6: Naïve CD8",
  # Cluster 8 (original, before merge)
  grepl("^8: ", as.character(TARA_sub$Manual_Annotation))
  ~ "8: TEMRA/CTL CD8",
  # Cluster 27
  grepl("^27: ", as.character(TARA_sub$Manual_Annotation))
  ~ "27: Tex CD8",
  # Cluster 12: leave as-is for now (will split below)
  grepl("^12: ", as.character(TARA_sub$Manual_Annotation))
  ~ "12: Mixed CD4/CD8",
  # All other clusters unchanged
  TRUE ~ as.character(TARA_sub$Manual_Annotation)
)

# ── Split Cluster 12 by CD8A expression ───────────────────────────────────────
# Cells in Cluster 12 with CD8A > 0 are classified as CD8; rest as CD4
DefaultAssay(TARA_sub) <- "RNA"
cd8a_expr <- GetAssayData(TARA_sub, slot = "data")["CD8A", ]

# Diagnostic: check CD8A in cluster 12
is_cl12 <- TARA_sub$Manual_Annotation_refined == "12: Mixed CD4/CD8"
cat("=== CD8A expression in Cluster 12 ===\n")
cat("  Cells in '12: Mixed CD4/CD8':", sum(is_cl12), "\n")
cat("  CD8A > 0:", sum(cd8a_expr[is_cl12] > 0, na.rm = TRUE), "\n")
cat("  CD8A == 0:", sum(cd8a_expr[is_cl12] == 0, na.rm = TRUE), "\n")
cat("  CD8A is NA:", sum(is.na(cd8a_expr[is_cl12])), "\n")
cat("  CD8A summary:\n")
print(summary(cd8a_expr[is_cl12]))
cat("\n")

# Split: CD8A > 0 → 12a, otherwise → 12b (covers 0 and NA)
TARA_sub$Manual_Annotation_refined[is_cl12 & cd8a_expr > 0] <- "12a: CD8"
TARA_sub$Manual_Annotation_refined[is_cl12 & (cd8a_expr <= 0 | is.na(cd8a_expr))] <- "12b: CD4"

TARA_sub$Manual_Annotation_refined <- factor(TARA_sub$Manual_Annotation_refined)
TARA_sub$Cluster_Number_refined <- gsub("^([0-9]+[ab]?):.*", "\\1",
                                        TARA_sub$Manual_Annotation_refined)

cat("=== Post-refinement cluster sizes (CD8-relevant clusters) ===\n")
print(table(TARA_sub$Manual_Annotation_refined)[
  grepl("^1:|^6:|^8:|^9:|^12|^27:", names(table(TARA_sub$Manual_Annotation_refined)))
])
cat("\n")

cat("=== Cluster 12 split results ===\n")
n_12a <- sum(TARA_sub$Manual_Annotation_refined == "12a: CD8")
n_12b <- sum(TARA_sub$Manual_Annotation_refined == "12b: CD4")
cat("  12a (CD8A+):", n_12a, "cells\n")
cat("  12b (CD4):  ", n_12b, "cells\n")

if (n_12a == 0 & n_12b == 0) {
  cat("\n  ⚠ WARNING: No cells assigned to 12a or 12b!\n")
  cat("  This means grepl('^12: ', Manual_Annotation) didn't match any cells.\n")
  cat("  Check your Manual_Annotation column — cluster 12's label may not start with '12: '\n")
  cat("  Unique annotations that contain '12':\n")
  print(unique(grep("12", as.character(TARA_sub$Manual_Annotation), value = TRUE)))
  cat("  Unique refined annotations that contain '12':\n")
  print(unique(grep("12", as.character(TARA_sub$Manual_Annotation_refined), value = TRUE)))
  cat("\n  You may need to adjust the grepl pattern in Step 1 to match your data.\n\n")
}
cat("\n")


################################################################################
# STEP 2: DEFINE CD8 CLUSTERS
#
# Cluster 12a (CD8A+) is included; 12b (CD4) is excluded from CD8 analyses.
# 12a will be further characterized by the violin/feature/avgexpr outputs.
################################################################################

cd8_cluster_names <- c(
  "1: Tcm CD8",
  "6: Naïve CD8",
  "8: TEMRA/CTL CD8",
  "9: γδ T cell",
  "12a: CD8",
  "27: Tex CD8"
)

cd8_short_labels <- setNames(
  c("Tcm CD8", "Naïve CD8", "TEMRA/CTL", "γδ T cell", "Cl.12a CD8", "Tex CD8"),
  cd8_cluster_names
)

# ── Subset CD8 clusters from the HEI/HEU/HUU subset ──────────────────────────
TARA_cd8 <- subset(TARA_sub,
                   subset = Manual_Annotation_refined %in% cd8_cluster_names)

# Drop unused factor levels so VlnPlot2 only shows the 6 CD8 clusters
TARA_cd8$Manual_Annotation_refined <- droplevels(TARA_cd8$Manual_Annotation_refined)

# Reorder factor levels to match cd8_cluster_names order
TARA_cd8$Manual_Annotation_refined <- factor(
  TARA_cd8$Manual_Annotation_refined,
  levels = cd8_cluster_names[cd8_cluster_names %in% levels(TARA_cd8$Manual_Annotation_refined)]
)

cat("=== CD8 subset cell counts ===\n")
print(table(TARA_cd8$Manual_Annotation_refined))
cat("\n")

# Verify 12a is present
stopifnot("12a: CD8 not found in TARA_cd8 — check Cluster 12 split logic" =
            "12a: CD8" %in% levels(TARA_cd8$Manual_Annotation_refined))

# ── Shared color palettes ──────────────────────────────────────────────────────
cond_cols <- c("HUU" = "#4E79A7", "HEU" = "#F28E2B", "HEI" = "#E15759")


################################################################################
# STEP 3: GENE / PROTEIN MARKER LISTS
#
# For VIOLIN + FEATURE PLOTS:
#   ADT → ALL 108 markers in the TotalSeq-C panel (no curated subset)
#   RNA → Expanded list including non-CD8 lineage markers
#
# For HEATMAPS / DOTPLOTS:
#   Curated CD8-specific subsets (defined later in those sections)
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# ADT: Use ALL markers in the ADT assay (all 108 TotalSeq-C)
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(TARA_cd8) <- "ADT"
all_adt_markers <- sort(rownames(TARA_cd8))

cat("=== Total ADT markers in assay:", length(all_adt_markers), "===\n")
cat(paste(all_adt_markers, collapse = ", "), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# RNA: Comprehensive list — CD8 markers + non-CD8 lineage markers
#
# Grouped for stacked violin plots; the flat list is used for feature plots
# ══════════════════════════════════════════════════════════════════════════════

rna_markers <- list(
  # ── CD8 Core Identity ──
  `CD8 Identity` = c(
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"
  ),
  # ── Naïve / Quiescence ──
  `Naive - Quiescence` = c(
    "CCR7", "SELL", "TCF7", "LEF1", "BACH2", "IL7R", "BCL2",
    "S1PR1", "KLF2", "FOXP1", "SATB1"
  ),
  # ── Effector / Cytotoxicity ──
  `Effector - Cytotoxicity` = c(
    "GZMB", "GZMA", "GZMH", "GZMK", "GZMM",
    "PRF1", "GNLY", "NKG7", "KLRG1", "CX3CR1",
    "FGFBP2", "FCGR3A"
  ),
  # ── Effector TFs ──
  `Effector TFs` = c(
    "TBX21", "EOMES", "RUNX3", "PRDM1", "ZEB2", "ID2"
  ),
  # ── Exhaustion / Inhibitory ──
  `Exhaustion` = c(
    "TOX", "TOX2", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
    "CTLA4", "ENTPD1", "CD160", "CD244"
  ),
  # ── Memory / Stemness ──
  `Memory - Stemness` = c(
    "IL7R", "CD27", "CD28", "BACH2", "TCF7", "BCL2",
    "BCL6", "CXCR5", "CCR7"
  ),
  # ── Activation / Proliferation ──
  `Activation - Proliferation` = c(
    "MKI67", "HLA-DRA", "HLA-DRB1", "CD38", "FAS",
    "ICOS", "TNFRSF9", "CXCR6", "CD69"
  ),
  # ── Tissue Residency ──
  `Tissue Residency` = c(
    "ITGAE", "ITGA1", "CD69", "CXCR6", "ZNF683", "PRDM1"
  ),
  # ── γδ TCR ──
  `GammaDelta TCR` = c(
    "TRDV1", "TRDV2", "TRGV9", "TRDC", "TRGC1", "TRGC2"
  ),
  # ── MAIT ──
  `MAIT` = c(
    "TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "RORC"
  ),
  # ── NK-like / Innate ──
  `NK-like - Innate` = c(
    "TYROBP", "KLRD1", "KIR2DL3", "KIR3DL1", "NCAM1",
    "KLRF1", "KLRC1", "KLRC2"
  ),
  # ── Type I IFN ──
  `Type I IFN` = c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44L", "MX1",
    "OAS1", "STAT1", "IRF7"
  ),
  # ── Stress / Heat Shock ──
  `Stress - Heat Shock` = c(
    "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "IFI27"
  ),
  # ── Chemokines / Cytokines ──
  `Chemokines - Cytokines` = c(
    "CCL3", "CCL4", "CCL4L2", "CCL5", "XCL1", "XCL2",
    "IFNG", "TNF", "IL2", "CSF2"
  ),
  # ── Terminal Differentiation / Senescence ──
  `Terminal Differentiation` = c(
    "B3GAT1", "KLRG1", "CX3CR1", "ZEB2", "GZMB"
  ),
  
  # ══════════════════════════════════════════════════════════════════════════
  # NON-CD8 LINEAGE MARKERS (for ruling out contamination / confirming purity)
  # ══════════════════════════════════════════════════════════════════════════
  
  # ── CD4 T cell ──
  `CD4 T cell` = c(
    "CD4", "IL7R", "FOXP3", "IL2RA", "CTLA4", "IKZF2",
    "CCR4", "RORC", "GATA3", "TBX21", "BCL6",
    "CXCR5", "MAF", "ICOS", "BATF"
  ),
  # ── Treg ──
  `Treg` = c(
    "FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT"
  ),
  # ── B cell ──
  `B cell` = c(
    "CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1",
    "BLK", "BLNK", "TCL1A", "IGHM", "IGHD", "IGKC", "IGLC2"
  ),
  # ── Plasma cell ──
  `Plasma cell` = c(
    "SDC1", "XBP1", "MZB1", "JCHAIN", "IGHG1", "IGHA1",
    "PRDM1", "IRF4"
  ),
  # ── NK cell ──
  `NK cell` = c(
    "NCAM1", "KLRF1", "KLRC1", "KLRB1", "NCR1", "NCR3",
    "FCGR3A", "TYROBP", "FCER1G", "SH2D1B",
    "SPON2", "CLIC3", "MYOM2"
  ),
  # ── Monocyte / Macrophage ──
  `Monocyte - Macrophage` = c(
    "CD14", "LYZ", "S100A8", "S100A9", "S100A12",
    "FCGR3A", "CSF1R", "CD68", "MARCO",
    "VCAN", "FCN1", "MNDA"
  ),
  # ── cDC (conventional dendritic cell) ──
  `cDC` = c(
    "FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1",
    "BATF3", "IRF8", "IRF4", "LAMP3"
  ),
  # ── pDC (plasmacytoid dendritic cell) ──
  `pDC` = c(
    "LILRA4", "CLEC4C", "IRF7", "TCF4", "GZMB",
    "IL3RA", "JCHAIN"
  ),
  # ── Platelet / Megakaryocyte ──
  `Platelet` = c(
    "PPBP", "PF4", "GP9", "ITGA2B", "TUBB1", "TREML1"
  ),
  # ── Erythrocyte (RBC contamination) ──
  `Erythrocyte` = c(
    "HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA"
  ),
  # ── HSC / Progenitor ──
  `HSC - Progenitor` = c(
    "CD34", "SPINK2", "AVP", "CRHBP", "CYTL1",
    "SOX4", "HOPX"
  ),
  # ── ILC (Innate Lymphoid Cell) ──
  `ILC` = c(
    "IL7R", "KIT", "GATA3", "RORC", "IL1RL1",
    "AREG", "IL22", "NCR2"
  )
)

# Flatten RNA markers
rna_markers_flat <- unique(unlist(rna_markers))

cat("Total unique RNA markers (CD8 + non-CD8):", length(rna_markers_flat), "\n")
cat("Total ADT markers (all in assay):", length(all_adt_markers), "\n\n")


################################################################################
# STEP 4: CHECK WHICH MARKERS EXIST IN THE OBJECT
################################################################################

DefaultAssay(TARA_cd8) <- "RNA"
rna_available <- rna_markers_flat[rna_markers_flat %in% rownames(TARA_cd8)]
rna_missing   <- setdiff(rna_markers_flat, rna_available)

cat("=== RNA markers: ", length(rna_available), " available, ",
    length(rna_missing), " missing ===\n")
if (length(rna_missing) > 0) cat("Missing:", paste(rna_missing, collapse = ", "), "\n\n")

# ADT: all markers are by definition in the assay, but double-check
DefaultAssay(TARA_cd8) <- "ADT"
adt_available <- all_adt_markers  # all 108
cat("=== ADT markers: ", length(adt_available), " (all in assay) ===\n\n")


################################################################################
# STEP 5: EXPORT AVERAGE EXPRESSION TABLES (for pasting back for annotation)
#
# Raw and min-max scaled, both RNA and ADT, per cluster
################################################################################

message("Exporting average expression tables...")

# ── RNA average expression ────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "RNA"

# Use all RNA genes available from our marker list for the avg expression table
avg_rna_raw <- AverageExpression(
  TARA_cd8,
  assays   = "RNA",
  features = rna_available,
  group.by = "Manual_Annotation_refined",
  slot     = "data"
)$RNA

# Also compute for ALL genes (full transcriptome) — useful for discovery
avg_rna_full <- AverageExpression(
  TARA_cd8,
  assays   = "RNA",
  group.by = "Manual_Annotation_refined",
  slot     = "data"
)$RNA

# Min-max scale (row-wise, 0–1)
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

# Clean column names
colnames(avg_rna_raw)    <- gsub("^g ", "", colnames(avg_rna_raw))
colnames(avg_rna_scaled) <- gsub("^g ", "", colnames(avg_rna_scaled))
colnames(avg_rna_full)   <- gsub("^g ", "", colnames(avg_rna_full))

write.csv(avg_rna_raw,
          file.path(out_avgexpr, "AvgExpr_RNA_per_cluster.csv"),
          row.names = TRUE)
write.csv(round(avg_rna_scaled, 4),
          file.path(out_avgexpr, "AvgExpr_RNA_per_cluster_scaled.csv"),
          row.names = TRUE)
write.csv(avg_rna_full,
          file.path(out_avgexpr, "AvgExpr_RNA_AllGenes_per_cluster.csv"),
          row.names = TRUE)

# ── ADT average expression ────────────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"

avg_adt_raw <- AverageExpression(
  TARA_cd8,
  assays   = "ADT",
  features = adt_available,
  group.by = "Manual_Annotation_refined",
  slot     = "data"
)$ADT

avg_adt_scaled <- scale_01(as.matrix(log10(avg_adt_raw + 1)))

colnames(avg_adt_raw)    <- gsub("^g ", "", colnames(avg_adt_raw))
colnames(avg_adt_scaled) <- gsub("^g ", "", colnames(avg_adt_scaled))

write.csv(avg_adt_raw,
          file.path(out_avgexpr, "AvgExpr_ADT_per_cluster.csv"),
          row.names = TRUE)
write.csv(round(avg_adt_scaled, 4),
          file.path(out_avgexpr, "AvgExpr_ADT_per_cluster_scaled.csv"),
          row.names = TRUE)

# ── Also export percent expression (% cells expressing each marker) ───────────
# This is very useful for annotation — a gene with high avg but low % is
# driven by a few outlier cells, not a real cluster signature.
DefaultAssay(TARA_cd8) <- "RNA"
pct_rna <- as.data.frame(
  t(sapply(levels(droplevels(TARA_cd8$Manual_Annotation_refined)), function(cl) {
    cells <- WhichCells(TARA_cd8, expression = Manual_Annotation_refined == cl)
    mat   <- GetAssayData(TARA_cd8, slot = "data")[rna_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_rna, 2),
          file.path(out_avgexpr, "PctExpr_RNA_per_cluster.csv"),
          row.names = TRUE)

DefaultAssay(TARA_cd8) <- "ADT"
pct_adt <- as.data.frame(
  t(sapply(levels(droplevels(TARA_cd8$Manual_Annotation_refined)), function(cl) {
    cells <- WhichCells(TARA_cd8, expression = Manual_Annotation_refined == cl)
    mat   <- GetAssayData(TARA_cd8, slot = "data")[adt_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_adt, 2),
          file.path(out_avgexpr, "PctExpr_ADT_per_cluster.csv"),
          row.names = TRUE)

message("✓ Average expression tables saved to: ", out_avgexpr)
cat("  Files:\n",
    "    AvgExpr_RNA_per_cluster.csv          (curated markers, raw)\n",
    "    AvgExpr_RNA_per_cluster_scaled.csv    (curated markers, 0-1 scaled)\n",
    "    AvgExpr_RNA_AllGenes_per_cluster.csv  (full transcriptome, raw)\n",
    "    AvgExpr_ADT_per_cluster.csv           (all 108 ADT, raw)\n",
    "    AvgExpr_ADT_per_cluster_scaled.csv    (all 108 ADT, 0-1 scaled)\n",
    "    PctExpr_RNA_per_cluster.csv           (% cells expressing, RNA)\n",
    "    PctExpr_ADT_per_cluster.csv           (% cells expressing, ADT)\n\n")


################################################################################
# STEP 6: VIOLIN PLOTS — ADT (ALL 108 markers)
#
# Using SeuratExtend VlnPlot2 with:
#   - cols = 'light' (light color scheme)

#   - show.mean = TRUE with mean_colors
#   - assay = "ADT"
#
# Individual per-marker violins for all 108 ADT markers
################################################################################

message("Generating ADT violin plots (all ", length(adt_available), " markers)...")

for (feat in adt_available) {
  
  p <- VlnPlot2(
    TARA_cd8,
    features     = feat,
    group.by     = "Manual_Annotation_refined",
    assay        = "ADT",
    cols         = "light",
    show.mean    = TRUE,
    mean_colors  = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(
    file.path(out_vln_adt, paste0("VlnADT_", safe_name, ".png")),
    p, width = 12, height = 5, dpi = 300, bg = "white"
  )
}

message("✓ All ADT violin plots saved (", length(adt_available), " markers)\n")


################################################################################
# STEP 7: VIOLIN PLOTS — RNA (CD8 + non-CD8 markers)
#
# Using SeuratExtend VlnPlot2 with:
#   - cols = 'light'
#   - show.mean = TRUE, mean_colors = c("red", "blue")
#   - assay = "RNA"
#
# Individual per-gene violins for all available RNA markers
################################################################################

message("Generating RNA violin plots (", length(rna_available), " markers)...")

for (gene in rna_available) {
  
  if (!gene %in% rownames(TARA_cd8[["RNA"]])) next
  
  p <- VlnPlot2(
    TARA_cd8,
    features     = gene,
    group.by     = "Manual_Annotation_refined",
    assay        = "RNA",
    cols         = "light",
    show.mean    = TRUE,
    mean_colors  = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(
    file.path(out_vln_rna, paste0("VlnRNA_", safe_name, ".png")),
    p, width = 12, height = 5, dpi = 300, bg = "white"
  )
}

message("✓ All RNA violin plots saved (", length(rna_available), " markers)\n")


################################################################################
# STEP 8: FEATURE PLOTS — ADT (ALL 108 markers, DimPlot2 on WNN UMAP)
#
# Using SeuratExtend DimPlot2 with features parameter
################################################################################

message("Generating ADT feature plots (all ", length(adt_available), " markers)...")

DefaultAssay(TARA_cd8) <- "ADT"

for (feat in adt_available) {
  
  if (!feat %in% rownames(TARA_cd8[["ADT"]])) next
  
  fp <- DimPlot2(
    TARA_cd8,
    features  = feat,
    reduction = "wnn.umap"
  ) + ggtitle(paste("ADT |", feat))
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(
    file.path(out_feat_adt, paste0("FeatADT_", safe_name, ".png")),
    fp, dpi = 500, width = 8, height = 6, bg = "white"
  )
}

message("✓ All ADT feature plots saved (", length(adt_available), " markers)\n")


################################################################################
# STEP 9: FEATURE PLOTS — RNA (CD8 + non-CD8, DimPlot2 on WNN UMAP)
#
# Using SeuratExtend DimPlot2 with features parameter
################################################################################

message("Generating RNA feature plots (", length(rna_available), " markers)...")

DefaultAssay(TARA_cd8) <- "RNA"

for (gene in rna_available) {
  
  if (!gene %in% rownames(TARA_cd8[["RNA"]])) next
  
  fp <- DimPlot2(
    TARA_cd8,
    features  = gene,
    reduction = "wnn.umap"
  ) + ggtitle(paste("RNA |", gene))
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(
    file.path(out_feat_rna, paste0("FeatRNA_", safe_name, ".png")),
    fp, dpi = 500, width = 8, height = 6, bg = "white"
  )
}

message("✓ All RNA feature plots saved (", length(rna_available), " markers)\n")


################################################################################
# STEP 10: COMBINED OVERVIEW — DotPlots (curated subsets for quick reference)
################################################################################

message("Generating DotPlots...")

# ── RNA DotPlot (curated CD8 identity markers) ───────────────────────────────
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
    "MKI67", "HLA-DRA", "CD69",
    "ISG15", "IFIT1", "IFI44L",
    "CCL3", "CCL4", "CCL5", "IFNG",
    # Non-CD8 exclusion markers
    "CD4", "FOXP3", "CD19", "MS4A1", "CD79A",
    "CD14", "LYZ", "S100A8",
    "NCAM1", "KLRF1",
    "PPBP", "PF4", "HBB"),
  rna_available
)

DefaultAssay(TARA_cd8) <- "RNA"

p_dot_rna <- DotPlot(
  TARA_cd8,
  features = rna_dotplot_genes,
  group.by = "Manual_Annotation_refined",
  cols     = c("lightgrey", "#D73027"),
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("RNA: CD8 Cluster Identity + Lineage Exclusion Markers") +
  theme(
    plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(out_annot, "DotPlot_RNA_CD8_Clusters.png"),
  p_dot_rna, width = 24, height = 7, dpi = 300, bg = "white"
)

# ── ADT DotPlot (all 108 markers) ────────────────────────────────────────────
DefaultAssay(TARA_cd8) <- "ADT"

p_dot_adt <- DotPlot(
  TARA_cd8,
  features = adt_available,
  group.by = "Manual_Annotation_refined",
  cols     = c("lightgrey", "#08519C"),
  dot.scale = 5
) +
  RotatedAxis() +
  ggtitle("ADT: All 108 Surface Protein Markers") +
  theme(
    plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(out_annot, "DotPlot_ADT_All108_CD8_Clusters.png"),
  p_dot_adt, width = 36, height = 7, dpi = 300, bg = "white"
)

message("✓ DotPlots saved\n")


################################################################################
# STEP 11: CLUSTER 12a (CD8A+) — ANNOTATION DIAGNOSTIC
#
# Cluster 12 was a mixed CD4/CD8 population. After splitting by CD8A,
# 12a contains the CD8+ cells. This section characterizes 12a to determine
# its identity (e.g., Tem, effector, activated, etc.)
################################################################################

message("Generating Cluster 12a annotation diagnostics...")

DefaultAssay(TARA_cd8) <- "RNA"

# ── DE: Cluster 12a vs all other CD8 clusters ─────────────────────────────────
TARA_cd8$is_cl12a <- ifelse(
  TARA_cd8$Manual_Annotation_refined == "12a: CD8",
  "Cluster 12a", "Other CD8"
)

cl12a_markers <- FindMarkers(
  TARA_cd8,
  ident.1         = "Cluster 12a",
  group.by        = "is_cl12a",
  test.use        = "MAST",
  min.pct         = 0.1,
  logfc.threshold = 0.25
)

write.csv(cl12a_markers,
          file.path(out_annot, "Cluster12a_vs_OtherCD8_DEgenes.csv"),
          row.names = TRUE)

top_up   <- head(rownames(cl12a_markers[cl12a_markers$avg_log2FC > 0, ]), 30)
top_down <- head(rownames(cl12a_markers[cl12a_markers$avg_log2FC < 0, ]), 30)

cat("=== Cluster 12a (CD8A+) — Top 30 upregulated genes ===\n")
print(top_up)
cat("\n=== Cluster 12a (CD8A+) — Top 30 downregulated genes ===\n")
print(top_down)
cat("\n")

# ── Feature plots: top 6 upregulated genes ────────────────────────────────────
top6 <- head(top_up, 6)
if (length(top6) > 0) {
  cl12a_feat_plots <- lapply(top6, function(g) {
    DimPlot2(TARA_cd8, features = g, reduction = "wnn.umap") +
      ggtitle(paste("RNA |", g))
  })
  p_cl12a_feat <- wrap_plots(cl12a_feat_plots, ncol = 3)
  
  ggsave(file.path(out_annot, "Cluster12a_Top6_FeaturePlot.png"),
         p_cl12a_feat, width = 20, height = 12, dpi = 300, bg = "white")
}

# ── Violin plots: key genes for CD8 subtype discrimination ────────────────────
cl12a_annotation_genes <- intersect(
  c(# Core identity
    "CD8A", "CD8B", "CD3D", "CD3E", "CD4",
    # Naïve
    "CCR7", "SELL", "TCF7", "LEF1", "BACH2", "S1PR1",
    # Memory
    "IL7R", "BCL2",
    # Effector / cytotoxic
    "GZMB", "GZMA", "GZMK", "GZMH", "PRF1", "GNLY", "NKG7",
    # Effector TFs
    "TBX21", "EOMES", "RUNX3", "ZEB2",
    # Exhaustion
    "TOX", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
    # Activation
    "MKI67", "HLA-DRA", "CD38", "FAS", "CD69",
    # Chemokines
    "CCL3", "CCL4", "CCL5", "IFNG",
    # Terminal diff
    "B3GAT1", "CX3CR1", "KLRG1",
    # IFN
    "ISG15", "IFIT1",
    # γδ exclusion
    "TRDV1", "TRDC"),
  rna_available
)

for (gene in cl12a_annotation_genes) {
  p_cl12a <- VlnPlot2(
    TARA_cd8,
    features     = gene,
    group.by     = "Manual_Annotation_refined",
    assay        = "RNA",
    cols         = "light",
    show.mean    = TRUE,
    mean_colors  = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(file.path(out_annot, paste0("Cluster12a_Vln_", safe_name, ".png")),
         p_cl12a, width = 12, height = 5, dpi = 300, bg = "white")
}

# ── Also run DE: 12a vs each individual CD8 cluster (pairwise) ────────────────
# This helps see what makes 12a distinct from each specific population
other_cd8_clusters <- setdiff(cd8_cluster_names, "12a: CD8")

for (ref_cl in other_cd8_clusters) {
  
  # Subset to just 12a + reference cluster
  pair_cells <- WhichCells(TARA_cd8, expression =
                             Manual_Annotation_refined %in% c("12a: CD8", ref_cl))
  TARA_pair <- subset(TARA_cd8, cells = pair_cells)
  
  de_pair <- tryCatch({
    FindMarkers(
      TARA_pair,
      ident.1         = "12a: CD8",
      group.by        = "Manual_Annotation_refined",
      test.use        = "MAST",
      min.pct         = 0.1,
      logfc.threshold = 0.25
    )
  }, error = function(e) {
    message("  ⚠ DE failed for 12a vs ", ref_cl, ": ", e$message)
    NULL
  })
  
  if (!is.null(de_pair)) {
    ref_safe <- gsub("[^A-Za-z0-9._-]", "_", ref_cl)
    write.csv(de_pair,
              file.path(out_annot, paste0("Cluster12a_vs_", ref_safe, "_DEgenes.csv")),
              row.names = TRUE)
    
    cat("=== 12a vs", ref_cl, ": top 10 UP ===\n")
    print(head(rownames(de_pair[de_pair$avg_log2FC > 0, ]), 10))
    cat("\n")
  }
}

TARA_cd8$is_cl12a <- NULL

message("✓ Cluster 12a diagnostics saved\n")


################################################################################
# FIG 1B — UMAP (HEI/HEU/HUU subset only)
################################################################################

message("Generating Fig 1B UMAP...")

p_umap_all <- DimPlot2(
  TARA_sub,
  reduction  = "wnn.umap",
  group.by   = "Cluster_Number_refined",
  cols       = "default",
  label      = TRUE,
  box        = TRUE,
  label.size = 6,
  repel      = TRUE,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
)

ggsave(file.path(out_umap, "Fig1B_UMAP_AllClusters.png"),
       p_umap_all, width = 6, height = 6, dpi = 300, bg = "white")

p_umap_split <- DimPlot2(
  TARA_sub,
  reduction  = "wnn.umap",
  group.by   = "Cluster_Number_refined",
  split.by   = "Condition",
  cols       = "default",
  label      = TRUE,
  label.size = 4,
  repel      = TRUE,
  pt.size    = 0.3,
  raster     = FALSE,
  ncol       = 3,
  theme      = list(NoLegend(), NoAxes())
)

ggsave(file.path(out_umap, "Fig1B_UMAP_AllClusters_SplitByCondition.png"),
       p_umap_split, width = 16, height = 6, dpi = 300, bg = "white")

message("✓ Fig 1B UMAP saved\n")


################################################################################
# FIG 1C — ADT + RNA Heatmap (curated markers, now with Cluster 12)
################################################################################

message("Generating Fig 1C heatmaps...")

adt_heatmap_markers <- c(
  "CD7"       = "Naïve",      "SELL"    = "Naïve",      "CD45RA"  = "Naïve",
  "TCR-vA7.2" = "TCR Identity", "TCR-AB" = "TCR Identity", "TCR-vD2" = "TCR Identity",
  "IL7R"      = "Memory/Stemness", "CD38" = "Memory/Stemness", "NT5E" = "Memory/Stemness",
  "ENTPD1"    = "Memory/Stemness", "CD45RO" = "Memory/Stemness",
  "CD27"      = "Co-stimulation", "CD28" = "Co-stimulation",
  "B3GAT1"    = "Effector/TEMRA", "KLRG1" = "Effector/TEMRA", "KIR3DL1" = "Effector/TEMRA",
  "NCAM1"     = "NK-like/Innate", "FCGR3A" = "NK-like/Innate", "SIGLEC7" = "NK-like/Innate",
  "TIGIT"     = "Exhaustion",    "PDCD1"  = "Exhaustion",    "LAG3"    = "Exhaustion"
)

rna_heatmap_markers <- c(
  "CCR7"    = "Naïve/Stemness", "SELL"  = "Naïve/Stemness", "TCF7"  = "Naïve/Stemness",
  "LEF1"    = "Naïve/Stemness",
  "TRDV1"   = "γδ TCR",   "TRGV9" = "γδ TCR",   "TRDC" = "γδ TCR",
  "IL7R"    = "Memory/Stemness", "BCL2" = "Memory/Stemness", "BACH2" = "Memory/Stemness",
  "GZMK"    = "Effector/TEMRA", "GZMB"  = "Effector/TEMRA", "GNLY" = "Effector/TEMRA",
  "PRF1"    = "Effector/TEMRA", "NKG7"  = "Effector/TEMRA",
  "TBX21"   = "Effector TFs",  "EOMES" = "Effector TFs",  "RUNX3" = "Effector TFs",
  "TOX"     = "Exhaustion",    "PDCD1" = "Exhaustion",    "TIGIT" = "Exhaustion",
  "HAVCR2"  = "Exhaustion",
  "TYROBP"  = "NK-like/Innate", "KLRD1" = "NK-like/Innate", "KIR2DL3" = "NK-like/Innate",
  "MKI67"   = "Activation",    "CXCR6" = "Activation",    "HLA-DRA" = "Activation",
  "CD8A"    = "CD8 Identity",  "CD8B"  = "CD8 Identity",  "CD3E"    = "CD8 Identity"
)

# Filter to available
adt_heatmap_feats <- intersect(names(adt_heatmap_markers), adt_available)
rna_heatmap_feats <- intersect(names(rna_heatmap_markers), rna_available)

group_colors <- c(
  "Naïve"           = "#74C2E1",
  "Naïve/Stemness"  = "#74C2E1",
  "TCR Identity"    = "#9C6FD6",
  "γδ TCR"          = "#9C6FD6",
  "Memory/Stemness" = "#66BB6A",
  "Co-stimulation"  = "#AED581",
  "Effector/TEMRA"  = "#EF5350",
  "Effector TFs"    = "#FFA726",
  "NK-like/Innate"  = "#F06292",
  "Exhaustion"      = "#A1887F",
  "Activation"      = "#FFD54F",
  "CD8 Identity"    = "#42A5F5"
)

DefaultAssay(TARA_cd8) <- "ADT"
avg_adt_hm <- AverageExpression(
  TARA_cd8, assays = "ADT", features = adt_heatmap_feats,
  group.by = "Manual_Annotation_refined", slot = "data"
)$ADT

DefaultAssay(TARA_cd8) <- "RNA"
avg_rna_hm <- AverageExpression(
  TARA_cd8, assays = "RNA", features = rna_heatmap_feats,
  group.by = "Manual_Annotation_refined", slot = "data"
)$RNA

avg_adt_hm_sc <- scale_01(as.matrix(log10(avg_adt_hm + 1)))
avg_rna_hm_sc <- scale_01(as.matrix(avg_rna_hm))

colnames(avg_adt_hm_sc) <- gsub("^g ", "", colnames(avg_adt_hm_sc))
colnames(avg_rna_hm_sc) <- gsub("^g ", "", colnames(avg_rna_hm_sc))

# Column ordering (dynamic: find columns matching cluster labels)
# Seurat AverageExpression column names are derived from factor levels,
# so they'll contain the annotation string. We match by the cluster prefix.
get_col_for_cluster <- function(colnames_vec, cluster_prefix) {
  # Try matching "_<prefix>$" pattern (e.g., "_12a$" or "_1$")
  pattern <- paste0("_", cluster_prefix, "$")
  matched <- grep(pattern, colnames_vec, value = TRUE)
  if (length(matched) == 1) return(matched)
  # Try matching "^<prefix>:" pattern
  pattern2 <- paste0("^", cluster_prefix, ":")
  matched2 <- grep(pattern2, colnames_vec, value = TRUE)
  if (length(matched2) == 1) return(matched2)
  # Try partial match anywhere
  pattern3 <- cluster_prefix
  matched3 <- grep(pattern3, colnames_vec, value = TRUE, fixed = TRUE)
  if (length(matched3) == 1) return(matched3)
  return(NA_character_)
}

target_clusters <- c("1", "6", "8", "9", "12a", "27")
col_order_adt <- na.omit(sapply(target_clusters, function(cl)
  get_col_for_cluster(colnames(avg_adt_hm_sc), cl)))
col_order_rna <- na.omit(sapply(target_clusters, function(cl)
  get_col_for_cluster(colnames(avg_rna_hm_sc), cl)))

avg_adt_hm_sc <- avg_adt_hm_sc[, col_order_adt, drop = FALSE]
avg_rna_hm_sc <- avg_rna_hm_sc[, col_order_rna, drop = FALSE]

make_display_label <- function(raw_name) {
  # Extract cluster ID: everything after the last underscore, or before the colon
  num <- sub(".*_([0-9a-z]+)$", "\\1", raw_name)
  if (num == raw_name) num <- sub("^([0-9a-z]+):.*", "\\1", raw_name)
  name_map <- c("1" = "Tcm CD8", "6" = "Naïve CD8", "8" = "TEMRA/CTL",
                "9" = "γδ T cell", "12a" = "Cl.12a CD8", "27" = "Tex CD8")
  short <- ifelse(num %in% names(name_map), name_map[num], paste0("Cl.", num))
  paste0(short, "\n(Cluster ", num, ")")
}

colnames(avg_adt_hm_sc) <- sapply(col_order_adt, make_display_label)
colnames(avg_rna_hm_sc) <- sapply(col_order_rna, make_display_label)

# Row annotations
adt_annot <- data.frame(
  Group = factor(adt_heatmap_markers[rownames(avg_adt_hm_sc)], levels = names(group_colors)),
  row.names = rownames(avg_adt_hm_sc)
)
rna_annot <- data.frame(
  Group = factor(rna_heatmap_markers[rownames(avg_rna_hm_sc)], levels = names(group_colors)),
  row.names = rownames(avg_rna_hm_sc)
)

cluster_within_groups <- function(mat, annot_df) {
  groups <- levels(droplevels(annot_df$Group))
  row_order <- character(0)
  for (grp in groups) {
    rows <- rownames(annot_df)[annot_df$Group == grp]
    if (length(rows) <= 1) { row_order <- c(row_order, rows)
    } else {
      hc <- hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2")
      row_order <- c(row_order, rows[hc$order])
    }
  }
  row_order
}

adt_row_order <- cluster_within_groups(avg_adt_hm_sc, adt_annot)
rna_row_order <- cluster_within_groups(avg_rna_hm_sc, rna_annot)

avg_adt_ordered   <- avg_adt_hm_sc[adt_row_order, ]
avg_rna_ordered   <- avg_rna_hm_sc[rna_row_order, ]
adt_annot_ordered <- adt_annot[adt_row_order, , drop = FALSE]
rna_annot_ordered <- rna_annot[rna_row_order, , drop = FALSE]

get_gaps <- function(annot_df) {
  grps <- as.character(annot_df$Group)
  which(grps[-length(grps)] != grps[-1])
}

annot_colors_adt <- list(Group = group_colors[levels(droplevels(adt_annot_ordered$Group))])
annot_colors_rna <- list(Group = group_colors[levels(droplevels(rna_annot_ordered$Group))])

heatmap_colors <- colorRampPalette(c("#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"))(100)

p_adt <- pheatmap(
  avg_adt_ordered, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
  color = heatmap_colors, border_color = "white",
  annotation_row = adt_annot_ordered, annotation_colors = annot_colors_adt,
  annotation_names_row = FALSE, gaps_row = get_gaps(adt_annot_ordered),
  cellwidth = 60, cellheight = 16, fontsize = 11, fontsize_row = 10,
  fontsize_col = 10, angle_col = 0, main = "Protein (ADT)", silent = TRUE
)

p_rna <- pheatmap(
  avg_rna_ordered, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
  color = heatmap_colors, border_color = "white",
  annotation_row = rna_annot_ordered, annotation_colors = annot_colors_rna,
  annotation_names_row = FALSE, gaps_row = get_gaps(rna_annot_ordered),
  cellwidth = 60, cellheight = 16, fontsize = 11, fontsize_row = 10,
  fontsize_col = 10, angle_col = 0, main = "mRNA (RNA)", silent = TRUE
)

png(file.path(out_heat, "Fig1C_ADT_RNA_Heatmap_CD8_Combined.png"),
    width = 24, height = 16, units = "in", res = 300, bg = "white")
grid.newpage()
grid.text("CD8 T Cell Cluster Validation: Protein (ADT) & mRNA",
          x = 0.5, y = 0.975, gp = gpar(fontsize = 18, fontface = "bold"))
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1)); grid.draw(p_adt$gtable); popViewport()
pushViewport(viewport(layout.pos.col = 2)); grid.draw(p_rna$gtable); popViewport()
dev.off()

png(file.path(out_heat, "Fig1C_ADT_Heatmap_CD8.png"),
    width = 12, height = 10, units = "in", res = 300, bg = "white")
grid.draw(p_adt$gtable)
dev.off()

png(file.path(out_heat, "Fig1C_RNA_Heatmap_CD8.png"),
    width = 12, height = 14, units = "in", res = 300, bg = "white")
grid.draw(p_rna$gtable)
dev.off()

message("✓ Fig 1C heatmaps saved\n")


################################################################################
# FIG 1D — CD8 Cluster Proportions (% of total PBMC)
################################################################################

message("Generating Fig 1D proportions...")

meta_sub <- as.data.frame(TARA_sub@meta.data)

sample_totals_all <- meta_sub %>%
  group_by(orig.ident, Condition) %>%
  summarise(total_pbmc = dplyr::n(), .groups = "drop")

cluster_freq <- meta_sub %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  group_by(orig.ident, Condition, Manual_Annotation_refined) %>%
  summarise(cluster_cells = dplyr::n(), .groups = "drop") %>%
  group_by(orig.ident, Condition) %>%
  tidyr::complete(
    Manual_Annotation_refined = cd8_cluster_names,
    fill = list(cluster_cells = 0)
  ) %>%
  ungroup() %>%
  left_join(sample_totals_all, by = c("orig.ident", "Condition")) %>%
  mutate(
    percent = 100 * cluster_cells / total_pbmc,
    Facet_Label = factor(
      cd8_short_labels[as.character(Manual_Annotation_refined)],
      levels = cd8_short_labels
    )
  )

facet_max <- cluster_freq %>%
  group_by(Facet_Label) %>%
  summarise(max_percent = max(percent, na.rm = TRUE), .groups = "drop")

stat_df <- cluster_freq %>%
  group_by(Facet_Label) %>%
  wilcox_test(percent ~ Condition) %>%
  left_join(facet_max, by = "Facet_Label") %>%
  mutate(
    xmin = case_when(
      group1 == "HUU" & group2 == "HEU" ~ 1,
      group1 == "HUU" & group2 == "HEI" ~ 1,
      group1 == "HEU" & group2 == "HEI" ~ 2
    ),
    xmax = case_when(
      group1 == "HUU" & group2 == "HEU" ~ 2,
      group1 == "HUU" & group2 == "HEI" ~ 3,
      group1 == "HEU" & group2 == "HEI" ~ 3
    ),
    y.position = case_when(
      group1 == "HUU" & group2 == "HEU" ~ max_percent * 1.10,
      group1 == "HUU" & group2 == "HEI" ~ max_percent * 1.22,
      group1 == "HEU" & group2 == "HEI" ~ max_percent * 1.16
    ),
    p.signif = case_when(
      p <= 0.0001 ~ "****", p <= 0.001 ~ "***",
      p <= 0.01   ~ "**",   p <= 0.05  ~ "*",
      TRUE        ~ "ns"
    )
  ) %>%
  ungroup() %>%
  filter(p.signif != "ns")

p_prop <- ggplot(cluster_freq, aes(x = Condition, y = percent,
                                   fill = Condition, color = Condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 2.3, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", nrow = 1) +
  stat_pvalue_manual(stat_df, label = "p.signif", tip.length = 0.01,
                     bracket.size = 0.6, size = 4.5) +
  scale_fill_manual(values  = cond_cols) +
  scale_color_manual(values = cond_cols) +
  labs(x = NULL, y = "% of total PBMC") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 10),
    axis.text.x     = element_text(face = "bold", size = 10),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(face = "bold", size = 12)
  )

ggsave(file.path(out_prop, "Fig1D_CD8_Cluster_Proportions.png"),
       p_prop, width = 22, height = 6, dpi = 600, bg = "white")

message("✓ Fig 1D proportions saved\n")


################################################################################
# FIG 1E — CD8 Cluster % vs HIV Viral Load (HEI only)
################################################################################

message("Generating Fig 1E viral load correlations...")

meta_hei_vl <- as.data.frame(TARA_sub@meta.data) %>%
  filter(Condition == "HEI", !is.na(Viral_Load_num), Viral_Load_num > 0)

sample_totals_hei <- meta_hei_vl %>%
  group_by(orig.ident) %>%
  summarise(total_pbmc = dplyr::n(), .groups = "drop")

sample_vl <- meta_hei_vl %>%
  group_by(orig.ident) %>%
  summarise(Viral_Load_num = dplyr::first(Viral_Load_num), .groups = "drop") %>%
  mutate(log10_VL = log10(Viral_Load_num))

cluster_vl_df <- meta_hei_vl %>%
  filter(Manual_Annotation_refined %in% cd8_cluster_names) %>%
  group_by(orig.ident, Manual_Annotation_refined) %>%
  summarise(cluster_cells = dplyr::n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  tidyr::complete(
    Manual_Annotation_refined = cd8_cluster_names,
    fill = list(cluster_cells = 0)
  ) %>%
  ungroup() %>%
  left_join(sample_totals_hei, by = "orig.ident") %>%
  left_join(sample_vl, by = "orig.ident") %>%
  mutate(
    percent_pbmc = 100 * cluster_cells / total_pbmc,
    Facet_Label  = factor(
      cd8_short_labels[as.character(Manual_Annotation_refined)],
      levels = cd8_short_labels
    )
  )

cor_stats <- cluster_vl_df %>%
  group_by(Facet_Label) %>%
  cor_test(log10_VL, percent_pbmc, method = "spearman") %>%
  mutate(
    stat_label = paste0("Spearman's \u03c1 = ", sprintf("%.2f", cor),
                        "\np-value = ", signif(p, 2)),
    x = -Inf, y = Inf
  ) %>%
  ungroup()

p_vl <- ggplot(cluster_vl_df, aes(x = log10_VL, y = percent_pbmc)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              color = "grey40", fill = "grey85", alpha = 0.3) +
  geom_point(shape = 21, size = 3.2, stroke = 0.4,
             fill = "#C44E52", color = "white") +
  geom_text(data = cor_stats,
            aes(x = x, y = y, label = stat_label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.2,
            size = 4.5, lineheight = 1.4, color = "grey15", fontface = "bold") +
  facet_wrap(~ Facet_Label, scales = "free", nrow = 1) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = expression(log[10]~"viral load"), y = "% of total PBMC") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = 10, color = "black"),
    axis.title.x     = element_text(face = "bold", size = 12, margin = margin(t = 6)),
    axis.title.y     = element_text(face = "bold", size = 12, margin = margin(r = 6)),
    legend.position  = "none",
    panel.spacing    = unit(0.8, "lines")
  )

ggsave(file.path(out_vl, "Fig1E_CD8_Clusters_ViralLoad_Correlation.png"),
       p_vl, width = 24, height = 5.5, dpi = 600, bg = "white")

cat("=== Viral Load Correlations ===\n")
print(cor_stats[, c("Facet_Label", "cor", "statistic", "p")])

message("✓ Fig 1E viral load correlations saved\n")


################################################################################
# SUMMARY
################################################################################

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " All Figure 1 analysis outputs saved to:\n",
        " ", fig1_dir, "/Analysis/\n",
        "══════════════════════════════════════════════════════════════\n\n",
        " WHAT CHANGED FROM PREVIOUS VERSION:\n",
        "   1. ADT violins + feature plots: ALL ", length(adt_available), " markers (not curated subset)\n",
        "      - Stacked violins in batches of 15 for readability\n",
        "      - Individual per-marker violins + feature plots for all 108\n",
        "   2. RNA violins + feature plots: expanded with NON-CD8 lineage markers\n",
        "      - CD4, Treg, B cell, Plasma, NK, Monocyte, cDC, pDC,\n",
        "        Platelet, Erythrocyte, HSC/Progenitor, ILC\n",
        "   3. Average expression CSVs exported to AvgExpression/:\n",
        "      - AvgExpr_RNA_per_cluster.csv           (curated markers)\n",
        "      - AvgExpr_RNA_per_cluster_scaled.csv     (0-1 scaled)\n",
        "      - AvgExpr_RNA_AllGenes_per_cluster.csv   (FULL transcriptome)\n",
        "      - AvgExpr_ADT_per_cluster.csv            (all 108 ADT)\n",
        "      - AvgExpr_ADT_per_cluster_scaled.csv     (0-1 scaled)\n",
        "      - PctExpr_RNA_per_cluster.csv            (% cells expressing)\n",
        "      - PctExpr_ADT_per_cluster.csv            (% cells expressing)\n",
        "\n",
        " → Paste the CSV contents back here and I'll help annotate clusters.\n"
)

