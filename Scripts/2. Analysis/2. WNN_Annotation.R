################################################################################
# TARA — Unified WNN Annotation Pipeline (Final — corrected DSB)
#
# All annotations verified against correct DSB-normalized ADT data.
# Run to CHECKPOINT for CSV review, then continue.
#
# OUTPUT: ~/Documents/CD8_Longitudinal/Analysis/Annotation/
################################################################################


# =============================================================================
# STEP 0: LIBRARIES
# =============================================================================

library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(Matrix)
library(scales)
library(patchwork)
library(ggrepel)
library(grid)
library(rstatix)
library(viridis)
library(scCustomize)
library(SeuratExtend)
library(Polychrome)
library(scRepertoire)
library(qs2)
library(enrichR)
library(openxlsx)


# =============================================================================
# STEP 1: PATHS
# =============================================================================

base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
out_dir   <- file.path(base_dir, "Analysis", "Annotation")

out_clusters <- file.path(out_dir, "Cluster_Plot")
out_annot    <- file.path(out_dir, "Annotation_Validation")
out_gating   <- file.path(out_annot, "Cluster1_Gating")
out_avgexpr  <- file.path(out_dir, "AvgExpression")
out_vln_rna  <- file.path(out_annot, "Violins_RNA")
out_vln_adt  <- file.path(out_annot, "Violins_ADT")
out_feat_rna <- file.path(out_annot, "FeaturePlots_RNA")
out_feat_adt <- file.path(out_annot, "FeaturePlots_ADT")
lineage_base <- file.path(out_dir, "Lineage_Exploration")
out_dge      <- file.path(base_dir, "Analysis", "Differential_Expression")
out_pathway  <- file.path(base_dir, "Analysis", "Pathway_Analysis")

for (d in c(out_clusters, out_annot, out_gating, out_avgexpr,
            out_vln_rna, out_vln_adt, out_feat_rna, out_feat_adt,
            lineage_base, out_dge, out_pathway)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# =============================================================================
# STEP 2: LOAD, REJOIN, CLEANUP
# =============================================================================

load(file.path(saved_dir, "TARA_TCR_Combined.RData"))

TARA_ALL[["RNA"]] <- JoinLayers(TARA_ALL[["RNA"]])
TARA_ALL[["ADT"]] <- JoinLayers(TARA_ALL[["ADT"]])

rm(TARA_ALL_TRB_0, TARA_ALL_TRB_1, TARA_ALL_TRB_2)
gc()

TARA_ALL <- subset(TARA_ALL, subset = orig.ident != "CP013_1m")


# =============================================================================
# STEP 3: FILTER SMALL CLUSTERS
# =============================================================================

min_cells <- 100
Idents(TARA_ALL) <- "snn.louvianmlr_1"
tara_keep <- names(which(table(Idents(TARA_ALL)) >= min_cells))
TARA_ALL  <- subset(TARA_ALL, idents = tara_keep)

cat("══ CLUSTERS AFTER FILTERING ══\n")
print(sort(table(Idents(TARA_ALL)), decreasing = TRUE))
cat("\n")


# =============================================================================
# STEP 4: APPLY ALL ANNOTATIONS
# =============================================================================

cat("══ STEP 4: ANNOTATIONS ══\n\n")

final_labels <- c(
  # --- CD4 T cells (3 naive kept separate) ---
  "0"  = "Naive 1 CD4",                   # CD71+CD24+TCRhi
  "2"  = "Naive 2 CD4",                   # CD27hi quiescent
  "41" = "Naive 2 CD4",                   # merge with 2
  "55" = "Naive 2 CD4",                   # merge with 2
  "5"  = "Naive 3 CD4",                   # IL4Rhi
  "7"  = "Transitional Memory CD4",       # CD45RA+RO+ TCF7hi
  "10" = "Th2/Th17 EM CD4",              # LAG3+PD1+CCR4+MAF+
  "20" = "Treg",                          # CD25+FOXP3+
  
  # --- Cluster 12: entirely removed (platelet RNA contamination in both fractions) ---
  "12" = "FLAG_Platelet_cl12",
  
  # --- CD8 T cells ---
  "1"  = "Memory CD8 pending",            # sub-gated in Step 5
  "6"  = "Naive 2 CD8",
  "42" = "Naive 2 CD8",
  "8"  = "TEMRA/CTL CD8",
  "46" = "TEMRA/CTL CD8",
  "52" = "TEMRA/CTL CD8",
  "9"  = "TRDV1+ pending",               # split by TCR in Step 5
  "43" = "TRDV1+ pending",
  "45" = "TRDV1+ pending",
  "50" = "TRDV1+ pending",
  "51" = "TRDV1+ pending",
  "27" = "Tex CD8",
  
  # --- NK cells ---
  "15" = "CD56dim CD122lo NK",
  "36" = "CD56dim CD122lo NK",
  "17" = "CD56dim CD122int NK",
  "3"  = "CD56dim CD122hi NK",
  "40" = "CD56dim CD122hi NK",
  "53" = "CD56dim CD122hi NK",
  "25" = "CD56bright NK",
  
  # --- B cells ---
  "4"  = "Follicular B cell",
  "44" = "Follicular B cell",
  "47" = "Follicular B cell",
  "11" = "Resting Naive B cell",
  "13" = "Transitional B cell",
  "22" = "Activated B cell",
  "35" = "TNF+ B cell",
  
  # --- Myeloid ---
  "14" = "Classical Monocyte",
  "19" = "Non-classical Monocyte",
  "24" = "Intermediate Monocyte",
  "33" = "pDC",
  
  # --- gd T cells ---
  "16" = "Vd1 gd T cell",
  "21" = "Vg9Vd2 gd T cell",
  
  # --- Reclassified former DN T cells ---
  "18" = "ISG+ CD4 T cell",              # CD4+(89%) with IFN/stress signature
  "23" = "CD4dim Naive T cell",           # CD4dim(76%), naive program
  "26" = "CD8+ MAIT cell",
  
  # --- Flagged for removal ---
  "28" = "FLAG_LowQuality",
  "29" = "FLAG_Doublet_T_B",
  "30" = "FLAG_Platelet_NK",
  "31" = "FLAG_Platelet_B",
  "32" = "FLAG_Mixed_Plasmablast",
  "34" = "FLAG_APC_Unclear",
  
  # --- Unresolved ---
  "37" = "T cell (unresolved)",
  "38" = "T cell (unresolved)",
  "39" = "T cell (unresolved)",
  "48" = "T cell (unresolved)",
  "49" = "T cell (unresolved)",
  "54" = "T cell (unresolved)"
)

cluster_ids <- as.character(TARA_ALL$snn.louvianmlr_1)
TARA_ALL$Annotation <- unname(final_labels[cluster_ids])

unmapped <- is.na(TARA_ALL$Annotation)
if (any(unmapped)) {
  TARA_ALL$Annotation[unmapped] <- paste0("Unmapped_cl", cluster_ids[unmapped])
  cat("WARNING: Unmapped clusters:", unique(cluster_ids[unmapped]), "\n")
}


# =============================================================================
# STEP 5: CD8 REFINEMENTS
# =============================================================================

cat("══ STEP 5: CD8 REFINEMENTS ══\n\n")

# --- 5a: Derive has_TCR ---
tcr_col <- intersect(c("ClonalFrequency", "clonalFrequency", "Frequency",
                       "CTstrict", "CTgene", "CTaa"),
                     colnames(TARA_ALL@meta.data))[1]
if (is.na(tcr_col)) stop("No TCR/clonal column found.")
TARA_ALL$has_TCR <- !is.na(TARA_ALL@meta.data[[tcr_col]])

# --- 5b: Cluster 9 split ---
is_trdv1 <- TARA_ALL$Annotation == "TRDV1+ pending"
TARA_ALL$Annotation[is_trdv1 & TARA_ALL$has_TCR == TRUE]  <- "TEMRA/CTL CD8"
TARA_ALL$Annotation[is_trdv1 & TARA_ALL$has_TCR == FALSE] <- "Cytotoxic gd T cell"

# --- 5c: Cluster 1 sub-gate ---
is_cl1    <- TARA_ALL$Annotation == "Memory CD8 pending"
cl1_cells <- colnames(TARA_ALL)[is_cl1]

DefaultAssay(TARA_ALL) <- "ADT"
adt_data   <- GetAssayData(TARA_ALL, slot = "data")
cd45ro_cl1 <- as.numeric(adt_data["CD45RO", cl1_cells])
fas_cl1    <- as.numeric(adt_data["FAS",    cl1_cells])
cd45ra_cl1 <- as.numeric(adt_data["CD45RA", cl1_cells])
cd62l_cl1  <- as.numeric(adt_data["SELL",   cl1_cells])

cd45ro_threshold <- 0
fas_threshold    <- 0

gate <- case_when(
  cd45ro_cl1 > cd45ro_threshold ~ "Naive Intermediate CD8",
  fas_cl1    > fas_threshold    ~ "Tscm CD8",
  TRUE                          ~ "Naive 1 CD8"
)
names(gate) <- cl1_cells
TARA_ALL$Annotation[match(names(gate), colnames(TARA_ALL))] <- gate

cat("Cluster 1 sub-gating:\n"); print(table(gate)); cat("\n")

# Gating diagnostic plots
gate_df <- data.frame(
  CD45RO = cd45ro_cl1, CD45RA = cd45ra_cl1,
  FAS = fas_cl1, CD62L = cd62l_cl1,
  Gate = gate[cl1_cells], stringsAsFactors = FALSE
)
gate_cols <- c("Naive 1 CD8" = "#4E79A7", "Tscm CD8" = "#59A14F",
               "Naive Intermediate CD8" = "#E15759")

p_g1 <- ggplot(gate_df, aes(CD45RA, CD45RO, color = Gate)) +
  geom_point(size = 0.5, alpha = 0.5) + scale_color_manual(values = gate_cols) +
  geom_hline(yintercept = cd45ro_threshold, linetype = "dashed") +
  labs(title = "Cluster 1: CD45RO vs CD45RA") + theme_minimal() +
  theme(legend.position = "bottom")
p_g2 <- ggplot(gate_df, aes(FAS, CD45RO, color = Gate)) +
  geom_point(size = 0.5, alpha = 0.5) + scale_color_manual(values = gate_cols) +
  geom_hline(yintercept = cd45ro_threshold, linetype = "dashed") +
  geom_vline(xintercept = fas_threshold, linetype = "dashed") +
  labs(title = "Cluster 1: CD45RO vs FAS") + theme_minimal() +
  theme(legend.position = "bottom")
p_hist <- (
  ggplot(gate_df, aes(CD45RO, fill = Gate)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = gate_cols) +
    geom_vline(xintercept = cd45ro_threshold, linetype = "dashed") +
    theme_minimal()
) / (
  ggplot(gate_df, aes(FAS, fill = Gate)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = gate_cols) +
    geom_vline(xintercept = fas_threshold, linetype = "dashed") +
    theme_minimal()
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(out_gating, "Gate_CD45RO_vs_CD45RA.png"), p_g1,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_gating, "Gate_CD45RO_vs_FAS.png"), p_g2,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_gating, "Gate_Histograms.png"), p_hist,
       width = 10, height = 10, dpi = 300, bg = "white")


# =============================================================================
# STEP 6: PLATELET FILTER + REMOVE FLAGGED
# =============================================================================

cat("══ STEP 6: PLATELET FILTER & REMOVAL ══\n\n")

DefaultAssay(TARA_ALL) <- "ADT"
cd41 <- GetAssayData(TARA_ALL, slot = "data")["ITGA2B", ]

cat("ITGA2B (CD41) — from data slot:\n")
cat("  Range:", range(as.numeric(cd41)), "\n")
print(quantile(as.numeric(cd41), probs = c(0.9, 0.95, 0.99, 0.995, 0.999, 1.0)))

# DSB threshold — positive signal above ~0.5
platelet_threshold <- 0.5
platelet_cells <- cd41 > platelet_threshold

cat("\nPlatelet cells (ITGA2B >", platelet_threshold, "):", sum(platelet_cells), "\n")
cat("By cluster:\n")
print(table(TARA_ALL$Annotation[platelet_cells]))

# Reference counts at other thresholds
cat("\nReference: threshold=0.25 →", sum(cd41 > 0.25), "cells\n")
cat("Reference: threshold=1.0  →", sum(cd41 > 1.0), "cells\n")

# Histogram
p_plat <- ggplot(data.frame(ITGA2B = as.numeric(cd41[cd41 > 0])),
                 aes(x = ITGA2B)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = platelet_threshold, linetype = "dashed",
             color = "red", linewidth = 1) +
  labs(title = "ITGA2B (CD41) — non-zero cells",
       x = "ITGA2B (DSB)", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_annot, "Platelet_Filter_Histogram.png"),
       p_plat, width = 10, height = 6, dpi = 300, bg = "white")

# Combined removal: flagged clusters + platelet cells
flagged <- grepl("^FLAG_", TARA_ALL$Annotation)
to_remove <- flagged | platelet_cells

cat("\nFlagged cluster cells:", sum(flagged), "\n")
cat("Platelet cells (non-flagged):", sum(platelet_cells & !flagged), "\n")
cat("Total removed:", sum(to_remove), "\n")
cat("Remaining:", sum(!to_remove), "\n\n")

TARA_ALL <- subset(TARA_ALL, cells = colnames(TARA_ALL)[!to_remove])
TARA_ALL$Annotation <- factor(TARA_ALL$Annotation)

cat("══ FINAL ANNOTATION COUNTS ══\n")
print(sort(table(TARA_ALL$Annotation), decreasing = TRUE))
cat("\n")


# =============================================================================
# STEP 7: METADATA CORRECTION
# =============================================================================

cat("══ STEP 7: METADATA CORRECTION ══\n\n")

sample_metadata <- list(
  "CP002_entry" = list(group = "PreART_Entry",        vl = 4284389,  age = 1),
  "CP003_entry" = list(group = "PreART_Entry",        vl = 656769,   age = 1),
  "CP003_V12"   = list(group = "PostART_Unsuppressed", vl = 503497,   age = 12),
  "CP003_V24"   = list(group = "PostART_Unsuppressed", vl = 489676,   age = 24),
  "CP006_entry" = list(group = "PreART_Entry",        vl = 10000000, age = 1),
  "CP006_12m"   = list(group = "PostART_Suppressed",   vl = 73,       age = 12),
  "CP006_V24"   = list(group = "PostART_Suppressed",   vl = 20,       age = 24),
  "CP011_entry" = list(group = "PreART_Entry",        vl = 36965,    age = 2),
  "CP013_entry" = list(group = "PreART_Entry",        vl = 3434,     age = 1),
  "CP013_12m"   = list(group = "PostART_Suppressed",   vl = 20,       age = 12),
  "CP013_24m"   = list(group = "PostART_Suppressed",   vl = 20,       age = 25),
  "CP016_entry" = list(group = "PreART_Entry",        vl = 1978332,  age = 2),
  "CP017_entry" = list(group = "PreART_Entry",        vl = 3167384,  age = 2),
  "CP018_entry" = list(group = "PreART_Entry",        vl = 176970,   age = 1),
  "CP018_V24"   = list(group = "PostART_Suppressed",   vl = 20,       age = 25),
  "CP018_42m"   = list(group = "PostART_Suppressed",   vl = 20,       age = 42),
  "CP020_V1"    = list(group = "PreART_Entry",        vl = 414409,   age = 1),
  "CP020_V12"   = list(group = "PostART_Unsuppressed", vl = 6293122,  age = 12),
  "CP020_V44"   = list(group = "PostART_Unsuppressed", vl = 534772,   age = 44),
  "CP022_entry" = list(group = "PreART_Entry",        vl = 5075764,  age = 1),
  "CP042_entry" = list(group = "PreART_Entry",        vl = 6113,     age = 1)
)

if (!"Timepoint_Group" %in% colnames(TARA_ALL@meta.data)) TARA_ALL$Timepoint_Group <- NA_character_
if (!"Viral_Load_num" %in% colnames(TARA_ALL@meta.data))  TARA_ALL$Viral_Load_num  <- NA_real_
if (!"Age" %in% colnames(TARA_ALL@meta.data))             TARA_ALL$Age             <- NA_real_

for (sname in names(sample_metadata)) {
  mask <- TARA_ALL$orig.ident == sname
  if (sum(mask) == 0) next
  meta <- sample_metadata[[sname]]
  TARA_ALL$Timepoint_Group[mask] <- meta$group
  TARA_ALL$Viral_Load_num[mask]  <- meta$vl
  TARA_ALL$Age[mask]             <- meta$age
}

TARA_ALL$orig.ident_raw <- TARA_ALL$orig.ident
ident_rename <- sapply(names(sample_metadata), function(s) {
  paste0(sub("_.*$", "", s), "_", sample_metadata[[s]]$age, "m")
})
names(ident_rename) <- names(sample_metadata)
new_ids <- ident_rename[TARA_ALL$orig.ident]
new_ids[is.na(new_ids)] <- TARA_ALL$orig.ident[is.na(new_ids)]
TARA_ALL$orig.ident <- unname(new_ids)

# Derive Condition column (HEI / HEU / HUU) if not present
if (!"Condition" %in% colnames(TARA_ALL@meta.data)) {
  # HEI samples have Timepoint_Group set; others are HEU or HUU
  # Attempt to derive from existing metadata
  TARA_ALL$Condition <- ifelse(!is.na(TARA_ALL$Timepoint_Group), "HEI",
                               TARA_ALL$Condition)
  cat("NOTE: Condition column derived. Verify HEU/HUU labels are correct.\n")
}
cat("Condition distribution:\n")
print(table(TARA_ALL$Condition, useNA = "ifany"))


# =============================================================================
# STEP 8: HELPERS
# =============================================================================

scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x)) / rng
  }))
  colnames(scaled) <- colnames(mat)
  scaled
}

export_avg_expression <- function(obj, group_by, rna_feats, adt_feats, out_d) {
  dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
  
  # RNA — AverageExpression is fine for log-normalized RNA
  DefaultAssay(obj) <- "RNA"
  rna_avail <- intersect(rna_feats, rownames(obj[["RNA"]]))
  avg_rna <- AverageExpression(obj, assays = "RNA", features = rna_avail,
                               group.by = group_by, slot = "data")$RNA
  colnames(avg_rna) <- gsub("^g ", "", colnames(avg_rna))
  write.csv(avg_rna, file.path(out_d, "AvgExpr_RNA.csv"))
  write.csv(round(scale_01(as.matrix(avg_rna)), 4), file.path(out_d, "AvgExpr_RNA_scaled.csv"))
  avg_rna_all <- AverageExpression(obj, assays = "RNA", group.by = group_by, slot = "data")$RNA
  colnames(avg_rna_all) <- gsub("^g ", "", colnames(avg_rna_all))
  write.csv(avg_rna_all, file.path(out_d, "AvgExpr_RNA_AllGenes.csv"))
  
  # ADT — manual averaging to avoid expm1 back-transformation
  DefaultAssay(obj) <- "ADT"
  adt_avail <- intersect(adt_feats, rownames(obj[["ADT"]]))
  adt_mat   <- GetAssayData(obj, slot = "data")[adt_avail, , drop = FALSE]
  clusters  <- obj@meta.data[[group_by]]
  avg_adt <- sapply(levels(factor(clusters)), function(cl) {
    cells <- which(clusters == cl)
    if (length(cells) == 0) return(rep(0, length(adt_avail)))
    rowMeans(adt_mat[, cells, drop = FALSE])
  })
  rownames(avg_adt) <- adt_avail
  write.csv(avg_adt, file.path(out_d, "AvgExpr_ADT.csv"))
  write.csv(round(scale_01(as.matrix(avg_adt)), 4), file.path(out_d, "AvgExpr_ADT_scaled.csv"))
  
  # Percent expression
  DefaultAssay(obj) <- "RNA"
  rna_mat <- GetAssayData(obj, slot = "data")[rna_avail, , drop = FALSE]
  pct_rna <- as.data.frame(t(sapply(levels(factor(clusters)), function(cl) {
    cells <- which(obj@meta.data[[group_by]] == cl)
    rowMeans(rna_mat[, cells, drop = FALSE] > 0) * 100
  })))
  write.csv(round(pct_rna, 2), file.path(out_d, "PctExpr_RNA.csv"))
  
  pct_adt <- as.data.frame(t(sapply(levels(factor(clusters)), function(cl) {
    cells <- which(obj@meta.data[[group_by]] == cl)
    rowMeans(adt_mat[, cells, drop = FALSE] > 0) * 100
  })))
  write.csv(round(pct_adt, 2), file.path(out_d, "PctExpr_ADT.csv"))
}

generate_violin_plots <- function(obj, features, assay, out_d, group_by) {
  dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
  DefaultAssay(obj) <- assay
  avail <- intersect(features, rownames(obj[[assay]]))
  for (feat in avail) {
    p <- VlnPlot2(obj, features = feat, group.by = group_by,
                  assay = assay, cols = "light",
                  show.mean = TRUE, mean_colors = c("red", "blue"))
    safe <- gsub("[^A-Za-z0-9._-]", "_", feat)
    ggsave(file.path(out_d, paste0("Vln_", safe, ".png")),
           p, width = 14, height = 5, dpi = 300, bg = "white")
  }
}

generate_feature_plots <- function(obj, features, assay, out_d) {
  dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
  DefaultAssay(obj) <- assay
  avail <- intersect(features, rownames(obj[[assay]]))
  for (feat in avail) {
    p <- DimPlot2(obj, features = feat, reduction = "wnn.umap") +
      ggtitle(paste(assay, "|", feat))
    safe <- gsub("[^A-Za-z0-9._-]", "_", feat)
    ggsave(file.path(out_d, paste0("Feat_", safe, ".png")),
           p, dpi = 500, width = 8, height = 6, bg = "white")
  }
}

rna_markers_flat <- unique(c(
  "CD8A","CD8B","CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2",
  "CCR7","SELL","TCF7","LEF1","BACH2","IL7R","BCL2","S1PR1","KLF2","FOXP1","SATB1",
  "GZMB","GZMA","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","KLRG1","CX3CR1","FGFBP2","FCGR3A",
  "TBX21","EOMES","RUNX3","PRDM1","ZEB2","ID2",
  "TOX","TOX2","PDCD1","HAVCR2","TIGIT","LAG3","CTLA4","ENTPD1","CD160","CD244",
  "CD27","CD28","BCL6","CXCR5",
  "MKI67","HLA-DRA","HLA-DRB1","CD38","FAS","ICOS","TNFRSF9","CXCR6","CD69",
  "ITGAE","ITGA1","ZNF683",
  "TRDV1","TRDV2","TRGV9","TRDC","TRGC1","TRGC2",
  "TRAV1-2","SLC4A10","KLRB1","ZBTB16","RORC",
  "TYROBP","KLRD1","KIR2DL3","KIR3DL1","NCAM1","KLRF1","KLRC1","KLRC2",
  "ISG15","IFIT1","IFIT2","IFIT3","IFI44L","MX1","OAS1","STAT1","IRF7",
  "HSPA1A","HSPA1B","HSP90AA1","DNAJB1","IFI27",
  "CCL3","CCL4","CCL4L2","CCL5","XCL1","XCL2","IFNG","TNF","IL2","CSF2",
  "B3GAT1",
  "CD4","FOXP3","IL2RA","IKZF2","CCR4","GATA3","MAF","BATF",
  "CD19","MS4A1","CD79A","CD79B","PAX5","BANK1","BLK","BLNK","TCL1A","IGHM","IGHD","IGKC","IGLC2",
  "SDC1","XBP1","MZB1","JCHAIN","IGHG1","IGHA1","IRF4",
  "NCR1","NCR3","FCER1G","SH2D1B","SPON2","CLIC3",
  "CD14","LYZ","S100A8","S100A9","S100A12","CSF1R","CD68","MARCO","VCAN","FCN1","MNDA",
  "FCER1A","CLEC10A","CD1C","CLEC9A","XCR1","BATF3","IRF8","LAMP3",
  "LILRA4","CLEC4C","IL3RA","TCF4",
  "PPBP","PF4","GP9","ITGA2B","TUBB1","TREML1",
  "HBB","HBA1","HBA2","ALAS2","SLC4A1","GYPA"
))

DefaultAssay(TARA_ALL) <- "ADT"
all_adt_markers <- sort(rownames(TARA_ALL[["ADT"]]))

# Publication-quality UMAP helper
plot_pub_umap <- function(obj, group_col, color_map, label_size = 12,
                          pt_size = 1.0, pt_alpha = 0.7, colored_labels = TRUE) {
  umap_df <- as.data.frame(Embeddings(obj, reduction = "wnn.umap"))
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Group <- obj@meta.data[[group_col]]
  
  centroids <- umap_df %>%
    group_by(Group) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")
  
  x_arrow <- min(umap_df$UMAP1, na.rm = TRUE) + 1
  y_arrow <- min(umap_df$UMAP2, na.rm = TRUE) + 1
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point(size = pt_size, alpha = pt_alpha, stroke = 0) +
    scale_color_manual(values = color_map)
  
  if (colored_labels) {
    p <- p + geom_label_repel(
      data = centroids,
      aes(x = UMAP1, y = UMAP2, label = Group, color = Group),
      fill = "white", size = label_size, fontface = "bold",
      label.size = 0.6, label.padding = unit(0.3, "lines"),
      box.padding = 1.0, point.padding = 0.5,
      segment.color = "grey40", segment.size = 0.6,
      min.segment.length = 0, max.overlaps = 40, show.legend = FALSE)
  } else {
    p <- p + geom_label_repel(
      data = centroids,
      aes(x = UMAP1, y = UMAP2, label = Group),
      fill = "white", color = "black", size = label_size, fontface = "bold",
      label.size = 0.6, label.padding = unit(0.3, "lines"),
      box.padding = 1.0, point.padding = 0.5,
      segment.color = "grey40", segment.size = 0.6,
      min.segment.length = 0, max.overlaps = 40, show.legend = FALSE)
  }
  
  p + annotate("segment", x = x_arrow, xend = x_arrow + 3.5,
               y = y_arrow, yend = y_arrow,
               arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
               linewidth = 1, color = "black") +
    annotate("text", x = x_arrow + 1.75, y = y_arrow - 1.2,
             label = "UMAP1", size = 6, fontface = "bold") +
    annotate("segment", x = x_arrow, xend = x_arrow,
             y = y_arrow, yend = y_arrow + 3.5,
             arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
             linewidth = 1, color = "black") +
    annotate("text", x = x_arrow - 1.2, y = y_arrow + 1.75,
             label = "UMAP2", size = 6, fontface = "bold", angle = 90) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10)) +
    coord_fixed()
}


################################################################################
#
#  ██████╗██╗  ██╗███████╗ ██████╗██╗  ██╗██████╗  ██████╗ ██╗███╗   ██╗████████╗
# ██╔════╝██║  ██║██╔════╝██╔════╝██║ ██╔╝██╔══██╗██╔═══██╗██║████╗  ██║╚══██╔══╝
# ██║     ███████║█████╗  ██║     █████╔╝ ██████╔╝██║   ██║██║██╔██╗ ██║   ██║
# ██║     ██╔══██║██╔══╝  ██║     ██╔═██╗ ██╔═══╝ ██║   ██║██║██║╚██╗██║   ██║
# ╚██████╗██║  ██║███████╗╚██████╗██║  ██╗██║     ╚██████╔╝██║██║ ╚████║   ██║
#  ╚═════╝╚═╝  ╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝      ╚═════╝ ╚═╝╚═╝  ╚═══╝   ╚═╝
#
# PAUSE HERE. Export CSVs, review, then continue from STEP 10.
#
################################################################################

# =============================================================================
# STEP 9: EXPORT CSVs FOR CONFIRMATION
# =============================================================================

export_avg_expression(TARA_ALL, "Annotation", rna_markers_flat, all_adt_markers, out_avgexpr)

# Quick checkpoint UMAP using publication helper
annot_lvl <- sort(unique(as.character(TARA_ALL$Annotation)))
n_a <- length(annot_lvl)
chk_pal <- unname(createPalette(n_a, c("#FF0000","#00FF00","#0000FF"), M = 1000))
names(chk_pal) <- annot_lvl

p_check <- plot_pub_umap(TARA_ALL, "Annotation", chk_pal,
                         label_size = 7, pt_size = 0.8, pt_alpha = 0.6)

ggsave(file.path(out_clusters, "TARA_Annotation_Checkpoint.png"),
       p_check, width = 16, height = 14, dpi = 300, bg = "white")

qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_checkpoint.qs2"))

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " CHECKPOINT: CSVs exported.\n",
        " Files: ", out_avgexpr, "\n",
        " Review, then continue from STEP 10.\n",
        "══════════════════════════════════════════════════════════════\n"
)

# ▼▼▼ UNCOMMENT TO STOP HERE ▼▼▼
# stop("CHECKPOINT — review CSVs before continuing")

# To resume: TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_checkpoint.qs2"))


# =============================================================================
# STEP 10: FINAL UMAPs (publication quality)
# =============================================================================

cat("══ STEP 10: FINAL UMAPs ══\n\n")

# --- 10a: Broad lineage UMAP ---
broad_labels <- case_when(
  TARA_ALL$Annotation %in% c("Naive 1 CD4", "Naive 2 CD4", "Naive 3 CD4",
                             "Transitional Memory CD4", "Th2/Th17 EM CD4", "Treg",
                             "ISG+ CD4 T cell", "CD4dim Naive T cell") ~ "CD4 T cells",
  TARA_ALL$Annotation %in% c("Naive 1 CD8", "Tscm CD8", "Naive Intermediate CD8",
                             "Naive 2 CD8", "TEMRA/CTL CD8", "Tex CD8") ~ "CD8 T cells",
  TARA_ALL$Annotation %in% c("Vd1 gd T cell", "Vg9Vd2 gd T cell",
                             "Cytotoxic gd T cell", "CD8+ MAIT cell") ~ "\u03b3\u03b4 / MAIT",
  TARA_ALL$Annotation %in% c("CD56dim CD122lo NK", "CD56dim CD122int NK",
                             "CD56dim CD122hi NK", "CD56bright NK") ~ "NK cells",
  TARA_ALL$Annotation %in% c("Follicular B cell", "Resting Naive B cell",
                             "Transitional B cell", "Activated B cell", "TNF+ B cell") ~ "B cells",
  TARA_ALL$Annotation %in% c("Classical Monocyte", "Non-classical Monocyte",
                             "Intermediate Monocyte") ~ "Monocytes",
  TARA_ALL$Annotation == "pDC" ~ "pDC",
  TRUE ~ "Other"
)
TARA_ALL$Broad_CellType <- factor(broad_labels)

broad_cols <- c(
  "CD4 T cells"  = "#1B9E77",
  "CD8 T cells"  = "#D95F02",
  "\u03b3\u03b4 / MAIT" = "#7570B3",
  "NK cells"     = "#66A61E",
  "B cells"      = "#E6AB02",
  "Monocytes"    = "#666666",
  "pDC"          = "#1F78B4",
  "Other"        = "#CCCCCC"
)

p_broad <- plot_pub_umap(TARA_ALL, "Broad_CellType", broad_cols,
                         label_size = 14, pt_size = 1.2)
ggsave(file.path(out_clusters, "TARA_Broad_UMAP.png"),
       p_broad, width = 12, height = 11, dpi = 300, bg = "white")

# --- 10c: Detailed annotation UMAP ---
annot_levels <- sort(unique(as.character(TARA_ALL$Annotation)))
n_annot <- length(annot_levels)
if (n_annot <= 36) {
  annot_pal <- unname(createPalette(n_annot, c("#FF0000","#00FF00","#0000FF"),
                                    M = 1000))
} else {
  annot_pal <- rep(hue_pal()(min(n_annot, 30)), length.out = n_annot)
}
names(annot_pal) <- annot_levels

p_detail <- plot_pub_umap(TARA_ALL, "Annotation", annot_pal,
                          label_size = 8, pt_size = 0.9, pt_alpha = 0.7)
ggsave(file.path(out_clusters, "TARA_Detailed_UMAP.png"),
       p_detail, width = 16, height = 14, dpi = 300, bg = "white")

# --- 10d: Black-label version of broad UMAP ---
p_broad_bk <- plot_pub_umap(TARA_ALL, "Broad_CellType", broad_cols,
                            label_size = 14, pt_size = 1.2, colored_labels = FALSE)
ggsave(file.path(out_clusters, "TARA_Broad_UMAP_blacklabels.png"),
       p_broad_bk, width = 12, height = 11, dpi = 300, bg = "white")

# --- 10e: Composition bar plot ---
p_comp <- ClusterDistrBar(
  TARA_ALL$orig.ident, TARA_ALL$Annotation,
  cols = "default", flip = FALSE, border = "black"
) + theme(axis.title.x = element_blank())
ggsave(file.path(out_clusters, "TARA_Cluster_Distribution.png"),
       p_comp, width = 25, height = 14, dpi = 300, bg = "white")

cat("✓ Publication UMAPs + composition plot saved\n\n")


# =============================================================================
# STEP 11: VIOLIN & FEATURE PLOTS (all clusters, single set)
# =============================================================================

cat("\n══ STEP 11: VIOLIN & FEATURE PLOTS ══\n\n")

DefaultAssay(TARA_ALL) <- "RNA"
rel_rna <- intersect(rna_markers_flat, rownames(TARA_ALL[["RNA"]]))

generate_violin_plots(TARA_ALL, rel_rna, "RNA", out_vln_rna, "Annotation")
generate_violin_plots(TARA_ALL, all_adt_markers, "ADT", out_vln_adt, "Annotation")
generate_feature_plots(TARA_ALL, rel_rna, "RNA", out_feat_rna)
generate_feature_plots(TARA_ALL, all_adt_markers, "ADT", out_feat_adt)

# DotPlots — all clusters
DefaultAssay(TARA_ALL) <- "RNA"
p_dot_rna <- DotPlot(TARA_ALL, features = rel_rna, group.by = "Annotation",
                     cols = c("lightgrey", "#D73027"), dot.scale = 5) +
  RotatedAxis() + ggtitle("All Clusters — RNA") +
  theme(axis.text.x = element_text(size = 6))
ggsave(file.path(out_annot, "DotPlot_RNA_AllClusters.png"),
       p_dot_rna, width = 28, height = 12, dpi = 300, bg = "white")

DefaultAssay(TARA_ALL) <- "ADT"
p_dot_adt <- DotPlot(TARA_ALL, features = all_adt_markers, group.by = "Annotation",
                     cols = c("lightgrey", "#08519C"), dot.scale = 4) +
  RotatedAxis() + ggtitle("All Clusters — ADT") +
  theme(axis.text.x = element_text(size = 5))
ggsave(file.path(out_annot, "DotPlot_ADT_AllClusters.png"),
       p_dot_adt, width = 38, height = 12, dpi = 300, bg = "white")

cat("✓ Violin, feature, and dot plots complete\n\n")


# =============================================================================
# STEP 12: DIFFERENTIAL GENE EXPRESSION
# =============================================================================

cat("══ STEP 12: DIFFERENTIAL GENE EXPRESSION ══\n\n")

dir.create(out_dge, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(TARA_ALL) <- "RNA"
Idents(TARA_ALL) <- "Annotation"
TARA_ALL$PID <- sub("_.*", "", TARA_ALL$orig.ident_raw, perl = TRUE)

# --- DGE helper ---
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
count_sig <- function(df, pcol = "p_val_adj", thr = 0.05) sum(df[[pcol]] < thr, na.rm = TRUE)

run_clusterwise_de <- function(obj, group_col, ident1, ident2, out_dir,
                               latent = c("nCount_RNA"),
                               min_cells_per_grp = 10,
                               title_prefix = "", file_stub = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  Idents(obj) <- "Annotation"
  DefaultAssay(obj) <- "RNA"
  
  de_list  <- list()
  count_tbl <- data.frame()
  cl_levels <- levels(droplevels(factor(obj$Annotation)))
  
  for (cl in cl_levels) {
    message("  DE: ", cl)
    cl_cells <- colnames(obj)[obj$Annotation == cl]
    sc <- subset(obj, cells = cl_cells)
    sc$.__grp <- sc@meta.data[[group_col]]
    
    tab <- table(sc$.__grp)
    if (!all(c(ident1, ident2) %in% names(tab))) next
    if (any(tab[c(ident1, ident2)] < min_cells_per_grp)) next
    
    # Check latent vars have >1 level
    use_latent <- latent
    for (lv in latent) {
      if (lv %in% colnames(sc@meta.data)) {
        if (length(unique(sc@meta.data[[lv]][sc$.__grp %in% c(ident1, ident2)])) < 2) {
          use_latent <- setdiff(use_latent, lv)
        }
      }
    }
    
    de <- tryCatch(
      FindMarkers(sc, ident.1 = ident1, ident.2 = ident2,
                  group.by = ".__grp", test.use = "MAST",
                  latent.vars = use_latent),
      error = function(e) { message("    Failed: ", e$message); NULL })
    if (is.null(de)) next
    
    comp_name <- if (is.null(file_stub)) paste0(ident1, "_vs_", ident2) else file_stub
    write.csv(de, file.path(out_dir, paste0(safe_name(cl), "_", comp_name, ".csv")))
    de_list[[cl]] <- de
    count_tbl <- rbind(count_tbl, data.frame(Cluster = cl, DE_Genes = count_sig(de)))
  }
  
  # Summary bar plot
  if (nrow(count_tbl) > 0) {
    p <- ggplot(count_tbl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      xlab("Cluster") + ylab("# DE Genes (p_adj < 0.05)") +
      ggtitle(paste0(title_prefix, ": ", ident1, " vs ", ident2))
    ggsave(file.path(out_dir, paste0("Summary_DE_Counts_", safe_name(comp_name), ".png")),
           p, width = 12, height = 6, dpi = 300, bg = "white")
  }
  
  invisible(list(de = de_list, counts = count_tbl))
}

# --- Comparison 1: HEI vs HEU (Pre-ART, Age ≤ 2) ---
cat("── HEI vs HEU (Pre-ART) ──\n")
TARA_hei_heu <- subset(TARA_ALL, subset = Age <= 2 & Condition %in% c("HEI", "HEU"))

run_clusterwise_de(
  obj = TARA_hei_heu, group_col = "Condition",
  ident1 = "HEI", ident2 = "HEU",
  out_dir = file.path(out_dge, "HEIvsHEU_PreART"),
  latent = c("nCount_RNA"),
  title_prefix = "Pre-ART (Age ≤ 2m)"
)

# --- Comparison 2: HEU vs HUU ---
cat("── HEU vs HUU ──\n")
TARA_heu_huu <- subset(TARA_ALL, subset = Condition %in% c("HEU", "HUU"))

run_clusterwise_de(
  obj = TARA_heu_huu, group_col = "Condition",
  ident1 = "HEU", ident2 = "HUU",
  out_dir = file.path(out_dge, "HEUvsHUU"),
  latent = c("nCount_RNA"),
  title_prefix = "TARA"
)

# --- Comparisons 3–5: Longitudinal HEI (control for PID) ---
cat("── Longitudinal HEI comparisons ──\n")
TARA_HEI <- subset(TARA_ALL, subset = !is.na(Timepoint_Group) &
                     Timepoint_Group %in% c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))

# 3: Suppressed vs PreART
run_clusterwise_de(
  obj = subset(TARA_HEI, subset = Timepoint_Group %in% c("PostART_Suppressed", "PreART_Entry")),
  group_col = "Timepoint_Group",
  ident1 = "PostART_Suppressed", ident2 = "PreART_Entry",
  out_dir = file.path(out_dge, "Suppressed_vs_PreART"),
  latent = c("nCount_RNA", "PID"),
  title_prefix = "HEI",
  file_stub = "Suppressed_vs_PreART"
)

# 4: Unsuppressed vs PreART
run_clusterwise_de(
  obj = subset(TARA_HEI, subset = Timepoint_Group %in% c("PostART_Unsuppressed", "PreART_Entry")),
  group_col = "Timepoint_Group",
  ident1 = "PostART_Unsuppressed", ident2 = "PreART_Entry",
  out_dir = file.path(out_dge, "Unsuppressed_vs_PreART"),
  latent = c("nCount_RNA", "PID"),
  title_prefix = "HEI",
  file_stub = "Unsuppressed_vs_PreART"
)

# 5: Suppressed vs Unsuppressed
run_clusterwise_de(
  obj = subset(TARA_HEI, subset = Timepoint_Group %in% c("PostART_Suppressed", "PostART_Unsuppressed")),
  group_col = "Timepoint_Group",
  ident1 = "PostART_Suppressed", ident2 = "PostART_Unsuppressed",
  out_dir = file.path(out_dge, "Suppressed_vs_Unsuppressed"),
  latent = c("nCount_RNA", "PID"),
  title_prefix = "HEI",
  file_stub = "Suppressed_vs_Unsuppressed"
)

cat("✓ DGE complete\n\n")


# =============================================================================
# STEP 13: PATHWAY ANALYSIS (EnrichR + lollipop plots)
# =============================================================================

cat("══ STEP 13: PATHWAY ANALYSIS ══\n\n")

dir.create(out_pathway, recursive = TRUE, showWarnings = FALSE)

enrichr_dbs <- c(
  "TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs",
  "KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023",
  "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019"
)
tf_dbs <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_dbs <- setdiff(enrichr_dbs, tf_dbs)

run_enrichment <- function(gene_df, label, out_base) {
  gene_list <- rownames(gene_df)
  if (length(gene_list) == 0) return(NULL)
  
  enrichment <- tryCatch(enrichr(gene_list, enrichr_dbs),
                         error = function(e) { message("  Enrichr failed: ", e$message); NULL })
  if (is.null(enrichment)) return(NULL)
  
  dir.create(file.path(out_base, "CSVs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_base, "Plots"), recursive = TRUE, showWarnings = FALSE)
  
  # Save Excel
  wb <- createWorkbook()
  for (db in names(enrichment)) {
    addWorksheet(wb, substr(db, 1, 31))
    writeData(wb, substr(db, 1, 31), enrichment[[db]])
  }
  saveWorkbook(wb, file.path(out_base, "CSVs", paste0(label, "_Enrichment.xlsx")), overwrite = TRUE)
  
  # Collect top hits
  collect_top <- function(db_names, n_top = 10) {
    out <- list()
    for (db_name in db_names) {
      res <- enrichment[[db_name]]
      if (is.null(res) || nrow(res) == 0) next
      if ("Combined Score" %in% colnames(res)) res <- rename(res, Combined.Score = `Combined Score`)
      if (!"Combined.Score" %in% colnames(res)) next
      sig <- filter(res, Adjusted.P.value < 0.05)
      if (nrow(sig) > 0)
        out[[db_name]] <- sig %>% arrange(desc(Combined.Score)) %>%
          slice_head(n = n_top) %>% mutate(Database = db_name)
    }
    bind_rows(out)
  }
  
  # Lollipop plot helper
  make_lollipop <- function(df, title, filename) {
    if (nrow(df) == 0) return(NULL)
    df <- df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 20)
    p <- ggplot(df, aes(x = reorder(Term, Combined.Score), y = Combined.Score, color = Database)) +
      geom_segment(aes(xend = reorder(Term, Combined.Score), y = 0, yend = Combined.Score),
                   linewidth = 0.8) +
      geom_point(aes(size = -log10(Adjusted.P.value)), alpha = 0.85) +
      scale_size_continuous(name = "-log10(adj.p)", range = c(2, 7)) +
      coord_flip() +
      labs(title = title, x = NULL, y = "Combined Score") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 10, color = "black"),
            plot.title = element_text(size = 14, hjust = 0.5))
    ggsave(file.path(out_base, "Plots", filename),
           p, width = 13, height = 10, dpi = 300, bg = "white")
  }
  
  tf_df      <- collect_top(tf_dbs)
  pathway_df <- collect_top(pathway_dbs)
  
  make_lollipop(tf_df, paste("Top TFs —", label),
                paste0(label, "_Transcription_Factors.png"))
  make_lollipop(pathway_df, paste("Top Pathways —", label),
                paste0(label, "_Pathways.png"))
}

# Loop over DGE output directories
dge_dirs <- list.dirs(out_dge, recursive = FALSE)

for (dge_d in dge_dirs) {
  comparison <- basename(dge_d)
  cat(sprintf("── Enrichment: %s ──\n", comparison))
  comp_out <- file.path(out_pathway, comparison)
  
  csv_files <- list.files(dge_d, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!grepl("Summary_", basename(csv_files))]  # skip summary PNGs
  
  for (csvf in csv_files) {
    dge <- tryCatch(read.csv(csvf, row.names = 1, check.names = FALSE),
                    error = function(e) { message("  Read failed: ", basename(csvf)); NULL })
    if (is.null(dge) || !"p_val_adj" %in% colnames(dge)) next
    
    sig <- dge %>% filter(!is.na(p_val_adj) & p_val_adj < 0.05)
    if (nrow(sig) == 0) next
    
    label <- tools::file_path_sans_ext(basename(csvf))
    message("  Enrichr: ", label)
    run_enrichment(sig, label, comp_out)
  }
}

cat("✓ Pathway analysis complete\n\n")


# =============================================================================
# STEP 14: SAVE FINAL OBJECTS
# =============================================================================

cat("══ STEP 14: SAVING ══\n\n")

qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_annotated_final.qs2"))
cat("✓ TARA_ALL_annotated_final.qs2\n")

# CD8 subset
TARA_cd8 <- subset(TARA_ALL, subset = Annotation %in% c(
  "Naive 1 CD8", "Tscm CD8", "Naive Intermediate CD8", "Naive 2 CD8",
  "TEMRA/CTL CD8", "Tex CD8"))
qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_final.qs2"))
cat("✓ TARA_cd8_final.qs2\n")

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " PIPELINE COMPLETE\n",
        "══════════════════════════════════════════════════════════════\n",
        " Annotation: ", out_dir, "\n",
        " DGE:        ", out_dge, "\n",
        " Pathways:   ", out_pathway, "\n",
        " Objects:    ", saved_dir, "\n",
        "══════════════════════════════════════════════════════════════\n"
)