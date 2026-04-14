################################################################################
# TARA — Unified WNN Annotation Pipeline
#
# All annotations (CD8 refined + non-CD8 finalized) are built into this script.
# No Excel file needed. Run end-to-end.
#
# OUTPUT:
#   ~/Documents/CD8_Longitudinal/Analysis/Annotation/
#     Cluster_Plot/                  — UMAPs (raw clusters + final annotation)
#     Annotation_Validation/         — Gating plots, DE tables
#     AvgExpression/                 — CSVs for all final clusters
#     Lineage_Exploration/<name>/    — Per-lineage plots & tables
#   ~/Documents/CD8_Longitudinal/saved_R_data/
#     TARA_ALL_annotated_final.qs2
#     TARA_<lineage>_subset.qs2
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

for (d in c(out_clusters, out_annot, out_gating, out_avgexpr,
            out_vln_rna, out_vln_adt, out_feat_rna, out_feat_adt,
            lineage_base)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# =============================================================================
# STEP 2: LOAD DATA, REJOIN LAYERS, CLEANUP
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

# Raw cluster UMAP (reference)
p_raw <- DimPlot_scCustom(
  TARA_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_1",
  label = TRUE, repel = TRUE, label.size = 5
) + ggtitle("TARA: Raw Cluster Numbers")
ggsave(file.path(out_clusters, "TARA_WNN_raw_clusters.png"),
       p_raw, width = 10, height = 8, dpi = 300)


# =============================================================================
# STEP 4: APPLY ALL ANNOTATIONS
#
# Rough labels from manual review, including cluster merges ("combine with X"),
# then overwrite with finalized non-CD8 renames/merges/reclassifications.
# CD8 refinements (sub-gating, splitting) applied in Step 5.
# =============================================================================

cat("══ STEP 4: APPLYING ANNOTATIONS ══\n\n")

# --- 4a: Map raw cluster numbers to FINAL labels ----------------------------
# Clusters that merge share the same final label.
# Suspicious clusters (Plasmablast, APC, DN4, DN5) flagged for removal.
final_labels <- c(
  # --- CD4 T cells ---
  "0"  = "Naive CD4",                     # was "CD4 T cell"
  "2"  = "Naive CD4",                     # was "Naive CD4 T cell_1" — merged
  "5"  = "Naive CD4",                     # was "Naive CD4 T cell_2" — merged
  "41" = "Naive CD4",                     # combine with 2
  "55" = "Naive CD4",                     # combine with 2
  "7"  = "Central Memory CD4",            # was "Memory CD4 T cell"
  "10" = "Th2/Th17 Effector Memory CD4",  # was "Memory CD4 T helper"
  "12" = "DNAM1+ CD4/CD8 mixed",          # split by CD8A in Step 5
  "20" = "Treg",                          # was "T reg"
  
  # --- CD8 T cells (refined in Step 5) ---
  "1"  = "Memory CD8 pending",            # sub-gated in Step 5
  "6"  = "Naive 2 CD8",
  "8"  = "TEMRA/CTL CD8",
  "46" = "TEMRA/CTL CD8",                 # combine with 8
  "52" = "TEMRA/CTL CD8",                 # combine with 8
  "9"  = "TRDV1+ pending",               # split by TCR in Step 5
  "43" = "TRDV1+ pending",               # combine with 9
  "45" = "TRDV1+ pending",
  "50" = "TRDV1+ pending",
  "51" = "TRDV1+ pending",
  "42" = "Naive 2 CD8",                   # combine with 6
  "27" = "Tex CD8",
  
  # --- NK cells ---
  "3"  = "CD56dim NK (CD122hi)",          # was "IL2RB+ NK cell"
  "40" = "CD56dim NK (CD122hi)",          # combine with 3
  "53" = "CD56dim NK (CD122hi)",          # combine with 3
  "15" = "CD56dim NK",                    # was "NK cell_1" — merged
  "36" = "CD56dim NK",                    # combine with 15
  "17" = "CD56dim NK",                    # was "NK cell_2" — merged with 15
  "25" = "CD56bright NK",
  
  # --- B cells ---
  "4"  = "Mature Naive B",                # was "IgD+IgM+ B cell" — merged
  "44" = "Mature Naive B",                # combine with 4
  "47" = "Mature Naive B",                # combine with 4
  "11" = "Mature Naive B",                # was "Naive B cell" — merged
  "31" = "Mature Naive B",                # was "B cell" — merged
  "13" = "Transitional B",                # was "Transitional B cell"
  "22" = "Activated B",                   # was "CD69+ B cell"
  "35" = "TNF+ B cell",
  
  # --- Myeloid ---
  "14" = "Classical Monocyte",            # was "CD14hi CD16+ Monocyte"
  "19" = "Non-classical Monocyte",        # was "CD14+ CD16hi Monocyte"
  "24" = "Intermediate Monocyte",         # was "CD14lo CD16lo Monocyte"
  "33" = "pDC",
  
  # --- γδ T cells ---
  "16" = "Vd1 gd T cell",                # was "Gamma Delta 1 T cells"
  "21" = "Vg9Vd2 gd T cell",             # was "Gamma Delta 2 T cells"
  
  # --- DN T cells (reclassified) ---
  "18" = "CD4dim Naive T cell",           # was "DN T cell_1"
  "23" = "Naive-like T cell",             # was "DN T cell_2"
  "26" = "MAIT / Innate-like T cell",     # was "DN T cell_3"
  "30" = "NKT cell",                      # was "DN T cell_6"
  
  # --- Flagged for removal ---
  "28" = "FLAG_LowQuality",              # DN T cell_4 — near-zero expression
  "29" = "FLAG_Doublet_T_B",             # DN T cell_5 — T+B doublet
  "32" = "FLAG_Mixed_Plasmablast",       # Plasmablast — NK/myeloid contamination
  "34" = "FLAG_APC_Unclear",             # APC — ambiguous identity
  
  # --- Generic T cell clusters (small, from merges) ---
  "37" = "T cell (unresolved)",
  "38" = "T cell (unresolved)",
  "39" = "T cell (unresolved)",
  "48" = "T cell (unresolved)",
  "49" = "T cell (unresolved)",
  "54" = "T cell (unresolved)"
)

cluster_ids <- as.character(TARA_ALL$snn.louvianmlr_1)
TARA_ALL$Annotation <- unname(final_labels[cluster_ids])

# Tag unmapped
unmapped <- is.na(TARA_ALL$Annotation)
if (any(unmapped)) {
  TARA_ALL$Annotation[unmapped] <- paste0("Unmapped_cl", cluster_ids[unmapped])
  cat("WARNING: Unmapped clusters:", unique(cluster_ids[unmapped]), "\n")
}

cat("Initial annotation counts:\n")
print(sort(table(TARA_ALL$Annotation), decreasing = TRUE))
cat("\n")


# =============================================================================
# STEP 5: CD8 REFINEMENTS (sub-gating + splitting)
# =============================================================================

cat("══ STEP 5: CD8 REFINEMENTS ══\n\n")

# --- 5a: Derive has_TCR from clonal data ------------------------------------
tcr_col <- intersect(c("ClonalFrequency", "clonalFrequency", "Frequency",
                       "CTstrict", "CTgene", "CTaa"),
                     colnames(TARA_ALL@meta.data))[1]
if (is.na(tcr_col)) stop("No TCR/clonal column found in metadata.")
cat("Deriving has_TCR from '", tcr_col, "'\n", sep = "")
TARA_ALL$has_TCR <- !is.na(TARA_ALL@meta.data[[tcr_col]])

# --- 5b: Cluster 9 + merges: split by αβ TCR --------------------------------
is_trdv1 <- TARA_ALL$Annotation == "TRDV1+ pending"
TARA_ALL$Annotation[is_trdv1 & TARA_ALL$has_TCR == TRUE]  <- "TEMRA/CTL CD8"
TARA_ALL$Annotation[is_trdv1 & TARA_ALL$has_TCR == FALSE] <- "gd T cell (from cl9)"

cat("Cluster 9 split → TEMRA/CTL CD8:",
    sum(TARA_ALL$Annotation == "TEMRA/CTL CD8" & is_trdv1),
    "| gd T cell:",
    sum(TARA_ALL$Annotation == "gd T cell (from cl9)"), "\n")

# --- 5c: Cluster 12: split by CD8A ------------------------------------------
DefaultAssay(TARA_ALL) <- "RNA"
cd8a_expr <- GetAssayData(TARA_ALL, slot = "data")["CD8A", ]
is_cl12   <- TARA_ALL$Annotation == "DNAM1+ CD4/CD8 mixed"

TARA_ALL$Annotation[is_cl12 &  cd8a_expr > 0]                    <- "Transitional CD8"
TARA_ALL$Annotation[is_cl12 & (cd8a_expr <= 0 | is.na(cd8a_expr))] <- "Activated CD4 (cl12)"

cat("Cluster 12 → Transitional CD8:",
    sum(TARA_ALL$Annotation == "Transitional CD8"),
    "| Activated CD4:",
    sum(TARA_ALL$Annotation == "Activated CD4 (cl12)"), "\n")

# --- 5d: Cluster 1: sub-gate by ADT CD45RO / FAS ----------------------------
is_cl1    <- TARA_ALL$Annotation == "Memory CD8 pending"
cl1_cells <- colnames(TARA_ALL)[is_cl1]
cat("Cluster 1 sub-gating:", length(cl1_cells), "cells\n")

DefaultAssay(TARA_ALL) <- "ADT"
adt_data <- GetAssayData(TARA_ALL, slot = "data")

cd45ro_cl1 <- as.numeric(adt_data["CD45RO", cl1_cells])
fas_cl1    <- as.numeric(adt_data["FAS",    cl1_cells])
cd45ra_cl1 <- as.numeric(adt_data["CD45RA", cl1_cells])
cd62l_cl1  <- as.numeric(adt_data["SELL",   cl1_cells])

# DSB-normalized thresholds (adjust if needed)
cd45ro_threshold <- 0
fas_threshold    <- 0

gate <- case_when(
  cd45ro_cl1 > cd45ro_threshold ~ "Naive Intermediate CD8",
  fas_cl1    > fas_threshold    ~ "Tscm CD8",
  TRUE                          ~ "Naive 1 CD8"
)
names(gate) <- cl1_cells

print(table(gate))
TARA_ALL$Annotation[match(names(gate), colnames(TARA_ALL))] <- gate

# --- 5e: Gating diagnostic plots --------------------------------------------
gate_df <- data.frame(
  CD45RO = cd45ro_cl1, CD45RA = cd45ra_cl1,
  FAS = fas_cl1, CD62L = cd62l_cl1,
  Gate = gate[cl1_cells], stringsAsFactors = FALSE
)

gate_cols <- c("Naive 1 CD8"            = "#4E79A7",
               "Tscm CD8"               = "#59A14F",
               "Naive Intermediate CD8"  = "#E15759")

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
    labs(title = "CD45RO") + theme_minimal()
) / (
  ggplot(gate_df, aes(FAS, fill = Gate)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = gate_cols) +
    geom_vline(xintercept = fas_threshold, linetype = "dashed") +
    labs(title = "FAS/CD95") + theme_minimal()
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(out_gating, "Gate_CD45RO_vs_CD45RA.png"), p_g1,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_gating, "Gate_CD45RO_vs_FAS.png"), p_g2,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_gating, "Gate_Histograms.png"), p_hist,
       width = 10, height = 10, dpi = 300, bg = "white")

cat("\n")


# =============================================================================
# STEP 6: REMOVE FLAGGED CLUSTERS
# =============================================================================

cat("══ STEP 6: REMOVING FLAGGED CLUSTERS ══\n\n")

flagged <- grepl("^FLAG_", TARA_ALL$Annotation)
cat("Removing", sum(flagged), "cells from flagged clusters:\n")
print(table(TARA_ALL$Annotation[flagged]))

TARA_ALL <- subset(TARA_ALL, cells = colnames(TARA_ALL)[!flagged])
TARA_ALL$Annotation <- factor(TARA_ALL$Annotation)

cat("\n══ FINAL ANNOTATION COUNTS ══\n")
print(sort(table(TARA_ALL$Annotation), decreasing = TRUE))
cat("\n")


# =============================================================================
# STEP 7: METADATA CORRECTION (HEI clinical ground truth)
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

# Initialize columns if they don't exist
if (!"Timepoint_Group" %in% colnames(TARA_ALL@meta.data)) {
  TARA_ALL$Timepoint_Group <- NA_character_
}
if (!"Viral_Load_num" %in% colnames(TARA_ALL@meta.data)) {
  TARA_ALL$Viral_Load_num <- NA_real_
}
if (!"Age" %in% colnames(TARA_ALL@meta.data)) {
  TARA_ALL$Age <- NA_real_
}

for (sname in names(sample_metadata)) {
  mask <- TARA_ALL$orig.ident == sname
  if (sum(mask) == 0) next
  meta <- sample_metadata[[sname]]
  TARA_ALL$Timepoint_Group[mask] <- meta$group
  TARA_ALL$Viral_Load_num[mask]  <- meta$vl
  TARA_ALL$Age[mask]             <- meta$age
}

# Rename orig.ident → PID_AgeM for HEI
TARA_ALL$orig.ident_raw <- TARA_ALL$orig.ident
ident_rename <- sapply(names(sample_metadata), function(s) {
  paste0(sub("_.*$", "", s), "_", sample_metadata[[s]]$age, "m")
})
names(ident_rename) <- names(sample_metadata)
new_ids <- ident_rename[TARA_ALL$orig.ident]
new_ids[is.na(new_ids)] <- TARA_ALL$orig.ident[is.na(new_ids)]
TARA_ALL$orig.ident <- unname(new_ids)

cat("Timepoint_Group:\n")
print(table(TARA_ALL$Timepoint_Group))
cat("\n")


# =============================================================================
# STEP 8: ANNOTATED UMAPs
# =============================================================================

p_final <- DimPlot2(
  TARA_ALL, reduction = "wnn.umap", group.by = "Annotation",
  label = TRUE, box = TRUE, repel = TRUE, label.color = "black",
  cols = "default", label.size = 3, pt.size = 0.5
) + ggtitle("TARA: Final Annotation")

ggsave(file.path(out_clusters, "TARA_Final_Annotation.png"),
       p_final, width = 15, height = 10, dpi = 300, bg = "white")

p_comp <- ClusterDistrBar(
  TARA_ALL$orig.ident, TARA_ALL$Annotation,
  cols = "default", flip = FALSE, border = "black"
) + theme(axis.title.x = element_blank())

ggsave(file.path(out_clusters, "TARA_Cluster_Distribution.png"),
       p_comp, width = 25, height = 14, dpi = 300, bg = "white")


# =============================================================================
# STEP 9: MARKER LISTS + HELPERS
# =============================================================================

rna_markers <- list(
  `CD8 Identity`        = c("CD8A","CD8B","CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2"),
  `Naive-Quiescence`    = c("CCR7","SELL","TCF7","LEF1","BACH2","IL7R","BCL2","S1PR1","KLF2","FOXP1","SATB1"),
  `Effector-Cytotox`    = c("GZMB","GZMA","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","KLRG1","CX3CR1","FGFBP2","FCGR3A"),
  `Effector TFs`        = c("TBX21","EOMES","RUNX3","PRDM1","ZEB2","ID2"),
  `Exhaustion`          = c("TOX","TOX2","PDCD1","HAVCR2","TIGIT","LAG3","CTLA4","ENTPD1","CD160","CD244"),
  `Memory-Stemness`     = c("IL7R","CD27","CD28","BACH2","TCF7","BCL2","BCL6","CXCR5","CCR7"),
  `Activation-Prolif`   = c("MKI67","HLA-DRA","HLA-DRB1","CD38","FAS","ICOS","TNFRSF9","CXCR6","CD69"),
  `GammaDelta TCR`      = c("TRDV1","TRDV2","TRGV9","TRDC","TRGC1","TRGC2"),
  `MAIT`                = c("TRAV1-2","SLC4A10","KLRB1","ZBTB16","RORC"),
  `NK-like`             = c("TYROBP","KLRD1","KIR2DL3","KIR3DL1","NCAM1","KLRF1","KLRC1","KLRC2"),
  `Type I IFN`          = c("ISG15","IFIT1","IFIT2","IFIT3","IFI44L","MX1","OAS1","STAT1","IRF7"),
  `Chemokines`          = c("CCL3","CCL4","CCL4L2","CCL5","XCL1","XCL2","IFNG","TNF","IL2","CSF2"),
  `CD4 T cell`          = c("CD4","IL7R","FOXP3","IL2RA","CTLA4","IKZF2","CCR4","RORC","GATA3","TBX21","BCL6","CXCR5","MAF","ICOS","BATF"),
  `B cell`              = c("CD19","MS4A1","CD79A","CD79B","PAX5","BANK1","BLK","BLNK","TCL1A","IGHM","IGHD","IGKC","IGLC2"),
  `Plasma`              = c("SDC1","XBP1","MZB1","JCHAIN","IGHG1","IGHA1","PRDM1","IRF4"),
  `NK cell`             = c("NCAM1","KLRF1","KLRC1","KLRB1","NCR1","NCR3","FCGR3A","TYROBP","FCER1G","SH2D1B","SPON2","CLIC3"),
  `Mono-Mac`            = c("CD14","LYZ","S100A8","S100A9","S100A12","FCGR3A","CSF1R","CD68","MARCO","VCAN","FCN1","MNDA"),
  `DC`                  = c("FCER1A","CLEC10A","CD1C","CLEC9A","XCR1","BATF3","IRF8","IRF4","LAMP3","LILRA4","CLEC4C","IL3RA"),
  `Platelet`            = c("PPBP","PF4","GP9","ITGA2B","TUBB1","TREML1"),
  `Erythrocyte`         = c("HBB","HBA1","HBA2","ALAS2","SLC4A1","GYPA")
)

rna_markers_flat <- unique(unlist(rna_markers))
DefaultAssay(TARA_ALL) <- "ADT"
all_adt_markers <- sort(rownames(TARA_ALL[["ADT"]]))

scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    rng <- max(x) - min(x)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x)) / rng
  }))
  colnames(scaled) <- colnames(mat)
  scaled
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

export_avg_expression <- function(obj, group_by, rna_feats, adt_feats, out_d) {
  dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
  DefaultAssay(obj) <- "RNA"
  rna_avail <- intersect(rna_feats, rownames(obj[["RNA"]]))
  avg_rna <- AverageExpression(obj, assays = "RNA", features = rna_avail,
                               group.by = group_by, slot = "data")$RNA
  colnames(avg_rna) <- gsub("^g ", "", colnames(avg_rna))
  write.csv(avg_rna, file.path(out_d, "AvgExpr_RNA.csv"))
  write.csv(round(scale_01(as.matrix(avg_rna)), 4), file.path(out_d, "AvgExpr_RNA_scaled.csv"))
  
  DefaultAssay(obj) <- "ADT"
  adt_avail <- intersect(adt_feats, rownames(obj[["ADT"]]))
  avg_adt <- AverageExpression(obj, assays = "ADT", features = adt_avail,
                               group.by = group_by, slot = "data")$ADT
  colnames(avg_adt) <- gsub("^g ", "", colnames(avg_adt))
  write.csv(avg_adt, file.path(out_d, "AvgExpr_ADT.csv"))
  write.csv(round(scale_01(as.matrix(log10(avg_adt + 1))), 4), file.path(out_d, "AvgExpr_ADT_scaled.csv"))
}


# =============================================================================
# STEP 10: EXPORT FULL-OBJECT AVG EXPRESSION
# =============================================================================

export_avg_expression(TARA_ALL, "Annotation", rna_markers_flat, all_adt_markers, out_avgexpr)
message("✓ Full-object CSVs → ", out_avgexpr)


# =============================================================================
# STEP 11: DEFINE LINEAGE SUBSETS
# =============================================================================

lineage_map <- list(
  CD8 = c("Naive 1 CD8", "Tscm CD8", "Naive Intermediate CD8", "Naive 2 CD8",
          "TEMRA/CTL CD8", "gd T cell (from cl9)", "Transitional CD8", "Tex CD8"),
  
  CD4 = c("Naive CD4", "Central Memory CD4", "Th2/Th17 Effector Memory CD4",
          "Activated CD4 (cl12)", "Treg"),
  
  NK  = c("CD56dim NK (CD122hi)", "CD56dim NK", "CD56bright NK", "NKT cell"),
  
  B_cell = c("Mature Naive B", "Transitional B", "Activated B", "TNF+ B cell"),
  
  Myeloid = c("Classical Monocyte", "Non-classical Monocyte",
              "Intermediate Monocyte", "pDC"),
  
  GammaDelta = c("Vd1 gd T cell", "Vg9Vd2 gd T cell",
                 "MAIT / Innate-like T cell"),
  
  Other_T = c("CD4dim Naive T cell", "Naive-like T cell", "T cell (unresolved)")
)

lineage_map <- lineage_map[sapply(lineage_map, length) > 0]


# =============================================================================
# STEP 12: LINEAGE EXPLORATION LOOP
# =============================================================================

cat("\n══ STEP 12: LINEAGE EXPLORATION ══\n\n")

lineage_objects <- list()

for (ln in names(lineage_map)) {
  cls <- lineage_map[[ln]]
  # Only keep clusters that actually exist in the data
  cls <- cls[cls %in% levels(TARA_ALL$Annotation)]
  if (length(cls) == 0) { cat(ln, ": no matching clusters, skipping\n"); next }
  
  cat(sprintf("── %s (%d clusters) ──\n", ln, length(cls)))
  
  obj <- subset(TARA_ALL, subset = Annotation %in% cls)
  obj$Annotation <- droplevels(factor(obj$Annotation, levels = cls))
  
  if (ncol(obj) < 50) { cat("   Skipping (< 50 cells)\n\n"); next }
  cat(sprintf("   %d cells\n", ncol(obj)))
  lineage_objects[[ln]] <- obj
  
  lin_dir <- file.path(lineage_base, ln)
  
  # Avg expression
  export_avg_expression(obj, "Annotation", rna_markers_flat, all_adt_markers,
                        file.path(lin_dir, "AvgExpression"))
  
  # Relevant RNA
  DefaultAssay(obj) <- "RNA"
  rel_rna <- intersect(rna_markers_flat, rownames(obj[["RNA"]]))
  
  # Violins
  generate_violin_plots(obj, rel_rna, "RNA", file.path(lin_dir, "Violins_RNA"), "Annotation")
  generate_violin_plots(obj, all_adt_markers, "ADT", file.path(lin_dir, "Violins_ADT"), "Annotation")
  
  # Feature plots
  generate_feature_plots(obj, rel_rna, "RNA", file.path(lin_dir, "FeaturePlots_RNA"))
  generate_feature_plots(obj, all_adt_markers, "ADT", file.path(lin_dir, "FeaturePlots_ADT"))
  
  # DotPlots
  DefaultAssay(obj) <- "RNA"
  p_dot <- DotPlot(obj, features = rel_rna, group.by = "Annotation",
                   cols = c("lightgrey", "#D73027"), dot.scale = 6) +
    RotatedAxis() + ggtitle(paste(ln, "— RNA")) +
    theme(axis.text.x = element_text(size = 7))
  ggsave(file.path(lin_dir, paste0("DotPlot_RNA_", ln, ".png")),
         p_dot, width = min(26, 6 + 0.15 * length(rel_rna)),
         height = max(6, 2 + 0.5 * length(cls)), dpi = 300, bg = "white")
  
  DefaultAssay(obj) <- "ADT"
  p_dot_a <- DotPlot(obj, features = all_adt_markers, group.by = "Annotation",
                     cols = c("lightgrey", "#08519C"), dot.scale = 5) +
    RotatedAxis() + ggtitle(paste(ln, "— ADT")) +
    theme(axis.text.x = element_text(size = 6))
  ggsave(file.path(lin_dir, paste0("DotPlot_ADT_", ln, ".png")),
         p_dot_a, width = min(38, 6 + 0.12 * length(all_adt_markers)),
         height = max(6, 2 + 0.5 * length(cls)), dpi = 300, bg = "white")
  
  # UMAP
  p_u <- DimPlot2(obj, reduction = "wnn.umap", group.by = "Annotation",
                  label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.8) +
    ggtitle(paste(ln, "Subset"))
  ggsave(file.path(lin_dir, paste0("UMAP_", ln, ".png")),
         p_u, width = 10, height = 8, dpi = 300, bg = "white")
  
  # Save subset
  qs_save(obj, file.path(saved_dir, paste0("TARA_", ln, "_subset.qs2")))
  cat(sprintf("   ✓ %s complete\n\n", ln))
}


# =============================================================================
# STEP 13: CD8 VALIDATION (pairwise DE)
# =============================================================================

if ("CD8" %in% names(lineage_objects)) {
  TARA_cd8 <- lineage_objects[["CD8"]]
  cat("══ STEP 13: CD8 VALIDATION ══\n\n")
  
  DefaultAssay(TARA_cd8) <- "RNA"
  
  run_de <- function(obj, id1, id2, outdir, label) {
    sub <- subset(obj, subset = Annotation %in% c(id1, id2))
    n1 <- sum(sub$Annotation == id1); n2 <- sum(sub$Annotation == id2)
    cat(sprintf("  %s: %d vs %d\n", label, n1, n2))
    if (n1 < 10 | n2 < 10) { cat("    Skipping\n"); return(NULL) }
    de <- tryCatch(
      FindMarkers(sub, ident.1 = id1, group.by = "Annotation",
                  test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.15),
      error = function(e) { message("    Failed: ", e$message); NULL })
    if (!is.null(de)) {
      write.csv(de, file.path(outdir, paste0("DE_", gsub("[^A-Za-z0-9._-]","_",label), ".csv")))
    }
    de
  }
  
  run_de(TARA_cd8, "Tscm CD8", "Naive 1 CD8", out_gating, "Tscm_vs_Naive1")
  run_de(TARA_cd8, "Naive Intermediate CD8", "Naive 1 CD8", out_gating, "NaiveInt_vs_Naive1")
  run_de(TARA_cd8, "Naive 1 CD8", "Naive 2 CD8", out_gating, "Naive1_vs_Naive2")
}


# =============================================================================
# STEP 14: SAVE FINAL OBJECTS
# =============================================================================

cat("\n══ STEP 14: SAVING ══\n\n")

qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_annotated_final.qs2"))
cat("✓ TARA_ALL_annotated_final.qs2\n")

if (exists("TARA_cd8")) {
  qs_save(TARA_cd8, file.path(saved_dir, "TARA_cd8_final.qs2"))
  cat("✓ TARA_cd8_final.qs2\n")
}


################################################################################
# DONE
################################################################################

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " ANNOTATION COMPLETE\n",
        "══════════════════════════════════════════════════════════════\n",
        "\n",
        " CD8:        Naive 1, Tscm, Naive Intermediate, Naive 2,\n",
        "             TEMRA/CTL, gd T cell, Transitional, Tex\n",
        " CD4:        Naive, Central Memory, Th2/Th17 EM, Activated, Treg\n",
        " NK:         CD56dim (CD122hi), CD56dim, CD56bright, NKT\n",
        " B:          Mature Naive, Transitional, Activated, TNF+\n",
        " Myeloid:    Classical Mono, Non-classical Mono, Intermediate Mono, pDC\n",
        " GammaDelta: Vd1, Vg9Vd2, MAIT/Innate-like\n",
        "\n",
        " REMOVED: FLAG_LowQuality, FLAG_Doublet_T_B,\n",
        "          FLAG_Mixed_Plasmablast, FLAG_APC_Unclear\n",
        "\n",
        " Output:  ", out_dir, "\n",
        " Objects: ", saved_dir, "\n",
        "══════════════════════════════════════════════════════════════\n"
)
