############################################################
# CP003: attach TCR (scRepertoire) to integrated Seurat object
# - Keep 001A and 001B separate in Seurat (recommended)
# - Keep TCR samples separate for combineExpression
############################################################

library(qs2)
library(Seurat)
library(scRepertoire)
library(dplyr)
library(stringr)
library(scCustomize)   # DimPlot_scCustom (optional, used below)
library(patchwork)
# -----------------------------
# Paths
# -----------------------------
save_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data"
out_qs2  <- file.path(save_dir, "CP003_RNA_integrated_CCA_MNN.qs2")

vdj_base <- "/home/akshay-iyer/CP003_multi_out/outs/per_sample_outs"

# -----------------------------
# 1) Load Seurat (qs2)
# -----------------------------
seu <- qs_read(out_qs2)

# Ensure orig.ident is your library id (should already be, but just to confirm)
# table(seu$orig.ident)

# -----------------------------
# 2) Rename Seurat cells to: orig.ident + "_" + raw_barcode
#    (matches scRepertoire convention best)
# -----------------------------
barcodes <- rownames(seu[[]])

# Rename cells: remove ONLY the leading "_" character
seu <- RenameCells(seu, new.names = sub("^_", "", rownames(seu[[]])))

barcodes <- rownames(seu[[]])


# -----------------------------
# 3) Read TCR contigs for ALL 4 libraries (unmerged)
# -----------------------------
libraries <- c("CP003_2m_001A","CP003_2m_001B","CP003_101m_003","CP003_101m_004")

contig_list <- setNames(
  lapply(libraries, function(lib) {
    read.csv(file.path(vdj_base, lib, "vdj_t", "filtered_contig_annotations.csv"),
             stringsAsFactors = FALSE)
  }),
  libraries
)

combined.TCR <- combineTCR(contig_list, samples = names(contig_list))

# -----------------------------
# 4) Add variables (sample-level metadata) to combined.TCR
# -----------------------------
samples <- names(combined.TCR)

months <- sub("^CP003_([^_]+)_.*$", "\\1", samples)    # "2m" / "101m"
ocm    <- sub("^CP003_[^_]+_(.*)$", "\\1", samples)    # "001A" "001B" "003" "004"
celltype_sort <- ifelse(ocm %in% c("001A","001B","003"), "CD8+", "CD8-")

combined.TCR <- addVariable(combined.TCR, "PID",           rep("CP003", length(samples)))
combined.TCR <- addVariable(combined.TCR, "Timepoint",     months)
combined.TCR <- addVariable(combined.TCR, "OCM_Barcode",   ocm)
combined.TCR <- addVariable(combined.TCR, "CellType_Sort", celltype_sort)

# -----------------------------
# 5) Attach TCR to Seurat
#    group.by MUST match the sample id used in Seurat cell names prefix.
#    Here: orig.ident should be "CP003_2m_001A" etc.
# -----------------------------
seu <- combineExpression(
  combined.TCR,
  seu,
  cloneCall = "strict",
  chain = "both",
  group.by = "sample",  
  cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500),
  proportion = FALSE
)

# -----------------------------
# Output dirs
# -----------------------------
base_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/VDJ/TCR/Seurat_Plots"
dir_hyper <- file.path(base_dir, "Hyperexpansion")
dir_occ   <- file.path(base_dir, "Clones_per_Cluster")
dir.create(dir_hyper, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_occ,   recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0) Prep metadata for your requested views
#    - "2m_merged" collapses 001A + 001B
#    - 101m stays split by CD8+ (003) vs CD8- (004)
# -----------------------------
# These columns should already exist from your earlier metadata creation:
# seu$Sample, seu$Timepoint, seu$OCM_Barcode, seu$CellType_Sort

seu$Timepoint_Merged <- ifelse(seu$Timepoint == "2m", "2m_merged", seu$Timepoint)

seu$Timepoint_Sort_View <- case_when(
  seu$Timepoint == "2m"   & seu$CellType_Sort == "CD8+" ~ "2m_merged_CD8+",
  seu$Timepoint == "101m" & seu$CellType_Sort == "CD8+" ~ "101m_CD8+",
  seu$Timepoint == "101m" & seu$CellType_Sort == "CD8-" ~ "101m_CD8-",
  TRUE ~ NA_character_
)

seu$Timepoint_Sort_View <- factor(
  seu$Timepoint_Sort_View,
  levels = c("2m_merged_CD8+", "101m_CD8+", "101m_CD8-")
)

# -----------------------------
# 1) Hyperexpansion / cloneSize on UMAP
# -----------------------------
# Choose the reduction you want to plot (change if needed)
# options you likely have: "umap.mnn.rna", "umap.cca.rna", "umap.pre"
red_use <- "umap.mnn.rna"

# Color palette (your style)
colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)
clone_cols <- rev(colorblind_vector[c(1,3,4,5,7)])

# A) All cells: cloneSize
p <- DimPlot_scCustom(
  seu,
  group.by = "cloneSize",
  reduction = red_use
) + scale_color_manual(values = clone_cols)

ggsave(file.path(dir_hyper, "CP003_cloneSize_AllCells.png"),
       plot = p, width = 8, height = 7, dpi = 400, bg = "white")

# B) Split by Sample (4 libraries) — helps check replicate consistency
p <- DimPlot_scCustom(
  seu,
  group.by = "cloneSize",
  reduction = red_use,
  split.by = "Sample",
  split_seurat = TRUE
) + scale_color_manual(values = clone_cols)

ggsave(file.path(dir_hyper, "CP003_cloneSize_SplitBy_Sample.png"),
       plot = p, width = 18, height = 7, dpi = 400, bg = "white")

# C) Your requested merged view: 2m merged + 101m split by CD8+/CD8-
p <- DimPlot_scCustom(
  seu,
  group.by = "cloneSize",
  reduction = red_use,
  split.by = "Timepoint_Sort_View",
  split_seurat = TRUE
) + scale_color_manual(values = clone_cols)

ggsave(file.path(dir_hyper, "CP003_cloneSize_SplitBy_2mMerged_vs_101m_CD8plus_CD8minus.png"),
       plot = p, width = 16, height = 7, dpi = 400, bg = "white")

# -----------------------------
# 2) Clonal occupancy (clonalOccupy)
#    - by clusters (mnn) and by predicted.celltype.l2
#    - both raw and proportion + export tables
# -----------------------------
# Choose your cluster column for occupancy
# You used mnn resolution 0.6. Use the factor you created earlier if present.
# Example expected: seu$mnn_clusters_rna
x_cluster <- "mnn_clusters_rna"
x_celltype <- "predicted.celltype.l2"

# A) Occupancy by clusters (raw)
p <- clonalOccupy(seu, x.axis = x_cluster, label = FALSE)
ggsave(file.path(dir_occ, "CP003_Clonal_Occupancy_by_MNNclusters_raw.png"),
       plot = p, width = 17, height = 11, dpi = 300, bg = "white")

# A2) Occupancy by clusters (proportion)
p <- clonalOccupy(seu, x.axis = x_cluster, proportion = TRUE, label = FALSE)
ggsave(file.path(dir_occ, "CP003_Clonal_Occupancy_by_MNNclusters_proportion.png"),
       plot = p, width = 17, height = 11, dpi = 300, bg = "white")

# A3) Export table
tbl <- clonalOccupy(seu, x.axis = x_cluster, exportTable = TRUE)
write.csv(tbl,
          file.path(dir_occ, "CP003_Clones_per_MNNcluster.csv"),
          row.names = FALSE)

# B) Occupancy by Azimuth cell types (raw)
p <- clonalOccupy(seu, x.axis = x_celltype, label = FALSE)
ggsave(file.path(dir_occ, "CP003_Clonal_Occupancy_by_Azimuth_raw.png"),
       plot = p, width = 27, height = 11, dpi = 300, bg = "white")

# B2) Occupancy by Azimuth cell types (proportion)
p <- clonalOccupy(seu, x.axis = x_celltype, proportion = TRUE, label = FALSE)
ggsave(file.path(dir_occ, "CP003_Clonal_Occupancy_by_Azimuth_proportion.png"),
       plot = p, width = 27, height = 11, dpi = 300, bg = "white")

# B3) Export table
tbl <- clonalOccupy(seu, x.axis = x_celltype, exportTable = TRUE)
write.csv(tbl,
          file.path(dir_occ, "CP003_Clones_per_AzimuthCellType.csv"),
          row.names = FALSE)

# -----------------------------
# 3) Same plots but restricted to your requested views
#    - keep only: 2m merged CD8+, 101m CD8+, 101m CD8-
# -----------------------------
seu_view <- subset(seu, subset = !is.na(Timepoint_Sort_View))

# cloneSize UMAP split (again but view-only; often cleaner)
p <- DimPlot_scCustom(
  seu_view,
  group.by = "cloneSize",
  reduction = red_use,
  split.by = "Timepoint_Sort_View",
  split_seurat = TRUE
) + scale_color_manual(values = clone_cols)

ggsave(file.path(dir_hyper, "CP003_VIEWONLY_cloneSize_SplitBy_2mMerged_vs_101m.png"),
       plot = p, width = 16, height = 7, dpi = 400, bg = "white")

# clonalOccupy view-only by cluster (raw + proportion + table)
p <- clonalOccupy(seu_view, x.axis = x_cluster, label = FALSE)
ggsave(file.path(dir_occ, "CP003_VIEWONLY_Clonal_Occupancy_by_MNNclusters_raw.png"),
       plot = p, width = 17, height = 11, dpi = 300, bg = "white")

p <- clonalOccupy(seu_view, x.axis = x_cluster, proportion = TRUE, label = FALSE)
ggsave(file.path(dir_occ, "CP003_VIEWONLY_Clonal_Occupancy_by_MNNclusters_proportion.png"),
       plot = p, width = 17, height = 11, dpi = 300, bg = "white")

tbl <- clonalOccupy(seu_view, x.axis = x_cluster, exportTable = TRUE)
write.csv(tbl,
          file.path(dir_occ, "CP003_VIEWONLY_Clones_per_MNNcluster.csv"),
          row.names = FALSE)

# clonalOccupy view-only by Azimuth (raw + proportion + table)
p <- clonalOccupy(seu_view, x.axis = x_celltype, label = FALSE)
ggsave(file.path(dir_occ, "CP003_VIEWONLY_Clonal_Occupancy_by_Azimuth_raw.png"),
       plot = p, width = 27, height = 11, dpi = 300, bg = "white")

p <- clonalOccupy(seu_view, x.axis = x_celltype, proportion = TRUE, label = FALSE)
ggsave(file.path(dir_occ, "CP003_VIEWONLY_Clonal_Occupancy_by_Azimuth_proportion.png"),
       plot = p, width = 27, height = 11, dpi = 300, bg = "white")

tbl <- clonalOccupy(seu_view, x.axis = x_celltype, exportTable = TRUE)
write.csv(tbl,
          file.path(dir_occ, "CP003_VIEWONLY_Clones_per_AzimuthCellType.csv"),
          row.names = FALSE)



############################################################
# CP003: Shared clones (2m CD8+ vs 101m CD8+) + UMAP highlight
############################################################
# -----------------------------
# Paths
# -----------------------------
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/VDJ/TCR/Shared_Clones"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Choose your UMAP reduction name here
red_use <- "umap.mnn.rna"

# -----------------------------
# 1) Define the two groups (cells)
# -----------------------------
cells_2m_cd8p <- WhichCells(seu, expression = Timepoint == "2m"   & CellType_Sort == "CD8+")
cells_101m_cd8p <- WhichCells(seu, expression = Timepoint == "101m" & CellType_Sort == "CD8+")

# Pull CTstrict values (clone IDs)
clones_2m <- seu$CTstrict[cells_2m_cd8p]
clones_101m <- seu$CTstrict[cells_101m_cd8p]

# Drop NA/empty
clones_2m <- clones_2m[!is.na(clones_2m) & clones_2m != ""]
clones_101m <- clones_101m[!is.na(clones_101m) & clones_101m != ""]

# Shared clone IDs (CTstrict)
shared_clones <- intersect(unique(clones_2m), unique(clones_101m))

# Save list
write.csv(
  data.frame(CTstrict = shared_clones),
  file.path(out_dir, "CP003_SharedClones_2mCD8plus_vs_101mCD8plus_CTstrict.csv"),
  row.names = FALSE
)

# -----------------------------
# 1b) Useful summary table: abundance / frequency in each timepoint
# -----------------------------
df_2m <- data.frame(
  CTstrict = seu$CTstrict[cells_2m_cd8p],
  clonalFrequency = seu$clonalFrequency[cells_2m_cd8p],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict) %>%
  summarise(
    nCells_2m = n(),
    mean_clonalFrequency_2m = mean(clonalFrequency, na.rm = TRUE),
    .groups = "drop"
  )

df_101m <- data.frame(
  CTstrict = seu$CTstrict[cells_101m_cd8p],
  clonalFrequency = seu$clonalFrequency[cells_101m_cd8p],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict) %>%
  summarise(
    nCells_101m = n(),
    mean_clonalFrequency_101m = mean(clonalFrequency, na.rm = TRUE),
    .groups = "drop"
  )

shared_tbl <- df_2m %>%
  inner_join(df_101m, by = "CTstrict") %>%
  arrange(desc(nCells_101m), desc(nCells_2m))

write.csv(
  shared_tbl,
  file.path(out_dir, "CP003_SharedClones_Summary_2m_vs_101m.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) Annotate each cell as Shared vs Unique vs NoTCR
# -----------------------------
seu$CloneShare_2m_vs_101m_CD8plus <- "No_TCR"

# Mark shared
shared_cells <- WhichCells(seu, expression = !is.na(CTstrict) & CTstrict %in% shared_clones)
seu$CloneShare_2m_vs_101m_CD8plus[shared_cells] <- "Shared_2m_and_101m"

# Mark unique within those CD8+ timepoints (optional but helpful)
uniq_2m <- setdiff(unique(clones_2m), shared_clones)
uniq_101m <- setdiff(unique(clones_101m), shared_clones)

uniq_2m_cells <- WhichCells(seu, expression = !is.na(CTstrict) & CTstrict %in% uniq_2m)
uniq_101m_cells <- WhichCells(seu, expression = !is.na(CTstrict) & CTstrict %in% uniq_101m)

seu$CloneShare_2m_vs_101m_CD8plus[uniq_2m_cells] <- "Unique_2m_CD8plus"
seu$CloneShare_2m_vs_101m_CD8plus[uniq_101m_cells] <- "Unique_101m_CD8plus"

seu$CloneShare_2m_vs_101m_CD8plus <- factor(
  seu$CloneShare_2m_vs_101m_CD8plus,
  levels = c("No_TCR", "Unique_2m_CD8plus", "Unique_101m_CD8plus", "Shared_2m_and_101m")
)

# -----------------------------
# 2a) UMAP: all cells, colored by shared/unique status
# -----------------------------
p1 <- DimPlot2(
  seu,
  reduction = red_use,
  group.by = "CloneShare_2m_vs_101m_CD8plus",
  pt.size = 0.25
)

ggsave(
  filename = file.path(out_dir, "CP003_UMAP_SharedCloneStatus_AllCells.png"),
  plot = p1,
  width = 10,
  height = 8,
  dpi = 400,
  bg = "white"
)

# -----------------------------
# 2b) UMAP split: 2m vs 101m (CD8+ only), shared clones highlighted
# -----------------------------
seu_cd8p_2m101m <- subset(
  seu,
  cells = union(cells_2m_cd8p, cells_101m_cd8p)
)

p2 <- DimPlot2(
  seu_cd8p_2m101m,
  reduction = red_use,
  group.by = "CloneShare_2m_vs_101m_CD8plus",
  split.by = "Timepoint",
  pt.size = 0.35
)

ggsave(
  filename = file.path(out_dir, "CP003_UMAP_SharedCloneStatus_CD8plusOnly_SplitBy_Timepoint.png"),
  plot = p2,
  width = 14,
  height = 6,
  dpi = 400,
  bg = "white"
)

# -----------------------------
# 2c) UMAP highlight (explicit): highlight only SHARED-clone cells
# -----------------------------
p3 <- DimPlot(
  seu_cd8p_2m101m,
  reduction = red_use,
  cells.highlight = WhichCells(seu_cd8p_2m101m, expression = CloneShare_2m_vs_101m_CD8plus == "Shared_2m_and_101m"),
  sizes.highlight = 0.6,
  pt.size = 0.2
)

ggsave(
  filename = file.path(out_dir, "CP003_UMAP_Highlight_SharedClonesOnly_CD8plusOnly.png"),
  plot = p3,
  width = 10,
  height = 8,
  dpi = 400,
  bg = "white"
)




# -----------------------------
# 6) Save updated Seurat (qs2 only)
# -----------------------------
out_qs2_tcr <- file.path(save_dir, "CP003_RNA_integrated_CCA_MNN_withTCR.qs2")
qs_save(seu, file = out_qs2_tcr)
