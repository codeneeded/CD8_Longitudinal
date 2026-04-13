############################################################
# CP003 scRepertoire (TCR) analysis from Cell Ranger VDJ outputs
# - Input folders:
#   /home/akshay-iyer/CP003_multi_out/outs/per_sample_outs/
#     CP003_2m_001A/
#     CP003_2m_001B/
#     CP003_101m_003/
#     CP003_101m_004/
# - Each has: vdj_t/filtered_contig_annotations.csv
# - Output saved under:
#   /home/akshay-iyer/Documents/CD8_Longitudinal/CP003/VDJ/TCR/...
############################################################


library(scRepertoire)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(readr)


# -----------------------------
# 0) Paths
# -----------------------------
in.path <- "/home/akshay-iyer/CP003_multi_out/outs/per_sample_outs/"

out.root <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/VDJ/TCR"
dir.create(out.root, recursive = TRUE, showWarnings = FALSE)

dir_clonal_viz <- file.path(out.root, "Clonal_Visualizations")
dir_cd3_comp   <- file.path(out.root, "CDR3_Composition")
dir_genes      <- file.path(out.root, "VJ_Gene_Usage")
dir_kmer       <- file.path(out.root, "Kmer")
dir_diversity  <- file.path(out.root, "Clonal_Diversity")
dir_overlap    <- file.path(out.root, "Clonal_Overlap")
dir_compare    <- file.path(out.root, "Clonal_Compare")
dir_tables     <- file.path(out.root, "Tables")

dirs <- c(dir_clonal_viz, dir_cd3_comp, dir_genes, dir_kmer, dir_diversity, dir_overlap, dir_compare, dir_tables)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# -----------------------------
# 1) Samples present (fixed names)
# -----------------------------
samples <- c("CP003_2m_001A", "CP003_2m_001B", "CP003_101m_003", "CP003_101m_004")

# -----------------------------
# 2) Read contigs into a named list (TCR only)
# -----------------------------
############################################################
# CP003: Merge 2m technical replicates at the contig level
# so combined.TCR has UNIQUE sample names (required for clonalOverlap)
############################################################

in.path <- "/home/akshay-iyer/CP003_multi_out/outs/per_sample_outs"

# Read each library contig
t_2m_001A  <- read.csv(file.path(in.path, "CP003_2m_001A",  "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)
t_2m_001B  <- read.csv(file.path(in.path, "CP003_2m_001B",  "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)
t_101m_003 <- read.csv(file.path(in.path, "CP003_101m_003", "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)
t_101m_004 <- read.csv(file.path(in.path, "CP003_101m_004", "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)

# Merge the 2m replicates into one contig table
t_2m_cd8p <- dplyr::bind_rows(t_2m_001A, t_2m_001B)

# Build contig list with UNIQUE sample names
contig_list_merged <- list(
  CP003_2m_CD8plus   = t_2m_cd8p,
  CP003_101m_CD8plus = t_101m_003,
  CP003_101m_CD8minus= t_101m_004
)


# -----------------------------
# 3) Combine for scRepertoire
# -----------------------------
# Combine (this is where duplicates get merged properly)
# Combine (now samples are unique, so clonalOverlap works)
combined.TCR <- combineTCR(
  contig_list_merged,
  samples = names(contig_list_merged)
)

# Sanity check
# -----------------------------
# 4) Add CP003 metadata (Timepoint / OCM / Sort)
# -----------------------------
# Sample names that exist AFTER you merged contigs + ran combineTCR()
samples <- names(combined.TCR)
# Parse: CP003_<months>_<CD8plus/CD8minus>
months <- sub("^CP003_([^_]+)_.*$", "\\1", samples)          # "2m" / "101m"
sort_tag <- sub("^CP003_[^_]+_(.*)$", "\\1", samples)        # "CD8plus" / "CD8minus"

celltype_sort <- ifelse(sort_tag == "CD8plus", "CD8+", "CD8-")

combined.TCR <- addVariable(combined.TCR, variable.name = "PID",           variables = rep("CP003", length(samples)))
combined.TCR <- addVariable(combined.TCR, variable.name = "Timepoint",     variables = months)
combined.TCR <- addVariable(combined.TCR, variable.name = "CellType_Sort", variables = celltype_sort)

# (Optional) keep the simplified sample label too
combined.TCR <- addVariable(combined.TCR, variable.name = "Sample_Merged", variables = samples)

# Convenience names
sample_2m_cd8p   <- "CP003_2m_CD8plus"
sample_101m_cd8p <- "CP003_101m_CD8plus"
sample_101m_cd8m <- "CP003_101m_CD8minus"

# Orders you’ll reuse
order_2m_vs_101m_cd8plus <- c(sample_2m_cd8p, sample_101m_cd8p)
order_all_merged         <- c(sample_2m_cd8p, sample_101m_cd8p, sample_101m_cd8m)
# -----------------------------
# 5) BASIC CLONAL VISUALIZATIONS
# -----------------------------
# Unique clones by Timepoint (strict, TRAB)
p <- clonalQuant(
  combined.TCR,
  cloneCall = "strict",
  chain = "both",
  group.by = "Timepoint",
  scale = FALSE
)
ggsave(file.path(dir_clonal_viz, "CP003_Unique_Clones_Strict_TRAB_by_Timepoint_raw.png"),
       plot = p, width = 18, height = 8, dpi = 300, bg = "white")

p <- clonalQuant(
  combined.TCR,
  cloneCall = "strict",
  chain = "both",
  group.by = "Timepoint",
  scale = TRUE
)
ggsave(file.path(dir_clonal_viz, "CP003_Unique_Clones_Strict_TRAB_by_Timepoint_scaled.png"),
       plot = p, width = 18, height = 8, dpi = 300, bg = "white")

# Unique clones by CellType_Sort
p <- clonalQuant(
  combined.TCR,
  cloneCall = "strict",
  chain = "both",
  group.by = "CellType_Sort",
  scale = TRUE
)

# -----------------------------
# 4) Add CP003 metadata (Timepoint / OCM / Sort)
# -----------------------------
# Parse: CP003_<months>_<OCM>
months <- sub("^CP003_([^_]+)_.*$", "\\1", samples)      # "2m" / "101m"
ocm    <- sub("^CP003_[^_]+_(.*)$", "\\1", samples)      # "001A" etc.

ggsave(file.path(dir_clonal_viz, "CP003_Unique_Clones_Strict_TRAB_by_Sort_scaled.png"),
       plot = p, width = 18, height = 8, dpi = 300, bg = "white")

# Clonal abundance
p <- clonalAbundance(combined.TCR, cloneCall = "strict", scale = FALSE)
ggsave(file.path(dir_clonal_viz, "CP003_MERGED_Clonal_Abundance_Strict_raw.png"),
       plot = p, width = 14, height = 10, dpi = 300, bg = "white")

p <- clonalAbundance(combined.TCR, cloneCall = "strict", scale = TRUE)
ggsave(file.path(dir_clonal_viz, "CP003_MERGED_Clonal_Abundance_Strict_scaled.png"),
       plot = p, width = 14, height = 10, dpi = 300, bg = "white")

# Clonal length
p <- clonalLength(combined.TCR, cloneCall = "aa", chain = "both", scale = FALSE)
ggsave(file.path(dir_clonal_viz, "CP003_MERGED_Clonal_Length_AA_TRAB_raw.png"),
       plot = p, width = 14, height = 10, dpi = 300, bg = "white")

p <- clonalLength(combined.TCR, cloneCall = "aa", chain = "both", scale = TRUE)
ggsave(file.path(dir_clonal_viz, "CP003_MERGED_Clonal_Length_AA_TRAB_scaled.png"),
       plot = p, width = 14, height = 10, dpi = 300, bg = "white")

# Homeostasis
p <- clonalHomeostasis(
  combined.TCR,
  cloneCall = "strict",
  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded = 1)
)
ggsave(file.path(dir_clonal_viz, "CP003_MERGED_Clonal_Homeostasis_Strict.png"),
       plot = p, width = 18, height = 8, dpi = 300, bg = "white")

# -----------------------------
# 6) CDR3 COMPOSITION
# -----------------------------
p <- percentAA(combined.TCR, chain = "TRA", aa.length = 20)
ggsave(file.path(dir_cd3_comp, "CP003_MERGED_Percent_AA_TRA.png"),
       plot = p, width = 20, height = 16, dpi = 300, bg = "white")

p <- percentAA(combined.TCR, chain = "TRB", aa.length = 20)
ggsave(file.path(dir_cd3_comp, "CP003_MERGED_Percent_AA_TRB.png"),
       plot = p, width = 20, height = 16, dpi = 300, bg = "white")

p <- positionalEntropy(combined.TCR, chain = "both", aa.length = 20)
ggsave(file.path(dir_cd3_comp, "CP003_MERGED_Positional_Entropy_TRAB.png"),
       plot = p, width = 16, height = 10, dpi = 300, bg = "white")

# -----------------------------
# 7) V/J gene usage heatmaps
# -----------------------------
p <- vizGenes(combined.TCR, x.axis = "TRAV", y.axis = NULL, plot = "heatmap", scale = TRUE)
ggsave(file.path(dir_genes, "CP003_MERGED_Heatmap_TRAV.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- vizGenes(combined.TCR, x.axis = "TRBV", y.axis = NULL, plot = "heatmap", scale = TRUE)
ggsave(file.path(dir_genes, "CP003_MERGED_Heatmap_TRBV.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- vizGenes(combined.TCR, x.axis = "TRAJ", y.axis = NULL, plot = "heatmap", scale = TRUE)
ggsave(file.path(dir_genes, "CP003_MERGED_Heatmap_TRAJ.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- vizGenes(combined.TCR, x.axis = "TRBJ", y.axis = NULL, plot = "heatmap", scale = TRUE)
ggsave(file.path(dir_genes, "CP003_MERGED_Heatmap_TRBJ.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- percentKmer(combined.TCR, cloneCall = "aa", chain = "TRA", motif.length = 3, top.motifs = 25)
ggsave(file.path(dir_kmer, "CP003_MERGED_Heatmap_TRA_kmer_len3_top25.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- percentKmer(combined.TCR, cloneCall = "aa", chain = "TRB", motif.length = 3, top.motifs = 25)
ggsave(file.path(dir_kmer, "CP003_MERGED_Heatmap_TRB_kmer_len3_top25.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")
# -----------------------------
# 8) Kmer motifs
# -----------------------------
p <- percentKmer(combined.TCR, cloneCall = "aa", chain = "TRA", motif.length = 3, top.motifs = 25)
ggsave(file.path(dir_kmer, "CP003_Heatmap_TRA_kmer_len3_top25.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- percentKmer(combined.TCR, cloneCall = "aa", chain = "TRB", motif.length = 3, top.motifs = 25)
ggsave(file.path(dir_kmer, "CP003_Heatmap_TRB_kmer_len3_top25.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

# -----------------------------
# 9) Diversity + overlap
# -----------------------------
p <- clonalDiversity(combined.TCR, cloneCall = "strict")
ggsave(file.path(dir_diversity, "CP003_Clonal_Diversity_Strict.png"),
       plot = p, width = 12, height = 9, dpi = 300, bg = "white")

p <- clonalOverlap(combined.TCR, cloneCall = "strict", method = "morisita")
ggsave(file.path(dir_overlap, "CP003_Clonal_Overlap_Strict_morisita.png"),
       plot = p, width = 16, height = 10, dpi = 300, bg = "white")

p <- clonalOverlap(combined.TCR, cloneCall = "strict", method = "raw")
ggsave(file.path(dir_overlap, "CP003_Clonal_Overlap_Strict_raw.png"),
       plot = p, width = 16, height = 10, dpi = 300, bg = "white")

# -----------------------------
# 10) Longitudinal clone tracking (CP003 only)
#     - Alluvial across your 4 libraries in a fixed order
# -----------------------------
# Make sure order matches EXACT sample names in combined.TCR
order_2m_vs_101m_cd8plus <- c("CP003_2m_CD8plus", "CP003_101m_CD8plus")

# -----------------------------
# 1) ALLUVIAL (counts)
# -----------------------------
p_alluvial <- clonalCompare(
  combined.TCR,
  top.clones = 20,
  samples = order_2m_vs_101m_cd8plus,
  order.by = order_2m_vs_101m_cd8plus,
  cloneCall = "strict",
  relabel.clones = TRUE,
  proportion = FALSE,
  graph = "alluvial"
)

ggsave(
  filename = file.path(dir_compare, "CP003_2m_vs_101m_CD8plus_MERGED_Top20_Strict_alluvial_counts.png"),
  plot = p_alluvial,
  width = 16,
  height = 10,
  dpi = 300,
  bg = "white"
)

tbl_alluvial <- clonalCompare(
  combined.TCR,
  top.clones = 50,
  samples = order_2m_vs_101m_cd8plus,
  order.by = order_2m_vs_101m_cd8plus,
  cloneCall = "strict",
  relabel.clones = TRUE,
  proportion = FALSE,
  graph = "alluvial",
  exportTable = TRUE
)

write.csv(
  tbl_alluvial,
  file.path(dir_tables, "CP003_2m_vs_101m_CD8plus_MERGED_Top50_Strict_alluvial_exportTable.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) STACKED (counts)
# -----------------------------
p_stacked <- clonalCompare(
  combined.TCR,
  top.clones = 20,
  samples = order_2m_vs_101m_cd8plus,
  order.by = order_2m_vs_101m_cd8plus,
  cloneCall = "strict",
  relabel.clones = TRUE,
  proportion = FALSE,
  graph = "area"
)

ggsave(
  filename = file.path(dir_compare, "CP003_2m_vs_101m_CD8plus_MERGED_Top20_Strict_stacked_counts.png"),
  plot = p_stacked,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

tbl_stacked <- clonalCompare(
  combined.TCR,
  top.clones = 50,
  samples = order_2m_vs_101m_cd8plus,
  order.by = order_2m_vs_101m_cd8plus,
  cloneCall = "strict",
  relabel.clones = TRUE,
  proportion = FALSE,
  graph = "area",
  exportTable = TRUE
)

write.csv(
  tbl_stacked,
  file.path(dir_tables, "CP003_2m_vs_101m_CD8plus_MERGED_Top50_Strict_stacked_exportTable.csv"),
  row.names = FALSE
)
