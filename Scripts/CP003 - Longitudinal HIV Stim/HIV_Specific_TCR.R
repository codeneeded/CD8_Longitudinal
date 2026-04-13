
library(qs2)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)
library(SeuratExtend)   # for DimPlot2 / VlnPlot2
library(scCustomize)
library(dplyr)
library(ggplot2)
library(qs2)

# -----------------------------
# 0) Paths
# -----------------------------
save_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data"
in_qs2  <- file.path(save_dir, "CP003_RNA_integrated_CCA_MNN_withTCR.qs2")



# -----------------------------
# 1) Load object
# -----------------------------
seu <- qs_read(in_qs2)
seu$mnn_clusters_rna <- factor(
  seu$mnn_snn_res.0.6,
  levels = as.character(sort(as.numeric(unique(as.character(seu$mnn_snn_res.0.6)))))
)
Idents(seu) <- "mnn_clusters_rna"
levels(Idents(seu))
seu$Sample
seu$CTstrict
#################### MKI67 ################
############################################################
# CP003 Proliferation / HIV-specific TCR labeling
############################################################

# ----------------------------- #
# 0) Paths
# ----------------------------- #
plot_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/Proliferation"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- #
# 1) Save requested plots
# ----------------------------- #

p_scatter <- FeatureScatter(seu, feature1 = "MKI67", feature2 = "CD8A")
ggsave(
  filename = file.path(plot_dir, "FeatureScatter_MKI67_vs_CD8A.png"),
  plot = p_scatter,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

p_vln1 <- VlnPlot2(seu, "MKI67")
ggsave(
  filename = file.path(plot_dir, "VlnPlot2_MKI67.png"),
  plot = p_vln1,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

p_vln2 <- VlnPlot2(
  seu,
  features = "MKI67",
  cols = "default",
  show.mean = TRUE,
  mean_colors = c("red", "blue")
)
ggsave(
  filename = file.path(plot_dir, "VlnPlot2_MKI67_withMean.png"),
  plot = p_vln2,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ----------------------------- #
# 2) Label HIV-specific TCR in metadata
# Definition:
# - MKI67 positive ~= expression >= 0.4
# - clonalFrequency > 1
# - ignore 0 / NA / non-expanded as not HIV-specific
# ----------------------------- #

freq_num <- suppressWarnings(as.numeric(as.character(seu$clonalFrequency)))
cluster_num <- seu$mnn_snn_res.0.6

seu$HIV_Specific_TCR <- ifelse(
  !is.na(cluster_num) & !(cluster_num %in% c(0, 1)) &
    !is.na(freq_num) & freq_num > 1,
  "HIV-Specific TCR",
  "Other"
)

seu$HIV_Specific_TCR <- factor(
  seu$HIV_Specific_TCR,
  levels = c("Other", "HIV-Specific TCR")
)
### Cell Cycle
# Run cell cycle scoring (uses built-in Seurat gene lists)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seu <- CellCycleScoring(
  seu,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# ---- Summary table: phase proportions per cluster ----
cc_cluster <- seu@meta.data %>%
  group_by(mnn_snn_res.0.6, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mnn_snn_res.0.6) %>%
  mutate(proportion = n / sum(n))

write.csv(
  cc_cluster,
  file = file.path(plot_dir, "CellCycle_by_Cluster.csv"),
  row.names = FALSE
)

# ---- Plot 1: stacked bar of phase proportions per cluster ----
p_cc_bar <- ggplot(cc_cluster,
                   aes(x = as.factor(mnn_snn_res.0.6),
                       y = proportion,
                       fill = Phase)) +
  geom_col() +
  scale_fill_manual(values = c("G1" = "#AEC6CF", "S" = "#FFD700", "G2M" = "#E63946")) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cluster",
    y = "Proportion of cells",
    fill = "Cell Cycle Phase",
    title = "Cell Cycle Phase by Cluster"
  )

ggsave(
  file.path(plot_dir, "CellCycle_Phase_by_Cluster.png"),
  p_cc_bar,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# ---- Plot 2: UMAP colored by Phase ----
p_cc_umap <- DimPlot2(
  seu,
  features = "Phase",
  reduction = "umap.mnn.rna",
  theme = NoAxes()
)

ggsave(
  file.path(plot_dir, "UMAP_CellCycle_Phase.png"),
  p_cc_umap,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ---- Plot 3: S and G2M scores per cluster as violin ----
p_cc_violin <- VlnPlot(
  seu,
  features = c("S.Score", "G2M.Score"),
  group.by = "mnn_snn_res.0.6",
  pt.size = 0,
  ncol = 2
)

ggsave(
  file.path(plot_dir, "CellCycle_Scores_Violin_by_Cluster.png"),
  p_cc_violin,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)
# ----------------------------- #
# 3) Extract HIV-specific clone list
# Includes:
# - clone
# - sample
# - cluster
# - clonalFrequency
# - number of cells
# - mean MKI67
# - where they came from
# ----------------------------- #

meta_hiv <- data.frame(
  Cell = colnames(seu),
  Sample = seu$Sample,
  Cluster = seu$mnn_snn_res.0.6,
  CTstrict = seu$CTstrict,
  clonalFrequency = seu$clonalFrequency,
  HIV_Specific_TCR = seu$HIV_Specific_TCR,
  MKI67 = mki67_expr,
  stringsAsFactors = FALSE
)

meta_hiv$clonalFrequency_num <- suppressWarnings(as.numeric(as.character(meta_hiv$clonalFrequency)))

hiv_cells <- meta_hiv %>%
  filter(
    HIV_Specific_TCR == "HIV-Specific TCR",
    !is.na(CTstrict),
    CTstrict != "",
    !is.na(clonalFrequency_num),
    clonalFrequency_num > 1
  )

# Per-cell table
write.csv(
  hiv_cells,
  file = file.path(plot_dir, "CP003_HIV_Specific_TCR_cells.csv"),
  row.names = FALSE
)

# Clone-level summary
hiv_clone_summary <- hiv_cells %>%
  group_by(CTstrict, Sample, Cluster) %>%
  summarise(
    n_cells = n(),
    max_clonalFrequency = max(clonalFrequency_num, na.rm = TRUE),
    mean_MKI67 = mean(MKI67, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_clonalFrequency), desc(n_cells))

# ----------------------------- #
# 2) Label HIV-specific TCR in metadata
# Definition:
# - NOT from cluster 0 or 1 (mnn_snn_res.0.6)
# - clonalFrequency > 1
# - ignore 0 / NA / non-expanded as not HIV-specific
# ----------------------------- #

freq_num <- suppressWarnings(as.numeric(as.character(seu$clonalFrequency)))
cluster_num <- seu$mnn_snn_res.0.6

seu$HIV_Specific_TCR <- ifelse(
  !is.na(cluster_num) & !(cluster_num %in% c(0, 1)) &
    !is.na(freq_num) & freq_num > 1,
  "HIV-Specific TCR",
  "Other"
)

seu$HIV_Specific_TCR <- factor(
  seu$HIV_Specific_TCR,
  levels = c("Other", "HIV-Specific TCR")
)

# Optional finer annotation if you want to see component logic later
seu$Not_Cluster0or1 <- ifelse(
  !is.na(cluster_num) & !(cluster_num %in% c(0, 1)),
  "Not_Cluster0or1", "Cluster0or1"
)
seu$Expanded_TCR_gt1 <- ifelse(!is.na(freq_num) & freq_num > 1, "Expanded", "Not_Expanded")

# ----------------------------- #
# 3) Extract HIV-specific clone list
# Includes:
# - clone
# - sample
# - cluster
# - clonalFrequency
# - number of cells
# - mean MKI67
# - where they came from
# ----------------------------- #

mki67_expr <- FetchData(seu, vars = "MKI67")[, 1]

meta_hiv <- data.frame(
  Cell = colnames(seu),
  Sample = seu$Sample,
  Cluster = seu$mnn_snn_res.0.6,
  CTstrict = seu$CTstrict,
  clonalFrequency = seu$clonalFrequency,
  HIV_Specific_TCR = seu$HIV_Specific_TCR,
  MKI67 = mki67_expr,
  stringsAsFactors = FALSE
)

meta_hiv$clonalFrequency_num <- suppressWarnings(as.numeric(as.character(meta_hiv$clonalFrequency)))

hiv_cells <- meta_hiv %>%
  filter(
    HIV_Specific_TCR == "HIV-Specific TCR",
    !is.na(CTstrict),
    CTstrict != "",
    !is.na(clonalFrequency_num),
    clonalFrequency_num > 1
  )

# Per-cell table
write.csv(
  hiv_cells,
  file = file.path(plot_dir, "CP003_HIV_Specific_TCR_cells.csv"),
  row.names = FALSE
)

# Clone-level summary
hiv_clone_summary <- hiv_cells %>%
  group_by(CTstrict, Sample, Cluster) %>%
  summarise(
    n_cells = n(),
    max_clonalFrequency = max(clonalFrequency_num, na.rm = TRUE),
    mean_MKI67 = mean(MKI67, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_clonalFrequency), desc(n_cells))

write.csv(
  hiv_clone_summary,
  file = file.path(plot_dir, "CP003_HIV_Specific_TCR_clone_summary.csv"),
  row.names = FALSE
)

# Clone -> all samples/clusters it appears in among HIV-specific cells
hiv_clone_list <- hiv_cells %>%
  group_by(CTstrict) %>%
  summarise(
    Samples = paste(sort(unique(as.character(Sample))), collapse = ";"),
    Clusters = paste(sort(unique(as.character(Cluster))), collapse = ";"),
    n_cells = n(),
    max_clonalFrequency = max(clonalFrequency_num, na.rm = TRUE),
    mean_MKI67 = mean(MKI67, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_clonalFrequency), desc(n_cells))

write.csv(
  hiv_clone_list,
  file = file.path(plot_dir, "CP003_HIV_Specific_TCR_clone_list.csv"),
  row.names = FALSE
)

# ----------------------------- #
# 3B) Quick check of specific target clones
# ----------------------------- #

target_clones <- c(
  "TRAV19.TRAJ39.TRAC;TGTGCTCTGAGTGTTTACCCCATTAATAATGCAGGCAACATGCTCACCTTT_TRBV2.NA.TRBJ2-2.TRBC2;TGTGCCAGCAGCCCATTTATACTAGCGGGGGCCGGGGAGCTGTTTTTT",
  "TRAV19.TRAJ48.TRAC;TGTGCTCTGAGTGAGGCGGGGATATCTAACTTTGGAAATGAGAAATTAACCTTT_TRBV14.TRBD1.TRBJ2-7.TRBC2;TGTGCCAGCAGCCAAGATCAAAGGGACAGGGGGCTGAACGAGCAGTACTTC",
  "TRAV3.TRAJ29.TRAC;TGTGCTGTGAGACCATTTTCAGGAAACACACCTCTTGTCTTT_TRBV24-1.TRBD1.TRBJ2-1.TRBC2;TGTGCCACCAGTGAGCCCAACAGGGGGCCTGAGGATGAGCAGTTCTTC",
  "NA;NA_TRBV2.NA.TRBJ2-2.TRBC2;TGTGCCAGCAGCCCATTTATACTAGCGGGGGCCGGGGAGCTGTTTTTT"
)

# Label cells carrying any of these target clones
seu$Target_Clone_Status <- ifelse(
  !is.na(seu$CTstrict) & seu$CTstrict %in% target_clones,
  "Target Clone",
  "Other"
)

seu$Target_Clone_Status <- factor(
  seu$Target_Clone_Status,
  levels = c("Other", "Target Clone")
)

# Combined status: target clone and whether it meets your HIV-specific rule
seu$Target_Clone_HIV_Status <- ifelse(
  !is.na(seu$CTstrict) & seu$CTstrict %in% target_clones & seu$HIV_Specific_TCR == "HIV-Specific TCR",
  "Target Clone + HIV-Specific TCR",
  ifelse(
    !is.na(seu$CTstrict) & seu$CTstrict %in% target_clones,
    "Target Clone only",
    "Other"
  )
)

seu$Target_Clone_HIV_Status <- factor(
  seu$Target_Clone_HIV_Status,
  levels = c("Other", "Target Clone only", "Target Clone + HIV-Specific TCR")
)

# Quick table
target_clone_table <- meta_hiv %>%
  filter(!is.na(CTstrict), CTstrict %in% target_clones) %>%
  group_by(CTstrict, Sample, Cluster, HIV_Specific_TCR) %>%
  summarise(
    n_cells = n(),
    max_clonalFrequency = max(clonalFrequency_num, na.rm = TRUE),
    mean_MKI67 = mean(MKI67, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(CTstrict, Sample, Cluster)

target_clone_table

# ----------------------------- #
# 4) Visualizations of HIV-specific TCRs
# ----------------------------- #

plot_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/Proliferation"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# 4A UMAP highlighting HIV-specific TCR cells (scCustom)
############################################################

p_umap_hiv <- DimPlot_scCustom(
  seu,
  group.by = "HIV_Specific_TCR",
  reduction = "umap.mnn.rna"
)

ggsave(
  file.path(plot_dir, "UMAP_HIV_Specific_TCR.png"),
  p_umap_hiv,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

p_umap_hiv_split <- DimPlot_scCustom(
  seu,
  group.by = "HIV_Specific_TCR",
  reduction = "umap.mnn.rna",
  split.by = "Sample",
  split_seurat = TRUE
)

ggsave(
  file.path(plot_dir, "UMAP_HIV_Specific_TCR_SplitBySample.png"),
  p_umap_hiv_split,
  width = 16,
  height = 10,
  dpi = 300,
  bg = "white"
)

############################################################
# 4A2 UMAP highlighting HIV-specific TCR cells (SeuratExtend DimPlot2)
############################################################

p_umap_hiv_dim2 <- DimPlot2(
  seu,
  features = "HIV_Specific_TCR",
  reduction = "umap.mnn.rna",
  cols = c("Other" = "grey85", "HIV-Specific TCR" = "#E63946"),
  theme = NoAxes()
)

ggsave(
  file.path(plot_dir, "UMAP_HIV_Specific_TCR_DimPlot2.png"),
  p_umap_hiv_dim2,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Split by sample using DimPlot2
p_umap_hiv_dim2_split <- DimPlot2(
  seu,
  features = "HIV_Specific_TCR",
  reduction = "umap.mnn.rna",
  split.by = "Sample",
  cols = c("Other" = "grey85", "HIV-Specific TCR" = "#E63946"),
  theme = NoAxes()
)

ggsave(
  file.path(plot_dir, "UMAP_HIV_Specific_TCR_DimPlot2_SplitBySample.png"),
  p_umap_hiv_dim2_split,
  width = 16,
  height = 10,
  dpi = 300,
  bg = "white"
)

############################################################
# 4B Number of HIV-specific cells per sample
############################################################

plot_df <- hiv_cells %>%
  count(Sample)

p_sample_bar <- ggplot(plot_df, aes(x = Sample, y = n)) +
  geom_col() +
  theme_classic(base_size = 14) +
  labs(
    y = "Number of HIV-specific cells",
    x = "Sample"
  )

ggsave(
  file.path(plot_dir, "HIV_Specific_Cells_by_Sample.png"),
  p_sample_bar,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

############################################################
# 4C HIV-specific cells by cluster
############################################################

plot_df <- hiv_cells %>%
  count(Cluster)

p_cluster_bar <- ggplot(plot_df, aes(x = reorder(as.character(Cluster), n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    y = "Number of HIV-specific cells",
    x = "Cluster"
  )

ggsave(
  file.path(plot_dir, "HIV_Specific_Cells_by_Cluster.png"),
  p_cluster_bar,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

############################################################
# 4D Cluster composition per sample
############################################################

plot_df <- hiv_cells %>%
  count(Sample, Cluster)

p_cluster_sample <- ggplot(plot_df,
                           aes(x = Sample, y = n, fill = as.factor(Cluster))) +
  geom_col(position = "fill") +
  theme_classic(base_size = 14) +
  labs(
    y = "Proportion of HIV-specific cells",
    x = "Sample",
    fill = "Cluster"
  )

ggsave(
  file.path(plot_dir, "Cluster_Composition_HIV_Specific_by_Sample.png"),
  p_cluster_sample,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

############################################################
# 4E Top expanded HIV-specific clones
############################################################

top_clones <- hiv_cells %>%
  count(CTstrict, sort = TRUE) %>%
  slice_head(n = 20)

p_top_clones <- ggplot(top_clones,
                       aes(x = reorder(CTstrict, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    y = "Number of HIV-specific cells",
    x = "Clone (CTstrict)"
  )

ggsave(
  file.path(plot_dir, "Top_HIV_Specific_Clones.png"),
  p_top_clones,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

############################################################
# 4H MKI67 vs clone expansion
############################################################

clone_plot_df <- hiv_cells %>%
  group_by(CTstrict) %>%
  summarise(
    n_cells = n(),
    mean_MKI67 = mean(MKI67, na.rm = TRUE),
    max_clonalFrequency = max(clonalFrequency_num, na.rm = TRUE),
    .groups = "drop"
  )

p_clone_prolif <- ggplot(
  clone_plot_df,
  aes(x = max_clonalFrequency,
      y = mean_MKI67,
      size = n_cells)
) +
  geom_point(alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    x = "Max clonal frequency",
    y = "Mean MKI67",
    size = "Cells per clone"
  )

ggsave(
  file.path(plot_dir, "Clone_Expansion_vs_MKI67.png"),
  p_clone_prolif,
  width = 8,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ----------------------------- #
# 5) Save Seurat object with new metadata
# ----------------------------- #
qs2::qs_save(
  seu,
  file = file.path(save_dir, "seu_CP003_HIVSpecificTCR_annotated.qs")
)
DimPlot2(
  seu,
  features = "HIV_Specific_TCR",
  reduction = "umap.mnn.rna",
  split.by = "mnn_snn_res.0.6",
  cols = c("Other" = "grey85", "HIV-Specific TCR" = "#E63946"),
  theme = NoAxes()
)

