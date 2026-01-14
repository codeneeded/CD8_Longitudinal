library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(SeuratExtend)
library(scCustomize)
library(tidyr)
library(clustree)
library(scRepertoire)
library(monocle3)
library(SeuratWrappers)

##############################3
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "tara_cdnk_annotated.rds")
tara_cdnk <- readRDS(rds_in)

############################# Pseudotime - CD8 ##################################################
DefaultAssay(tara_cdnk) <-'RNA'
tara_cdnk$T_Cell_A
Idents(tara_cdnk) <- "T_Cell_A"
########## CD8 Pseudotime ##############

cd8_annots <- c(
  "CD8+ T_resting",
  "CD8+ T_naive_IL7R+",
  "CD8+ T_naive",
  "CD8+ T_memory",
  "CD8+ T_B-like",
  "CD8+ T_memory_costim+integrin",
  "CD8+ T_effector_prog",
  "NK-like_CD8+ T_PD1+",
  "NK-like_CD8+ T_CCR6+",
  "CTL",
  "γδ2_CTL",
  "NK-like_CD8+ T_KLRG1+",
  "NK-like_CD8+ T_KIR3DL1+"
)

cd8 <- subset(tara_cdnk, idents = cd8_annots)
cd8$annot <- Idents(cd8)  # keep a stable label column

####### Psuedotime


cd8_cds <- as.cell_data_set(cd8)
# Inject Seurat's WNN UMAP into Monocle
cd8_cds@int_colData@listData$reducedDims$UMAP <- cd8@reductions$wnn.umap@cell.embeddings

# Cluster and learn principal graph
cd8_cds <- cluster_cells(cd8_cds, reduction_method = "UMAP")
cd8_cds <- learn_graph(cd8_cds, use_partition = TRUE)

######### Define root nodes for pseudotime
naive_genes    <- intersect(c("TCF7","CCR7","LEF1","IL7R","SELL","LTB"), rownames(cd8[["RNA"]]))
effector_genes <- intersect(c("PRF1","GZMB","GZMH","NKG7","KLRG1","ZEB2","CX3CR1"), rownames(cd8[["RNA"]]))

cd8 <- AddModuleScore(cd8, features = list(naive_genes),    name = "naive_")
cd8 <- AddModuleScore(cd8, features = list(effector_genes), name = "eff_")

clust_scores <- cd8@meta.data |>
  transmute(annot = cd8$annot, naive = naive_1, eff = eff_1) |>
  group_by(annot) |>
  summarise(
    n     = n(),
    naive = mean(naive, na.rm=TRUE),
    eff   = mean(eff,   na.rm=TRUE),
    delta = naive - eff,
    .groups = "drop"
  ) |>
  arrange(desc(delta))

print(clust_scores)   # inspect: who looks earliest?

# Candidate root annotations (prefer clear naïve/resting labels, filter by delta>0)
# 2) Pick *all* annotations that look early by simple rules
min_cells <- 50                  # avoid tiny clusters
thr_delta <- 0                  # >=0 means naïve≥effector; raise to 0.05 if you want stricter
thr_naive <- quantile(clust_scores$naive, 0.50, na.rm=TRUE)  # above-median naïve score

########## Cluster based roots
root_annots <- c("CD8+ T_resting","CD8+ T_naive","CD8+ T_naive_IL7R+")
root_cells_A <- WhichCells(cd8, idents = root_annots)

cd8_cds_A <- cd8_cds
cd8_cds_A <- order_cells(cd8_cds_A, root_cells = root_cells_A)
pt_A <- pseudotime(cd8_cds_A)


# --- Pseudotime plot ---
p_pseudo <- plot_cells(
  cd8_cds_A,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2.5,
  group_label_size = 3.5,
  label_cell_groups = FALSE
)

ggsave(
  filename = "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Trajectory/CD8_pseudotime_plot.png",
  plot = p_pseudo,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white",# ensure PCA model exists for downstream plots (does not touch your WNN UMAP)
)
cd8_cds_A <- estimate_size_factors(cd8_cds_A)
cd8_cds_A <- preprocess_cds(cd8_cds_A, method = "PCA", num_dim = 50, norm_method = "log")



# --- Timepoint_Group plot ---
p_group <- plot_cells(
  cd8_cds_A,
  color_cells_by = "Timepoint_Group",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2.5,
  group_label_size = 3.5,
  label_cell_groups = FALSE
)

ggsave(
  filename = "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Trajectory/CD8_TimepointGroup_plot.png",
  plot = p_group,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white"
)

# --- Cluster plot ---
p_group <- plot_cells(
  cd8_cds_A,
  color_cells_by = "annot",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2.5,
  group_label_size = 3.5,
  label_cell_groups = FALSE
)

ggsave(
  filename = "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Trajectory/CD8_Cluster_plot.png",
  plot = p_group,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white"
)

# --- ident plot ---
p_group <- plot_cells(
  cd8_cds_A,
  color_cells_by = "orig.ident",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2.5,
  group_label_size = 3.5,
  label_cell_groups = FALSE
)

ggsave(
  filename = "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Trajectory/CD8_ident_plot.png",
  plot = p_group,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white"
)


########### Detect genes across pseudotime ########
# ensure PCA model exists for downstream plots (does not touch your WNN UMAP)
cd8_cds_A <- estimate_size_factors(cd8_cds_A)
cd8_cds_A <- preprocess_cds(cd8_cds_A, method = "PCA", num_dim = 50, norm_method = "log")

# Fit a smooth model for each gene along pseudotime
deg_pseudotime <- graph_test(
  cd8_cds_A,
  neighbor_graph = "principal_graph",   # Monocle3 pseudotime graph
  cores = 8
)

# Filter significant trajectory genes
deg_sig <- deg_pseudotime %>%
  arrange(q_value) %>%
  filter(q_value < 0.05)
head(deg_sig)

pt <- pseudotime(cd8_cds_A)
mat <- as.matrix(SingleCellExperiment::logcounts(cd8_cds_A))[rownames(deg_sig), names(pt)]

# Add gene column explicitly
deg_sig <- deg_sig %>%
  mutate(gene = rownames(.))

# Spearman correlation of each gene with pseudotime
rho <- apply(mat, 1, function(x) suppressWarnings(cor(x, pt, method = "spearman", use = "complete.obs")))
deg_sig$rho <- rho

# keep strongest, well-supported changers (tweak n as you like)
n_up <- 20; n_down <- 20

top_up <- deg_sig %>%
  filter(rho > 0) %>%
  arrange(q_value, desc(abs(rho))) %>%
  slice_head(n = n_up) %>%
  pull(gene)

top_down <- deg_sig %>%
  filter(rho < 0) %>%
  arrange(q_value, desc(abs(rho))) %>%
  slice_head(n = n_down) %>%
  pull(gene)

top_genes <- unique(c(top_up, top_down))
top_genes

### Add shortname metadata for moncole
# Add gene identifiers to rowData so Monocle can label facets
rowData(cd8_cds_A)$gene_id <- rownames(cd8_cds_A)
if (is.null(rowData(cd8_cds_A)$gene_short_name)) {
  rowData(cd8_cds_A)$gene_short_name <- rownames(cd8_cds_A)
}

#can use the default labeling
top_genes <- intersect(top_genes, rownames(cd8_cds_A))
top_genes <- unique(c(top_genes, "FASLG"))

# colored by pseudotime (continuous)
p_pt <- plot_genes_in_pseudotime(
  cd8_cds_A[top_genes, ],
  color_cells_by = "pseudotime",
  ncol = 4,
  min_expr = 0.1
)

# colored by Group  (continuous)
p_tg <- plot_genes_in_pseudotime(
  cd8_cds_A[top_genes, ],
  color_cells_by = "Timepoint_Group",
  ncol = 4,
  min_expr = 0.1
)

# colored by Cluster  (continuous)
p_clust <- plot_genes_in_pseudotime(
  cd8_cds_A[top_genes, ],
  color_cells_by = "annot",
  ncol = 4,
  min_expr = 0.1
)

# colored by ident  (continuous)
p_id <- plot_genes_in_pseudotime(
  cd8_cds_A[top_genes, ],
  color_cells_by = "orig.ident",
  ncol = 4,
  min_expr = 0.1
)

setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Trajectory")
ggsave("CD8_topGenes_pseudotime.png", p_pt,  width = 14, height = 10, dpi = 300, bg = "white")
ggsave("CD8_topGenes_pseudotime_Timepoint_Group.png", p_tg,  width = 14, height = 10, dpi = 300, bg = "white")
ggsave("CD8_topGenes_pseudotime_Cluster.png", p_clust,  width = 14, height = 10, dpi = 300, bg = "white")
ggsave("CD8_topGenes_pseudotime_Orig.png", p_id,  width = 14, height = 10, dpi = 300, bg = "white")

### Plots

