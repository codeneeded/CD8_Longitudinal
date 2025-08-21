# Load libraries
library(Seurat)
library(dplyr)
library(SeuratExtend)
library(slingshot)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

# Define paths
load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"
file.name <- "TARA_ALL_post_annotation_comparisonsplit.rds"
output.dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Trajectory_Analysis"

# Load the Seurat object
tara_file <- file.path(load.path, file.name)
TARA_ALL <- readRDS(tara_file)

# Make sure output directory exists
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

################ MONOCOLE  ###########################
# HEI only
HEI <- subset(TARA_ALL, subset = Condition == "HEI")

# Define CD8+ lineage
cd8_clusters <- c(
  "1: Memory CD8 T cell", 
  "6: Naïve CD8 T cell", 
  "8: CTL-like", 
  "9: TRDV1+ CTL-like", 
  "27: GZMK+ CD8 T cell"
)

# Subset CD8+ cells
HEI_CD8 <- subset(HEI, subset = Manual_Annotation %in% cd8_clusters)
DefaultAssay(HEI_CD8) <- "RNA"
HEI_CD8 <- RunPCA(HEI_CD8, verbose = FALSE)
HEI_CD8 <- RunUMAP(HEI_CD8, dims = 1:30, verbose = FALSE)


# Convert to CDS
cds <- as.cell_data_set(HEI_CD8)

# Required by Monocle: gene_short_name in rowData
rowData(cds)$gene_short_name <- rownames(rowData(cds))

# Partition (only one partition for now)
partition_vec <- factor(rep(1, ncol(cds)))
names(partition_vec) <- colnames(cds)
cds@clusters$UMAP$partitions <- partition_vec

# Use Seurat cluster labels (Manual_Annotation)
cds@clusters$UMAP$clusters <- HEI_CD8$Manual_Annotation

# This assumes your UMAP is stored in the reduction slot as wnn.umap
reducedDims(cds)$UMAP <- Embeddings(HEI_CD8, "wnn.umap")

plot_cells(
  cds,
  color_cells_by = "Manual_Annotation",
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  group_label_size = 4
)

cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(
  cds,
  color_cells_by = "Manual_Annotation",
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  group_label_size = 4
)



root_cells <- colnames(cds)[cds@colData$Manual_Annotation == "6: Naïve CD8 T cell"]
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_branch_points = TRUE,
  label_roots = TRUE,
  label_leaves = TRUE
)

############ Bar plot
# Save pseudotime into cds metadata
cds$monocle3_pseudotime <- pseudotime(cds)

# Convert metadata to data frame
data.pseudo <- as.data.frame(colData(cds))

# Boxplot of pseudotime across Manual_Annotation (Seurat clusters)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Manual_Annotation, monocle3_pseudotime, median), fill = Manual_Annotation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Monocle3 Pseudotime", y = "CD8 Cluster (Manual_Annotation)") +
  guides(fill = "none")


########### Genes That Change Along Pseudotime
# Run test on the principal graph
deg_cd8 <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

# Get top changing genes
deg_cd8 %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head(10)

plot_genes_in_pseudotime(
  cds_subset = cds[c("TNFRSF18", "EFHD2", "RCAN3"), ],
  min_expr = 0.1,
  cell_size = 0.5,
  ncol = 3,
  color_cells_by = "Manual_Annotation",  # or "pseudotime"
  trend_formula = "~ splines::ns(pseudotime, df=3)"
)
