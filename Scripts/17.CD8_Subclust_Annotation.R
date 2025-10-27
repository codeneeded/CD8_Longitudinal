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

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "tara_cdnk.rds")

tara_cdnk <- readRDS(rds_in) ### START FROM PSEUDOTIME

#################### LOAD VDJ DATA  ##########################

#Filenames

in.path <- "~/Documents/10x_Genomics/"

# Function to get all folder names within a specified path
get_folder_names <- function(in.path) {
  # List all directories within the specified path without recursion
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  
  # Filter out the root path itself if present
  folder_names <- folder_names[folder_names != ""]
  
  return(folder_names)
}

f_names <- get_folder_names(in.path)

# Split into TARA and EARTH names
TARA_names <- f_names[grepl("^C", f_names)]


for (name in f_names) {
  
  # Construct the file paths
  t_file <- paste0(in.path, name, '/per_sample_outs/TCR/filtered_contig_annotations.csv')
  b_file <- paste0(in.path, name, '/per_sample_outs/BCR/filtered_contig_annotations.csv')
  
  # Read the files
  t_data <- read.csv(t_file)
  b_data <- read.csv(b_file)
  
  # Create dynamically named variables
  assign(paste0(name, ".TCR"), t_data)
  assign(paste0(name, ".BCR"), b_data)
}

# Create Contig list

TARA.TCR <- paste(TARA_names, ".TCR", sep="")


TARA.contig_list.TCR <- as.list(mget(TARA.TCR))


#Combine For downstream Analysis

combined.TCR.TARA <- combineTCR(TARA.contig_list.TCR,samples = TARA_names)

# Merge with Seurat
tara_cdnk <- combineExpression(combined.TCR.TARA, 
                              tara_cdnk, 
                              cloneCall="strict",
                              chain='both',
                              group.by = 'sample',
                              cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                              500),
                              proportion = FALSE)

### Save output
#saveRDS(
#  tara_cdnk,
#  file = "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data/tara_cdnk_vdj.rds"
#)

######## Add metadata

tara_cdnk@meta.data <- tara_cdnk@meta.data %>%
  mutate(
    Viral_Load_Category = case_when(
      Condition == "HEI" & Viral_Load >= 100000 ~ "High",
      Condition == "HEI" & Viral_Load < 100000  ~ "Low",
      Condition == "HEU" ~ "HEU",
      Condition == "HUU" ~ "HUU",
      TRUE ~ NA_character_
    )
  )

tara_cdnk@meta.data <- tara_cdnk@meta.data %>%
  mutate(
    Timepoint_Group = case_when(
      Condition == "HEI" & Age <= 2 ~ "PreART_Entry",
      Condition == "HEI" & Viral_Load < 200 ~ "PostART_Suppressed",
      Condition == "HEI" & !is.na(Viral_Load) & Viral_Load >= 200 ~ "PostART_Unsuppressed",
      Condition %in% c("HEU", "HUU") ~ Condition,
      TRUE ~ NA_character_
    )
  )

##### Annotate ##################
Idents(tara_cdnk) <- 'wsnn_res.0.4'

tara_cdnk <- subset(tara_cdnk, idents = setdiff(levels(Idents(tara_cdnk)), c("23", "24", "25", "26")))



######################################33


### Plot Output
setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/VDJ")


#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "plasma", fixup = TRUE)

DimPlot_scCustom(tara_cdnk, group.by = "cloneSize", reduction = 'wnn.umap') +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave('tara_cdnk_Hyperexpansion.png',width=10,height=7)

DimPlot_scCustom(tara_cdnk, group.by = "cloneSize", reduction = 'wnn.umap',split.by = 'Condition',split_seurat = T) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave('tara_cdnk_Hyperexpansion_split_by_condition.png',width=14,height=7)





clonalOccupy(tara_cdnk, 
             x.axis = "wsnn_res.0.4",
             label = F)
ggsave('tara_cdnk_Clonal_Occupancy.png',width=13,height=11)


clonalOverlay(tara_cdnk, 
              reduction = "wnn.umap", 
              cutpoint = 10, 
              bins = 25, 
              facet.by = "orig.ident") + 
  guides(color = "none")
ggsave('TARA_cdnk_Clonal_Overlay_By_Sample.png',width=22,height=17)


clonalOverlay(tara_cdnk, 
              reduction = "wnn.umap", 
              cutpoint = 10, 
              bins = 25, 
              facet.by = "Condition") + 
  guides(color = "none")
ggsave('TARA_cdnk_Clonal_Overlay_By_Condition.png',width=17,height=9)


############# Heatmap Comparisons ##############################


################ Cluster Comparison Heatmaps ##################################
Idents(tara_cdnk)
### Function
make_heatmap <- function(seu, clusters, outfile, n_top_genes = 10, width = 10, height = 6) {
  # Subset
  seu_sub <- subset(seu, subset = wsnn_res.0.4 %in% clusters)
  
  # Get variable features
  genes <- VariableFeatures(seu_sub)
  
  # Calculate z-scores
  toplot <- CalcStats(seu_sub, features = genes, method = "zscore", order = "p", n = n_top_genes)
  
  # Heatmap
  p <- Heatmap(toplot, lab_fill = "zscore")
  
  # Save with custom size
  ggsave(outfile, p, width = width, height = height)
}

### RNA

DefaultAssay(tara_cdnk) <-'RNA'

setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Heatmaps/RNA")

make_heatmap(
  tara_cdnk,
  clusters = c(19, 0, 13, 18, 5, 7, 14, 20),
  outfile = "Resting_Naive_T_Cells_tara_cdnk.png",
  width = 10,
  height = 10
)

make_heatmap(
  tara_cdnk,
  clusters = c(15, 12, 6, 16, 22, 2, 10),
  outfile = "NK_and_NK_like_T_Cells_tara_cdnk.png",
  width = 10,
  height = 10
)
make_heatmap(
  tara_cdnk,
  clusters = c(4, 8),
  outfile = "CD8_T_Memory_tara_cdnk.png",
  width = 6,
  height = 5
)
make_heatmap(
  tara_cdnk,
  clusters = c(9, 1, 3, 6, 11, 16, 17, 21, 22),
  outfile = "NK_and_NK_like_Cells_tara_cdnk.png",
  width = 10,
  height = 10
)

# Set default assay to ADT
DefaultAssay(tara_cdnk) <- "ADT"

# Set working directory
setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Heatmaps/ADT")

# --- Resting / Naïve T Cells ---
make_heatmap(
  tara_cdnk,
  clusters = c(19, 0, 13, 18, 5, 7, 14, 20),
  outfile = "Resting_Naive_T_Cells_tara_cdnk_ADT.png",
  width = 10,
  height = 10
)

# --- NK and NK-like T Cells ---
make_heatmap(
  tara_cdnk,
  clusters = c(15, 12, 6, 16, 22, 2, 10),
  outfile = "NK_and_NK_like_T_Cells_tara_cdnk_ADT.png",
  width = 10,
  height = 10
)

# --- CD8 T Memory Cells ---
make_heatmap(
  tara_cdnk,
  clusters = c(4, 8),
  outfile = "CD8_T_Memory_tara_cdnk_ADT.png",
  width = 6,
  height = 5
)

# --- NK and NK-like Cells ---
make_heatmap(
  tara_cdnk,
  clusters = c(9, 1, 3, 6, 11, 16, 17, 21, 22),
  outfile = "NK_and_NK_like_Cells_tara_cdnk_ADT.png",
  width = 10,
  height = 10
)

######### Cluster Distribution #########33
########### Cluster Distribution Plots ##############
setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets")


ClusterDistrPlot(
  origin = tara_cdnk$orig.ident,
  cluster = tara_cdnk$wsnn_res.0.4,
  condition = tara_cdnk$Condition
)
ggsave('Clust_Distribution_tara_cdnk_By_Condition.png', width = 15, height = 13)

ClusterDistrPlot(
  origin = tara_cdnk$orig.ident, ### make tara entry first which is in the next section, needs reordering)
  cluster = tara_cdnk$wsnn_res.0.4,
  condition = tara_cdnk$Viral_Load_Category
)
ggsave('Clust_Distribution_tara_cdnk_By_Viral_Load.png', width = 15, height = 13)

ClusterDistrPlot(
  origin = tara_cdnk$orig.ident, ### make tara entry first which is in the next section, needs reordering)
  cluster = tara_cdnk$wsnn_res.0.4,
  condition = tara_cdnk$Timepoint_Group
)
ggsave('Clust_Distribution_Ttara_cdnk_By_ART.png', width = 15, height = 13)

##### Annotate ##################
Idents(tara_cdnk) <- 'wsnn_res.0.4'

tara_cdnk <- subset(tara_cdnk, idents = setdiff(levels(Idents(tara_cdnk)), c("23", "24", "25", "26")))

T_Cell_A_map <- c(
  "19" = "CD4+ T",
  "0"  = "CD8+ T_resting",
  "5"  = "CD8+ T_naive_IL7R+",
  "7"  = "CD8+ T_naive",
  "13" = "CD8+ T_B-like",
  "14" = "γδ1 T",
  "18" = "CD8+ T_effector_prog",
  "20" = "γδ2 T",
  "4"  = "CD8+ T_memory",
  "8"  = "CD8+ T_memory_costim+integrin",
  "2"  = "NK-T_CD64+",
  "6"  = "CTL",
  "10" = "NK-T",
  "12" = "NK_56low",
  "1"  = "NK_resting",
  "3"  = "NK_activated",
  "9"  = "NK_bright",
  "11" = "NK-like_CD8+ T_PD1+",
  "17" = "NK-like_CD8+ T_CCR6+",
  "21" = "NK_resting",
  "15" = "γδ2_CTL",
  "16" = "NK-like_CD8+ T_KLRG1+",
  "22" = "NK-like_CD8+ T_KIR3DL1+"
)

# 3. Map annotations to each cell based on remaining clusters
cluster_ids <- as.character(Idents(tara_cdnk))
tara_cdnk$T_Cell_A <- unname(T_Cell_A_map[cluster_ids])

# 4. Preserve cluster order
cluster_levels <- levels(Idents(tara_cdnk))
tara_cdnk$T_Cell_A <- factor(tara_cdnk$T_Cell_A, levels = unique(T_Cell_A_map[cluster_levels]))

DimPlot2(
  tara_cdnk,
  reduction = "wnn.umap",
  group.by = "T_Cell_A",
  cols = 'light',
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
)

ggsave('cd8nk_subcluster_annotated.png',
       dpi = 300,
       width = 14,
       bg='white')

saveRDS(
  tara_cdnk,
  file = "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data/tara_cdnk_annotated.rds"
)


##############################3
rds_in     <- file.path(saved_dir, "tara_cdnk_annotated.rds")
tara_cdnk <- readRDS(rds_in)

############################# Pseudotime ##################################################
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


plot_cells(
  cd8_cds_A,
  color_cells_by = "pseudotime",       # your annotation column
  label_groups_by_cluster = F,
  label_leaves = T,
  label_branch_points = T,
  graph_label_size = 2.5,
  group_label_size=3.5,
  label_cell_groups=FALSE
)

########### Detect genes across pseudotime ########
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


# colored by pseudotime (continuous)
p_pt <- plot_genes_in_pseudotime(
  cd8_cds_A[top_genes, ],
  color_cells_by = "pseudotime",
  ncol = 4,
  min_expr = 0.1
)

# colored by your clinical/timepoint group (HEU/HUU/Pre-/Post-ART)
p_grp <- plot_genes_in_pseudotime(
  cd8_cds_A[top_genes, ],
  color_cells_by = "Timepoint_Group",
  ncol = 4,
  min_expr = 0.1
)



setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Trajectory")
ggsave("CD8_topGenes_pseudotime.png", p_pt,  width = 14, height = 10, dpi = 300, bg = "white")
ggsave("CD8_topGenes_byGroup.png",  p_grp, width = 14, height = 10, dpi = 300, bg = "white")

### Plots
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
  bg = "white"
)


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

