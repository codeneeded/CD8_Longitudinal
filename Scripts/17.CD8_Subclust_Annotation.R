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

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "tara_cdnk.rds")

tara_cdnk <- readRDS(rds_in)

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

T_Cell_A_map <- c(
  "19" = "CD4 T cell",
  "0"  = "CD8 T cell",
  "13" = "CD8 T cell",
  "18" = "CD8 T cell",
  "4"  = "CD8 T memory",
  "8"  = "CD8 T memory",
  "5"  = "CD8 T Naïve",
  "7"  = "CD8 T Naïve",
  "14" = "gamma delta 1",
  "20" = "gamma delta 2",
  "15" = "gamma delta 2",
  "12" = "NK 56 low",
  "9"  = "NK bright",
  "1"  = "NK classic",
  "3"  = "NK classic_CD161+",
  "11" = "NK-like CD8 T cell",
  "17" = "NK-like CD8 T cell",
  "21" = "NK-like CD8 T cell",
  "6"  = "NK-like CD8 T cell",
  "16" = "NK-like CD8 T cell",
  "22" = "NK-like CD8 T cell",
  "2"  = "NK-T",
  "10" = "NK-T"
)
# 3. Map annotations to each cell based on remaining clusters
cluster_ids <- as.character(Idents(tara_cdnk))
tara_cdnk$T_Cell_A <- unname(T_Cell_A_map[cluster_ids])

# 4. Preserve cluster order
cluster_levels <- levels(Idents(tara_cdnk))
tara_cdnk$T_Cell_A <- factor(tara_cdnk$T_Cell_A, levels = unique(T_Cell_A_map[cluster_levels]))
tara_cdnk$wsnn_res.0.4
DimPlot2(
  tara_cdnk,
  reduction = "wnn.umap",
  group.by = "wsnn_res.0.4",
  cols = 'light',
  label = TRUE, box = TRUE, repel = TRUE,
  label.color = "black", pt.size = 1
)

ggsave('cd8nk_subcluster_res0.4.png',
       dpi = 300,
       bg='white')


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

