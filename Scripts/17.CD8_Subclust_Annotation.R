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
tara_cdnk$Viral_Load <- as.numeric(tara_cdnk$Viral_Load)

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


orig_timepoint_vl_table <- tara_cdnk@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(orig.ident) %>%
  summarise(
    PID = sub("_.*$", "", unique(orig.ident)[1]),
    Timepoint_Group = unique(Timepoint_Group)[1],
    Condition = unique(Condition)[1],
    Age_months = median(Age, na.rm = TRUE),
    Viral_Load = unique(Viral_Load)[1],
    .groups = "drop"
  ) %>%
  arrange(PID, Age_months)

orig_timepoint_vl_table
##### Annotate ##################
Idents(tara_cdnk) <- 'wsnn_res.0.4'

tara_cdnk <- subset(tara_cdnk, idents = setdiff(levels(Idents(tara_cdnk)), c("23", "24", "25", "26")))



######################################33


### Plot Output
setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/VDJ/Seurat_Overlays")


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

setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Pre-Annotation/Heatmaps/RNA")

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

############################################################
# Clonal succession genes heatmaps (RNA)
# 1) Z-score heatmap (like your make_heatmap workflow, but ONLY these genes)
# 2) Seurat DoHeatmap (cells shown; grouped by wsnn_res.0.4)
#
# Output dir:
# /home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Pre-Annotation/Heatmaps/RNA
############################################################
############################################################
# Clonal succession genes heatmaps (RNA)
# Restricted to clusters: 11,17,6,2,10,16,22
############################################################

setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Pre-Annotation/Heatmaps/RNA")

DefaultAssay(tara_cdnk) <- "RNA"

clusters_use <- c(11, 17, 6, 2, 10, 16, 22)

clonal_genes <- c(
  "GZMK","IL7R","TCF7","NELL2","CD28","TC2N","FOS","FOSB","JUN",
  "PRF1","GZMB","GZMH","FGFBP2","SLA","SLA2","PRDM1","ADGRG1","NR3C1",
  "PLAC8","HAVCR2","ISG15","IFIT3","IFI6","IFI27","IFI44L","IFI35","LY6E"
)

# ---------------------------- #
# Subset object to clusters
# ---------------------------- #
tara_sub <- subset(tara_cdnk, subset = wsnn_res.0.4 %in% clusters_use)

Idents(tara_sub) <- tara_sub$wsnn_res.0.4

# Keep only genes present
genes_use <- intersect(clonal_genes, rownames(tara_sub[["RNA"]]))
missing   <- setdiff(clonal_genes, genes_use)

message("Genes found (RNA): ", length(genes_use), "/", length(clonal_genes))
if (length(missing) > 0) message("Missing genes: ", paste(missing, collapse = ", "))

############################################################
# 1) Z-score cluster-level heatmap
############################################################

toplot <- CalcStats(
  tara_sub,
  features = genes_use,
  method   = "zscore",
  order    = "p",
  n        = length(genes_use)
)

p_z <- Heatmap(toplot, lab_fill = "zscore") +
  ggtitle("Clonal succession genes (RNA) — clusters 11,17,6,2,10,16,22")

ggsave("ClonalSuccession_Genes_RNA_zscore_selectedClusters.png",
       p_z, width = 10, height = 8, dpi = 300)

############################################################
# 2) Cell-level DoHeatmap
############################################################

p_cells <- DoHeatmap(
  tara_sub,
  features = genes_use,
  group.by = "wsnn_res.0.4",
  assay    = "RNA",
  slot     = "data"
) +
  ggtitle("Clonal succession genes (RNA) — DoHeatmap (selected clusters)")

ggsave("ClonalSuccession_Genes_RNA_DoHeatmap_cells_selectedClusters.png",
       p_cells, width = 14, height = 8, dpi = 300)


# Set default assay to ADT
DefaultAssay(tara_cdnk) <- "ADT"

# Set working directory
setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Pre-Annotation/Heatmaps/ADT")

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

# --- PAVE ---
# ---------------------------- #
# PAVE (force-include proteins)
# ---------------------------- #

# Make sure ADT assay is active
DefaultAssay(tara_cdnk) <- "ADT"

# ADT feature names
adt_feats <- rownames(tara_cdnk[["ADT"]])

# Synonyms / patterns to try (ADT naming varies a lot)
syn_map <- list(
  CD3  = c("^CD3$", "CD3[DEGZ]?", "CD3\\s*\\(", "CD3-"),
  CD8  = c("^CD8$", "^CD8A$", "^CD8B$", "CD8A", "CD8B", "CD8\\s*\\(", "CD8-"),
  CD56 = c("^CD56$", "CD56\\s*\\(", "CD56-", "NCAM1", "NCAM1\\s*\\("),
  `PD-1` = c("^PD-1$", "^PD1$", "PD-1\\s*\\(", "PD1\\s*\\(", "PDCD1", "CD279", "CD279\\s*\\(")
)

# Helper to find ADT hits for a given target
find_adt_hits <- function(target, patterns, feats) {
  hits <- unique(unlist(lapply(patterns, function(pat) {
    grep(pat, feats, value = TRUE, ignore.case = TRUE)
  })))
  hits
}

# Get hits
hits_list <- lapply(names(syn_map), function(k) find_adt_hits(k, syn_map[[k]], adt_feats))
names(hits_list) <- names(syn_map)

# Print what we found (so you can sanity check)
for (nm in names(hits_list)) {
  message(nm, " hits: ", ifelse(length(hits_list[[nm]])==0, "<none>", paste(hits_list[[nm]], collapse = " | ")))
}

# Choose what to include:
# - If multiple hits exist, keep them all (safe; guarantees inclusion)
pave_hits_extra <- unique(unlist(hits_list))

# Your original hits (already found)
pave_hits_existing <- c("FCGR3A", "KIR2DL1", "KLRG1", "B3GAT1", "CCR5")

# Final set to force-include
pave_force <- unique(c(pave_hits_existing, pave_hits_extra))
pave_force <- intersect(pave_force, adt_feats)  # ensure valid

message("FORCING these ADT features into heatmap: ", paste(pave_force, collapse = ", "))

seu_pave <- subset(tara_cdnk, subset = wsnn_res.0.4 %in% c(11, 17, 6, 2, 10, 16, 22))
DefaultAssay(seu_pave) <- "ADT"

vf_old <- VariableFeatures(seu_pave)
VariableFeatures(seu_pave) <- unique(c(vf_old, pave_force))

make_heatmap(
  seu_pave,
  clusters = c(11, 17, 6, 2, 10, 16, 22),
  outfile = "PAVE_talk_specified_clusters_tara_cdnk_ADT.png",
  width = 10,
  height = 10
)


######### Cluster Distribution #########33
########### Cluster Distribution Plots ##############
setwd("/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/Pre-Annotation")


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


