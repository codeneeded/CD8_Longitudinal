#Load Required Libraries
library(Seurat)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scCustomize)
library(Polychrome)
library(viridis)
library(readxl)
library(Seurat)
library(dplyr)
library(ggplot2)

# Load annotated Seurat object
TARA_ALL <- readRDS("~/Documents/CD8_Longitudinal/saved_R_data/TARA_ALL_post_annotation.rds")

# Define output directories
dir_hei_heu <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/HEIvsHEU"
dir_heu_huu <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/HEUvsHUU"
dir.create(dir_hei_heu, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_heu_huu, recursive = TRUE, showWarnings = FALSE)

# Set identity to manual annotation
Idents(TARA_ALL) <- "Manual_Annotation"

# Make sure group label is correct
TARA_ALL$Comparison_Group <- TARA_ALL$Condition  # values: HEI, HEU, HUU
TARA_ALL$PID <- sub("_.*", "", TARA_ALL$orig.ident)
DefaultAssay(TARA_ALL) <- "RNA"

# Initialize containers
de_stats_hei_heu <- list()
de_stats_heu_huu <- list()
cluster_counts <- data.frame()

# Loop over clusters
for (cluster in levels(TARA_ALL)) {
  message("Processing cluster: ", cluster)
  
  # Subset cluster
  cells_in_cluster <- WhichCells(TARA_ALL, idents = cluster)
  subset_cluster <- subset(TARA_ALL, cells = cells_in_cluster)
  
  # Check that required groups exist
  comp_table <- table(subset_cluster$Comparison_Group)
  
  # ---- HEI vs HEU ----
  if (all(c("HEI", "HEU") %in% names(comp_table)) && all(comp_table[c("HEI", "HEU")] >= 10)) {
    de <- FindMarkers(
      subset_cluster,
      ident.1 = "HEI", ident.2 = "HEU",
      group.by = "Comparison_Group",
      test.use = "MAST",
      latent.vars = c("nCount_RNA", "PID")
    )
    write.csv(de, file = file.path(dir_hei_heu, paste0(gsub("[:/ ]", "_", cluster), "_HEIvsHEU.csv")))
    de_stats_hei_heu[[cluster]] <- de
    
    # Count DE genes (adj. p < 0.05)
    num_sig <- sum(de$p_val_adj < 0.05, na.rm = TRUE)
    cluster_counts <- rbind(cluster_counts, data.frame(Cluster = cluster, DE_Genes = num_sig))
  }
  
  # ---- HEU vs HUU ----
  if (all(c("HEU", "HUU") %in% names(comp_table)) && all(comp_table[c("HEU", "HUU")] >= 10)) {
    de <- FindMarkers(
      subset_cluster,
      ident.1 = "HEU", ident.2 = "HUU",
      group.by = "Comparison_Group",
      test.use = "MAST",
      latent.vars = c("nCount_RNA", "PID")
    )
    write.csv(de, file = file.path(dir_heu_huu, paste0(gsub("[:/ ]", "_", cluster), "_HEUvsHUU.csv")))
    de_stats_heu_huu[[cluster]] <- de
  }
}

# ---- Plot: HEI vs HEU ----
ggplot(cluster_counts, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  xlab("Cluster") + ylab("# of DE Genes (adj. p < 0.05)") +
  ggtitle("TARA: Clusters Ranked by DE Genes (HEI vs HEU)") +
  ggsave(filename = file.path(dir_hei_heu, "Cluster_DE_Gene_Counts_HEIvsHEU.png"), width = 10, height = 6)

# ---- Plot: HEU vs HUU ----
# Build cluster count table for HEU vs HUU
cluster_counts_heuhuu <- data.frame(
  Cluster = names(de_stats_heu_huu),
  DE_Genes = sapply(de_stats_heu_huu, function(x) sum(x$p_val_adj < 0.05, na.rm = TRUE))
)

ggplot(cluster_counts_heuhuu, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  xlab("Cluster") + ylab("# of DE Genes (adj. p < 0.05)") +
  ggtitle("TARA: Clusters Ranked by DE Genes (HEU vs HUU)") +
  ggsave(filename = file.path(dir_heu_huu, "Cluster_DE_Gene_Counts_HEUvsHUU.png"), width = 10, height = 6)

###################################### High vs Low Viral Load ############################################################################################
TARA_Entry <- subset(TARA_ALL, subset = Age %in% c(1, 2))
TARA_Entry <- subset(TARA_Entry, subset = !(Condition %in% c("HEU", "HUU")))

TARA_Entry$Viral_Load_Category <- ifelse(
  as.numeric(as.character(TARA_Entry$Viral_Load)) >= 100000,
  "High",
  "Low"
)

# Setup: Output directory
# -----------------------------
dir_vl_highlow <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_Entry/HighvsLowVL"
dir.create(dir_vl_highlow, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Assign identities and group labels
# -----------------------------

# Set cluster identity to manual annotations
Idents(TARA_Entry) <- "Manual_Annotation"

# Define comparison group based on Viral Load Category ("High", "Low")
TARA_Entry$Comparison_Group <- TARA_Entry$Viral_Load_Category

# Set default assay to RNA
DefaultAssay(TARA_Entry) <- "RNA"

# -----------------------------
# Initialize containers to store results
# -----------------------------
de_stats_vl <- list()           # Stores DE results for each cluster
cluster_counts_vl <- data.frame()  # Tracks number of DE genes per cluster

# -----------------------------
# Loop through clusters and run DE
# -----------------------------
for (cluster in levels(TARA_Entry)) {
  message("Processing cluster: ", cluster)
  
  # Subset cells in the current cluster
  cells_in_cluster <- WhichCells(TARA_Entry, idents = cluster)
  subset_cluster <- subset(TARA_Entry, cells = cells_in_cluster)
  
  # Check if both High and Low groups exist with at least 10 cells each
  comp_table <- table(subset_cluster$Comparison_Group)
  
  if (all(c("High", "Low") %in% names(comp_table)) && all(comp_table[c("High", "Low")] >= 10)) {
    
    # Perform differential expression using MAST, adjusting for RNA content and individual (PID)
    de <- FindMarkers(
      subset_cluster,
      ident.1 = "High", ident.2 = "Low",
      group.by = "Comparison_Group",
      test.use = "MAST",
      latent.vars = c("nCount_RNA")
    )
    
    # Save DE results to CSV
    out_file <- file.path(dir_vl_highlow, paste0(gsub("[:/ ]", "_", cluster), "_HighvsLowVL.csv"))
    write.csv(de, file = out_file)
    
    # Store DE result in list
    de_stats_vl[[cluster]] <- de
    
    # Count number of significant genes (adj p < 0.05)
    num_sig <- sum(de$p_val_adj < 0.05, na.rm = TRUE)
    cluster_counts_vl <- rbind(cluster_counts_vl, data.frame(Cluster = cluster, DE_Genes = num_sig))
  }
}

# -----------------------------
# Plot: Number of DE genes per cluster
# -----------------------------
plot_file <- file.path(dir_vl_highlow, "Cluster_DE_Gene_Counts_HighvsLowVL.png")

ggplot(cluster_counts_vl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  xlab("Cluster") + 
  ylab("# of DE Genes (adj. p < 0.05)") +
  ggtitle("TARA Entry: Clusters Ranked by DE Genes (High vs Low Viral Load)")

ggsave(filename = plot_file, width = 10, height = 6)


################# Pre-ART vs Post-ART ###########################
# -----------------------------
# Step 1: Subset to only HIV-infected infants (HEI)
# -----------------------------
TARA_HEI <- subset(TARA_ALL, subset = Condition == "HEI")

# -----------------------------
# Step 2: Ensure Viral_Load is numeric
# -----------------------------
TARA_HEI$Viral_Load <- as.numeric(as.character(TARA_HEI$Viral_Load))

# -----------------------------
# Step 3: Create Timepoint_Group label
# PreART_Entry: Age <= 2 months
# PostART_Suppressed: Age > 2 AND VL < 200
# -----------------------------
TARA_HEI$Timepoint_Group <- ifelse(
  TARA_HEI$Age <= 2, 
  "PreART_Entry",
  ifelse(TARA_HEI$Viral_Load < 200, "PostART_Suppressed", NA)
)
# -----------------------------
# Step 4: Filter to just Pre vs Post comparison cells
# -----------------------------
TARA_HEI_compare <- subset(TARA_HEI, subset = !is.na(Timepoint_Group))

# -----------------------------
# Set output directory
# -----------------------------
dir_post_pre <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/PostARTvsPreART"
dir.create(dir_post_pre, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Set identities and metadata
# -----------------------------
Idents(TARA_HEI_compare) <- "Manual_Annotation"
TARA_HEI_compare$Comparison_Group <- TARA_HEI_compare$Timepoint_Group  # PreART_Entry / PostART_Suppressed
TARA_HEI_compare$PID <- sub("_.*", "", TARA_HEI_compare$orig.ident)  # Extract PID
DefaultAssay(TARA_HEI_compare) <- "RNA"

# -----------------------------
# Initialize storage containers
# -----------------------------
de_stats_post_pre <- list()
cluster_counts_post_pre <- data.frame()

# -----------------------------
# Loop through each cluster and run DE with MAST
# -----------------------------
for (cluster in levels(TARA_HEI_compare)) {
  message("Processing cluster: ", cluster)
  
  # Subset cells from the current cluster
  cells_in_cluster <- WhichCells(TARA_HEI_compare, idents = cluster)
  subset_cluster <- subset(TARA_HEI_compare, cells = cells_in_cluster)
  
  # Ensure both groups exist with enough cells
  comp_table <- table(subset_cluster$Comparison_Group)
  
  if (all(c("PreART_Entry", "PostART_Suppressed") %in% names(comp_table)) &&
      all(comp_table[c("PreART_Entry", "PostART_Suppressed")] >= 10)) {
    
    # Run differential expression using MAST with PID correction
    de <- FindMarkers(
      subset_cluster,
      ident.1 = "PostART_Suppressed", ident.2 = "PreART_Entry",
      group.by = "Comparison_Group",
      test.use = "MAST",
      latent.vars = c("nCount_RNA", "PID")
    )
    
    # Save results
    out_file <- file.path(dir_post_pre, paste0(gsub("[:/ ]", "_", cluster), "_PostARTvsPreART.csv"))
    write.csv(de, file = out_file)
    
    # Store results in list
    de_stats_post_pre[[cluster]] <- de
    
    # Track number of significant DE genes
    num_sig <- sum(de$p_val_adj < 0.05, na.rm = TRUE)
    cluster_counts_post_pre <- rbind(cluster_counts_post_pre, data.frame(Cluster = cluster, DE_Genes = num_sig))
  }
}

# -----------------------------
# Plot: Number of DE genes per cluster
# -----------------------------
plot_file <- file.path(dir_post_pre, "Cluster_DE_Gene_Counts_PostARTvsPreART.png")

ggplot(cluster_counts_post_pre, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  xlab("Cluster") +
  ylab("# of DE Genes (adj. p < 0.05)") +
  ggtitle("TARA HEI: Clusters Ranked by DE Genes (Post-ART vs Pre-ART)")

ggsave(filename = plot_file, width = 10, height = 6)
