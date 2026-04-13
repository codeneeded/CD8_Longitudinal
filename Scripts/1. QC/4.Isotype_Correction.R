#Load Required Libraries
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scCustomize)
library(SeuratWrappers)

library(future)
library(future.apply)

# Set up a parallel plan with multiple workers
plan(multisession, workers = 30)  # Adjust the number of workers based on your CPU cores

##set path to load data


##set path to load data
setwd('~/Documents/CD8_Longitudinal/Integration/Isotype')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"


load(paste0(load.path,'SeuratV5_SplitSeurat_Preintegration.RData'))

#### Delete all assays that are not RNA and ADT (including the sketched assay- this is to control size)
### we also get rid of the azimuth assays, however we keep the azimuth annotations
# Loop through each Seurat object in the list and keep only RNA and ADT assays
seurat_layer <- lapply(split_seurat, function(seurat_obj) {
  # Get the names of assays in the Seurat object
  current_assays <- names(seurat_obj@assays)
  
  # Define assays to keep
  assays_to_keep <- c("RNA", "ADT")
  
  # Identify assays to drop
  assays_to_drop <- setdiff(current_assays, assays_to_keep)
  
  # Drop unwanted assays
  for (assay in assays_to_drop) {
    seurat_obj[[assay]] <- NULL  # Remove the unwanted assay
  }
  
  return(seurat_obj)
})

rm(split_seurat)
gc()

seurat_layer <- Merge_Seurat_List(seurat_layer,merge.data = TRUE)

#####
# Start with the first object in the list
seurat_layer <- split_seurat[[1]]

# Sequentially merge the remaining objects
for (i in 2:length(split_seurat)) {
  seurat_layer <- merge(seurat_layer, y = split_seurat[[i]])
}

gc()
########################Calculate thresholds for Isotype Controls ##################################

#### Define the isotype controls and extract data for plotting:
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b', 'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')

# Extract the data for isotype controls
isotype_data <- seurat_layer[["ADT"]]$data[isotype_genes, ]
seurat_layer[["ADT"]]$data['Mouse-IgG1',]
rownames(seurat_layer[["ADT"]]$data)

# Convert the matrix to a data frame, retain cell barcodes, and reshape to long format
plot_data <- as.data.frame(t(isotype_data)) %>%
  rownames_to_column("CellBarcode") %>%
  gather(key = "Isotype", value = "Expression", -CellBarcode)

# Add the orig.ident information using the cell barcodes
plot_data$orig.ident <- seurat_layer@meta.data[plot_data$CellBarcode, "orig.ident", drop = TRUE]


#### Calculate the 99% threshold for each isotype control, split by sample:
threshold_data_list <- lapply(isotype_genes, function(isotype) {
  thresholds <- tapply(seurat_layer[["ADT"]]$data[isotype, ], seurat_layer@meta.data$orig.ident, function(x) quantile(x, 0.99))
  data.frame(Isotype = isotype, orig.ident = names(thresholds), Threshold = as.numeric(thresholds))
})

threshold_data <- do.call(rbind, threshold_data_list)

seurat_layer@assays$ADT@data@Dimnames[[1]]
#### Plot using ggplot2 with the added threshold data:

ggplot(plot_data, aes(x = Isotype, y = Expression)) +
  geom_violin(scale = "width", fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) + 
  geom_line(data = threshold_data, aes(y = Threshold, group = orig.ident), color = "red", linetype = "solid", size = 0.5) +
  facet_wrap(~ orig.ident) +  # Split the plot by sample
  theme_bw() +
  labs(title = "Isotype Control Expression with 99% Thresholds", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave('Isotype_Threshold.png',dpi=500, width = 17)

########################################### Isotype Correction #####################################
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b', 'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')

# Read the Excel file
isotype_df <- read_excel(paste0(load.path,"Isotype.xlsx"))

# Create a named vector where names are proteins and values are isotypes
protein_to_isotype <- setNames(isotype_df$isotype, isotype_df$names)

# Calculate the 99% threshold for each isotype
names(isotype_genes) <- isotype_genes
thresholds_list <- lapply(isotype_genes, function(isotype) {
  tapply(seurat_layer[["ADT"]]@data[isotype, ], seurat_layer@meta.data$orig.ident, function(x) quantile(x, 0.99))
})

# Create a copy of the Seurat object to avoid overwriting the original data
seurat_isotype <- seurat_layer

# Apply the thresholds
for (sample in unique(seurat_layer@meta.data$orig.ident)) {
  for (protein in names(protein_to_isotype)) {
    isotype <- protein_to_isotype[protein]
    threshold <- thresholds_list[[isotype]][sample]
    
    # Subtract the threshold
    seurat_isotype[["ADT"]]@data[protein, seurat_layer@meta.data$orig.ident == sample] <- 
      seurat_isotype[["ADT"]]@data[protein, seurat_layer@meta.data$orig.ident == sample] - threshold
    
    # Set negative values to zero
    seurat_isotype[["ADT"]]@data[protein, seurat_isotype[["ADT"]]@data[protein, ] < 0] <- 0
  }
}

save(seurat_isotype, file=paste0(load.path,"Seuratv5_isotype_Assay3.RData"))

