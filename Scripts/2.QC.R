
# Quality Control Visualizations
##Ref -> https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html
##Ref -> http://bioconductor.org/books/3.15/OSCA.intro/getting-scrna-seq-datasets.html
#Ref -> https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html

library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(ggplot2)
library(biomaRt)
library(scCustomize)
library(patchwork)
############################################ Global Variables #######################################################
## Set path to load data
setwd("~/Documents/CD8_Longitudinal/QC")  # Set working directory

# Define input and output paths
in.path <- "~/Documents/10x_Genomics/"
out.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"  # Update as needed

load(paste0(out.path,"Seuratv5_CITEseq_dsbnorm_merged_Seurat.RData"))

############################################## Cell Level QC #########################################################

#Post dsb QC, pre-final filteration

colnames(merged_seurat@meta.data)[4] <- "nCounts_ADT"
colnames(merged_seurat@meta.data)[5] <- "nFeatures_ADT"
colnames(merged_seurat@meta.data)

## Number of Cells per Sample
#The cell numbers can also vary by protocol, producing cell numbers that are much higher than what we loaded. 
#For example, during the inDrops protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets 
#with a single cell and lysis/reaction mixture. While each hydrogel should have a single cellular barcode associated with it, 
#occasionally a hydrogel can have more than one cellular barcode. 
#Similarly, with the 10X protocol there is a chance of obtaining only a barcoded bead in the emulsion droplet (GEM) and no actual cell.

### Create Folder as needed
# Define the folder name
folder_name <- "Pre-QC"

# Check if the folder exists, if not, create it
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  message("Folder '", folder_name, "' created.")
} else {
  message("Folder '", folder_name, "' already exists.")
}

# Change the working directory to the new folder
setwd(folder_name)

metadata <- merged_seurat@meta.data

png(file="Cells_per_sample.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

png(file="Cells_per_Condition.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(x=Condition, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# QC Features rough look

feats.1 <- c("nCount_RNA", "nFeature_RNA","nCounts_ADT", "nFeatures_ADT",
             "percent_mito","percent_ribo","percent_hb","percent_plat")



png(file="Pre-QC_features_grouped.png", width = 1800, height = 1200)
VlnPlot(merged_seurat, group.by = "orig.ident", features = feats.1, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

#UMI Count
#Should generally be above 500, that is the low end of what we expect.
#If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.

png(file="UMI_Count.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()

#nGenes
#For high quality data, the proportional histogram should contain a single large peak that represents cells that were encapsulated. 
#If we see a small shoulder to the left of the major peak (not present in our data), or a bimodal distribution of the cells, 
#that can indicate a couple of things. It might be that there are a set of cells that failed for some reason. 
#It could also be that there are biologically different types of cells (i.e. quiescent cell populations, less complex cells of interest), 
#and/or one type is much smaller than the other (i.e. cells with high counts may be cells that are larger in size).

png(file="nGenes.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250)
dev.off()

#Complexity Score
## Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
#We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score. 
#The novelty score is computed by taking the ratio of nGenes over nUMI. If there are many captured transcripts (high nUMI) 
#and a low number of genes detected in a cell, this likely means that you only captured a low number of genes and simply sequenced transcripts 
#from those lower number of genes over and over again. These low complexity (low novelty) cells could represent a specific cell type 
#(i.e. red blood cells which lack a typical transcriptome), or could be due to an artifact or contamination. 
#Generally, we expect the novelty score to be above 0.80 for good quality cells.

png(file="Complexity_Score.png", width = 1800, height = 1200)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

#Mito Ratio
#This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. 
#We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark, 
#unless of course you are expecting this in your sample.
# Visualize the distribution of mitochondrial gene expression detected per cell

png(file="Mito_Ratio.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_mito, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 15)
dev.off()

# Ribo Ratio

png(file="Ribo_Ratio.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 5)
dev.off()

#Heme Ratio

png(file="Heme_Ratio.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_hb, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)
dev.off()

#Platelet Ratio

png(file="Platlet_Ratio.png", width = 1800, height = 1200)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_plat, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2)
dev.off()

######## SC_Customize Output
# All functions contain
p1 <- QC_Plots_Genes(seurat_object = merged_seurat, low_cutoff = 600, high_cutoff = 5500)
p2 <- QC_Plots_UMIs(seurat_object = merged_seurat, low_cutoff = 1200, high_cutoff = 45000)
p3 <- QC_Plots_Mito(seurat_object = merged_seurat, high_cutoff = 20)
p4 <- QC_Plots_Complexity(seurat_object = merged_seurat, high_cutoff = 0.8)

png(file="Grouped_Cuttoff.png", width = 1800, height = 1200)
wrap_plots(p1, p2, p3, p4, ncol = 4)
dev.off()


### Scatter QC
# All functions contain
png(file="UMIvsGene.png", width = 1800, height = 1200)
QC_Plot_UMIvsGene(seurat_object = merged_seurat, low_cutoff_gene = 600, high_cutoff_gene = 5500, low_cutoff_UMI = 500,
                  high_cutoff_UMI = 50000,group.by = 'orig.ident')
dev.off()

png(file="MitovsGene.png", width = 1800, height = 1200)
QC_Plot_GenevsFeature(seurat_object = merged_seurat, feature1 = "percent_mito", low_cutoff_gene = 600,
                      high_cutoff_gene = 5500, high_cutoff_feature = 20,group.by = 'orig.ident')
dev.off()

png(file="MitovsGene_gradient.png", width = 1800, height = 1200)
QC_Plot_UMIvsGene(seurat_object = merged_seurat, meta_gradient_name = "percent_mito", low_cutoff_gene = 600,
                  high_cutoff_gene = 5500, high_cutoff_UMI = 45000)
dev.off()


#################################################### FILTERING ######################################################

merged_seurat <- JoinLayers(merged_seurat)


DefaultAssay(merged_seurat) <- 'RNA'



# Filter out low quality cells using selected thresholds -> these will change with experiment
##Cell Level Filtering
#nUMI > 500
#nGene > 600
#log10GenesPerUMI > 0.8
#mitoRatio < 20%
#riboratio > 5%
#Heme <20%
#Platlet <2%


filtered_seurat <- subset(x = filtered_seurat, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 600) & 
                            (log10GenesPerUMI > 0.80) & 
                            (percent_mito < 15) &
                            (percent_ribo > 5) &
                            (percent_hb < 20) &
                            (percent_plat < 2)
)


#Within our data we will have many genes with zero counts. 
#These genes can dramatically reduce the average expression for a cell and so we will remove them from our data. 
# Identify genes expressed in at least 10 cells
genes_in_10_cells <- rowSums(filtered_seurat@assays$RNA@layers$counts > 0) >= 10

# Subset the Seurat object to keep only these genes
filtered_seurat <- subset(filtered_seurat, features = names(genes_in_10_cells[genes_in_10_cells]))


### Make Post-QC Folder
setwd("~/Documents/CD8_Longitudinal/QC") 
folder_name <- "Post-QC"

# Check if the folder exists, if not, create it
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  message("Folder '", folder_name, "' created.")
} else {
  message("Folder '", folder_name, "' already exists.")
}

# Change the working directory to the new folder
setwd(folder_name)


# Cells Post QC
png(file="Post-QC_Cells_per_sample.png", width = 1800, height = 1200)
filtered_seurat@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

png(file="Post-QC_features_grouped.png", width = 1800, height = 1200)
VlnPlot(filtered_seurat, group.by = "orig.ident", features = feats.1, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()



f_metadata <- filtered_seurat@meta.data

save(filtered_seurat, file="Seuratv5_filtered_seurat.RData")

