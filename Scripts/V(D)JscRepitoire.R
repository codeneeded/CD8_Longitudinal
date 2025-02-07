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
library(SeuratWrappers)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(scRepertoire)

##set path to load data

setwd('~/Documents/CD8_Longitudinal/VDJ')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"


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
EARTH_names <- f_names[grepl("^S", f_names)]



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
EARTH.TCR <- paste(EARTH_names, ".TCR", sep="")
TARA.BCR <- paste(TARA_names, ".BCR", sep="")
EARTH.BCR <- paste(EARTH_names, ".BCR", sep="")

TARA.contig_list.TCR <- as.list(mget(TARA.TCR))
EARTH.contig_list.TCR <- as.list(mget(EARTH.TCR))
TARA.contig_list.BCR <- as.list(mget(TARA.BCR))
EARTH.contig_list.BCR <- as.list(mget(EARTH.BCR))

#Combine For downstream Analysis

combined.TCR.TARA <- combineTCR(TARA.contig_list.TCR,samples = TARA_names)
combined.TCR.EARTH <- combineTCR(EARTH.contig_list.TCR,samples = EARTH_names)
combined.BCR.TARA <- combineBCR(TARA.contig_list.BCR,samples = TARA_names)
combined.BCR.EARTH <- combineBCR(EARTH.contig_list.BCR,samples = EARTH_names)

### Definitions ###
#we will use clone and define this as the cells with shared/trackable complementarity-determining region 3 (CDR3) sequences. 
#Within this definition, one might use amino acid (aa) sequences of one or both chains to define a clone. 
#Alternatively, we could use nucleotide (nt) or the V(D)JC genes (genes) to define a clone. 
#The latter genes would be a more permissive definition of “clones”, as multiple amino acid or nucleotide sequences can result from the same gene combination. #
#Another option to define clone is the use of the V(D)JC and nucleotide sequence (strict). 

################################################### Basic Clonal Visualizations ################################################
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Clonal_Visualizations'
      )

clonalQuant(combined.TCR.TARA, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)
ggsave('TARA_Unique_Clones_Strict_TRAB_raw.png',width=24,height=12)


clonalQuant(combined.TCR.TARA, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)
ggsave('TARA_Unique_Clones_Strict_TRAB_scaled.png',width=24,height=12)

clonalQuant(combined.TCR.EARTH, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)
ggsave('EARTH_Unique_Clones_Strict_TRAB_raw.png',width=24,height=12)


clonalQuant(combined.TCR.EARTH, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)
ggsave('EARTH_Unique_Clones_Strict_TRAB_scaled.png',width=24,height=12)

clonalQuant(combined.TCR.TARA, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)
ggsave('TARA_Unique_Clones_Strict_TRAB_scaled.png',width=24,height=12)

### Clonal Abundance
clonalAbundance(combined.TCR.TARA, 
                cloneCall = "strict", 
                scale = FALSE)
ggsave('TARA_Clonal_Abundance_TRAB_raw.png',width=16,height=12)


clonalAbundance(combined.TCR.TARA, 
                cloneCall = "strict", 
                scale = TRUE)
ggsave('TARA_Clonal_Abundance_TRAB_scaled.png',width=16,height=12)

clonalAbundance(combined.TCR.EARTH, 
                cloneCall = "strict", 
                scale = FALSE)
ggsave('EARTH_Clonal_Abundance_TRAB_raw.png',width=16,height=12)


clonalAbundance(combined.TCR.EARTH, 
                cloneCall = "strict", 
                scale = TRUE)
ggsave('EARTH_Clonal_Abundance_TRAB_scaled.png',width=16,height=12)

### Clonal Length
clonalLength(combined.TCR.TARA, 
             cloneCall="aa", 
             chain = "both") 
ggsave('TARA_Clonal_Length_TRAB_raw.png',width=16,height=12)


clonalLength(combined.TCR.TARA, 
             cloneCall="aa", 
             chain = "both", 
             scale = TRUE)
ggsave('TARA_Clonal_Abundance_TRAB_scaled.png',width=16,height=12)

clonalLength(combined.TCR.EARTH, 
             cloneCall="aa", 
             chain = "both") 
ggsave('EARTH_Clonal_Length_TRAB_raw.png',width=16,height=12)


clonalLength(combined.TCR.EARTH, 
             cloneCall="aa", 
             chain = "both", 
             scale = TRUE)
ggsave('EARTH_Clonal_Abundance_TRAB_scaled.png',width=16,height=12)

TARA

clonalHomeostasis(combined.TCR.TARA, 
                  cloneCall = "strict",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =
                                  1))
ggsave('TARA_Clonal_Homeostasis_TRAB_scaled.png',width=24,height=12)

clonalHomeostasis(combined.TCR.EARTH, 
                  cloneCall = "strict",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =
                                  1))
ggsave('EARTH_Clonal_Homeostasis_TRAB_scaled.png',width=20,height=12)

######################################## CD3 Composition ##################################
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/CD3_Composition')

percentAA(combined.TCR.TARA, 
          chain = "TRA", 
          aa.length = 20)
ggsave('TARA_Percent_AA_TRA.png',width=26,height=24)

percentAA(combined.TCR.TARA, 
          chain = "TRB", 
          aa.length = 20)
ggsave('TARA_Percent_AA_TRAB.png',width=26,height=24)

percentAA(combined.TCR.EARTH, 
          chain = "TRA", 
          aa.length = 20)
ggsave('EARTH_Percent_AA_TRA.png',width=26,height=24)

percentAA(combined.TCR.EARTH, 
          chain = "TRB", 
          aa.length = 20)
ggsave('EARTH_Percent_AA_TRAB.png',width=26,height=24)

###
positionalEntropy(combined.TCR.TARA, 
                  chain = "both", 
                  aa.length = 20)
ggsave('TARA_Positional_Entropy_TRAB.png',width=20,height=15)

positionalEntropy(combined.TCR.EARTH, 
                  chain = "both", 
                  aa.length = 20)
ggsave('EARTH_Positional_Entropy_TRAB.png',width=20,height=15)

####

vizGenes(combined.TCR.TARA, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('TARA_Heatmap_TRA_V_gene.png',width=12,height=9)

vizGenes(combined.TCR.TARA, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('TARA_Heatmap_TRB_V_gene.png',width=12,height=9)

vizGenes(combined.TCR.TARA, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('TARA_Heatmap_TRA_J_gene.png',width=12,height=9)

vizGenes(combined.TCR.TARA, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('TARA_Heatmap_TRB_J_gene.png',width=12,height=9)

vizGenes(combined.TCR.EARTH, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('EARTH_Heatmap_TRA_V_gene.png',width=12,height=9)

vizGenes(combined.TCR.EARTH, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('EARTH_Heatmap_TRB_V_gene.png',width=12,height=9)

vizGenes(combined.TCR.EARTH, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('EARTH_Heatmap_TRA_J_gene.png',width=12,height=9)

vizGenes(combined.TCR.EARTH, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('EARTH_Heatmap_TRB_J_gene.png',width=12,height=9)



percentKmer(combined.TCR.TARA, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)
ggsave('TARA_Heatmap_TRA_kmer.png',width=12,height=9)

percentKmer(combined.TCR.TARA, 
            cloneCall = "aa",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)
ggsave('TARA_Heatmap_TRB_kmer.png',width=12,height=9)


percentKmer(combined.TCR.EARTH, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)
ggsave('EARTH_Heatmap_TRA_kmer.png',width=12,height=9)

percentKmer(combined.TCR.EARTH, 
            cloneCall = "aa",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)
ggsave('EARTH_Heatmap_TRB_kmer.png',width=12,height=9)


############################### Clonal Diversity and Overlap #####################################################
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Clonal_Diversity')


clonalDiversity(combined.TCR.TARA, 
                cloneCall = "strict")
ggsave('TARA_Clonal_Diversity_Strict.png',width=12,height=9)
clonalDiversity(combined.TCR.EARTH, 
                cloneCall = "strict")
ggsave('EARTH_Clonal_Diversity_Strict.png',width=12,height=9)


clonalOverlap(combined.TCR.TARA, 
              cloneCall = "strict", 
              method = "morisita")
ggsave('TARA_Clonal_Overlap_Strict_morisita.png',width=24,height=11)

clonalOverlap(combined.TCR.TARA, 
              cloneCall = "strict", 
              method = "raw")
ggsave('TARA_Clonal_Overlap_Strict_raw.png',width=24,height=11)

clonalOverlap(combined.TCR.EARTH, 
              cloneCall = "strict", 
              method = "morisita")
ggsave('EARTH_Clonal_Overlap_Strict_morisita.png',width=15,height=11)

clonalOverlap(combined.TCR.EARTH, 
              cloneCall = "strict", 
              method = "raw")
ggsave('EARTH_Clonal_Overlap_Strict_raw.png',width=15,height=11)

setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Clonal_Diversity/By_Sample')

#### Tara - CP003, CP006, CP013, CP018, CP20

clonalCompare(combined.TCR.TARA, 
              top.clones = 20, 
              samples = c("CP003_entry", "CP003_V12","CP003_V24"), 
              order.by = c("CP003_entry", "CP003_V12","CP003_V24"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('TARA_Clonal_Comparison_CP003.png',width=15,height=11)

x <- clonalCompare(combined.TCR.TARA, 
              top.clones = 20, 
              samples = c("CP003_entry", "CP003_V12","CP003_V24"), 
              order.by = c("CP003_entry", "CP003_V12","CP003_V24"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial",
              exportTable = T)
write.csv(x,'TARA_Clonal_Comparison_CP003.csv',row.names = F)

#
clonalCompare(combined.TCR.TARA, 
              top.clones = 20, 
              samples = c("CP006_entry", "CP006_12m","CP006_V24"), 
              order.by = c("CP006_entry", "CP006_12m","CP006_V24"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('TARA_Clonal_Comparison_CP006.png',width=15,height=11)

x <- clonalCompare(combined.TCR.TARA, 
                   top.clones = 20, 
                   samples = c("CP006_entry", "CP006_12m","CP006_V24"), 
                   order.by = c("CP006_entry", "CP006_12m","CP006_V24"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   proportion = F,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'TARA_Clonal_Comparison_CP006.csv',row.names = F)


#
clonalCompare(combined.TCR.TARA, 
              top.clones = 20, 
              samples = c("CP013_1m", "CP013_12m","CP013_24m"), 
              order.by = c("CP013_1m", "CP013_12m","CP013_24m"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('TARA_Clonal_Comparison_CP013.png',width=15,height=11)

x <- clonalCompare(combined.TCR.TARA, 
                   top.clones = 20, 
                   samples = c("CP013_1m", "CP013_12m","CP013_24m"), 
                   order.by = c("CP013_1m", "CP013_12m","CP013_24m"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   proportion = F,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'TARA_Clonal_Comparison_CP013.csv',row.names = F)

#
clonalCompare(combined.TCR.TARA, 
              top.clones = 20, 
              samples = c("CP018_entry", "CP018_V24","CP018_42m"), 
              order.by = c("CP018_entry", "CP018_V24","CP018_42m"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('TARA_Clonal_Comparison_CP018.png',width=15,height=11)

x <- clonalCompare(combined.TCR.TARA, 
                   top.clones = 20, 
                   samples = c("CP018_entry", "CP018_V24","CP018_42m"), 
                   order.by = c("CP018_entry", "CP018_V24","CP018_42m"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   proportion = F,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'TARA_Clonal_Comparison_CP018.csv',row.names = F)


#
clonalCompare(combined.TCR.TARA, 
              top.clones = 20, 
              samples = c("CP020_V1", "CP020_V12", "CP020_V44"), 
              order.by = c("CP020_V1", "CP020_V12", "CP020_V44"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('TARA_Clonal_Comparison_CP020.png',width=15,height=11)

x <- clonalCompare(combined.TCR.TARA, 
                   top.clones = 20, 
                   samples = c("CP020_V1", "CP020_V12", "CP020_V44"), 
                   order.by = c("CP020_V1", "CP020_V12", "CP020_V44"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   proportion = F,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'TARA_Clonal_Comparison_CP020.csv',row.names = F)

### EARTH - SA_AH_004, SA_TY_026, SA_CH_009

clonalCompare(combined.TCR.EARTH, 
              top.clones = 20, 
              samples = c("SA_AH_004_V0", "SA_AH_004_V1"), 
              order.by = c("SA_AH_004_V0", "SA_AH_004_V1"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('EARTH_Clonal_Comparison_SA_AH_004.png',width=15,height=11)

x <- clonalCompare(combined.TCR.EARTH, 
              top.clones = 20, 
              samples = c("SA_AH_004_V0", "SA_AH_004_V1"), 
              order.by = c("SA_AH_004_V0", "SA_AH_004_V1"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial",
              exportTable = T)
write.csv(x,'EARTH_Clonal_Comparison_SA_AH_004.csv',row.names = F)

#
clonalCompare(combined.TCR.EARTH, 
              top.clones = 20, 
              samples = c("SA_TY_026_V2", "SA_TY_026_V5","SA_TY_026_V6","SA_TY_026_V7"), 
              order.by = c("SA_TY_026_V2", "SA_TY_026_V5","SA_TY_026_V6","SA_TY_026_V7"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('EARTH_Clonal_Comparison_SA_TY_026.png',width=15,height=11)

x <- clonalCompare(combined.TCR.EARTH, 
                   top.clones = 20, 
                   samples = c("SA_TY_026_V2", "SA_TY_026_V5","SA_TY_026_V6","SA_TY_026_V7"), 
                   order.by = c("SA_TY_026_V2", "SA_TY_026_V5","SA_TY_026_V6","SA_TY_026_V7"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   proportion = F,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'EARTH_Clonal_Comparison_SA_TY_026.csv',row.names = F)

#
clonalCompare(combined.TCR.EARTH, 
              top.clones = 20, 
              samples = c("SA_CH_009_V1", "SA_CH_009_V5","SA_CH_009_V7","SA_CH_009_V9"), 
              order.by = c("SA_CH_009_V1", "SA_CH_009_V5","SA_CH_009_V7","SA_CH_009_V9"),
              cloneCall="strict", 
              relabel.clones = T,
              proportion = F,
              graph = "alluvial")
ggsave('EARTH_Clonal_Comparison_SA_CH_009.png',width=15,height=11)

x <- clonalCompare(combined.TCR.EARTH, 
                   top.clones = 20, 
                   samples = c("SA_CH_009_V1", "SA_CH_009_V5","SA_CH_009_V7","SA_CH_009_V9"), 
                   order.by = c("SA_CH_009_V1", "SA_CH_009_V5","SA_CH_009_V7","SA_CH_009_V9"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   proportion = F,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'EARTH_Clonal_Comparison_SA_CH_009.csv',row.names = F)

################# Clonal Clustering ############################3




  ############ Merge Seurat #####################
# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(seurat_isotype[[]])

# Use gsub to modify the barcodes, removing everything before and including the fourth '_'
modified_barcodes <- gsub(".*_.*_.*_.*_(.*)", "\\1", barcodes)
modified_barcodes <- paste0(seurat_isotype$orig.ident, "_", modified_barcodes)

# Assign the modified barcodes back to the Seurat object
seurat_isotype <- RenameCells(seurat_isotype, new.names = modified_barcodes)


seurat.tcr <- combineExpression(combined.TCR, 
                                seurat_isotype, 
                                cloneCall="strict",
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)

seurat.bcr <- combineExpression(combined.BCR, 
                                seurat_isotype, 
                                cloneCall="strict",
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)

###################################################### Hyper Expansion Plots ##################################

cols <- c("Hyperexpanded (100 < X <= 500)"="#F0F921","Large (20 < X <= 100)" = "#F69441",
          "Medium (5 < X <= 20)" = "#CA4778" ,"Small (1 < X <= 5)" = "#7D06A5",
          "Single (0 < X <= 1)" = "#0D0887")

## first ordering the Clone Size as a factor, this prevents the coloring from being in alphabetical order. 
slot(seurat.tcr, "meta.data")$cloneSize <- factor(slot(seurat.tcr, "meta.data")$cloneSize, 
                                                  levels = c("Hyperexpanded (100 < X <= 500)", 
                                                             "Large (20 < X <= 100)", 
                                                             "Medium (5 < X <= 20)", 
                                                             "Small (1 < X <= 5)", 
                                                             "Single (0 < X <= 1)", NA))
slot(seurat.bcr, "meta.data")$cloneSize <- factor(slot(seurat.bcr, "meta.data")$cloneSize, 
                                                  levels = c("Hyperexpanded (100 < X <= 500)", 
                                                             "Large (20 < X <= 100)", 
                                                             "Medium (5 < X <= 20)", 
                                                             "Small (1 < X <= 5)", 
                                                             "Single (0 < X <= 1)", NA))


### TCR
setwd("C:/Users/axi313/Documents/TARA_Entry/WNN/VDJ/Seurat_Plots")

tcr.condition <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'Condition', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')

tcr.VL <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'CTLGrouping', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')

tcr.CTL <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'Viral_Load_Category', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')

ggsave('TCR_Seurat_Condition.png', width = 15, dpi = 500, tcr.condition)
ggsave('TCR_Seurat_VL.png', width = 12, dpi = 500, tcr.VL)
ggsave('TCR_Seurat_CTL.png', width = 12, dpi = 500, tcr.CTL)


DimPlot(subset(seurat.tcr, ident='20'), group.by = "cloneSize", reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')



DimPlot(seurat_isotype, label=T, reduction='wnn.umap')
DimPlot(subset(seurat_isotype, ident='singleton'), label=T, reduction='wnn.umap')

### BCR
bcr.condition <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'condition', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')
bcr.stim <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'stim', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')
#bcr.sample <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'sample', reduction = 'wnn.umap') +
#  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')

ggsave('BCR_Seurat_Condition.png', width = 12, dpi = 500, bcr.condition)
ggsave('BCR_Seurat_Stim.png', width = 12, dpi = 500, bcr.stim)
#ggsave('BCR_Seurat_Sample.png', width = 18, dpi = 500, bcr.sample)

### Recluster

seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", 
                               group.singletons=F,
                               algorithm = 3, 
                               resolution = 0.8,
                               random.seed = 1990)

save(seurat_isotype, file=paste0(load.path,"Seuratv5_WNN_Complete_2.RData"))



##############
load(paste0(load.path,'Seuratv5_WNN_labled.RData'))
load(paste0(load.path,'Seuratv5_WNN_labled.RData'))
load(paste0(load.path,'Seuratv5_WNN_labled.RData'))


save(TARA_ALL, file = paste0(load.path, "TARA_ALL_WNN.Rdata"))
save(TARA_HEI, file = paste0(load.path, "TARA_HEI_WNN.Rdata"))
save(EARTH, file = paste0(load.path, "EARTH_WNN.Rdata"))


#### VDJ ####
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

