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
library(igraph)
library(Cairo)
library(RColorBrewer)
library(Polychrome)


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


### For TARA Entry
# --- Step 1: Define VL groups ---
high_VL <- c("CP006","CP022","CP002","CP017","CP016","CP003","CP018")
low_VL  <- c("CP011","SAAH29","CP042","CP013","SATY021")

# --- Step 2: Filter to only _entry samples ---
# TARA: only CP* entry
TARA_entry_CP_names <- TARA_names[grepl("^CP.*_entry$", TARA_names)]
TARA_entry_CP_names

# EARTH: only entry

EARTH_entry_names   <- EARTH_names[grepl("_entry$", EARTH_names)]
EARTH_entry_names

# --- Step 3: Build vector of names ---
Entry_names <- c(TARA_entry_CP_names, EARTH_entry_names)

# --- Step 4: Create contig list (standard method) ---
EntryTCR <- paste(Entry_names, ".TCR", sep = "")
Entry.contig_list.TCR <- as.list(mget(EntryTCR))

# --- Step 5: Combine ---
combined.TCR.entry <- combineTCR(
  Entry.contig_list.TCR,
  samples = Entry_names
)

# --- Step 6: Add Viral Load group as a new variable ---
pids <- sub("_.*", "", Entry_names)
VL_group <- ifelse(pids %in% high_VL, "High_VL",
                   ifelse(pids %in% low_VL,  "Low_VL", "Unknown"))


combined.TCR.entry <- addVariable(
  combined.TCR.entry,
  variable.name = "VL_Group",
  variables = VL_group
)

# (Optional) also add PID as its own variable
combined.TCR.entry <- addVariable(
  combined.TCR.entry,
  variable.name = "PID",
  variables = pids
)



### Definitions ###
#we will use clone and define this as the cells with shared/trackable complementarity-determining region 3 (CDR3) sequences. 
#Within this definition, one might use amino acid (aa) sequences of one or both chains to define a clone. 
#Alternatively, we could use nucleotide (nt) or the V(D)JC genes (genes) to define a clone. 
#The latter genes would be a more permissive definition of “clones”, as multiple amino acid or nucleotide sequences can result from the same gene combination. #
#Another option to define clone is the use of the V(D)JC and nucleotide sequence (strict). 

################################################### Basic Clonal Visualizations ################################################
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Clonal_Visualizations'
      )


clonalQuant(combined.TCR.entry, 
            cloneCall="strict", 
            chain = "both", 
            group.by = 'VL_Group',
            scale = FALSE)
ggsave('TARA_Entry_Unique_Clone_TRAB_HighvsLow.png',width=24,height=12)





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
#raph representation of TCR sequences that are grouped based on their similarity (edit distance threshold of 0.85). 
#The network structure helps identify clonal expansions and relationships between TCR sequences.
#Clustering 
#Nodes (Vertices): Represent unique TCR sequences (or clonotypes).
#Edges (Links): Indicate that two clonotypes are similar (i.e., have an edit distance below the set threshold of 0.85).
#Clusters (Communities): Groups of connected nodes represent TCR sequences that are considered part of the same clonal family.
#Highly connected clusters → Groups of similar TCRs (shared clonal families).
#Larger nodes → Higher clonal expansion.
#Different colors → Different samples, allowing comparison of TCR sharing across conditions.
#Sparse or isolated nodes → Unique or sample-specific clonotypes.
library(igraph)


# Clustering samples for TRA (threshold 0.90)
igraph.TRA <- clonalCluster(
  combined.TCR.EARTH[c("SA_TY_026_V2", "SA_TY_026_V5", "SA_TY_026_V6", "SA_TY_026_V7")],
  chain = "TRA",
  sequence = "aa",
  group.by = "sample",
  threshold = 0.90, 
  exportGraph = TRUE
)

# Clustering samples for TRB (threshold 0.85)
igraph.TRB <- clonalCluster(
  combined.TCR.EARTH[c("SA_TY_026_V2", "SA_TY_026_V5", "SA_TY_026_V6", "SA_TY_026_V7")],
  chain = "TRB",
  sequence = "aa",
  group.by = "sample",
  threshold = 0.85, 
  exportGraph = TRUE
)

# Setting color scheme
custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F", "#E5C494"
)

# Function to set colors and order legend dynamically
set_colors_and_legend <- function(igraph.object) {
  col_legend <- factor(igraph::V(igraph.object)$group)
  color.legend <- unique(igraph::V(igraph.object)$group)
  
  # Order legend from longest to shortest sample names
  ordered_indices <- order(nchar(color.legend), decreasing = TRUE)
  color.legend <- color.legend[ordered_indices]
  
  # Assign distinct colors dynamically based on this order
  num_colors <- length(color.legend)  # Using ordered legend length
  col_samples <- custom_colors[seq_len(num_colors)][ordered_indices]  # Ensure color order matches legend
  
  return(list(col_legend = col_legend, color.legend = color.legend, col_samples = col_samples))
}

# Apply function for TRA and TRB
TRA_colors <- set_colors_and_legend(igraph.TRA)
TRB_colors <- set_colors_and_legend(igraph.TRB)

# Function to plot graph without edge arrows
plot_igraph <- function(igraph.object, col_samples, color.legend, title_text) {
  plot(
    igraph.object,
    vertex.size     = sqrt(igraph::V(igraph.object)$size),
    vertex.label    = NA,
    edge.arrow.size = 0,   # Removes edge arrows (bidirectional)
    edge.curved     = 0.3,
    vertex.color    = col_samples
  )
  
  # Add ordered legend
  legend(
    "topleft", 
    legend = color.legend, 
    pch = 16, 
    col = unique(col_samples), 
    bty = "n"
  )
  
  # Add title
  title(title_text, cex.main = 1.5, font.main = 2)
}

# Set working directory
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Network_Analysis/EARTH')
paste( c("CP020_V1", "CP020_V12", "CP020_V44"), collapse = "_")
# Save TRA plot as a high-quality PNG using Cairo
CairoPNG("SA_TY_026_TRA.png", width = 7200, height = 5500, res = 500)  # Higher resolution
plot_igraph(igraph.TRA, TRA_colors$col_samples, TRA_colors$color.legend, "SA_TY_026 TRA Sequences (Threshold 0.90)")
dev.off()  # Close Cairo device

# Save TRB plot as a high-quality PNG using Cairo
CairoPNG("SA_TY_026_TRB.png", width = 7200, height = 5500, res = 500)
plot_igraph(igraph.TRB, TRB_colors$col_samples, TRB_colors$color.legend, "SA_TY_026 TRB Sequences (Threshold 0.85)")
dev.off()  # Close Cairo device

############################################################ Clonal Clustering #################################################################


# Function to set colors and order legend dynamically
set_colors_and_legend <- function(igraph.object) {
  col_legend <- factor(igraph::V(igraph.object)$group)
  color.legend <- unique(igraph::V(igraph.object)$group)
  
  # Order legend from longest to shortest sample names
  ordered_indices <- order(nchar(color.legend), decreasing = TRUE)
  color.legend <- color.legend[ordered_indices]
  
  # Assign distinct colors dynamically based on this order
  num_colors <- length(color.legend)  # Using ordered legend length
  col_samples <- custom_colors[seq_len(num_colors)][ordered_indices]  # Ensure color order matches legend
  
  return(list(col_legend = col_legend, color.legend = color.legend, col_samples = col_samples))
}

# Function to plot graph without edge arrows
plot_igraph <- function(igraph.object, col_samples, color.legend, title_text) {
  plot(
    igraph.object,
    vertex.size     = sqrt(igraph::V(igraph.object)$size),
    vertex.label    = NA,
    edge.arrow.size = 0,   # Removes edge arrows (bidirectional)
    edge.curved     = 0.3,
    vertex.color    = col_samples
  )
  
  # Add ordered legend
  legend(
    "topleft", 
    legend = color.legend, 
    pch = 16, 
    col = unique(col_samples), 
    bty = "n"
  )
  
  # Add title
  title(title_text, cex.main = 1.5, font.main = 2)
}

# Define color scheme (same for all runs)
custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F", "#E5C494"
)

# Define sample groupings
earth_cohort <- list(
  c("SA_CH_009_V1", "SA_CH_009_V5", "SA_CH_009_V7", "SA_CH_009_V9"),
  c("SA_TY_026_V2", "SA_TY_026_V5", "SA_TY_026_V6", "SA_TY_026_V7"),
  c("SA_AH_004_V0", "SA_AH_004_V1")
)

tara_cohort <- list(
  c("CP020_V1", "CP020_V12", "CP020_V44"),
  c("CP018_entry", "CP018_V24", "CP018_42m"),
  c("CP013_1m", "CP013_12m", "CP013_24m"),
  c("CP006_entry", "CP006_12m", "CP006_V24"),
  c("CP003_entry", "CP003_V12", "CP003_V24")
)

# Function to process and save TRA and TRB graphs for each cohort
process_cohort <- function(cohort_samples, dataset, save_dir) {
  
  # Set working directory for cohort
  setwd(save_dir)
  
  for (samples in cohort_samples) {
    
    # Generate file name prefix from sample names
    sample_label <- unique(sub("(_[A-Za-z0-9]+)?$", "", samples)) # Removes last _ and whatever follows and picks common remnant
    
    # Clustering for TRA (threshold 0.90)
    igraph.TRA <- clonalCluster(
      dataset[samples],
      chain = "TRA",
      sequence = "aa",
      group.by = "sample",
      threshold = 0.90, 
      exportGraph = TRUE
    )
    
    # Clustering for TRB (threshold 0.85)
    igraph.TRB <- clonalCluster(
      dataset[samples],
      chain = "TRB",
      sequence = "aa",
      group.by = "sample",
      threshold = 0.85, 
      exportGraph = TRUE
    )
    
    # Assign colors
    TRA_colors <- set_colors_and_legend(igraph.TRA)
    TRB_colors <- set_colors_and_legend(igraph.TRB)
    
    # Save TRA plot
    CairoPNG(paste0(sample_label, "_TRA.png"), width = 7200, height = 5500, res = 500)
    plot_igraph(igraph.TRA, TRA_colors$col_samples, TRA_colors$color.legend, paste0(sample_label, " TRA Sequences (Threshold 0.90)"))
    dev.off()
    
    # Save TRB plot
    CairoPNG(paste0(sample_label, "_TRB.png"), width = 7200, height = 5500, res = 500)
    plot_igraph(igraph.TRB, TRB_colors$col_samples, TRB_colors$color.legend, paste0(sample_label, " TRB Sequences (Threshold 0.85)"))
    dev.off()
    
    print(paste0("Saved: ", sample_label, " in ", save_dir))
  }
}

# Process EARTH cohort
process_cohort(earth_cohort, combined.TCR.EARTH, "~/Documents/CD8_Longitudinal/VDJ/TCR/Network_Analysis/EARTH")

# Process TARA cohort
process_cohort(tara_cohort, combined.TCR.TARA, "~/Documents/CD8_Longitudinal/VDJ/TCR/Network_Analysis/TARA")

############## Merge with seurat object ###################################

load(paste0(load.path,'TARA_ALL_WNN.Rdata'))
load(paste0(load.path,'TARA_HEI_WNN.Rdata'))
load(paste0(load.path,'EARTH_WNN.Rdata'))

############ Merge Seurat #####################

# TARA All
# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(TARA_ALL[[]])

# Use gsub to modify the barcodes, removing everything before and including the third '_'
modified_barcodes <- gsub(".*_.*_.*_(.*)", "\\1", barcodes)
modified_barcodes <- paste0(TARA_ALL$orig.ident, "_", modified_barcodes)

# Assign the modified barcodes back to the Seurat object
TARA_ALL <- RenameCells(TARA_ALL, new.names = modified_barcodes)


# TARA HEI
# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(TARA_HEI[[]])

# Use gsub to modify the barcodes, removing everything before and including the third '_'
modified_barcodes <- gsub(".*_.*_.*_(.*)", "\\1", barcodes)
modified_barcodes <- paste0(TARA_HEI$orig.ident, "_", modified_barcodes)

# Assign the modified barcodes back to the Seurat object
TARA_HEI <- RenameCells(TARA_HEI, new.names = modified_barcodes)


# EARTH
# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(EARTH[[]])

# Use gsub to modify the barcodes, removing everything before and including the third '_'
modified_barcodes <- gsub(".*_.*_(.*)", "\\1", barcodes)
modified_barcodes <- paste0(EARTH$orig.ident, "_", modified_barcodes)

# Assign the modified barcodes back to the Seurat object
EARTH <- RenameCells(EARTH, new.names = modified_barcodes)


########## Combine TCR Expression with Seurat Object

TARA_ALL <- combineExpression(combined.TCR.TARA, 
                              TARA_ALL, 
                                cloneCall="strict",
                                chain='both',
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)

TARA_HEI <- combineExpression(combined.TCR.TARA, 
                              TARA_HEI, 
                              cloneCall="strict",
                              chain='both',
                              group.by = 'sample',
                              cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                              500),
                              proportion = FALSE)

EARTH <- combineExpression(combined.TCR.EARTH, 
                              EARTH, 
                              cloneCall="strict",
                              chain='both',
                              group.by = 'sample',
                              cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                              500),
                              proportion = FALSE)


######################################################################### Seurat Hyperexpansion Visualisations #################################################

setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Seurat_Plots/Hyperexpansion')


#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "plasma", fixup = TRUE)

DimPlot_scCustom(TARA_ALL, group.by = "cloneSize", reduction = 'wnn.umap') +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave('TARA_ALL_Hyperexpansion.png',width=8,height=7)


DimPlot_scCustom(TARA_ALL, group.by = "cloneSize", reduction = 'wnn.umap',split.by = 'Condition',split_seurat = T) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave('EARTH_Clonal_Comparison_SA_AH_004.png',width=15,height=11)


DimPlot_scCustom(TARA_HEI, group.by = "cloneSize", reduction = 'wnn.umap') +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave('TARA_HEI_Hyperexpansion.png',width=8,height=7)


DimPlot_scCustom(EARTH, group.by = "cloneSize", reduction = 'wnn.umap') +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave('EARTH_Hyperexpansion.png',width=8,height=7)

### Clones per cluster
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Seurat_Plots/Cloness_per_Cluster')

### TARA ALL
clonalOccupy(TARA_ALL, 
             x.axis = "snn.louvianmlr_1",
             label = F)
ggsave('TARA_All_Clonal_Occupancy.png',width=17,height=11)
clonalOccupy(TARA_ALL, 
             x.axis = "snn.louvianmlr_1",
             proportion = T,
             label = F)
ggsave('TARA_All_Clonal_Occupancy_proportion.png',width=17,height=11)

table <- clonalOccupy(TARA_ALL, x.axis = "snn.louvianmlr_1", exportTable = TRUE)
write.csv(table,'TARA_All_Clones_per_Cluster.csv',row.names = F)

clonalOccupy(TARA_ALL, 
             x.axis = "predicted.celltype.l2",
             label = F)
ggsave('TARA_All_Clonal_Occupancy_Azimuth.png',width=27,height=11)
clonalOccupy(TARA_ALL, 
             x.axis = "predicted.celltype.l2",
             proportion = T,
             label = F)
ggsave('TARA_All_Clonal_Occupancy_proportion_Azimuth.png',width=27,height=11)

table <- clonalOccupy(TARA_ALL, x.axis = "predicted.celltype.l2", exportTable = TRUE)
write.csv(table,'TARA_All_Clones_per_Cluster_Azimuth.csv',row.names = F)


### TARA HEI
clonalOccupy(TARA_HEI, 
             x.axis = "snn.louvianmlr_1",
             label = F)
ggsave('TARA_HEI_Clonal_Occupancy.png',width=17,height=11)
clonalOccupy(TARA_HEI, 
             x.axis = "snn.louvianmlr_1",
             proportion = T,
             label = F)
ggsave('TARA_HEI_Clonal_Occupancy_proportion.png',width=17,height=11)

table <- clonalOccupy(TARA_HEI, x.axis = "snn.louvianmlr_1", exportTable = TRUE)
write.csv(table,'TARA_HEI_Clones_per_Cluster.csv',row.names = F)

clonalOccupy(TARA_HEI, 
             x.axis = "predicted.celltype.l2",
             label = F)
ggsave('TARA_HEI_Clonal_Occupancy_Azimuth.png',width=27,height=11)
clonalOccupy(TARA_HEI, 
             x.axis = "predicted.celltype.l2",
             proportion = T,
             label = F)
ggsave('TARA_HEI_Clonal_Occupancy_proportion_Azimuth.png',width=27,height=11)

table <- clonalOccupy(TARA_HEI, x.axis = "predicted.celltype.l2", exportTable = TRUE)
write.csv(table,'TARA_HEIClones_per_Cluster_Azimuth.csv',row.names = F)

### EARTH
clonalOccupy(EARTH, 
             x.axis = "snn.louvianmlr_1",
             label = F)
ggsave('EARTH_Clonal_Occupancy.png',width=17,height=11)
clonalOccupy(EARTH, 
             x.axis = "snn.louvianmlr_1",
             proportion = T,
             label = F)
ggsave('EARTH_Clonal_Occupancy_proportion.png',width=17,height=11)

table <- clonalOccupy(EARTH, x.axis = "snn.louvianmlr_1", exportTable = TRUE)
write.csv(table,'TARA_All_Clones_per_Cluster.csv',row.names = F)

clonalOccupy(EARTH, 
             x.axis = "predicted.celltype.l2",
             label = F)
ggsave('EARTH_Clonal_Occupancy_Azimuth.png',width=27,height=11)
clonalOccupy(EARTH, 
             x.axis = "predicted.celltype.l2",
             proportion = T,
             label = F)
ggsave('EARTH_Clonal_Occupancy_proportion_Azimuth.png',width=27,height=11)

table <- clonalOccupy(EARTH, x.axis = "predicted.celltype.l2", exportTable = TRUE)
write.csv(table,'EARTH_Clones_per_Cluster_Azimuth.csv',row.names = F)

#### Clonal Overlay ####

setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Seurat_Plots/Clonal_Overlay')
#Dr. Francesco Mazziotta and inspired by Drs. Carmona and Andreatta and their work with ProjectTIL,
clonalOverlay(TARA_ALL, 
              reduction = "wnn.umap", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "orig.ident") + 
  guides(color = "none")
ggsave('TARA_Clonal_Overlay_By_Sample.png',width=22,height=17)


clonalOverlay(TARA_ALL, 
              reduction = "wnn.umap", 
              cutpoint = 10, 
              bins = 25, 
              facet.by = "Condition") + 
  guides(color = "none")
ggsave('TARA_Clonal_Overlay_By_Condition.png',width=17,height=9)


clonalOverlay(EARTH, 
              reduction = "wnn.umap", 
              cutpoint = 10, 
              bins = 25, 
              facet.by = "orig.ident") + 
  guides(color = "none")
ggsave('EARTH_Clonal_Overlay_By_Sample.png',width=18,height=13)


#### TCRX ####
library (Trex)
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Trex')

TARA_ALL_TRB_0 <- annotateDB(TARA_ALL, 
                           chains = "TRB")

TARA_ALL_TRB_1 <- annotateDB(TARA_ALL, 
                           chains = "TRB", edit.distance = 1)

TARA_ALL_TRB_2 <- annotateDB(TARA_ALL, 
                           chains = "TRB", edit.distance = 2)


# Extract metadata to dataframe with reordered columns
TARA_TRB_df_0 <- TARA_ALL_TRB_0@meta.data %>%
  dplyr::select(
    cells,
    PID = orig.ident,
    Age,
    Viral_Load,
    cluster = snn.louvianmlr_1,
    CTstrict,
    clonalFrequency,
    TRB_Epitope.target,
    TRB_Epitope.sequence,
    TRB_Epitope.species,
    TRB_Tissue,
    TRB_Cell.type,
    TRB_Database
  ) %>%
  dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)

# For edit distance 1
TARA_TRB_df_1 <- TARA_ALL_TRB_1@meta.data %>%
  dplyr::select(
    cells,
    PID = orig.ident,
    Age,
    Viral_Load,
    cluster = snn.louvianmlr_1,
    CTstrict,
    clonalFrequency,
    TRB_Epitope.target,
    TRB_Epitope.sequence,
    TRB_Epitope.species,
    TRB_Tissue,
    TRB_Cell.type,
    TRB_Database
  ) %>%
  dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)

# For edit distance 2
TARA_TRB_df_2 <- TARA_ALL_TRB_2@meta.data %>%
  dplyr::select(
    cells,
    PID = orig.ident,
    Age,
    Viral_Load,
    cluster = snn.louvianmlr_1,
    CTstrict,
    clonalFrequency,
    TRB_Epitope.target,
    TRB_Epitope.sequence,
    TRB_Epitope.species,
    TRB_Tissue,
    TRB_Cell.type,
    TRB_Database
  ) %>%
  dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)

# Rename TRB columns by edit distance (keeping other columns untouched)
TARA_TRB_df_0_labeled <- TARA_TRB_df_0 %>%
  rename_with(~ paste0(., "_ED0"), starts_with("TRB_"))

TARA_TRB_df_1_labeled <- TARA_TRB_df_1 %>%
  rename_with(~ paste0(., "_ED1"), starts_with("TRB_"))

TARA_TRB_df_2_labeled <- TARA_TRB_df_2 %>%
  rename_with(~ paste0(., "_ED2"), starts_with("TRB_"))

# Define the **key columns** for joining
key_cols <- c("cells", "PID", "Age", "Viral_Load", "cluster", "CTstrict", "clonalFrequency")

# Now **full join them by the key columns**
TARA_TRB_combined <- TARA_TRB_df_0_labeled %>%
  full_join(TARA_TRB_df_1_labeled, by = key_cols) %>%
  full_join(TARA_TRB_df_2_labeled, by = key_cols)


TARA_TRB_combined <- TARA_TRB_combined %>%
  mutate(PID = sub("_.*$", "", PID))



# Pivot TRB_Epitope.target columns to long format
TARA_TRB_long_fixed <- TARA_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  ) %>%
  # Ensure Epitope_Target is labeled
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target)) %>%
  # Remove duplicate entries by ensuring 1 row per CTstrict per Edit Distance
  distinct(CTstrict, Edit_Distance, Epitope_Target, clonalFrequency)

############### Pivot both species and target
# Pivot target
epitope_target_long <- TARA_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  )

# Pivot species
epitope_species_long <- TARA_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.species"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.species_(.*)",
    values_to = "Epitope_Species"
  )

# Join them by cell ID, clone, and edit distance
TARA_TRB_long_species <- epitope_target_long %>%
  select(cells, PID, Age, CTstrict, clonalFrequency, Edit_Distance, Epitope_Target) %>%
  left_join(
    epitope_species_long %>%
      select(cells, CTstrict, Edit_Distance, Epitope_Species),
    by = c("cells", "CTstrict", "Edit_Distance")
  ) %>%
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target),
         Epitope_Species = ifelse(is.na(Epitope_Species), "Unknown", Epitope_Species)) %>%
  distinct()

######

# Plot with CTstrict fill, no legend
ggplot(TARA_TRB_long_fixed, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Remove legend
  ) +
  ylab("Clonal Frequency") +
  xlab("Epitope Target") +
  ggtitle("Clonal Frequency per Epitope Target split by Clone and Edit Distance")


# Plot
# Full Plot 1 – Including All Epitope Targets
ggplot(TARA_TRB_long_fixed, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance, ncol = 1) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 20, l = 30)
  ) +
  ylab("Clonal Frequency") +
  xlab("Epitope Target") +
  ggtitle("Clonal Frequency per Epitope Target split by Clone and Edit Distance")

ggsave('Tara_TRB_ClonalFreq_vs_Target.png', width = 18, height = 13, dpi = 300, bg = 'white')

# Remove Unknown Targets
TARA_TRB_long_no_unknown <- TARA_TRB_long_fixed %>%
  filter(Epitope_Target != "Unknown")

# Full Plot 2 – Excluding Unknown Epitope Targets
ggplot(TARA_TRB_long_no_unknown, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance, ncol = 1) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 20, l = 30)
  ) +
  ylab("Clonal Frequency") +
  xlab("Epitope Target") +
  ggtitle("Clonal Frequency per Epitope Target split by Clone and Edit Distance (Excluding Unknown)")

ggsave('Tara_TRB_ClonalFreq_vs_Target_Unkownexclude.png', width = 20, height = 13, dpi = 300, bg = 'white')

write_csv(TARA_TRB_combined,'Trex_TRB_Epitope_Database.csv')

# Redo the long format while retaining PID, Age and other metadata
TARA_TRB_long_fixed <- TARA_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  ) %>%
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target)) %>%
  distinct(PID, Age, CTstrict, Edit_Distance, Epitope_Target, clonalFrequency)

# Now summarize by Edit Distance, PID, Age, and Epitope Target
TARA_TRB_summary_by_PID_Age <- TARA_TRB_long_fixed %>%
  group_by(Edit_Distance, PID, Age, Epitope_Target) %>%
  summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE), .groups = "drop") %>%
  arrange(Edit_Distance, PID, Age, desc(Total_Clonal_Frequency))

# View the table
write_csv(TARA_TRB_summary_by_PID_Age,'Trex_TARA_Epitope_Specificity_By_Sample.csv')

### Pie charts #######
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Trex/PID_Pie_Chart')

species_pie_data_full <- TARA_TRB_long_species %>%
  group_by(PID, Age, Edit_Distance, Epitope_Species) %>%
  summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE), .groups = "drop")

# Normalize clonal frequency within each PID x Age x Edit_Distance group
species_pie_data_normalized <- species_pie_data_full %>%
  group_by(PID, Age, Edit_Distance) %>%
  mutate(Percent = Total_Clonal_Frequency / sum(Total_Clonal_Frequency)) %>%
  ungroup()

species_pie_data_normalized <- species_pie_data_normalized %>%
  mutate(Age = paste0(Age, ifelse(Age == 1, " Month", " Months")))

# Set output directory
output_dir <- "~/Documents/CD8_Longitudinal/VDJ/TCR/Trex/PID_Pie_Chart"
setwd(output_dir)


# Get unique PIDs
pid_list <- levels(as.factor(species_pie_data_normalized$PID))



# Generate short, single-line labels
species_pie_data_normalized <- species_pie_data_normalized %>%
  mutate(Epitope_Species_short = sapply(strsplit(Epitope_Species, ";"), function(x) {
    if (length(x) > 2) paste0(x[1], ";", x[2], "...")
    else paste(x, collapse = ";")
  }))

# Use short label for fill
all_labels <- unique(species_pie_data_normalized$Epitope_Species_short)
n_colors <- length(all_labels)
palette <- Polychrome::createPalette(
  N = n_colors,
  seedcolors = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
)

names(palette) <- all_labels

# PIDs
pid_list <- levels(as.factor(species_pie_data_normalized$PID))

for (pid_to_plot in pid_list) {
  
  plot_data <- species_pie_data_normalized %>%
    filter(PID == pid_to_plot) %>%
    droplevels()
  
  if (nrow(plot_data) == 0) next
  
  n_rows <- length(unique(plot_data$Age))
  n_cols <- length(unique(plot_data$Edit_Distance))
  plot_width <- n_cols * 5
  plot_height <- n_rows * 5
  
  p <- ggplot(plot_data, aes(x = "", y = Percent, fill = Epitope_Species_short)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    facet_grid(rows = vars(Age), cols = vars(Edit_Distance), switch = "y") +
    scale_fill_manual(values = palette) +
    theme_void(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),                 # Increased from 10 → 12
      legend.key.height = unit(0.6, "lines"),                # Slightly taller blocks
      strip.text.y.left = element_text(angle = 0, size = 14, face = "bold"),
      strip.text.x = element_text(size = 14, face = "bold", margin = margin(t = 6)),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      plot.margin = margin(t = 12, r = 10, b = 24, l = 10)    # More space below for legend
    ) +
    guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
    ggtitle(paste("Epitope Species Specificity for", pid_to_plot))
  
  ggsave(
    filename = paste0(pid_to_plot, "_pie.png"),
    plot = p,
    bg = "white",
    width = plot_width,
    height = plot_height,
    dpi = 300
  )
}

### EDIT PLOTS FOR EARTH COHORT TOO! (Dont Forget)
##############################
### Earth
setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Trex')

EARTH_TRB_0 <- annotateDB(EARTH, 
                             chains = "TRB")

EARTH_TRB_1 <- annotateDB(EARTH, 
                             chains = "TRB", edit.distance = 1)

EARTH_TRB_2 <- annotateDB(EARTH, 
                             chains = "TRB", edit.distance = 2)


# Extract metadata to dataframe with reordered columns
EARTH_TRB_df_0 <- EARTH_TRB_0@meta.data %>%
  dplyr::select(
    cells,
    PID = orig.ident,
    Age,
    Viral_Load,
    cluster = snn.louvianmlr_1,
    CTstrict,
    clonalFrequency,
    TRB_Epitope.target,
    TRB_Epitope.sequence,
    TRB_Epitope.species,
    TRB_Tissue,
    TRB_Cell.type,
    TRB_Database
  ) %>%
  dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)

# For edit distance 1
EARTH_TRB_df_1 <- EARTH_TRB_1@meta.data %>%
  dplyr::select(
    cells,
    PID = orig.ident,
    Age,
    Viral_Load,
    cluster = snn.louvianmlr_1,
    CTstrict,
    clonalFrequency,
    TRB_Epitope.target,
    TRB_Epitope.sequence,
    TRB_Epitope.species,
    TRB_Tissue,
    TRB_Cell.type,
    TRB_Database
  ) %>%
  dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)

# For edit distance 2
EARTH_TRB_df_2 <- EARTH_TRB_2@meta.data %>%
  dplyr::select(
    cells,
    PID = orig.ident,
    Age,
    Viral_Load,
    cluster = snn.louvianmlr_1,
    CTstrict,
    clonalFrequency,
    TRB_Epitope.target,
    TRB_Epitope.sequence,
    TRB_Epitope.species,
    TRB_Tissue,
    TRB_Cell.type,
    TRB_Database
  ) %>%
  dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)

# Rename TRB columns by edit distance (keeping other columns untouched)
EARTH_TRB_df_0_labeled <- EARTH_TRB_df_0 %>%
  rename_with(~ paste0(., "_ED0"), starts_with("TRB_"))

EARTH_TRB_df_1_labeled <- EARTH_TRB_df_1 %>%
  rename_with(~ paste0(., "_ED1"), starts_with("TRB_"))

EARTH_TRB_df_2_labeled <- EARTH_TRB_df_2 %>%
  rename_with(~ paste0(., "_ED2"), starts_with("TRB_"))

# Define the **key columns** for joining
key_cols <- c("cells", "PID", "Age", "Viral_Load", "cluster", "CTstrict", "clonalFrequency")

# Now **full join them by the key columns**
EARTH_TRB_combined <- EARTH_TRB_df_0_labeled %>%
  full_join(EARTH_TRB_df_1_labeled, by = key_cols) %>%
  full_join(EARTH_TRB_df_2_labeled, by = key_cols)

levels(as.factor(EARTH_TRB_combined$PID))

EARTH_TRB_combined <- EARTH_TRB_combined %>%
  mutate(PID = sub("_(V\\d+|entry)$", "", PID))

# Pivot TRB_Epitope.target columns to long format
EARTH_TRB_long_fixed <- EARTH_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  ) %>%
  # Ensure Epitope_Target is labeled
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target)) %>%
  # Remove duplicate entries by ensuring 1 row per CTstrict per Edit Distance
  distinct(CTstrict, Edit_Distance, Epitope_Target, clonalFrequency)


# Plot with CTstrict fill, no legend
ggplot(EARTH_TRB_long_fixed, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Remove legend
  ) +
  ylab("Clonal Frequency") +
  xlab("Epitope Target") +
  ggtitle("Clonal Frequency per Epitope Target split by Clone and Edit Distance")


# Plot
ggplot(EARTH_TRB_long_fixed, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance, ncol = 1) +  # Stack facets vertically
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  ylab("Clonal Frequency") +
  xlab("Epitope Target") +
  ggtitle("Clonal Frequency per Epitope Target split by Clone and Edit Distance")
ggsave('EARTH_TRB_ClonalFreq_vs_Target.png',width=18,height=13)

EARTH_TRB_long_no_unknown <- EARTH_TRB_long_fixed %>%
  filter(Epitope_Target != "Unknown")

ggplot(EARTH_TRB_long_no_unknown, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance, ncol = 1) +  # Stack facets vertically
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  ylab("Clonal Frequency") +
  xlab("Epitope Target") +
  ggtitle("Clonal Frequency per Epitope Target split by Clone and Edit Distance (Excluding Unknown)")
ggsave('EARTH_TRB_ClonalFreq_vs_Target_Unkownexclude.png',width=18,height=13)

write_csv(EARTH_TRB_combined,'Earth_Trex_TRB_Epitope_Database.csv')

# Redo the long format while retaining PID, Age and other metadata
EARTH_TRB_long_fixed <- EARTH_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  ) %>%
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target)) %>%
  distinct(PID, Age, CTstrict, Edit_Distance, Epitope_Target, clonalFrequency)

# Now summarize by Edit Distance, PID, Age, and Epitope Target
EARTH_TRB_summary_by_PID_Age <- EARTH_TRB_long_fixed %>%
  group_by(Edit_Distance, PID, Age, Epitope_Target) %>%
  summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE), .groups = "drop") %>%
  arrange(Edit_Distance, PID, Age, desc(Total_Clonal_Frequency))

# View the table
write_csv(EARTH_TRB_summary_by_PID_Age,'Trex_EARTH_Epitope_Specificity_By_Sample.csv')

#### Per i ndividual pie chart with and without unkown, specificity to different disease pathogens
### Maintained clones per timepoint known specificity
### Unkown specificity
### Future experiment summary
### Pie charts #######


############### Pivot both species and target
# Pivot target
epitope_target_long <- EARTH_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  )

# Pivot species
epitope_species_long <- EARTH_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.species"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.species_(.*)",
    values_to = "Epitope_Species"
  )

# Join them by cell ID, clone, and edit distance
EARTH_TRB_long_species <- epitope_target_long %>%
  select(cells, PID, Age, CTstrict, clonalFrequency, Edit_Distance, Epitope_Target) %>%
  left_join(
    epitope_species_long %>%
      select(cells, CTstrict, Edit_Distance, Epitope_Species),
    by = c("cells", "CTstrict", "Edit_Distance")
  ) %>%
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target),
         Epitope_Species = ifelse(is.na(Epitope_Species), "Unknown", Epitope_Species)) %>%
  distinct()



setwd('~/Documents/CD8_Longitudinal/VDJ/TCR/Trex/PID_Pie_Chart')

species_pie_data_full <- EARTH_TRB_long_species %>%
  group_by(PID, Age, Edit_Distance, Epitope_Species) %>%
  summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE), .groups = "drop")

# Normalize clonal frequency within each PID x Age x Edit_Distance group
species_pie_data_normalized <- species_pie_data_full %>%
  group_by(PID, Age, Edit_Distance) %>%
  mutate(Percent = Total_Clonal_Frequency / sum(Total_Clonal_Frequency)) %>%
  ungroup()

species_pie_data_normalized <- species_pie_data_normalized %>%
  mutate(Age = paste0(Age, ifelse(Age == 1, " Month", " Months")))

# Set output directory
output_dir <- "~/Documents/CD8_Longitudinal/VDJ/TCR/Trex/PID_Pie_Chart"
setwd(output_dir)


# Get unique PIDs
pid_list <- levels(as.factor(species_pie_data_normalized$PID))

# Loop through each PID and save as PNG
for (pid_to_plot in pid_list) {
  
  plot_data <- species_pie_data_normalized %>%
    filter(PID == pid_to_plot) %>%
    droplevels()
  
  # Skip if no data
  if (nrow(plot_data) == 0) next
  
  # Generate plot
  p <- ggplot(plot_data, aes(x = "", y = Percent, fill = Epitope_Species)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    facet_grid(rows = vars(Age), cols = vars(Edit_Distance), switch = "y") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(paste("Epitope Species Specificity for", pid_to_plot))
  
  # Save as PNG
  ggsave(
    filename = paste0(pid_to_plot, "_pie.png"),
    plot = p,
    bg = "white",
    width = 21,
    height = 14
  )
}


####################################### Viral Load ####################################################

TARA_clonal_expansion_summary <- TARA_TRB_combined %>%
  group_by(PID, Age) %>%
  summarise(
    Viral_Load = unique(Viral_Load),  # One per timepoint assumed
    TARA_Expanded_Clone_Count = sum(clonalFrequency > 1, na.rm = TRUE),
    TARA_Total_Clone_Count = n(),
    .groups = "drop"
  )
TARA_clonal_expansion_summary_2 <- TARA_clonal_expansion_summary %>%
  mutate(
    Viral_Load_Numeric = as.numeric(Viral_Load),
    Viral_Load_Log = log10(Viral_Load_Numeric + 1),
    Viral_Load_Binned = cut(
      Viral_Load_Log,
      breaks = c(-Inf, 1, 2, 3, 4, Inf),
      labels = c("Very Low", "Low", "Medium", "High", "Very High")
    ),
    Age = paste0(Age, ifelse(Age == 1, " Month", " Months")),
    Age = factor(Age, levels = c("1 Month", "12 Months", "24 Months"))
  ) %>%
  filter(!is.na(Viral_Load_Binned), !is.na(TARA_Expanded_Clone_Count))

EARTH_clonal_expansion_summary <- EARTH_TRB_combined %>%
  group_by(PID, Age) %>%
  summarise(
    Viral_Load = unique(Viral_Load),
    EARTH_Expanded_Clone_Count = sum(clonalFrequency > 1, na.rm = TRUE),
    EARTH_Total_Clone_Count = n(),
    .groups = "drop"
  )

ggplot(TARA_clonal_expansion_summary, aes(x = Viral_Load, y = TARA_Expanded_Clone_Count)) +
  geom_point(aes(color = Age)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~PID) +
  theme_minimal() +
  labs(
    title = "TARA: Clonal Expansion vs Viral Load",
    x = "Viral Load",
    y = "Expanded Clones (clonalFrequency > 1)"
  )
ggplot(EARTH_clonal_expansion_summary, aes(x = Viral_Load, y = EARTH_Expanded_Clone_Count)) +
  geom_point(aes(color = Age)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~PID) +
  theme_minimal() +
  labs(
    title = "EARTH: Clonal Expansion vs Viral Load",
    x = "Viral Load",
    y = "Expanded Clones (clonalFrequency > 1)"
  )

ggplot(TARA_clonal_expansion_summary_2, aes(x = Viral_Load_Log, y = TARA_Expanded_Clone_Count, color = PID, shape = Age)) +
  geom_point(size = 3) +
  geom_line(aes(group = PID), linewidth = 1) +
  theme_minimal() +
  labs(
    title = "TARA: Log Viral Load vs Clonal Expansion (Age as Shape)",
    x = "log10(Viral Load + 1)",
    y = "Expanded Clones (clonalFrequency > 1)"
  ) +
  theme(legend.title = element_blank())

ggplot(TARA_clonal_expansion_summary_2, aes(x = Age, y = TARA_Expanded_Clone_Count, color = PID, shape = Viral_Load_Binned)) +
  geom_point(size = 3) +
  geom_line(aes(group = PID), linewidth = 1) +
  theme_minimal() +
  labs(
    title = "TARA: Clonal Expansion Over Time (Binned Viral Load as Shape)",
    x = "Age",
    y = "Expanded Clones (clonalFrequency > 1)",
    shape = "log10(Viral Load + 1)"
  ) +
  theme(
    legend.title = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#### Save TCR Combined Seurat ####

save(TARA_ALL, TARA_ALL_TRB_0, TARA_ALL_TRB_1, TARA_ALL_TRB_2,
     file = "/home/akshay-iyer/Documents/CD8_Longitudinal/TARA_TCR_Combined.RData")

save(EARTH, EARTH_TRB_0, EARTH_TRB_1, EARTH_TRB_2,
     file = "/home/akshay-iyer/Documents/CD8_Longitudinal/EARTH_TCR_Combined.RData")
