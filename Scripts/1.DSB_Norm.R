###Citeseq Pipeline
#Cite seurat and ds packages
#Load Required Libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(ggplot2)
library(hdf5r)
###############################################Path + GLobal Variables ##################################################
## Set path to load data
setwd("~/Documents/10x_Genomics")  # Set working directory

# Define input and output paths
in.path <- "~/Documents/10x_Genomics/"
out.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"  # Update as needed
r.str <- "/multi/GEX_FB/raw_feature_bc_matrix.h5"
f.str <- "/per_sample_outs/GEX_FB/sample_filtered_feature_bc_matrix.h5"

# Function to get all folder names within a specified path
get_folder_names <- function(in.path) {
  # List all directories within the specified path without recursion
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  
  # Filter out the root path itself if present
  folder_names <- folder_names[folder_names != ""]
  
  return(folder_names)
}

# List folder names

f_names <- get_folder_names(in.path)

# Print the folder names
print(f_names)

####################################### Load Data + Basic Pre-processing ##############################################################
# Standard protein names (desired output)
standard_protein_names <- c(
  "CD86", "CD274", "TNFRSF14", "PVR", "NECTIN2", "CD47", "CD48", "CD40", "CD40LG", "CD52", 
  "CD3D", "CD8A", "NCAM1", "CD19", "CD33", "ITGAX", "HLA-A", "CD45RA", "IL3RA", "CD7", 
  "ENG", "ITGA6", "CCR4", "CD4", "CD44", "CD14", "FCGR3A", "IL2RA", "CD45RO", "PDCD1", 
  "TIGIT", "Mouse-IgG1", "Mouse-IgG2a", "Mouse-IgG2b", "Rat-IgG2b", "MS4A1", "NCR1", 
  "PECAM1", "MCAM", "IGHM", "CD5", "CXCR3", "CCR5", "FCGR2A", "CCR6", "CXCR5", "ITGAE", 
  "CD69", "SELL", "KLRB1", "CTLA4", "LAG3", "KLRG1", "CD27", "LAMP1", "FAS", "TNFRSF4", 
  "HLA-DRA", "CD1C", "ITGAM", "FCGR1A", "THBD", "CD1D", "KLRK1", "CR1", "B3GAT1", "BTLA", 
  "ICOS", "CD58", "ENTPD1", "CX3CR1", "CD24", "CR2", "ITGAL", "CD79B", "CD244", "SIGLEC1", 
  "ITGB7", "TNFRSF13C", "GP1BB", "ICAM1", "SELP", "IFNGR1", "TCR-AB", "Rat-IgG1", 
  "Rat-IgG2a", "Hamster-IgG", "IL2RB", "TNFRSF13B", "FCER1A", "ITGA2B", "TNFRSF9", 
  "CD163", "CD83", "IL4R", "ANPEP", "CD2", "CD226", "ITGB1", "CLEC4C", "ITGA2", "CD81", 
  "IGHD", "ITGB2", "CD28", "CD38", "IL7R", "CD45", "CD22", "TFRC", "DPP4", "CD36", 
  "KIR2DL1", "ITGA1", "ITGA4", "NT5E", "TCR-vA7.2", "TCR-vD2", "OLR1", "KIR2DL3", 
  "KIR3DL1", "SLAMF7", "CD99", "CLEC12A", "SLAMF6", "KLRD1", "IGKC", "LILRB1", "FCER2", 
  "Human-Ig", "SIGLEC7", "ADGRG1", "HLA-E", "CD82", "CD101", "C5AR1", "GGT1"
)

# Loop through the files and standardize protein names
for (i in f_names) {
  # Read data from 10x outputs
  raw <- Read10X_h5(paste0(in.path, i, r.str))
  cells <- Read10X_h5(paste0(in.path, i, f.str))
  
  # Replace antibody capture names with standardized names based on position
  raw$`Antibody Capture`@Dimnames[[1]] <- standard_protein_names
  cells$`Antibody Capture`@Dimnames[[1]] <- standard_protein_names
  
  # define a vector of cell-containing barcodes and remove them from unfiltered data 
  stained_cells <- colnames(cells$`Gene Expression`)
  background <- setdiff(colnames(raw$`Gene Expression`), stained_cells)
  
  # Assign final processed data
  prot <- raw$`Antibody Capture`
  rna <- raw$`Gene Expression`
  
  # Create metadata
  rna.size <- log10(Matrix::colSums(rna))
  prot.size <- log10(Matrix::colSums(prot))
  nCount_RNA <- Matrix::colSums(rna) 
  nCount_ADT <- Matrix::colSums(prot) 
  nFeature_RNA <- Matrix::colSums(rna > 0) 
  nFeature_ADT <- Matrix::colSums(prot > 0) 
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  mt.prop <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  
  # Combine metadata into a data frame
  md <- as.data.frame(cbind(nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, rna.size, prot.size, mt.prop))
  
  # add indicator for barcodes Cell Ranger called as cells
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # remove barcodes with no evidence of capture in the experiment
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  # Save metadata and matrices
  assign(paste0(i, ".md"), md)
  assign(paste0(i, ".prot"), prot)
  assign(paste0(i, ".rna"), rna)
}

# Example: Check protein names
print(`SATY021_entry.prot`@Dimnames[[1]])


###  Droplet Settings
# Define the base path
base_path <- "~/Documents/CD8_Longitudinal"

# Define the directories to create
qc_dir <- file.path(base_path, "QC")
droplet_settings_dir <- file.path(qc_dir, "Droplet_Settings")

# Create directories if they do not exist
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
  cat("Created directory:", qc_dir, "\n")
}

if (!dir.exists(droplet_settings_dir)) {
  dir.create(droplet_settings_dir, recursive = TRUE)
  cat("Created directory:", droplet_settings_dir, "\n")
}

for (i in f_names) {
  # Set metadata object and protein object
  md <- eval(parse(text = paste0(i, ".md")))
  prot <- eval(parse(text = paste0(i, ".prot")))
  rna <- eval(parse(text = paste0(i, ".rna")))
  
  # Output plot for detected genes vs protein library size
  png(file.path(droplet_settings_dir, paste0(i, "_genevsprotlibsize.png")),
      width = 800, height = 600)
  
  p <- ggplot(md, aes(x = log10(nFeature_RNA), y = prot.size)) + 
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~drop.class)
  
  print(p)
  dev.off()
}

############################################### Droplet Thresholds for DSB Norm###################################################

# Define the function to process each entry
process_entry <- function(entry_name) {
  # Dynamically create the variable name for the data frame
  md_var_name <- paste0(entry_name, ".md")
  
  # Get the data frame from the variable name
  md <- get(md_var_name)
  
  # Filter the rows based on the given conditions
  filtered_md <- md[md$prot.size > 0.5 & md$prot.size < 4 & md$rna.size < 2.5, ]
  
  # Extract the row names of the filtered rows
  background_drops <- rownames(filtered_md)
  
  # Dynamically assign the result to the new variable
  assign(paste0(entry_name, ".background_drops"), background_drops, envir = .GlobalEnv)
}


# Loop through each entry and process it
for (entry_name in f_names) {
  process_entry(entry_name)
}


###################################Initial QC + DSB Normalization###################################################


for (i in f_names) {
  #Set metadata object and protein object to prevent constant eval parse calling  
  md <- eval(parse(text= paste0(i,'.md')))
  prot <- eval(parse(text = paste0(i,'.prot')))
  rna <- eval(parse(text = paste0(i,'.rna')))
  background_drops <- eval(parse(text = paste0(i,'.background_drops')))
  
  background.adt.mtx = as.matrix(prot[ , background_drops])
  
  cellmd = md[md$drop.class == 'cell', ]  #Define Cell Metadata on only cells not background  
  
  # filter drops with + / - 3 median absolute deviations from the median library size
  rna.mult = (3*mad(cellmd$rna.size))
  prot.mult = (3*mad(cellmd$prot.size))
  rna.lower = median(cellmd$rna.size) - rna.mult
  rna.upper = median(cellmd$rna.size) + rna.mult
  prot.lower = median(cellmd$prot.size) - prot.mult
  prot.upper = median(cellmd$prot.size) + prot.mult
  
  # filter rows based on droplet qualty control metrics
  qc_cells = rownames(
    cellmd[cellmd$prot.size > prot.lower & 
             cellmd$prot.size < prot.upper & 
             cellmd$rna.size > rna.lower & 
             cellmd$rna.size < rna.upper & 
             cellmd$mt.prop < 0.25, ]
  )
  
  # Output thresholds for quality control metrics as in any standard scRNAseq analysis
  png(paste0(i,'_qc_thresholds.png'),width = 800, height = 600)
  plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
  p1 = ggplot(cellmd, aes(x = rna.size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
  p2 = ggplot(cellmd, aes(x = mt.prop)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
  p3 = ggplot(cellmd, aes(x = log10(nFeature_RNA), y = rna.size, fill = mt.prop )) + plot_aes
  p4 = ggplot(cellmd, aes(x = nFeature_RNA, y = prot.size, fill = mt.prop )) + plot_aes
  print (p1+p2+p3+p4)
  dev.off()
  cell.adt.raw = as.matrix(prot[ , qc_cells])
  cell.rna.raw = rna[ ,qc_cells]
  cellmd = cellmd[qc_cells, ]
  
  #Proteins without Staining
  pm = sort(apply(cell.adt.raw, 1, max))
  pm2 = apply(background.adt.mtx, 1, max)
  head(pm2)
  
  #Assign it to your final output
  assign(paste0(i,'.cell.adt.raw'),cell.adt.raw)
  assign(paste0(i,'.cell.rna.raw'),cell.rna.raw)
  assign(paste0(i,'.background.adt.mtx'),background.adt.mtx)
  assign(paste0(i,'.cellmd'),cellmd)
  assign(paste0(i,'.pm'),pm)
}


# Check if you need to remove proteins without staining
#https://www.rdocumentation.org/packages/dsb/versions/0.3.0
#prot.expres.total <- rbindlist(adt.list)

#In this case we do not

### DSB Normalisation
#Set isotype control

isotype.controls <- c('Mouse-IgG1', 'Mouse-IgG2a','Mouse-IgG2b', 'Rat-IgG2b'
                      ,'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')


#normalize protein data for the cell containing droplets with the dsb method. 

for (i in f_names) {
  cell.adt.raw <- eval(parse(text = paste0(i,'.cell.adt.raw')))
  background.adt.mtx <- eval(parse(text = paste0(i,'.background.adt.mtx')))
  # normalize and denoise with dsb with 
  cells.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.raw, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls
  )
  cells.dsb.norm <- Matrix(as.matrix(cells.dsb.norm),sparse=TRUE)
  assign(paste0(i,'.cells.dsb.norm'),cells.dsb.norm)
}


########################################Create + Merge Seurat Object (norm dsb+RNA)########################################

# Create Seurat Object

for (i in f_names) {
  cellmd <- eval(parse(text= paste0(i,'.cellmd')))
  cell.adt.raw <- eval(parse(text= paste0(i,'.cell.adt.raw')))
  cells.dsb.norm <- eval(parse(text= paste0(i,'.cells.dsb.norm')))
  cell.rna.raw <- eval(parse(text= paste0(i,'.cell.rna.raw')))
  
  # integrating with Seurat
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))
  
  # create Seurat object note: min.cells is a gene filter, not a cell filter
  s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                                 meta.data = cellmd,
                                 assay = "RNA", 
                                 min.cells = 20)
  
  # add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
  s[["ADT"]] <- CreateAssayObject(data = cells.dsb.norm)
  s$orig.ident <- i
  #Assign
  assign(paste0(i,'.seurat'),s)
}


#Merge Seurat Objects

# Dynamically construct the merge() function
merged_seurat <- merge(
  x = eval(parse(text = paste0(f_names[1], ".seurat"))), 
  y = lapply(f_names[-1], function(name) eval(parse(text = paste0(name, ".seurat")))), 
  add.cell.id = f_names
)

###################################################### Edit Seurat Metadata ###############################################


### Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

### Compute percent mitochondrial genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent_mito")

### Compute percentage of ribosomal genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
merged_seurat <- PercentageFeatureSet(merged_seurat, "^HB[^(P)]", col.name = "percent_hb")

#Same for platelets
merged_seurat <- PercentageFeatureSet(merged_seurat, "PECAM1|PF4", col.name = "percent_plat")


## Create a metadata Object from seurat to work on

### Create metadata dataframe
metadata <- merged_seurat@meta.data

### Add cell IDs to metadata
metadata$cells <- rownames(metadata)

### Create Condition Column
metadata$Condition <- NA
metadata$Condition[which(str_detect(metadata$cells, "CE"))] <- "HEU"
metadata$Condition[which(str_detect(metadata$cells, "CP"))] <- "HEI"
metadata$Condition[which(str_detect(metadata$cells, "CS"))] <- "HUU"
metadata$Condition[which(str_detect(metadata$cells, "SA"))] <- "HEI"


### Create Timepoint Column
# Initialize Age column in metadata
metadata$Age <- NA

# Assign age for each entry based on f_names
metadata$Age[which(str_detect(metadata$cells, "SA_CH_009_V1"))] <- 1.5
metadata$Age[which(str_detect(metadata$cells, "SA_CH_009_V5"))] <- 11.6
metadata$Age[which(str_detect(metadata$cells, "SA_CH_009_V7"))] <- 22.7
metadata$Age[which(str_detect(metadata$cells, "SA_CH_009_V9"))] <- 33.9
metadata$Age[which(str_detect(metadata$cells, "SA_TY_026_V2"))] <- 3.6
metadata$Age[which(str_detect(metadata$cells, "SA_TY_026_V6"))] <- 20
metadata$Age[which(str_detect(metadata$cells, "CP003_V12"))] <- 12
metadata$Age[which(str_detect(metadata$cells, "CP003_V24"))] <- 24
metadata$Age[which(str_detect(metadata$cells, "CP006_V24"))] <- 24
metadata$Age[which(str_detect(metadata$cells, "CP018_V24"))] <- 24
metadata$Age[which(str_detect(metadata$cells, "CP020_V1"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP020_V12"))] <- 12
metadata$Age[which(str_detect(metadata$cells, "CP020_V44"))] <- 44
metadata$Age[which(str_detect(metadata$cells, "SA_AH_004_V0"))] <- 2.3
metadata$Age[which(str_detect(metadata$cells, "SA_AH_004_V1"))] <- 2.8
metadata$Age[which(str_detect(metadata$cells, "SA_TY_026_V5"))] <- 14
metadata$Age[which(str_detect(metadata$cells, "SA_TY_026_V7"))] <- 26.9
metadata$Age[which(str_detect(metadata$cells, "CP013_1m"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP013_12m"))] <- 12
metadata$Age[which(str_detect(metadata$cells, "CP013_24m"))] <- 24
metadata$Age[which(str_detect(metadata$cells, "CP018_42m"))] <- 42
metadata$Age[which(str_detect(metadata$cells, "CE021_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CE025_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CE037_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP002_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP003_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP006_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP006_12m"))] <- 12
metadata$Age[which(str_detect(metadata$cells, "CP011_entry"))] <- 2
metadata$Age[which(str_detect(metadata$cells, "CP013_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP016_entry"))] <- 2
metadata$Age[which(str_detect(metadata$cells, "CP017_entry"))] <- 2
metadata$Age[which(str_detect(metadata$cells, "CP018_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP022_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CP042_entry"))] <- 1
metadata$Age[which(str_detect(metadata$cells, "CS005_entry"))] <- 2
metadata$Age[which(str_detect(metadata$cells, "CS015_entry"))] <- 2
metadata$Age[which(str_detect(metadata$cells, "SAAH29_entry"))] <- 2.5
metadata$Age[which(str_detect(metadata$cells, "SATY021_entry"))] <- 1.3

### Create Cohort Column
metadata$Cohort <- NA
metadata$Cohort <- ifelse(str_detect(metadata$cells, "SA"), "EARTH", "TARA")

### Add Viral Load Column
# Initialize Viral_Load column in metadata
metadata$Viral_Load <- NA

# Assign viral loads for the first set
metadata$Viral_Load[which(str_detect(metadata$cells, "CE021_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CE025_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CE037_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP002_entry"))] <- "4284389"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP003_entry"))] <- "656769"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP006_12m"))] <- "73"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP006_entry"))] <- "10000000"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP011_entry"))] <- "36965"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP013_entry"))] <- "3434"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP016_entry"))] <- "1978332"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP017_entry"))] <- "3167384"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP018_entry"))] <- "176970"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP022_entry"))] <- "5075764"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP042_entry"))] <- "6113"
metadata$Viral_Load[which(str_detect(metadata$cells, "CS005_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CS015_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "SAAH29_entry"))] <- "24769"
metadata$Viral_Load[which(str_detect(metadata$cells, "SATY021_entry"))] <- "488"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_CH_009_V1"))] <- NA
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_CH_009_V5"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_CH_009_V7"))] <- "128"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_CH_009_V9"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_TY_026_V2"))] <- NA
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_TY_026_V6"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP003_V12"))] <- "503497"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP003_V24"))] <- "489676"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP006_V24"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP018_V24"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP020_V1"))] <- "414409"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP020_V12"))] <- "6293122.5"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP020_V44"))] <- "158"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_AH_004_V0"))] <- "48700"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_AH_004_V1"))] <- "4279"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_TY_026_V5"))] <- "187"
metadata$Viral_Load[which(str_detect(metadata$cells, "SA_TY_026_V7"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP013_1m"))] <- "3434"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP013_12m"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP013_24m"))] <- "20"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP018_42m"))] <- "20"


### After verifying created metadata is correct, add metadata back to seurat object
merged_seurat@meta.data <- metadata


# Create .RData object to load at any time
save(merged_seurat, file=paste0(out.path,"Seuratv5_CITEseq_dsbnorm_merged_Seurat.RData"))

########################################################################################################################