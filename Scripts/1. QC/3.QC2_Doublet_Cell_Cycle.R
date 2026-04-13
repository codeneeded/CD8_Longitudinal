###Citeseq Pipeline
#Cite seurat and ds packages
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
library(DoubletFinder)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
#library(Nebulosa)
#library(tricycle)
#library(org.Hs.eg.db)
library(Azimuth)
library(scDblFinder)

# Data Input, Creation of Merged Seurat Object
## https://samuel-marsh.github.io/scCustomize/articles/Sequencing_QC_Plots.html -> Check this out, very cool
##set path to load data

setwd("~/Documents/CD8_Longitudinal/QC")  # Set working directory
load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/" 
# Load Data

load(paste0(load.path,'Seuratv5_filtered_seurat.RData'))
load(paste0(load.path,'Cell_Cycle_Genes.RData'))

rownames(filtered_seurat@assays$ADT@data)

###################################Complete Initialisation Steps ######################################################
DefaultAssay(filtered_seurat) <-'RNA'
### Azimuth Analysis
seurat_phase <- RunAzimuth(filtered_seurat, reference = "pbmcref")

### Basic Transformations
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst")

#### Use Sketch to Create a subset of 50k cells

seurat_phase <- SketchData(
  object = seurat_phase,
  ncells = 150000,
  method = "LeverageScore",
  sketched.assay = "sketch.rna"
)

# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(seurat_phase) <- "sketch.rna"
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst")
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "pca")
seurat_phase <- FindClusters(seurat_phase, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated",return.model=T)
DimPlot(seurat_phase, group.by = 'predicted.celltype.l2', label = T,repel = T)

# Project data back to original datast
seurat_phase <- ProjectData(
  object = seurat_phase,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch.rna",
  sketched.reduction = "pca",
  umap.model = "umap.unintegrated",
  dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters")
)
DefaultAssay(seurat_phase) <- "RNA"
DimPlot(seurat_phase, label = T, label.size = 3, reduction = "full.umap.unintegrated", group.by = "predicted.celltype.l2", alpha = 0.2,repel = T) + NoLegend()

#save(seurat_phase, file = paste0(load.path,'SeuratV5_SeuratPhase.RData'))
#load(paste0(load.path,'SeuratV5_SeuratPhase.RData'))

############################################## EXPLORE CELL CYCLE ########################################################
## https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

setwd("~/Documents/CD8_Longitudinal/QC/Cell_Cycle")  # Set working directory

s.genes <- s_genes
g2m.genes <- g2m_genes

#Evaluate effects of Cell Cycle on Data
FeaturePlot(object = seurat_phase, features = 'NCAM1')

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes,
                                 set.ident = TRUE)


#seurat_phase <- RunPCA(seurat_phase, features = c(s.genes, g2m.genes))
DimPlot(seurat_phase)
# Visualize the distribution of cell cycle markers across
RidgePlot(seurat_phase, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# View cell cycle scores and phases assigned to cells 
png(file="Cell_Cycle_Phase_PerSample.png", width = 2100, height = 1200)
VlnPlot(seurat_phase, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", pt.size = 0.1)
dev.off()

# Plot the PCA colored by cell cycle phase
png(file="PCA_by_Cell_Cycle_Split.png", width = 2100, height = 1200)
DimPlot(seurat_phase,
        reduction = "pca.full",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

png(file="UMAP_by_Cell_Cycle_Split.png", width = 2100, height = 1200)
DimPlot(seurat_phase,
        reduction = "full.umap.unintegrated",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

#Check effect of Condition
png(file="UMAP_by_Cell_Cycle_Split_by_phasae_groupby_Condition.png", width = 2100, height = 1200)
DimPlot(seurat_phase,
        reduction = "full.umap.unintegrated",
        group.by= "Condition",
        split.by = "Phase")
dev.off()

png(file="UMAP_by_ident_and_condition.png", width = 2100, height = 1200)
DimPlot(seurat_phase, label = T, label.size = 3, reduction = "full.umap.unintegrated", split.by = "Condition", group.by='orig.ident', alpha = 0.2,repel = T) + NoLegend()
dev.off()

#####################################################################################################################
##As an alternative, we suggest regressing out the difference between the G2M and S phase scores. 
##This means that signals separating non-cycling cells and cycling cells will be maintained, 
##but differences in cell cycle phase among proliferating cells (which are often uninteresting), 
##will be regressed out of the data

############################################Doublet Finder###########################################################
### 
#load(paste0(load.path,'SeuratV5_SeuratPhase.RData'))
# Split Seurat Object For Doublet Filtration

split_seurat <- SplitObject(seurat_phase, split.by = "orig.ident")

samples <-  levels(as.factor(seurat_phase$orig.ident))

setwd("~/Documents/CD8_Longitudinal/QC/Doublets")

for (i in samples) {
  sce <- scDblFinder(GetAssayData(split_seurat[[i]], slot="counts"))
  split_seurat[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seurat[[i]]$scDblFinder.class <- sce$scDblFinder.class
  plot <- DimPlot(split_seurat[[i]],group.by = 'scDblFinder.class') + ggtitle(paste0(i,' Doublets'))
  ggsave(paste0(i,'_doublets.png'),device = 'png', width = 6, 
         height = 6, dpi = 600,plot)
}

#save(split_seurat, file = paste0(load.path,'SeuratV5_SplitSeurat_Preintegration_WithDoublets.RData')) -> DID NOT USE THIS TIME
### Remove Doublets

for (i in samples) {
  split_seurat[[i]]<- subset(split_seurat[[i]],subset= `scDblFinder.class`== 'singlet')
}

save(split_seurat, file = paste0(load.path,'SeuratV5_SplitSeurat_Preintegration.RData'))

