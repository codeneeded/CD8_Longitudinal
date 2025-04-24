#Load Required Libraries
library(Nebulosa)
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
library(pheatmap)

##set path to load data

setwd('~/Documents/CD8_Longitudinal/Grant_Analysis')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"


############## Merge with seurat object ###################################

load(paste0(load.path,'TARA_ALL_WNN.Rdata'))
load(paste0(load.path,'TARA_HEI_WNN.Rdata'))
load(paste0(load.path,'EARTH_WNN.Rdata'))
#### Select individuals at entry
entry_samples <- grep("_entry$", levels(as.factor(TARA_ALL$orig.ident)), value = TRUE)
TARA_entry <- subset(TARA_ALL, subset = orig.ident %in% entry_samples)

levels(as.factor(TARA_entry$Condition))
TARA_entry <- JoinLayers(TARA_entry,assay = 'RNA')
?JoinLayers
################################################# Markers #############################################################


GeneralInflammation <- c(
  "CRP",      # C-reactive protein
  "SAA1",     # Serum amyloid A1
  "SAA2",     # Serum amyloid A2
  "HP",       # Haptoglobin
  "GCH1",     # Linked to neopterin production
  "IL6",      # Key inflammatory cytokine
  "TNF",      # Tumor necrosis factor
  "CXCL10",   # Interferon-induced chemokine
  "IL1B",     # IL-1β
  "PTX3"      # Pentraxin 3 (inflammatory pattern recognition molecule)
)

TNFSuperfamily <- c(
  # Ligands
  "TNF",        # TNF-alpha
  "TNFSF10",    # TRAIL
  "TNFSF13B",   # BAFF
  "TNFSF13",    # APRIL
  "TNFSF11",    # RANKL
  "TNFSF14",    # LIGHT
  "TNFSF15",    # TL1A
  "FASLG",      # Fas Ligand
  
  # Receptors
  "TNFRSF1A",   # TNFR1
  "TNFRSF1B",   # TNFR2
  "TNFRSF10A",  # TRAIL-R1
  "TNFRSF10B",  # TRAIL-R2
  "TNFRSF13B",  # TACI
  "TNFRSF13C",  # BAFF-R
  "TNFRSF14",   # HVEM
  "TNFRSF25",   # DR3
  "FAS"         # Fas receptor (CD95)
)

Cytokines <- c(
  "IL6", "CXCL8", "IL10", "TNF", "IL18", "IFNA1", "IFNG", "CXCL10", "CXCL11",
  "IL1B", "IL12A", "IL12B", "IL17A", "IL23A", "IL4", "IL5", "IL13",
  "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL20", "CXCL9", "CXCL12",
  "TGFB1"  # Transforming growth factor beta 1
)
ImmuneActivation <- c(
  "TNFRSF1A", "TNFRSF1B", "IL2RA", "CRP", "GCH1", "FGA", "FGB", "FGG",
  "SAA1", "SAA2", "HP", "SERPINA1", "IL1RN", "IDO1", "CD274", "PDCD1LG2",
  "CXCL13", "B2M", "LAG3", "HAVCR2"
)

AdhesionMolecules <- c(
  "VCAM1", "ICAM1", "SELL", "SELE",
  "ITGAL", "ITGAM", "ITGB2", "PECAM1", "CDH5", "SPN", "CD44", "JAM3"
)

GutTranslocation <- c(
  "LBP", "CD14", "CD163", "CCL2",
  "REG3A", "DEFA5", "DEFB1", "TJP1", "OCLN", "CLDN2", "CLDN3", "MUC2",
  "S100A8", "S100A9", "TLR2", "TLR4", "NOD2", "MYD88"
)

SenescenceMarkers <- c(
  "CDKN2A",   # p16^INK4a
  "CDKN1A",   # p21^CIP1
  "TP53",     # p53
  "GLB1",     # β-galactosidase (SA-β-gal gene)
  "SERPINE1", # PAI-1 (SASP component)
  "IL6",      # SASP cytokine
  "MMP3",     # SASP protease
  "MMP9",     # SASP protease
  "CXCL8",    # IL-8, SASP chemokine
  "CCL2"      # SASP chemokine
)

ExhaustionMarkers <- c(
  "PDCD1",    # PD-1
  "CTLA4",    # CTLA-4
  "LAG3",     # LAG-3
  "HAVCR2",   # TIM-3
  "TIGIT",    # TIGIT
  "ENTPD1",   # CD39
  "TOX",      # Exhaustion-driving transcription factor
  "EOMES",    # Transcriptional regulator (context-dependent, up in exhausted)
  "BATF"      # Cooperates with TOX in exhaustion programs
)


################################ Global Comparison across all PBMC #######################################
gene_sets <- list(
  Cytokines = Cytokines,
  ImmuneActivation = ImmuneActivation,
  AdhesionMolecules = AdhesionMolecules,
  GutTranslocation = GutTranslocation,
  SenescenceMarkers = SenescenceMarkers,
  ExhaustionMarkers = ExhaustionMarkers,
  GeneralInflammation = GeneralInflammation,
  TNFSuperfamily = TNFSuperfamily
)

avg_expr_list <- list()

for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  
  avg_exp <- tryCatch({
    AverageExpression(TARA_entry,
                      group.by = "Condition",
                      features = genes,
                      assays = "RNA",
                      return.seurat = FALSE)
  }, error = function(e) NULL)
  
  avg_expr_list[[set_name]] <- avg_exp
}

### DGE for HEI vs HEU and HEU vs HUU #########################
degs_stars <- list()

get_significance_labels <- function(seurat_obj, genes, group1, group2) {
  # Safety check: skip if gene list is empty
  if (length(genes) == 0) return(rep("", 0))
  
  markers <- tryCatch({
    FindMarkers(
      seurat_obj,
      ident.1 = group1,
      ident.2 = group2,
      group.by = "Condition",
      features = genes,
      assay = "RNA",
      logfc.threshold = 0,
      min.pct = 0
    )
  }, error = function(e) {
    message(paste("Skipping", group1, "vs", group2, "- error:", e$message))
    return(data.frame())
  })
  
  stars <- rep("", length(genes))
  names(stars) <- genes
  
  if (nrow(markers) == 0) return(stars)
  
  sig_genes <- intersect(rownames(markers), genes)
  pvals <- markers[sig_genes, "p_val_adj"]
  
  stars[sig_genes[pvals < 0.05]] <- "*"
  stars[sig_genes[pvals < 0.01]] <- "**"
  stars[sig_genes[pvals < 0.001]] <- "***"
  return(stars)
}

for (set_name in names(gene_sets)) {
  genes <- intersect(gene_sets[[set_name]], rownames(TARA_entry[["RNA"]]))
  
  stars_list <- list()
  stars_list$HEI_vs_HEU <- get_significance_labels(TARA_entry, genes, "HEI", "HEU")
  stars_list$HEU_vs_HUU <- get_significance_labels(TARA_entry, genes, "HEU", "HUU")
  
  degs_stars[[set_name]] <- stars_list
}
###################

# Custom color palette (no red)
my_colors <- colorRampPalette(c("navy", "white", "goldenrod"))(100)


# Mapping from original Seurat condition labels to display names
condition_labels <- c("HEI" = "pHIV", "HEU" = "pHEU", "HUU" = "pHUU")

# Create output folder
dir.create("PBMC_Heatmaps", showWarnings = FALSE)

# Store cleaned, scaled, renamed matrices
clean_matrices <- list()


for (set_name in names(avg_expr_list)) {
  avg_exp <- avg_expr_list[[set_name]]
  if (is.null(avg_exp) || !"RNA" %in% names(avg_exp)) {
    message("Skipping ", set_name, ": no valid RNA entry.")
    next
  }
  
  mat <- avg_exp$RNA
  
  # Force coercion if not a matrix
  if (!is.matrix(mat)) {
    mat <- tryCatch(as.matrix(mat), error = function(e) {
      message("Skipping ", set_name, ": failed matrix coercion.")
      return(NULL)
    })
    if (is.null(mat)) next
  }
  
  # Ensure sample columns exist and rename
  present_conditions <- intersect(colnames(mat), names(condition_labels))
  if (length(present_conditions) < 2) {
    message("Skipping ", set_name, ": not enough condition columns.")
    next
  }
  
  colnames(mat) <- condition_labels[colnames(mat)]
  
  # Scale rows and clean
  scaled <- t(scale(t(mat)))
  scaled <- scaled[apply(scaled, 1, function(x) all(is.finite(x))), , drop = FALSE]
  
  if (nrow(scaled) == 0) {
    message("Skipping ", set_name, ": no valid rows after cleaning.")
    next
  }
  
  clean_matrices[[set_name]] <- scaled
}


#### 2 add stars to rownames
# Store row label vectors with stars
label_rows <- list()

for (set_name in names(clean_matrices)) {
  mat <- clean_matrices[[set_name]]
  rows <- rownames(mat)
  
  label_list <- list()
  
  if (!is.null(degs_stars[[set_name]]$HEI_vs_HEU)) {
    stars <- degs_stars[[set_name]]$HEI_vs_HEU
    label_list$HEI_vs_HEU <- paste0(rows, " ", stars[rows])
  }
  
  if (!is.null(degs_stars[[set_name]]$HEU_vs_HUU)) {
    stars <- degs_stars[[set_name]]$HEU_vs_HUU
    label_list$HEU_vs_HUU <- paste0(rows, " ", stars[rows])
  }
  
  label_rows[[set_name]] <- label_list
}

######### Plot

for (set_name in names(clean_matrices)) {
  mat <- clean_matrices[[set_name]]
  labels <- label_rows[[set_name]]
  
  # HEI vs HEU
  if (all(c("pHIV", "pHEU") %in% colnames(mat))) {
    png(paste0("PBMC_Heatmaps/", set_name, "_pHIV_vs_pHEU.png"),
        width = 2000, height = 3000, res = 300)
    pheatmap(mat[, c("pHIV", "pHEU")],
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = my_colors,
             labels_row = labels$HEI_vs_HEU,
             main = paste("pHIV vs pHEU -", set_name))
    dev.off()
  }
  
  # HEU vs HUU
  if (all(c("pHEU", "pHUU") %in% colnames(mat))) {
    png(paste0("PBMC_Heatmaps/", set_name, "_pHEU_vs_pHUU.png"),
        width = 2000, height = 3000, res = 300)
    pheatmap(mat[, c("pHEU", "pHUU")],
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = my_colors,
             labels_row = labels$HEU_vs_HUU,
             main = paste("pHEU vs pHUU -", set_name))
    dev.off()
  }
  
  # Full comparison
  if (all(c("pHIV", "pHEU", "pHUU") %in% colnames(mat))) {
    png(paste0("PBMC_Heatmaps/", set_name, "_pHIV_pHEU_pHUU.png"),
        width = 2500, height = 3000, res = 300)
    pheatmap(mat,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = my_colors,
             main = paste("Bulk PBMC Expression -", set_name))
    dev.off()
  }
}














################
# Custom color palette: blue → white → gold (avoiding red)
my_colors <- colorRampPalette(c("navy", "white", "goldenrod"))(100)

# Output directory
dir.create("PBMC_Heatmaps", showWarnings = FALSE)

for (set_name in names(gene_sets)) {
  genes <- gene_sets[set_name]
  
  avg_exp <- AverageExpression(TARA_entry,
                               group.by = "Condition",
                               features = genes,
                               assays = 'RNA',
                               return.seurat = FALSE)
  
  scaled_matrix <- t(scale(t(avg_exp$RNA)))
  clean_scaled_matrix <- scaled_matrix[apply(scaled_matrix, 1, function(x) all(is.finite(x))), ]
  
  if (nrow(clean_scaled_matrix) == 0) next  # Skip if no usable data
  
  # Determine how many comparisons are possible
  has_HEI_HEU <- all(c("HEI", "HEU") %in% colnames(clean_scaled_matrix))
  has_HEU_HUU <- all(c("HEU", "HUU") %in% colnames(clean_scaled_matrix))
  
  # Full comparison: HEI, HEU, HUU
  if (has_HEI_HEU && has_HEU_HUU) {
    png(paste0("PBMC_Heatmaps/", set_name, "_HEI_HEU_HUU.png"),
        width = 2500, height = 3000, res = 300)
    pheatmap(clean_scaled_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = my_colors,
             main = paste("Bulk PBMC Expression -", set_name))
    dev.off()
  }
  
  # Pairwise: HEI vs HEU
  if (has_HEI_HEU) {
    png(paste0("PBMC_Heatmaps/", set_name, "_HEI_vs_HEU.png"),
        width = 2000, height = 3000, res = 300)
    pheatmap(clean_scaled_matrix[, c("HEI", "HEU")],
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = my_colors,
             main = paste("HEI vs HEU -", set_name))
    dev.off()
  }
  
  # Pairwise: HEU vs HUU
  if (has_HEU_HUU) {
    png(paste0("PBMC_Heatmaps/", set_name, "_HEU_vs_HUU.png"),
        width = 2000, height = 3000, res = 300)
    pheatmap(clean_scaled_matrix[, c("HEU", "HUU")],
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = my_colors,
             main = paste("HEU vs HUU -", set_name))
    dev.off()
  }
}
