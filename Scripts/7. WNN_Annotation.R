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

### Read Files in

setwd('~/Documents/CD8_Longitudinal/Annotation')

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"



### TARA
load(paste0(load.path,'TARA_TCR_Combined.RData'))
load(paste0(load.path,'EARTH_TCR_Combined.RData'))


######## Rejoin layers ####
TARA_ALL[["RNA"]] <- JoinLayers(TARA_ALL[["RNA"]])
TARA_ALL[["ADT"]] <- JoinLayers(TARA_ALL[["ADT"]])

rm(TARA_ALL_TRB_0, TARA_ALL_TRB_1, TARA_ALL_TRB_2)
gc()


EARTH[["RNA"]] <- JoinLayers(EARTH[["RNA"]])
EARTH[["ADT"]] <- JoinLayers(EARTH[["ADT"]])

rm(EARTH_TRB_0, EARTH_TRB_1, EARTH_TRB_2)
gc()
#################################### RNA and Protein Features of interest ##############################################library(dplyr)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)

# Define IFN gene list
IFN_genes_all <- c(
  "CCR5", "MX1", "IFIT1", "IFIT2", "IFIT3", "ISG15", "OAS1", 
  "STAT1", "STAT2", "IRF7", "IFI6", "IFI44", "IFI44L", 
  "ISG20", "MX2", "OAS2", "OAS3", "IRF9", "RSAD2", 
  "BST2", "USP18", "TRIM22", "SIGLEC1"
)

# Define comparisons and paths
comparisons <- list(
  HEIvsHEU = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/HEIvsHEU",
  HEUvsHUU = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/HEUvsHUU",
  PostARTvsPreART = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_ALL/PostARTvsPreART",
  HighvsLowVL = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Expression/TARA_Entry/HighvsLowVL"
)

base_output <- "/home/akshay-iyer/Documents/CD8_Longitudinal/IFNa_Analysis"

# Loop through each comparison
for (comp_name in names(comparisons)) {
  
  de_dir <- comparisons[[comp_name]]
  output_dir <- file.path(base_output, comp_name)
  volcano_dir <- file.path(output_dir, "Volcano_Plots")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
  
  de_files <- list.files(de_dir, pattern = "\\.csv$", full.names = TRUE)
  summary_list <- list()
  
  # Loop through each cluster DE file
  for (file in de_files) {
    df <- read.csv(file, row.names = 1)
    df <- df %>%
      mutate(
        gene = rownames(.),
        pval = p_val,
        padj = p_val_adj,
        logfc = avg_log2FC
      ) %>%
      filter(gene %in% IFN_genes_all & !is.na(padj) & padj < 0.05)
    
    up <- sum(df$logfc > 0, na.rm = TRUE)
    down <- sum(df$logfc < 0, na.rm = TRUE)
    total <- nrow(df)
    
    cluster <- str_remove(basename(file), "\\.csv$")
    cell_type <- str_split(cluster, "__", simplify = TRUE)[, 2]
    
    summary_list[[cluster]] <- data.frame(
      comparison = comp_name,
      cluster_id = cluster,
      cell_type = cell_type,
      total_ifn_genes = total,
      upregulated = up,
      downregulated = down
    )
    
    ##### Volcano Plot for this cluster #####
    full_df <- read.csv(file, row.names = 1)
    full_df <- full_df %>%
      mutate(
        gene = rownames(.),
        logfc = avg_log2FC,
        pval = p_val,
        padj = p_val_adj
      )
    
    keyvals <- ifelse(
      abs(full_df$logfc) > 1.5 & full_df$padj < 0.01, "#CD0BBC",
      ifelse(full_df$padj < 0.01, "#28E2E5", "gray30")
    )
    keyvals[full_df$gene %in% IFN_genes_all] <- "#FF9900"
    
    names(keyvals)[keyvals == "gray30"] <- "NS"
    names(keyvals)[keyvals == "#28E2E5"] <- "adj(p-value) < 0.01"
    names(keyvals)[keyvals == "#CD0BBC"] <- "FC > 1.5 & p < 0.01"
    names(keyvals)[keyvals == "#FF9900"] <- "IFN Gene"
    
    vp <- EnhancedVolcano(full_df,
                          lab = full_df$gene,
                          x = "logfc",
                          y = "pval",
                          pCutoffCol = "padj",
                          pCutoff = 0.01,
                          FCcutoff = 1.5,
                          xlab = bquote(~Log[2]~ 'fold change'),
                          ylab = bquote(-Log[10]~italic(P)),
                          pointSize = 2.0,
                          labSize = 2.5,
                          colCustom = keyvals,
                          colAlpha = 0.8,
                          labCol = "black",
                          labFace = "bold",
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = "gray40",
                          legendPosition = "right",
                          legendLabSize = 10,
                          legendIconSize = 3.5,
                          title = paste0(cluster),
                          subtitle = paste0("Volcano Plot - ", comp_name, " (IFNα genes highlighted)")
    )
    
    ggsave(
      filename = file.path(volcano_dir, paste0(cluster, "_volcano.png")),
      plot = vp + guides(color = guide_legend(reverse = TRUE)),
      width = 10, height = 7, dpi = 400
    )
  }
  
  # Save cluster-wise IFN summary
  ifn_summary <- bind_rows(summary_list)
  write.csv(ifn_summary, file.path(output_dir, "IFN_gene_summary_by_cluster.csv"), row.names = FALSE)
  
  # Plot total IFN genes per cluster
  ifn_plot <- ggplot(ifn_summary, aes(x = reorder(cluster_id, -total_ifn_genes), y = total_ifn_genes)) +
    geom_col(fill = "steelblue") +
    theme_minimal(base_size = 14) +
    xlab("Cluster") + ylab("Significant IFN Genes (adj.p < 0.05)") +
    ggtitle(paste("IFN Gene Expression by Cluster -", comp_name)) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    filename = file.path(output_dir, "IFN_gene_barplot_by_cluster.png"),
    plot = ifn_plot,
    width = 13,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}


rna.features <-  c('CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
                   'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
                   'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
                   'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
                   'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
                   'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
                   'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
                   'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67')

prots <- rownames(TARA_ALL@assays$ADT)

Idents(TARA_ALL) <- 'snn.louvianmlr_1'
Idents(EARTH) <- 'snn.louvianmlr_1'

#Remove clusters with <50 cells
# Identify clusters with fewer than 20 cells
cluster_sizes <- table(Idents(TARA_ALL))
small_clusters <- names(cluster_sizes[cluster_sizes < 100])
large_clusters <- names(cluster_sizes[cluster_sizes > 20])
TARA_ALL <- subset(TARA_ALL, idents = large_clusters)

cluster_sizes <- table(Idents(EARTH))
small_clusters <- names(cluster_sizes[cluster_sizes < 100])
large_clusters <- names(cluster_sizes[cluster_sizes > 20])
EARTH <- subset(EARTH, idents = large_clusters)

############# Plots ###############################
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Cluster_Plot')
saveRDS(TARA_ALL, file = file.path(load.path, "TARA_ALL_post_annotation.rds"))

p1 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

ggsave(
  filename = "TARA_WNN_lmlr1.png",
  plot=p1,
  width = 10,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300
)

p2 <- DimPlot_scCustom(
  EARTH,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

ggsave(
  filename = "EARTH_WNN_lmlr1.png",
  plot=p2,
  width = 10,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300 
)

#### Average Expression ##############


avg_TARA <- AverageExpression(TARA_ALL, return.seurat = FALSE)
write.csv(avg_TARA$RNA, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/avg_rna_TARA.csv")
write.csv(avg_TARA$ADT, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/avg_adt_TARA.csv")

avg_EARTH <- AverageExpression(EARTH, return.seurat = FALSE)
write.csv(avg_EARTH$RNA, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/EARTH/avg_rna_EARTH.csv")
write.csv(avg_EARTH$ADT, "/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/EARTH/avg_adt_EARTH.csv")



########################################## Feature Plots and VLN Plots ###############################################
### TARA
### ADT

DefaultAssay(TARA_ALL) <-'ADT'

# VLN Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Violin_Plot/ADT')

for (i in prots) {
  
  vln.pl <-VlnPlot_scCustom(TARA_ALL,features = i)
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 14, vln.pl)
}

# Feature Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Feature_Plot/ADT')

for (i in prots) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(TARA_ALL, reduction = 'wnn.umap', features = i, 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_Magma.png'),dpi=500, width = 8, fea.pl)
}

### RNA

DefaultAssay(TARA_ALL) <-'RNA'

# VLN Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Violin_Plot/RNA')

for (i in rna.features) {
  
  vln.pl <-VlnPlot_scCustom(TARA_ALL,features = i)
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 14, vln.pl)
}

# Feature Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Feature_Plot/RNA')

for (i in rna.features) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(TARA_ALL, reduction = 'wnn.umap', features = i, 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_Magma.png'),dpi=500, width = 8, fea.pl)
}

### Earth
### ADT

DefaultAssay(EARTH) <-'ADT'

# VLN Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/EARTH/Violin_Plot/ADT')

for (i in prots) {
  
  vln.pl <-VlnPlot_scCustom(EARTH,features = i)
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 14, vln.pl)
}

# Feature Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/EARTH/Feature_Plot/ADT')

for (i in prots) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(EARTH, reduction = 'wnn.umap', features = i, 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_Magma.png'),dpi=500, width = 8, fea.pl)
}

### RNA

DefaultAssay(EARTH) <-'RNA'

# VLN Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/EARTH/Violin_Plot/RNA')

for (i in rna.features) {
  
  vln.pl <-VlnPlot_scCustom(EARTH,features = i)
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 14, vln.pl)
}

# Feature Plots
setwd('~/Documents/CD8_Longitudinal/Annotation/EARTH/Feature_Plot/RNA')

for (i in rna.features) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(EARTH, reduction = 'wnn.umap', features = i, 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_Magma.png'),dpi=500, width = 8, fea.pl)
}

############################################################ Annotated Plots ####################################################
setwd('~/Documents/CD8_Longitudinal/Annotation/TARA_ALL')

TARA_ALL$Manual_Annotation <- NULL
# Set cluster identity to Louvain clustering
Idents(TARA_ALL) <- 'snn.louvianmlr_1'

# Load and process annotation file
TARA_annotation_df <- read_excel("TARA_ALL_Annotation.xlsx") %>%
  rename(cluster = `Cluster Number`, annotation = `Cell Type`) %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(cluster %in% as.character(0:35)) %>%
  mutate(label = paste0(cluster, ": ", annotation))

# Get cluster order from Seurat object
cluster_levels <- levels(Idents(TARA_ALL))

# Create named vector for mapping cluster -> label
TARA_cluster_map <- setNames(TARA_annotation_df$label, TARA_annotation_df$cluster)

# Map annotations to cells based on cluster identity
cluster_ids <- as.character(Idents(TARA_ALL))
TARA_ALL$Manual_Annotation <- unname(TARA_cluster_map[cluster_ids])

# Ensure factor levels follow cluster order
TARA_ALL$Manual_Annotation <- factor(TARA_ALL$Manual_Annotation, levels = TARA_cluster_map[cluster_levels])

# Set identity class to Manual_Annotation
Idents(TARA_ALL) <- 'Manual_Annotation'

Idents(TARA_ALL)

p1 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE,
  label.box = TRUE,
  label.size = 3.5,
  pt.size = 1
)

ggsave(
  filename = "TARA_Annotated_1.png",
  plot=p1,
  width = 13,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300
)

p2 <- DimPlot_scCustom(
  TARA_ALL,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE,
  label.box = TRUE,
  label.size = 3.5,
  pt.size = 0.1
)

ggsave(
  filename = "TARA_Annotated_0.1.png",
  plot=p2,
  width = 13,  # Adjust width as needed
  height = 8,  # Adjust height as needed
  dpi = 300
)

# Create frequency table
bar_data <- as.data.frame(table(TARA_ALL$Manual_Annotation))
colnames(bar_data) <- c("Cluster", "Cell_Count")

# Get colors matching cluster identities
cluster_colors <- DiscretePalette_scCustomize(
  num_colors = length(levels(TARA_ALL$Manual_Annotation)),
  palette = "polychrome"
)

names(cluster_colors) <- levels(TARA_ALL$Manual_Annotation)

# Bar plot
p <- ggplot(bar_data, aes(x = Cluster, y = Cell_Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  xlab("Manual Cluster Annotation") +
  ylab("Number of Cells") +
  ggtitle("TARA: Cell Count per Annotated Cluster")

ggsave("TARA_Cluster_Barplot.png", plot = p, width = 10, height = 6, bg='white')


saveRDS(TARA_ALL, file = file.path(load.path, "TARA_ALL_post_annotation.rds"))


### Differential Expression ###





### EARTH


