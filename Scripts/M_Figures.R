################# Libraries ##################################################
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
library(qs2)

############################## Read files ######################################
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")

TARA_ALL <- qs_read(
  file = file.path(saved_dir, "TARA_ALL_sorted.qs2")
)

rds_in     <- file.path(saved_dir, "tara_cdnk_annotated.rds")
tara_cdnk <- readRDS(rds_in)

######################################################### Counts Flowchart #########################################
exp_levels <- c(
  "Medium (5 < X <= 20)",
  "Large (20 < X <= 100)",
  "Hyperexpanded (100 < X <= 500)"
)

tara_cdnk$Expanding <- tara_cdnk$cloneSize %in% exp_levels

# total cells
ncol(TARA_ALL)
ncol(tara_cdnk)

table(TARA_ALL$Timepoint_Group)
table(tara_cdnk$Expanding) 
#################################################### Figure 1 ######################################################
DefaultAssay(TARA_ALL) <- 'RNA'

############# 1A UMAP
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 1/01_UMAP"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep_groups <- c("PreART_Entry", "HEU", "HUU")

TARA_fig1 <- subset(
  TARA_ALL,
  subset = Timepoint_Group %in% keep_groups
)

TARA_fig1$Timepoint_Group <- factor(
  TARA_fig1$Timepoint_Group,
  levels = c("HUU", "HEU", "PreART_Entry")
)
p_celltypes <- DimPlot2(TARA_fig1,
  reduction = "wnn.umap",
  group.by = "Manual_Annotation",
  cols='default',
  label = F,
  repel = TRUE,
  raster = FALSE,
  , theme = NoLegend()
)

ggsave(
  filename = file.path(out_dir, "Fig1A_DimPlot2_CellTypes.png"),
  plot = p_celltypes,
  width = 6,
  height = 6,
  bg='white'
)
TARA_fig1$predicted.celltype.l2
DimPlot2(TARA_ALL,
                        reduction = "wnn.umap",
                        group.by = "Manual_Annotation",
                        cols='default',
                        label = T,
                        repel = TRUE,
                        raster = FALSE,
)

############### 1B - Cluster Distribution #################

table(TARA_ALL$Manual_Annotation)
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 1/02_Cluster_Proportions"

pbmc_sig <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(
    Timepoint_Group %in% c("HUU","HEU", "PreART_Entry"),
    Manual_Annotation %in% c("1: Memory CD8 T cell", "2: Naïve CD4 T cell_1", "8: CTL-like")
  ) %>%
  transmute(
    sample_id = orig.ident,
    cluster   = Manual_Annotation,
    condition = Timepoint_Group
  )

pbmc_sig$condition <- factor(pbmc_sig$condition, levels = c("HUU", "HEU", "PreART_Entry"))
pbmc_sig$cluster   <- factor(pbmc_sig$cluster, levels = c("1: Memory CD8 T cell", "2: Naïve CD4 T cell_1", "8: CTL-like"))

p_sig <- ClusterDistrPlot(
  origin    = pbmc_sig$sample_id,
  cluster   = pbmc_sig$cluster,
  condition = pbmc_sig$condition,
  stat.method = "wilcox.test",
  plot = TRUE
)
ggsave(
  filename = file.path(out_dir, "Fig1B_ClusterProportions_PreART_vs_HEU_Clusters1_2_8.png"),
  plot = p_sig,
  width = 7,
  height = 4,
  dpi = 600,
  bg = "white"
)
############### 1E

library(EnhancedVolcano)



in_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/HEIvsHEU_PreART"
in_csv <- file.path(in_dir, "1_Memory_CD8_T_cell_HEI_vs_HEU.csv")

df <- read_csv(in_csv, show_col_types = FALSE)

df2 <- df
colnames(df2)[1] <- "gene"

df2$p_val_adj <- ifelse(
  is.na(df2$p_val_adj) | df2$p_val_adj <= 0,
  1e-300,
  df2$p_val_adj
)

head(df2)

label_genes <- c(
  "APOO","GBP5","AC105052.4","RBFOX2","ANK3","GBP1","STAT1","AC079793.1","CCR7",
  "SGK1","AL138963.4","GSTM3","SPINK2","HDLBP","MTRNR2L12","AC004448.2","FBLN2",
  "IFIT1","ISG15","MX1","GZMK","NKG7","FGFBP2","CCR7","IL7R","LAG3","TCF7"
)

p <- EnhancedVolcano(
  df2,
  lab = df2$gene,
  x = "avg_log2FC",
  y = "p_val_adj",
  title = "Memory CD8 T cell: HEI vs HEU (PreART)",
  subtitle = NULL,
  pCutoff = 0.05,
  FCcutoff = 0.25,
  selectLab = intersect(label_genes, df2$gene),
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  boxedLabels = TRUE,
  labSize = 4,
  pointSize = 2.2,
  legendPosition = "none"
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black")
  )

p

getwd()
# Save PNG
manuscript_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 1"
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
volcano_dir <- file.path(manuscript_dir, "04_Volcano")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(volcano_dir, "Fig1_MemoryCD8_HEI_vs_HEU_Volcano.png")

ggsave(out_png,p, width = 7.5, height = 6.5, dpi = 600, bg = "white")
