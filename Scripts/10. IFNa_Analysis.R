# Install (if not already) and load required packages
library(msigdbr)
library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(ggplot2)
library(Seurat)
library(ggpubr)
library(tidyr)

# 1. Load MSigDB gene sets for human
msig_human <- msigdbr(species = "Homo sapiens")

# 2. Extract gene sets of interest
hallmark_ifna <- msig_human %>%
  filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene_symbol)

go_ifna <- msig_human %>%
  filter(gs_name == "GOBP_RESPONSE_TO_INTERFERON_ALPHA") %>%
  pull(gene_symbol)

reactome_ifna <- msig_human %>%
  filter(gs_name == "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING") %>%
  pull(gene_symbol)

# 3. Manual ISG list (from literature)
manual_isg_combined <- c(
  "BST2", "CCR5", "CD274", "CMPK2", "CXCL10", "CXCL11", "DDX58", "DHX58", "EPSTI1", 
  "GBP1", "GBP2", "GBP4", "HERC5", "IFI27", "IFI35", "IFI44", "IFI44L", "IFI6", 
  "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IRF7", "IRF9", "ISG15", "ISG20", 
  "LY6E", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "PLSCR1", "RSAD2", "SIGLEC1", 
  "SP100", "STAT1", "STAT2", "TAP1", "TRIM22", "TRIM5", "USP18", "XAF1"
)

# 4. Combine all lists and remove duplicates
IFN_genes_all <- unique(c(hallmark_ifna, go_ifna, reactome_ifna, manual_isg_combined))

# Optional: Sort alphabetically
IFN_genes_all <- sort(IFN_genes_all)

# 5. Output length and preview
length(IFN_genes_all)
head(IFN_genes_all, 10)
##############

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

################# IFN a Pathway Comparisons #########################################


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
    
    clean_title <- gsub("__", ": ", cluster)              # turn '0__CD4_T_cell' into '0: CD4_T_cell'
    clean_title <- gsub("_", " ", clean_title)            # turn 'CD4_T_cell' into 'CD4 T cell'
    
    
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
                          title = clean_title,
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

######################################################## IDNa Module Scoring ##################################################################

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"
TARA_ALL <- readRDS(file = file.path(load.path, "TARA_ALL_post_annotation.rds"))
DefaultAssay(TARA_ALL) <-'RNA'

# Exclude CCR5 if you're using it for correlation
IFN_genes_for_score <- setdiff(IFN_genes_all, "CCR5")

TARA_ALL <- AddModuleScore(
  object = TARA_ALL,
  features = list(IFN_genes_for_score),
  name = "IFN_score"
)

# For ADT (dsb-normalized CCR5 protein expression)
DefaultAssay(TARA_ALL) <-'RNA'
TARA_ALL$RNA_CCR5 <- FetchData(TARA_ALL, vars = "CCR5", layer = "data", assay = "ADT")[, 1]
DefaultAssay(TARA_ALL) <-'ADT'
TARA_ALL$Protein_CCR5 <- FetchData(TARA_ALL, vars = "CCR5", layer = "data", assay = "ADT")[, 1]

### Functions for plotting 



plot_module_score_unpaired <- function(seurat_obj, group_var, score_var = "IFN_score1", output_path, plot_title) {
  meta <- seurat_obj@meta.data
  gg <- ggplot(meta, aes_string(x = group_var, y = score_var, fill = group_var)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.3) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(~Manual_Annotation, scales = "free_y") +
    theme_minimal(base_size = 14) +
    xlab("") + ylab("IFN Module Score") +
    ggtitle(plot_title) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.margin = margin(10, 10, 15, 10)
    )
  
  ggsave(output_path, plot = gg, width = 20, height = 18, dpi = 300, bg='white')
}


# Ensure PID is extracted
TARA_ALL$PID <- sub("_.*", "", TARA_ALL$orig.ident)

# Output directory
module_score_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/IFNa_Analysis/Module_scores"
dir.create(module_score_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------
## 1. HEI vs HEU
## ------------------------
subset1 <- subset(TARA_ALL, (Condition == "HEI" & Age <= 2) | Condition == "HEU")

plot_module_score_unpaired(
  seurat_obj = subset1,
  group_var = "Condition",
  output_path = file.path(module_score_dir, "HEI_vs_HEU_IFNscore.png"),
  plot_title = "IFN Pathway Activity: HEI (Entry, <2 Months) vs HEU"
)

## ------------------------
## 2. HEU vs HUU
## ------------------------
subset2 <- subset(TARA_ALL, Condition %in% c("HEU", "HUU"))

plot_module_score_unpaired(
  seurat_obj = subset2,
  group_var = "Condition",
  output_path = file.path(module_score_dir, "HEU_vs_HUU_IFNscore.png"),
  plot_title = "IFN Pathway Activity: HEU vs HUU"
)

## ------------------------
## 3. High vs Low Viral Load (within HEI only)
## ------------------------
subset3 <- subset(TARA_ALL, Condition == "HEI"& Age <= 2)
subset3$VL_Group <- ifelse(as.numeric(as.character(subset3$Viral_Load)) >= 100000, "High", "Low") 
plot_module_score_unpaired(
  seurat_obj = subset3,
  group_var = "VL_Group",
  output_path = file.path(module_score_dir, "HighVL_vs_LowVL_IFNscore.png"),
  plot_title = "IFN Pathway Activity: High vs Low Viral Load (HEI, Entry)"
)

## ------------------------
## 4. Pre vs Post ART (within HEI)
## ------------------------
subset4 <- subset(TARA_ALL, Condition == "HEI")
subset4$ART_Group <- ifelse(
  subset4$Age <= 2, "PreART_Entry",
  ifelse(subset4$Viral_Load < 200 & subset4$Age > 2, "PostART_Suppressed", NA)
)
# Subset again based on metadata condition
subset4 <- subset(subset4, !is.na(ART_Group))

plot_module_score_unpaired(
  seurat_obj = subset4,
  group_var = "ART_Group",
  output_path = file.path(module_score_dir, "PreART_vs_PostART_IFNscore.png"),
  plot_title = "IFN Pathway Activity: Pre vs Post ART"
)

## ------------------------
## 5. PreART vs PostART Suppressed vs PostART Unsuppressed
## ------------------------
subset5 <- subset(TARA_ALL, Condition == "HEI")
subset5$ART_3Group <- ifelse(
  subset5$Age <= 2, "PreART_Entry",
  ifelse(subset5$Viral_Load < 200 & subset5$Age > 2, "PostART_Suppressed",
         ifelse(subset5$Viral_Load >= 200 & subset5$Age > 2, "PostART_Unsuppressed", NA))
)

# Drop NAs
subset5 <- subset(subset5, !is.na(ART_3Group))

subset5$ART_3Group <- factor(
  subset5$ART_3Group,
  levels = c("PreART_Entry",  "PostART_Unsuppressed","PostART_Suppressed")
)

# Get metadata
meta <- subset5@meta.data

# Define comparisons
comparisons <- list(
  c("PreART_Entry", "PostART_Suppressed"),
  c("PreART_Entry", "PostART_Unsuppressed"),
  c("PostART_Suppressed", "PostART_Unsuppressed")
)


# Plot
gg <- ggplot(meta, aes(x = ART_3Group, y = IFN_score1, fill = ART_3Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.2, size = 0.4, alpha = 0.3) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.signif",       # use significance stars
    tip.length = 0.03,
    hide.ns = TRUE,           # hide non-significant comparisons
    size = 5
  ) +
  facet_wrap(~Manual_Annotation, scales = "free_y") +
  theme_minimal(base_size = 18) +
  xlab("") + ylab("IFN Module Score") +
  ggtitle("IFN Pathway Activity: Pre vs Post ART (Suppressed vs Unsuppressed)") +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1.5, "lines")  # increase space between facets
    
  )

gg

ggsave(
  filename = "/home/akshay-iyer/Documents/CD8_Longitudinal/IFNa_Analysis/Module_scores/PreART_vs_PostART_3groups_IFNscore.png",
  plot = gg,
  width = 30, height = 29, dpi = 400, bg = "white"
)

######################## CCR5 Corelation ########################


correlation_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/IFNa_Analysis/Correlation_CCR5"
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

#### Function


plot_ifn_vs_ccr5_by_cluster_grouped <- function(seurat_obj, score_var, gene_var, group_var, output_path, plot_title) {
  meta <- seurat_obj@meta.data
  expr <- FetchData(seurat_obj, vars = c(score_var, gene_var, group_var, "Manual_Annotation"))
  expr <- na.omit(expr)
  
  p <- ggplot(expr, aes_string(x = score_var, y = gene_var, color = group_var)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = FALSE, size = 0.6) +
    stat_cor(aes_string(color = group_var), method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    facet_wrap(~Manual_Annotation, scales = "free_y") +
    labs(title = plot_title, x = "IFN Module Score", y = paste0(gene_var, " Expression")) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(output_path, plot = p, width = 20, height = 18, dpi = 300, bg = "white")
}
### Plotting ###################
# mRNA
plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset1,
  score_var = "IFN_score1",
  gene_var = "RNA_CCR5",
  group_var = "Condition",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_mRNA_HEIvsHEU.png"),
  plot_title = "IFN Score vs CCR5 (mRNA): HEI vs HEU"
)

plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset2,
  score_var = "IFN_score1",
  gene_var = "RNA_CCR5",
  group_var = "Condition",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_mRNA_HEUvsHUU.png"),
  plot_title = "IFN Score vs CCR5 (mRNA): HEU vs HUU"
)

plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset3,
  score_var = "IFN_score1",
  gene_var = "RNA_CCR5",
  group_var = "VL_Group",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_mRNA_HighVLvsLowVL.png"),
  plot_title = "IFN Score vs CCR5 (mRNA): High vs Low Viral Load (Entry, HEI)"
)

plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset4,
  score_var = "IFN_score1",
  gene_var = "RNA_CCR5",
  group_var = "ART_Group",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_mRNA_PreARTvsPostART.png"),
  plot_title = "IFN Score vs CCR5 (mRNA): Pre vs Post ART (HEI)"
)

# Protein
plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset1,
  score_var = "IFN_score1",
  gene_var = "Protein_CCR5",
  group_var = "Condition",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_protein_HEIvsHEU.png"),
  plot_title = "IFN Score vs CCR5 (Protein): HEI vs HEU"
)

plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset2,
  score_var = "IFN_score1",
  gene_var = "Protein_CCR5",
  group_var = "Condition",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_protein_HEUvsHUU.png"),
  plot_title = "IFN Score vs CCR5 (Protein): HEU vs HUU"
)

plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset3,
  score_var = "IFN_score1",
  gene_var = "Protein_CCR5",
  group_var = "VL_Group",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_protein_HighVLvsLowVL.png"),
  plot_title = "IFN Score vs CCR5 (Protein): High vs Low Viral Load (Entry, HEI)"
)

plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset4,
  score_var = "IFN_score1",
  gene_var = "Protein_CCR5",
  group_var = "ART_Group",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_protein_PreARTvsPostART.png"),
  plot_title = "IFN Score vs CCR5 (Protein): Pre vs Post ART (HEI)"
)

### Subset 5 (3)
# mRNA
plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset5,
  score_var = "IFN_score1",
  gene_var = "RNA_CCR5",
  group_var = "ART_3Group",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_mRNA_PrePostSuppressedUnsuppressed.png"),
  plot_title = "IFN Score vs CCR5 (mRNA): Pre vs Post ART (Suppressed vs Unsuppressed)"
)

# Protein
plot_ifn_vs_ccr5_by_cluster_grouped(
  seurat_obj = subset5,
  score_var = "IFN_score1",
  gene_var = "Protein_CCR5",
  group_var = "ART_3Group",
  output_path = file.path(correlation_dir, "IFN_vs_CCR5_protein_PrePostSuppressedUnsuppressed.png"),
  plot_title = "IFN Score vs CCR5 (Protein): Pre vs Post ART (Suppressed vs Unsuppressed)"
)

load.path <- "~/Documents/CD8_Longitudinal/saved_R_data/"
saveRDS(TARA_ALL, file = file.path(load.path, "TARA_ALL_post_annotation_comparisonsplit.rds"))
