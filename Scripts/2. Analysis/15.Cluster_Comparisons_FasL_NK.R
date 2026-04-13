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

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "TARA_ALL_post_annotation.rds")

TARA_ALL <- readRDS(rds_in)
DefaultAssay(TARA_ALL) <- 'RNA'
assay_use   <- DefaultAssay(TARA_ALL)




################ Cluster Comparison Heatmaps ##################################

### Function
make_heatmap <- function(seu, clusters, outfile, n_top_genes = 10) {
  # Subset
  seu_sub <- subset(seu, subset = Manual_Annotation %in% clusters)
  
  # Get variable features
  genes <- VariableFeatures(seu_sub)
  
  # Calculate z-scores
  toplot <- CalcStats(seu_sub, features = genes, method = "zscore", order = "p", n = n_top_genes)
  
  # Heatmap
  p <- Heatmap(toplot, lab_fill = "zscore")
  
  # Save
  ggsave(outfile, p, width = 10, height = 6)
}

make_heatmap_adt <- function(seu, clusters, outfile, n_top_genes = 10) {
  # Subset
  seu_sub <- subset(seu, subset = Manual_Annotation %in% clusters)
  
  # Use ADT assay
  DefaultAssay(seu_sub) <- "ADT"
  
  # Use all proteins
  genes <- rownames(seu_sub[["ADT"]])
  
  # Calculate z-scores
  toplot <- CalcStats(
    seu_sub,
    features = genes,
    method   = "zscore",
    order    = "p",
    n        = n_top_genes
  )
  
  # Heatmap
  p <- Heatmap(toplot, lab_fill = "zscore")
  
  # Save
  ggsave(outfile, p, width = 10, height = 6)
}
### RNA

DefaultAssay(TARA_ALL) <-'RNA'

setwd ('/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Cluster_Comparison_Heatmaps/RNA')

# DN T Cells
make_heatmap(
  TARA_ALL,
  clusters = c("18: DN T cell_1","23: DN T cell_2","26: DN T cell_3","28: DN T cell_4"),
  outfile = "DN_T_Cells_Heatmap_TARA_ALL.png"
)

# CD8 
make_heatmap(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_All_Heatmap_TARA_ALL.png"
)

make_heatmap(
  TARA_ALL,
  clusters = c("8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_Only_Heatmap_TARA_ALL.png"
)

make_heatmap(
  TARA_ALL,
  clusters = c("6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_and_Naive_Heatmap_TARA_ALL.png"
)

make_heatmap(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell"),
  outfile = "CD8_MemoryvsNaive_Heatmap_TARA_ALL.png"
)

# NK ALL
make_heatmap(
  TARA_ALL,
  clusters = c("3: IL2RB+ NK cell","15: NK cell_1","17: NK cell_2","25: CD56bright NK"),
  outfile = "NK_All_Heatmap_TARA_ALL.png"
)

# B cells
make_heatmap(
  TARA_ALL,
  clusters = c("4: IgD+IgM+ B cell","11: Naïve B cell","13: Transitional B cell","22: CD69+ B cell","31: B cell","32: Plasmablast","35: TNF+ B cell"),
  outfile = "B_Cells_All_Heatmap_TARA_ALL.png"
)



### ADT
DefaultAssay(TARA_ALL) <-'ADT'


setwd ('/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/Cluster_Comparison_Heatmaps/ADT')

# DN T Cells
make_heatmap_adt(
  TARA_ALL,
  clusters = c("18: DN T cell_1","23: DN T cell_2","26: DN T cell_3","28: DN T cell_4"),
  outfile = "DN_T_Cells_Heatmap_TARA_ALL.png"
)

# CD8 
make_heatmap_adt(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_All_Heatmap_TARA_ALL.png"
)

make_heatmap_adt(
  TARA_ALL,
  clusters = c("8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_Only_Heatmap_TARA_ALL.png"
)

make_heatmap_adt(
  TARA_ALL,
  clusters = c("6: Naïve CD8 T cell","8: CTL-like","9: TRDV1+ CTL-like","27: GZMK+ CD8 T cell"),
  outfile = "CD8_CTL_and_Naive_Heatmap_TARA_ALL.png"
)

make_heatmap_adt(
  TARA_ALL,
  clusters = c("1: Memory CD8 T cell","6: Naïve CD8 T cell"),
  outfile = "CD8_MemoryvsNaive_Heatmap_TARA_ALL.png"
)

# NK ALL
make_heatmap_adt(
  TARA_ALL,
  clusters = c("3: IL2RB+ NK cell","15: NK cell_1","17: NK cell_2","25: CD56bright NK"),
  outfile = "NK_All_Heatmap_TARA_ALL.png"
)

# B cells
make_heatmap_adt(
  TARA_ALL,
  clusters = c("4: IgD+IgM+ B cell","11: Naïve B cell","13: Transitional B cell","22: CD69+ B cell","31: B cell","32: Plasmablast","35: TNF+ B cell"),
  outfile = "B_Cells_All_Heatmap_TARA_ALL.png"
)


########### Cluster Distribution Plots ##############
setwd ('/home/akshay-iyer/Documents/CD8_Longitudinal/Annotation/TARA_ALL/')

ClusterDistrPlot(
  origin = TARA_ALL$orig.ident,
  cluster = TARA_ALL$Manual_Annotation,
  condition = TARA_ALL$Condition
)
ggsave('Clust_Distribution_TARA_ALL_By_Condition.png', width = 15, height = 13)

ClusterDistrPlot(
  origin = TARA_entry$orig.ident, ### make tara entry first which is in the next section, needs reordering)
  cluster = TARA_entry$Manual_Annotation,
  condition = TARA_entry$Viral_Load_Category
)
ggsave('Clust_Distribution_TARA_Entry_HEI_By_Viral_Load.png', width = 15, height = 13)



######################### NK Cell FasL Analysis #################################

setwd ('/home/akshay-iyer/Documents/CD8_Longitudinal/FASL')

DimPlot2(TARA_ALL, features = c("FAS", "FASLG"),reduction = 'wnn.umap')
ggsave('FAS_FASLG_Featureplot.png', width = 12, height = 6)


VlnPlot2(TARA_ALL,'FAS',show.mean = T,cols = 'default',angle = 90)
ggsave('FAS_VLN_Plot_ALL.png', width = 12, height =7)

VlnPlot2(TARA_ALL,'FASLG',show.mean = T,cols = 'default',angle = 90)
ggsave('FASLG_VLN_Plot_ALL.png', width = 12, height = 7)

### Set viral load status
TARA_entry <- subset(TARA_ALL, subset = Age %in% c(1,2) & Condition == "HEI")
TARA_entry$Viral_Load_Category <- ifelse(TARA_entry$Viral_Load >= 100000, "High", "Low")

VlnPlot2(TARA_entry,c('FAS','FASLG'),show.mean = T,split.by = 'Viral_Load_Category',
         stat.method = "wilcox.test")
ggsave('FASL_FASLG_VLN_Plot_HEI_ENTRY_HIGHvsLOW.png', width = 17, height = 11)


# ---- SETTINGS ----------------------------------------------------------------
nk_clusters <- c("3: IL2RB+ NK cell","15: NK cell_1","17: NK cell_2","25: CD56bright NK")

assay_to_use <- "RNA"; fasl_feature <- "FASLG"   

# ---- 1) LABEL CELLS AS FASL+ / FASL- -----------------------------------------
DefaultAssay(TARA_ALL) <- assay_to_use

# Use normalized values (slot='data'). Threshold: > 0 (for RNA log-normalized or ADT CLR)
vals <- as.numeric(GetAssayData(TARA_ALL, assay = assay_to_use, slot = "data")[fasl_feature, ])
TARA_ALL$FASL_expr <- vals
TARA_ALL$FASL_status <- ifelse(TARA_ALL$FASL_expr > 0, "FASL_pos", "FASL_neg")


fasl_pct_all <- TARA_ALL@meta.data %>%
  count(Manual_Annotation, FASL_status) %>%
  group_by(Manual_Annotation) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# Save as CSV
write.csv(fasl_pct_all, "FASL_Positive_Fractions_byCluster.csv", row.names = FALSE)

# order clusters by % FASL_pos for readability
order_df <- fasl_pct_all %>%
  filter(FASL_status == "FASL_pos") %>%
  arrange(pct)

fasl_pct_all$Manual_Annotation <- factor(
  fasl_pct_all$Manual_Annotation,
  levels = order_df$Manual_Annotation
)

fasl_pct_all <- fasl_pct_all %>%
  filter(!is.na(Manual_Annotation), FASL_status %in% c("FASL_neg","FASL_pos")) %>%
  complete(Manual_Annotation, FASL_status = c("FASL_neg","FASL_pos"),
           fill = list(n = 0, pct = 0)) %>%
  droplevels()
p_stack <- ggplot(fasl_pct_all,
                  aes(x = Manual_Annotation, y = pct/100, fill = FASL_status)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  labs(x = NULL, y = "Fraction of cells", fill = "FASL status",
       title = "FASL+ / FASL− fractions by cluster") +
  theme_bw()


ggsave("FASL_fraction_byCluster_stacked.png", p_stack, width = 8, height = 10, dpi = 300)


########################### DGE and DPE ########################################
library(dplyr)
library(ggplot2)
# Optional: label a few points nicely if you have ggrepel installed
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# ----------------------- SETTINGS --------------------------------------------
clusters_fasl <- c(
  "3: IL2RB+ NK cell",
  "21: Gamma Delta 2 T cells",
  "30: DN T cell_6",
  "25: CD56bright NK",
  "15: NK cell_1",
  "17: NK cell_2",
  "9: TRDV1+ CTL-like",
  "27: GZMK+ CD8 T cell",
  "8: CTL-like"
)

base_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/FASL"
dge_dir  <- file.path(base_dir, "DGE")  # RNA
dpe_dir  <- file.path(base_dir, "DPE")  # ADT
dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dpe_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------- HELPERS ---------------------------------------------
run_de <- function(seurat_object, clusters, assay = "RNA", latent_var = "nCount_RNA",
                   pos_label = "FASL_pos", neg_label = "FASL_neg") {
  DefaultAssay(seurat_object) <- assay
  de_list <- list()
  
  for (cl in clusters) {
    de <- FindMarkers(
      subset(seurat_object, subset = Manual_Annotation == cl),
      ident.1        = pos_label,
      ident.2        = neg_label,
      test.use       = "MAST",
      group.by       = "FASL_status",
      latent.vars    = latent_var,
      logfc.threshold= 0
    )
    de$cluster <- cl
    de_list[[cl]] <- de
  }
  
  do.call(rbind, de_list)
}
write_de <- function(de_table, outdir, prefix) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Preserve original rownames (which currently look like "<cluster>.<gene>")
  rn <- rownames(de_table)
  df <- de_table
  df$._rn <- rn  # keep for debugging if needed
  
  # Vectorized extraction of gene by removing the "<cluster>." prefix ONLY if present
  # Works even when cluster has spaces/colons/plus signs, and gene may contain dots.
  gene <- ifelse(
    startsWith(rn, paste0(df$cluster, ".")),
    substring(rn, nchar(df$cluster) + 2L),  # +1 for the dot and +1 for 1-based indexing
    rn
  )
  df$gene <- gene
  
  # Reorder columns: gene first, then the stats, then cluster
  # Keep a clean table for output
  keep_order <- c("gene", setdiff(colnames(df), c("gene", "._rn")))
  df <- df[, keep_order, drop = FALSE]
  
  # Sort for readability
  df <- df[order(df$cluster, df$p_val_adj, -abs(df$avg_log2FC)), ]
  
  # Combined CSV (no row names)
  write.csv(df, file.path(outdir, paste0(prefix, "_ALL.csv")), row.names = FALSE)
  
  # Per-cluster CSVs
  split(df, df$cluster) |>
    lapply(function(dd) {
      fn <- file.path(outdir, paste0(prefix, "_", make.names(unique(dd$cluster)), ".csv"))
      write.csv(dd, fn, row.names = FALSE)
      invisible(NULL)
    })
}


plot_volcano <- function(de_table, outdir, padj_thr = 0.05, lfc_thr = 0.25, top_n = 10) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Ensure we have a gene column (if write_de ran first, it’s there already)
  df_all <- de_table
  if (!"gene" %in% colnames(df_all)) {
    rn <- rownames(df_all)
    df_all$gene <- ifelse(
      startsWith(rn, paste0(df_all$cluster, ".")),
      substring(rn, nchar(df_all$cluster) + 2L),
      rn
    )
  }
  
  for (cl in unique(df_all$cluster)) {
    df <- subset(df_all, cluster == cl)
    df$neglog10_padj <- -log10(df$p_val_adj + 1e-300)
    df$sig <- ifelse(df$p_val_adj < padj_thr & df$avg_log2FC >  lfc_thr, "Up in FASL+",
                     ifelse(df$p_val_adj < padj_thr & df$avg_log2FC < -lfc_thr, "Down in FASL+", "NS"))
    
    # pick top labels on each side by significance
    lab_up <- subset(df, p_val_adj < padj_thr & avg_log2FC >  lfc_thr)
    lab_dn <- subset(df, p_val_adj < padj_thr & avg_log2FC < -lfc_thr)
    lab_up <- lab_up[order(-lab_up$neglog10_padj), ][seq_len(min(nrow(lab_up), top_n)), , drop = FALSE]
    lab_dn <- lab_dn[order(-lab_dn$neglog10_padj), ][seq_len(min(nrow(lab_dn), top_n)), , drop = FALSE]
    lab_genes <- unique(c(lab_up$gene, lab_dn$gene))
    
    p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj, color = sig)) +
      geom_point(size = 1.5) +
      geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
      geom_hline(yintercept = -log10(padj_thr), linetype = 2) +
      scale_color_manual(values = c("Up in FASL+" = "#D55E00", "Down in FASL+" = "#0072B2", "NS" = "grey")) +
      labs(title = paste("Volcano:", cl),
           x = "avg_log2FC (FASL+ vs FASL-)",
           y = "-log10(adj p)") +
      theme_bw() +
      geom_text_repel(
        data = subset(df, gene %in% lab_genes),
        aes(label = gene),
        max.overlaps = Inf,
        size = 3
      )
    
    ggsave(file.path(outdir, paste0("Volcano_", make.names(cl), ".png")),
           p, width = 7, height = 5, dpi = 300)
  }
}

    
# RNA
de_rna <- run_de(TARA_ALL, clusters_fasl, assay = "RNA", latent_var = "nCount_RNA")
write_de(de_rna, dge_dir, "DGE_RNA")
plot_volcano(de_rna, file.path(dge_dir, "Volcano"))

# ADT
de_adt <- run_de(TARA_ALL, clusters_fasl, assay = "ADT", latent_var = "nCount_ADT")
write_de(de_adt, dpe_dir, "DPE_ADT")
plot_volcano(de_adt, file.path(dpe_dir, "Volcano"))

colnames(TARA_ALL@assays$ADT@layers$data[1])
###################3 Plot specific volcano plots ##################
# ---- Setup ----
library(tidyverse)
library(ggrepel)

# Input file
infile <- "/home/akshay-iyer/Documents/CD8_Longitudinal/FASL/DPE/DPE_ADT_X3..IL2RB..NK.cell.csv"

# Output dir (will be created if it doesn't exist)
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/FASL/DPE/Volcano"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Thresholds
padj_thr <- 0.05
lfc_thr  <- 0.25   # change if you prefer another log2FC cutoff

# Labels to highlight (case-insensitive match; aliases handled below)
lab_genes <- c(
  "FCGR3A",
  "C5AR1",
  "GGT1",
  "KLRD1",
  "SLAMF7",
  "IL2RB",
  "KLRB1",
  "CD81",
  "CD226",
  "TIGIT",
  "CX3CR1",
  "KIR3DL1",
  "TNFRSF14",
  "KIR2DL3",
  "PECAM1",
  "ICOS",
  "ITGAM",
  "FAS"
)



# ---- Load & prep data ----
df <- readr::read_csv(infile, show_col_types = FALSE)

# Ensure required columns exist
req_cols <- c("gene","p_val_adj","avg_log2FC","cluster")
stopifnot(all(req_cols %in% names(df)))

# Clean/augment
df <- df %>%
  mutate(
    neglog10_padj = -log10(p_val_adj),
    sig = case_when(
      p_val_adj < padj_thr & avg_log2FC >=  lfc_thr ~ "Up in FASL+",
      p_val_adj < padj_thr & avg_log2FC <= -lfc_thr ~ "Down in FASL+",
      TRUE ~ "NS"
    )
  )
df <- df %>% filter(between(avg_log2FC, -10, 10))
message("Post-filter FC range: ", paste(range(df$avg_log2FC, na.rm = TRUE), collapse = " to "))

# Cluster name for titles/filenames (use unique cluster in file)
cl <- df %>% dplyr::pull(cluster) %>% unique() %>% paste(collapse = "; ")

# Case-insensitive matching for labels, also check alias names
gene_lower <- tolower(df$gene)
want_lower <- tolower(lab_genes)
label_idx  <- gene_lower %in% want_lower

# Report any requested labels not found (informative message only)
not_found <- setdiff(tolower(lab_genes), gene_lower)
if (length(not_found) > 0) {
  message("Note: Some requested labels not found in 'gene' column (case-insensitive): ",
          paste(sort(unique(not_found)), collapse = ", "))
}

# ---- Plot ----
p <- ggplot(df, aes(x = avg_log2FC, y = neglog10_padj, color = sig)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
  geom_hline(yintercept = -log10(padj_thr), linetype = 2) +
  scale_color_manual(values = c("Up in FASL+" = "#D55E00",
                                "Down in FASL+" = "#0072B2",
                                "NS" = "grey70")) +
  labs(
    title = paste0("Volcano: ", cl),
    x = "avg_log2FC (FASL+ vs FASL-)",
    y = "-log10(adj p)",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  geom_label_repel(
    data = df[label_idx, ],
    aes(label = gene),
    max.overlaps = Inf,
    size = 3,
    min.segment.length = 0,
    box.padding = 0.4,
    label.size = 0.3,          # thickness of the box border
    fill = "white",            # background color of box
    alpha = 0.9                # slight transparency if you want
  )

'#36A374'
# ---- Save ----


ggsave(paste0(out_dir,'/IL2RB_NK_volcano_specified_prot.png'), p, width =11, height = 7, dpi = 300)

