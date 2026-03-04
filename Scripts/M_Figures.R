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
library(tidyr)
library(Matrix)
library(ggrepel)
library(pheatmap)
library(patchwork)
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

clusters_use <- c("1: Memory CD8 T cell",
                  "2: Naïve CD4 T cell_1",
                  "8: CTL-like")

meta_df <- as.data.frame(TARA_ALL@meta.data)

meta_sub <- meta_df %>%
  filter(as.character(Timepoint_Group) %in% c("HUU", "HEU", "PreART_Entry")) %>%
  mutate(
    Condition = case_when(
      as.character(Timepoint_Group) == "PreART_Entry" ~ "HEI",
      TRUE ~ as.character(Timepoint_Group)
    ),
    Condition = factor(Condition, levels = c("HUU", "HEU", "HEI"))
  )

cluster_freq <- meta_sub %>%
  group_by(orig.ident, Condition, Manual_Annotation) %>%
  summarise(cluster_cells = dplyr::n(), .groups = "drop") %>%
  group_by(orig.ident, Condition) %>%
  mutate(total_cells = sum(cluster_cells)) %>%
  ungroup() %>%
  mutate(percent = 100 * cluster_cells / total_cells) %>%
  filter(Manual_Annotation %in% clusters_use)

# Make sure Condition order is set (you already did this, but safe)
# Ensure correct factor order
cluster_freq <- cluster_freq %>%
  mutate(Condition = factor(Condition, levels = c("HUU","HEU","HEI")))

# Publication colors
cond_cols <- c("HUU" = "#4E79A7",
               "HEU" = "#F28E2B",
               "HEI" = "#E15759")

# Pairwise Wilcoxon (RAW p-values)
stat_df <- cluster_freq %>%
  group_by(Manual_Annotation) %>%
  wilcox_test(percent ~ Condition) %>%
  add_xy_position(
    x = "Condition",
    step.increase = 0.05   
  ) %>%
  mutate(p.signif = case_when(
    p <= 0.0001 ~ "****",
    p <= 0.001  ~ "***",
    p <= 0.01   ~ "**",
    p <= 0.05   ~ "*",
    TRUE        ~ "ns"
  )) %>%
  filter(p.signif != "ns")
ggplot(cluster_freq,
       aes(x = Condition, y = percent,
           fill = Condition, color = Condition)) +
  
  geom_boxplot(width = 0.65,
               outlier.shape = NA,
               alpha = 0.25,
               linewidth = 0.7) +
  
  geom_jitter(width = 0.12,
              size = 2.3,
              alpha = 0.8) +
  
  facet_wrap(~ Manual_Annotation, scales = "free_y") +
  
  stat_pvalue_manual(stat_df,
                     label = "p.signif",
                     tip.length = 0.01,
                     bracket.size = 0.6,
                     size = 4.5) +
  
  scale_fill_manual(values = cond_cols) +
  scale_color_manual(values = cond_cols) +
  
  labs(x = NULL,
       y = "% of total PBMC") +
  
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    panel.spacing = unit(1.2, "lines")
  )
output_path_png <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 1/02_Cluster_Proportions/Cluster_Proportions_HEI_HEU_HUU.png"

ggsave(
  filename = output_path_png,
  plot = last_plot(),
  width = 11,
  height = 6,
  dpi = 600
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

############################################################
#### FIGURE 2 — Viral Load vs Cluster Abundance (PreART_Entry) ##########
############################################################

# -----------------------------
# Output directory + structure
# -----------------------------
fig2_root <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 2"
fig2_dir  <- file.path(fig2_root, "Cluster_Size_Vs_Viral_Load")
dir.create(fig2_dir, showWarnings = FALSE, recursive = TRUE)

dir_data   <- file.path(fig2_dir, "01_DataTables")
dir_corr   <- file.path(fig2_dir, "02_Correlation")
dir_plots  <- file.path(fig2_dir, "03_ScatterPlots")
dir_allpct <- file.path(dir_plots, "All_Clusters", "Percent")
dir_allcnt <- file.path(dir_plots, "All_Clusters", "Counts")
dir_sig    <- file.path(dir_plots, "Significant_Only")
dir_panel  <- file.path(fig2_dir, "04_FinalPanels")

dir.create(dir_data,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_corr,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_allpct, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_allcnt, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_sig,    showWarnings = FALSE, recursive = TRUE)
dir.create(dir_panel,  showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Subset to PreART_Entry only
# -----------------------------
seu_pre <- subset(
  TARA_ALL,
  subset = Timepoint_Group == "PreART_Entry"
)

# Make viral load numeric
seu_pre$Viral_Load_numeric <- as.numeric(seu_pre$Viral_Load)

# Extract metadata
meta_pre <- seu_pre@meta.data

# Compute cluster counts per sample
#  Viral load per sample
vl_df <- meta_pre %>%
  group_by(orig.ident) %>%
  summarise(
    Viral_Load = mean(Viral_Load_numeric, na.rm = TRUE),
    .groups = "drop"
  )

#  Cluster counts
cluster_counts <- meta_pre %>%
  group_by(orig.ident, Manual_Annotation) %>%
  summarise(cluster_cells = n(), .groups = "drop")

#  Add totals + percent
cluster_freq <- cluster_counts %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells = sum(cluster_cells),
    percent = (cluster_cells / total_cells) * 100
  ) %>%
  ungroup() %>%
  left_join(vl_df, by = "orig.ident")
# Add log10 viral load (avoid -Inf if any VL == 0)
cluster_freq <- cluster_freq %>%
  mutate(
    Viral_Load_log10 = ifelse(is.finite(Viral_Load) & Viral_Load > 0,
                              log10(Viral_Load),
                              NA_real_)
  )

readr::write_csv(
  cluster_freq,
  file.path(dir_data, "PreART_Entry_Cluster_Frequencies_with_ViralLoad.csv")
)



##Correlation Per Cluster

safe_spearman <- function(x, y, min_n = 5) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  
  n_ok <- length(x2)
  if (n_ok < min_n) {
    return(list(rho = NA_real_, p = NA_real_, n = n_ok))
  }
  # If no variation in x or y, correlation test is undefined
  if (length(unique(x2)) < 2 || length(unique(y2)) < 2) {
    return(list(rho = NA_real_, p = NA_real_, n = n_ok))
  }
  
  ct <- suppressWarnings(cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = n_ok)
}

cor_results <- cluster_freq %>%
  group_by(Manual_Annotation) %>%
  summarise(
    tmp_percent = list(safe_spearman(percent, Viral_Load_log10, min_n = 5)),
    tmp_counts  = list(safe_spearman(cluster_cells, Viral_Load_log10, min_n = 5)),
    cor_percent = tmp_percent[[1]]$rho,
    p_percent   = tmp_percent[[1]]$p,
    n_percent   = tmp_percent[[1]]$n,
    cor_counts  = tmp_counts[[1]]$rho,
    p_counts    = tmp_counts[[1]]$p,
    n_counts    = tmp_counts[[1]]$n,
    .groups = "drop"
  ) %>%
  select(-tmp_percent, -tmp_counts) %>%
  arrange(p_percent)

readr::write_csv(
  cor_results,
  file.path(dir_corr, "PreART_Entry_ViralLoad_LOG10_Correlation_by_Cluster.csv")
)

###Scatter Plots
### Percent of PBMC vs Viral Load
for (cl in unique(cluster_freq$Manual_Annotation)) {
  
  df_plot <- cluster_freq %>% filter(Manual_Annotation == cl)
  
  p <- ggplot(df_plot, aes(x = Viral_Load_log10, y = percent)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(
      title = paste0(cl, " (% of PBMC)"),
      x = "log10(Viral Load)",
      y = "% of Total PBMC"
    ) +
    theme_classic(base_size = 16)
  
  ggsave(
    filename = file.path(dir_allpct, paste0("PreART_Entry_log10VL_vs_PERCENT_", cl, ".png")),
    plot = p, width = 6, height = 5, dpi = 300
  )
}


### Raw Counts vs Viral Load
for (cl in unique(cluster_freq$Manual_Annotation)) {
  
  df_plot <- cluster_freq %>% filter(Manual_Annotation == cl)
  
  p <- ggplot(df_plot, aes(x = Viral_Load_log10, y = cluster_cells)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(
      title = paste0(cl, " (Raw Count)"),
      x = "log10(Viral Load)",
      y = "Cell Count"
    ) +
    theme_classic(base_size = 16)
  
  ggsave(
    filename = file.path(dir_allcnt, paste0("PreART_Entry_log10VL_vs_COUNTS_", cl, ".png")),
    plot = p, width = 6, height = 5, dpi = 300
  )
}

#########Significant-only panel uses log10 viral load
alpha <- 0.05
min_n <- 5
top_n <- 6

sig_tbl <- cor_results %>%
  filter(!is.na(p_percent), n_percent >= min_n, p_percent < alpha) %>%
  arrange(p_percent)

readr::write_csv(sig_tbl, file.path(dir_corr, "SignificantClusters_PERCENT_log10VL.csv"))

sig_clusters <- head(sig_tbl$Manual_Annotation, top_n)

sig_plot_list <- list()

for (cl in sig_clusters) {
  
  df_plot <- cluster_freq %>% filter(Manual_Annotation == cl)
  
  p <- ggplot(df_plot, aes(x = Viral_Load_log10, y = percent)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE) +
    ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(
      title = cl,
      x = "log10(Viral Load)",
      y = "% of Total PBMC"
    ) +
    theme_classic(base_size = 16)
  
  ggsave(
    filename = file.path(dir_sig, paste0("SIG_PreART_Entry_log10VL_vs_PERCENT_", cl, ".png")),
    plot = p, width = 6, height = 5, dpi = 300
  )
  
  sig_plot_list[[cl]] <- p
}

if (length(sig_plot_list) > 0) {
  panel <- patchwork::wrap_plots(sig_plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = "PreART_Entry: Cluster frequency vs log10(Viral Load) (Significant clusters)"
    )
  
  ggsave(
    filename = file.path(dir_panel, "Fig2_PreART_Entry_SignificantClusters_PERCENT_log10VL_Panel.png"),
    plot = panel, width = 12, height = 10, dpi = 300
  )
  
  ggsave(
    filename = file.path(dir_panel, "Fig2_PreART_Entry_SignificantClusters_PERCENT_log10VL_Panel.pdf"),
    plot = panel, width = 12, height = 10, device = cairo_pdf
  )
}


ClusterDistrPlot(
  origin = seu_pre$orig.ident,
  cluster = seu_pre$Manual_Annotation,
  condition = seu_pre$Gender
)
############################ CLECL4C Correlation ############
# PBMC-level folder
dir_pbmc_level <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 2/Cluster_Size_Vs_Viral_Load/05_VL_Associated_Features/01_PBMC_Level"
dir_pbmc_plots <- file.path(dir_pbmc_level, "03_Plots")
dir_pbmc_res   <- file.path(dir_pbmc_level, "02_Results")
dir.create(dir_pbmc_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pbmc_res,   showWarnings = FALSE, recursive = TRUE)

# Match sample IDs: rna_avg/adt_avg columns are like "CP002-entry"
# Ensure vl_df_use uses hyphen too
vl_df_use_match <- vl_df_use %>%
  mutate(orig.ident = gsub("_", "-", orig.ident)) %>%
  filter(is.finite(Viral_Load_log10))

# ---- find CLEC4C row robustly ----
if (!"CLEC4C" %in% rownames(adt_avg)) {
  hit <- grep("^CLEC4C$|CLEC4C", rownames(adt_avg), value = TRUE)
  if (length(hit) == 0) stop("CLEC4C not found in rownames(adt_avg). Check ADT feature names.")
  clec_name <- hit[1]
} else {
  clec_name <- "CLEC4C"
}

# ---- build plotting df (align samples) ----
samp <- intersect(colnames(adt_avg), vl_df_use_match$orig.ident)

df_clec <- tibble(
  orig.ident = samp,
  Viral_Load_log10 = vl_df_use_match$Viral_Load_log10[match(samp, vl_df_use_match$orig.ident)],
  CLEC4C_ADT = as.numeric(adt_avg[clec_name, samp])
) %>% filter(is.finite(Viral_Load_log10), is.finite(CLEC4C_ADT))

# ---- plot ----
p_clec <- ggplot(df_clec, aes(x = Viral_Load_log10, y = CLEC4C_ADT)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = paste0("PBMC-level ADT: ", clec_name, " vs log10(Viral Load)"),
    x = "log10(Viral Load)",
    y = paste0(clec_name, " (ADT, DSB normalized)")
  ) +
  theme_classic(base_size = 16)

ggsave(
  filename = file.path(dir_pbmc_plots, "PBMC_ADT_CLEC4C_vs_log10VL.png"),
  plot = p_clec, width = 6.5, height = 5.5, dpi = 300
)

# save df used
readr::write_csv(df_clec, file.path(dir_pbmc_res, "PBMC_ADT_CLEC4C_vs_log10VL_data.csv"))

# Curated Type I IFN/ISG set (edit if you have your own)
ifn1_genes <- c(
  "ISG15","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM3",
  "MX1","MX2","OAS1","OAS2","OASL","RSAD2","USP18",
  "IRF7","STAT1","DDX58","IFIH1","HERC5","BST2"
)

genes_present <- intersect(rownames(rna_avg), ifn1_genes)
if (length(genes_present) < 5) {
  stop(paste0("Too few IFN genes found in rna_avg (found ", length(genes_present), "). Check gene symbols / rownames(rna_avg)."))
}

samp_rna <- intersect(colnames(rna_avg), vl_df_use_match$orig.ident)

# Option: z-score each gene across samples, then average (recommended)
rna_sub <- as.matrix(rna_avg[genes_present, samp_rna, drop = FALSE])
rna_z   <- t(scale(t(rna_sub)))  # gene-wise z-score across samples

ifn1_score <- colMeans(rna_z, na.rm = TRUE)

df_ifn <- tibble(
  orig.ident = samp_rna,
  Viral_Load_log10 = vl_df_use_match$Viral_Load_log10[match(samp_rna, vl_df_use_match$orig.ident)],
  IFN1_Module_Zmean = as.numeric(ifn1_score)
) %>% filter(is.finite(Viral_Load_log10), is.finite(IFN1_Module_Zmean))

readr::write_csv(df_ifn, file.path(dir_pbmc_res, "PBMC_RNA_IFN1_Module_vs_log10VL_data.csv"))

p_ifn <- ggplot(df_ifn, aes(x = Viral_Load_log10, y = IFN1_Module_Zmean)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = paste0("PBMC-level RNA: Type I IFN Module vs log10(Viral Load)\n(n genes=", length(genes_present), ")"),
    x = "log10(Viral Load)",
    y = "Type I IFN module score (mean z-scored expression)"
  ) +
  theme_classic(base_size = 16)

ggsave(
  filename = file.path(dir_pbmc_plots, "PBMC_RNA_IFN1_Module_vs_log10VL.png"),
  plot = p_ifn, width = 6.5, height = 5.5, dpi = 300
)
#### Per cluster
dir_cluster_level <- file.path(dir_pbmc_level, "04_Cluster_Level")
dir_cluster_plots <- file.path(dir_cluster_level, "03_Plots")
dir_cluster_res   <- file.path(dir_cluster_level, "02_Results")

dir.create(dir_cluster_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_cluster_res, showWarnings = FALSE, recursive = TRUE)


DefaultAssay(seu_pre2) <- "ADT"

clec_name <- if ("CLEC4C" %in% rownames(seu_pre2)) "CLEC4C" else grep("CLEC4C", rownames(seu_pre2), value = TRUE)[1]

# Extract per-cell data
adt_df <- FetchData(seu_pre2, vars = c("orig.ident", "Manual_Annotation", clec_name))

colnames(adt_df)[3] <- "CLEC4C_ADT"

# Aggregate per cluster per sample
adt_cluster_avg <- adt_df %>%
  group_by(orig.ident, Manual_Annotation) %>%
  summarise(
    CLEC4C_mean = mean(CLEC4C_ADT, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# Merge VL
adt_cluster_avg <- adt_cluster_avg %>%
  mutate(orig.ident = gsub("_", "-", orig.ident)) %>%
  left_join(vl_df_use_match, by = "orig.ident")

write.csv(adt_cluster_avg,
          file.path(dir_cluster_res, "CLEC4C_ADT_ClusterLevel_Data.csv"),
          row.names = FALSE)


min_n_samples <- 5

clusters_to_plot <- adt_cluster_avg %>%
  filter(is.finite(CLEC4C_mean), is.finite(Viral_Load_log10)) %>%
  group_by(Manual_Annotation) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  filter(n_samples >= min_n_samples) %>%
  pull(Manual_Annotation)

for (cl in clusters_to_plot) {
  
  df_plot <- adt_cluster_avg %>%
    filter(Manual_Annotation == cl,
           is.finite(CLEC4C_mean),
           is.finite(Viral_Load_log10))
  
  p <- ggplot(df_plot, aes(x = Viral_Load_log10, y = CLEC4C_mean)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(
      title = paste0("CLEC4C ADT vs log10(VL) — ", cl),
      x = "log10(Viral Load)",
      y = "Mean CLEC4C (ADT)"
    ) +
    theme_classic(base_size = 15)
  
  ggsave(
    file.path(dir_cluster_plots,
              paste0("CLEC4C_ADT_vs_log10VL_", gsub("[/: ]+", "_", cl), ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}






################## Build sample-level expression matrices (RNA + ADT) at PreART_Entry #########

####### Create Dir 

dir_expr <- file.path(fig2_dir, "05_VL_Associated_Features")
dir_pbmc <- file.path(dir_expr, "01_PBMC_Level")
dir_clst <- file.path(dir_expr, "02_Cluster_Within_Sample")

dir.create(dir_expr, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pbmc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_clst, showWarnings = FALSE, recursive = TRUE)

dir_pbmc_data  <- file.path(dir_pbmc, "01_Data")
dir_pbmc_res   <- file.path(dir_pbmc, "02_Results")
dir_pbmc_plots <- file.path(dir_pbmc, "03_Plots")
dir_clst_res   <- file.path(dir_clst, "02_Results")
dir_clst_plots <- file.path(dir_clst, "03_Plots")

dir.create(dir_pbmc_data,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pbmc_res,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pbmc_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_clst_res,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_clst_plots, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Cluster-within-sample unified folders
# -----------------------------
dir_clst_rna <- file.path(dir_clst, "RNA")
dir_clst_adt <- file.path(dir_clst, "ADT")

dir_rna_avg   <- file.path(dir_clst_rna, "01_AvgExpr")
dir_rna_res   <- file.path(dir_clst_rna, "02_Results")
dir_rna_plots <- file.path(dir_clst_rna, "03_Plots")

dir_adt_avg   <- file.path(dir_clst_adt, "01_AvgExpr")
dir_adt_res   <- file.path(dir_clst_adt, "02_Results")
dir_adt_plots <- file.path(dir_clst_adt, "03_Plots")

dir.create(dir_rna_avg,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rna_res,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rna_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_adt_avg,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_adt_res,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_adt_plots, showWarnings = FALSE, recursive = TRUE)

### Func - Spearman per feature + volcano-like plot


safe_spearman_vec <- function(x, y, min_n = 5) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]; y2 <- y[ok]
  n_ok <- length(x2)
  if (n_ok < min_n) return(c(rho = NA_real_, p = NA_real_, n = n_ok))
  if (length(unique(x2)) < 2 || length(unique(y2)) < 2) return(c(rho = NA_real_, p = NA_real_, n = n_ok))
  ct <- suppressWarnings(cor.test(x2, y2, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = n_ok)
}

correlate_matrix_vs_vl <- function(mat_features_by_sample, vl_df, vl_col = "Viral_Load_log10", min_n = 5) {
  
  samp <- intersect(colnames(mat_features_by_sample), vl_df$orig.ident)
  mat  <- mat_features_by_sample[, samp, drop = FALSE]
  y    <- vl_df[[vl_col]][match(samp, vl_df$orig.ident)]
  
  res <- t(apply(mat, 1, function(x) safe_spearman_vec(as.numeric(x), y, min_n = min_n)))
  res <- as.data.frame(res)
  res$feature <- rownames(res)
  
  res <- res %>%
    mutate(
      rho = as.numeric(rho),
      p = as.numeric(p),
      n = as.integer(n),
      neglog10p = -log10(p)
    ) %>%
    arrange(p)
  
  res
}

plot_corr_volcano <- function(df, title, out_png, label_top = 12, alpha = 0.05) {
  df2 <- df %>% filter(is.finite(rho), is.finite(p))
  if (nrow(df2) == 0) return(NULL)
  
  # label top by p-value
  lab <- df2 %>% arrange(p) %>% slice_head(n = label_top)
  
  p1 <- ggplot(df2, aes(x = rho, y = -log10(p))) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    ggrepel::geom_text_repel(data = lab, aes(label = feature), size = 3, max.overlaps = Inf) +
    labs(title = title, x = "Spearman rho", y = "-log10(p)") +
    theme_classic(base_size = 16)
  
  ggsave(out_png, p1, width = 8.5, height = 6.5, dpi = 300)
  p1
}


###### PBMC-level (sample aggregated) RNA + ADT correlations
# Ensure we use only samples with finite log10(VL)
vl_df_use <- vl_df %>% filter(is.finite(Viral_Load_log10))
cells_keep <- rownames(seu_pre@meta.data)[seu_pre$orig.ident %in% vl_df_use$orig.ident]
seu_pre2 <- subset(seu_pre, cells = cells_keep)

# RNA (log-normalized data slot)
DefaultAssay(seu_pre2) <- "RNA"
rna_avg <- AverageExpression(
  seu_pre2,
  group.by = "orig.ident",
  assays = "RNA",
  slot = "data",
  verbose = FALSE
)$RNA

write.csv(as.data.frame(rna_avg),
          file.path(dir_pbmc_data, "PBMC_RNA_AverageExpression_bySample.csv"))

# ADT (normalized data slot, e.g. CLR/DSB)

adt_test <- AverageExpression(
  seu_pre2,
  group.by = "orig.ident",
  assays = "ADT",
  slot = "data",
  verbose = FALSE
)

adt_avg <- adt_test$ADT

write.csv(
  as.data.frame(adt_avg),
  file.path(dir_pbmc_data, "PBMC_ADT_AverageExpression_bySample.csv")
)

vl_df_use <- vl_df_use %>%
  mutate(orig.ident = gsub("_", "-", orig.ident))

# RNA correlations
rna_cor <- correlate_matrix_vs_vl(rna_avg, vl_df_use, vl_col = "Viral_Load_log10", min_n = 5)
write_csv(rna_cor, file.path(dir_pbmc_res, "PBMC_RNA_FeatureCorrelation_vs_log10VL.csv"))
plot_corr_volcano(
  rna_cor,
  title = "PBMC-level RNA features vs log10(VL) (PreART_Entry)",
  out_png = file.path(dir_pbmc_plots, "PBMC_RNA_CorrVolcano_vs_log10VL.png"),
  label_top = 15,
  alpha = 0.05
)

# ADT correlations
if (exists("adt_avg")) {
  adt_cor <- correlate_matrix_vs_vl(adt_avg, vl_df_use, vl_col = "Viral_Load_log10", min_n = 5)
  write_csv(adt_cor, file.path(dir_pbmc_res, "PBMC_ADT_FeatureCorrelation_vs_log10VL.csv"))
  plot_corr_volcano(
    adt_cor,
    title = "PBMC-level ADT features vs log10(VL) (PreART_Entry)",
    out_png = file.path(dir_pbmc_plots, "PBMC_ADT_CorrVolcano_vs_log10VL.png"),
    label_top = 15,
    alpha = 0.05
  )
}


#### Cluster-within-sample correlations (RNA + ADT)
run_cluster_within_sample <- function(seu, cluster_col = "Manual_Annotation", cl, vl_df_use,
                                      assay_name = "RNA", slot_name = "data",
                                      min_cells_per_sample = 30, min_n_samples = 5,
                                      out_avg_dir, out_res_dir, out_plot_dir) {
  
  # sanitize cluster name for filenames
  cl_tag <- gsub("[^A-Za-z0-9]+", "_", cl)
  
  # keep cells in this cluster
  cells_cl <- rownames(seu@meta.data)[seu@meta.data[[cluster_col]] == cl]
  if (length(cells_cl) == 0) return(NULL)
  
  seu_cl <- subset(seu, cells = cells_cl)
  DefaultAssay(seu_cl) <- assay_name
  
  # cell counts per sample within cluster (explicit dplyr::count to avoid masking)
  cc <- seu_cl@meta.data %>%
    dplyr::count(orig.ident, name = "n_cells_cluster")
  
  keep_samples <- cc %>%
    filter(n_cells_cluster >= min_cells_per_sample) %>%
    pull(orig.ident)
  
  if (length(keep_samples) < min_n_samples) {
    return(NULL)
  }
  
  seu_cl <- subset(seu_cl, subset = orig.ident %in% keep_samples)
  
  # average expression within cluster per sample
  avg <- AverageExpression(
    seu_cl,
    group.by = "orig.ident",
    assays = assay_name,
    slot = slot_name,
    verbose = FALSE
  )[[assay_name]]
  
  # correlate features vs VL (only those samples)
  vl_sub <- vl_df_use %>% filter(orig.ident %in% colnames(avg))
  res <- correlate_matrix_vs_vl(avg, vl_sub, vl_col = "Viral_Load_log10", min_n = min_n_samples)
  
  # save (unique filenames per cluster)
  write.csv(as.data.frame(avg), file.path(out_avg_dir, paste0(cl_tag, "_", assay_name, "_AvgExpr_bySample.csv")))
  readr::write_csv(res,          file.path(out_res_dir, paste0(cl_tag, "_", assay_name, "_FeatureCorrelation_vs_log10VL.csv")))
  
  # plot
  plot_corr_volcano(
    res,
    title = paste0(cl, ": ", assay_name, " features vs log10(VL)"),
    out_png = file.path(out_plot_dir, paste0(cl_tag, "_", assay_name, "_CorrVolcano_vs_log10VL.png")),
    label_top = 12,
    alpha = 0.05
  )
  
  # return summary-style table for merging
  res %>%
    mutate(cluster = cl, assay = assay_name) %>%
    select(cluster, assay, feature, rho, p, n, neglog10p)
}


### Run 
clusters_use <- sort(unique(seu_pre2$Manual_Annotation))
min_cells_per_sample <- 30
min_n_samples <- 5
all_cluster_hits <- list()

for (cl in clusters_use) {
  
  # RNA
  rna_res <- run_cluster_within_sample(
    seu = seu_pre2,
    cluster_col = "Manual_Annotation",
    cl = cl,
    vl_df_use = vl_df_use,
    assay_name = "RNA",
    slot_name = "data",
    min_cells_per_sample = min_cells_per_sample,
    min_n_samples = min_n_samples,
    out_avg_dir  = dir_rna_avg,
    out_res_dir  = dir_rna_res,
    out_plot_dir = dir_rna_plots
  )
  if (is.data.frame(rna_res)) all_cluster_hits[[paste0(cl, "_RNA")]] <- rna_res
  
  # ADT (no if; safely skip if ADT not available)
  adt_res <- tryCatch(
    run_cluster_within_sample(
      seu = seu_pre2,
      cluster_col = "Manual_Annotation",
      cl = cl,
      vl_df_use = vl_df_use,
      assay_name = "ADT",
      slot_name = "data",
      min_cells_per_sample = min_cells_per_sample,
      min_n_samples = min_n_samples,
      out_avg_dir  = dir_adt_avg,
      out_res_dir  = dir_adt_res,
      out_plot_dir = dir_adt_plots
    ),
    error = function(e) NULL
  )
  if (is.data.frame(adt_res)) all_cluster_hits[[paste0(cl, "_ADT")]] <- adt_res
}

# Combine and save global summaries
dir.create(dir_clst_res, showWarnings = FALSE, recursive = TRUE)

if (length(all_cluster_hits) > 0) {
  summary_hits <- dplyr::bind_rows(all_cluster_hits) %>% arrange(p)
  readr::write_csv(summary_hits, file.path(dir_clst_res, "ClusterWithinSample_AllFeatures_vs_log10VL_Summary.csv"))
  
  top_hits <- summary_hits %>%
    group_by(cluster, assay) %>%
    slice_min(order_by = p, n = 20, with_ties = FALSE) %>%
    ungroup()
  
  readr::write_csv(top_hits, file.path(dir_clst_res, "ClusterWithinSample_Top20Hits_perClusterAssay.csv"))
}
levels(as.factor(TARA_ALL$Manual_Annotation))

###### Feature Selection VL corelations ################

summary_csv <- file.path(dir_clst_res, "ClusterWithinSample_AllFeatures_vs_log10VL_Summary.csv")
summary_hits <- readr::read_csv(summary_csv, show_col_types = FALSE)

# -----------------------------
# Parameters
# -----------------------------
min_n <- 8
top_per_cluster_rna <- 5
top_per_cluster_adt <- 3
max_total_features <- 40

# -----------------------------
# Exclude DN T cell clusters
# -----------------------------
summary_hits_noDN <- summary_hits %>%
  filter(!grepl("DN T cell", cluster))

# -----------------------------
# Pure math ranking
# -----------------------------
hits_clean <- summary_hits_noDN %>%
  filter(is.finite(rho), is.finite(n)) %>%
  filter(n >= min_n) %>%
  mutate(abs_rho = abs(rho)) %>%
  arrange(desc(abs_rho))

# Save full ranked table
write_csv(hits_clean,
          file.path(dir_heatmap, "MathRanked_NoDN_AllFeatures_by_absRho.csv"))

# -----------------------------
# Top per cluster
# -----------------------------
top_rna <- hits_clean %>%
  filter(assay == "RNA") %>%
  group_by(cluster) %>%
  slice_max(order_by = abs_rho,
            n = top_per_cluster_rna,
            with_ties = FALSE) %>%
  ungroup()

top_adt <- hits_clean %>%
  filter(assay == "ADT") %>%
  group_by(cluster) %>%
  slice_max(order_by = abs_rho,
            n = top_per_cluster_adt,
            with_ties = FALSE) %>%
  ungroup()

shortlist <- bind_rows(top_rna, top_adt) %>%
  arrange(desc(abs_rho)) %>%
  slice_head(n = max_total_features)

write_csv(shortlist,
          file.path(dir_heatmap, "MathRanked_NoDN_Shortlist_forHeatmap.csv"))

shortlist

feats <- unique(shortlist$feature)
clust <- unique(shortlist$cluster)

mat <- matrix(NA_real_,
              nrow = length(feats),
              ncol = length(clust),
              dimnames = list(feats, clust))

for (i in seq_len(nrow(shortlist))) {
  mat[shortlist$feature[i], shortlist$cluster[i]] <- shortlist$rho[i]
}

# fill missing with 0 so clustering works
mat2 <- mat
mat2[is.na(mat2)] <- 0

png(file.path(dir_heatmap, "Ranked_VL_Associated_Features_Heatmap.png"),
    width = 6000, height = 9000, res = 300)

pheatmap::pheatmap(
  mat2,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  # --- Make tiles big so fonts + borders are visible
  cellwidth  = 55,
  cellheight = 30,
  
  # --- Fonts
  fontsize_row = 12,
  fontsize_col = 14,
  angle_col = 90,
  
  # --- Grid boxes (dark + visible)
  border_color = "black",
  
  # --- Dendrogram space
  treeheight_row = 60,
  treeheight_col = 60,
  
  main = "Top Cluster-Level Associations with log10(Viral Load) (Spearman ρ)"
)

dev.off()

