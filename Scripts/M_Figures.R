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