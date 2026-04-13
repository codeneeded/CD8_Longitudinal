############################################################
# Trex epitope analysis (CMV vs HIV) — TARA CD/NK object
# - Uses Trex::annotateDB
# - Keeps ALL cells for UMAP context (wnn.umap)
# - Uses SeuratExtend::DimPlot2, SeuratExtend::CalcStats + SeuratExtend::Heatmap
# - DE via MAST + volcano
# - Summaries by orig.ident and Timepoint_Group
# - Step-wise output folders
############################################################

# -----------------------------
# Libraries (no suppressMessages)
# -----------------------------
library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(Trex)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "C:/Users/ammas/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "tara_cdnk_annotated.rds")

tara_cdnk <- readRDS(rds_in)
DefaultAssay(tara_cdnk) <- "RNA"

# -----------------------------
# Output folders
# -----------------------------
out_root <- file.path(base_dir, "Trex_Epitope_CMV_HIV")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

step01_dir <- file.path(out_root, "Step_01_Annotate_and_Flags")
step02_dir <- file.path(out_root, "Step_02_UMAPs")
step03_dir <- file.path(out_root, "Step_03_Composition_Summaries")
step04_dir <- file.path(out_root, "Step_04_DE_MAST_and_Volcano")
step05_dir <- file.path(out_root, "Step_05_Heatmaps_SeuratExtend")

dir.create(step01_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step02_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step03_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step04_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step05_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# User settings
# -----------------------------
reduction_use <- "wnn.umap"

# Choose your main cluster label here:
# - "seurat_clusters" OR your annotation column (e.g., "Manual_Annotation", "T_Cell_A", etc.)
cluster_col <- "seurat_clusters"

# Clone expansion definition
exp_levels <- c(
  "Medium (5 < X <= 20)",
  "Large (20 < X <= 100)",
  "Hyperexpanded (100 < X <= 500)"
)

# -----------------------------
# Helper: safe ggplot saving
# -----------------------------
save_plot <- function(p, out_png, w = 7, h = 6, dpi = 300) {
  ggsave(
    filename = out_png,
    plot = p,
    width = w,
    height = h,
    dpi = dpi,
    bg = "white"   # ← THIS fixes transparent background
  )
}
# Helper: save SeuratExtend heatmap (not always a ggplot)
save_seuratextend_heatmap <- function(hm_obj, out_png, w_px = 1400, h_px = 900, res = 150) {
  png(out_png, width = w_px, height = h_px, res = res)
  print(hm_obj)
  dev.off()
}

# -----------------------------
# STEP 01: Trex annotate + flags
# -----------------------------
tara_cdnk_TRB_2 <- annotateDB(
  tara_cdnk,
  chains = "TRB",
  edit.distance = 2
)

# Save raw table of predicted species labels
tab_species <- sort(table(tara_cdnk_TRB_2$TRB_Epitope.species), decreasing = TRUE)
write.csv(
  data.frame(TRB_Epitope.species = names(tab_species), n = as.integer(tab_species)),
  file.path(step01_dir, "Table_TRB_Epitope_species_counts.csv"),
  row.names = FALSE
)

# Epitope grouping (multi-label strings)
tara_cdnk_TRB_2$Ep_CMV <- grepl("CMV", tara_cdnk_TRB_2$TRB_Epitope.species %||% "")
tara_cdnk_TRB_2$Ep_HIV <- grepl("HIV", tara_cdnk_TRB_2$TRB_Epitope.species %||% "")

tara_cdnk_TRB_2$Epitope_Group <- "Other"
tara_cdnk_TRB_2$Epitope_Group[tara_cdnk_TRB_2$Ep_CMV & !tara_cdnk_TRB_2$Ep_HIV] <- "CMV"
tara_cdnk_TRB_2$Epitope_Group[tara_cdnk_TRB_2$Ep_HIV & !tara_cdnk_TRB_2$Ep_CMV] <- "HIV"
tara_cdnk_TRB_2$Epitope_Group[tara_cdnk_TRB_2$Ep_HIV &  tara_cdnk_TRB_2$Ep_CMV] <- "CMV_HIV"

tara_cdnk_TRB_2$Epitope_Group <- factor(
  tara_cdnk_TRB_2$Epitope_Group,
  levels = c("HIV","CMV_HIV","CMV","Other")
)

# Expansion flag
tara_cdnk_TRB_2$Expanding <- tara_cdnk_TRB_2$cloneSize %in% exp_levels

# Convenience combined label for plotting
tara_cdnk_TRB_2$Epitope_Exp <- "Other"
tara_cdnk_TRB_2$Epitope_Exp[tara_cdnk_TRB_2$Expanding & tara_cdnk_TRB_2$Epitope_Group == "CMV"] <- "CMV_exp"
tara_cdnk_TRB_2$Epitope_Exp[tara_cdnk_TRB_2$Expanding & tara_cdnk_TRB_2$Epitope_Group == "HIV"] <- "HIV_exp"
tara_cdnk_TRB_2$Epitope_Exp[tara_cdnk_TRB_2$Expanding & tara_cdnk_TRB_2$Epitope_Group == "CMV_HIV"] <- "CMV_HIV_exp"
tara_cdnk_TRB_2$Epitope_Exp <- factor(
  tara_cdnk_TRB_2$Epitope_Exp,
  levels = c("HIV_exp","CMV_HIV_exp","CMV_exp","Other")
)

# Save quick counts
write.csv(
  as.data.frame(table(Epitope_Group = tara_cdnk_TRB_2$Epitope_Group)),
  file.path(step01_dir, "Table_Epitope_Group_counts_ALL.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(Expanding = tara_cdnk_TRB_2$Expanding)),
  file.path(step01_dir, "Table_Expanding_counts_ALL.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(Epitope_Exp = tara_cdnk_TRB_2$Epitope_Exp)),
  file.path(step01_dir, "Table_EpitopeExp_counts_ALL.csv"),
  row.names = FALSE
)

# -----------------------------
# STEP 02: UMAP plots (ALL cells)
# -----------------------------
p_umap_expand <- DimPlot2(
  tara_cdnk_TRB_2,
  group.by  = "Expanding",
  reduction = reduction_use
)
save_plot(p_umap_expand, file.path(step02_dir, "UMAP_Expanding.png"), w = 7, h = 6)

p_umap_epitopeexp <- DimPlot2(
  tara_cdnk_TRB_2,
  group.by  = "Epitope_Exp",
  reduction = reduction_use
)
save_plot(p_umap_epitopeexp, file.path(step02_dir, "UMAP_EpitopeExp.png"), w = 7, h = 6)

p_umap_highlight_any_epitope <- DimPlot2(
  tara_cdnk_TRB_2,
  reduction = reduction_use,
  cells.highlight = WhichCells(
    tara_cdnk_TRB_2,
    expression = Epitope_Group %in% c("CMV","HIV","CMV_HIV")
  ),
  cols = "grey85"
)
save_plot(p_umap_highlight_any_epitope, file.path(step02_dir, "UMAP_Highlight_EpitopeCells.png"), w = 7, h = 6)

# Split UMAP by Timepoint_Group (can be wide; increase width)
p_umap_split_tp <- DimPlot2(
  tara_cdnk_TRB_2,
  group.by  = "Epitope_Exp",
  split.by  = "Timepoint_Group",
  reduction = reduction_use
)
save_plot(p_umap_split_tp, file.path(step02_dir, "UMAP_EpitopeExp_splitBy_TimepointGroup.png"), w = 14, h = 6)

# -----------------------------
# STEP 03: Composition summaries
# - Focus on expanding epitope groups (CMV/HIV/CMV_HIV) but retain ALL context in object
# -----------------------------
meta <- tara_cdnk_TRB_2@meta.data %>%
  mutate(
    Epitope_Group = as.character(Epitope_Group),
    Epitope_Exp   = as.character(Epitope_Exp),
    Expanding     = as.logical(Expanding),
    cluster_label = .data[[cluster_col]]
  )

meta_exp_epi <- meta %>%
  filter(Expanding, Epitope_Group %in% c("CMV","HIV","CMV_HIV"))

# A) orig.ident breakdown (fraction within each orig.ident)
tab_orig <- meta_exp_epi %>%
  count(orig.ident, Epitope_Group, name = "n") %>%
  group_by(orig.ident) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
write.csv(tab_orig, file.path(step03_dir, "Expanded_Epitopes_byOrigIdent_counts_and_fraction.csv"), row.names = FALSE)

p_orig_frac <- ggplot(tab_orig, aes(x = orig.ident, y = frac, fill = Epitope_Group)) +
  geom_col() +
  coord_flip()
save_plot(p_orig_frac, file.path(step03_dir, "BarFrac_ExpandedEpitopes_byOrigIdent.png"), w = 10, h = 9)

# B) Timepoint_Group breakdown (fraction within each Timepoint_Group)
tab_tp <- meta_exp_epi %>%
  count(Timepoint_Group, Epitope_Group, name = "n") %>%
  group_by(Timepoint_Group) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
write.csv(tab_tp, file.path(step03_dir, "Expanded_Epitopes_byTimepointGroup_counts_and_fraction.csv"), row.names = FALSE)

p_tp_frac <- ggplot(tab_tp, aes(x = Timepoint_Group, y = frac, fill = Epitope_Group)) +
  geom_col()
save_plot(p_tp_frac, file.path(step03_dir, "BarFrac_ExpandedEpitopes_byTimepointGroup.png"), w = 8, h = 5)

# C) Cluster composition (fraction within cluster) for expanding epitope groups
tab_cluster <- meta_exp_epi %>%
  count(cluster_label, Epitope_Group, name = "n") %>%
  group_by(cluster_label) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

write.csv(tab_cluster, file.path(step03_dir, paste0("Expanded_Epitopes_by_", cluster_col, "_counts_and_fraction.csv")), row.names = FALSE)

p_cluster_frac <- ggplot(tab_cluster, aes(x = cluster_label, y = frac, fill = Epitope_Group)) +
  geom_col() +
  coord_flip()
save_plot(p_cluster_frac, file.path(step03_dir, paste0("BarFrac_ExpandedEpitopes_by_", cluster_col, ".png")), w = 9, h = 10)

# -----------------------------
# STEP 04: Differential expression (MAST) + Volcano
# - Expanded only
# - CMV vs HIV (exclude CMV_HIV)
# - Keep the big object; subset only inside analysis
# -----------------------------
cells_use <- WhichCells(
  tara_cdnk_TRB_2,
  expression = Expanding == TRUE & Epitope_Group %in% c("CMV","HIV")
)

obj_use <- subset(tara_cdnk_TRB_2, cells = cells_use)
Idents(obj_use) <- "Epitope_Group"

markers_cmv_vs_hiv <- FindMarkers(
  obj_use,
  ident.1 = "CMV",
  ident.2 = "HIV",
  test.use = "MAST",
  latent.vars = "nCount_RNA",
  logfc.threshold = 0.0,
  min.pct = 0.05
)

markers_cmv_vs_hiv$gene <- rownames(markers_cmv_vs_hiv)
write.csv(markers_cmv_vs_hiv, file.path(step04_dir, "MAST_Expanded_CMV_vs_HIV.csv"), row.names = TRUE)

# Volcano
pval_col <- if ("p_val_adj" %in% colnames(markers_cmv_vs_hiv)) "p_val_adj" else "p.adjust"
lfc_col  <- if ("avg_log2FC" %in% colnames(markers_cmv_vs_hiv)) "avg_log2FC" else "avg_logFC"

vol <- EnhancedVolcano(
  markers_cmv_vs_hiv,
  lab = markers_cmv_vs_hiv$gene,
  x = lfc_col,
  y = pval_col,
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 1.8,
  labSize = 3.0,
  title = "Expanded clones: CMV vs HIV (MAST)",
  subtitle = "Idents = Epitope_Group; latent.vars = nCount_RNA"
)
save_plot(vol, file.path(step04_dir, "Volcano_Expanded_CMV_vs_HIV_MAST.png"), w = 7.5, h = 6.5)

# Optional: CMV vs HIV within each Timepoint_Group (skips if too few cells)
tp_levels <- intersect(
  c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"),
  unique(as.character(tara_cdnk_TRB_2$Timepoint_Group))
)

dge_by_tp <- list()

for (tp in tp_levels) {
  cells_tp <- WhichCells(
    tara_cdnk_TRB_2,
    expression = Timepoint_Group == tp & Expanding == TRUE & Epitope_Group %in% c("CMV","HIV")
  )
  
  if (length(cells_tp) < 50) next
  
  obj_tp <- subset(tara_cdnk_TRB_2, cells = cells_tp)
  Idents(obj_tp) <- "Epitope_Group"
  
  if (!all(c("CMV","HIV") %in% levels(Idents(obj_tp)))) next
  if (min(table(Idents(obj_tp))[c("CMV","HIV")]) < 20) next
  
  dge <- FindMarkers(
    obj_tp,
    ident.1 = "CMV",
    ident.2 = "HIV",
    test.use = "MAST",
    latent.vars = "nCount_RNA",
    logfc.threshold = 0.0,
    min.pct = 0.05
  )
  dge$gene <- rownames(dge)
  dge_by_tp[[tp]] <- dge
  
  write.csv(dge, file.path(step04_dir, paste0("MAST_Expanded_CMV_vs_HIV__", tp, ".csv")), row.names = TRUE)
  
  pval_col_tp <- if ("p_val_adj" %in% colnames(dge)) "p_val_adj" else "p.adjust"
  lfc_col_tp  <- if ("avg_log2FC" %in% colnames(dge)) "avg_log2FC" else "avg_logFC"
  
  vol_tp <- EnhancedVolcano(
    dge,
    lab = dge$gene,
    x = lfc_col_tp,
    y = pval_col_tp,
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 1.8,
    labSize = 3.0,
    title = paste0("Expanded CMV vs HIV (MAST) — ", tp)
  )
  save_plot(vol_tp, file.path(step04_dir, paste0("Volcano_Expanded_CMV_vs_HIV__", tp, ".png")), w = 7.5, h = 6.5)
}

# -----------------------------
# STEP 05: Heatmaps (SeuratExtend CalcStats + Heatmap)
# You asked for auto-selected genes like:
# genes <- VariableFeatures(pbmc)
# toplot <- CalcStats(pbmc, features=genes, method="zscore", order="p", n=5)
# Heatmap(toplot, lab_fill="zscore", color_scheme="BuGn")
# -----------------------------

# 5A) Heatmap on ALL cells using VariableFeatures
# Make sure VariableFeatures exist; if empty, run FindVariableFeatures quickly
if (length(VariableFeatures(tara_cdnk_TRB_2)) == 0) {
  tara_cdnk_TRB_2 <- FindVariableFeatures(tara_cdnk_TRB_2, selection.method = "vst", nfeatures = 2000)
}

genes_all <- VariableFeatures(tara_cdnk_TRB_2)

toplot_all <- CalcStats(
  tara_cdnk_TRB_2,
  features = genes_all,
  method   = "zscore",
  order    = "p",
  n        = 5
)
Heatmap(
  toplot_all,
  lab_fill     = "zscore",
  feature_text_subset = genes_all[1:20],
  color_scheme = "BuGn"
)

hm_all <- Heatmap(
  toplot_all,
  lab_fill     = "zscore",
  feature_text_subset = genes[1:20],
  color_scheme = "BuGn"
)

save_seuratextend_heatmap(hm_all, file.path(step05_dir, "Heatmap_ALLcells_VariableFeatures_zscore_top5_byP.png"))

# 5B) Heatmap for EXPANDED epitope comparison using top DE genes (autoselect from MAST)
# Select top N up/down genes from MAST results (robust to col names)
pval_col <- if ("p_val_adj" %in% colnames(markers_cmv_vs_hiv)) "p_val_adj" else "p.adjust"
lfc_col  <- if ("avg_log2FC" %in% colnames(markers_cmv_vs_hiv)) "avg_log2FC" else "avg_logFC"

# Remove NA and rank
mk <- markers_cmv_vs_hiv %>%
  filter(!is.na(.data[[pval_col]]), !is.na(.data[[lfc_col]]))

topN <- 25
top_up   <- mk %>% arrange(.data[[pval_col]], desc(.data[[lfc_col]])) %>% head(topN) %>% pull(gene)
top_down <- mk %>% arrange(.data[[pval_col]], .data[[lfc_col]])        %>% head(topN) %>% pull(gene)
genes_de <- unique(c(top_up, top_down))
genes_de <- genes_de[genes_de %in% rownames(tara_cdnk_TRB_2)]

# Build a small object focused on expanded CMV/HIV for heatmap interpretability
obj_hm <- obj_use
Idents(obj_hm) <- "Epitope_Group"

# CalcStats over the DE gene set
toplot_de <- CalcStats(
  obj_hm,
  features = genes_de,
  method   = "zscore",
  order    = "p",
  n        = 5
)

hm_de <- Heatmap(
  toplot_de,
  lab_fill     = "zscore",
  color_scheme = "BuGn"
)

save_seuratextend_heatmap(hm_de, file.path(step05_dir, "Heatmap_EXPANDED_CMVvsHIV_TopDEgenes_zscore_top5_byP.png"))

# 5C) Optional: Cluster-based heatmap (using your chosen cluster_col)
# This uses CalcStats ordering; interpret as "top features differentiating groups" within this object
if (cluster_col %in% colnames(tara_cdnk_TRB_2@meta.data)) {
  obj_cl <- tara_cdnk_TRB_2
  Idents(obj_cl) <- cluster_col
  
  # Use variable features again for cluster separation
  genes_cl <- VariableFeatures(obj_cl)
  
  toplot_cl <- CalcStats(
    obj_cl,
    features = genes_cl,
    method   = "zscore",
    order    = "p",
    n        = 5
  )
  
  hm_cl <- Heatmap(
    toplot_cl,
    lab_fill     = "zscore",
    color_scheme = "BuGn"
  )
  
  save_seuratextend_heatmap(hm_cl, file.path(step05_dir, paste0("Heatmap_", cluster_col, "_VariableFeatures_zscore_top5_byP.png")))
}

# -----------------------------
# Done
# -----------------------------
message("DONE. Outputs in: ", out_root)

