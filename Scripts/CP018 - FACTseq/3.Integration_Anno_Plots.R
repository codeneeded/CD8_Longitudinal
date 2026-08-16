############################################################
# CP018 FACT-Seq — integration (CCA + FastMNN) and annotation plots
#
# Adapted from "CP003 - Longitudinal HIV Stim/2.Integration_Anno_Plots.R".
# Differences: batch = 4 libraries (2 timepoints x 2 replicates); grouping by
# Replicate / Timepoint_Rep instead of CellType_Sort.
#
# Input : QS_SINGLETS
# Output: QS_INTEGRATED + INT_ROOT / ANNO_ROOT plots
############################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)
library(SeuratExtend)
library(scCustomize)

mkdirs(INT_ROOT)

require_file(QS_SINGLETS, "post-QC singlets object")
seu <- qs_read(QS_SINGLETS)
DefaultAssay(seu) <- "RNA"

# Batch = library. Replicate A/B is exactly the technical effect integration
# should absorb; Timepoint is biological and must survive it.
seu$batch <- seu$Sample_Subfolder

# -----------------------------
# 1) Unintegrated
# -----------------------------
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$batch)

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca", verbose = FALSE)
seu <- FindClusters(seu, resolution = 2, cluster.name = "unintegrated_clusters", verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
               reduction.name = "umap.unintegrated", verbose = FALSE)

p_un <- wrap_plots(
  DimPlot2(seu, reduction = "umap.unintegrated", group.by = "batch"),
  DimPlot2(seu, reduction = "umap.unintegrated", group.by = "Timepoint"),
  DimPlot2(seu, reduction = "umap.unintegrated", group.by = "Replicate"),
  if ("predicted.celltype.l2" %in% colnames(seu@meta.data)) {
    DimPlot2(seu, reduction = "umap.unintegrated", group.by = "predicted.celltype.l2",
             label = TRUE, repel = TRUE, box = TRUE)
  } else {
    DimPlot2(seu, reduction = "umap.unintegrated", group.by = "unintegrated_clusters",
             label = TRUE, repel = TRUE, box = TRUE)
  },
  ncol = 2
)
ggsave(file.path(INT_ROOT, "CP018_RNA_Unintegrated_UMAPs.png"), p_un,
       width = 14, height = 10, dpi = 350, bg = "white")

# -----------------------------
# 2) Integrations
# -----------------------------
seu <- IntegrateLayers(object = seu, method = CCAIntegration,
                       orig.reduction = "pca", assay = "RNA",
                       new.reduction = "integrated.cca.rna", verbose = FALSE)

seu <- IntegrateLayers(object = seu, method = FastMNNIntegration, assay = "RNA",
                       new.reduction = "integrated.mnn.rna", verbose = FALSE)

# -----------------------------
# 3) Clustering grids + clustree
# -----------------------------
res_grid <- seq(0.2, 3.0, by = 0.2)

seu <- FindNeighbors(seu, reduction = "integrated.cca.rna", dims = 1:30,
                     graph.name = "cca_snn", verbose = FALSE)
for (r in res_grid) {
  seu <- FindClusters(seu, graph.name = "cca_snn", resolution = r,
                      algorithm = 1, verbose = FALSE)
}
ggsave(file.path(INT_ROOT, "CP018_clustree_CCA.png"),
       clustree(seu, prefix = "cca_snn_res."),
       width = 14, height = 12, dpi = 350, bg = "white")

seu <- FindNeighbors(seu, reduction = "integrated.mnn.rna", dims = 1:30,
                     graph.name = "mnn_snn", verbose = FALSE)
for (r in res_grid) {
  seu <- FindClusters(seu, graph.name = "mnn_snn", resolution = r,
                      algorithm = 1, verbose = FALSE)
}
ggsave(file.path(INT_ROOT, "CP018_clustree_MNN.png"),
       clustree(seu, prefix = "mnn_snn_res."),
       width = 14, height = 12, dpi = 350, bg = "white")

# Finalise resolutions (values live in config.R — revisit after clustree)
cca_col <- paste0("cca_snn_res.", CCA_RES_FINAL)
mnn_col <- paste0("mnn_snn_res.", MNN_RES_FINAL)
if (!cca_col %in% colnames(seu@meta.data)) stop("Missing CCA resolution column: ", cca_col)
if (!mnn_col %in% colnames(seu@meta.data)) stop("Missing MNN resolution column: ", mnn_col)

num_levels <- function(x) as.character(sort(as.numeric(unique(as.character(x)))))
seu$cca_clusters_rna <- factor(seu[[cca_col]][, 1], levels = num_levels(seu[[cca_col]][, 1]))
seu$mnn_clusters_rna <- factor(seu[[mnn_col]][, 1], levels = num_levels(seu[[mnn_col]][, 1]))

# Drop singleton clusters (<10 cells). At res 0.4 CP018 has one 2-cell cluster;
# it cannot support a marker test and would appear as a spurious annotation.
.tiny <- names(which(table(seu$mnn_clusters_rna) < MIN_CLUSTER_SIZE))
if (length(.tiny)) {
  message("Dropping tiny cluster(s): ", paste(.tiny, collapse = ", "),
          " (", sum(seu$mnn_clusters_rna %in% .tiny), " cells)")
  seu <- seu[, !seu$mnn_clusters_rna %in% .tiny]
  seu$mnn_clusters_rna <- droplevels(seu$mnn_clusters_rna)
}

# ── Apply the annotation: labels replace cluster numbers everywhere downstream ─
# Cluster NUMBERS are not portable (they change with resolution and between
# datasets); labels are. Fig_Annotation becomes the grouping variable for all
# plots from here on, and scripts 4-8 use it rather than mnn_clusters_rna.
if (!is.null(CP018_CLUSTER_ANNOTATION)) {
  present <- levels(seu$mnn_clusters_rna)
  missing <- setdiff(present, names(CP018_CLUSTER_ANNOTATION))
  if (length(missing)) {
    szs <- table(seu$mnn_clusters_rna)[missing]
    stop("Clusters in the object have no entry in CP018_CLUSTER_ANNOTATION.\n",
         "  unmapped cluster(s): ",
         paste(sprintf("%s (n=%d)", missing, as.integer(szs)), collapse = ", "), "\n",
         "  object has ", length(present), " clusters at ", mnn_col, "\n",
         "  annotation covers ", length(CP018_CLUSTER_ANNOTATION), ": ",
         paste(names(CP018_CLUSTER_ANNOTATION), collapse = ","), "\n",
         "  -> clustering produced a different number of clusters than the ",
         "annotation was built on.\n",
         "     Either raise the singleton cutoff, or re-derive the annotation ",
         "with 3b.Annotation_Evidence.R.", call. = FALSE)
  }
  # unname(): a named-vector lookup carries the LOOKUP's names, and Seurat then
  # tries to match those against cell barcodes and errors on no overlap.
  lab <- unname(CP018_CLUSTER_ANNOTATION[as.character(seu$mnn_clusters_rna)])
  bad <- setdiff(unique(lab), ANNO_ORDER)
  if (length(bad)) stop("Label(s) not in ANNO_ORDER: ", paste(bad, collapse = ", "),
                        call. = FALSE)
  stopifnot(length(lab) == ncol(seu), !any(is.na(lab)))
  fig_anno <- factor(lab, levels = intersect(ANNO_ORDER, unique(lab)))
  names(fig_anno) <- colnames(seu)
  seu$Fig_Annotation <- fig_anno
  Idents(seu) <- "Fig_Annotation"
  message("Annotation applied: ", paste(levels(seu$Fig_Annotation), collapse = " | "))
  print(table(seu$Fig_Annotation))
} else {
  warning("CP018_CLUSTER_ANNOTATION is NULL - plots will use cluster numbers.")
  seu$Fig_Annotation <- seu$mnn_clusters_rna
}

seu <- RunUMAP(seu, reduction = "integrated.cca.rna", dims = 1:30,
               reduction.name = "umap.cca.rna", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.mnn.rna", dims = 1:30,
               reduction.name = "umap.mnn.rna", verbose = FALSE)

# -----------------------------
# 4) Integrated embeddings
# -----------------------------
azi_or <- function(red, fallback) {
  if ("predicted.celltype.l2" %in% colnames(seu@meta.data)) {
    DimPlot2(seu, reduction = red, group.by = "predicted.celltype.l2",
             label = TRUE, repel = TRUE, box = TRUE)
  } else {
    DimPlot2(seu, reduction = red, group.by = fallback,
             label = TRUE, repel = TRUE, box = TRUE)
  }
}

ggsave(file.path(INT_ROOT, "CP018_RNA_Integration_CCA_MNN.png"),
       wrap_plots(
         azi_or("umap.cca.rna", "cca_clusters_rna"),
         DimPlot2(seu, reduction = "umap.cca.rna", group.by = "cca_clusters_rna",
                  label = TRUE, repel = TRUE, box = TRUE),
         azi_or("umap.mnn.rna", "mnn_clusters_rna"),
         DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "mnn_clusters_rna",
                  label = TRUE, repel = TRUE, box = TRUE),
         ncol = 2, byrow = TRUE),
       width = 18, height = 12, dpi = 350, bg = "white")

# Batch-effect checks. Replicate A/B SHOULD overlap after integration; if they
# separate, integration has not removed the technical effect.
for (nm in c("cca", "mnn")) {
  red <- paste0("umap.", nm, ".rna")
  ggsave(file.path(INT_ROOT, paste0("CP018_", toupper(nm), "_UMAP_checks.png")),
         wrap_plots(
           DimPlot2(seu, reduction = red, group.by = "batch"),
           DimPlot2(seu, reduction = red, group.by = "Timepoint"),
           DimPlot2(seu, reduction = red, group.by = "Replicate"),
           DimPlot2(seu, reduction = red, group.by = "Phase"),
           ncol = 4),
         width = 24, height = 6, dpi = 350, bg = "white")
}

ggsave(file.path(INT_ROOT, "CP018_MNN_UMAP_split_by_Timepoint.png"),
       DimPlot2(seu, reduction = "umap.mnn.rna", split.by = "Timepoint"),
       width = 11, height = 5.5, dpi = 350, bg = "white")

ggsave(file.path(INT_ROOT, "CP018_MNN_UMAP_split_by_Replicate.png"),
       DimPlot2(seu, reduction = "umap.mnn.rna", split.by = "Timepoint_Rep"),
       width = 20, height = 5.5, dpi = 350, bg = "white")

# -----------------------------
# CSV exports (source data for every plot above)
# -----------------------------
csv_dir <- file.path(INT_ROOT, "tables"); mkdirs(csv_dir)

# cluster sizes and composition
save_csv(as.data.frame(table(cluster = seu$mnn_clusters_rna)),
         file.path(csv_dir, "CP018_cluster_sizes.csv"))
for (v in c("Sample_Subfolder", "Timepoint", "Replicate", "Phase")) {
  if (v %in% colnames(seu@meta.data))
    save_csv(composition_table(seu, "mnn_clusters_rna", v),
             file.path(csv_dir, paste0("CP018_composition_cluster_by_", v, ".csv")))
}

# Azimuth label x cluster
if ("predicted.celltype.l2" %in% colnames(seu@meta.data)) {
  save_csv(as.data.frame.matrix(table(seu$mnn_clusters_rna, seu$predicted.celltype.l2)),
           file.path(csv_dir, "CP018_cluster_by_Azimuth_l2.csv"), row.names = TRUE)
}

# differential markers for the chosen resolution
Idents(seu) <- "mnn_clusters_rna"
markers_all <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25, verbose = FALSE)
save_csv(markers_all, file.path(csv_dir, "CP018_markers_all_clusters.csv"))
save_csv(markers_all %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 25) %>% ungroup(),
         file.path(csv_dir, "CP018_markers_top25_per_cluster.csv"))

# UMAP coordinates + per-cell metadata (lets anyone re-plot without Seurat)
umap_df <- as.data.frame(Embeddings(seu, "umap.mnn.rna"))
names(umap_df) <- c("UMAP1", "UMAP2")
umap_df$cell_barcode <- rownames(umap_df)
umap_df$cluster      <- as.character(seu$mnn_clusters_rna)
umap_df$Timepoint    <- seu$Timepoint
umap_df$Replicate    <- seu$Replicate
save_csv(umap_df, file.path(csv_dir, "CP018_UMAP_coordinates.csv"))
export_metadata(seu, file.path(csv_dir, "CP018_cell_metadata.csv"))

# -----------------------------
# Annotated plots (labels, not cluster numbers)
# -----------------------------
anno_dir <- file.path(ANNO_ROOT, "Annotated"); mkdirs(anno_dir)
pal <- ANNO_COLS[levels(seu$Fig_Annotation)]

ggsave(file.path(anno_dir, "CP018_UMAP_annotated.png"),
       DimPlot(seu, reduction = "umap.mnn.rna", group.by = "Fig_Annotation",
               cols = pal, label = TRUE, repel = TRUE, label.size = 4.5) +
         ggtitle("CP018 annotated clusters") + theme_cp018(),
       width = 9, height = 7, dpi = 350, bg = "white")

ggsave(file.path(anno_dir, "CP018_UMAP_annotated_by_Timepoint.png"),
       DimPlot(seu, reduction = "umap.mnn.rna", group.by = "Fig_Annotation",
               split.by = "Timepoint", cols = pal) +
         ggtitle("CP018 annotated clusters by timepoint") + theme_cp018(),
       width = 14, height = 6.5, dpi = 350, bg = "white")

ggsave(file.path(anno_dir, "CP018_UMAP_annotated_by_Replicate.png"),
       DimPlot(seu, reduction = "umap.mnn.rna", group.by = "Fig_Annotation",
               split.by = "Timepoint_Rep", cols = pal, ncol = 2) +
         ggtitle("CP018 annotated clusters by library") + theme_cp018(),
       width = 13, height = 11, dpi = 350, bg = "white")

# composition: stacked bars + table
comp <- composition_table(seu, "Fig_Annotation", "Sample_Subfolder")
save_csv(comp, file.path(csv_dir, "CP018_composition_annotation_by_library.csv"))
comp_tp <- composition_table(seu, "Fig_Annotation", "Timepoint")
save_csv(comp_tp, file.path(csv_dir, "CP018_composition_annotation_by_timepoint.csv"))

ggsave(file.path(anno_dir, "CP018_composition_by_library.png"),
       ggplot(comp, aes(x = Sample_Subfolder, y = pct_of_sample, fill = Fig_Annotation)) +
         geom_col(colour = "black", linewidth = 0.25, width = 0.7) +
         scale_fill_manual(values = pal, name = NULL) +
         labs(x = NULL, y = "% of library", title = "Cell-type composition per library") +
         theme_cp018() +
         theme(axis.text.x = element_text(angle = 30, hjust = 1)),
       width = 9, height = 6, dpi = 350, bg = "white")

ggsave(file.path(anno_dir, "CP018_composition_by_timepoint.png"),
       ggplot(comp_tp, aes(x = Timepoint, y = pct_of_sample, fill = Fig_Annotation)) +
         geom_col(colour = "black", linewidth = 0.25, width = 0.6) +
         geom_text(aes(label = ifelse(pct_of_sample >= 3,
                                      sprintf("%.0f%%", pct_of_sample), "")),
                   position = position_stack(vjust = 0.5), size = 4) +
         scale_fill_manual(values = pal, name = NULL) +
         labs(x = NULL, y = "% of timepoint", title = "Composition: 2m vs 90m") +
         theme_cp018(),
       width = 7.5, height = 6.5, dpi = 350, bg = "white")

# marker dot plot by LABEL
mk_feats <- intersect(c("SELL","CCR7","IL7R","TCF7","LEF1","ID3","CD27","CD38",
                        "GZMK","GZMB","GNLY","PRF1","FGFBP2","NKG7",
                        "MKI67","TOP2A","TOX","PDCD1","LAG3","HAVCR2",
                        "KLRB1","SLC4A10","TRAV1-2","BNIP3","P4HA1","MIR210HG"),
                      rownames(seu))
ggsave(file.path(anno_dir, "CP018_marker_dotplot_annotated.png"),
       DotPlot(seu, features = mk_feats, group.by = "Fig_Annotation") +
         scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B") +
         labs(x = NULL, y = NULL, title = "Marker expression by annotation") +
         theme_cp018(base_size = 13) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       width = 13, height = 6, dpi = 350, bg = "white")

save_csv(as.data.frame(table(Annotation = seu$Fig_Annotation)),
         file.path(csv_dir, "CP018_annotation_sizes.csv"))

qs_save(seu, file = QS_INTEGRATED)
message("Saved: ", QS_INTEGRATED)

############################################################
# Annotation helper plots (MNN clusters)
############################################################

dir_vln  <- file.path(ANNO_ROOT, "VlnPlot2_RNA")
dir_ftr  <- file.path(ANNO_ROOT, "FeaturePlot2_RNA")
dir_dist <- file.path(ANNO_ROOT, "Cluster_Distribution")
mkdirs(dir_vln, dir_ftr, dir_dist)

Idents(seu) <- "mnn_clusters_rna"

# Marker panel: CP003 FACT-Seq panel plus the module genes used to annotate the
# six functional populations in Figure 3 (naive/stemness, exhaustion,
# cytotoxicity, Tpex, cycling) and the Table 2 IFN-memory module.
rna.features <- c(
  # naive / stemness
  'TCF7','LEF1','SELL','CCR7','IL7R','BACH2','BCL2','S1PR1','KLF2',
  # exhaustion
  'TOX','PDCD1','HAVCR2','LAG3','TIGIT','CTLA4','ENTPD1',
  # cytotoxicity
  'GZMB','GNLY','PRF1','NKG7','GZMA','GZMH','GZMK','FGFBP2',
  # terminal differentiation
  'ZEB2','PRDM1','TBX21','CX3CR1','S1PR5','ID2','KLRG1',
  # type I IFN memory (Table 2)
  'IFIT1','IFIT3','ISG15','MX1',
  # acute stress / IFN
  'HSPA1A','HSPA1B','DNAJB1','IFI27','IFI44L',
  # inflammatory chemokines
  'CCL3','CCL4','CCL5','CCL4L2',
  # proliferation
  'MKI67','TOP2A','PCNA','MCM6','HIST1H4C',
  # lineage / identity + contamination checks
  'CD3D','CD3E','CD8A','CD8B','CD4','TRAC','TRBC1','TRDC','TRGC1','KLRB1',
  'FCGR3A','NCAM1','KLRD1','KLRC1','CD14','CD19','MS4A1','FOXP3',
  # activation
  'FAS','HLA-DRA','CD27','CD28','TNFRSF9','ICOS','IFNG','IL2RA'
)

rna.features.use <- intersect(rna.features, rownames(seu))
missing_genes <- setdiff(rna.features, rna.features.use)
if (length(missing_genes) > 0)
  message("Not in object: ", paste(missing_genes, collapse = ", "))

for (g in rna.features.use) {
  ggsave(file.path(dir_vln, paste0(g, "_VLNplot.png")),
         VlnPlot2(seu, features = g, cols = "default", show.mean = TRUE,
                  mean_colors = c("red", "blue")) + ggtitle(paste("RNA |", g)),
         dpi = 500, width = 14, height = 8, bg = "white")

  ggsave(file.path(dir_ftr, paste0(g, "_FeaturePlot.png")),
         FeaturePlot_scCustom(seu, features = g, reduction = "umap.mnn.rna") +
           ggtitle(paste("RNA |", g)),
         dpi = 500, width = 10, height = 8, bg = "white")
}

for (pc in c(FALSE, TRUE)) {
  tag <- if (pc) "percent" else "counts"
  ggsave(file.path(dir_dist, paste0("CP018_ClusterDistrBar_", tag, ".png")),
         ClusterDistrBar(origin = seu$Sample_Subfolder, cluster = seu$mnn_clusters_rna,
                         cols = "default", flip = FALSE, border = "black",
                         percent = pc) +
           theme(axis.title.x = element_blank()) +
           ggtitle(paste0("Cluster distribution | mnn_clusters_rna (", tag, ")")),
         dpi = 400, width = 16, height = 8, bg = "white")
}

# Cluster x replicate table — a cluster driven by one replicate is suspect
clust_tab <- as.data.frame.matrix(table(seu$mnn_clusters_rna, seu$Timepoint_Rep))
clust_tab$Cluster <- rownames(clust_tab)
write.csv(clust_tab, file.path(dir_dist, "CP018_cluster_by_library.csv"), row.names = FALSE)
print(clust_tab)
