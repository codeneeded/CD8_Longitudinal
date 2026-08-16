################################################################################
# CP018 — set the final resolution and apply the annotation (NO re-integration)
#
# WHY THIS EXISTS
#   Script 3 does the expensive work: CCA + FastMNN integration, clustering at
#   many resolutions, UMAPs. That work is ALREADY DONE and saved in
#   QS_INTEGRATED - the object carries mnn_snn_res.0.2 ... mnn_snn_res.2.0 as
#   separate metadata columns. The only thing wrong with it was that
#   mnn_clusters_rna pointed at resolution 0.6 while the annotation was derived
#   at 0.4.
#
#   Re-running script 3 to fix that re-runs IntegrateLayers on a 5,131-cell
#   object, which exhausted memory on this machine (the run died at the layer
#   split with no R traceback - the OOM killer's signature).
#
#   This script does the cheap part only: point mnn_clusters_rna at the right
#   resolution, drop sub-threshold clusters, attach Fig_Annotation, regenerate
#   the annotated plots, and re-save. One object load, no integration.
#
# SAFE TO RE-RUN. Reads QS_INTEGRATED, writes QS_INTEGRATED.
#
# Input : QS_INTEGRATED (from a previous successful script 3)
# Output: QS_INTEGRATED (mnn_clusters_rna + Fig_Annotation) + annotated plots
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(dplyr)
library(ggplot2)

require_file(QS_INTEGRATED, "integrated object (script 3)")
message("Loading ", basename(QS_INTEGRATED), " ...")
seu <- qs_read(QS_INTEGRATED)
message("  ", ncol(seu), " cells x ", nrow(seu), " features")

################################################################################
# 1. Point mnn_clusters_rna at the resolution the annotation was built on
################################################################################
mnn_col <- paste0("mnn_snn_res.", MNN_RES_FINAL)
cca_col <- paste0("cca_snn_res.", CCA_RES_FINAL)

if (!mnn_col %in% colnames(seu@meta.data)) {
  stop("Object has no column '", mnn_col, "'.\n",
       "  available: ",
       paste(grep("^mnn_snn_res", colnames(seu@meta.data), value = TRUE), collapse = ", "),
       "\n  -> MNN_RES_FINAL in config.R does not match a clustering that was run.",
       call. = FALSE)
}

num_levels <- function(x) as.character(sort(as.numeric(unique(as.character(x)))))
before <- if ("mnn_clusters_rna" %in% colnames(seu@meta.data))
            nlevels(factor(seu$mnn_clusters_rna)) else NA

seu$mnn_clusters_rna <- factor(seu@meta.data[[mnn_col]],
                               levels = num_levels(seu@meta.data[[mnn_col]]))
if (cca_col %in% colnames(seu@meta.data)) {
  seu$cca_clusters_rna <- factor(seu@meta.data[[cca_col]],
                                 levels = num_levels(seu@meta.data[[cca_col]]))
}
message("mnn_clusters_rna: ", before, " clusters -> ",
        nlevels(seu$mnn_clusters_rna), " (from ", mnn_col, ")")

################################################################################
# 2. Drop sub-threshold clusters
################################################################################
.tiny <- names(which(table(seu$mnn_clusters_rna) < MIN_CLUSTER_SIZE))
if (length(.tiny)) {
  message("Dropping cluster(s) < ", MIN_CLUSTER_SIZE, " cells: ",
          paste(.tiny, collapse = ", "),
          " (", sum(seu$mnn_clusters_rna %in% .tiny), " cells)")
  seu <- seu[, !seu$mnn_clusters_rna %in% .tiny]
  seu$mnn_clusters_rna <- droplevels(seu$mnn_clusters_rna)
}
print(table(seu$mnn_clusters_rna))

################################################################################
# 3. Apply the annotation
################################################################################
if (is.null(CP018_CLUSTER_ANNOTATION))
  stop("CP018_CLUSTER_ANNOTATION is NULL in config.R.", call. = FALSE)

present <- levels(seu$mnn_clusters_rna)
missing <- setdiff(present, names(CP018_CLUSTER_ANNOTATION))
if (length(missing)) {
  szs <- table(seu$mnn_clusters_rna)[missing]
  stop("Clusters with no entry in CP018_CLUSTER_ANNOTATION.\n",
       "  unmapped: ", paste(sprintf("%s (n=%d)", missing, as.integer(szs)),
                             collapse = ", "), "\n",
       "  object has ", length(present), " clusters at ", mnn_col, "\n",
       "  annotation covers ", length(CP018_CLUSTER_ANNOTATION), ": ",
       paste(names(CP018_CLUSTER_ANNOTATION), collapse = ","), "\n",
       "  -> raise MIN_CLUSTER_SIZE, or re-derive labels with 3b.", call. = FALSE)
}
# Indexing a named vector returns names from the LOOKUP table, not the cells.
# Seurat matches metadata by cell name, so those must be stripped and replaced
# with the barcodes or the assignment fails with "No cell overlap".
lab <- unname(CP018_CLUSTER_ANNOTATION[as.character(seu$mnn_clusters_rna)])
bad <- setdiff(unique(lab), ANNO_ORDER)
if (length(bad))
  stop("Label(s) not in ANNO_ORDER: ", paste(bad, collapse = ", "), call. = FALSE)
stopifnot(length(lab) == ncol(seu), !any(is.na(lab)))

fig_anno <- factor(lab, levels = intersect(ANNO_ORDER, unique(lab)))
names(fig_anno) <- colnames(seu)
seu$Fig_Annotation <- fig_anno
Idents(seu) <- "Fig_Annotation"
message("Annotation applied: ", paste(levels(seu$Fig_Annotation), collapse = " | "))
print(table(seu$Fig_Annotation))

################################################################################
# 4. Annotated plots — labels, never cluster numbers
################################################################################
anno_dir <- file.path(ANNO_ROOT, "Annotated")
csv_dir  <- file.path(INT_ROOT, "tables")
mkdirs(anno_dir, csv_dir)

pal <- ANNO_COLS[levels(seu$Fig_Annotation)]
umap_key <- if ("umap.mnn.rna" %in% Reductions(seu)) "umap.mnn.rna" else Reductions(seu)[1]
message("UMAP reduction: ", umap_key)

# Two versions: labelled (standalone) and legend-only (for multi-panel figures,
# where a repeated legend wastes space).
p_um <- DimPlot(seu, reduction = umap_key, group.by = "Fig_Annotation",
                cols = pal, pt.size = 0.35) + ggtitle(NULL)

ggsave(file.path(anno_dir, "CP018_UMAP_annotated.png"),
       tidy_umap(umap_labels(p_um), title = "CP018 annotated clusters",
                 subtitle = paste0(ncol(seu), " cells | resolution ", MNN_RES_FINAL),
                 legend = FALSE),
       width = 8, height = 7, dpi = 350, bg = "white")

ggsave(file.path(anno_dir, "CP018_UMAP_annotated_legend.png"),
       tidy_umap(p_um, title = "CP018 annotated clusters"),
       width = 9.5, height = 7, dpi = 350, bg = "white")

ggsave(file.path(anno_dir, "CP018_UMAP_annotated_by_Timepoint.png"),
       tidy_umap(DimPlot(seu, reduction = umap_key, group.by = "Fig_Annotation",
                         split.by = "Timepoint", cols = pal, pt.size = 0.35),
                 title = "CP018 annotated clusters by timepoint"),
       width = 14, height = 6.5, dpi = 350, bg = "white")

if ("Timepoint_Rep" %in% colnames(seu@meta.data)) {
  ggsave(file.path(anno_dir, "CP018_UMAP_annotated_by_library.png"),
         tidy_umap(DimPlot(seu, reduction = umap_key, group.by = "Fig_Annotation",
                           split.by = "Timepoint_Rep", cols = pal, ncol = 2,
                           pt.size = 0.35),
                   title = "CP018 annotated clusters by library"),
         width = 13, height = 11, dpi = 350, bg = "white")
}

comp    <- composition_table(seu, "Fig_Annotation", "Sample_Subfolder")
comp_tp <- composition_table(seu, "Fig_Annotation", "Timepoint")
save_csv(comp,    file.path(csv_dir, "CP018_composition_annotation_by_library.csv"))
save_csv(comp_tp, file.path(csv_dir, "CP018_composition_annotation_by_timepoint.csv"))

ggsave(file.path(anno_dir, "CP018_composition_by_library.png"),
       ggplot(comp, aes(x = Sample_Subfolder, y = pct_of_sample, fill = Fig_Annotation)) +
         geom_col(colour = "black", linewidth = 0.25, width = 0.7) +
         scale_fill_manual(values = pal, name = NULL) +
         labs(x = NULL, y = "% of library", title = "Composition per library") +
         theme_cp018() + theme(axis.text.x = element_text(angle = 30, hjust = 1)),
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
save_csv(as.data.frame(table(cluster = seu$mnn_clusters_rna)),
         file.path(csv_dir, "CP018_cluster_sizes.csv"))
export_metadata(seu, file.path(csv_dir, "CP018_cell_metadata.csv"))

################################################################################
# 5. Save
################################################################################
qs_save(seu, file = QS_INTEGRATED)
message("\nSaved: ", QS_INTEGRATED)
cat("\nmnn_clusters_rna is now resolution ", MNN_RES_FINAL,
    " and Fig_Annotation is attached.\nScripts 3b / 4 / 4b / 5 / 8 can run.\n", sep = "")
