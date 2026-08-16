############################################################
# CP018 FACT-Seq — attach TCR (scRepertoire) + repertoire analysis
#
# Adapted from CP003 scripts 3.VDJ.R and 4.VDJ_Seurat.R.
#
# Two views of the repertoire, deliberately kept separate:
#   (a) PER-LIBRARY (4 samples) -> attached to Seurat via combineExpression,
#       because Seurat cell barcodes are prefixed per library.
#   (b) REPLICATE-MERGED (2 samples: 2m, 90m) -> for clonalOverlap and
#       diversity, which need unique sample names per biological group.
#       CP003 did the same thing when merging its 2m replicates.
#
# Input : QS_INTEGRATED, CR_OUT/<lib>/vdj_t/filtered_contig_annotations.csv
# Output: QS_WITH_TCR + VDJ_ROOT plots/tables
############################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(scRepertoire)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(SeuratExtend)

dir_clonal    <- file.path(VDJ_ROOT, "Clonal_Visualizations")
dir_diversity <- file.path(VDJ_ROOT, "Clonal_Diversity")
dir_overlap   <- file.path(VDJ_ROOT, "Clonal_Overlap")
dir_tables    <- file.path(VDJ_ROOT, "Tables")
dir_seurat    <- file.path(VDJ_ROOT, "Seurat_Plots")
mkdirs(dir_clonal, dir_diversity, dir_overlap, dir_tables, dir_seurat)

libraries <- CP018_SAMPLES$Sample_Subfolder

# -----------------------------
# 1) Read contigs
# -----------------------------
contig_paths <- file.path(CR_OUT, libraries, "vdj_t", "filtered_contig_annotations.csv")
names(contig_paths) <- libraries
for (i in seq_along(contig_paths)) require_file(contig_paths[i], paste0("TCR contigs for ", libraries[i]))

contig_list <- lapply(contig_paths, read.csv, stringsAsFactors = FALSE)
names(contig_list) <- libraries

contig_summary <- data.frame(
  Sample_Subfolder = libraries,
  n_contigs = sapply(contig_list, nrow),
  n_barcodes = sapply(contig_list, function(d) length(unique(d$barcode)))
)
write.csv(contig_summary, file.path(dir_tables, "CP018_contig_summary.csv"), row.names = FALSE)
print(contig_summary)

# -----------------------------
# 2) (a) Per-library combineTCR — for Seurat attachment
# -----------------------------
combined.TCR <- combineTCR(contig_list, samples = names(contig_list))

samples <- names(combined.TCR)
meta_by_sample <- CP018_SAMPLES[match(samples, CP018_SAMPLES$Sample_Subfolder), ]

combined.TCR <- addVariable(combined.TCR, "PID",         meta_by_sample$PID)
combined.TCR <- addVariable(combined.TCR, "Timepoint",   meta_by_sample$Timepoint)
combined.TCR <- addVariable(combined.TCR, "Replicate",   meta_by_sample$Replicate)
combined.TCR <- addVariable(combined.TCR, "OCM_Barcode", meta_by_sample$OCM_Barcode)

# -----------------------------
# 3) (b) Replicate-merged combineTCR — for overlap/diversity
# -----------------------------
merged_contigs <- list(
  CP018_2m  = dplyr::bind_rows(contig_list[["CP018_2m_A"]],  contig_list[["CP018_2m_B"]]),
  CP018_90m = dplyr::bind_rows(contig_list[["CP018_90m_A"]], contig_list[["CP018_90m_B"]])
)
combined.TCR.merged <- combineTCR(merged_contigs, samples = names(merged_contigs))
combined.TCR.merged <- addVariable(combined.TCR.merged, "PID", rep("CP018", 2))
combined.TCR.merged <- addVariable(combined.TCR.merged, "Timepoint", c("2m", "90m"))

# -----------------------------
# 4) Repertoire visualisations
# -----------------------------
# Per-replicate first: does clonal structure reproduce across A/B?
ggsave(file.path(dir_clonal, "CP018_clonalQuant_by_library.png"),
       clonalQuant(combined.TCR, cloneCall = "strict", chain = "both", scale = TRUE),
       width = 7, height = 5, dpi = 400, bg = "white")

ggsave(file.path(dir_clonal, "CP018_clonalHomeostasis_by_library.png"),
       clonalHomeostasis(combined.TCR, cloneCall = "strict"),
       width = 8, height = 5, dpi = 400, bg = "white")

ggsave(file.path(dir_clonal, "CP018_clonalAbundance_by_library.png"),
       clonalAbundance(combined.TCR, cloneCall = "strict", scale = FALSE),
       width = 7, height = 5, dpi = 400, bg = "white")

# Replicate-merged: timepoint comparison
ggsave(file.path(dir_clonal, "CP018_clonalQuant_by_timepoint.png"),
       clonalQuant(combined.TCR.merged, cloneCall = "strict", chain = "both", scale = TRUE),
       width = 6, height = 5, dpi = 400, bg = "white")

ggsave(file.path(dir_diversity, "CP018_clonalDiversity_by_timepoint.png"),
       clonalDiversity(combined.TCR.merged, cloneCall = "strict", group.by = "Timepoint"),
       width = 8, height = 5, dpi = 400, bg = "white")

# Overlap: 2m vs 90m clonal sharing — the persistence question, and the direct
# analogue of the CP003 1m->101m tracking in Figure 3C.
ggsave(file.path(dir_overlap, "CP018_clonalOverlap_Jaccard.png"),
       clonalOverlap(combined.TCR.merged, cloneCall = "strict", method = "jaccard"),
       width = 6, height = 5, dpi = 400, bg = "white")

# cloneCall = "strict" concatenates V/J genes AND the full nucleotide sequence,
# so each legend entry is ~150 characters and the legend crowds the panel out
# entirely. "aa" gives the CDR3 amino-acid pair, which is the readable
# identifier and the one used in the manuscript.
p_cmp <- clonalCompare(combined.TCR.merged, top.clones = 20,
                       cloneCall = "aa", graph = "alluvial") +
  guides(fill = guide_legend(ncol = 1, keyheight = unit(0.42, "cm"))) +
  theme(legend.text  = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"))
ggsave(file.path(dir_overlap, "CP018_clonalCompare_top20.png"), p_cmp,
       width = 13, height = 7.5, dpi = 400, bg = "white", limitsize = FALSE)

# Same comparison without the legend - the flows are the point, and the
# clonotype identities are in CP018_clonotype_table.csv.
ggsave(file.path(dir_overlap, "CP018_clonalCompare_top20_nolegend.png"),
       p_cmp + theme(legend.position = "none") +
         labs(title = "Top 20 clonotypes, 2m vs 90m",
              subtitle = "Identities in CP018_clonotype_table.csv"),
       width = 7, height = 6.5, dpi = 400, bg = "white")

# Shared-clonotype table (2m vs 90m)
cl_2m  <- unique(na.omit(combined.TCR.merged[["CP018_2m"]]$CTstrict))
cl_90m <- unique(na.omit(combined.TCR.merged[["CP018_90m"]]$CTstrict))
shared <- intersect(cl_2m, cl_90m)

overlap_stats <- data.frame(
  n_clonotypes_2m  = length(cl_2m),
  n_clonotypes_90m = length(cl_90m),
  n_shared         = length(shared),
  jaccard          = round(length(shared) / length(union(cl_2m, cl_90m)), 4),
  pct_of_2m_shared = round(100 * length(shared) / length(cl_2m), 2),
  pct_of_90m_shared= round(100 * length(shared) / length(cl_90m), 2)
)
write.csv(overlap_stats, file.path(dir_tables, "CP018_2m_vs_90m_overlap_stats.csv"),
          row.names = FALSE)
write.csv(data.frame(CTstrict = shared),
          file.path(dir_tables, "CP018_shared_clonotypes_2m_90m.csv"), row.names = FALSE)
cat("\n2m vs 90m clonal overlap:\n"); print(overlap_stats)

# Replicate-level clonotype overlap — a technical reproducibility measure
rep_overlap <- do.call(rbind, lapply(unique(CP018_SAMPLES$Timepoint), function(tp) {
  a <- paste0("CP018_", tp, "_A"); b <- paste0("CP018_", tp, "_B")
  ca <- unique(na.omit(combined.TCR[[a]]$CTstrict))
  cb <- unique(na.omit(combined.TCR[[b]]$CTstrict))
  data.frame(Timepoint = tp, n_A = length(ca), n_B = length(cb),
             n_shared = length(intersect(ca, cb)),
             jaccard  = round(length(intersect(ca, cb)) / length(union(ca, cb)), 4))
}))
write.csv(rep_overlap, file.path(dir_tables, "CP018_replicate_clonotype_overlap.csv"),
          row.names = FALSE)
cat("\nReplicate clonotype overlap (A vs B):\n"); print(rep_overlap)

# -----------------------------
# 5) Attach TCR to Seurat
# -----------------------------
require_file(QS_INTEGRATED, "integrated object")
seu <- qs_read(QS_INTEGRATED)

# Fig_Annotation is created in script 3. Plots use LABELS, not cluster numbers.
if (!"Fig_Annotation" %in% colnames(seu@meta.data)) {
  stop("Fig_Annotation missing - re-run 3.Integration_Anno_Plots.R first ",
       "(it applies CP018_CLUSTER_ANNOTATION).", call. = FALSE)
}
seu$Fig_Annotation <- factor(as.character(seu$Fig_Annotation),
                             levels = intersect(ANNO_ORDER,
                                                unique(as.character(seu$Fig_Annotation))))


# Strip any leading underscore introduced by RenameCells (as in CP003)
seu <- RenameCells(seu, new.names = sub("^_", "", rownames(seu[[]])))

seu <- combineExpression(
  combined.TCR, seu,
  cloneCall  = "strict",
  chain      = "both",
  group.by   = "sample",
  cloneSize  = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500),
  proportion = FALSE
)

# How many cells actually received a clonotype
tcr_recovery <- seu@meta.data %>%
  group_by(Sample_Subfolder) %>%
  summarise(n_cells = n(),
            n_with_TCR = sum(!is.na(CTstrict)),
            pct_with_TCR = round(100 * mean(!is.na(CTstrict)), 1),
            .groups = "drop")
write.csv(tcr_recovery, file.path(dir_tables, "CP018_TCR_recovery_per_library.csv"),
          row.names = FALSE)
cat("\nTCR recovery:\n"); print(tcr_recovery)

# -----------------------------
# 6) Clonal-expansion plots on the UMAP
# -----------------------------
ggsave(file.path(dir_seurat, "CP018_UMAP_cloneSize.png"),
       DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "cloneSize"),
       width = 8, height = 6, dpi = 400, bg = "white")

ggsave(file.path(dir_seurat, "CP018_UMAP_cloneSize_by_Timepoint.png"),
       DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "cloneSize",
                split.by = "Timepoint"),
       width = 13, height = 6, dpi = 400, bg = "white")

# proportion view answers "which subsets are clonally expanded?" - the count
# view is dominated by the 2,885-cell naive compartment and hides everything else
ggsave(file.path(dir_seurat, "CP018_clonalOccupy_cluster.png"),
       clonalOccupy(seu, x.axis = "Fig_Annotation") +
         labs(x = NULL, title = "Clonal expansion by annotation (cell counts)") +
         rotate_x(35) + theme(plot.title = element_text(face = "bold")),
       width = 11, height = 6, dpi = 400, bg = "white")

ggsave(file.path(dir_seurat, "CP018_clonalOccupy_cluster_proportion.png"),
       clonalOccupy(seu, x.axis = "Fig_Annotation", proportion = TRUE) +
         labs(x = NULL, y = "Proportion of cells",
              title = "Clonal expansion by annotation (proportion)") +
         rotate_x(35) + theme(plot.title = element_text(face = "bold")),
       width = 11, height = 6, dpi = 400, bg = "white")

# source data for both
save_csv(as.data.frame(table(Annotation = seu$Fig_Annotation,
                             cloneSize  = seu$cloneSize)),
         file.path(dir_tables, "CP018_cloneSize_by_annotation.csv"))

# -----------------------------
# CSV exports — repertoire source data
# -----------------------------
csv_dir <- file.path(VDJ_ROOT, "tables"); mkdirs(csv_dir)

md <- seu@meta.data
md$cell_barcode <- rownames(md)

# every clonotype with size, timepoint spread and cluster distribution
clono <- md %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict) %>%
  summarise(n_cells        = n(),
            CTaa           = dplyr::first(CTaa),
            CTgene         = dplyr::first(CTgene),
            n_2m           = sum(Timepoint == "2m"),
            n_90m          = sum(Timepoint == "90m"),
            n_timepoints   = dplyr::n_distinct(Timepoint),
            n_replicates   = dplyr::n_distinct(Timepoint_Rep),
            in_both_reps_2m  = sum(Timepoint == "2m"  & Replicate == "A") > 0 &
                               sum(Timepoint == "2m"  & Replicate == "B") > 0,
            in_both_reps_90m = sum(Timepoint == "90m" & Replicate == "A") > 0 &
                               sum(Timepoint == "90m" & Replicate == "B") > 0,
            clusters       = paste(sort(unique(as.character(mnn_clusters_rna))), collapse = ";"),
            .groups = "drop") %>%
  arrange(desc(n_cells))
save_csv(clono, file.path(csv_dir, "CP018_clonotype_table.csv"))

# expanded clonotypes only
save_csv(clono %>% filter(n_cells > 1),
         file.path(csv_dir, "CP018_clonotypes_expanded.csv"))

# clonotypes shared across the two timepoints — the persistence question
save_csv(clono %>% filter(n_2m > 0, n_90m > 0) %>% arrange(desc(n_cells)),
         file.path(csv_dir, "CP018_clonotypes_shared_2m_90m.csv"))

# per-library TCR recovery
save_csv(md %>% group_by(Sample_Subfolder, Timepoint, Replicate) %>%
           summarise(n_cells = n(),
                     n_with_TCR = sum(!is.na(CTstrict)),
                     pct_with_TCR = round(100 * mean(!is.na(CTstrict)), 1),
                     n_clonotypes = dplyr::n_distinct(CTstrict[!is.na(CTstrict)]),
                     n_expanded_cells = sum(clonalFrequency > 1, na.rm = TRUE),
                     .groups = "drop"),
         file.path(csv_dir, "CP018_TCR_recovery_by_library.csv"))

save_csv(composition_table(seu, "cloneSize", "Timepoint"),
         file.path(csv_dir, "CP018_cloneSize_by_Timepoint.csv"))
export_metadata(seu, file.path(csv_dir, "CP018_cell_metadata_withTCR.csv"))

qs_save(seu, file = QS_WITH_TCR)
message("Saved: ", QS_WITH_TCR)
