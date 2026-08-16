############################################################
# CP018 FACT-Seq — central configuration
#
# source() this at the top of every CP018 script so paths live in ONE place.
# (The CP003 scripts hardcode /home/akshay-iyer/... in ~40 places; this avoids
#  repeating that when the project moves machines.)
#
#   source("Scripts/CP018 - FACTseq/config.R")
############################################################

# ── Root paths ────────────────────────────────────────────────────────────────
PROJ_ROOT   <- path.expand("~/Desktop/Projects")
REPO_DIR    <- file.path(PROJ_ROOT, "Repositories", "CD8_Longitudinal")
SAVED_DIR   <- file.path(REPO_DIR, "saved_R_data")          # gitignored
CR_OUT      <- file.path(PROJ_ROOT, "cellranger_out", "CP018_multi_out",
                         "outs", "per_sample_outs")

# ── Analysis output roots ────────────────────────────────────────────────────
CP018_ANALYSIS <- file.path(REPO_DIR, "Analysis", "CP018")
QC_ROOT        <- file.path(CP018_ANALYSIS, "QC")
CONC_ROOT      <- file.path(CP018_ANALYSIS, "Replicate_Concordance")
INT_ROOT       <- file.path(CP018_ANALYSIS, "Integration")
ANNO_ROOT      <- file.path(CP018_ANALYSIS, "Annotation_Plots")
VDJ_ROOT       <- file.path(CP018_ANALYSIS, "VDJ", "TCR")
HIV_ROOT       <- file.path(CP018_ANALYSIS, "HIV_Specific")

# ── Object checkpoints (qs2) ─────────────────────────────────────────────────
QS_SINGLETS   <- file.path(SAVED_DIR, "CP018_postQC_cellcycle_scDblFinder_singlets.qs2")
QS_INTEGRATED <- file.path(SAVED_DIR, "CP018_RNA_integrated_CCA_MNN.qs2")
QS_WITH_TCR   <- file.path(SAVED_DIR, "CP018_RNA_integrated_CCA_MNN_withTCR.qs2")
QS_ANNOTATED  <- file.path(SAVED_DIR, "CP018_HIVSpecificTCR_annotated.qs2")

# CP003 / TARA objects for cross-referencing (already in saved_R_data)
QS_CP003_ANNOT_FIG3 <- file.path(SAVED_DIR, "cp003_fig3_annotated.qs2")
QS_CP003_ANNOT <- file.path(SAVED_DIR, "seu_CP003_HIVSpecificTCR_annotated.qs2")
QS_TARA_CD8    <- file.path(SAVED_DIR, "TARA_cd8_HEI_annotated_final.qs2")
RDATA_TARA_TCR <- file.path(SAVED_DIR, "TARA_TCR_Combined.RData")

# ── Sample sheet ─────────────────────────────────────────────────────────────
# From L. de Armas (2026-08-07). OCM = on-chip multiplexing barcode.
# Both timepoints run in duplicate: 2m = lanes 1,2 / 90m = lanes 3,4.
#
# NOTE: CellType_Sort is constant "CD8+" here — unlike CP003, where OCM 004 was
# a CD8- fraction. FACT-Seq sorts proliferating (CellTrace-low) CD8+ T cells,
# and A/B are technical replicates of the SAME sort. If Lesley confirms a
# different sort for any lane, change it here only.
CP018_SAMPLES <- data.frame(
  Sample_Subfolder = c("CP018_2m_A", "CP018_2m_B", "CP018_90m_A", "CP018_90m_B"),
  OCM_Barcode      = c("OB1",        "OB2",        "OB3",         "OB4"),
  Tube_Label       = c("1",          "2",          "3",           "4"),
  Replicate        = c("A",          "B",          "A",           "B"),
  Months           = c(2L,           2L,           90L,           90L),
  CellType_Sort    = "CD8+",
  PID              = "CP018",
  stringsAsFactors = FALSE
)
CP018_SAMPLES$Timepoint <- paste0(CP018_SAMPLES$Months, "m")

# ── Shared aesthetics ────────────────────────────────────────────────────────
TP_COLORS  <- c("2m" = "#4A90D9", "90m" = "#E76F51")
REP_COLORS <- c("A"  = "#52B788", "B"   = "#B07AA1")

# ── QC thresholds (identical to CP003 for comparability) ─────────────────────
QC_MIN_COUNT    <- 500
QC_MIN_FEATURE  <- 600
QC_MIN_COMPLEX  <- 0.80
QC_MAX_MITO     <- 15
QC_MIN_RIBO     <- 5
QC_MAX_HB       <- 20
QC_MAX_PLAT     <- 2
QC_MIN_CELLS_PER_GENE <- 10

# ── HIV-specific definition (CP018) ──────────────────────────────────────────
# FACT-Seq sorts and stimulates with HIV peptide, then reads out which
# clonotypes PROLIFERATED. Clonal expansion after antigen stimulation IS the
# assay readout, so expansion alone defines specificity here.
#
# WHY THIS DIFFERS FROM CP003. The published CP003 definition was
#   clonalFrequency > 1 AND cluster not in {bystander clusters}
# where the bystander clusters were hand-picked cluster NUMBERS. Cluster
# numbering is not portable across resolutions or datasets - and it silently
# broke here: objects built at res 0.6 were filtered with cluster numbers
# derived at res 0.4, which excluded the Cycling Effector population as
# "bystander". Defining specificity by expansion alone removes that dependency
# entirely and makes the call reproducible.
#
# THRESHOLD. n > 2 (i.e. >= 3 cells) rather than > 1. A 2-cell "clone" is one
# division, or a barcode collision; 3+ cells is a population that responded.
HIV_MIN_CLONE_SIZE   <- 2      # clonalFrequency STRICTLY GREATER than this
HIV_PER_TIMEPOINT    <- TRUE   # count expansion WITHIN a timepoint, not pooled
                               # (pooling lets 2m and 90m cells sum to a "clone"
                               #  that never expanded at either timepoint)
HIV_EXCLUDE_MAIT     <- TRUE   # MAIT TCRs are semi-invariant: unrelated cells
                               # share a clonotype, so expansion is guaranteed
                               # by receptor biology, not antigen response.
                               # Applied by ANNOTATION LABEL, never by number.

# ── Clustering resolution (set after inspecting clustree) ────────────────────
# Clusters smaller than this are dropped before annotation - they cannot
# support a marker test. Raise it if clustering yields extra small clusters.
MIN_CLUSTER_SIZE <- 10

MNN_RES_FINAL <- 0.4   # CP018: chosen from clustree + marker profiling (CP003 used 0.6)
CCA_RES_FINAL <- 0.4

# ── Figure 3 annotation vocabulary (matches the published CP003 figure) ──────
# Keep these names and colours identical to Scripts/3. Manuscript/Figure_3.R so
# a CP018 panel can sit beside the CP003 one without a legend mismatch.
ANNO_ORDER <- c("Naive/Bystander", "Naive (hypoxic)", "Transitional Tem",
                "Activated Stem-like", "Tpex", "TEMRA/Effector",
                "Cycling Effector", "MAIT", "Tex")
# NOTE: CP018 has NO Tex cluster (see below). "Tex" is retained in the vocabulary
# and palette so a CP018 panel can sit beside the published CP003 figure without
# a legend mismatch; it will simply be unused here.
# "MAIT" is NEW for CP018 - CP003 had no distinct MAIT cluster.

ANNO_COLS <- c(
  "Naive/Bystander"  = "#2166AC",
  "Naive (hypoxic)"  = "#6BAED6",   # naive surface phenotype + hypoxia response
  "Transitional Tem" = "#1A9850",
  "Activated Stem-like" = "#F4A261",   # CP018 cluster 2 (see note below)
  "Tpex"             = "#E9C46A",
  "TEMRA/Effector"   = "#E63946",
  "Cycling Effector" = "#9D4EDD",
  "MAIT"             = "#00838F",
  "Tex"              = "#6D4C41"
)

HIV_COLS   <- c("Other" = "#E0E0E0", "HIV-Specific TCR" = "#E63946")
PHASE_COLS <- c("G1" = "#AEC6CF", "S" = "#FFD700", "G2M" = "#E63946")

# ══════════════════════════════════════════════════════════════════════════════
# MANUAL STEP — map CP018 MNN clusters onto the annotation vocabulary above.
# Cluster numbering is dataset-specific, so this CANNOT be inherited from CP003.
# Fill in after reviewing script 3's marker plots and script 5's cluster profile
# (CP018_cluster_profile_for_bystander_call.csv).
#
#   Naive/Bystander   high TCF7/LEF1/SELL/CCR7/IL7R, no clonal expansion
#   Transitional Tem  intermediate; GZMK+, some CCR7
#   Tpex              TCF7+ with TOX/PDCD1 — stem-like exhausted
#   TEMRA/Effector    GZMB/GNLY/PRF1/FGFBP2/CX3CR1 high, TCF7 low
#   Cycling Effector  MKI67/TOP2A/PCNA high, S/G2M phase
#   Tex               TOX/PDCD1/HAVCR2/LAG3/TIGIT high, low stemness
#
# Example once known:
#   CP018_CLUSTER_ANNOTATION <- c("0" = "Naive/Bystander", "1" = "Naive/Bystander",
#                                 "2" = "TEMRA/Effector",  "3" = "Cycling Effector")
CP018_CLUSTER_ANNOTATION <- c(
  "0" = "Naive/Bystander",   # SELL 12.4  CCR7 6.6  TCF7 5.2  IL7R 14.1  GZMB 0.02   n=2404
  "1" = "Transitional Tem",  # GZMK 8.1   XCL1 3.7  TOX 1.2   SELL 2.3   CCR7 1.0    n=496
  "2" = "Activated Stem-like",
  #   TCF7 10.0 (highest) ID3 5.66  CD27 6.02  SELL 13.1  CCR7 8.5  LEF1 3.9
  #   CD38 2.81 (10x naive) IL7R 2.50 (naive 14.1) -> recently TCR-activated
  #   PDCD1 0.09  TOX 0.62  SLAMF6 0.63 -> exhaustion programme NOT engaged
  #   NOT called "Tpex": Tpex requires the exhaustion programme (PDCD1/TOX/SLAMF6
  #   high). CP003's equivalent stem-like cluster (c3) has ENTPD1 0.98, LAG3 4.07,
  #   HAVCR2 1.59 and IS a genuine Tpex. CP018's does not.
  #   NOT called "Tscm" either: classic Tscm is IL7R-HIGH; here IL7R is low, which
  #   is expected after peptide stimulation. Whether these are Tscm at rest is
  #   testable in the unstimulated TARA data (script 6) - revisit the label then.
  #
  #   EXPANDED-PANEL CHECK (script 3b): Stem_TF z=+1.5 (highest), Tscm panel
  #   z=+1.0, Activation_late z=+1.0, Inhibitory z=+1.1, Exhaustion_TF z=+0.7,
  #   OXPHOS z=+1.3, Glycolysis z=+0.4. Top DE vs naive: XCL2, XCL1, ENTPD1,
  #   CD70, KLRC1, DUSP4, HOPX, TNFRSF18.
  #   XCL1/XCL2 + TCF7 + ID3 IS the Tpex hallmark, and ENTPD1 (CD39) is a
  #   chronic-activation marker - so this cluster sits BETWEEN activated
  #   stem-like and Tpex. It is called "Activated Stem-like" because PDCD1
  #   (0.09) and TOX (0.62) remain low: the exhaustion programme is not
  #   engaged. If script 6 shows these clonotypes are exhausted at rest, Tpex
  #   becomes the better label.
  #   CAVEAT: mean nFeature 5218 and nCount 29928 are ~2x other non-cycling
  #   clusters and the mean doublet score (0.236) is the highest outside the
  #   cycling cluster. CD4 is absent (0.05) and CD3E/CD8A/CD8B are normal, so
  #   this is NOT a CD4/CD8 or T/myeloid doublet - high transcriptional output
  #   is expected in recently activated blasts. Flagged in
  #   CP018_cluster_QC_confounders.csv; worth a second look if this cluster
  #   drives a key result.
  "3" = "TEMRA/Effector",    # GNLY 76.5  FGFBP2 4.8 GZMB 5.2  TCF7 1.3              n=426
  "4" = "Cycling Effector",  # MKI67 3.4  GZMB 13.1 LAG3 3.2  97% S/G2M              n=408
  "5" = "Naive/Bystander",   # SELL 12.4  CCR7 6.2  IL7R 14.5                        n=250
  "6" = "MAIT",              # KLRB1 15.3 TRAV1-2 2.7 SLC4A10 1.1 ZBTB16 0.54        n=246
  "7" = "Naive (hypoxic)",
  #   Naive surface intact (SELL 9.5, CCR7 5.2, IL7R 13.7) BUT a dominant
  #   hypoxia/glycolytic response: hypoxia panel z=+2.6, glycolysis z=+1.6,
  #   OXPHOS z=-1.1. Top DE vs other naive: MIR210HG +4.6, DDIT4 +3.9,
  #   P4HA1 +3.2, BNIP3 +2.9, ALDOC +2.3, ANKRD37 +2.1, EGLN1 +1.9 (all
  #   canonical HIF1A targets). QC is otherwise normal (mito 2.6%,
  #   complexity 0.865, doublet score 0.042) so this is a CELL STATE, not
  #   low-quality debris - most likely cells that experienced hypoxic stress
  #   during handling/culture. Kept as a bystander for the HIV-specific call
  #   (phenotypically naive, not antigen-expanded) but labelled separately so
  #   the hypoxia signal is not silently read as biology.
  "8" = "Naive/Bystander"    # SELL 12.6  CCR7 6.3  IL7R 14.5                        n=231
)
# Cluster 9 (2 cells) is dropped in script 3 as a singleton artefact.
# ══════════════════════════════════════════════════════════════════════════════

# ── Namespace conflict guard ─────────────────────────────────────────────────
# Bioconductor packages (SingleCellExperiment/S4Vectors/IRanges, pulled in by
# scDblFinder, batchelor, scRepertoire) mask several dplyr/tidyr verbs. If you
# run the scripts sequentially in ONE R session, script 1's attachments are
# still on the search path when script 2 runs, and dplyr::count resolves to the
# wrong function ("Argument 'x' is not a vector: list").
#
# Binding the dplyr/tidyr versions into the global environment fixes this for
# every downstream script: globalenv is searched before any attached package,
# so these win regardless of what gets attached later.
.prefer_verbs <- function(pkg, fns) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(invisible(NULL))
  ns <- asNamespace(pkg)
  for (f in fns) {
    if (exists(f, envir = ns, inherits = FALSE)) {
      assign(f, get(f, envir = ns), envir = globalenv())
    }
  }
}

.prefer_verbs("dplyr", c("count","filter","select","rename","slice","mutate",
                         "summarise","summarize","arrange","group_by","ungroup",
                         "first","last","desc","n","n_distinct","bind_rows",
                         "left_join","inner_join","full_join","transmute",
                         "recode","pull","distinct","across","all_of","case_when",
                         "if_else","lag","lead","between","setdiff","union",
                         "intersect","collapse"))
.prefer_verbs("tidyr", c("pivot_longer","pivot_wider","expand","separate","unite",
                         "replace_na","drop_na","fill","complete"))

# ── Helpers ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(ggrepel))

mkdirs <- function(...) invisible(lapply(c(...), dir.create,
                                        recursive = TRUE, showWarnings = FALSE))

# ── Plot saving ──────────────────────────────────────────────────────────────
# Publication defaults. The old base-png() device at res=200 produced crowded
# axis labels and small text; ggsave with an explicit inch canvas + 350 dpi
# matches the existing TARA figures (5400x3600 @ 300 dpi).
FIG_DPI <- 350

# Shared theme: larger type, breathing room, no chart junk.
theme_cp018 <- function(base_size = 16) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = base_size * 1.25, face = "bold",
                                              hjust = 0.5, margin = ggplot2::margin(b = 8)),
      plot.subtitle   = ggplot2::element_text(size = base_size * 0.95, hjust = 0.5,
                                              colour = "grey30"),
      axis.title      = ggplot2::element_text(size = base_size * 1.05, face = "bold"),
      axis.text       = ggplot2::element_text(size = base_size * 0.9, colour = "black"),
      axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5),
      strip.text      = ggplot2::element_text(size = base_size * 0.95, face = "bold"),
      strip.background= ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.title    = ggplot2::element_text(size = base_size * 0.95, face = "bold"),
      legend.text     = ggplot2::element_text(size = base_size * 0.9),
      plot.margin     = ggplot2::margin(12, 14, 10, 12)
    )
}

# width/height are INCHES (not pixels). Defaults suit a single panel; pass
# larger values for multi-panel patchworks.
save_png <- function(filename, plot_obj, width = 9, height = 6.5,
                     dpi = FIG_DPI, limitsize = FALSE) {
  # Back-compat: earlier calls passed pixels (e.g. width = 1800). Detect and convert.
  if (width  > 60) width  <- width  / 200
  if (height > 60) height <- height / 200
  ggplot2::ggsave(filename = filename, plot = plot_obj,
                  width = width, height = height, dpi = dpi,
                  bg = "white", limitsize = limitsize)
}

# Size a multi-panel figure by its panel count so nothing gets crushed.
panel_size <- function(n_panels, ncol = 2, panel_w = 5.0, panel_h = 3.6) {
  nrow <- ceiling(n_panels / ncol)
  c(width = ncol * panel_w, height = nrow * panel_h)
}

# ── CSV export ───────────────────────────────────────────────────────────────
# Every figure should ship the numbers behind it. save_csv() writes to the same
# directory as the plot and echoes the path, so nothing is silently skipped.
save_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = row.names)
  message("  [csv] ", basename(path), "  (", nrow(x), " rows)")
  invisible(path)
}

# Cell-level metadata dump: the master table every downstream question can be
# answered from without reloading the Seurat object.
export_metadata <- function(obj, path, extra_cols = NULL) {
  md <- obj@meta.data
  md$cell_barcode <- rownames(md)
  keep <- c("cell_barcode", intersect(
    c("Sample_Subfolder","Sample","Timepoint","Replicate","Timepoint_Rep","OCM_Barcode",
      "nCount_RNA","nFeature_RNA","percent_mito","percent_ribo","log10GenesPerUMI",
      "S.Score","G2M.Score","Phase","scDblFinder.class",
      "predicted.celltype.l1","predicted.celltype.l2",
      "mnn_clusters_rna","cca_clusters_rna","Fig_Annotation",
      "CTgene","CTaa","CTnt","CTstrict","clonalFrequency","clonalProportion","cloneSize",
      "HIV_Specific_TCR","Stemness_score1","Cytotoxicity_score1","Exhaustion_score1",
      "IFN_score1"),
    colnames(md)), extra_cols)
  save_csv(md[, unique(keep), drop = FALSE], path)
}

# Composition table: counts + within-group percentages for any two factors.
composition_table <- function(obj, group, split) {
  md <- obj@meta.data
  tb <- as.data.frame(table(group = md[[group]], split = md[[split]]))
  names(tb) <- c(group, split, "n_cells")
  tb <- tb[tb$n_cells > 0 | TRUE, ]
  tot <- tapply(tb$n_cells, tb[[split]], sum)
  tb$pct_of_sample <- round(100 * tb$n_cells / tot[as.character(tb[[split]])], 2)
  tb
}

# Annotation labels are words ("Activated Stem-like"), not digits, so any plot
# with annotation on the x-axis needs rotated tick labels or they collide.
rotate_x <- function(angle = 35, size = NULL) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = angle, hjust = 1, vjust = 1, size = size))
}

# ── UMAP presentation ────────────────────────────────────────────────────────
# Seurat names UMAP axes after the reduction key, giving "umapmnnrna_1" on the
# axis. Also, on-plot cluster labels drawn straight onto dark points are
# unreadable, and duplicate the legend. tidy_umap() fixes both.
#
#   tidy_umap(DimPlot(...))                      # legend only (default)
#   tidy_umap(DimPlot(..., label = TRUE))        # keeps labels, adds halo
tidy_umap <- function(p, title = NULL, subtitle = NULL,
                      xlab = "UMAP 1", ylab = "UMAP 2", legend = TRUE) {
  p <- p + ggplot2::labs(x = xlab, y = ylab)
  if (!is.null(title))    p <- p + ggplot2::ggtitle(title)
  if (!is.null(subtitle)) p <- p + ggplot2::labs(subtitle = subtitle)
  p <- p + theme_cp018() +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      axis.title    = ggplot2::element_text(face = "plain"),
      legend.key    = ggplot2::element_blank())
  if (!legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}

# Cluster labels with a white halo so they stay legible over dense points.
# Use INSTEAD of DimPlot(label = TRUE), which has no contrast control.
umap_labels <- function(p, size = 4.2) {
  d   <- p$data
  xy  <- names(d)[1:2]
  grp <- if ("ident" %in% names(d)) "ident" else names(d)[3]
  cen <- stats::aggregate(d[, xy], by = list(lab = d[[grp]]), FUN = stats::median)
  p + ggrepel::geom_label_repel(
        data = cen, ggplot2::aes(x = .data[[xy[1]]], y = .data[[xy[2]]], label = lab),
        inherit.aes = FALSE, size = size, fontface = "bold",
        fill = grDevices::adjustcolor("white", alpha.f = 0.78),
        label.size = 0.15, box.padding = 0.35, seed = 1,
        min.segment.length = 0.2, segment.colour = "grey35")
}

# Fail loudly and early if an input object is missing.
require_file <- function(path, what = "input") {
  if (!file.exists(path)) {
    stop(what, " not found:\n  ", path,
         "\nCheck config.R paths, or that the upstream script has been run.",
         call. = FALSE)
  }
  invisible(path)
}

message("CP018 config loaded | ", nrow(CP018_SAMPLES), " libraries | CR_OUT: ", CR_OUT)
