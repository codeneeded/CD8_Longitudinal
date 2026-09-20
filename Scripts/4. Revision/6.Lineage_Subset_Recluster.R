################################################################################
# Revision script 6 — subset, re-integrate and re-cluster each lineage
#
# ANSWERS
#   Rev 1 major 2  "The focus on CD8 T cells is briefly explained ... but it
#                   still seems odd that the investigators did not perform any
#                   substantial analyses of any other cell lineages (CD4 T
#                   cells, B cells, monocytes, DCs, NK cells etc). Given that
#                   this is all sitting there in the dataset it seems like it
#                   would be fairly easy to look for patterns in these cells."
#
# WHY RE-DERIVE RATHER THAN USE THE EXISTING SUBSETS
#   saved_R_data/ already holds TARA_CD4_subset.qs2 (51,660 cells, 23 clusters),
#   TARA_NK_subset.qs2 (16,877, 10) and TARA_B_cell_subset.qs2 (17,671, 9) - each
#   with its own wsnn graph and wnn.umap, so each WAS re-clustered. But they are
#   STALE relative to the published object:
#     - their parent labels are from a superseded annotation round
#       ("Naive CD4" vs the published "Naive 1/2/3 CD4"; "CD56dim NK" vs
#        "CD56dim CD122hi/int/lo NK")
#     - only 94.4%, 95.1% and 97.3% of their barcodes still exist in
#       TARA_ALL_annotated_final.qs2
#     - none is traceable to a script in Scripts/
#   Annotating them would produce per-lineage counts that DISAGREE with Figure 1,
#   which is exactly what Rev 2 major 1a (per-participant cell counts) exposes.
#
#   NOTE the CD8 case is different and already answered: the published
#   TARA_cd8_HEI_annotated_final.qs2 (18,843 cells, 10 populations, res 0.4) IS a
#   dedicated re-clustering and backs Figures 4-5. TARA_CD8_subset.qs2 (29,978
#   cells, 24 clusters) is the same superseded round as the others - do not use.
#
# METHOD — deliberately identical to 4.CD8_Subseting.R so the lineages are
# directly comparable to the published CD8 analysis:
#   subset by published Annotation -> re-normalise RNA -> HVG -> PCA -> MNN
#   -> ADT scale -> aPCA -> MNN -> FindMultiModalNeighbors -> WNN UMAP
#   -> FindClusters over 0.2-2.0 -> clustree
#
# OUTPUT
#   Analysis/<LIN>_subset/01_integration_diagnostics/*.png   QC + integration
#   Analysis/<LIN>_subset/02_clustering_exploration/*.png     clustree, res sweep
#   saved_R_data/rev_<lin>_reclustered.qs2                    checkpoint
#
#   NOTHING here goes to Reviewer_Comments/ - this is QC and integration, not a
#   reviewer-facing result. Script 7 annotates; script 8 produces the panels
#   that answer the comment.
#
# USAGE — one lineage per invocation so a failure never loses the others:
#   Rscript "Scripts/4. Revision/6.Lineage_Subset_Recluster.R" CD4
#   Rscript "Scripts/4. Revision/6.Lineage_Subset_Recluster.R" NK
#   Rscript "Scripts/4. Revision/6.Lineage_Subset_Recluster.R" B
#   Rscript "Scripts/4. Revision/6.Lineage_Subset_Recluster.R" Myeloid
################################################################################

.cfg <- Sys.getenv("CD8_REV_CONFIG", "")
if (!nzchar(.cfg)) {
  for (b in c("~/Documents/CD8_Longitudinal",
              "~/Desktop/Projects/Repositories/CD8_Longitudinal")) {
    p <- file.path(path.expand(b), "Scripts", "4. Revision", "config_revision.R")
    if (file.exists(p)) { .cfg <- p; break }
  }
}
source(.cfg)
suppressPackageStartupMessages(library(clustree))

# INTEGRATION METHOD - deliberately explicit, never silently substituted.
#   The published CD8 analysis (4.CD8_Subseting.R) used FastMNNIntegration, which
#   lives in SeuratWrappers and needs batchelor. If those are absent the script
#   STOPS rather than falling back to CCA: a different integration method would
#   make these lineages non-comparable to the published CD8 populations, and a
#   silent substitution is exactly the kind of thing a reviewer would catch.
#
#   Set INTEGRATION=cca to override deliberately (documented in the output), e.g.
#   for a quick look on a machine without batchelor.
INTEGRATION <- Sys.getenv("INTEGRATION", "mnn")

# ── PREFLIGHT: prove this interpreter can finish BEFORE anything is destroyed ──
# A caller that clears checkpoints and then discovers the packages are missing
# has destroyed work for nothing. Every package the script needs end-to-end is
# checked here, at the top, before a single file is written or removed.
.need <- c("Seurat", "SeuratObject", "qs2", "clustree", "dplyr", "tidyr", "ggplot2")
if (INTEGRATION == "mnn") .need <- c(.need, "SeuratWrappers", "batchelor")
.absent <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.absent)) {
  stop("PREFLIGHT FAILED - nothing has been written or deleted.\n",
       "  missing: ", paste(.absent, collapse = ", "), "\n",
       "  running: ", R.home(), " (", R.version.string, ")\n",
       "  This interpreter cannot complete the run. Do NOT clear checkpoints\n",
       "  before re-running with an interpreter where these load.\n",
       "  If SeuratWrappers/batchelor are the only gap, INTEGRATION=cca works\n",
       "  but is NOT comparable to the published FastMNN CD8 analysis.",
       call. = FALSE)
}
message("Preflight OK: ", length(.need), " packages available")

if (INTEGRATION == "mnn") {
  missing <- c("SeuratWrappers", "batchelor")[
    !vapply(c("SeuratWrappers", "batchelor"), requireNamespace, logical(1),
            quietly = TRUE)]
  if (length(missing)) {
    stop("FastMNN integration needs: ", paste(missing, collapse = ", "),
         "\n  These ARE installed in your host R (SeuratWrappers at ",
         "~/R/x86_64-pc-linux-gnu-library/4.6/, batchelor in ",
         "/usr/lib/R/site-library/), so run this script with your own Rscript:",
         "\n    Rscript \"Scripts/4. Revision/6.Lineage_Subset_Recluster.R\" CD4",
         "\n  Or set INTEGRATION=cca to use Seurat's built-in CCA instead - but ",
         "note\n  that makes these lineages NOT directly comparable to the ",
         "published CD8\n  analysis, which used FastMNN. State it in Methods if ",
         "you do.", call. = FALSE)
  }
  suppressPackageStartupMessages({library(SeuratWrappers); library(batchelor)})
  INT_METHOD <- FastMNNIntegration
  INT_LABEL  <- "FastMNN (matches the published CD8 analysis)"
} else if (INTEGRATION == "cca") {
  INT_METHOD <- CCAIntegration
  INT_LABEL  <- "CCA (Seurat built-in) - NOT the published CD8 method"
  warning("Using CCA, not FastMNN. These lineages are not directly comparable ",
          "to the published CD8 analysis.", call. = FALSE)
} else {
  stop("INTEGRATION must be 'mnn' or 'cca', got: ", INTEGRATION, call. = FALSE)
}
message("Integration method: ", INT_LABEL)

SEED <- 42
set.seed(SEED)

args <- commandArgs(trailingOnly = TRUE)
LIN  <- if (length(args)) args[1] else stop(
  "Usage: Rscript 6.Lineage_Subset_Recluster.R <CD4|NK|B|Myeloid>", call. = FALSE)

# Which published Annotation values belong to each lineage. Matched by regex on
# the real level names so a renamed population cannot be silently dropped.
LIN_PATTERNS <- list(
  CD4     = "CD4|Treg",
  NK      = "NK",
  B       = "B cell|Naive B|Transitional B|Activated B|Follicular B|Resting Naive B|TNF\\+ B",
  Myeloid = "Monocyte|DC|pDC"
)
if (!LIN %in% names(LIN_PATTERNS))
  stop("LIN must be one of: ", paste(names(LIN_PATTERNS), collapse = ", "), call. = FALSE)

OUT  <- file.path(BASE, "Analysis", paste0(LIN, "_subset"))
D_INT <- file.path(OUT, "01_integration_diagnostics")
D_CLU <- file.path(OUT, "02_clustering_exploration")
mkdirs(D_INT, D_CLU)
CKPT <- paste0(tolower(LIN), "_reclustered")

cat(strrep("=", 70), "\n  LINEAGE:", LIN, "\n", strrep("=", 70), "\n\n")

if (has_qs(CKPT)) {
  # A checkpoint is only reusable if it was built with the SAME integration
  # method. Otherwise a CCA object would silently satisfy a FastMNN run.
  prev <- tryCatch({
    mf <- file.path(REV_TAB, "rev_object_manifest.csv")
    if (file.exists(mf)) {
      m <- read.csv(mf, stringsAsFactors = FALSE)
      n <- m$note[m$name == CKPT]
      if (length(n)) tail(n, 1) else ""
    } else ""
  }, error = function(e) "")
  prev_method <- if (grepl("FastMNN", prev, fixed = TRUE)) "mnn"
                 else if (grepl("CCA", prev, fixed = TRUE)) "cca" else "unknown"

  if (prev_method == INTEGRATION || prev_method == "unknown" &&
      !nzchar(Sys.getenv("INTEGRATION"))) {
    cat("Checkpoint exists:", rev_qs_path(CKPT), "\n")
    cat("  built with:", prev_method, "- matches the requested method.\n")
    cat("Delete it to recompute. Nothing to do.\n")
    quit(save = "no", status = 0)
  }
  cat("Checkpoint exists but was built with a DIFFERENT integration method:\n")
  cat("    on disk:  ", prev_method, "\n")
  cat("    requested:", INTEGRATION, "\n")
  cat("  Recomputing - a", prev_method, "object must not be passed off as",
      paste0(INTEGRATION, "."), "\n\n")
}

################################################################################
# 1. Subset from the PUBLISHED object, so counts reconcile with Figure 1
################################################################################
require_obj(QS_ALL)
all_obj <- qs_read(QS_ALL)
stopifnot("Annotation" %in% colnames(all_obj@meta.data))

lvls  <- unique(as.character(all_obj$Annotation))
keep  <- grep(LIN_PATTERNS[[LIN]], lvls, value = TRUE)
# CD4 must not swallow "CD4dim Naive T cell" if that is not a conventional CD4;
# report what matched so the selection is auditable rather than implicit.
cat("Populations selected for", LIN, ":\n")
for (k in sort(keep)) cat("   ", k, " (", sum(all_obj$Annotation == k), " cells )\n", sep = "")
if (!length(keep)) stop("No Annotation levels matched for ", LIN, call. = FALSE)

cells <- colnames(all_obj)[all_obj$Annotation %in% keep]
obj   <- subset(all_obj, cells = cells)
rm(all_obj); gc(verbose = FALSE)

n0 <- ncol(obj)
cat("\nSubset:", format(n0, big.mark = ","), "cells |",
    length(unique(obj$orig.ident)), "samples |",
    length(unique(obj$PID)), "participants\n")

# Integration dimensions cannot exceed the cells in the SMALLEST batch, and
# very small batches make an unstable anchor set regardless of method. Myeloid
# has samples with 19 and 26 cells, which fails outright at 30 dims.
#   - drop samples below MIN_CELLS_PER_SAMPLE from the integration
#   - cap dims at what the smallest RETAINED sample supports
# Both are reported and written to a CSV so the exclusion is never silent.
MIN_CELLS_PER_SAMPLE <- as.integer(Sys.getenv("MIN_CELLS_PER_SAMPLE", "50"))
per_samp <- table(obj$orig.ident)
too_small <- names(per_samp)[per_samp < MIN_CELLS_PER_SAMPLE]

if (length(too_small)) {
  cat("\nSamples below", MIN_CELLS_PER_SAMPLE, "cells, EXCLUDED from",
      LIN, "integration:\n")
  for (s in too_small) cat("   ", s, ":", per_samp[[s]], "cells\n")
  write.csv(data.frame(sample = too_small,
                       n_cells = as.integer(per_samp[too_small]),
                       reason = paste0("fewer than ", MIN_CELLS_PER_SAMPLE,
                                       " cells; too few to anchor reliably")),
            file.path(OUT, paste0(LIN, "_samples_excluded.csv")), row.names = FALSE)
  obj <- subset(obj, cells = colnames(obj)[!obj$orig.ident %in% too_small])
  obj$orig.ident <- droplevels(factor(obj$orig.ident))
  cat("  retained:", format(ncol(obj), big.mark = ","), "cells in",
      length(unique(obj$orig.ident)), "samples\n")
}

smallest <- min(table(obj$orig.ident))
DIMS_RNA <- min(30, smallest - 1L)
DIMS_ADT <- min(20, smallest - 1L)
# CCA anchor weighting also cannot exceed the smallest batch (default k.weight
# = 100 fails on a 61-cell sample). MNN ignores this argument, so it is passed
# only when it applies.
K_WEIGHT <- min(100L, max(5L, smallest - 1L))   # CCA only; MNN ignores it
cat("\nsmallest retained sample:", smallest, "cells ->",
    "RNA dims 1:", DIMS_RNA, " ADT dims 1:", DIMS_ADT,
    if (INTEGRATION == "cca") paste0(" k.weight ", K_WEIGHT) else "",
    "\n", sep = "")
if (DIMS_RNA < 10)
  stop("Smallest sample has only ", smallest, " cells, leaving ", DIMS_RNA,
       " dimensions - too few for a meaningful integration. Raise ",
       "MIN_CELLS_PER_SAMPLE.", call. = FALSE)

sel <- data.frame(lineage = LIN, population = sort(keep),
                  n_cells = vapply(sort(keep),
                                   function(k) sum(obj$Annotation == k), integer(1)))
write.csv(sel, file.path(OUT, paste0(LIN, "_populations_selected.csv")), row.names = FALSE)

################################################################################
# 2. QC before anything else - the reviewer asks for cell counts that passed QC
################################################################################
qcs <- obj@meta.data |>
  dplyr::group_by(Sample = orig.ident, PID, Timepoint_Group, Condition) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   median_genes = median(nFeature_RNA),
                   median_umi   = median(nCount_RNA),
                   median_pct_mt = if ("percent.mt" %in% names(obj@meta.data))
                     round(median(percent.mt), 2) else NA_real_,
                   .groups = "drop")
write.csv(qcs, file.path(D_INT, paste0(LIN, "_qc_per_sample.csv")), row.names = FALSE)

p_qc <- ggplot(qcs, aes(x = reorder(Sample, -n_cells), y = n_cells,
                        fill = Timepoint_Group)) +
  geom_col() +
  scale_fill_manual(values = ART_COLORS, labels = ART_LABELS, name = NULL) +
  labs(title = paste0(LIN, " cells per sample after QC"),
       subtitle = paste0(format(n0, big.mark = ","), " cells across ",
                         nrow(qcs), " samples"),
       x = NULL, y = "Cells") +
  theme_rev(10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7.5))
ggsave(file.path(D_INT, paste0(LIN, "_cells_per_sample.png")), p_qc,
       width = 11, height = 5.5, dpi = 400, bg = "white")

################################################################################
# 3. RNA: re-normalise, HVG, PCA, MNN integration by sample
################################################################################
cat("\n-- RNA re-processing --\n")
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]])) > 3) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
obj <- FindVariableFeatures(obj, assay = "RNA", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, assay = "RNA", verbose = FALSE)
obj <- RunPCA(obj, assay = "RNA", npcs = 50, reduction.name = "pca",
              seed.use = SEED, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:DIMS_RNA, seed.use = SEED,
               reduction.name = "umap.unintegrated", verbose = FALSE)

# FastMNN takes ONLY assay + new.reduction, exactly as 4.CD8_Subseting.R does.
# Passing dims/orig.reduction sends them through ... into batchelor, which reads
# `dims` as another batch and dies with "number of rows is not the same across
# batches (see batch \"dims\")". CCA needs dims and k.weight; MNN must not get them.
int_args_rna <- if (INTEGRATION == "mnn") {
  list(object = obj, method = INT_METHOD, assay = "RNA",
       new.reduction = "integrated.mnn.rna")
} else {
  list(object = obj, method = INT_METHOD, orig.reduction = "pca",
       dims = 1:DIMS_RNA, k.weight = K_WEIGHT,
       new.reduction = "integrated.mnn.rna", verbose = FALSE)
}
obj <- do.call(IntegrateLayers, int_args_rna)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- RunUMAP(obj, reduction = "integrated.mnn.rna", dims = 1:DIMS_RNA, seed.use = SEED,
               reduction.name = "umap.mnn.rna", verbose = FALSE)

for (nm in c("umap.unintegrated", "umap.mnn.rna")) {
  pp <- DimPlot(obj, reduction = nm, group.by = "orig.ident", raster = FALSE) +
    NoLegend() + ggtitle(paste0(LIN, " - ", nm, " (by sample)"))
  ggsave(file.path(D_INT, paste0(LIN, "_", nm, "_by_sample.png")), pp,
         width = 7.5, height = 7, dpi = 350, bg = "white")
}

################################################################################
# 4. ADT: scale, aPCA, MNN
################################################################################
cat("-- ADT re-processing --\n")
DefaultAssay(obj) <- "ADT"
if (length(Layers(obj[["ADT"]])) > 3) obj[["ADT"]] <- JoinLayers(obj[["ADT"]])
VariableFeatures(obj) <- rownames(obj[["ADT"]])
obj[["ADT"]] <- split(obj[["ADT"]], f = obj$orig.ident)

obj <- ScaleData(obj, assay = "ADT", verbose = FALSE)
obj <- RunPCA(obj, assay = "ADT", npcs = 30, reduction.name = "apca",
              reduction.key = "aPC_", seed.use = SEED, verbose = FALSE)
int_args_adt <- if (INTEGRATION == "mnn") {
  list(object = obj, method = INT_METHOD, assay = "ADT",
       new.reduction = "integrated.mnn.adt")
} else {
  list(object = obj, method = INT_METHOD, orig.reduction = "apca",
       assay = "ADT", dims = 1:DIMS_ADT, k.weight = K_WEIGHT,
       new.reduction = "integrated.mnn.adt", verbose = FALSE)
}
obj <- do.call(IntegrateLayers, int_args_adt)
obj[["ADT"]] <- JoinLayers(obj[["ADT"]])
obj <- RunUMAP(obj, reduction = "integrated.mnn.adt", dims = 1:DIMS_ADT, seed.use = SEED,
               reduction.name = "umap.mnn.adt", verbose = FALSE)

################################################################################
# 5. WNN, UMAP, clustering sweep
################################################################################
cat("-- WNN + clustering --\n")
obj <- FindMultiModalNeighbors(obj,
         reduction.list = list("integrated.mnn.rna", "integrated.mnn.adt"),
         dims.list = list(1:DIMS_RNA, 1:DIMS_ADT), verbose = FALSE)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", seed.use = SEED, verbose = FALSE)

for (res in seq(0.2, 2.0, by = 0.2)) {
  obj <- FindClusters(obj, graph.name = "wsnn", resolution = res,
                      algorithm = 3, random.seed = SEED, verbose = FALSE)
}
for (col in grep("^wsnn_res\\.", colnames(obj@meta.data), value = TRUE)) {
  l <- as.character(sort(as.numeric(levels(obj[[col]][, 1]))))
  obj[[col]] <- factor(obj[[col]][, 1], levels = l)
}

p_ct <- clustree(obj, prefix = "wsnn_res.") +
  theme(plot.margin = margin(10, 10, 22, 10))
ggsave(file.path(D_CLU, paste0(LIN, "_clustree.png")), p_ct,
       width = 14, height = 13, dpi = 350, bg = "white")

res_tab <- do.call(rbind, lapply(
  grep("^wsnn_res\\.", colnames(obj@meta.data), value = TRUE), function(cc) {
    tb <- table(obj@meta.data[[cc]])
    data.frame(resolution = as.numeric(sub("wsnn_res\\.", "", cc)),
               n_clusters = length(tb), min_cluster = min(tb),
               median_cluster = median(tb), n_below_50 = sum(tb < 50))
  }))
write.csv(res_tab, file.path(D_CLU, paste0(LIN, "_resolution_sweep.csv")),
          row.names = FALSE)
cat("\nresolution sweep:\n"); print(res_tab, row.names = FALSE)

# Provisional resolution: 0.4, matching the published CD8 analysis. Script 7
# revisits it against the marker profiles before annotation.
CHOSEN <- 0.4
obj$seurat_clusters <- obj[[paste0("wsnn_res.", CHOSEN)]][, 1]
Idents(obj) <- "seurat_clusters"
cat("\nprovisional resolution", CHOSEN, "->",
    nlevels(obj$seurat_clusters), "clusters\n")

p_um <- DimPlot(obj, reduction = "wnn.umap", group.by = "seurat_clusters",
                label = TRUE, raster = FALSE) +
  NoLegend() + ggtitle(paste0(LIN, " WNN UMAP, res ", CHOSEN))
ggsave(file.path(D_CLU, paste0(LIN, "_wnn_umap_res", CHOSEN, ".png")), p_um,
       width = 8, height = 7.5, dpi = 400, bg = "white")

################################################################################
# 6. Checkpoint
################################################################################
save_qs(obj, CKPT,
        note = paste0(LIN, " subset re-integrated (", INT_LABEL,
                      ", RNA+ADT) and WNN clustered 0.2-2.0; provisional res ",
                      CHOSEN))

cat("\n", strrep("=", 70), "\n")
cat("  ", LIN, "done:", format(ncol(obj), big.mark = ","), "cells,",
    nlevels(obj$seurat_clusters), "clusters at res", CHOSEN, "\n")
cat("  diagnostics:", OUT, "\n")
cat("  next: Rscript \"Scripts/4. Revision/7.Lineage_Annotation.R\" ", LIN, "\n", sep = "")
cat(strrep("=", 70), "\n")
