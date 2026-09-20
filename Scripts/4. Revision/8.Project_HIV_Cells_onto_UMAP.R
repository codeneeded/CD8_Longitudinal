################################################################################
# Revision script 8 — place the QC-excluded HIV+ cells onto the PUBLISHED UMAP
#
# ANSWERS  Rev 1 major 3, follow-up: where do the infected cells sit relative to
#          the cohort, including the ones our QC removed?
#
# WHY NOT ProjectUMAP / MapQuery
#   Neither works here. Every UMAP in TARA_ALL_annotated_final.qs2 was built
#   without return.model=TRUE, so no uwot model is stored (verified: misc is
#   empty for all 7 UMAP reductions). ProjectUMAP and MapQuery's UMAP projection
#   both require that model. Re-running RunUMAP(return.model=TRUE) would produce
#   a DIFFERENT embedding and invalidate every published coordinate.
#
# WHAT THIS DOES INSTEAD — k-nearest-neighbour placement
#   The published coordinates are left untouched. Each query cell is placed at
#   the similarity-weighted centroid of its k nearest transcriptional neighbours
#   among reference cells, using the neighbours' existing UMAP coordinates.
#
#   The key design choice: candidate neighbours are restricted to cells from THE
#   SAME LIBRARY (CP006_1m / CP006_entry). The query cells come from a different
#   alignment (Cell Ranger 10.1.0 + HIV reference) than the published object, so
#   a cohort-wide neighbour search would confound alignment and batch effects.
#   Restricting to the same library matches batch by construction and removes
#   that confound entirely - there are 6,233 same-library reference cells, which
#   is ample.
#
# WHAT IS AND IS NOT CLAIMED
#   A projected point is NOT a cell that was integrated. It is a statement that
#   this cell's transcriptome most closely resembles reference cells sitting
#   HERE. The figure marks projected cells distinctly and the legend says so.
#
# VALIDATION (this is the point of the script, not the plot)
#   1. Annotation consensus - the neighbours' published labels must agree with
#      the INDEPENDENT Cell Ranger reference annotation of the query cell. Two
#      unrelated methods agreeing is real evidence; a dot on a plot is not.
#   2. Similarity floor - a query whose best neighbours are no more similar than
#      a random reference pair is unplaceable and is reported as such, not drawn.
#   3. Leave-one-out control - reference cells of known annotation are projected
#      by the same procedure with themselves excluded, and we report how often
#      the recovered label matches the truth. That number bounds how much the
#      placement of the query cells can be trusted.
#
# Input : QS_ALL (published) + the pilot filtered matrix
# Output: Reviewer_Comments/Figures/R1M3_04_HIV_projected_UMAP.png + source CSVs
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

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(Matrix); library(ggplot2); library(dplyr)
})

PILOT <- Sys.getenv("PILOT_OUTS",
  path.expand("~/Desktop/Projects/cellranger_out/pilot_hiv_CP006_1m/outs"))
LIB_IN_PUB <- "CP006_entry"     # how this library is named in the published object
K          <- 30                # neighbours used for placement
MIN_VIRAL  <- 3                 # viral UMI floor for "infected" (dose-response cut)

################################################################################
# 1. Load the published reference and restrict to the same library
################################################################################
ref <- qs_read(QS_ALL)
stopifnot("wnn.umap" %in% Reductions(ref))

# Record, rather than assume, that no uwot model exists
.has_model <- vapply(grep("umap", Reductions(ref), value = TRUE, ignore.case = TRUE),
                     function(r) length(slot(ref[[r]], "misc")) > 0, logical(1))
if (any(.has_model)) {
  message("NOTE: a UMAP model IS stored (", paste(names(.has_model)[.has_model],
          collapse = ", "), "). ProjectUMAP would be preferable to kNN placement.")
} else {
  message("Confirmed: no uwot model stored in any UMAP - using kNN placement.")
}

same_lib <- grepl(LIB_IN_PUB, colnames(ref))
if (!any(same_lib))
  stop("No reference cells match '", LIB_IN_PUB, "'.", call. = FALSE)
cat("\nReference cells from", LIB_IN_PUB, ":", sum(same_lib), "\n")

################################################################################
# 2. Load the pilot matrix and find the viral cells
################################################################################
mat  <- Read10X(file.path(PILOT, "filtered_feature_bc_matrix"))
if (is.list(mat)) mat <- mat[["Gene Expression"]]
vgen <- grep("^HIV1_", rownames(mat), value = TRUE)
stopifnot(length(vgen) == 2)

vumi <- Matrix::colSums(mat[vgen, , drop = FALSE])
qbc  <- names(vumi)[vumi >= MIN_VIRAL]
cat("Pilot cells with >=", MIN_VIRAL, "viral UMI:", length(qbc), "\n")

ct <- read.csv(file.path(PILOT, "cell_types", "Azimuth", "cell_types.csv"),
               stringsAsFactors = FALSE)
bc_col <- intersect(c("barcode", "Barcode"), names(ct))[1]
cr_lab <- setNames(ct$fine_cell_type, ct[[bc_col]])

################################################################################
# 3. Shared feature space, normalised independently on each side
#    The two matrices come from different Cell Ranger versions, so the published
#    normalisation must not be assumed to transfer.
################################################################################
ref_counts <- LayerData(ref, assay = "RNA", layer = "counts")[, same_lib, drop = FALSE]
shared <- setdiff(intersect(rownames(ref_counts), rownames(mat)), vgen)
cat("Shared features:", length(shared), "\n")

use <- intersect(VariableFeatures(ref), shared)
if (length(use) < 500) use <- shared
cat("Features used for neighbour search:", length(use), "\n")

lognorm <- function(m) {
  cs <- Matrix::colSums(m); cs[cs == 0] <- 1
  m <- m %*% Matrix::Diagonal(x = 1e4 / cs)
  m@x <- log1p(m@x)
  m
}
R <- lognorm(ref_counts[use, , drop = FALSE])
Q <- lognorm(mat[use, qbc, drop = FALSE])

mu  <- Matrix::rowMeans(R)
sdv <- sqrt(Matrix::rowMeans((R - mu)^2)); sdv[sdv == 0] <- 1
zs  <- function(m) { m <- as.matrix(m); m <- (m - mu) / sdv
                     sweep(m, 2, pmax(sqrt(colSums(m^2)), 1e-9), "/") }
Rz <- zs(R); Qz <- zs(Q)

################################################################################
# 4. kNN placement on the published wnn.umap coordinates
################################################################################
emb <- Embeddings(ref, "wnn.umap")[same_lib, , drop = FALSE]
ann <- as.character(ref@meta.data[same_lib, "Annotation"])

sim <- t(Rz) %*% Qz
set.seed(1)
i1 <- sample(ncol(Rz), 4000, replace = TRUE); i2 <- sample(ncol(Rz), 4000, replace = TRUE)
keep <- i1 != i2
null_sim <- colSums(Rz[, i1[keep]] * Rz[, i2[keep]])
SIM_FLOOR <- as.numeric(quantile(null_sim, 0.99))
cat(sprintf("\nNull similarity (random reference pairs): median %.3f, 99th pct %.3f\n",
            median(null_sim), SIM_FLOOR))

place <- function(j) {
  s <- sim[, j]
  o <- order(s, decreasing = TRUE)[seq_len(K)]
  w <- pmax(s[o], 0); if (sum(w) == 0) w <- rep(1, length(o))
  tb <- sort(tapply(w, ann[o], sum), decreasing = TRUE)
  list(x = sum(w * emb[o, 1]) / sum(w), y = sum(w * emb[o, 2]) / sum(w),
       top_sim = max(s), mean_sim_k = mean(s[o]),
       consensus = names(tb)[1], consensus_frac = as.numeric(tb[1] / sum(tb)),
       placeable = max(s) > SIM_FLOOR)
}
res <- lapply(seq_along(qbc), place)

pub_core <- sub(paste0("^.*", LIB_IN_PUB, "_"), "", colnames(ref)[same_lib])
proj <- data.frame(
  barcode          = qbc,
  viral_UMI        = as.integer(vumi[qbc]),
  cellranger_label = unname(cr_lab[qbc]),
  knn_consensus    = vapply(res, `[[`, character(1), "consensus"),
  consensus_frac   = round(vapply(res, `[[`, numeric(1), "consensus_frac"), 3),
  top_similarity   = round(vapply(res, `[[`, numeric(1), "top_sim"), 3),
  mean_sim_k       = round(vapply(res, `[[`, numeric(1), "mean_sim_k"), 3),
  placeable        = vapply(res, `[[`, logical(1), "placeable"),
  umap_1           = round(vapply(res, `[[`, numeric(1), "x"), 3),
  umap_2           = round(vapply(res, `[[`, numeric(1), "y"), 3),
  in_published     = qbc %in% pub_core,
  stringsAsFactors = FALSE
)
proj <- proj[order(-proj$viral_UMI), ]

################################################################################
# 5. VALIDATION 1 — do two independent annotations agree?
################################################################################
lineage <- function(x) ifelse(grepl("CD8|TEMRA|Tex|Tscm|CTL", x), "CD8",
                       ifelse(grepl("CD4|Treg|Th1|Th2|Th17|Tfh|Naive [0-9]", x), "CD4",
                       ifelse(grepl("NK", x), "NK",
                       ifelse(grepl("B cell|Plasma", x), "B",
                       ifelse(grepl("Mono|DC|Macro", x), "Myeloid", "other")))))
proj$cr_lineage    <- lineage(proj$cellranger_label)
proj$knn_lineage   <- lineage(proj$knn_consensus)
proj$lineage_agree <- proj$cr_lineage == proj$knn_lineage

################################################################################
# 6. VALIDATION 3 — leave-one-out control, the ceiling on trust
################################################################################
set.seed(42)
ctrl_idx <- sample(ncol(Rz), min(300, ncol(Rz)))
sim_rr <- t(Rz) %*% Rz[, ctrl_idx, drop = FALSE]
loo_ok <- vapply(seq_along(ctrl_idx), function(j) {
  s <- sim_rr[, j]; s[ctrl_idx[j]] <- -Inf
  o <- order(s, decreasing = TRUE)[seq_len(K)]
  w <- pmax(s[o], 0); if (sum(w) == 0) w <- rep(1, length(o))
  tb <- sort(tapply(w, ann[o], sum), decreasing = TRUE)
  lineage(names(tb)[1]) == lineage(ann[ctrl_idx[j]])
}, logical(1))
LOO <- mean(loo_ok)
cat(sprintf("\nLeave-one-out control: correct LINEAGE recovered for %.1f%% of %d cells\n",
            100 * LOO, length(ctrl_idx)))
cat("  That is the ceiling on how far any single projected point can be trusted.\n")

save_tab("HIV_projection_validation", data.frame(
  metric = c("k", "reference_cells_same_library", "features_used",
             "null_similarity_99pct", "leave_one_out_lineage_accuracy",
             "query_cells", "query_placeable", "lineage_agreement"),
  value  = c(K, sum(same_lib), length(use), round(SIM_FLOOR, 3), round(LOO, 3),
             nrow(proj), sum(proj$placeable), round(mean(proj$lineage_agree), 3))))
save_tab("HIV_projected_cells", proj)

################################################################################
# 7. Figure
################################################################################
# All-grey reference makes the panel unreadable: the reader cannot tell WHICH
# compartment the projected cells land in, which is the entire point. Colour the
# reference by broad lineage instead, and label each lineage at its centroid.
ref_df <- data.frame(umap_1 = emb[, 1], umap_2 = emb[, 2],
                     Annotation = ann, Lineage = lineage(ann))
pl <- proj[proj$placeable, ]

lin_cols <- c(CD4 = "#4C9F70", CD8 = "#7FA8D9", NK = "#C9A227",
              B = "#B07AA1", Myeloid = "#E08B4F", other = "grey75")
lab_df <- aggregate(cbind(umap_1, umap_2) ~ Lineage, data = ref_df, FUN = median)

# Plain red dots. The earlier version encoded viral burden in point SIZE and
# QC status in point SHAPE (circle vs triangle), which needed two extra legends
# and made the panel busy without adding anything the table does not carry
# better. One mark per cell, one size, one colour; the per-cell burden and
# identity live in HIV_projected_cells.csv and in the accompanying table.
p <- ggplot(ref_df, aes(umap_1, umap_2)) +
  geom_point(aes(colour = Lineage), size = 0.7, alpha = 0.55, stroke = 0) +
  scale_colour_manual(values = lin_cols, name = "reference lineage",
                      guide = guide_legend(override.aes = list(size = 3.4,
                                                               alpha = 1))) +
  geom_text(data = lab_df, aes(label = Lineage), size = 4.1, fontface = "bold",
            colour = "grey15") +
  geom_point(data = pl, colour = "#C00000", size = 3.4, alpha = 0.95) +
  labs(title = "HIV-1 transcript-positive cells on the published UMAP",
       subtitle = paste(strwrap(paste0(
                         "Red dots mark the ", nrow(pl),
                         " HIV-1 transcript-positive cells; ",
                         sum(pl$knn_lineage == "CD4"),
                         " of them fall in the CD4 compartment. Positions are ",
                         "inferred by ", K, "-nearest-neighbour placement against ",
                         sum(same_lib), " same-library reference cells, not by ",
                         "integration, and the published embedding is unchanged. ",
                         "Leave-one-out lineage accuracy of the placement: ",
                         sprintf("%.0f%%", 100 * LOO), ". ",
                         nrow(proj) - nrow(pl),
                         " cell below the similarity floor is not drawn. Per-cell ",
                         "viral burden and identity are in the accompanying table."),
                         width = 96), collapse = "\n"),
       x = "wnnUMAP 1", y = "wnnUMAP 2") +
  theme_rev() + theme(legend.position = "right",
                      legend.box.margin = margin(0, 6, 0, 10),
                      plot.margin = margin(12, 26, 12, 12))

save_fig("R1M3_04_HIV_projected_UMAP", p, width = 11.5, height = 7.4)

cat("\n=== projected cells ===\n")
print(proj[, c("viral_UMI","cellranger_label","knn_consensus","consensus_frac",
               "top_similarity","lineage_agree","in_published","placeable")],
      row.names = FALSE)
