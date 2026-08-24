################################################################################
# CP003 + CP018 — joint integration of the two FACT-Seq datasets
#
# SCOPE: the two STIMULATED (FACT-Seq) objects only. Both are HIV-peptide
# stimulated and therefore directly comparable. The unstimulated eCITE-seq data
# is deliberately excluded: it differs in assay content (carries ADT), in
# stimulation state, and in library chemistry, so folding it in would confound
# participant with experiment.
#
# THE CONFOUND, STATED UP FRONT: n = 1 participant per suppression group. Any
# axis separating CP003 from CP018 is inseparable from individual identity, and
# the two were aligned with different Cell Ranger versions. Integration cannot
# fix this. Its purpose is a SHARED coordinate system for visualisation and for
# testing whether the Tex asymmetry survives batch correction - not to create
# statistical power that the design does not have.
#
# Section 3 is the load-bearing check: batch correction can erase a genuine
# population difference by forcing clusters to mix. Pseudobulk correlation
# (script 3b) does not require that assumption, so it stays the primary
# evidence; integration is corroboration.
#
# Input : QS_ANNOTATED (CP018, script 5), QS_CP003_ANNOT_FIG3 (CP003)
# Output: INT_JOINT — UMAPs, composition, Tex-survival check, source CSVs
#         checkpoint: saved_R_data/CP003_CP018_integrated.qs2
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(scales); library(patchwork)
})

INT_JOINT <- file.path(CP018_ANALYSIS, "CP003_vs_CP018", "Integrated"); mkdirs(INT_JOINT)
csv_j     <- file.path(INT_JOINT, "tables");                            mkdirs(csv_j)
QS_JOINT  <- file.path(SAVED_DIR, "CP003_CP018_integrated.qs2")

require_file(QS_ANNOTATED,        "CP018 annotated stim object")
require_file(QS_CP003_ANNOT_FIG3, "CP003 Fig3-annotated object")

################################################################################
# 1. Load, harmonise metadata, merge
################################################################################
cp18 <- qs_read(QS_ANNOTATED)
cp03 <- qs_read(QS_CP003_ANNOT_FIG3)
DefaultAssay(cp18) <- "RNA"; DefaultAssay(cp03) <- "RNA"

# Keep RNA only. The prediction-score and mnn.reconstructed assays are
# per-dataset artefacts and must not travel into a joint object.
cp18 <- DietSeurat(cp18, assays = "RNA")
cp03 <- DietSeurat(cp03, assays = "RNA")

# Each object still carries the per-library split layers from its OWN
# integration (CP018 has counts.CP018_2m_A etc.). Those must be joined before
# the joint object can be re-split by participant, or split() refuses.
join_if_split <- function(o, lab) {
  n <- length(Layers(o[["RNA"]]))
  if (n > 3) { o <- JoinLayers(o, assay = "RNA")
               cat("  joined", lab, "-", n, "->", length(Layers(o[["RNA"]])), "layers\n") }
  o
}
cat("Layer state before merge:\n")
cp18 <- join_if_split(cp18, "CP018")
cp03 <- join_if_split(cp03, "CP003")

cp18$participant <- "CP018"; cp03$participant <- "CP003"
cp18$ART_group   <- "Suppressed"; cp03$ART_group <- "Unsuppressed"
cp18$cellranger  <- "10.1.0";     cp03$cellranger <- "9.0.1"

# One annotation column spanning both, using each dataset's own published labels.
cp18$joint_annotation <- as.character(cp18$Fig_Annotation)
cp03$joint_annotation <- as.character(cp03$Fig3_Annotation)
cp18$joint_timepoint  <- paste0("CP018_", as.character(cp18$Timepoint))
cp03$joint_timepoint  <- paste0("CP003_", as.character(cp03$Timepoint))

shared <- intersect(rownames(cp18), rownames(cp03))
cat("Shared genes:", length(shared), "\n")
stopifnot(length(shared) > 5000)

comb <- merge(cp18[shared, ], cp03[shared, ],
              add.cell.ids = c("CP018", "CP003"))
rm(cp18, cp03); gc(verbose = FALSE)
cat("Merged:", ncol(comb), "cells |", nrow(comb), "genes\n")
print(table(comb$participant, comb$ART_group))

################################################################################
# 2. Integrate with participant as batch
################################################################################
# Same method as the within-participant pipeline (FastMNN) so cluster structure
# is comparable to the per-participant results rather than method-dependent.
# merge() already leaves one layer PER OBJECT (counts.1/counts.2), which IS the
# participant split - re-splitting is both unnecessary and an error. Join first
# only if the layer count does not match the number of participants.
n_lay  <- length(grep("^counts", Layers(comb[["RNA"]]), value = TRUE))
n_part <- n_distinct(comb$participant)
cat("post-merge count layers:", n_lay, "| participants:", n_part, "\n")
if (n_lay != n_part) {
  comb <- JoinLayers(comb, assay = "RNA")
  comb[["RNA"]] <- split(comb[["RNA"]], f = comb$participant)
  cat("re-split by participant ->", length(Layers(comb[["RNA"]])), "layers\n")
}

comb <- NormalizeData(comb, verbose = FALSE)
comb <- FindVariableFeatures(comb, nfeatures = 2000, verbose = FALSE)
comb <- ScaleData(comb, verbose = FALSE)
comb <- RunPCA(comb, npcs = 30, verbose = FALSE)

# CCAIntegration ships with Seurat itself (no extra dependency) and script 3
# already uses it for the within-participant integration, so cluster structure
# stays comparable rather than method-dependent. FastMNN would need batchelor,
# which is not installed here.
# FastMNNIntegration is exported by SeuratWrappers, NOT by batchelor - checking
# only for batchelor selects the FastMNN branch and then fails on an unresolved
# symbol, which is what happened on the analysis machine. Require BOTH, and
# resolve the function object explicitly rather than by bare name.
.fastmnn_ok <- requireNamespace("batchelor", quietly = TRUE) &&
               requireNamespace("SeuratWrappers", quietly = TRUE) &&
               exists("FastMNNIntegration",
                      where = asNamespace("SeuratWrappers"), inherits = FALSE)
INT_METHOD <- if (.fastmnn_ok) "FastMNN" else "CCA"
cat("Integration method:", INT_METHOD,
    if (!.fastmnn_ok) "(FastMNN unavailable: needs batchelor + SeuratWrappers)" else "",
    "\n")
comb <- if (INT_METHOD == "FastMNN") {
  IntegrateLayers(comb, method = get("FastMNNIntegration",
                                     envir = asNamespace("SeuratWrappers")),
                  assay = "RNA",
                  orig.reduction = "pca", new.reduction = "integrated",
                  verbose = FALSE)
} else {
  IntegrateLayers(comb, method = CCAIntegration, assay = "RNA",
                  orig.reduction = "pca", new.reduction = "integrated",
                  verbose = FALSE)
}
comb <- JoinLayers(comb)
comb$integration_method <- INT_METHOD

comb <- FindNeighbors(comb, reduction = "integrated", dims = 1:30, verbose = FALSE)
comb <- FindClusters(comb, resolution = MNN_RES_FINAL, cluster.name = "joint_clusters",
                     verbose = FALSE)
comb <- RunUMAP(comb, reduction = "integrated", dims = 1:30,
                reduction.name = "umap.joint", verbose = FALSE)

# Newly created object - stamp the canonical clonotype column before saving so
# anything reading this checkpoint gets a consistent identifier.
if ("CTstrict" %in% colnames(comb@meta.data)) comb <- add_ct_canon(comb)
qs_save(comb, file = QS_JOINT)
message("Saved checkpoint: ", QS_JOINT)

################################################################################
# 3. Did integration erase the biology, or preserve it?
################################################################################
# A cluster containing cells from only ONE participant is either a real
# participant-specific population or a batch-correction failure. Quantify rather
# than eyeball: for every joint cluster, report the participant mix.
mix <- comb@meta.data %>%
  count(joint_clusters, participant, name = "cells") %>%
  group_by(joint_clusters) %>%
  mutate(pct = round(100 * cells / sum(cells), 1),
         cluster_total = sum(cells)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(joint_clusters, cluster_total),
              names_from = participant, values_from = pct, values_fill = 0) %>%
  mutate(dominance = pmax(CP003, CP018),
         verdict = case_when(
           dominance >= 95 ~ "participant-EXCLUSIVE",
           dominance >= 80 ~ "participant-skewed",
           TRUE            ~ "mixed"))
save_csv(mix, file.path(csv_j, "joint_cluster_participant_mix.csv"))
cat("\n=== Participant composition of each joint cluster ===\n")
print(as.data.frame(mix), row.names = FALSE)

# Where do CP003's Tex cells land? If Tex is genuinely absent from CP018, the
# cluster containing them should be strongly CP003-dominated. If integration
# smoothed it away, Tex cells will be spread across mixed clusters.
if (any(comb$joint_annotation == "Tex", na.rm = TRUE)) {
  tex_cells <- comb$joint_annotation == "Tex" & !is.na(comb$joint_annotation)
  tex_home <- comb@meta.data[tex_cells, ] %>%
    count(joint_clusters, name = "tex_cells") %>%
    left_join(mix %>% select(joint_clusters, cluster_total, CP003, CP018, verdict),
              by = "joint_clusters") %>%
    mutate(tex_pct_of_cluster = round(100 * tex_cells / cluster_total, 1)) %>%
    arrange(desc(tex_cells))
  save_csv(tex_home, file.path(csv_j, "CP003_Tex_destination_after_integration.csv"))
  cat("\n=== Where CP003's Tex cells sit after integration ===\n")
  print(as.data.frame(tex_home), row.names = FALSE)

  top <- tex_home[1, ]
  cat("\nLargest Tex-containing cluster:", as.character(top$joint_clusters),
      "| CP003", top$CP003, "% / CP018", top$CP018, "% |", top$verdict, "\n")
  if (top$CP018 > 20) {
    cat("CAUTION: that cluster is >20% CP018, so CP018 cells DO occupy the\n",
        "region CP003's Tex cells map to. Either the population is not fully\n",
        "absent in CP018, or integration has over-corrected. The pseudobulk\n",
        "correlation (script 3b, Tex r = 0.695 vs >=0.87 elsewhere) does not\n",
        "depend on batch correction and remains the primary evidence.\n", sep = "")
  } else {
    cat("Tex maps to a CP003-dominated cluster, consistent with the pseudobulk\n",
        "correlation result. Integration corroborates the asymmetry.\n", sep = "")
  }
}

################################################################################
# 3b. Is a "mixed" Tex cluster over-correction, or a real rare subpopulation?
################################################################################
# A Tex cluster containing both participants has two readings, and they differ
# in what we would report:
#   (a) integration over-corrected and forced CP018 cells into a CP003 region
#   (b) CP018 genuinely has a small exhausted subpopulation that whole-cluster
#       pseudobulk correlation averaged away
# These are separable: under (a) the CP018 cells in that cluster look like their
# own participant's neighbours; under (b) they look exhausted on their own terms.
if (exists("tex_home") && nrow(tex_home) > 0) {
  tex_cl <- as.character(tex_home$joint_clusters[1])
  DefaultAssay(comb) <- "RNA"
  if (length(Layers(comb[["RNA"]])) > 3) comb <- JoinLayers(comb, assay = "RNA")
  Xj <- GetAssayData(comb, layer = "data")

  EX <- intersect(c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1"), rownames(Xj))
  ST <- intersect(c("TCF7","IL7R","SELL","CCR7","ID3","LEF1","BACH2"),        rownames(Xj))
  in_cl <- comb$joint_clusters == tex_cl

  score_grp <- function(mask, label) {
    if (sum(mask) < 5) return(NULL)
    ex <- mean(Matrix::colMeans(Xj[EX, mask, drop = FALSE]))
    st <- mean(Matrix::colMeans(Xj[ST, mask, drop = FALSE]))
    data.frame(group = label, n_cells = sum(mask),
               exhaustion = round(ex, 3), stemness = round(st, 3),
               exh_stem_ratio = round(ex / (st + 0.01), 2),
               stringsAsFactors = FALSE)
  }
  # The decisive contrast: CP018 cells INSIDE the Tex cluster vs CP018 cells
  # carrying the SAME published label elsewhere.
  lab_here <- comb$joint_annotation
  cmp <- do.call(rbind, list(
    score_grp(in_cl & comb$participant == "CP003", paste0("CP003 in cluster ", tex_cl)),
    score_grp(in_cl & comb$participant == "CP018", paste0("CP018 in cluster ", tex_cl)),
    score_grp(!in_cl & comb$participant == "CP018" & lab_here == "Cycling Effector",
              "CP018 Cycling Effector elsewhere"),
    score_grp(!in_cl, "all cells outside the cluster")))
  save_csv(cmp, file.path(csv_j, "Tex_cluster_exhaustion_check.csv"))
  cat("\n=== Is the shared Tex cluster over-correction or real? ===\n")
  print(cmp, row.names = FALSE)

  gene_tab <- data.frame(
    gene  = c(EX, ST),
    panel = c(rep("exhaustion", length(EX)), rep("stemness", length(ST))),
    CP003_in_cluster = round(Matrix::rowMeans(
      Xj[c(EX, ST), in_cl & comb$participant == "CP003", drop = FALSE]), 3),
    CP018_in_cluster = round(Matrix::rowMeans(
      Xj[c(EX, ST), in_cl & comb$participant == "CP018", drop = FALSE]), 3),
    row.names = NULL)
  save_csv(gene_tab, file.path(csv_j, "Tex_cluster_per_gene_by_participant.csv"))
  cat("\n=== Per-gene inside the Tex cluster ===\n")
  print(gene_tab, row.names = FALSE)

  # Prevalence, normalised within participant - the difference is magnitude.
  prev <- comb@meta.data %>% group_by(participant) %>%
    summarise(cells_in_tex_cluster = sum(joint_clusters == tex_cl),
              participant_total = n(),
              pct_of_participant = round(100 * sum(joint_clusters == tex_cl) / n(), 2),
              .groups = "drop")
  save_csv(prev, file.path(csv_j, "Tex_cluster_prevalence.csv"))
  cat("\n=== Prevalence of the exhausted cluster, within participant ===\n")
  print(as.data.frame(prev), row.names = FALSE)

  cp18_ratio <- cmp$exh_stem_ratio[grepl("CP018 in cluster", cmp$group)]
  cp18_else  <- cmp$exh_stem_ratio[grepl("elsewhere", cmp$group)]
  if (length(cp18_ratio) && length(cp18_else) && cp18_ratio > 1.5 * cp18_else) {
    cat("\nVERDICT: the CP018 cells in this cluster are exhausted ON THEIR OWN\n",
        "TERMS (ratio ", cp18_ratio, " vs ", cp18_else, " for the same label\n",
        "elsewhere), so this is a REAL rare exhausted subpopulation, not\n",
        "over-correction. Report the CP003/CP018 difference as one of\n",
        "MAGNITUDE, not absence. Whole-cluster pseudobulk correlation\n",
        "averages a subpopulation this small away.\n", sep = "")
  } else {
    cat("\nVERDICT: CP018 cells here resemble their own participant's cells of\n",
        "the same label, consistent with batch over-correction. Treat the\n",
        "shared cluster as an integration artefact and rely on the pseudobulk\n",
        "correlation.\n", sep = "")
  }

  p_ex <- ggplot(cmp, aes(x = reorder(group, exh_stem_ratio), y = exh_stem_ratio,
                          fill = exh_stem_ratio > 1.5)) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.68) +
    geom_text(aes(label = paste0(exh_stem_ratio, "  (n=", n_cells, ")")),
              hjust = -0.12, size = 3.6) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    coord_flip() +
    labs(x = NULL, y = "Exhaustion : stemness ratio",
         title = paste0("Exhausted phenotype inside joint cluster ", tex_cl),
         subtitle = "Ratio > 1.5 = inhibitory receptors high with stemness collapsed") +
    theme_cp018(base_size = 13) +
    theme(plot.title.position = "plot",
          plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(hjust = 0))
  save_png(file.path(INT_JOINT, "Tex_cluster_exhaustion_check.png"), p_ex,
           width = 10, height = 4.6)
}

# Reciprocal view: which joint cluster does each published label occupy?
lab_map <- comb@meta.data %>%
  filter(!is.na(joint_annotation)) %>%
  count(participant, joint_annotation, joint_clusters, name = "cells") %>%
  group_by(participant, joint_annotation) %>%
  slice_max(cells, n = 1, with_ties = FALSE) %>% ungroup() %>%
  arrange(joint_annotation, participant)
save_csv(lab_map, file.path(csv_j, "joint_label_to_cluster_map.csv"))

################################################################################
# 4. Plots
################################################################################
pal_part <- c(CP003 = "#B2182B", CP018 = "#2166AC")

p1 <- tidy_umap(DimPlot(comb, reduction = "umap.joint", group.by = "participant",
                        cols = pal_part, pt.size = 0.3, shuffle = TRUE),
                title = "Joint integration: participant",
                subtitle = paste0(ncol(comb), " cells | ", INT_METHOD,
                                  " integration, participant as batch"))
save_png(file.path(INT_JOINT, "Joint_UMAP_participant.png"), p1, width = 8.5, height = 6.5)

p2 <- tidy_umap(DimPlot(comb, reduction = "umap.joint", group.by = "joint_clusters",
                        label = TRUE, pt.size = 0.3) + NoLegend(),
                title = paste0("Joint clusters (resolution ", MNN_RES_FINAL, ")"))
save_png(file.path(INT_JOINT, "Joint_UMAP_clusters.png"), p2, width = 8, height = 6.5)

p3 <- tidy_umap(DimPlot(comb, reduction = "umap.joint", group.by = "joint_annotation",
                        split.by = "participant", pt.size = 0.35),
                title = "Published labels in the joint space")
save_png(file.path(INT_JOINT, "Joint_UMAP_annotation_by_participant.png"),
         p3, width = 13, height = 6)

# Composition: proportion, not count - the two participants differ in cell number.
comp <- comb@meta.data %>%
  count(participant, joint_clusters, name = "cells") %>%
  group_by(participant) %>% mutate(pct = 100 * cells / sum(cells)) %>% ungroup()
save_csv(comp, file.path(csv_j, "joint_cluster_composition.csv"))

p4 <- ggplot(comp, aes(x = joint_clusters, y = pct, fill = participant)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = pal_part, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Joint cluster", y = "% of participant's cells",
       title = "Cluster composition, normalised within participant",
       subtitle = "Proportions, not counts: CP018 has 5,129 cells vs CP003's 4,032") +
  theme_cp018(base_size = 13) +
  theme(plot.title.position = "plot",
        plot.title = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0))
save_png(file.path(INT_JOINT, "Joint_cluster_composition.png"), p4, width = 10, height = 5.5)

cat("\nOutputs in:\n  ", INT_JOINT, "\n")
cat("\nREMINDER: n = 1 participant per group. Integration provides a shared\n",
    "coordinate system, NOT statistical power. Any CP003/CP018 difference\n",
    "remains confounded with individual identity and Cell Ranger version.\n", sep = "")
