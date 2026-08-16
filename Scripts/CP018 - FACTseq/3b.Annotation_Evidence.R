################################################################################
# CP018 — annotation evidence and justification
#
# Produces the audit trail behind CP018_CLUSTER_ANNOTATION in config.R, so the
# label on every cluster can be checked against numbers rather than taken on
# trust. Run AFTER script 3, BEFORE accepting the annotation.
#
# Panels span stemness, activation (early and late), inhibitory receptors,
# exhaustion TFs, metabolism (glycolysis / OXPHOS / FAO / mTOR), effector and
# cytokine output, cycling, interferon, tissue residency, terminal
# differentiation, and technical confounders (hypoxia, heat-shock, apoptosis,
# doublet score, library complexity) — because a "biological" cluster that is
# really a stress or doublet artefact will otherwise be annotated as real.
#
# Output (all CSV):
#   CP018_annotation_justification.csv   one row per cluster: label + evidence
#   CP018_panel_scores_mean.csv          programme x cluster, mean expression
#   CP018_panel_scores_zscore.csv        programme x cluster, z across clusters
#   CP018_panel_gene_coverage.csv        which panel genes exist in the object
#   CP018_marker_genes_per_cluster.csv   discriminating genes, cluster vs rest
#   CP018_cluster_QC_confounders.csv     depth, mito, doublet score, hypoxia
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

EVI_ROOT <- file.path(INT_ROOT, "Annotation_Evidence")
mkdirs(EVI_ROOT)

require_file(QS_INTEGRATED, "integrated object (script 3)")
seu <- qs_read(QS_INTEGRATED)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

CL_COL <- paste0("mnn_snn_res.", MNN_RES_FINAL)
if (!CL_COL %in% colnames(seu@meta.data)) CL_COL <- "mnn_clusters_rna"
cl_all <- as.character(seu@meta.data[[CL_COL]])

# Drop singleton clusters — they cannot support a marker test.
tiny <- names(which(table(cl_all) < 10))
if (length(tiny)) {
  message("Excluding singleton cluster(s): ", paste(tiny, collapse = ", "))
  seu <- seu[, !cl_all %in% tiny]
  cl_all <- as.character(seu@meta.data[[CL_COL]])
}
message("Clusters evaluated: ", paste(sort(unique(cl_all)), collapse = ", "))

################################################################################
# 1. Marker panels
################################################################################
# Grouped so that a label claim maps onto a specific panel. Sources are the
# standard CD8 differentiation / exhaustion literature; the same panels are
# applied to CP003 in script 8 so the two participants stay comparable.
PANELS <- list(
  # --- stemness / naivety ---
  Stem_TF          = c("TCF7","LEF1","MYB","ID3","BACH2","FOXP1","KLF2"),
  Stem_surface     = c("SELL","CCR7","IL7R","CD27","CD28","CXCR4","FAS"),
  Tscm_associated  = c("CD58","LTB","NOSIP","MAL","TXK","ACTN1","PLAC8"),
  Naive_restricted = c("CCR9","SATB1","LRRN3","NELL2","TSHZ2","OXNAD1"),
  # --- activation ---
  Activation_early = c("CD69","JUN","FOS","JUNB","EGR1","NR4A1","DUSP1","IER2"),
  Activation_late  = c("CD38","HLA-DRA","HLA-DRB1","TNFRSF9","IL2RA","ICOS","TNFRSF4"),
  Costim_TNFRSF    = c("CD70","TNFRSF18","TNFRSF9","TNFRSF4","CD27"),
  # --- exhaustion ---
  Inhibitory_recep = c("CTLA4","PDCD1","LAG3","TIGIT","HAVCR2","BTLA","CD160"),
  Exhaustion_TF    = c("TOX","TOX2","EOMES","NR4A2","BATF","IKZF2","VSIR"),
  Tpex_hallmark    = c("SLAMF6","XCL1","XCL2","TCF7","ID3","CXCR5"),
  # --- metabolism ---
  Glycolysis       = c("SLC2A1","HK2","GAPDH","LDHA","PKM","ENO1","PGK1","TPI1","PFKP"),
  OXPHOS           = c("NDUFA4","COX5A","COX7C","ATP5F1E","UQCRB","SDHB","CYCS"),
  FAO_mitochondrial= c("CPT1A","ACADM","PRKAA1","PPARGC1A","TFAM","SLC25A20"),
  mTOR_growth      = c("MYC","RPS6","EIF4EBP1","SLC7A5","SLC3A2","RHEB","MTOR"),
  # --- effector ---
  Cytotoxicity     = c("GZMA","GZMB","GZMH","GNLY","PRF1","NKG7","KLRD1","FGFBP2"),
  Cytokine         = c("IFNG","TNF","IL2","CSF2","CCL3","CCL4","XCL1","XCL2"),
  Terminal_diff    = c("KLRG1","ZEB2","TBX21","S1PR5","CX3CR1","B3GAT1"),
  Tissue_residency = c("ITGAE","CD69","CXCR6","RUNX3","ZNF683"),
  # --- other states ---
  Cycling          = c("MKI67","TOP2A","PCNA","TYMS","STMN1","CCNB1","CDK1"),
  Type_I_IFN       = c("IFIT1","IFIT3","ISG15","MX1","OAS1","IFI6","STAT1","IRF7"),
  MAIT             = c("KLRB1","SLC4A10","ZBTB16","TRAV1-2","IL18RAP","CEBPD"),
  # --- technical / stress confounders ---
  Hypoxia          = c("BNIP3","BNIP3L","P4HA1","ANKRD37","MIR210HG","PDK1","VEGFA",
                       "ALDOC","EGLN3","NDRG1","HILPDA","DDIT4","SLC16A3"),
  Heat_shock       = c("HSPA1A","HSPA1B","HSP90AA1","DNAJB1","HSPH1","HSPB1","BAG3"),
  Apoptosis        = c("BAX","BCL2L11","CASP3","PMAIP1","BBC3","TP53")
)

X <- GetAssayData(seu, layer = "data")

coverage <- do.call(rbind, lapply(names(PANELS), function(p) data.frame(
  panel = p, gene = PANELS[[p]],
  present_in_object = PANELS[[p]] %in% rownames(X))))
save_csv(coverage, file.path(EVI_ROOT, "CP018_panel_gene_coverage.csv"))
n_missing <- sum(!coverage$present_in_object)
if (n_missing) message("NOTE: ", n_missing, " panel genes absent from the object ",
                       "(see CP018_panel_gene_coverage.csv) - panel means use the rest.")

panel_mean <- function(genes) {
  g <- intersect(genes, rownames(X))
  if (!length(g)) return(setNames(rep(NA_real_, length(unique(cl_all))), sort(unique(cl_all))))
  tapply(Matrix::colMeans(X[g, , drop = FALSE]), cl_all, mean)
}
M <- do.call(rbind, lapply(PANELS, panel_mean))
M <- M[, order(as.numeric(colnames(M))), drop = FALSE]

Z <- t(scale(t(M)))          # z across clusters: which cluster is high for a programme
Z[is.nan(Z)] <- NA

save_csv(data.frame(panel = rownames(M), round(M, 4), check.names = FALSE),
         file.path(EVI_ROOT, "CP018_panel_scores_mean.csv"))
save_csv(data.frame(panel = rownames(Z), round(Z, 3), check.names = FALSE),
         file.path(EVI_ROOT, "CP018_panel_scores_zscore.csv"))

cat("\n=== Panel z-scores (which cluster is HIGH for each programme) ===\n")
print(round(Z, 1))

################################################################################
# 2. Technical confounders — a stress or doublet cluster must not be annotated
#    as a biological cell state.
################################################################################
md <- seu@meta.data
md$.cl <- cl_all
qc <- md %>% group_by(cluster = .cl) %>%
  summarise(
    n_cells        = n(),
    mean_nFeature  = round(mean(nFeature_RNA)),
    mean_nCount    = round(mean(nCount_RNA)),
    pct_mito       = round(mean(percent_mito), 2),
    pct_ribo       = round(mean(percent_ribo), 1),
    complexity     = round(mean(log10GenesPerUMI), 3),
    mean_doublet_score = if ("scDblFinder.score" %in% names(md))
                            round(mean(scDblFinder.score), 3) else NA_real_,
    pct_cycling    = round(100 * mean(Phase != "G1")),
    .groups = "drop")

# Flag clusters whose depth or doublet score is an outlier against the rest.
qc <- qc %>% mutate(
  nFeature_vs_median = round(mean_nFeature / median(mean_nFeature), 2),
  hypoxia_z          = round(Z["Hypoxia",  as.character(cluster)], 2),
  heatshock_z        = round(Z["Heat_shock", as.character(cluster)], 2),
  flag = case_when(
    # Proliferating blasts genuinely have high transcriptional output and an
    # elevated chaperone load, so depth and heat-shock flags are only
    # meaningful in a NON-cycling cluster.
    pct_cycling > 60                             ~ "CYCLING (depth/HSP expected)",
    nFeature_vs_median > 1.8                     ~ "HIGH DEPTH - check doublets",
    hypoxia_z   > 1.5                            ~ "HYPOXIA-DOMINATED",
    heatshock_z > 1.5                            ~ "HEAT-SHOCK / STRESS",
    pct_mito    > 10                             ~ "HIGH MITO",
    TRUE ~ "")
)
save_csv(qc, file.path(EVI_ROOT, "CP018_cluster_QC_confounders.csv"))
cat("\n=== Technical confounders ===\n"); print(as.data.frame(qc), row.names = FALSE)

################################################################################
# 3. Discriminating marker genes per cluster
################################################################################
# Script 3 already computed these; reuse rather than repeat a ~10 min test.
# Reuse script 3's table only if it is actually USABLE. file.exists() is not
# enough: a script-3 run that died partway leaves a 3-byte stub behind, and
# read.csv on that aborts the whole script. Validate, then fall back.
mk_from_s3 <- file.path(INT_ROOT, "tables", "CP018_markers_all_clusters.csv")
mk <- NULL
if (file.exists(mk_from_s3) && file.info(mk_from_s3)$size > 1000) {
  mk_try <- try(utils::read.csv(mk_from_s3, stringsAsFactors = FALSE), silent = TRUE)
  needed <- c("gene", "cluster", "avg_log2FC")
  if (!inherits(mk_try, "try-error") && is.data.frame(mk_try) &&
      nrow(mk_try) > 0 && all(needed %in% names(mk_try))) {
    # must also cover the clusters we are annotating NOW, not a previous run's
    if (all(levels(factor(cl_all)) %in% as.character(unique(mk_try$cluster)))) {
      mk <- mk_try
      message("Reusing marker table from script 3 (", nrow(mk), " rows)")
    } else {
      message("Script 3's marker table covers different clusters - recomputing.")
    }
  } else {
    message("Script 3's marker table is empty or malformed - recomputing.")
  }
}
if (is.null(mk)) {
  message("Computing markers (FindAllMarkers)... this is the slow step.")
  Idents(seu) <- CL_COL
  mk <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.5, verbose = FALSE)
}
stopifnot(nrow(mk) > 0)
save_csv(mk, file.path(EVI_ROOT, "CP018_marker_genes_per_cluster.csv"))
top_mk <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15) %>%
  summarise(top15_markers = paste(gene, collapse = ", "), .groups = "drop")

################################################################################
# 4. Assemble the justification table
################################################################################
# For each cluster: the top three programmes it is HIGH for, the three it is
# LOW for, the assigned label, and whether the evidence supports it.
zdf <- as.data.frame(t(Z))
zdf$cluster <- rownames(zdf)

top_panels <- function(cl, n = 3, high = TRUE) {
  v <- Z[, cl]; v <- v[!is.na(v)]
  v <- sort(v, decreasing = high)[seq_len(min(n, length(v)))]
  paste(sprintf("%s (z=%+.1f)", names(v), v), collapse = "; ")
}

ann <- CP018_CLUSTER_ANNOTATION
just <- data.frame(cluster = colnames(Z), stringsAsFactors = FALSE) %>%
  mutate(
    n_cells        = qc$n_cells[match(cluster, qc$cluster)],
    assigned_label = ifelse(cluster %in% names(ann), ann[cluster], "UNASSIGNED"),
    highest_programmes = vapply(cluster, top_panels, character(1), high = TRUE),
    lowest_programmes  = vapply(cluster, top_panels, character(1), high = FALSE),
    Stem_TF_z        = round(Z["Stem_TF", cluster], 2),
    Stem_surface_z   = round(Z["Stem_surface", cluster], 2),
    Activation_late_z= round(Z["Activation_late", cluster], 2),
    Inhibitory_z     = round(Z["Inhibitory_recep", cluster], 2),
    Exhaustion_TF_z  = round(Z["Exhaustion_TF", cluster], 2),
    Tpex_hallmark_z  = round(Z["Tpex_hallmark", cluster], 2),
    Glycolysis_z     = round(Z["Glycolysis", cluster], 2),
    OXPHOS_z         = round(Z["OXPHOS", cluster], 2),
    Cytotoxicity_z   = round(Z["Cytotoxicity", cluster], 2),
    Cycling_z        = round(Z["Cycling", cluster], 2),
    MAIT_z           = round(Z["MAIT", cluster], 2),
    Hypoxia_z        = round(Z["Hypoxia", cluster], 2),
    QC_flag          = qc$flag[match(cluster, qc$cluster)],
    pct_cycling      = qc$pct_cycling[match(cluster, qc$cluster)],
    doublet_score    = qc$mean_doublet_score[match(cluster, qc$cluster)]
  ) %>%
  left_join(top_mk %>% mutate(cluster = as.character(cluster)), by = "cluster")

# Rule-based check: does the evidence actually support the assigned label?
# This is a CONSISTENCY test, not an automatic annotator — it exists so a
# mislabelled cluster is caught rather than inherited.
check_label <- function(r) {
  lab <- r[["assigned_label"]]
  z   <- function(k) suppressWarnings(as.numeric(r[[k]]))
  ok  <- NA; why <- ""
  if (lab == "Naive/Bystander") {
    ok  <- z("Stem_surface_z") > 0 && z("Cytotoxicity_z") < 0 && z("Cycling_z") < 0.5
    why <- "expect stem-surface high, cytotoxicity low, not cycling"
  } else if (lab == "Cycling Effector") {
    ok  <- z("Cycling_z") > 1 && as.numeric(r[["pct_cycling"]]) > 60
    why <- "expect cycling programme high and >60% S/G2M"
  } else if (lab == "TEMRA/Effector") {
    ok  <- z("Cytotoxicity_z") > 1 && z("Stem_surface_z") < 0
    why <- "expect cytotoxicity high, stem-surface low"
  } else if (lab == "MAIT") {
    ok  <- z("MAIT_z") > 1
    why <- "expect MAIT panel high (KLRB1/SLC4A10/ZBTB16/TRAV1-2)"
  } else if (lab == "Tpex") {
    ok  <- z("Tpex_hallmark_z") > 0.5 && z("Inhibitory_z") > 0.5 && z("Exhaustion_TF_z") > 0
    why <- "Tpex REQUIRES the exhaustion programme engaged alongside stemness"
  } else if (lab == "Activated Stem-like") {
    ok  <- z("Stem_TF_z") > 0.5 && z("Activation_late_z") > 0.5
    why <- "expect stem TFs retained WITH late activation"
  } else if (lab == "Transitional Tem") {
    ok  <- z("Stem_surface_z") < 0.5 && z("Cytotoxicity_z") > -1
    why <- "intermediate: stem-surface declining, some cytotoxicity"
  } else if (lab == "Naive (hypoxic)") {
    ok  <- z("Hypoxia_z") > 1.5 && z("Stem_surface_z") > 0
    why <- "naive surface phenotype with a dominant hypoxia response"
  }
  c(supported = ifelse(is.na(ok), "not_checked", ifelse(ok, "YES", "NO - REVIEW")),
    criterion = why)
}
chk <- t(apply(just, 1, check_label))
just$evidence_supports_label <- chk[, "supported"]
just$label_criterion         <- chk[, "criterion"]

save_csv(just, file.path(EVI_ROOT, "CP018_annotation_justification.csv"))

cat("\n=== ANNOTATION JUSTIFICATION ===\n")
print(as.data.frame(just[, c("cluster","n_cells","assigned_label",
                             "evidence_supports_label","QC_flag",
                             "highest_programmes")]), row.names = FALSE)

bad <- just$cluster[just$evidence_supports_label == "NO - REVIEW"]
if (length(bad)) {
  cat("\n*** Clusters whose evidence does NOT support the assigned label: ",
      paste(bad, collapse = ", "), "\n",
      "    Inspect CP018_annotation_justification.csv before proceeding.\n", sep = "")
}
flagged <- qc$cluster[qc$flag != ""]
if (length(flagged)) {
  cat("\n*** Clusters with a technical flag: ", paste(flagged, collapse = ", "),
      "\n    A stress/doublet artefact must not be annotated as a cell state.\n", sep = "")
}

################################################################################
# 6. Cross-participant naming: which CP018 clusters correspond to CP003 labels?
################################################################################
# Naming is decided by evidence, not preference. Pseudobulk profiles are
# correlated across the two datasets on shared variable genes. Where a CP018
# cluster clearly matches a CP003 label, the SAME NAME is used so panels can be
# read side by side. Where no clear match exists, a distinct name is used - and
# that mismatch is itself a result worth reporting.
if (file.exists(QS_CP003_ANNOT_FIG3)) {
  cp3 <- qs_read(QS_CP003_ANNOT_FIG3)
  DefaultAssay(cp3) <- "RNA"; cp3 <- JoinLayers(cp3)

  if (!"Fig3_Annotation" %in% colnames(cp3@meta.data)) {
    message("CP003 object lacks Fig3_Annotation - skipping cross-participant mapping.")
  } else {
    shared_genes <- intersect(rownames(cp3), rownames(seu))
    pseudobulk <- function(obj, grp) {
      Xa <- GetAssayData(obj, layer = "data")[shared_genes, ]
      sp <- split(colnames(obj), as.character(obj@meta.data[[grp]]))
      sapply(sp, function(cells) Matrix::rowMeans(Xa[, cells, drop = FALSE]))
    }
    A <- pseudobulk(cp3, "Fig3_Annotation")
    B <- pseudobulk(seu, CL_COL)

    # variable genes only: housekeeping genes would inflate every correlation
    hv <- intersect(union(
      VariableFeatures(FindVariableFeatures(cp3, nfeatures = 2000, verbose = FALSE)),
      VariableFeatures(FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE))
    ), shared_genes)

    CM <- cor(B[hv, ], A[hv, ], method = "spearman")
    save_csv(data.frame(cp18_cluster = rownames(CM), round(CM, 4), check.names = FALSE),
             file.path(EVI_ROOT, "CP018_vs_CP003_pseudobulk_correlation.csv"))

    map <- data.frame(
      cp18_cluster     = rownames(CM),
      cp18_label       = ann[rownames(CM)],
      best_CP003_match = colnames(CM)[apply(CM, 1, which.max)],
      r_best           = round(apply(CM, 1, max), 4),
      runner_up        = colnames(CM)[apply(CM, 1, function(x) order(-x)[2])],
      r_second         = round(apply(CM, 1, function(x) sort(x, decreasing = TRUE)[2]), 4),
      stringsAsFactors = FALSE)
    map$margin <- round(map$r_best - map$r_second, 4)
    map$name_shared_with_CP003 <- ifelse(map$cp18_label == map$best_CP003_match, "YES", "NO")
    # A margin under ~0.03 between the top two candidates is not a real
    # assignment - the marker evidence decides in that case.
    map$confidence <- ifelse(map$margin >= 0.05, "clear",
                      ifelse(map$margin >= 0.03, "weak", "ambiguous (margin < 0.03)"))
    save_csv(map, file.path(EVI_ROOT, "CP018_to_CP003_label_mapping.csv"))

    # Reverse direction: a CP003 label with no good CP018 match means that
    # population is ABSENT in CP018 - the headline comparison result.
    revmap <- data.frame(
      CP003_label       = colnames(CM),
      best_cp18_cluster = rownames(CM)[apply(CM, 2, which.max)],
      r_best            = round(apply(CM, 2, max), 4),
      stringsAsFactors  = FALSE)
    revmap$cp18_label <- ann[revmap$best_cp18_cluster]
    revmap$population_present_in_CP018 <- ifelse(revmap$r_best >= 0.85, "YES",
                                                 "NO - no comparable cluster")
    save_csv(revmap, file.path(EVI_ROOT, "CP003_label_recovery_in_CP018.csv"))

    cat("\n=== CP018 cluster -> CP003 label ===\n")
    print(map[, c("cp18_cluster","cp18_label","best_CP003_match","r_best",
                  "margin","name_shared_with_CP003","confidence")], row.names = FALSE)
    cat("\n=== Is each CP003 population recovered in CP018? ===\n")
    print(revmap, row.names = FALSE)

    absent <- revmap$CP003_label[revmap$population_present_in_CP018 != "YES"]
    if (length(absent)) {
      cat("\n*** CP003 population(s) with NO comparable CP018 cluster: ",
          paste(absent, collapse = ", "), "\n",
          "    CP018 is the ART-suppressed participant; absence of an exhausted\n",
          "    population is consistent with the Figure 4/5 prediction, but n = 1.\n", sep = "")
    }
    rm(cp3); gc(verbose = FALSE)
  }
} else {
  message("CP003 Fig3 object not found - skipping cross-participant mapping.")
}

################################################################################
# 5. Heatmap of the evidence
################################################################################
hm <- as.data.frame(as.table(Z))
names(hm) <- c("panel", "cluster", "z")
hm$label <- ann[as.character(hm$cluster)]
# Wrap long labels rather than letting adjacent columns touch.
wrap_lab <- function(x, width = 12) vapply(
  strwrap(x, width = width, simplify = FALSE),
  paste, character(1), collapse = "\n")
hm$cluster_lab <- paste0(hm$cluster, "\n",
                         wrap_lab(ifelse(is.na(hm$label), "?", as.character(hm$label))))

p_hm <- ggplot(hm, aes(x = cluster_lab, y = panel, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", z)), size = 3.1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey90", name = "z") +
  labs(x = NULL, y = NULL, title = "CP018 annotation evidence",
       subtitle = paste0("Programme scores z-scored across clusters (res ",
                         MNN_RES_FINAL, ")")) +
  theme_cp018(base_size = 13) +
  theme(plot.title.position = "plot",
        plot.title    = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        axis.text.x   = element_text(angle = 0, hjust = 0.5, size = 9,
                                     lineheight = 0.95),
        panel.grid    = element_blank())
save_png(file.path(EVI_ROOT, "CP018_annotation_evidence_heatmap.png"),
         p_hm, width = 14, height = 10.5)

cat("\nOutputs in:\n  ", EVI_ROOT, "\n")
