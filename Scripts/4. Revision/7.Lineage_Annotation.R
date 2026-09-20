################################################################################
# Revision script 7 — annotate each re-clustered lineage
#
# ANSWERS  Rev 1 major 2 (analyse the other lineages)
#
# WHAT THIS DOES
#   Takes the re-clustered checkpoint from script 6 and produces the evidence
#   needed to LABEL its clusters: marker profiles at RNA and protein, canonical
#   panel scores, reference-mapping agreement, and the same technical-confounder
#   checks used for CD8 in script 5 (depth, doublets, participant concentration).
#
# IT DOES NOT INVENT LABELS
#   The script writes an evidence table and a PROPOSED label per cluster from a
#   transparent rule (highest-scoring canonical panel, with a margin test). Any
#   cluster whose top panel does not clear the margin is left as
#   "UNRESOLVED - inspect" rather than given a confident name. The proposed
#   labels must be reviewed against the marker plots before they go in a figure,
#   exactly as the CD8 annotation was.
#
# OUTPUT
#   Analysis/<LIN>_subset/03_annotation/*.png            marker evidence
#   Analysis/<LIN>_subset/03_annotation/<LIN>_cluster_evidence.csv
#   saved_R_data/rev_<lin>_annotated.qs2                 checkpoint
#
#   Reviewer-facing panels come from script 8, not here. This is annotation
#   evidence, which belongs in Analysis/ per the repo convention.
#
# USAGE
#   Rscript "Scripts/4. Revision/7.Lineage_Annotation.R" CD4
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

args <- commandArgs(trailingOnly = TRUE)
LIN  <- if (length(args)) args[1] else stop(
  "Usage: Rscript 7.Lineage_Annotation.R <CD4|NK|B|Myeloid>", call. = FALSE)

OUT <- file.path(BASE, "Analysis", paste0(LIN, "_subset"), "03_annotation")
mkdirs(OUT)
CK_IN  <- paste0(tolower(LIN), "_reclustered")
CK_OUT <- paste0(tolower(LIN), "_annotated")

################################################################################
# Canonical panels per lineage. Deliberately conservative: only genes whose
# expression is restricted enough to discriminate WITHIN the lineage.
################################################################################
PANELS <- list(
  CD4 = list(
    Naive              = c("CCR7","SELL","LEF1","TCF7","NOSIP","MAL"),
    Central_memory     = c("CCR7","SELL","IL7R","CD27","GPR183","AQP3"),
    Effector_memory    = c("IL7R","LTB","AQP3","S100A4","KLRB1","CCL5"),
    Th1                = c("TBX21","CXCR3","IFNG","IL12RB2","CCR5"),
    Th2                = c("GATA3","CCR4","IL4R","PTGDR2","IL17RB"),
    Th17               = c("RORC","CCR6","IL23R","KLRB1","IL17A"),
    Treg               = c("FOXP3","IL2RA","IKZF2","CTLA4","TNFRSF4","RTKN2"),
    Tfh_like           = c("CXCR5","BCL6","ICOS","PDCD1","TOX2"),
    Cytotoxic_CD4      = c("GZMB","GNLY","PRF1","NKG7","KLRD1"),
    Proliferating      = c("MKI67","TOP2A","STMN1","TUBB"),
    IFN_stimulated     = c("ISG15","IFI6","MX1","OAS1","IFI44L","STAT1")
  ),
  NK = list(
    CD56bright         = c("NCAM1","GZMK","SELL","XCL1","IL7R","KLRC1"),
    CD56dim_mature     = c("FCGR3A","FGFBP2","GZMB","PRF1","CX3CR1","S1PR5"),
    Adaptive_NK        = c("KLRC2","IL32","CD3E","LAG3","CD52"),
    Cytokine_activated = c("IFNG","CCL3","CCL4","TNF","XCL2"),
    Proliferating      = c("MKI67","TOP2A","STMN1"),
    IFN_stimulated     = c("ISG15","IFI6","MX1","OAS1","IFI44L")
  ),
  B = list(
    Naive_B            = c("TCL1A","IGHD","IGHM","FCER2","CD72","YBX3"),
    Transitional_B     = c("VPREB3","MME","CD9","SOX4","IGLL5"),
    Memory_B           = c("CD27","TNFRSF13B","AIM2","IGHG1","IGHA1"),
    Activated_B        = c("CD69","CD83","NFKBIA","JUN","FOS","EGR1"),
    Atypical_B         = c("ITGAX","TBX21","FCRL5","ZEB2","SIGLEC6"),
    Plasmablast        = c("MZB1","XBP1","JCHAIN","PRDM1","DERL3"),
    Proliferating      = c("MKI67","TOP2A","STMN1"),
    IFN_stimulated     = c("ISG15","IFI6","MX1","OAS1","IFI44L")
  ),
  Myeloid = list(
    Classical_mono     = c("CD14","S100A8","S100A9","LYZ","VCAN","FCN1"),
    Nonclassical_mono  = c("FCGR3A","MS4A7","CDKN1C","LST1","CSF1R"),
    Intermediate_mono  = c("CD14","FCGR3A","HLA-DRA","CD74"),
    cDC1               = c("CLEC9A","XCR1","IDO1","BATF3"),
    cDC2               = c("CD1C","FCER1A","CLEC10A","HLA-DQA1"),
    pDC                = c("LILRA4","IL3RA","CLEC4C","TCF4","IRF7"),
    IFN_stimulated     = c("ISG15","IFI6","MX1","OAS1","IFI44L")
  )
)
if (!LIN %in% names(PANELS))
  stop("no panel defined for ", LIN, call. = FALSE)

# Technical confounders - the same checks applied to CD8 in script 5.
TECHNICAL <- list(
  Stress_HeatShock = c("HSPA1A","HSPA1B","HSPB1","DNAJB1","JUN","FOS"),
  Hypoxia_Glyco    = c("HILPDA","IGFBP2","PFKFB4","VEGFA","NDRG1","SLC2A1"),
  Mito_high        = c("MT-CO1","MT-CO2","MT-ND1","MT-ATP6")
)
# Cross-lineage panels - a cluster scoring high here is contamination/doublets.
# NOTE on gene choice: TRBC2 is deliberately EXCLUDED from the T-cell panel.
# It is detected in ~72% of cells in every NK cluster of this dataset (ambient
# transcript or mismapping), so including it dilutes the panel to the point
# where a genuinely CD3-high cluster escapes the flag.
CROSS <- list(
  T_cell   = c("CD3E","CD3D","CD3G"),
  B_cell   = c("MS4A1","CD79A","CD79B"),
  NK_cell  = c("NKG7","KLRD1","GNLY","NCR1"),
  Myeloid  = c("LYZ","CD14","S100A8","FCN1"),
  Platelet = c("PPBP","PF4","GP9")
)
# A cross-lineage call needs the marker genuinely DETECTED, not merely ranked
# highest by z-score - z-scoring across clusters guarantees some cluster is top
# at everything. Require this fraction of cells positive for the panel's genes.
CROSS_MIN_PCT <- 20

################################################################################
# Load
################################################################################
if (!has_qs(CK_IN))
  stop("Run script 6 first: no checkpoint ", rev_qs_path(CK_IN), call. = FALSE)
# Report which integration produced the input, so an annotation built on a CCA
# object is never mistaken for one built on the published FastMNN method.
.mf <- file.path(REV_TAB, "rev_object_manifest.csv")
if (file.exists(.mf)) {
  .m <- read.csv(.mf, stringsAsFactors = FALSE)
  .n <- .m$note[.m$name == CK_IN]
  if (length(.n)) {
    .meth <- if (grepl("FastMNN", tail(.n, 1), fixed = TRUE)) "FastMNN"
             else if (grepl("CCA", tail(.n, 1), fixed = TRUE)) "CCA" else "unrecorded"
    message("Input checkpoint integration method: ", .meth)
    if (.meth == "CCA")
      warning("Annotating a CCA-integrated object. The published CD8 analysis used ",
              "FastMNN, so these labels are NOT directly comparable to it. ",
              "Re-run script 6 without INTEGRATION=cca before using any figure.",
              call. = FALSE)
  }
}
obj <- load_qs(CK_IN)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]])) > 3) obj <- JoinLayers(obj, assay = "RNA")
md <- obj@meta.data
stopifnot("seurat_clusters" %in% names(md))
cl <- factor(md$seurat_clusters)
cat("\n", LIN, ":", ncol(obj), "cells,", nlevels(cl), "clusters\n\n")

################################################################################
# 1. Module scores for every panel
################################################################################
ALL_PANELS <- c(PANELS[[LIN]], TECHNICAL, CROSS)
for (nm in names(ALL_PANELS)) {
  g <- intersect(ALL_PANELS[[nm]], rownames(obj))
  if (length(g) < 2) { message("  skip ", nm, " (<2 genes present)"); next }
  obj <- AddModuleScore(obj, features = list(g), name = paste0("P_", nm),
                        nbin = 20, ctrl = min(50, length(g) * 10), seed = 42)
}
for (s in grep("^P_", colnames(obj@meta.data), value = TRUE)) {
  clean <- sub("1$", "", s)
  if (clean != s) { obj@meta.data[[clean]] <- obj@meta.data[[s]]
                    obj@meta.data[[s]] <- NULL }
}
md <- obj@meta.data

score_cols <- grep("^P_", names(md), value = TRUE)
per_cl <- md[, c("seurat_clusters", score_cols)] |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(dplyr::across(all_of(score_cols), mean), .groups = "drop")
write.csv(per_cl, file.path(OUT, paste0(LIN, "_panel_scores_raw.csv")), row.names = FALSE)

# z-score across clusters so panels of different baseline are comparable
z <- per_cl
for (cc in score_cols) {
  v <- z[[cc]]; z[[cc]] <- if (sd(v) > 0) (v - mean(v)) / sd(v) else 0
}
zl <- tidyr::pivot_longer(z, cols = all_of(score_cols),
                          names_to = "panel", values_to = "zscore")
zl$panel <- sub("^P_", "", zl$panel)
zl$kind <- ifelse(zl$panel %in% names(PANELS[[LIN]]), "identity",
           ifelse(zl$panel %in% names(TECHNICAL), "technical", "cross-lineage"))
write.csv(zl, file.path(OUT, paste0(LIN, "_panel_scores_z.csv")), row.names = FALSE)

p_hm <- ggplot(zl, aes(x = panel, y = seurat_clusters, fill = zscore)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(zscore, 1)), size = 2.4) +
  facet_grid(~ kind, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#4A90D9", mid = "grey95", high = "#C0392B",
                       midpoint = 0, name = "z across\nclusters") +
  labs(title = paste0(LIN, " cluster identity evidence"),
       subtitle = paste0("z-scored across clusters; ", n_cells_label(obj)),
       x = NULL, y = "WNN cluster",
       caption = paste0("technical and cross-lineage panels are shown on ",
                        "purpose - a cluster driven by either must not be ",
                        "given a biological identity.")) +
  theme_rev(9.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, paste0(LIN, "_panel_heatmap.png")), p_hm,
       width = 15, height = 1.4 + 0.42 * nlevels(cl), dpi = 350,
       bg = "white", limitsize = FALSE)

################################################################################
# 2. Reference-mapping agreement (independent of our clustering)
################################################################################
if ("predicted.celltype.l2" %in% names(md)) {
  ct <- as.data.frame(table(cluster = md$seurat_clusters,
                            azimuth = md$predicted.celltype.l2))
  ct <- ct[ct$Freq > 0, ]
  ct$pct <- round(100 * ct$Freq / ave(ct$Freq, ct$cluster, FUN = sum), 1)
  ct <- ct[ct$pct >= 5, ]
  write.csv(ct[order(ct$cluster, -ct$pct), ],
            file.path(OUT, paste0(LIN, "_azimuth_agreement.csv")), row.names = FALSE)

  p_az <- ggplot(ct, aes(x = azimuth, y = cluster, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = pct), size = 2.3) +
    scale_fill_gradient(low = "grey95", high = "#2C6E49", name = "% of\ncluster") +
    labs(title = paste0(LIN, " clusters vs Azimuth reference mapping"),
         subtitle = "Percentage of each cluster assigned to each reference cell type",
         x = NULL, y = "WNN cluster",
         caption = "Only assignments covering at least 5% of a cluster are shown.") +
    theme_rev(9.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT, paste0(LIN, "_azimuth_agreement.png")), p_az,
         width = 12, height = 1.4 + 0.42 * nlevels(cl), dpi = 350,
         bg = "white", limitsize = FALSE)
}

################################################################################
# 3. Differential markers per cluster
################################################################################
Idents(obj) <- "seurat_clusters"
mk <- tryCatch(
  FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, min.pct = 0.25,
                 logfc.threshold = 0.4, verbose = FALSE),
  error = function(e) { message("FindAllMarkers failed: ", conditionMessage(e)); NULL })
if (!is.null(mk) && nrow(mk)) {
  write.csv(mk, file.path(OUT, paste0(LIN, "_markers_all.csv")), row.names = FALSE)
  top <- mk |> dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = 8, with_ties = FALSE) |> dplyr::ungroup()
  write.csv(top, file.path(OUT, paste0(LIN, "_markers_top8.csv")), row.names = FALSE)
}

################################################################################
# 4. Technical QC per cluster
################################################################################
dbl <- grep("scDblFinder.score|doublet_score", names(md), value = TRUE)[1]
qc <- data.frame(cluster = md$seurat_clusters, genes = md$nFeature_RNA,
                 umi = md$nCount_RNA, PID = md$PID,
                 dbl = if (!is.na(dbl)) md[[dbl]] else NA_real_) |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   median_genes = median(genes), median_umi = median(umi),
                   doublet_score = round(mean(dbl, na.rm = TRUE), 3),
                   n_participants = dplyr::n_distinct(PID),
                   max_single_participant_pct =
                     round(100 * max(table(PID)) / dplyr::n(), 1),
                   .groups = "drop")
write.csv(qc, file.path(OUT, paste0(LIN, "_cluster_qc.csv")), row.names = FALSE)

################################################################################
# 5. PROPOSED labels — transparent rule, margin-tested, never asserted
################################################################################
ident_panels <- names(PANELS[[LIN]])
zi <- zl[zl$panel %in% ident_panels, ]

# Detection rate per cluster for each cross-lineage panel: the fraction of
# cells expressing ANY of that panel's genes. Used to gate the flag above.
cross_pct <- lapply(names(CROSS), function(nm) {
  g <- intersect(CROSS[[nm]], rownames(obj))
  if (!length(g)) return(setNames(rep(NA_real_, nlevels(cl)), levels(cl)))
  X <- GetAssayData(obj, assay = "RNA", layer = "data")[g, , drop = FALSE]
  pos <- Matrix::colSums(X > 0) > 0
  setNames(vapply(levels(cl), function(k)
    round(100 * mean(pos[cl == k]), 1), numeric(1)), levels(cl))
})
names(cross_pct) <- names(CROSS)
write.csv(data.frame(cluster = levels(cl),
                     as.data.frame(lapply(cross_pct, as.numeric))),
          file.path(OUT, paste0(LIN, "_cross_lineage_detection.csv")),
          row.names = FALSE)

prop <- lapply(levels(cl), function(k) {
  s <- zi[zi$seurat_clusters == k, ]
  s <- s[order(-s$zscore), ]
  q <- qc[qc$cluster == k, ]
  tech <- zl[zl$seurat_clusters == k & zl$kind == "technical", ]
  xl   <- zl[zl$seurat_clusters == k & zl$kind == "cross-lineage", ]

  top1 <- s$panel[1]; z1 <- s$zscore[1]
  top2 <- if (nrow(s) > 1) s$panel[2] else NA
  z2   <- if (nrow(s) > 1) s$zscore[2] else NA
  margin <- if (!is.na(z2)) round(z1 - z2, 2) else NA

  # Three-way verdict. A single "z AND margin" rule conflated two different
  # situations and refused clusters with strong absolute signal merely because a
  # RELATED panel scored close behind (e.g. Memory_B z=3.36 with Atypical_B 0.32
  # behind - those panels share genes, so a close runner-up is expected).
  #   confident  - strong signal, clear winner
  #   ambiguous  - strong signal but two panels within the margin: the label is
  #                one of a NAMED pair, which is far more useful than "unknown"
  #   unresolved - no panel stands out at all
  STRONG <- 1.5; WEAK <- 0.5; MARGIN <- 0.4
  verdict <- if (is.na(z1)) "UNRESOLVED - inspect"
    else if (z1 >= STRONG && (is.na(margin) || margin > MARGIN)) "confident"
    else if (z1 >= STRONG) "ambiguous"
    else if (z1 > WEAK && !is.na(margin) && margin > MARGIN) "confident"
    else "UNRESOLVED - inspect"
  ok <- verdict == "confident"
  worst_tech <- if (nrow(tech)) tech$panel[which.max(tech$zscore)] else NA
  worst_tz   <- if (nrow(tech)) round(max(tech$zscore), 2) else NA
  worst_xl   <- if (nrow(xl)) xl$panel[which.max(xl$zscore)] else NA
  worst_xz   <- if (nrow(xl)) round(max(xl$zscore), 2) else NA

  flags <- c()
  if (!is.na(worst_tz) && worst_tz > 1.5)
    flags <- c(flags, paste0("technical: ", worst_tech, " z=", worst_tz))
  # Cross-lineage needs BOTH a high z AND real detection of that panel's genes.
  if (!is.na(worst_xz) && worst_xz > 1.5 && !is.na(worst_xl)) {
    pctd <- cross_pct[[worst_xl]][[k]]
    if (!is.na(pctd) && pctd >= CROSS_MIN_PCT) {
      flags <- c(flags, paste0("cross-lineage: ", worst_xl, " z=", worst_xz,
                               ", ", pctd, "% of cells positive"))
    } else {
      flags <- c(flags, paste0("weak ", worst_xl, " signal (z=", worst_xz,
                               " but only ", pctd, "% positive - not called)"))
    }
  }
  if (!is.na(q$doublet_score[1]) && q$doublet_score[1] > 0.3)
    flags <- c(flags, paste0("doublet score ", q$doublet_score[1]))
  if (q$max_single_participant_pct[1] > 60)
    flags <- c(flags, paste0(q$max_single_participant_pct[1], "% one participant"))
  if (q$median_genes[1] < 700)
    flags <- c(flags, paste0("low complexity, ", q$median_genes[1], " genes"))

  data.frame(
    cluster = k, n_cells = q$n_cells[1],
    proposed_label = if (verdict == "confident") top1
                     else if (verdict == "ambiguous") paste0(top1, " or ", top2)
                     else "UNRESOLVED - inspect",
    call_quality = verdict,
    top_panel = top1, top_z = round(z1, 2),
    runner_up = top2, runner_up_z = round(z2, 2), margin = margin,
    median_genes = q$median_genes[1], doublet_score = q$doublet_score[1],
    n_participants = q$n_participants[1],
    max_single_participant_pct = q$max_single_participant_pct[1],
    top_marker = if (!is.null(mk) && nrow(mk))
      paste(head(mk$gene[mk$cluster == k], 5), collapse = ", ") else "",
    flags = if (length(flags)) paste(flags, collapse = "; ") else "")
})
ev <- do.call(rbind, prop)
write.csv(ev, file.path(OUT, paste0(LIN, "_cluster_evidence.csv")), row.names = FALSE)

# Attach the proposal, clearly named so it can never be mistaken for a final call
prop_vec <- setNames(ev$proposed_label, ev$cluster)[as.character(md$seurat_clusters)]
names(prop_vec) <- colnames(obj)
obj$Proposed_Annotation <- factor(prop_vec)

save_qs(obj, CK_OUT,
        note = paste0(LIN, " re-clustered + panel-scored; Proposed_Annotation ",
                      "is RULE-BASED and needs review"))

cat("\n=== PROPOSED LABELS (review before use) ===\n")
print(ev[, c("cluster","n_cells","proposed_label","call_quality","top_z","margin","flags")],
      row.names = FALSE)

cat("\n  confident:", sum(ev$call_quality == "confident"),
    "| ambiguous (named pair):", sum(ev$call_quality == "ambiguous"),
    "| unresolved:", sum(ev$call_quality == "UNRESOLVED - inspect"),
    "of", nrow(ev), "clusters\n")
fl <- ev[nzchar(ev$flags), ]
if (nrow(fl)) cat("  clusters carrying a technical flag:", nrow(fl), "\n")
cat("\n  Evidence in:", OUT, "\n")
cat("  Review ", LIN, "_panel_heatmap.png and ", LIN,
    "_azimuth_agreement.png, then set final labels.\n", sep = "")
