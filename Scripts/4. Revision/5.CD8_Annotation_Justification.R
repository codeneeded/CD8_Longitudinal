################################################################################
# NAMING: panels are written under their reviewer-keyed names (R1M2_0n_*)
# directly. They were originally written as RevFig_5* and renamed by a one-off
# pass, which meant a re-run produced a fresh RevFig_* file while the cited
# R1M2_* copy stayed stale - two versions of one panel, with the figure PDF
# pointing at the old one. Writing the final name removes that failure mode.
# Revision script 5 — annotation justification for every CD8 population
#
# ANSWERS
#   Rev 2 major 1b  "is 'naive 2' really CD8+ T cells? What is the CD4
#                    expression in this cluster?"
#   Rev 1 minor 4   CD8 Tn cells clustering with CD4 T cells
#   Rev 1 minor 5   the term "Tpex" is not explained
#   Rev 1 minor 1   central memory cells are not mentioned
#   Rev 2 major 5c  CD16+ effector CD8 described in congenital CMV (PMID 39531313)
#
# WHY THIS EXISTS BEYOND THOSE ITEMS
#   Script 4's identity panel surfaced two populations whose CD8 protein is an
#   order of magnitude below the rest:
#       Naive CD8 4          ADT CD8A  1.70,  CD4 detected in 11% of cells
#       CD16+ Effector CD8   ADT CD8A  4.38,  CD4 detected in  6.9%
#   against ~15-21 in every other population. Sequencing depth does NOT explain
#   it (median genes 1,861 and 1,757; both mid-range). Both are load-bearing:
#   Naive CD8 4 is one of the four naive subsets in the abstract, and CD16+
#   Effector is a named novel population in the discussion.
#
#   So this script tests each population against competing explanations rather
#   than asserting a label:
#     1. lineage       CD8A/CD8B/CD4 at RNA and protein
#     2. dropout       is low CD8 explained by depth?
#     3. doublet       scDblFinder score, and CD4/B/mono/NK cross-lineage signal
#     4. participant   is it one participant's batch artefact?
#     5. programme     naive / memory / effector / exhaustion / IFN / stress
#   and flags any population whose evidence contradicts its label.
#
# OUTPUT
#   R1M2_02_lineage_identity_RNA_ADT.png   both modalities side by side
#   R1M2_03_qc_vs_CD8.png                  dropout test
#   R1M2_04_crosslineage_doublet.png       doublet / contamination test
#   R1M2_05_participant_composition.png    batch test
#   R1M2_06_programme_heatmap.png          programme scores
#   Tables/R2M1d_01_CD8_annotation_justification.csv  the evidence table + verdict
#   saved_R_data/rev_cd8_annotated_scored.qs2   checkpoint WITH module scores
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

################################################################################
# Programme definitions. Technical panels are included deliberately - a cluster
# driven by stress or hypoxia is a real cell state, but it must not be read as
# a differentiation stage.
################################################################################
PROGRAMMES <- list(
  Naive           = c("CCR7","SELL","LEF1","TCF7","IL7R","BACH2","FOXP1","MYB","NOSIP"),
  Tscm            = c("CD27","CD28","FAS","IL7R","CXCR3","TCF7","BCL2"),
  Tcm             = c("CCR7","SELL","IL7R","CD27","CD28","GPR183"),
  Effector        = c("GZMB","GNLY","PRF1","NKG7","KLRG1","FGFBP2","FCGR3A"),
  Cytotoxicity    = c("GZMA","GZMH","GZMK","CTSW","KLRD1","CST7"),
  Exhaustion      = c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1","CD160"),
  Innate_like     = c("KLRC1","KLRC2","KIR2DL3","KIR3DL1","NCAM1","TYROBP","FCER1G"),
  MAIT            = c("SLC4A10","KLRB1","NCR3","ZBTB16","RORC","IL18R1"),
  Proliferation   = c("MKI67","TOP2A","STMN1","TUBB","CCNB1","CDK1"),
  IFN_I           = c("ISG15","IFI6","MX1","OAS1","STAT1","IRF7","IFI44L"),
  Stress_HeatShock = c("HSPA1A","HSPA1B","HSPB1","DNAJB1","HSPH1","JUN","FOS"),
  Hypoxia_Glyco   = c("HILPDA","IGFBP2","PFKFB4","VEGFA","NDRG1","SLC2A1","LDHA")
)
LINEAGE_CHECK <- list(
  CD8_lineage = c("CD8A","CD8B"),
  # CD4-RESTRICTED ONLY. IL7R and MAL are expressed by naive/memory CD8 cells
  # too, so including them made this panel score high on every naive population
  # - it flagged 6 of 10 as CD4-contaminated when they are not.
  CD4_lineage = c("CD4","CD40LG"),
  B_lineage   = c("MS4A1","CD79A","CD79B","IGHM"),
  Mono_lineage = c("LYZ","CD14","S100A8","S100A9","FCN1"),
  NK_lineage  = c("NCAM1","NKG7","KLRF1","NCR1"),
  Platelet    = c("PPBP","PF4")
)

################################################################################
# Load, or reuse the scored checkpoint
################################################################################
CKPT <- "cd8_annotated_scored"
if (has_qs(CKPT)) {
  message("Reusing scored checkpoint - delete it to rescore.")
  cd8 <- load_qs(CKPT)
} else {
  require_obj(QS_CD8)
  cd8 <- qs_read(QS_CD8)
  DefaultAssay(cd8) <- "RNA"
  if (length(Layers(cd8[["RNA"]])) > 3) cd8 <- JoinLayers(cd8, assay = "RNA")

  for (nm in names(c(PROGRAMMES, LINEAGE_CHECK))) {
    g <- c(PROGRAMMES, LINEAGE_CHECK)[[nm]]
    g <- intersect(g, rownames(cd8))
    if (length(g) < 2) { message("skipping ", nm, " - <2 genes present"); next }
    cd8 <- AddModuleScore(cd8, features = list(g), name = paste0("SC_", nm),
                          nbin = 20, ctrl = min(50, length(g) * 10), seed = 42)
  }
  # AddModuleScore appends "1"; drop it so column names match the panel names.
  sc <- grep("^SC_", colnames(cd8@meta.data), value = TRUE)
  for (s in sc) {
    clean <- sub("1$", "", s)
    if (clean != s) {
      cd8@meta.data[[clean]] <- cd8@meta.data[[s]]
      cd8@meta.data[[s]] <- NULL
    }
  }
  save_qs(cd8, CKPT, note = "CD8 object + programme and lineage module scores")
}

md <- cd8@meta.data
md$CD8_Annotation <- factor(md$CD8_Annotation,
  levels = intersect(names(CD8_COLORS), unique(md$CD8_Annotation)))
pops <- levels(md$CD8_Annotation)
cat("\nPopulations:", length(pops), "|", ncol(cd8), "cells\n")

################################################################################
# PANEL A — lineage identity, RNA and protein together
################################################################################
DefaultAssay(cd8) <- "RNA"
rna_id <- intersect(c("CD8A","CD8B","CD4","CD3E","CD3D"), rownames(cd8))
rna_e  <- FetchData(cd8, vars = rna_id, layer = "data")
rna_sum <- do.call(rbind, lapply(rna_id, function(g) data.frame(
  population = md$CD8_Annotation, gene = g, expr = rna_e[[g]])) ) |>
  dplyr::group_by(population, gene) |>
  dplyr::summarise(pct_pos = 100 * mean(expr > 0),
                   mean_expr = mean(expr), .groups = "drop")
rna_sum$modality <- "RNA"

DefaultAssay(cd8) <- "ADT"
adt_id <- intersect(c("CD8A","CD4","CD3D","CD45RA","CD45RO","CD27"), rownames(cd8))
adt_e  <- FetchData(cd8, vars = adt_id, layer = "data")
adt_sum <- do.call(rbind, lapply(adt_id, function(g) data.frame(
  population = md$CD8_Annotation, gene = g, expr = adt_e[[g]])) ) |>
  dplyr::group_by(population, gene) |>
  dplyr::summarise(pct_pos = 100 * mean(expr > 0),
                   mean_expr = mean(expr), .groups = "drop")
adt_sum$modality <- "Protein (ADT, DSB)"
DefaultAssay(cd8) <- "RNA"

id_all <- rbind(rna_sum, adt_sum)
save_src("R1M2_02_lineage_identity_RNA_ADT", id_all)

pA <- ggplot(id_all, aes(x = gene, y = population)) +
  geom_point(aes(size = pct_pos, colour = mean_expr)) +
  geom_text(aes(label = round(mean_expr, 1)), size = 2.3, nudge_x = 0.32, hjust = 0) +
  facet_wrap(~ modality, scales = "free_x") +
  scale_colour_gradient(low = "grey90", high = "#C0392B", name = "Mean") +
  scale_size_area(max_size = 9, name = "% positive") +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.0))) +
  labs(title = "Lineage identity of every CD8 population, RNA and protein",
       subtitle = paste0("Surface protein and transcript for each population, side by side; ",
                         "", n_cells_label(cd8)),
       x = NULL, y = NULL,
       caption = paste0("Protein CD8A around 15-21 indicates a confidently ",
                        "CD8+ population. Values far below that are flagged in ",
                        "the justification table.")) +
  theme_rev(10)
save_fig("R1M2_02_lineage_identity_RNA_ADT", pA, width = 12, height = 6.5)

################################################################################
# PANEL B — dropout test: is low CD8 just low depth?
################################################################################
adt_cd8 <- if ("CD8A" %in% adt_id) adt_e$CD8A else rep(NA_real_, nrow(md))
qc <- data.frame(population = md$CD8_Annotation,
                 genes = md$nFeature_RNA, umis = md$nCount_RNA,
                 pct_mt = if ("percent.mt" %in% names(md)) md$percent.mt else NA,
                 adt_CD8A = adt_cd8) |>
  dplyr::group_by(population) |>
  dplyr::summarise(n = dplyr::n(),
                   median_genes = median(genes), median_umi = median(umis),
                   median_mt = round(median(pct_mt, na.rm = TRUE), 2),
                   adt_CD8A = round(mean(adt_CD8A, na.rm = TRUE), 2),
                   .groups = "drop")
save_src("R1M2_03_qc_vs_CD8", qc)

pB <- ggplot(qc, aes(x = median_genes, y = adt_CD8A)) +
  geom_point(aes(size = n, colour = population)) +
  ggrepel::geom_text_repel(aes(label = population), size = 3,
                           min.segment.length = 0, box.padding = 0.5) +
  scale_colour_manual(values = CD8_COLORS, guide = "none") +
  scale_size_area(max_size = 9, name = "Cells") +
  labs(title = "Surface CD8A against sequencing depth, per population",
       subtitle = paste0("If low CD8A tracked depth these points would sit on ",
                         "a rising line; ", n_cells_label(cd8)),
       x = "Median genes detected per cell", y = "Mean protein CD8A (DSB)",
       caption = paste0("A population low on the y-axis but mid-range on the ",
                        "x-axis is NOT a dropout artefact - its CD8 protein is ",
                        "genuinely low.")) +
  theme_rev(11)
save_fig("R1M2_03_qc_vs_CD8", pB, width = 9, height = 6.5)

################################################################################
# PANEL C — cross-lineage and doublet test
################################################################################
lin_cols <- paste0("SC_", names(LINEAGE_CHECK))
lin_cols <- intersect(lin_cols, names(md))
dbl_col  <- grep("scDblFinder.score|doublet_score", names(md), value = TRUE)[1]

cl <- md[, c("CD8_Annotation", lin_cols), drop = FALSE]
cl_long <- tidyr::pivot_longer(cl, cols = all_of(lin_cols),
                               names_to = "lineage", values_to = "score")
cl_long$lineage <- sub("^SC_", "", sub("_lineage$", "", cl_long$lineage))
cl_sum <- cl_long |> dplyr::group_by(CD8_Annotation, lineage) |>
  dplyr::summarise(mean_score = mean(score), .groups = "drop")
if (!is.na(dbl_col)) {
  d <- md |> dplyr::group_by(CD8_Annotation) |>
    dplyr::summarise(doublet_score = round(mean(.data[[dbl_col]], na.rm = TRUE), 3),
                     .groups = "drop")
  cl_sum <- dplyr::left_join(cl_sum, d, by = "CD8_Annotation")
}
save_src("R1M2_04_crosslineage_doublet", cl_sum)

pC <- ggplot(cl_sum, aes(x = lineage, y = CD8_Annotation, fill = mean_score)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(mean_score, 2)), size = 2.7) +
  scale_fill_gradient2(low = "#4A90D9", mid = "grey95", high = "#C0392B",
                       midpoint = 0, name = "Module\nscore") +
  labs(title = "Cross-lineage marker scores per CD8 population",
       subtitle = paste0("A genuinely CD8 population scores high on CD8 and ",
                         "near zero elsewhere; ", n_cells_label(cd8)),
       x = NULL, y = NULL,
       caption = if (!is.na(dbl_col))
         paste0("Doublet scores in the source CSV (column ", dbl_col, ").")
         else "No doublet score column in this object.") +
  theme_rev(10)
save_fig("R1M2_04_crosslineage_doublet", pC, width = 8.5, height = 6)

################################################################################
# PANEL D — participant composition: batch artefact test
################################################################################
pc <- md |> dplyr::group_by(CD8_Annotation, PID) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") |>
  dplyr::mutate(pct = 100 * n / sum(n)) |> dplyr::ungroup()
conc <- pc |> dplyr::group_by(CD8_Annotation) |>
  dplyr::summarise(n_participants = dplyr::n_distinct(PID),
                   max_single_pct = round(max(pct), 1), .groups = "drop")
save_src("R1M2_05_participant_composition", pc)
save_tab("CD8_population_participant_concentration", conc)

pD <- ggplot(pc, aes(x = CD8_Annotation, y = pct, fill = PID)) +
  geom_col(width = 0.75) +
  geom_text(data = conc, inherit.aes = FALSE,
            aes(x = CD8_Annotation, y = 103,
                label = paste0(n_participants, " ppts")), size = 2.7) +
  labs(title = "Participant composition of each CD8 population",
       subtitle = paste0("A population from one participant only is a batch ",
                         "artefact, not a cell state; ", n_cells_label(cd8)),
       x = NULL, y = "% of the population's cells",
       caption = paste0("Bars are participant shares within each population. ",
                        "Label is how many participants contribute.")) +
  scale_y_continuous(limits = c(0, 107), breaks = seq(0, 100, 25)) +
  theme_rev(10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right", legend.key.size = unit(0.32, "cm"))
save_fig("R1M2_05_participant_composition", pD, width = 11, height = 6.5)

################################################################################
# PANEL E — programme scores, z-scored across populations
################################################################################
prog_cols <- intersect(paste0("SC_", names(PROGRAMMES)), names(md))
pr <- md[, c("CD8_Annotation", prog_cols), drop = FALSE] |>
  dplyr::group_by(CD8_Annotation) |>
  dplyr::summarise(dplyr::across(all_of(prog_cols), mean), .groups = "drop")
prz <- pr
for (cc in prog_cols) {
  v <- prz[[cc]]
  prz[[cc]] <- if (sd(v) > 0) (v - mean(v)) / sd(v) else 0
}
prz_long <- tidyr::pivot_longer(prz, cols = all_of(prog_cols),
                                names_to = "programme", values_to = "z")
prz_long$programme <- sub("^SC_", "", prz_long$programme)
save_src("R1M2_06_programme_heatmap", prz_long)
save_tab("CD8_programme_scores_raw", pr)

pE <- ggplot(prz_long, aes(x = programme, y = CD8_Annotation, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(z, 1)), size = 2.6) +
  scale_fill_gradient2(low = "#4A90D9", mid = "grey95", high = "#C0392B",
                       midpoint = 0, name = "z across\npopulations") +
  labs(title = "Transcriptional programme per CD8 population",
       subtitle = paste0("z-scored across populations, so colour is relative ",
                         "not absolute; ", n_cells_label(cd8)),
       x = NULL, y = NULL,
       caption = paste0("Stress and hypoxia panels are included on purpose: a ",
                        "population driven by either is a real state but must ",
                        "not be read as a differentiation stage.")) +
  theme_rev(10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig("R1M2_06_programme_heatmap", pE, width = 11, height = 6)

################################################################################
# THE JUSTIFICATION TABLE — evidence plus a rule-based verdict
################################################################################
wide_id <- id_all |>
  dplyr::mutate(key = paste0(modality, "_", gene)) |>
  dplyr::select(population, key, mean_expr, pct_pos)

getv <- function(pop, modality, gene, what = "mean_expr") {
  r <- id_all[id_all$population == pop & id_all$modality == modality &
              id_all$gene == gene, ]
  if (!nrow(r)) return(NA_real_)
  round(r[[what]][1], 2)
}
ADT <- "Protein (ADT, DSB)"
adt_cd8_all <- vapply(pops, function(p) getv(p, ADT, "CD8A"), numeric(1))
adt_ok <- median(adt_cd8_all, na.rm = TRUE)

rows <- lapply(pops, function(p) {
  q <- qc[qc$population == p, ]
  cc <- cl_sum[cl_sum$CD8_Annotation == p, ]
  pz <- prz_long[prz_long$CD8_Annotation == p, ]
  # z is RELATIVE across populations, so every population is "top" at something.
  # Keep the relative winner but note whether its RAW score is actually positive
  # - a population can lead on a programme it barely expresses.
  top <- pz$programme[order(-pz$z)][1:2]
  raw_top <- pr[[paste0("SC_", top[1])]][pr$CD8_Annotation == p]
  raw_top <- if (length(raw_top)) round(raw_top[1], 3) else NA
  a_cd8 <- getv(p, ADT, "CD8A"); a_cd4 <- getv(p, ADT, "CD4")
  r_cd8 <- getv(p, "RNA", "CD8A", "pct_pos"); r_cd4 <- getv(p, "RNA", "CD4", "pct_pos")
  dbl <- if ("doublet_score" %in% names(cc)) cc$doublet_score[1] else NA
  cd4_lin <- cc$mean_score[cc$lineage == "CD4"]
  cd4_lin <- if (length(cd4_lin)) round(cd4_lin[1], 3) else NA

  flags <- c()
  if (!is.na(a_cd8) && a_cd8 < adt_ok / 3)
    flags <- c(flags, sprintf("CD8 protein %.2f vs cohort median %.2f", a_cd8, adt_ok))
  if (!is.na(r_cd4) && r_cd4 > 5)
    flags <- c(flags, sprintf("CD4 RNA in %.1f%% of cells", r_cd4))
  # Only flag CD4 signal when BOTH the restricted module and CD4 RNA detection
  # are elevated - either alone is noise at this panel size.
  if (!is.na(cd4_lin) && cd4_lin > 0.10 && !is.na(r_cd4) && r_cd4 > 5)
    flags <- c(flags, sprintf("CD4 lineage module %.3f with CD4 RNA %.1f%%",
                              cd4_lin, r_cd4))
  if (!is.na(dbl) && dbl > 0.3)
    flags <- c(flags, sprintf("doublet score %.3f", dbl))
  cn <- conc[conc$CD8_Annotation == p, ]
  if (nrow(cn) && cn$max_single_pct[1] > 60)
    flags <- c(flags, sprintf("%.0f%% from one participant", cn$max_single_pct[1]))
  if (!is.na(q$median_genes[1]) && q$median_genes[1] < 700)
    flags <- c(flags, sprintf("low complexity, %d genes", q$median_genes[1]))

  data.frame(
    population = p, n_cells = q$n[1],
    RNA_CD8A_pct = r_cd8, RNA_CD4_pct = r_cd4,
    ADT_CD8A = a_cd8, ADT_CD4 = a_cd4,
    CD4_lineage_module = cd4_lin, doublet_score = dbl,
    median_genes = q$median_genes[1], median_UMI = q$median_umi[1],
    n_participants = if (nrow(cn)) cn$n_participants[1] else NA,
    max_single_participant_pct = if (nrow(cn)) cn$max_single_pct[1] else NA,
    top_programmes = paste(top, collapse = " / "),
    top_programme_raw_score = raw_top,
    verdict = if (length(flags)) "REVIEW" else "supported",
    flags = if (length(flags)) paste(flags, collapse = "; ") else "")
})
just <- do.call(rbind, rows)
save_tab("R2M1d_01_CD8_annotation_justification", just)

cat("\n=== ANNOTATION JUSTIFICATION ===\n")
print(just[, c("population","n_cells","ADT_CD8A","RNA_CD4_pct",
               "n_participants","top_programmes","verdict")], row.names = FALSE)

flagged <- just[just$verdict == "REVIEW", ]
cat("\n=== POPULATIONS NEEDING A DECISION BEFORE THE RESPONSE LETTER ===\n")
if (nrow(flagged)) {
  for (i in seq_len(nrow(flagged))) {
    cat("\n  ", flagged$population[i], " (", flagged$n_cells[i], " cells)\n", sep = "")
    cat("    ", flagged$flags[i], "\n", sep = "")
    cat("     top programmes: ", flagged$top_programmes[i], "\n", sep = "")
  }
  cat("\n  These are NOT automatically wrong - a low-CD8 population can be a real\n")
  cat("  state. But each needs an explicit call in the manuscript, because a\n")
  cat("  reviewer running this same check will ask.\n")
} else {
  cat("  none - every label is consistent with its evidence.\n")
}

cat("\nFigures: ", REV_FIG, "\nTables:  ", REV_TAB, "\n")
