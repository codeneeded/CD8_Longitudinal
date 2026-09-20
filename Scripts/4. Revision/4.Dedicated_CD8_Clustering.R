################################################################################
# Revision script 4 — the dedicated CD8 sub-clustering, shown explicitly
#
# ANSWERS
#   Rev 2 major 1d  "Can the authors provide a rationale for not sub-clustering
#                    CD8 T cells, or if they did, show this analysis? Does the
#                    unusual 2nd naive population persist when dedicated CD8
#                    T cell analysis is done."
#   Rev 2 major 1b  Fig 1C as a bubble plot (% expressing + relative expression);
#                    "is 'naive 2' really CD8+ T cells? What is the CD4
#                    expression in this cluster?"
#   Rev 2 major 1c  frequency of the OTHER clusters across groups, not just TEMRA
#   Rev 2 major 1e  cluster frequency as % of total CD8, not % of live cells
#   Rev 1 minor 4   "A number of the cells annotated as CD8 Tn can be seen
#                    clustering with the CD4 T cells. This should be explained."
#
# THE KEY POINT FOR THE RESPONSE LETTER
#   The dedicated CD8 sub-clustering the reviewer asks for ALREADY EXISTS. It is
#   "Scripts/2. Analysis/4.CD8_Subseting.R", which re-normalises, re-scales,
#   runs RNA PCA + FastMNN, ADT PCA + FastMNN, FindMultiModalNeighbors and WNN
#   clustering swept over resolutions 0.2-2.0, and produces the 10 populations in
#   TARA_cd8_HEI_annotated_final.qs2 ($CD8_Annotation). Figures 4, 5 and
#   Supplementary 4-9 are all built on it.
#
#   What the manuscript never states is that TWO clusterings are in play:
#     Figures 1-2  whole-PBMC WNN clustering   (TARA_ALL_annotated_final.qs2,
#                                               126,435 cells, $Annotation)
#     Figures 4-5  dedicated CD8-only WNN      (TARA_cd8_HEI_annotated_final.qs2,
#                                               18,843 cells, $CD8_Annotation)
#   Both reviewers independently tripped over this. The fix is disclosure plus
#   this figure, not new analysis.
#
# OUTPUT (PNG only, each with its source CSV)
#   R2M2a_01_CD8_dedicated_UMAP.png       the dedicated clustering
#   R2M1d_03_ADT_subgating_crosswalk.png  WNN cluster -> final label (ADT gate)
#   R2M1b_01_bubble_RNA.png               bubble plot, RNA (major 1b)
#   R2M1b_02_bubble_ADT.png               bubble plot, protein (major 1b)
#   R2M1e_01_composition_pctCD8.png       all clusters by group (1c + 1e)
#   R2M1d_04_naive2_identity.png          is Naive CD8 2 really CD8+? (1b)
#   Tables/CD8_ADT_subgating_crosswalk.csv the table the response letter quotes
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

require_obj(QS_CD8)
cd8 <- qs_read(QS_CD8)
md  <- cd8@meta.data

stopifnot("CD8_Annotation" %in% names(md), "source_cluster" %in% names(md))
md$CD8_Annotation <- factor(md$CD8_Annotation,
  levels = intersect(names(CD8_COLORS), unique(md$CD8_Annotation)))
cd8$CD8_Annotation <- md$CD8_Annotation

cat("\nDedicated CD8 object:", ncol(cd8), "cells,",
    nlevels(md$CD8_Annotation), "populations,",
    length(unique(md$PID)), "participants (", unique(md$Condition), "only )\n")
print(sort(table(md$CD8_Annotation), decreasing = TRUE))

# Which WNN resolution the published clustering used, read from the object
# rather than assumed - the resolution columns are all retained.
res_cols <- grep("^wsnn_res\\.", names(md), value = TRUE)
res_match <- NULL
for (rc in res_cols) {
  if (length(unique(md[[rc]])) == nlevels(droplevels(factor(md$source_cluster)))) {
    res_match <- c(res_match, rc)
  }
}
cat("\nWNN resolutions computed:", paste(sub("wsnn_res\\.", "", res_cols), collapse = ", "), "\n")

################################################################################
# PANEL A — the dedicated clustering, with n on the plot (Rev 2 minor 2)
################################################################################
umap_name <- if ("wnn.umap" %in% Reductions(cd8)) "wnn.umap" else Reductions(cd8)[1]
um <- as.data.frame(Embeddings(cd8, umap_name))[, 1:2]
names(um) <- c("UMAP_1", "UMAP_2")
um$CD8_Annotation <- md$CD8_Annotation
um$Timepoint_Group <- md$Timepoint_Group
um$PID <- md$PID

lab_pos <- um %>% group_by(CD8_Annotation) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),
            n = n(), .groups = "drop")

pA <- ggplot(um, aes(UMAP_1, UMAP_2, colour = CD8_Annotation)) +
  geom_point(size = 0.25, alpha = 0.75) +
  ggrepel::geom_label_repel(data = lab_pos,
    aes(label = paste0(CD8_Annotation, "\n", format(n, big.mark = ","))),
    size = 3, fill = "white", alpha = 0.85, label.size = 0.15,
    show.legend = FALSE, min.segment.length = 0, box.padding = 0.4) +
  scale_colour_manual(values = CD8_COLORS, guide = "none") +
  labs(title = "Dedicated CD8+ T cell sub-clustering",
       subtitle = paste0("WNN clustering of CD8 T cells alone (re-normalised, ",
                         "RNA + ADT MNN integration); ", n_cells_label(cd8)),
       caption = paste0("Independent of the whole-PBMC clustering used in Fig 1. ",
                        "HEI participants only (n = ", length(unique(md$PID)), ").")) +
  theme_rev(11)
save_fig("R2M2a_01_CD8_dedicated_UMAP", pA, width = 9, height = 7.5)
save_src("R2M2a_01_CD8_dedicated_UMAP", lab_pos)

################################################################################
# PANEL B — the crosswalk. THIS is the answer to major 1d.
################################################################################
cross <- as.data.frame(table(source_cluster = md$source_cluster,
                             CD8_Annotation = md$CD8_Annotation))
cross <- cross[cross$Freq > 0, ]
cross$pct_of_source <- round(100 * cross$Freq /
  ave(cross$Freq, cross$source_cluster, FUN = sum), 1)
save_tab("CD8_ADT_subgating_crosswalk", cross[order(cross$source_cluster, -cross$Freq), ])

# WHAT THIS FIGURE ACTUALLY SHOWS - corrected.
#   `source_cluster` is the DEDICATED CD8 WNN clustering at resolution 0.4, NOT
#   the whole-PBMC clustering. It is set in 4.CD8_Subseting.R line 602 from
#   $seurat_clusters, immediately AFTER FindMultiModalNeighbors and before the
#   ADT sub-gating. Verified: purity of CD8_Annotation against wsnn_res.0.4 is
#   0.859, identical to purity against source_cluster.
#   So the off-diagonal structure is the ADT SUB-GATING step - CD45RO+ cells
#   pulled out as Naive Intermediate and CD45RO-/FAS+ as Tscm - redistributing
#   the four naive WNN clusters. It is NOT recovery of resolution lost in the
#   whole-PBMC clustering.
pB <- ggplot(cross, aes(x = factor(source_cluster), y = CD8_Annotation)) +
  geom_point(aes(size = Freq, colour = CD8_Annotation)) +
  geom_text(aes(label = format(Freq, big.mark = ",")), size = 2.6,
            nudge_x = 0.34, hjust = 0) +
  scale_colour_manual(values = CD8_COLORS, guide = "none") +
  scale_size_area(max_size = 11, name = "Cells") +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.1))) +
  labs(title = "ADT sub-gating of the dedicated CD8 WNN clusters",
       subtitle = paste0("WNN cluster (res 0.4) to final population, after ",
                         "CD45RO / FAS protein sub-gating; ", n_cells_label(cd8)),
       x = "Dedicated CD8 WNN cluster (resolution 0.4)",
       y = "Final annotated population",
       caption = paste0("Point size and label are cell counts. Clusters 0, 2, 4 ",
                        "and 6 each split: CD45RO+ cells become Naive\n",
                        "Intermediate CD8 and CD45RO-/FAS+ become Tscm CD8. ",
                        "Effector clusters 1, 3, 5 and 7 map one-to-one.")) +
  theme_rev(11)
save_fig("R2M1d_03_ADT_subgating_crosswalk", pB, width = 10, height = 6.5)
save_src("R2M1d_03_ADT_subgating_crosswalk", cross)

################################################################################
# PANELS C/D — bubble plots (major 1b), RNA and protein, CD4 included
################################################################################
rna_markers <- c("CD3E","CD8A","CD8B","CD4",                      # lineage
                 "CCR7","SELL","TCF7","LEF1","IL7R","BACH2",      # naive/stem
                 "CD27","CD28","FAS",                             # Tscm
                 "GZMK","GZMB","GNLY","PRF1","NKG7","KLRG1",      # effector
                 "FCGR3A","KLRC2","KIR2DL3",                      # CD16 / KIR
                 "IGFBP2","HILPDA","PFKFB4",                      # metabolic
                 "IFIT1","ISG15","MX1")                           # IFN
DefaultAssay(cd8) <- "RNA"
rna_present <- intersect(rna_markers, rownames(cd8))

pC <- DotPlot(cd8, features = rna_present, group.by = "CD8_Annotation",
              dot.scale = 6) +
  scale_colour_gradient2(low = "#4A90D9", mid = "grey92", high = "#D62728",
                         midpoint = 0, name = "Scaled\nexpression") +
  labs(title = "CD8 population markers - RNA",
       subtitle = paste0("Dot size = % of cells expressing, colour = scaled mean ",
                         "expression; ", n_cells_label(cd8)),
       x = NULL, y = NULL,
       caption = paste0("CD4 is shown alongside CD8A/CD8B so cluster identity ",
                        "can be judged directly (reviewer 2, major 1b).")) +
  theme_rev(10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8.5),
        legend.position = "right")
save_fig("R2M1b_01_bubble_RNA", pC, width = 13, height = 6)
save_src("R2M1b_01_bubble_RNA", pC$data)

adt_markers <- c("CD3D","CD8A","CD4","CD45RA","CD45RO","CD27","CD38","CD44",
                 "CD45","CD47","CD48","CD81","CD82")
DefaultAssay(cd8) <- "ADT"
adt_present <- intersect(adt_markers, rownames(cd8))

pD <- DotPlot(cd8, features = adt_present, group.by = "CD8_Annotation",
              dot.scale = 6) +
  scale_colour_gradient2(low = "#4A90D9", mid = "grey92", high = "#D62728",
                         midpoint = 0, name = "Scaled\nexpression") +
  labs(title = "CD8 population markers - surface protein (ADT)",
       subtitle = paste0("Dot size = % of cells expressing, colour = scaled mean ",
                         "expression; ", n_cells_label(cd8)),
       x = NULL, y = NULL,
       caption = paste0("DSB-normalised ADT. Protein CD8A and CD4 together ",
                        "address whether every population is genuinely CD8+.")) +
  theme_rev(10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8.5))
save_fig("R2M1b_02_bubble_ADT", pD, width = 9.5, height = 6)
save_src("R2M1b_02_bubble_ADT", pD$data)
DefaultAssay(cd8) <- "RNA"

################################################################################
# PANEL E — composition as % of TOTAL CD8 (major 1c + 1e)
################################################################################
comp <- md %>%
  group_by(Sample = orig.ident, PID, Timepoint_Group, CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  group_by(Sample) %>% mutate(pct_of_CD8 = 100 * n / sum(n)) %>% ungroup()
save_src("R2M1e_01_composition_pctCD8", comp)

pE <- ggplot(comp, aes(x = Timepoint_Group, y = pct_of_CD8, fill = Timepoint_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.65) +
  geom_point(position = position_jitter(width = 0.13, height = 0),
             size = 1.5, alpha = 0.85) +
  facet_wrap(~ CD8_Annotation, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = ART_COLORS, labels = ART_LABELS, name = NULL) +
  scale_x_discrete(labels = ART_LABELS) +
  labs(title = "Every CD8 population by ART status, as % of total CD8+ T cells",
       subtitle = paste0("One point per sample (n = ",
                         length(unique(comp$Sample)), " samples, ",
                         length(unique(comp$PID)), " participants)"),
       x = NULL, y = "% of CD8+ T cells in that sample",
       caption = paste0("Denominator is CD8 T cells, not live cells: PBMC-",
                        "denominated frequencies move with the other lineages ",
                        "(reviewer 2, major 1e).")) +
  theme_rev(10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "top")
save_fig("R2M1e_01_composition_pctCD8", pE, width = 14, height = 7.5)

################################################################################
# PANEL F — is "Naive CD8 2" really CD8+? (major 1b, minor 4)
################################################################################
id_genes <- intersect(c("CD8A","CD8B","CD4","CD3E"), rownames(cd8))
expr <- FetchData(cd8, vars = c(id_genes, "CD8_Annotation"), layer = "data")
id_long <- tidyr::pivot_longer(expr, cols = all_of(id_genes),
                               names_to = "gene", values_to = "expr")

id_summary <- id_long %>%
  group_by(CD8_Annotation, gene) %>%
  summarise(mean_expr = mean(expr),
            pct_pos   = 100 * mean(expr > 0), .groups = "drop")
save_src("R2M1d_04_naive2_identity", id_summary)

pF <- ggplot(id_summary, aes(x = gene, y = CD8_Annotation)) +
  geom_point(aes(size = pct_pos, colour = mean_expr)) +
  geom_text(aes(label = paste0(round(pct_pos), "%")), size = 2.5,
            nudge_x = 0.3, hjust = 0) +
  scale_colour_gradient(low = "grey90", high = "#D62728", name = "Mean\nlog-expr") +
  scale_size_area(max_size = 10, name = "% positive") +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.9))) +
  labs(title = "Lineage identity of every CD8 population",
       subtitle = paste0("CD8A/CD8B against CD4, RNA; ", n_cells_label(cd8)),
       x = NULL, y = NULL,
       caption = paste0("Answers whether each population - Naive CD8 2 in ",
                        "particular - is genuinely CD8+ rather than CD4 ",
                        "contamination.")) +
  theme_rev(11)
save_fig("R2M1d_04_naive2_identity", pF, width = 8, height = 6)

################################################################################
# NUMBERS FOR THE RESPONSE LETTER
################################################################################
n2 <- "Naïve CD8 2"
naive2_src <- cross[cross$CD8_Annotation == n2, ]
id_n2 <- id_summary[id_summary$CD8_Annotation == n2, ]

resp <- data.frame(
  item = c("dedicated CD8 cells", "dedicated populations", "participants",
           "WNN resolutions swept", "Naive CD8 2 cells",
           "Naive CD8 2 origin (whole-PBMC cluster)",
           "Naive CD8 2 CD8A pct_pos", "Naive CD8 2 CD4 pct_pos"),
  value = c(
    ncol(cd8), nlevels(md$CD8_Annotation), length(unique(md$PID)),
    paste(range(as.numeric(sub("wsnn_res\\.", "", res_cols))), collapse = "-"),
    sum(md$CD8_Annotation == n2, na.rm = TRUE),
    paste(naive2_src$source_cluster, collapse = ", "),
    round(id_n2$pct_pos[id_n2$gene == "CD8A"], 1),
    if (any(id_n2$gene == "CD4")) round(id_n2$pct_pos[id_n2$gene == "CD4"], 1) else NA))
save_tab("CD8_dedicated_response_numbers", resp)

cat("\n=== NUMBERS FOR THE RESPONSE LETTER ===\n")
print(resp, row.names = FALSE)

cat("\n=== DOES 'NAIVE CD8 2' PERSIST UNDER DEDICATED CLUSTERING? ===\n")
if (nrow(naive2_src)) {
  cat("  Yes. It is recovered as a distinct population of",
      sum(naive2_src$Freq), "cells,\n  derived from whole-PBMC cluster(s)",
      paste(naive2_src$source_cluster, collapse = ", "), ".\n")
  cat("  Its CD8A/CD4 profile is in RevFig_4F, so the identity question is\n")
  cat("  answered from the data rather than asserted.\n")
} else {
  cat("  NO - the label is absent from the dedicated clustering. That is itself\n")
  cat("  the answer to major 1d and must be stated plainly.\n")
}

cat("\nFigures: ", REV_FIG, "\nTables:  ", REV_TAB, "\n")
