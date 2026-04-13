################################################################################
# SUPPLEMENTARY FIGURES S4–S9: CD8 sub-cluster supporting data
#
# EDITS:
#   - Removed panel numbers (S4A, S5B, etc.) from all plot titles
#   - S7: common legend saved separately, removed from individual volcanos
#   - S8: all three panels combined in one row with shared legend underneath
#   - S9: increased all font sizes
#   - CSV outputs added for all supplementary data
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(pheatmap)
library(grid)
library(qs2)
library(ggrepel)
library(ggpubr)
library(SeuratExtend)
library(patchwork)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir     <- "~/Documents/CD8_Longitudinal"
saved_dir    <- file.path(base_dir, "saved_R_data")
manuscript   <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Manuscript/Fig 4-5"
supp_base    <- file.path(manuscript, "Supplementary")
analysis_dir <- file.path(manuscript, "analysis")

s4_dir <- file.path(supp_base, "S4_Annotation_Validation")
s5_dir <- file.path(supp_base, "S5_Cluster_Composition")
s6_dir <- file.path(supp_base, "S6_Effector_DGE")
s7_dir <- file.path(supp_base, "S7_Naive_Volcanos")
s8_dir <- file.path(supp_base, "S8_Pseudotime_Modules")
s9_dir <- file.path(supp_base, "S9_KIR_Innatelike")
csv_dir <- file.path(supp_base, "CSV_outputs")

for (d in c(s4_dir, s5_dir, s6_dir, s7_dir, s8_dir, s9_dir, csv_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Load object ──────────────────────────────────────────────────────────────
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("Loaded:", ncol(TARA_cd8), "cells\n")

# ── Shared aesthetics ────────────────────────────────────────────────────────
art_colors <- c(
  "PreART_Entry" = "#4A90D9", "PostART_Suppressed" = "#52B788", "PostART_Unsuppressed" = "#E76F51"
)

col_order_cd8 <- c(
  "Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
  "Tscm CD8", "Naïve Intermediate CD8",
  "Transitional Tem CD8", "TEM CD8", "TEMRA CD8",
  "KIR+ innate-like CD8", "MAIT-like Trm",
  "γδ1 T cell", "Naïve γδ1 T cell", "γδ2 T cell"
)

effector_clusters <- c("TEM CD8", "TEMRA CD8", "Transitional Tem CD8")

# Volcano gene categories (shared across S6, S7)
exhaustion_genes   <- c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1")
stemness_genes     <- c("TCF7","SELL","CCR7","LEF1","BACH2","BCL2","IL7R","S1PR1","KLF2")
cytotoxic_genes    <- c("GZMB","GNLY","PRF1","NKG7","GZMA","GZMH","FGFBP2")
stress_genes       <- c("HSPA1A","HSPA1B","DNAJB1","IFI27","IFI44L")
ifn_memory_genes   <- c("IFIT1","IFIT3","ISG15","MX1")
chemokine_genes    <- c("CCL3","CCL3L3","CCL4","CCL4L2","CCL5")
terminal_genes     <- c("ZEB2","PRDM1","TBX21","CX3CR1","S1PR5","ID2","EOMES")
activation_genes   <- c("HLA-DRA","CD38","FAS","CD69")

highlight_cols <- c(
  "Exhaustion"="#A32D2D","Stemness"="#185FA5","Cytotoxicity"="#3B6D11",
  "Stress response"="#BA7517","Type I IFN memory"="#0F6E56","Chemokines"="#993556",
  "Terminal diff."="#534AB7","Activation"="#D85A30","Other"="grey80"
)

hm_cols <- colorRampPalette(c("#F7FCF5","#C7E9C0","#74C476","#31A354","#006D2C"))(100)
scale_01 <- function(mat) {
  t(apply(mat, 1, function(x) { r <- max(x)-min(x); if(r==0) rep(0.5,length(x)) else (x-min(x))/r }))
}

# Helper: assign highlight category to DGE genes
assign_highlight <- function(dge) {
  dge$highlight <- "Other"
  dge$highlight[dge$gene %in% exhaustion_genes]  <- "Exhaustion"
  dge$highlight[dge$gene %in% stemness_genes]    <- "Stemness"
  dge$highlight[dge$gene %in% cytotoxic_genes]   <- "Cytotoxicity"
  dge$highlight[dge$gene %in% stress_genes]      <- "Stress response"
  dge$highlight[dge$gene %in% ifn_memory_genes]  <- "Type I IFN memory"
  dge$highlight[dge$gene %in% chemokine_genes]   <- "Chemokines"
  dge$highlight[dge$gene %in% terminal_genes]    <- "Terminal diff."
  dge$highlight[dge$gene %in% activation_genes]  <- "Activation"
  dge$highlight <- factor(dge$highlight, levels = names(highlight_cols))
  dge
}


################################################################################
# S4: ANNOTATION VALIDATION
################################################################################
cat("\n=== S4: Annotation validation ===\n")

# ── Shared group colors (same order as Fig 1) ────────────────────────────────
s4_group_colors <- c(
  "Naive/Quiescence"      = "#74C2E1",
  "Memory/Survival"       = "#66BB6A",
  "Co-stimulation/Homing" = "#AED581",
  "Effector/Cytotoxicity" = "#EF5350",
  "Terminal Diff."        = "#AB47BC",
  "Exhaustion"            = "#A1887F",
  "Activation"            = "#FFD54F",
  "Innate-like/NK"        = "#F06292",
  "Tissue Residency"      = "#9C755F",
  "TCR Identity"          = "#9C6FD6",
  "MAIT"                  = "#26A69A"
)

get_gaps <- function(annot_df) {
  grps <- as.character(annot_df$Group)
  which(grps[-length(grps)] != grps[-1])
}

annot_colors <- list(Group = s4_group_colors)

# ── S4A: ADT heatmap ────────────────────────────────────────────────────────
cat("  S4A: ADT heatmap...\n")

adt_heatmap_markers <- c(
  "CD45RA"    = "Naive/Quiescence",     "SELL"      = "Naive/Quiescence",
  "CD7"       = "Naive/Quiescence",     "IL7R"      = "Naive/Quiescence",
  "CD45RO"    = "Memory/Survival",      "CD44"      = "Memory/Survival",
  "CD27"      = "Co-stimulation/Homing","CD28"      = "Co-stimulation/Homing",
  "ICOS"      = "Co-stimulation/Homing","ITGB7"     = "Co-stimulation/Homing",
  "B3GAT1"    = "Effector/Cytotoxicity","KLRG1"     = "Effector/Cytotoxicity",
  "CX3CR1"    = "Effector/Cytotoxicity",
  "TIGIT"     = "Exhaustion",           "PDCD1"     = "Exhaustion",
  "LAG3"      = "Exhaustion",           "ENTPD1"    = "Exhaustion",
  "NT5E"      = "Exhaustion",
  "CD38"      = "Activation",           "CD69"      = "Activation",
  "FAS"       = "Activation",
  "CXCR3"     = "Tissue Residency",     "ITGB7"     = "Tissue Residency",
  "ITGA1"     = "Tissue Residency",
  "KIR3DL1"   = "Innate-like/NK",       "NCAM1"     = "Innate-like/NK",
  "FCGR3A"    = "Innate-like/NK",       "SIGLEC7"   = "Innate-like/NK",
  "KLRD1"     = "Innate-like/NK",       "KLRB1"     = "Innate-like/NK",
  "TCR-AB"    = "TCR Identity",         "TCR-vA7.2" = "TCR Identity",
  "TCR-vD2"   = "TCR Identity"
)

adt_heatmap_markers <- adt_heatmap_markers[!duplicated(names(adt_heatmap_markers))]

DefaultAssay(TARA_cd8) <- "ADT"
adt_feats <- names(adt_heatmap_markers)[names(adt_heatmap_markers) %in% rownames(TARA_cd8)]

avg_adt <- AverageExpression(TARA_cd8, assays = "ADT", features = adt_feats,
                             group.by = "CD8_Annotation", slot = "data")$ADT
avg_adt_s <- scale_01(as.matrix(log10(avg_adt + 1)))
colnames(avg_adt_s) <- gsub("^g ", "", colnames(avg_adt_s))
co <- col_order_cd8[col_order_cd8 %in% colnames(avg_adt_s)]
avg_adt_s <- avg_adt_s[, co]

adt_row_order <- adt_feats[adt_feats %in% rownames(avg_adt_s)]
avg_adt_s <- avg_adt_s[adt_row_order, , drop = FALSE]

adt_annot <- data.frame(
  Group = factor(adt_heatmap_markers[rownames(avg_adt_s)], levels = names(s4_group_colors)),
  row.names = rownames(avg_adt_s)
)

adt_mat <- matrix(as.numeric(avg_adt_s), nrow = nrow(avg_adt_s), dimnames = dimnames(avg_adt_s))

p_adt <- pheatmap(
  adt_mat, cluster_rows = FALSE, cluster_cols = FALSE, color = hm_cols,
  border_color = "white", cellwidth = 48, cellheight = 22, fontsize = 16,
  fontsize_row = 14, fontsize_col = 13, angle_col = 45,
  annotation_row = adt_annot, annotation_colors = annot_colors,
  annotation_names_row = FALSE, annotation_legend = TRUE,
  gaps_row = get_gaps(adt_annot),
  main = "Surface protein (ADT) expression — 13 CD8 sub-clusters",
  silent = TRUE
)

png(file.path(s4_dir, "S4A_ADT_heatmap.png"), width = 26, height = 20, units = "in", res = 300, bg = "white")
grid.draw(p_adt$gtable)
dev.off()

# CSV: ADT heatmap values
write.csv(as.data.frame(adt_mat), file.path(csv_dir, "S4A_ADT_heatmap_scaled_values.csv"))

# ── S4B: RNA heatmap ────────────────────────────────────────────────────────
cat("  S4B: RNA heatmap...\n")

rna_heatmap_markers <- c(
  "CCR7"    = "Naive/Quiescence",      "SELL"    = "Naive/Quiescence",
  "TCF7"    = "Naive/Quiescence",      "LEF1"    = "Naive/Quiescence",
  "KLF2"    = "Naive/Quiescence",      "S1PR1"   = "Naive/Quiescence",
  "BACH2"   = "Memory/Survival",       "BCL2"    = "Memory/Survival",
  "IL7R"    = "Memory/Survival",       "ID3"     = "Memory/Survival",
  "HLA-DRA" = "Activation",            "CD38"    = "Activation",
  "MKI67"   = "Activation",
  "GZMB"    = "Effector/Cytotoxicity", "GNLY"    = "Effector/Cytotoxicity",
  "PRF1"    = "Effector/Cytotoxicity", "NKG7"    = "Effector/Cytotoxicity",
  "GZMA"    = "Effector/Cytotoxicity", "FGFBP2"  = "Effector/Cytotoxicity",
  "TBX21"   = "Terminal Diff.",        "ZEB2"    = "Terminal Diff.",
  "PRDM1"   = "Terminal Diff.",        "CX3CR1"  = "Terminal Diff.",
  "EOMES"   = "Terminal Diff.",
  "TOX"     = "Exhaustion",            "PDCD1"   = "Exhaustion",
  "TIGIT"   = "Exhaustion",            "HAVCR2"  = "Exhaustion",
  "LAG3"    = "Exhaustion",
  "TRDV1"   = "TCR Identity",          "TRGV9"   = "TCR Identity",
  "TRDC"    = "TCR Identity",
  "SLC4A10" = "MAIT",                  "ZBTB16"  = "MAIT"
)

DefaultAssay(TARA_cd8) <- "RNA"
rna_feats <- names(rna_heatmap_markers)[names(rna_heatmap_markers) %in% rownames(TARA_cd8)]

avg_rna <- AverageExpression(TARA_cd8, assays = "RNA", features = rna_feats,
                             group.by = "CD8_Annotation", slot = "data")$RNA
avg_rna_s <- scale_01(as.matrix(avg_rna))
colnames(avg_rna_s) <- gsub("^g ", "", colnames(avg_rna_s))
avg_rna_s <- avg_rna_s[, col_order_cd8[col_order_cd8 %in% colnames(avg_rna_s)]]

rna_row_order <- rna_feats[rna_feats %in% rownames(avg_rna_s)]
avg_rna_s <- avg_rna_s[rna_row_order, , drop = FALSE]

rna_annot <- data.frame(
  Group = factor(rna_heatmap_markers[rownames(avg_rna_s)], levels = names(s4_group_colors)),
  row.names = rownames(avg_rna_s)
)

rna_mat <- matrix(as.numeric(avg_rna_s), nrow = nrow(avg_rna_s), dimnames = dimnames(avg_rna_s))

p_rna <- pheatmap(
  rna_mat, cluster_rows = FALSE, cluster_cols = FALSE, color = hm_cols,
  border_color = "white", cellwidth = 48, cellheight = 22, fontsize = 16,
  fontsize_row = 14, fontsize_col = 13, angle_col = 45,
  annotation_row = rna_annot, annotation_colors = annot_colors,
  annotation_names_row = FALSE, annotation_legend = TRUE,
  gaps_row = get_gaps(rna_annot),
  main = "mRNA expression — 13 CD8 sub-clusters",
  silent = TRUE
)

png(file.path(s4_dir, "S4B_RNA_heatmap.png"), width = 26, height = 20, units = "in", res = 300, bg = "white")
grid.draw(p_rna$gtable)
dev.off()

# CSV: RNA heatmap values
write.csv(as.data.frame(rna_mat), file.path(csv_dir, "S4B_RNA_heatmap_scaled_values.csv"))

# ── Standalone horizontal group legend ──
cat("  Shared group legend...\n")

legend_labels <- c(
  "Naive/\nQuiescence", "Memory/\nSurvival", "Co-stim./\nHoming",
  "Effector/\nCytotoxicity", "Terminal\nDiff.", "Exhaustion",
  "Activation", "Innate-like/\nNK", "Tissue\nResidency",
  "TCR\nIdentity", "MAIT"
)

legend_df <- data.frame(
  Group = factor(names(s4_group_colors), levels = names(s4_group_colors)),
  y     = seq(1, by = 0.7, length.out = length(s4_group_colors)),
  label = legend_labels
)

p_legend_only <- ggplot(legend_df, aes(x = 1, y = y)) +
  geom_point(aes(fill = Group), shape = 22, size = 10, stroke = 0) +
  geom_text(aes(x = 0.82, label = label), size = 5, fontface = "bold",
            lineheight = 0.8, vjust = 0.5) +
  scale_fill_manual(values = s4_group_colors) +
  scale_x_continuous(limits = c(0.6, 1.1)) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(0, 2, 0, 2))

ggsave(file.path(s4_dir, "S4_Group_Legend.png"),
       p_legend_only, width = 16, height = 1.6, dpi = 300, bg = "white")

# ── S4C: ADT gating scatter plots ───────────────────────────────────────────
cat("  S4C: ADT gating...\n")

naive_clusters <- c("Naïve CD8", "Naïve CD8 2", "Naïve CD8 3",
                    "Tscm CD8", "Naïve Intermediate CD8")
naive_mask <- TARA_cd8$CD8_Annotation %in% naive_clusters
cat("    Naïve-lineage cells:", sum(naive_mask), "\n")

if (sum(naive_mask) > 0) {
  naive_cells <- colnames(TARA_cd8)[naive_mask]
  DefaultAssay(TARA_cd8) <- "ADT"
  adt_data <- GetAssayData(TARA_cd8, slot = "data")
  
  gate_plot_df <- data.frame(
    CD45RO     = as.numeric(adt_data["CD45RO", naive_cells]),
    FAS        = as.numeric(adt_data["FAS", naive_cells]),
    annotation = TARA_cd8$CD8_Annotation[naive_cells],
    stringsAsFactors = FALSE
  )
  
  gate_plot_df$plot_group <- dplyr::case_when(
    gate_plot_df$annotation == "Tscm CD8"               ~ "Tscm CD8",
    gate_plot_df$annotation == "Naïve Intermediate CD8"  ~ "Naïve Intermediate CD8",
    TRUE                                                 ~ "Naïve (ungated)"
  )
  gate_plot_df$plot_group <- factor(gate_plot_df$plot_group,
                                    levels = c("Naïve (ungated)", "Tscm CD8", "Naïve Intermediate CD8"))
  
  if ("source_cluster" %in% colnames(TARA_cd8@meta.data)) {
    gate_plot_df$source <- TARA_cd8$source_cluster[naive_cells]
    cat("    source_cluster values:", paste(sort(unique(gate_plot_df$source)), collapse = ", "), "\n")
    gate_plot_df$facet <- dplyr::case_when(
      !is.na(gate_plot_df$source) ~ paste0("Source cluster ", gate_plot_df$source),
      gate_plot_df$annotation == "Naïve CD8"   ~ "Naïve CD8",
      gate_plot_df$annotation == "Naïve CD8 2" ~ "Naïve CD8 2",
      gate_plot_df$annotation == "Naïve CD8 3" ~ "Naïve CD8 3",
      TRUE ~ "Other"
    )
  } else {
    gate_plot_df$facet <- "All naïve clusters"
  }
  gate_plot_df <- gate_plot_df %>% filter(facet != "Other")
  
  gate_cols <- c(
    "Tscm CD8"               = "#E31A1C",
    "Naïve Intermediate CD8" = "#1F78B4"
  )
  
  fg_df <- gate_plot_df %>% filter(plot_group != "Naïve (ungated)")
  cat("    Gated cells to plot:", nrow(fg_df), "\n")
  
  p_s4c <- ggplot(fg_df, aes(x = FAS, y = CD45RO, color = plot_group)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = gate_cols, name = "Population",
                       labels = c("Tscm CD8" = "Tscm (CD45RO\u2212 FAS+)",
                                  "Naïve Intermediate CD8" = "Intermediate (CD45RO+)")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.7) +
    facet_wrap(~ facet, nrow = 1) +
    labs(title = "ADT sub-gating of naïve clusters",
         subtitle = "Tscm = CD45RO\u2212 FAS+  |  Intermediate = CD45RO+  |  True naïve = CD45RO\u2212 FAS\u2212",
         x = "FAS/CD95 (DSB)", y = "CD45RO (DSB)") +
    theme_cowplot(font_size = 22) +
    theme(
      strip.background  = element_rect(fill = "grey95", color = NA),
      strip.text         = element_text(size = 20, face = "bold"),
      axis.text          = element_text(size = 18),
      axis.title         = element_text(size = 20),
      plot.title         = element_text(size = 24, face = "bold"),
      plot.subtitle      = element_text(size = 18),
      legend.position    = "bottom",
      legend.text        = element_text(size = 20),
      legend.title       = element_text(size = 22, face = "bold"),
      legend.key.size    = unit(1.5, "cm"),
      plot.background    = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 8, alpha = 1), nrow = 1))
  
  ggsave(file.path(s4_dir, "S4C_ADT_gating.png"), plot = p_s4c,
         width = 24, height = 10, dpi = 300, bg = "white")
  
  # CSV: gating data
  write.csv(fg_df, file.path(csv_dir, "S4C_ADT_gating_data.csv"), row.names = FALSE)
}
cat("  S4 done.\n")

################################################################################
# S5: CLUSTER COMPOSITION
################################################################################
cat("\n=== S5: Cluster composition ===\n")

comp_mat <- as.data.frame.matrix(
  prop.table(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 1))
comp_mat <- comp_mat[col_order_cd8[col_order_cd8 %in% rownames(comp_mat)], , drop = FALSE]
colnames(comp_mat) <- c("Pre-ART","Suppressed","Unsuppressed")[
  match(colnames(comp_mat), c("PreART_Entry","PostART_Suppressed","PostART_Unsuppressed"))]

display_mat <- matrix(paste0(round(as.matrix(comp_mat)*100,1),"%"),
                      nrow=nrow(comp_mat), dimnames=dimnames(comp_mat))

prop_colors <- colorRampPalette(c("#FFFFFF","#D6E8F7","#7FB3D8","#2171B5","#08306B"))(100)

png(file.path(s5_dir, "S5A_composition_by_ART.png"), width = 14, height = 16, units = "in", res = 300, bg = "white")
pheatmap(as.matrix(comp_mat), cluster_rows=FALSE, cluster_cols=FALSE, color=prop_colors,
         border_color="white", cellwidth=100, cellheight=35, fontsize=16, fontsize_row=14,
         fontsize_col=16, angle_col=0, display_numbers=display_mat, fontsize_number=14,
         number_color=ifelse(as.matrix(comp_mat)>0.45,"white","black"),
         main="Cluster composition by ART status\n(% of each cluster from each condition)")
dev.off()

# CSV: composition by ART
write.csv(comp_mat, file.path(csv_dir, "S5A_composition_by_ART.csv"))

comp_col <- as.data.frame.matrix(
  prop.table(table(TARA_cd8$CD8_Annotation, TARA_cd8$Timepoint_Group), margin = 2))
comp_col <- comp_col[col_order_cd8[col_order_cd8 %in% rownames(comp_col)], , drop = FALSE]
colnames(comp_col) <- c("Pre-ART","Suppressed","Unsuppressed")[
  match(colnames(comp_col), c("PreART_Entry","PostART_Suppressed","PostART_Unsuppressed"))]

display_col <- matrix(paste0(round(as.matrix(comp_col)*100,1),"%"),
                      nrow=nrow(comp_col), dimnames=dimnames(comp_col))

png(file.path(s5_dir, "S5B_condition_distribution.png"), width = 14, height = 16, units = "in", res = 300, bg = "white")
pheatmap(as.matrix(comp_col), cluster_rows=FALSE, cluster_cols=FALSE, color=prop_colors,
         border_color="white", cellwidth=100, cellheight=35, fontsize=16, fontsize_row=14,
         fontsize_col=16, angle_col=0, display_numbers=display_col, fontsize_number=14,
         number_color=ifelse(as.matrix(comp_col)>0.15,"white","black"),
         main="Distribution of each condition across clusters")
dev.off()

# CSV: condition distribution
write.csv(comp_col, file.path(csv_dir, "S5B_condition_distribution.csv"))

cat("  S5 done.\n")

################################################################################
# S6: EFFECTOR DGE — TEM + TEMRA volcanos + ADT dot plot
################################################################################
cat("\n=== S6: Effector DGE ===\n")

genes_to_label_all <- c(exhaustion_genes, stemness_genes, cytotoxic_genes, stress_genes,
                        ifn_memory_genes, chemokine_genes, terminal_genes, activation_genes)

volc_files <- list(
  list(label="S6A", title="TEM CD8", file="DGE_MAST_Expanding_TEMCD8_Sup_vs_Unsup.csv"),
  list(label="S6B", title="TEMRA CD8", file="DGE_MAST_Expanding_TEMRACD8_Sup_vs_Unsup.csv")
)

for (vf in volc_files) {
  fpath <- file.path(analysis_dir, "03_DGE", vf$file)
  if (!file.exists(fpath)) { cat("    SKIP", vf$label, "\n"); next }
  cat("    ", vf$label, ":", vf$title, "...\n")
  
  dge <- read.csv(fpath, row.names = 1)
  dge$gene <- rownames(dge)
  dge$neg_log10_padj <- -log10(dge$p_val_adj + 1e-300)
  dge <- assign_highlight(dge)
  
  gp <- genes_to_label_all[genes_to_label_all %in% dge$gene]
  dge$label <- ""
  dge$label[dge$gene %in% gp & dge$p_val_adj < 0.05] <- dge$gene[dge$gene %in% gp & dge$p_val_adj < 0.05]
  
  x_lim <- quantile(abs(dge$avg_log2FC[is.finite(dge$avg_log2FC)]), 0.995) * 1.15
  y_cap <- quantile(dge$neg_log10_padj[is.finite(dge$neg_log10_padj)], 0.995) * 1.10
  
  # TITLE WITHOUT panel number
  p <- ggplot(dge, aes(x = pmax(pmin(avg_log2FC, x_lim*0.98), -x_lim*0.98),
                       y = pmin(neg_log10_padj, y_cap), color = highlight)) +
    geom_point(data = dge %>% filter(highlight == "Other"), size = 2, alpha = 0.3) +
    geom_point(data = dge %>% filter(highlight != "Other"), size = 7, alpha = 0.85) +
    geom_label_repel(data = dge %>% filter(label != ""),
                     aes(label = label, fill = highlight), color = "white",
                     size = 8, fontface = "bold.italic", max.overlaps = 50,
                     segment.size = 0.4, segment.color = "grey40",
                     box.padding = 0.5, label.size = 0.2, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = highlight_cols, name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = highlight_cols, guide = "none") +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~italic(p)),
         title = paste0(vf$title, " — expanding clones: suppressed vs unsuppressed")) +
    theme_cowplot(font_size = 20) +
    theme(legend.position = "right", legend.text = element_text(size = 20),
          legend.title = element_text(size = 22, face = "bold"), legend.key.size = unit(1.3, "cm"),
          axis.text = element_text(size = 16), axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(10, 10, 80, 10)) +
    guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
  
  arrow_right <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
              gp = gpar(col = "#52B788", lwd = 8)),
    textGrob("Higher in suppressed", x = unit(0.74, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#52B788", fontsize = 24, fontface = "bold"))
  )
  arrow_left <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 8)),
    textGrob("Higher in unsuppressed", x = unit(0.26, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 24, fontface = "bold"))
  )
  
  p <- p +
    annotation_custom(grob = arrow_right, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.24, ymax = -y_cap * 0.10) +
    annotation_custom(grob = arrow_left, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.24, ymax = -y_cap * 0.10) +
    coord_cartesian(ylim = c(0, y_cap), xlim = c(-x_lim, x_lim), clip = "off")
  
  ggsave(file.path(s6_dir, paste0(vf$label, "_Volcano_", gsub(" ", "_", vf$title), ".png")),
         plot = p, width = 14, height = 11, dpi = 300, bg = "white")
  
  # CSV: copy DGE source file
  file.copy(fpath, file.path(csv_dir, paste0(vf$label, "_DGE_", gsub(" ", "_", vf$title), ".csv")),
            overwrite = TRUE)
}

# ── S6C: ADT dot plot ────────────────────────────────────────────────────────
cat("    S6C: ADT dot plot...\n")
adt_prots <- c("TIGIT","PDCD1","LAG3","FAS","CD38","ICOS","CD27","CD28",
               "B3GAT1","KLRG1","CX3CR1","CD45RA","CD45RO","SELL","KIR3DL1","KLRD1")
adt_prots <- adt_prots[adt_prots %in% rownames(TARA_cd8[["ADT"]])]

expand_eff <- subset(TARA_cd8, CD8_Annotation %in% effector_clusters &
                       !is.na(cloneSize) & cloneSize != "Single (0 < X <= 1)")
DefaultAssay(expand_eff) <- "ADT"

adt_list <- lapply(adt_prots, function(prot) {
  expr <- FetchData(expand_eff, vars = prot)[, 1]
  data.frame(protein = prot, Timepoint_Group = expand_eff$Timepoint_Group, expression = expr) %>%
    group_by(protein, Timepoint_Group) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE),
              pct_pos = mean(expression > 0, na.rm = TRUE) * 100, .groups = "drop")
})
adt_df <- bind_rows(adt_list)
adt_df$protein <- factor(adt_df$protein, levels = adt_prots)
adt_df$Condition <- dplyr::recode(adt_df$Timepoint_Group, "PreART_Entry"="Pre-ART",
                                  "PostART_Suppressed"="Suppressed", "PostART_Unsuppressed"="Unsuppressed")
adt_df$Condition <- factor(adt_df$Condition, levels = c("Pre-ART","Suppressed","Unsuppressed"))

# TITLE WITHOUT panel number
p_s6c <- ggplot(adt_df, aes(x = protein, y = Condition)) +
  geom_point(aes(size = pct_pos, color = mean_expr)) +
  scale_size_continuous(name = "% expressing", range = c(2, 12), limits = c(0, 100)) +
  scale_color_viridis_c(option = "magma", name = "Mean expression\n(DSB)", direction = -1) +
  labs(x = NULL, y = NULL, title = "Surface protein — expanding effector clones") +
  theme_cowplot(font_size = 18) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

ggsave(file.path(s6_dir, "S6C_ADT_dotplot.png"), plot = p_s6c,
       width = 14, height = 6, dpi = 300, bg = "white")

# CSV: ADT dot plot data
write.csv(adt_df, file.path(csv_dir, "S6C_ADT_dotplot_data.csv"), row.names = FALSE)

cat("  S6 done.\n")

################################################################################
# S7: NAÏVE PAIRWISE VOLCANOS — common legend saved separately
################################################################################
cat("\n=== S7: Naïve pairwise volcanos ===\n")

naive_cols <- c("Quiescence / homing"="#185FA5","Hypoxia / metabolism"="#A32D2D",
                "Type I IFN"="#0F6E56","Survival / quiescence"="#7B68EE",
                "Activation"="#D85A30","Other"="grey80")

quiescence_genes  <- c("KLF2","S1PR1","SELL","BACH2","BCL2","IL7R","FOXP1","CCR7","LEF1","TCF7")
hypoxia_genes     <- c("IGFBP2","HILPDA","PFKFB4","EGLN3","GBP5","AK4","HIF1A")
ifn_genes_naive   <- c("MX1","IFIT1","ISG15","IFI44L","STAT1","IRF7","IFIT3")
survival_genes    <- c("BACH2","BCL2","FOXP1","SATB1","ID3")
activation_naive  <- c("HLA-DRA","CD38","ICOS","CD69","FAS")

naive_label_genes <- unique(c(quiescence_genes, hypoxia_genes, ifn_genes_naive,
                              survival_genes, activation_naive))

naive_comps <- list(
  list(label="S7A", title="Naïve CD8 vs Naïve CD8 2", id1="Naïve CD8", id2="Naïve CD8 2",
       file="DGE_MAST_RNA_NaveCD8_vs_NaveCD82.csv"),
  list(label="S7B", title="Naïve CD8 vs Naïve CD8 3", id1="Naïve CD8", id2="Naïve CD8 3",
       file="DGE_MAST_RNA_NaveCD8_vs_NaveCD83.csv"),
  list(label="S7C", title="Naïve CD8 2 vs Naïve CD8 3", id1="Naïve CD8 2", id2="Naïve CD8 3",
       file="DGE_MAST_RNA_NaveCD82_vs_NaveCD83.csv")
)

# Build all S7 plots in a list (legend removed from individual plots)
s7_plots <- list()

for (nc in naive_comps) {
  fpath <- file.path(analysis_dir, "03_DGE", nc$file)
  if (!file.exists(fpath)) { cat("    SKIP", nc$label, "\n"); next }
  cat("    ", nc$label, ":", nc$title, "...\n")
  
  dge <- read.csv(fpath, row.names = 1)
  dge$gene <- rownames(dge)
  dge$neg_log10_padj <- -log10(dge$p_val_adj + 1e-300)
  
  dge$highlight <- "Other"
  dge$highlight[dge$gene %in% quiescence_genes] <- "Quiescence / homing"
  dge$highlight[dge$gene %in% hypoxia_genes]    <- "Hypoxia / metabolism"
  dge$highlight[dge$gene %in% ifn_genes_naive]  <- "Type I IFN"
  dge$highlight[dge$gene %in% survival_genes]   <- "Survival / quiescence"
  dge$highlight[dge$gene %in% activation_naive] <- "Activation"
  dge$highlight <- factor(dge$highlight, levels = names(naive_cols))
  
  gp <- naive_label_genes[naive_label_genes %in% dge$gene]
  dge$label <- ""
  dge$label[dge$gene %in% gp & dge$p_val_adj < 0.05] <- dge$gene[dge$gene %in% gp & dge$p_val_adj < 0.05]
  dge$label_ns <- ""
  dge$label_ns[dge$gene %in% gp & dge$p_val_adj >= 0.05] <- dge$gene[dge$gene %in% gp & dge$p_val_adj >= 0.05]
  
  x_lim <- quantile(abs(dge$avg_log2FC[is.finite(dge$avg_log2FC)]), 0.995) * 1.15
  y_cap <- quantile(dge$neg_log10_padj[is.finite(dge$neg_log10_padj)], 0.995) * 1.10
  
  # TITLE WITHOUT panel number; LEGEND REMOVED (saved separately below)
  p <- ggplot(dge, aes(x = pmax(pmin(avg_log2FC, x_lim*0.98), -x_lim*0.98),
                       y = pmin(neg_log10_padj, y_cap), color = highlight)) +
    geom_point(data = dge %>% filter(highlight == "Other"), size = 2, alpha = 0.3) +
    geom_point(data = dge %>% filter(highlight != "Other"), size = 7, alpha = 0.85) +
    geom_label_repel(data = dge %>% filter(label != ""),
                     aes(label = label, fill = highlight), color = "white",
                     size = 8, fontface = "bold.italic", max.overlaps = 50,
                     segment.size = 0.4, segment.color = "grey40",
                     box.padding = 0.5, label.size = 0.2, show.legend = FALSE) +
    ggrepel::geom_text_repel(data = dge %>% filter(label_ns != ""),
                             aes(label = label_ns), size = 5, fontface = "italic", color = "grey40",
                             max.overlaps = 20, segment.size = 0.2, segment.color = "grey60",
                             segment.linetype = "dashed", box.padding = 0.3, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = naive_cols, name = "Gene category", drop = FALSE) +
    scale_fill_manual(values = naive_cols, guide = "none") +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~italic(p)),
         title = nc$title) +
    theme_cowplot(font_size = 20) +
    theme(legend.position = "none",
          axis.text = element_text(size = 16), axis.title = element_text(size = 18),
          plot.title = element_text(size = 20, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(10, 10, 80, 10)) +
    guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
  
  # Directional arrows
  arrow_right_s7 <- grobTree(
    linesGrob(x = unit(c(0.52, 0.95), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
              gp = gpar(col = "#185FA5", lwd = 8)),
    textGrob(paste0("Higher in ", nc$id1), x = unit(0.74, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#185FA5", fontsize = 22, fontface = "bold"))
  )
  arrow_left_s7 <- grobTree(
    linesGrob(x = unit(c(0.48, 0.05), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
              gp = gpar(col = "#A32D2D", lwd = 8)),
    textGrob(paste0("Higher in ", nc$id2), x = unit(0.26, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#A32D2D", fontsize = 22, fontface = "bold"))
  )
  
  p <- p +
    annotation_custom(grob = arrow_right_s7, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.24, ymax = -y_cap * 0.10) +
    annotation_custom(grob = arrow_left_s7, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.24, ymax = -y_cap * 0.10) +
    coord_cartesian(ylim = c(0, y_cap), xlim = c(-x_lim, x_lim), clip = "off")
  
  # Save individual plot (no legend)
  ggsave(file.path(s7_dir, paste0(nc$label, "_Volcano.png")),
         plot = p, width = 14, height = 11, dpi = 300, bg = "white")
  
  s7_plots[[nc$label]] <- p
  
  # CSV: copy DGE source file
  file.copy(fpath, file.path(csv_dir, paste0(nc$label, "_DGE_Naive_pairwise.csv")),
            overwrite = TRUE)
}

# ── S7 COMMON LEGEND (saved as separate file, horizontal, for placing under plots) ──
cat("  S7: Saving common legend...\n")

p_s7_legend_src <- ggplot(data.frame(x = 1:5, y = 1:5,
                                     cat = factor(names(naive_cols)[1:5], levels = names(naive_cols))),
                          aes(x = x, y = y, color = cat)) +
  geom_point(size = 8) +
  scale_color_manual(values = naive_cols, name = "Gene category", drop = FALSE) +
  theme_void() +
  theme(
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 22),
    legend.title      = element_text(size = 24, face = "bold"),
    legend.key.size   = unit(1.5, "cm"),
    legend.spacing.x  = unit(0.8, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 10, alpha = 1), nrow = 1))

s7_legend <- cowplot::get_legend(p_s7_legend_src)
ggsave(file.path(s7_dir, "S7_common_legend.png"),
       plot = as_ggplot(s7_legend), width = 20, height = 2.5, dpi = 300, bg = "white")

cat("  S7 done.\n")

################################################################################
# S8: PSEUDOTIME MODULES — all 3 in one row, common legend underneath
################################################################################
cat("\n=== S8: Additional pseudotime modules ===\n")

m3_pt_path <- file.path(analysis_dir, "07_trajectory", "monocle3", "Monocle3_pseudotime_per_cell.csv")
module_path <- file.path(analysis_dir, "10_module_scores", "ModuleScores_per_cell.csv")

if (file.exists(m3_pt_path) & file.exists(module_path)) {
  m3_pt <- read.csv(m3_pt_path)
  if (!"cell" %in% colnames(m3_pt)) m3_pt$cell <- m3_pt[,1]
  mod_scores <- read.csv(module_path, row.names = 1)
  
  common <- intersect(m3_pt$cell, rownames(mod_scores))
  pt_mod <- data.frame(
    pseudotime = m3_pt$monocle3_pseudotime[match(common, m3_pt$cell)],
    Timepoint  = mod_scores[common, "Timepoint_Group"]
  )
  if ("Cytotoxicity" %in% colnames(mod_scores)) pt_mod$Cytotoxicity <- mod_scores[common, "Cytotoxicity"]
  if ("Inflammatory_Chemokines" %in% colnames(mod_scores)) pt_mod$Inflammatory_Chemokines <- mod_scores[common, "Inflammatory_Chemokines"]
  if ("TypeI_IFN_Memory" %in% colnames(mod_scores)) pt_mod$TypeI_IFN_Memory <- mod_scores[common, "TypeI_IFN_Memory"]
  pt_mod <- pt_mod %>% filter(!is.na(pseudotime) & is.finite(pseudotime) & !is.na(Timepoint))
  pt_mod$Timepoint <- factor(pt_mod$Timepoint, levels = names(art_colors))
  
  # CSV: pseudotime + module scores
  write.csv(pt_mod, file.path(csv_dir, "S8_pseudotime_module_scores.csv"), row.names = FALSE)
  
  s8_panels <- list(
    list(col = "Cytotoxicity", ylab = "Cytotoxicity score",
         title = "Cytotoxicity along differentiation"),
    list(col = "Inflammatory_Chemokines", ylab = "Inflammatory chemokines score",
         title = "Inflammatory chemokines along differentiation"),
    list(col = "TypeI_IFN_Memory", ylab = "Type I IFN memory score",
         title = "Type I IFN memory along differentiation")
  )
  
  s8_plot_list <- list()
  
  for (i in seq_along(s8_panels)) {
    sp <- s8_panels[[i]]
    if (!sp$col %in% colnames(pt_mod)) { cat("    SKIP", sp$col, "\n"); next }
    cat("    ", sp$title, "...\n")
    
    p_s8 <- ggplot(pt_mod, aes(x = pseudotime, y = .data[[sp$col]], color = Timepoint)) +
      geom_point(size = 0.15, alpha = 0.1) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.8, alpha = 0.2, span = 0.4) +
      scale_color_manual(values = art_colors, name = "ART Status",
                         labels = c("Pre-ART","Suppressed","Unsuppressed")) +
      labs(x = "Monocle3 pseudotime", y = sp$ylab, title = sp$title) +
      theme_cowplot(font_size = 20) +
      theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
            plot.title = element_text(size = 18, face = "bold"),
            legend.position = "none",
            plot.background = element_rect(fill = "white", color = NA))
    
    s8_plot_list[[i]] <- p_s8
  }
  
  # Build shared legend
  p_s8_legend_src <- ggplot(pt_mod, aes(x = pseudotime, y = Cytotoxicity, color = Timepoint)) +
    geom_point(size = 5) +
    scale_color_manual(values = art_colors, name = "ART Status",
                       labels = c("Pre-ART","Suppressed","Unsuppressed")) +
    theme_void() +
    theme(
      legend.direction  = "horizontal",
      legend.text       = element_text(size = 22),
      legend.title      = element_text(size = 24, face = "bold"),
      legend.key.size   = unit(1.5, "cm"),
      legend.spacing.x  = unit(0.8, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 10, alpha = 1), nrow = 1))
  
  s8_legend <- cowplot::get_legend(p_s8_legend_src)
  
  # Combine: 3 plots in one row + legend underneath
  if (length(s8_plot_list) == 3) {
    combined_s8 <- (s8_plot_list[[1]] | s8_plot_list[[2]] | s8_plot_list[[3]]) /
      wrap_elements(full = s8_legend) +
      plot_layout(heights = c(1, 0.08))
    
    ggsave(file.path(s8_dir, "S8_combined_pseudotime_modules.png"),
           plot = combined_s8, width = 27, height = 8, dpi = 300, bg = "white")
  }
  
  # Also save legend separately in case needed
  ggsave(file.path(s8_dir, "S8_common_legend.png"),
         plot = as_ggplot(s8_legend), width = 14, height = 2, dpi = 300, bg = "white")
}
cat("  S8 done.\n")

################################################################################
# S9: KIR+ INNATE-LIKE CD8 — increased fonts throughout
################################################################################
cat("\n=== S9: KIR+ innate-like CD8 ===\n")

if ("KIR+ innate-like CD8" %in% unique(TARA_cd8$CD8_Annotation)) {
  kir_cells <- subset(TARA_cd8, CD8_Annotation == "KIR+ innate-like CD8")
  kir_cells$Timepoint_Group <- factor(kir_cells$Timepoint_Group, levels = names(art_colors))
  
  # S9A: Proportion — INCREASED FONTS
  prop_df <- as.data.frame(table(kir_cells$Timepoint_Group))
  colnames(prop_df) <- c("ART_Status", "Count")
  prop_df$Pct <- round(100 * prop_df$Count / sum(prop_df$Count), 1)
  
  p_s9a <- ggplot(prop_df, aes(x = ART_Status, y = Count, fill = ART_Status)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(Pct, "%")), vjust = -0.5, size = 9, fontface = "bold") +
    scale_fill_manual(values = art_colors, guide = "none") +
    scale_x_discrete(labels = c("Pre-ART","Suppressed","Unsuppressed")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "", y = "Number of KIR+ cells", title = "KIR+ innate-like CD8 by ART status") +
    theme_cowplot(font_size = 24) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 24, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(s9_dir, "S9A_KIR_proportion.png"), plot = p_s9a,
         width = 9, height = 8, dpi = 300, bg = "white")
  
  # CSV: KIR proportion
  write.csv(prop_df, file.path(csv_dir, "S9A_KIR_proportion.csv"), row.names = FALSE)
  
  # S9B: Clonal expansion — INCREASED FONTS
  kir_clone <- kir_cells@meta.data %>%
    filter(!is.na(cloneSize)) %>%
    mutate(cloneSize = factor(cloneSize, levels = c("Single (0 < X <= 1)","Small (1 < X <= 5)",
                                                    "Medium (5 < X <= 20)","Large (20 < X <= 100)",
                                                    "Hyperexpanded (100 < X <= 500)")))
  
  p_s9b <- ggplot(kir_clone, aes(x = Timepoint_Group, fill = cloneSize)) +
    geom_bar(position = "fill", width = 0.6) +
    scale_fill_viridis_d(name = "Clone size", option = "mako", direction = -1) +
    scale_x_discrete(labels = c("Pre-ART","Suppressed","Unsuppressed")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "", y = "Proportion", title = "Clonal expansion in KIR+ CD8") +
    theme_cowplot(font_size = 24) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(s9_dir, "S9B_KIR_clonal_expansion.png"), plot = p_s9b,
         width = 11, height = 8, dpi = 300, bg = "white")
  
  # CSV: KIR clonal expansion counts
  kir_clone_summary <- kir_clone %>%
    group_by(Timepoint_Group, cloneSize) %>%
    summarise(n = n(), .groups = "drop")
  write.csv(kir_clone_summary, file.path(csv_dir, "S9B_KIR_clonal_expansion.csv"), row.names = FALSE)
  
  # S9C: Volcano — INCREASED FONTS
  cat("    S9C: KIR+ volcano...\n")
  kir_dge_path <- file.path(analysis_dir, "03_DGE", "DGE_MAST_KIR_Sup_vs_Unsup.csv")
  if (!file.exists(kir_dge_path)) {
    kir_post <- subset(kir_cells, Timepoint_Group %in% c("PostART_Suppressed","PostART_Unsuppressed"))
    if (min(table(kir_post$Timepoint_Group)) >= 10) {
      Idents(kir_post) <- "Timepoint_Group"
      DefaultAssay(kir_post) <- "RNA"
      kir_post <- JoinLayers(kir_post)
      kir_dge <- FindMarkers(kir_post, ident.1 = "PostART_Suppressed",
                             ident.2 = "PostART_Unsuppressed", test.use = "MAST",
                             min.pct = 0.1, logfc.threshold = 0.1)
      write.csv(kir_dge, file.path(s9_dir, "DGE_MAST_KIR_Sup_vs_Unsup.csv"))
      kir_dge_path <- file.path(s9_dir, "DGE_MAST_KIR_Sup_vs_Unsup.csv")
    }
  }
  
  if (file.exists(kir_dge_path)) {
    kir_dge <- read.csv(kir_dge_path, row.names = 1)
    kir_dge$gene <- rownames(kir_dge)
    kir_dge$neg_log10_padj <- -log10(kir_dge$p_val_adj + 1e-300)
    
    kir_dge <- kir_dge %>% filter(is.finite(avg_log2FC) & is.finite(neg_log10_padj))
    
    if (nrow(kir_dge) > 0) {
      kir_dge <- assign_highlight(kir_dge)
      
      kir_label <- genes_to_label_all[genes_to_label_all %in% kir_dge$gene]
      kir_dge$label <- ""
      kir_dge$label[kir_dge$gene %in% kir_label & kir_dge$p_val_adj < 0.05] <-
        kir_dge$gene[kir_dge$gene %in% kir_label & kir_dge$p_val_adj < 0.05]
      
      kx <- max(quantile(abs(kir_dge$avg_log2FC), 0.995) * 1.15, 0.5)
      ky <- max(quantile(kir_dge$neg_log10_padj, 0.995) * 1.10, 2)
      
      kir_dge$x_plot <- pmax(pmin(kir_dge$avg_log2FC, kx * 0.98), -kx * 0.98)
      kir_dge$y_plot <- pmin(kir_dge$neg_log10_padj, ky)
      
      p_s9c <- ggplot(kir_dge, aes(x = x_plot, y = y_plot, color = highlight)) +
        geom_point(data = kir_dge %>% filter(highlight == "Other"), size = 2, alpha = 0.3) +
        geom_point(data = kir_dge %>% filter(highlight != "Other"), size = 6, alpha = 0.85) +
        geom_label_repel(data = kir_dge %>% filter(label != ""),
                         aes(label = label, fill = highlight), color = "white",
                         size = 7, fontface = "bold.italic", max.overlaps = 40,
                         segment.size = 0.4, box.padding = 0.4, label.size = 0.2, show.legend = FALSE) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.4) +
        scale_color_manual(values = highlight_cols, name = "Gene category", drop = FALSE) +
        scale_fill_manual(values = highlight_cols, guide = "none") +
        scale_x_continuous(limits = c(-kx, kx)) +
        scale_y_continuous(limits = c(0, ky)) +
        labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~italic(p)),
             title = "KIR+ innate-like CD8 — suppressed vs unsuppressed") +
        theme_cowplot(font_size = 24) +
        theme(plot.title = element_text(size = 22, face = "bold"),
              axis.text = element_text(size = 18),
              axis.title = element_text(size = 20),
              legend.position = "right",
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20, face = "bold"),
              legend.key.size = unit(1.3, "cm"),
              plot.background = element_rect(fill = "white", color = NA)) +
        guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))
      
      ggsave(file.path(s9_dir, "S9C_KIR_volcano.png"), plot = p_s9c,
             width = 14, height = 10, dpi = 300, bg = "white")
      
      # CSV: copy KIR DGE
      file.copy(kir_dge_path, file.path(csv_dir, "S9C_KIR_DGE_Sup_vs_Unsup.csv"),
                overwrite = TRUE)
    } else {
      cat("    WARNING: No finite DGE values for KIR+ volcano.\n")
    }
  }
  
  # S9D: Shared clonotypes — INCREASED FONTS
  cat("    S9D: Shared clonotypes...\n")
  tcr_col <- NULL
  for (candidate in c("CTaa", "CTstrict", "CTgene", "CTnt", "clonotype_id")) {
    if (candidate %in% colnames(kir_cells@meta.data)) {
      tcr_col <- candidate
      break
    }
  }
  
  if (!is.null(tcr_col)) {
    cat("      Using TCR column:", tcr_col, "\n")
    top_clones <- kir_cells@meta.data %>%
      filter(!is.na(.data[[tcr_col]])) %>%
      group_by(.data[[tcr_col]]) %>% summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n)) %>% head(5) %>% pull(.data[[tcr_col]])
    
    if (length(top_clones) > 0) {
      shared <- TARA_cd8@meta.data %>%
        filter(.data[[tcr_col]] %in% top_clones) %>%
        group_by(.data[[tcr_col]], CD8_Annotation) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(Clone_short = paste0("Clone ", match(.data[[tcr_col]], top_clones)))
      
      p_s9d <- ggplot(shared, aes(x = Clone_short, y = n, fill = CD8_Annotation)) +
        geom_col(width = 0.6) +
        scale_fill_brewer(palette = "Set2", name = "CD8 cluster") +
        labs(x = "", y = "Number of cells", title = "Top KIR+ clonotypes across clusters") +
        theme_cowplot(font_size = 24) +
        theme(plot.title = element_text(size = 22, face = "bold"),
              axis.text = element_text(size = 18),
              axis.title = element_text(size = 20),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20, face = "bold"),
              legend.position = "right",
              plot.background = element_rect(fill = "white", color = NA)) +
        coord_flip()
      
      ggsave(file.path(s9_dir, "S9D_KIR_shared_clonotypes.png"), plot = p_s9d,
             width = 14, height = 8, dpi = 300, bg = "white")
      
      # CSV: shared clonotypes
      write.csv(shared, file.path(csv_dir, "S9D_KIR_shared_clonotypes.csv"), row.names = FALSE)
    } else {
      cat("      No expanded KIR+ clonotypes found.\n")
    }
  } else {
    cat("      No TCR clonotype column found. Available columns:\n")
    cat("      ", paste(grep("CT|clone|tcr", colnames(kir_cells@meta.data),
                             value = TRUE, ignore.case = TRUE), collapse = ", "), "\n")
  }
  cat("  S9 done.\n")
}

################################################################################
# SUMMARY
################################################################################
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUPPLEMENTARY COMPLETE\n")
cat("Output:", supp_base, "\n\n")
for (d in list(c("S4", s4_dir), c("S5", s5_dir), c("S6", s6_dir),
               c("S7", s7_dir), c("S8", s8_dir), c("S9", s9_dir),
               c("CSV", csv_dir))) {
  n <- length(list.files(d[2], "\\.png$|\\.csv$"))
  cat("  ", d[1], ":", n, "files\n")
}