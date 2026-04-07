################################################################################
# FIGURE 3 — MAIN PANELS + SUPPLEMENTARY S3 (FINAL)
#
# MAIN FIGURE 3:
#   A  — Experimental schematic (Illustrator)
#   B1 — UMAP: annotated clusters (no legend, labels on plot)
#   B2 — UMAP: HIV-specific TCR overlay
#   C  — Alluvial: clonotype tracking across timepoints
#   D  — Functional gene expression: HIV-validated cells in TARA
#   E  — Persistence histogram
#   F  — TARA cluster composition of HIV-validated clones
#
# SUPPLEMENTARY S3:
#   S3A — Cell cycle by annotation
#   S3B — Module scores: HIV-specific vs Other (stimulation)
#   S3C — Dot plot of key markers per annotation
#   S3D — HIV-specific % per annotation
#   S3E — Timepoint composition per annotation
#   S3F — Module scores within HIV-specific cells, 2m vs 101m
#
# REQUIRES: cp003_annotated_fig3.qs2, tara_cp003_fig3_annotated.qs2
################################################################################

library(Seurat)
library(SeuratExtend)
library(scRepertoire)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)
library(RColorBrewer)
library(qs2)

base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
panel_dir <- file.path(base_dir, "Manuscript/Fig 3/panels")
s3_dir    <- file.path(base_dir, "Manuscript/Supplementary 3")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(s3_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading...\n")
cp003_clean <- qs_read(file.path(saved_dir, "cp003_annotated_fig3.qs2"))
tara_cp003  <- qs_read(file.path(saved_dir, "tara_cp003_fig3_annotated.qs2"))

red_use <- "umap.mnn.rna"

# ═══════════════════════════════════════════════════════════════
# PALETTES
# ═══════════════════════════════════════════════════════════════

anno_order <- c("Naive/Bystander", "Transitional Tem", "Tpex",
                "TEMRA/Effector", "Cycling Effector", "Tex")

anno_cols <- c(
  "Naive/Bystander"  = "#4E79A7",
  "Transitional Tem" = "#59A14F",
  "Tpex"             = "#EDC948",
  "TEMRA/Effector"   = "#E15759",
  "Cycling Effector"  = "#FF9DA7",
  "Tex"              = "#9C755F"
)

cp003_clean$Fig3_Annotation <- factor(cp003_clean$Fig3_Annotation, levels = anno_order)

hiv_cols   <- c("Other" = "grey85", "HIV-Specific TCR" = "#E63946")
phase_cols <- c("G1" = "#AEC6CF", "S" = "#FFD700", "G2M" = "#E63946")
tp_cols    <- c("2 months" = "#E63946", "101 months" = "#4E79A7")


################################################################################
#                          MAIN FIGURE 3 PANELS
################################################################################

# ═══════════════════════════════════════════════════════════════
# PANEL B1: UMAP — Annotated clusters (no legend, large labels)
# ═══════════════════════════════════════════════════════════════
cat("Panel B1: Annotated UMAP...\n")

panelB1 <- DimPlot2(
  cp003_clean, reduction = red_use, group.by = "Fig3_Annotation",
  cols = 'default', label = TRUE, box = TRUE, label.size = 13,
  repel = TRUE, pt.size = 0.6, raster = FALSE,
  theme = list(NoAxes(), NoLegend(), theme_umap_arrows())
) + ggtitle("CP003 stimulation") +
  theme(plot.title = element_text(size = 32, face = "bold"))

ggsave(file.path(panel_dir, "Fig3B1_UMAP_Annotated.png"),
       panelB1, width = 9, height = 9, dpi = 400, bg = "white")


# ═══════════════════════════════════════════════════════════════
# PANEL B2: UMAP — HIV-Specific TCR overlay
# ═══════════════════════════════════════════════════════════════
cat("Panel B2: HIV-Specific UMAP...\n")

panelB2 <- DimPlot(
  cp003_clean, reduction = red_use, group.by = "HIV_Specific_TCR",
  cols = hiv_cols, pt.size = 0.6, raster = FALSE
) + NoAxes() + NoLegend() +
  ggtitle("HIV-specific clones (red)") +
  theme(plot.title = element_text(size = 32, face = "bold", color = "#E63946"))

ggsave(file.path(panel_dir, "Fig3B2_UMAP_HIV_Specific.png"),
       panelB2, width = 9, height = 9, dpi = 400, bg = "white")
# ═══════════════════════════════════════════════════════════════
# PANEL C: Alluvial
# ═══════════════════════════════════════════════════════════════
cat("Panel C: Alluvial...\n")

stim_alluvial <- cp003_clean@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
  mutate(timepoint_label = ifelse(Months == 2, "2 months\n(HIV-stim)", "101 months\n(HIV-stim)")) %>%
  group_by(CTstrict, timepoint_label) %>%
  summarise(n_cells = n(), .groups = "drop")

tara_alluvial <- tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  mutate(timepoint_label = paste0(Age, " months")) %>%
  group_by(CTstrict, timepoint_label) %>%
  summarise(n_cells = n(), .groups = "drop")

alluvial_combined <- bind_rows(stim_alluvial, tara_alluvial)
tp_order <- c("2 months\n(HIV-stim)", "1 months", "12 months",
              "24 months", "101 months\n(HIV-stim)")
alluvial_combined$timepoint_label <- factor(alluvial_combined$timepoint_label, levels = tp_order)

multi_tp <- alluvial_combined %>%
  group_by(CTstrict) %>% filter(n_distinct(timepoint_label) >= 2) %>% ungroup()

if (nrow(multi_tp) > 0) {
  clone_order <- multi_tp %>%
    group_by(CTstrict) %>% summarise(total = sum(n_cells), .groups = "drop") %>%
    arrange(desc(total))
  n_clones <- nrow(clone_order)
  clone_pal <- colorRampPalette(brewer.pal(min(n_clones, 12), "Set3"))(n_clones)
  names(clone_pal) <- clone_order$CTstrict
  
  panelC <- ggplot(multi_tp,
                   aes(x = timepoint_label, y = n_cells,
                       stratum = CTstrict, alluvium = CTstrict, fill = CTstrict)) +
    geom_stratum(width = 0.3, color = "grey30", linewidth = 0.3) +
    geom_flow(alpha = 0.5, width = 0.3) +
    scale_fill_manual(values = clone_pal, guide = "none") +
    scale_x_discrete(expand = c(0.08, 0.08), drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic(base_size = 28) +
    theme(
      axis.text.x   = element_text(size = 24, hjust = 0.5, lineheight = 0.9, face = "bold"),
      axis.text.y   = element_text(size = 22),
      axis.title    = element_text(size = 26, face = "bold"),
      plot.title    = element_text(size = 32, face = "bold"),
      plot.subtitle = element_text(size = 22, color = "grey40")
    ) +
    labs(x = NULL, y = "Number of cells",
         title = "HIV-specific clonotype persistence",
         subtitle = paste0(n_clones, " clonotypes at \u22652 timepoints"))
  
  ggsave(file.path(panel_dir, "Fig3C_Alluvial.png"),
         panelC, width = 12, height = 9, dpi = 400, bg = "white")
}


# ═══════════════════════════════════════════════════════════════
# PANEL D: Functional gene expression — HIV-validated in TARA
# ═══════════════════════════════════════════════════════════════
cat("Panel D: Functional profile...\n")

DefaultAssay(tara_cp003) <- "RNA"
hiv_val_barcodes <- rownames(tara_cp003@meta.data[tara_cp003$HIV_validated_status == "HIV-Specific (validated)", ])

func_genes <- c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7",
                "IFNG", "CCL3", "CCL4", "CCL5",
                "TBX21", "EOMES", "PRDM1",
                "CX3CR1", "KLRG1", "FAS",
                "PDCD1", "LAG3", "HAVCR2", "TOX")

gene_stats <- do.call(rbind, lapply(func_genes, function(g) {
  expr <- GetAssayData(tara_cp003, assay = "RNA", layer = "data")[g, hiv_val_barcodes]
  data.frame(gene = g,
             pct_expressing = round(mean(expr > 0) * 100, 1),
             mean_expr = round(mean(expr), 3))
}))

gene_stats$category <- case_when(
  gene_stats$gene %in% c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7") ~ "Cytotoxicity",
  gene_stats$gene %in% c("IFNG", "CCL3", "CCL4", "CCL5") ~ "Cytokines/Chemokines",
  gene_stats$gene %in% c("TBX21", "EOMES", "PRDM1") ~ "Effector TFs",
  gene_stats$gene %in% c("CX3CR1", "KLRG1", "FAS") ~ "Effector Markers",
  gene_stats$gene %in% c("PDCD1", "LAG3", "HAVCR2", "TOX") ~ "Exhaustion"
)

gene_stats$category <- factor(gene_stats$category,
                              levels = c("Cytotoxicity", "Cytokines/Chemokines",
                                         "Effector TFs", "Effector Markers", "Exhaustion"))
gene_stats$gene <- factor(gene_stats$gene, levels = rev(func_genes))

cat_cols <- c(
  "Cytotoxicity" = "#E15759", "Cytokines/Chemokines" = "#FF9DA7",
  "Effector TFs" = "#F28E2B", "Effector Markers" = "#EDC948",
  "Exhaustion" = "#9C755F"
)

panelD <- ggplot(gene_stats, aes(x = gene, y = pct_expressing, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0(pct_expressing, "%")),
            hjust = -0.08, size = 7, fontface = "bold") +
  scale_fill_manual(values = cat_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_flip() +
  theme_classic(base_size = 28) +
  theme(
    axis.text.x   = element_text(size = 22),
    axis.text.y   = element_text(size = 22, face = "italic"),
    axis.title    = element_text(size = 26, face = "bold"),
    legend.position = "bottom",
    legend.text   = element_text(size = 20),
    plot.title    = element_text(size = 32, face = "bold"),
    plot.subtitle = element_text(size = 22, color = "grey40")
  ) +
  labs(x = NULL, y = "% cells expressing",
       title = "Functional profile of HIV-specific clones",
       subtitle = "127 validated cells in TARA unstimulated data") +
  guides(fill = guide_legend(nrow = 2))

ggsave(file.path(panel_dir, "Fig3D_Functional_Profile.png"),
       panelD, width = 12, height = 10, dpi = 400, bg = "white")


# ═══════════════════════════════════════════════════════════════
# PANEL E: Persistence histogram
# ═══════════════════════════════════════════════════════════════
cat("Panel E: Persistence...\n")

stim_tp_per_clone <- cp003_clean@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
  mutate(tp = paste0("stim_", Months, "m")) %>%
  group_by(CTstrict) %>%
  summarise(stim_tps = list(unique(tp)), .groups = "drop")

tara_tp_per_clone <- tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  mutate(tp = paste0("tara_", Age, "m")) %>%
  group_by(CTstrict) %>%
  summarise(tara_tps = list(unique(tp)), .groups = "drop")

all_hiv_clones <- data.frame(CTstrict = unique(
  cp003_clean@meta.data$CTstrict[cp003_clean$HIV_Specific_TCR == "HIV-Specific TCR"]))

persist_full <- all_hiv_clones %>%
  left_join(stim_tp_per_clone, by = "CTstrict") %>%
  left_join(tara_tp_per_clone, by = "CTstrict") %>%
  rowwise() %>%
  mutate(n_timepoints = length(unique(c(unlist(stim_tps), unlist(tara_tps))))) %>%
  ungroup()

persist_summary <- persist_full %>% count(n_timepoints) %>%
  mutate(n_timepoints = factor(n_timepoints))

panelE <- ggplot(persist_summary, aes(x = n_timepoints, y = n)) +
  geom_col(fill = "#E63946", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 10, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 28) +
  theme(
    axis.text     = element_text(size = 24),
    axis.title    = element_text(size = 26, face = "bold"),
    plot.title    = element_text(size = 32, face = "bold"),
    plot.subtitle = element_text(size = 22, color = "grey40")
  ) +
  labs(x = "Timepoints detected", y = "N clonotypes",
       title = "Clonotype persistence",
       subtitle = paste0(nrow(persist_full), " HIV-specific clonotypes"))

ggsave(file.path(panel_dir, "Fig3E_Persistence.png"),
       panelE, width = 8, height = 8, dpi = 400, bg = "white")


# ═══════════════════════════════════════════════════════════════
# PANEL F: TARA cluster composition
# ═══════════════════════════════════════════════════════════════
cat("Panel F: TARA clusters...\n")

tara_hiv_clusters <- tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  group_by(Manual_Annotation_refined) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = round(n / sum(n) * 100, 1),
         Manual_Annotation_refined = reorder(Manual_Annotation_refined, n))

panelF <- ggplot(tara_hiv_clusters,
                 aes(x = Manual_Annotation_refined, y = n, fill = Manual_Annotation_refined)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", pct, "%)")),
            hjust = -0.08, size = 7, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  scale_fill_manual(values = c("TEMRA/CTL" = "#E15759", "Tpex CD8" = "#EDC948",
                               "DNAM1+ CD4 T cell" = "#4E79A7", "DN T cell_3" = "#76B7B2",
                               "Gamma Delta 1 T cells" = "#B07AA1", "Plasmablast" = "#FF9DA7")) +
  coord_flip() +
  theme_classic(base_size = 28) +
  theme(
    axis.text.x   = element_text(size = 22),
    axis.text.y   = element_text(size = 22),
    axis.title    = element_text(size = 26, face = "bold"),
    plot.title    = element_text(size = 32, face = "bold"),
    plot.subtitle = element_text(size = 22, color = "grey40")
  ) +
  labs(x = NULL, y = "N cells",
       title = "Cluster distribution in TARA",
       subtitle = "127 cells from 7 validated clonotypes")

ggsave(file.path(panel_dir, "Fig3F_TARA_Clusters.png"),
       panelF, width = 10, height = 6, dpi = 400, bg = "white")


################################################################################
#                        SUPPLEMENTARY FIGURE S3
################################################################################
cat("\nGenerating Supplementary S3...\n")

# ═══════════════════════════════════════════════════════════════
# S3A: Cell cycle by annotation
# ═══════════════════════════════════════════════════════════════
cat("  S3A: Cell cycle...\n")

cc_data <- cp003_clean@meta.data %>%
  filter(!is.na(Phase), !is.na(Fig3_Annotation)) %>%
  group_by(Fig3_Annotation, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Fig3_Annotation) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

p_s3a <- ggplot(cc_data, aes(x = Fig3_Annotation, y = proportion, fill = Phase)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = phase_cols) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x  = element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 18),
    axis.title   = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title   = element_text(size = 24, face = "bold")
  ) +
  labs(x = NULL, y = "Proportion", fill = "Cell cycle",
       title = "Cell cycle by cluster annotation")

ggsave(file.path(s3_dir, "S3A_CellCycle.png"),
       p_s3a, width = 9, height = 6, dpi = 300, bg = "white")


# ═══════════════════════════════════════════════════════════════
# S3B: Module scores — HIV-specific vs Other (stimulation)
# ═══════════════════════════════════════════════════════════════
cat("  S3B: Module scores HIV vs Other (stim)...\n")

modules <- list(
  Cytotoxicity = c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7", "FASLG"),
  Exhaustion   = c("TOX", "PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "CD244", "ENTPD1"),
  Stemness     = c("TCF7", "LEF1", "SELL", "CCR7", "IL7R", "BACH2", "BCL2", "S1PR1"),
  `IFN Memory` = c("IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "STAT1")
)

for (mod_name in names(modules)) {
  col_name <- paste0(gsub(" ", "_", mod_name), "_s3b_")
  genes_present <- modules[[mod_name]][modules[[mod_name]] %in% rownames(cp003_clean)]
  if (length(genes_present) >= 2) {
    cp003_clean <- AddModuleScore(cp003_clean, features = list(genes_present),
                                  name = col_name, seed = 42)
  }
}

s3b_cols <- grep("_s3b_1$", colnames(cp003_clean@meta.data), value = TRUE)

violin_stim <- cp003_clean@meta.data %>%
  select(HIV_Specific_TCR, all_of(s3b_cols)) %>%
  pivot_longer(cols = all_of(s3b_cols), names_to = "Module", values_to = "Score") %>%
  mutate(
    Module = gsub("_s3b_1$", "", Module) %>% gsub("_", " ", .),
    Module = factor(Module, levels = c("Cytotoxicity", "Exhaustion", "Stemness", "IFN Memory")),
    HIV_Specific_TCR = factor(HIV_Specific_TCR, levels = c("Other", "HIV-Specific TCR"))
  )

p_s3b <- ggplot(violin_stim, aes(x = Module, y = Score, fill = HIV_Specific_TCR)) +
  geom_violin(scale = "width", width = 0.8, alpha = 0.8,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Other" = "grey70", "HIV-Specific TCR" = "#E63946"), name = NULL) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x  = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 18),
    axis.title   = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 24, face = "bold")
  ) +
  labs(x = NULL, y = "Module score",
       title = "HIV-specific vs bystander (stimulation data)")

ggsave(file.path(s3_dir, "S3B_ModuleScores_HIVvsOther.png"),
       p_s3b, width = 10, height = 7, dpi = 300, bg = "white")


# ═══════════════════════════════════════════════════════════════
# S3C: Dot plot of key markers per annotation
# ═══════════════════════════════════════════════════════════════
cat("  S3C: Dot plot...\n")

dot_genes <- c("TCF7", "LEF1", "SELL", "IL7R", "BACH2",
               "GZMK", "EOMES",
               "GZMB", "GNLY", "PRF1", "NKG7", "CX3CR1", "KLRG1",
               "TOX", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
               "MKI67", "TOP2A",
               "CCL3", "CCL4", "CCL5",
               "IFNG", "CD69", "HLA-DRA")
dot_genes <- dot_genes[dot_genes %in% rownames(cp003_clean)]
Idents(cp003_clean) <- "Fig3_Annotation"

p_s3c <- DotPlot(cp003_clean, features = dot_genes, group.by = "Fig3_Annotation",
                 cols = c("grey90", "#E15759"), dot.scale = 8) +
  coord_flip() +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x  = element_text(size = 16, angle = 40, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title   = element_text(size = 22, face = "bold")
  ) +
  labs(title = "Key marker expression by annotation")

ggsave(file.path(s3_dir, "S3C_DotPlot.png"),
       p_s3c, width = 11, height = 11, dpi = 300, bg = "white")


# ═══════════════════════════════════════════════════════════════
# S3D: HIV-specific % per annotation
# ═══════════════════════════════════════════════════════════════
cat("  S3D: HIV % per annotation...\n")

hiv_pct <- cp003_clean@meta.data %>%
  group_by(Fig3_Annotation) %>%
  summarise(pct_hiv = round(sum(HIV_Specific_TCR == "HIV-Specific TCR", na.rm = TRUE) / n() * 100, 1),
            .groups = "drop")

p_s3d <- ggplot(hiv_pct, aes(x = Fig3_Annotation, y = pct_hiv, fill = Fig3_Annotation)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct_hiv, "%")), vjust = -0.5, size = 7, fontface = "bold") +
  scale_fill_manual(values = anno_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x  = element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 18),
    axis.title   = element_text(size = 20, face = "bold"),
    plot.title   = element_text(size = 24, face = "bold")
  ) +
  labs(x = NULL, y = "% HIV-Specific TCR",
       title = "HIV-specific cell enrichment by cluster")

ggsave(file.path(s3_dir, "S3D_HIV_Pct.png"),
       p_s3d, width = 8, height = 6, dpi = 300, bg = "white")


# ═══════════════════════════════════════════════════════════════
# S3E: Timepoint composition per annotation
# ═══════════════════════════════════════════════════════════════
cat("  S3E: Timepoint composition...\n")

tp_comp <- cp003_clean@meta.data %>%
  mutate(Timepoint = ifelse(Months == 2, "2 months", "101 months")) %>%
  group_by(Fig3_Annotation, Timepoint) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Fig3_Annotation) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

p_s3e <- ggplot(tp_comp, aes(x = Fig3_Annotation, y = pct, fill = Timepoint)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tp_cols) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x  = element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 18),
    axis.title   = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 24, face = "bold")
  ) +
  labs(x = NULL, y = "Proportion", fill = NULL,
       title = "Timepoint composition by cluster")

ggsave(file.path(s3_dir, "S3E_Timepoint.png"),
       p_s3e, width = 9, height = 6, dpi = 300, bg = "white")


# ═══════════════════════════════════════════════════════════════
# S3F: Module scores within HIV-specific cells, 2m vs 101m
# ═══════════════════════════════════════════════════════════════
cat("  S3F: Module scores 2m vs 101m...\n")

for (mod_name in names(modules)) {
  col_name <- paste0(gsub(" ", "_", mod_name), "_s3f_")
  genes_present <- modules[[mod_name]][modules[[mod_name]] %in% rownames(cp003_clean)]
  if (length(genes_present) >= 2) {
    cp003_clean <- AddModuleScore(cp003_clean, features = list(genes_present),
                                  name = col_name, seed = 42)
  }
}

s3f_cols <- grep("_s3f_1$", colnames(cp003_clean@meta.data), value = TRUE)

hiv_tp_violin <- cp003_clean@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
  mutate(Timepoint = ifelse(Months == 2, "2 months", "101 months"),
         Timepoint = factor(Timepoint, levels = c("2 months", "101 months"))) %>%
  select(Timepoint, all_of(s3f_cols)) %>%
  pivot_longer(cols = all_of(s3f_cols), names_to = "Module", values_to = "Score") %>%
  mutate(Module = gsub("_s3f_1$", "", Module) %>% gsub("_", " ", .),
         Module = factor(Module, levels = c("Cytotoxicity", "Exhaustion", "Stemness", "IFN Memory")))

p_s3f <- ggplot(hiv_tp_violin, aes(x = Module, y = Score, fill = Timepoint)) +
  geom_violin(scale = "width", width = 0.8, alpha = 0.8,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = tp_cols, name = NULL) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x  = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 18),
    axis.title   = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 24, face = "bold"),
    plot.subtitle = element_text(size = 18, color = "grey40")
  ) +
  labs(x = NULL, y = "Module score",
       title = "HIV-specific cells: 2 months vs 101 months",
       subtitle = paste0("n = ", sum(cp003_clean$HIV_Specific_TCR == "HIV-Specific TCR" & cp003_clean$Months == 2),
                         " (2m), n = ", sum(cp003_clean$HIV_Specific_TCR == "HIV-Specific TCR" & cp003_clean$Months == 101),
                         " (101m)"))

ggsave(file.path(s3_dir, "S3F_ModuleScores_2m_vs_101m.png"),
       p_s3f, width = 10, height = 7, dpi = 300, bg = "white")


################################################################################
# DONE
################################################################################
cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(" MAIN FIGURE 3:", panel_dir, "\n")
cat("  B1: UMAP annotated (no legend, labels on plot)\n")
cat("  B2: UMAP HIV-specific overlay\n")
cat("  C:  Alluvial tracking\n")
cat("  D:  Functional profile (% expressing, TARA unstimulated)\n")
cat("  E:  Persistence histogram\n")
cat("  F:  TARA cluster composition\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(" SUPPLEMENTARY S3:", s3_dir, "\n")
cat("  S3A: Cell cycle by annotation\n")
cat("  S3B: Module scores HIV vs Other (stimulation)\n")
cat("  S3C: Dot plot markers\n")
cat("  S3D: HIV % per annotation\n")
cat("  S3E: Timepoint composition\n")
cat("  S3F: Module scores 2m vs 101m (HIV-specific)\n")
cat("══════════════════════════════════════════════════════════════\n")