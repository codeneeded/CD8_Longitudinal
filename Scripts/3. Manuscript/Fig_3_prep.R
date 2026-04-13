################################################################################
# FIGURE 3 — ANALYSIS & ANNOTATION (FINAL)
#
# Steps: Load → Clean → Score → Annotate → Verify → TARA cross-ref → Save
# Uses mnn_clusters_rna (12 clusters)
# Annotations: Naive/Bystander, Transitional Tem, Tpex, TEMRA/Effector,
#              Cycling Effector, Tex
################################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(qs2)

base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
ana_dir   <- file.path(base_dir, "Manuscript/Fig 3/analysis")
dir.create(ana_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading objects...\n")
cp003      <- qs_read(file.path(saved_dir, "cp003_with_cc_epitope.qs2"))
tara_cp003 <- qs_read(file.path(saved_dir, "tara_cp003_fig3_annotated.qs2"))

cluster_col <- "mnn_clusters_rna"


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: REMOVE CD8- SORTED CELLS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n################################################################\n")
cat("# STEP 1: REMOVE CD8- CELLS\n")
cat("################################################################\n\n")

cat("Total cells:", ncol(cp003), "\n")
print(table(cp003$CellType_Sort, useNA = "ifany"))

cp003_clean <- subset(cp003, CellType_Sort == "CD8+")
cat("After removing CD8-:", ncol(cp003_clean), "\n")
cat("Clusters:\n")
print(table(cp003_clean@meta.data[[cluster_col]]))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: CLUSTER SCORING
# ══════════════════════════════════════════════════════════════════════════════

cat("\n################################################################\n")
cat("# STEP 2: CLUSTER SCORING\n")
cat("################################################################\n\n")

gene_sets <- list(
  naive_stem   = c("TCF7", "LEF1", "SELL", "IL7R", "BACH2", "BCL2"),
  cytotoxicity = c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7"),
  exhaustion   = c("TOX", "PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT"),
  temra_markers = c("CX3CR1", "KLRG1", "B3GAT1"),
  cycling      = c("MKI67", "TOP2A"),
  tpex         = c("TCF7", "TOX"),
  transitional = c("GZMK", "EOMES", "CCL5")
)

avg_all <- AverageExpression(cp003_clean, features = unique(unlist(gene_sets)),
                             group.by = cluster_col, assays = "RNA")
avg_mat <- as.data.frame(avg_all$RNA)

scores <- data.frame(cluster = colnames(avg_mat))
for (gs in names(gene_sets)) {
  genes_use <- gene_sets[[gs]][gene_sets[[gs]] %in% rownames(avg_mat)]
  if (length(genes_use) > 0) scores[[gs]] <- round(colMeans(avg_mat[genes_use, , drop = FALSE]), 2)
}

cell_features <- cp003_clean@meta.data %>%
  group_by(!!sym(cluster_col)) %>%
  summarise(
    n_cells      = n(),
    pct_G2M_S    = round(sum(Phase %in% c("G2M", "S")) / n() * 100, 1),
    pct_HIV      = round(sum(HIV_Specific_TCR == "HIV-Specific TCR", na.rm = TRUE) / n() * 100, 1),
    pct_expanded = round(sum(clonalFrequency > 1, na.rm = TRUE) / sum(!is.na(clonalFrequency)) * 100, 1),
    .groups = "drop"
  )
cell_features$cluster <- paste0("g", cell_features[[cluster_col]])
scores <- merge(scores, cell_features[, c("cluster", "n_cells", "pct_G2M_S", "pct_HIV", "pct_expanded")],
                by = "cluster", all.x = TRUE)

print(scores)


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: ANNOTATE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n################################################################\n")
cat("# STEP 3: ANNOTATION\n")
cat("################################################################\n\n")

scores$annotation <- NA
for (i in seq_len(nrow(scores))) {
  r <- scores[i, ]
  if (r$pct_G2M_S > 30 & r$cytotoxicity > 5)
    scores$annotation[i] <- "Cycling Effector"
  else if (r$naive_stem > 5 & r$pct_HIV < 15 & r$pct_expanded < 15)
    scores$annotation[i] <- "Naive/Bystander"
  else if (r$tpex > 2.5 & r$exhaustion > 0.5 & r$cytotoxicity < 10)
    scores$annotation[i] <- "Tpex"
  else if (r$naive_stem > 3 & r$transitional > 5 & r$cytotoxicity < 15)
    scores$annotation[i] <- "Transitional Tem"
  else if (r$exhaustion > 2 & r$cytotoxicity < 10)
    scores$annotation[i] <- "Tex"
  else if (r$cytotoxicity > 8 & r$pct_G2M_S < 30)
    scores$annotation[i] <- "TEMRA/Effector"
  else if (r$pct_HIV > 30 & r$cytotoxicity > 3)
    scores$annotation[i] <- "TEMRA/Effector"
  else
    scores$annotation[i] <- "Unassigned"
}

cat("=== Annotations ===\n")
scores %>%
  select(cluster, n_cells, annotation, naive_stem, cytotoxicity, exhaustion, tpex,
         pct_G2M_S, pct_HIV, pct_expanded) %>%
  arrange(annotation, cluster) %>%
  as.data.frame() %>%
  print()

cat("\n=== Summary ===\n")
scores %>%
  group_by(annotation) %>%
  summarise(n_clusters = n(), total_cells = sum(n_cells),
            clusters = paste(cluster, collapse = ", "),
            mean_pct_HIV = round(mean(pct_HIV), 1), .groups = "drop") %>%
  as.data.frame() %>%
  print()


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: ASSIGN TO OBJECT + VERIFY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n################################################################\n")
cat("# STEP 4: VERIFICATION\n")
cat("################################################################\n\n")

cluster_to_anno <- setNames(scores$annotation, gsub("g", "", scores$cluster))
anno_vec <- cluster_to_anno[as.character(cp003_clean@meta.data[[cluster_col]])]
names(anno_vec) <- colnames(cp003_clean)
cp003_clean <- AddMetaData(cp003_clean, metadata = anno_vec, col.name = "Fig3_Annotation")

cat("=== Cells per annotation ===\n")
print(table(cp003_clean$Fig3_Annotation, useNA = "ifany"))

cat("\n=== Average expression per annotation ===\n")
verify_genes <- c("TCF7", "LEF1", "SELL", "IL7R", "GZMB", "GNLY", "PRF1", "NKG7",
                  "GZMK", "EOMES", "TOX", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
                  "CX3CR1", "KLRG1", "MKI67", "TOP2A", "CCL3", "CCL4", "CCL5",
                  "CD69", "HLA-DRA", "IFNG", "CD8A")
verify_genes <- verify_genes[verify_genes %in% rownames(cp003_clean)]
avg_by_anno <- AverageExpression(cp003_clean, features = verify_genes,
                                 group.by = "Fig3_Annotation", assays = "RNA")
print(round(as.data.frame(avg_by_anno$RNA), 2))

cat("\n=== Cell cycle by annotation ===\n")
cp003_clean@meta.data %>%
  group_by(Fig3_Annotation, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Fig3_Annotation) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  select(-n) %>%
  pivot_wider(names_from = Phase, values_from = pct, values_fill = 0) %>%
  as.data.frame() %>% print()

cat("\n=== HIV-specific % by annotation ===\n")
cp003_clean@meta.data %>%
  group_by(Fig3_Annotation) %>%
  summarise(n = n(), pct_hiv = round(sum(HIV_Specific_TCR == "HIV-Specific TCR", na.rm = TRUE) / n() * 100, 1),
            .groups = "drop") %>%
  as.data.frame() %>% print()

cat("\n=== HIV-specific cells by annotation AND timepoint ===\n")
cp003_clean@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
  group_by(Fig3_Annotation, Months) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Months, values_from = n, values_fill = 0, names_prefix = "m") %>%
  as.data.frame() %>% print()


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: TARA CROSS-REFERENCE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n################################################################\n")
cat("# STEP 5: TARA CROSS-REFERENCE\n")
cat("################################################################\n\n")

tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  group_by(orig.ident, Age) %>%
  summarise(n_cells = n(), n_clones = n_distinct(CTstrict), .groups = "drop") %>%
  as.data.frame() %>% print()

cat("\n")
tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  group_by(Manual_Annotation_refined) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  as.data.frame() %>% print()

hiv_val <- tara_cp003@meta.data %>% filter(HIV_validated_status == "HIV-Specific (validated)")
cat("\nTotal HIV-validated cells in TARA:", nrow(hiv_val), "\n")
cat("Unique clonotypes:", n_distinct(hiv_val$CTstrict), "\n")

stim_hiv_clones <- unique(cp003_clean@meta.data$CTstrict[cp003_clean$HIV_Specific_TCR == "HIV-Specific TCR"])
tara_clones <- unique(tara_cp003@meta.data$CTstrict[!is.na(tara_cp003$CTstrict)])
shared_clones <- intersect(stim_hiv_clones, tara_clones)
cat("Stim HIV-specific clonotypes:", length(stim_hiv_clones), "\n")
cat("Shared with TARA:", length(shared_clones), "\n")

if (length(shared_clones) > 0) {
  tara_tp <- tara_cp003@meta.data %>% filter(CTstrict %in% shared_clones) %>%
    group_by(CTstrict) %>%
    summarise(tara_tps = list(sort(unique(orig.ident))), n_tara = n_distinct(orig.ident),
              tara_cells = n(), .groups = "drop")
  stim_tp <- cp003_clean@meta.data %>%
    filter(CTstrict %in% shared_clones & HIV_Specific_TCR == "HIV-Specific TCR") %>%
    group_by(CTstrict) %>%
    summarise(n_stim = n_distinct(Months), stim_cells = n(), .groups = "drop")
  combined <- merge(tara_tp, stim_tp, by = "CTstrict", all = TRUE) %>%
    mutate(total_tp = n_tara + n_stim, total_cells = tara_cells + stim_cells) %>%
    arrange(desc(total_tp))
  cat("\nPer-clone persistence:\n")
  for (j in seq_len(nrow(combined))) {
    r <- combined[j, ]
    cat(sprintf("  Clone %d: %d TARA tp (%s) + %d stim tp -> %d total, %d cells\n",
                j, r$n_tara, paste(r$tara_tps[[1]], collapse = ","),
                r$n_stim, r$total_tp, r$total_cells))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════════════

qs_save(cp003_clean, file.path(saved_dir, "cp003_annotated_fig3.qs2"))
cat("\n✓ Saved:", file.path(saved_dir, "cp003_annotated_fig3.qs2"), "\n")
cat("  Cells:", ncol(cp003_clean), "| Clusters:", length(unique(cp003_clean@meta.data[[cluster_col]])), "\n")
cat("  Annotations:\n")
print(table(cp003_clean$Fig3_Annotation))

