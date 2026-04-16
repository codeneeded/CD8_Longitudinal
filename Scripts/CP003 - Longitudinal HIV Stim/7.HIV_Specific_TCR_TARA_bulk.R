################################################################################
# CP003 — COMBINED ANALYSIS PIPELINE (v5)
#
# Merges: clone matching + cell cycle + Trex + cluster annotation
# Outputs to: ~/Documents/CD8_Longitudinal/Analysis/CP003/
#
# INPUT:
#   - seu_CP003_HIVSpecificTCR_annotated.qs2  (CP003 stim data)
#   - TARA_ALL_annotated_final.qs2            (full TARA, new annotations)
#
# OUTPUT (saved to Analysis/CP003/):
#   - cp003_fig3_annotated.qs2    → CP003 stim: cell cycle, Trex, cluster annotations
#   - tara_cp003_fig3.qs2         → TARA CP003 subset: HIV_validated_status
#   - Master_Clone_Match.csv      → clone-level match table
#   - HIV_Specific_Clones.csv     → HIV-specific subset
#
# CHANGES FROM v4:
#   - Uses TARA_ALL_annotated_final.qs2 (new annotations)
#   - Annotation column instead of Manual_Annotation_refined / CD8_Annotation
#   - Combined from two separate scripts
#   - Tscm gating removed (separate analysis, not Figure 3)
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)
library(Trex)
library(ggpubr)
library(rstatix)
library(SeuratExtend)
library(patchwork)


# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════

base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
out_dir   <- file.path(base_dir, "Analysis", "CP003")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: LOAD OBJECTS
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("  STEP 1: LOADING OBJECTS\n")
cat("================================================================\n\n")

cp003    <- qs_read(file.path(saved_dir, "seu_CP003_HIVSpecificTCR_annotated.qs2"))
TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_annotated_final.qs2"))

# CHANGED: Subset TARA to CP003 participant (all timepoints)
TARA_ALL$PID <- sub("_.*$", "", TARA_ALL$orig.ident)
tara_cp003 <- subset(TARA_ALL, subset = PID == "CP003")

cat("CP003 stim cells:", ncol(cp003), "\n")
cat("TARA CP003 cells:", ncol(tara_cp003), "\n")
cat("TARA CP003 timepoints:\n")
print(table(tara_cp003$orig.ident))
cat("\n")

rm(TARA_ALL); gc()


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: CELL CYCLE SCORING — CP003
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("  STEP 2: CELL CYCLE SCORING\n")
cat("================================================================\n\n")

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cp003 <- JoinLayers(cp003)
cp003 <- CellCycleScoring(
  cp003,
  s.features   = s.genes,
  g2m.features = g2m.genes,
  set.ident    = FALSE
)
cat("Cell cycle phases:\n")
print(table(cp003$Phase))
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: TREX EPITOPE ANNOTATION
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("  STEP 3: TREX EPITOPE ANNOTATION\n")
cat("================================================================\n\n")

pull_trex_meta <- function(obj, ed_label) {
  obj@meta.data %>%
    transmute(
      Cell         = rownames(.),
      CTstrict     = CTstrict,
      !!paste0("Epitope.target_", ed_label) := TRB_Epitope.target,
      !!paste0("Epitope.sequence_", ed_label) := TRB_Epitope.sequence,
      !!paste0("Epitope.species_", ed_label) := TRB_Epitope.species,
      !!paste0("Epitope.database_", ed_label) := TRB_Database
    )
}

# ── CP003 Trex ───────────────────────────────────────────────────────────────
cat("Running Trex on CP003 (ED0, ED1, ED2)...\n")
cp003_trex_0 <- annotateDB(cp003, chains = "TRB")
cp003_trex_1 <- annotateDB(cp003, chains = "TRB", edit.distance = 1)
cp003_trex_2 <- annotateDB(cp003, chains = "TRB", edit.distance = 2)

cp003_trex_df <- pull_trex_meta(cp003_trex_0, "ED0") %>%
  left_join(pull_trex_meta(cp003_trex_1, "ED1"), by = c("Cell", "CTstrict")) %>%
  left_join(pull_trex_meta(cp003_trex_2, "ED2"), by = c("Cell", "CTstrict"))

cp003 <- AddMetaData(cp003, cp003_trex_df %>% tibble::column_to_rownames("Cell"))
cat("CP003 Trex complete.\n\n")

# ── TARA CP003 Trex ──────────────────────────────────────────────────────────
cat("Running Trex on TARA CP003 subset (ED0, ED1, ED2)...\n")
tara_trex_0 <- annotateDB(tara_cp003, chains = "TRB")
tara_trex_1 <- annotateDB(tara_cp003, chains = "TRB", edit.distance = 1)
tara_trex_2 <- annotateDB(tara_cp003, chains = "TRB", edit.distance = 2)

tara_trex_df <- pull_trex_meta(tara_trex_0, "ED0") %>%
  left_join(pull_trex_meta(tara_trex_1, "ED1"), by = c("Cell", "CTstrict")) %>%
  left_join(pull_trex_meta(tara_trex_2, "ED2"), by = c("Cell", "CTstrict"))

tara_cp003 <- AddMetaData(tara_cp003, tara_trex_df %>% tibble::column_to_rownames("Cell"))
cat("TARA CP003 Trex complete.\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: HIV-SPECIFIC CLONE MATCHING
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("  STEP 4: CLONE MATCHING\n")
cat("================================================================\n\n")

# ── CP003 clone lists ────────────────────────────────────────────────────────
cp003_meta <- cp003@meta.data %>%
  transmute(
    Cell             = rownames(.),
    CTstrict         = CTstrict,
    HIV_Specific_TCR = HIV_Specific_TCR,
    Sample           = Sample,
    CP003_Cluster    = mnn_snn_res.0.6,
    clonalFrequency  = clonalFrequency,
    CP003_Phase      = Phase,
    CP003_S_Score    = S.Score,
    CP003_G2M_Score  = G2M.Score,
    CP003_Epitope_ED0 = Epitope.target_ED0,
    CP003_Epitope_ED1 = Epitope.target_ED1,
    CP003_Epitope_ED2 = Epitope.target_ED2,
    CP003_Species_ED0 = Epitope.species_ED0,
    CP003_Species_ED1 = Epitope.species_ED1,
    CP003_Species_ED2 = Epitope.species_ED2
  )

cp003_hiv_clones <- cp003_meta %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict), CTstrict != "") %>%
  distinct(CTstrict) %>% pull(CTstrict)

cp003_all_clones <- cp003_meta %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  distinct(CTstrict) %>% pull(CTstrict)

cat("Unique HIV-specific CP003 clones:", length(cp003_hiv_clones), "\n")
cat("Total CP003 clones:", length(cp003_all_clones), "\n")

# ── CP003 clone-level summary ────────────────────────────────────────────────
cp003_clone_origin <- cp003_meta %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict) %>%
  summarise(
    CP003_HIV_Specific   = paste(sort(unique(HIV_Specific_TCR)), collapse = ";"),
    CP003_Samples        = paste(sort(unique(Sample)), collapse = ";"),
    CP003_Clusters       = paste(sort(unique(CP003_Cluster)), collapse = ";"),
    CP003_n_cells        = n(),
    CP003_max_clonalFreq = max(suppressWarnings(as.numeric(as.character(clonalFrequency))), na.rm = TRUE),
    CP003_Phase_breakdown = paste(sort(unique(CP003_Phase)), collapse = ";"),
    CP003_pct_S          = round(mean(CP003_Phase == "S", na.rm = TRUE) * 100, 1),
    CP003_pct_G2M        = round(mean(CP003_Phase == "G2M", na.rm = TRUE) * 100, 1),
    CP003_pct_G1         = round(mean(CP003_Phase == "G1", na.rm = TRUE) * 100, 1),
    CP003_Epitope_ED0    = paste(sort(unique(na.omit(CP003_Epitope_ED0))), collapse = ";"),
    CP003_Epitope_ED2    = paste(sort(unique(na.omit(CP003_Epitope_ED2))), collapse = ";"),
    CP003_Species_ED0    = paste(sort(unique(na.omit(CP003_Species_ED0))), collapse = ";"),
    CP003_Species_ED2    = paste(sort(unique(na.omit(CP003_Species_ED2))), collapse = ";"),
    .groups = "drop"
  )

# ── Label matching clones in TARA CP003 ──────────────────────────────────────
# CHANGED: uses Annotation column instead of CD8_Annotation
tara_meta <- tara_cp003@meta.data %>%
  transmute(
    Cell             = rownames(.),
    CTstrict         = CTstrict,
    Sample           = orig.ident,
    Tara_Cluster     = Annotation,
    clonalFrequency  = clonalFrequency,
    Tara_Epitope_ED0 = Epitope.target_ED0,
    Tara_Epitope_ED2 = Epitope.target_ED2,
    Tara_Species_ED0 = Epitope.species_ED0,
    Tara_Species_ED2 = Epitope.species_ED2
  ) %>%
  mutate(
    CP003_clone_status = case_when(
      !is.na(CTstrict) & CTstrict %in% cp003_hiv_clones ~ "HIV-Specific Match",
      !is.na(CTstrict) & CTstrict %in% cp003_all_clones ~ "Non-HIV Match",
      TRUE ~ "Other"
    ),
    CP003_clone_status = factor(CP003_clone_status,
                                levels = c("Other", "Non-HIV Match", "HIV-Specific Match"))
  )

tara_cp003$CP003_clone_status <- tara_meta$CP003_clone_status

# ── Build master match table ─────────────────────────────────────────────────
tara_clone_summary <- tara_meta %>%
  filter(CP003_clone_status != "Other", !is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict, CP003_clone_status) %>%
  summarise(
    Tara_Samples        = paste(sort(unique(trimws(as.character(Sample)))), collapse = ";"),
    Tara_n_samples      = n_distinct(Sample),
    Tara_Clusters       = paste(sort(unique(Tara_Cluster)), collapse = ";"),
    Tara_n_cells        = n(),
    Tara_max_clonalFreq = max(suppressWarnings(as.numeric(as.character(clonalFrequency))), na.rm = TRUE),
    Tara_Epitope_ED0    = paste(sort(unique(na.omit(Tara_Epitope_ED0))), collapse = ";"),
    Tara_Epitope_ED2    = paste(sort(unique(na.omit(Tara_Epitope_ED2))), collapse = ";"),
    Tara_Species_ED0    = paste(sort(unique(na.omit(Tara_Species_ED0))), collapse = ";"),
    Tara_Species_ED2    = paste(sort(unique(na.omit(Tara_Species_ED2))), collapse = ";"),
    .groups = "drop"
  )

master_match_table <- tara_clone_summary %>%
  left_join(cp003_clone_origin, by = "CTstrict") %>%
  select(-CP003_HIV_Specific) %>%
  mutate(
    All_Samples = apply(
      select(., CP003_Samples, Tara_Samples), 1,
      function(x) paste(sort(unique(unlist(strsplit(paste(na.omit(x), collapse = ";"), ";")))), collapse = ";")
    )
  ) %>%
  select(-CP003_Samples, -Tara_Samples) %>%
  relocate(All_Samples, .after = CTstrict)

# ── NOTE: Species refinement removed ─────────────────────────────────────────
# The previous pipeline used Trex species predictions to promote "Non-HIV Match"
# clones to "HIV-Specific Match" if no non-HIV species was predicted. This is
# speculative — a clone that was expanded in the stim assay but NOT classified
# as HIV-Specific TCR (e.g., in Naive/Bystander cluster) should not be assumed
# HIV-reactive just because no database match was found.
#
# HIV_validated_status is now set ONLY from clones that were functionally
# classified as HIV-Specific TCR in the peptide stimulation assay.

master_match_table <- master_match_table %>%
  arrange(desc(CP003_clone_status), desc(Tara_n_cells))

# ── Create HIV_validated_status in TARA object ───────────────────────────────
# Use ONLY stim-confirmed HIV-specific clones (cp003_hiv_clones from Step 4)
tara_cp003$HIV_validated_status <- ifelse(
  !is.na(tara_cp003$CTstrict) & tara_cp003$CTstrict %in% cp003_hiv_clones,
  "HIV-Specific (validated)",
  "Other"
)

cat("\nHIV-validated cells in TARA CP003:\n")
print(table(tara_cp003$HIV_validated_status))
cat("\n")

# ── Save match tables ────────────────────────────────────────────────────────
write.csv(master_match_table,
          file.path(out_dir, "Master_Clone_Match.csv"), row.names = FALSE)
write.csv(master_match_table %>% filter(CP003_clone_status == "HIV-Specific Match"),
          file.path(out_dir, "HIV_Specific_Clones.csv"), row.names = FALSE)
write.csv(cp003_clone_origin,
          file.path(out_dir, "CP003_Clone_Origin_Summary.csv"), row.names = FALSE)

cat("Match tables saved to:", out_dir, "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: REMOVE CD8- CELLS FROM CP003
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("  STEP 5: REMOVE CD8- CELLS\n")
cat("================================================================\n\n")

cat("Total CP003 cells:", ncol(cp003), "\n")
print(table(cp003$CellType_Sort, useNA = "ifany"))

cp003_clean <- subset(cp003, CellType_Sort == "CD8+")
cat("After removing CD8-:", ncol(cp003_clean), "\n\n")

cluster_col <- "mnn_clusters_rna"


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: CLUSTER SCORING + ANNOTATION
# ══════════════════════════════════════════════════════════════════════════════

cat("================================================================\n")
cat("  STEP 6: CLUSTER SCORING + ANNOTATION\n")
cat("================================================================\n\n")

gene_sets <- list(
  naive_stem    = c("TCF7", "LEF1", "SELL", "IL7R", "BACH2", "BCL2"),
  cytotoxicity  = c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7"),
  exhaustion    = c("TOX", "PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT"),
  temra_markers = c("CX3CR1", "KLRG1", "B3GAT1"),
  cycling       = c("MKI67", "TOP2A"),
  tpex          = c("TCF7", "TOX"),
  transitional  = c("GZMK", "EOMES", "CCL5")
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

cat("Cluster scores:\n")
print(scores)

# ── Annotation rules ─────────────────────────────────────────────────────────
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

cat("\nAnnotation summary:\n")
scores %>%
  select(cluster, n_cells, annotation, naive_stem, cytotoxicity, exhaustion,
         pct_G2M_S, pct_HIV, pct_expanded) %>%
  arrange(annotation) %>%
  as.data.frame() %>% print()

# ── Assign to object ─────────────────────────────────────────────────────────
cluster_to_anno <- setNames(scores$annotation, gsub("g", "", scores$cluster))
anno_vec <- cluster_to_anno[as.character(cp003_clean@meta.data[[cluster_col]])]
names(anno_vec) <- colnames(cp003_clean)
cp003_clean <- AddMetaData(cp003_clean, metadata = anno_vec, col.name = "Fig3_Annotation")

cat("\nCells per annotation:\n")
print(table(cp003_clean$Fig3_Annotation, useNA = "ifany"))


# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: VERIFICATION
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  STEP 7: VERIFICATION\n")
cat("================================================================\n\n")

cat("Cell cycle by annotation:\n")
cp003_clean@meta.data %>%
  group_by(Fig3_Annotation, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Fig3_Annotation) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  select(-n) %>%
  pivot_wider(names_from = Phase, values_from = pct, values_fill = 0) %>%
  as.data.frame() %>% print()

cat("\nHIV-specific % by annotation:\n")
cp003_clean@meta.data %>%
  group_by(Fig3_Annotation) %>%
  summarise(n = n(),
            pct_hiv = round(sum(HIV_Specific_TCR == "HIV-Specific TCR", na.rm = TRUE) / n() * 100, 1),
            .groups = "drop") %>%
  as.data.frame() %>% print()

cat("\nHIV-specific cells by annotation AND timepoint:\n")
cp003_clean@meta.data %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
  group_by(Fig3_Annotation, Months) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Months, values_from = n, values_fill = 0, names_prefix = "m") %>%
  as.data.frame() %>% print()


# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: TARA CROSS-REFERENCE VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  STEP 8: TARA CROSS-REFERENCE\n")
cat("================================================================\n\n")

cat("HIV-validated cells per timepoint:\n")
tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  group_by(orig.ident, Age) %>%
  summarise(n_cells = n(), n_clones = n_distinct(CTstrict), .groups = "drop") %>%
  as.data.frame() %>% print()

# CHANGED: uses Annotation column
cat("\nHIV-validated cells per cluster:\n")
tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)") %>%
  group_by(Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  as.data.frame() %>% print()

hiv_val <- tara_cp003@meta.data %>%
  filter(HIV_validated_status == "HIV-Specific (validated)")
cat("\nTotal HIV-validated cells in TARA:", nrow(hiv_val), "\n")
cat("Unique validated clonotypes:", n_distinct(hiv_val$CTstrict), "\n")

# ── Clone persistence tracking ───────────────────────────────────────────────
stim_hiv_clones <- unique(cp003_clean@meta.data$CTstrict[
  cp003_clean$HIV_Specific_TCR == "HIV-Specific TCR"])
tara_clones <- unique(tara_cp003@meta.data$CTstrict[!is.na(tara_cp003$CTstrict)])
shared_clones <- intersect(stim_hiv_clones, tara_clones)

cat("Stim HIV-specific clonotypes:", length(stim_hiv_clones), "\n")
cat("Shared with TARA:", length(shared_clones), "\n")

if (length(shared_clones) > 0) {
  tara_tp <- tara_cp003@meta.data %>%
    filter(CTstrict %in% shared_clones) %>%
    group_by(CTstrict) %>%
    summarise(tara_tps = list(sort(unique(orig.ident))),
              n_tara = n_distinct(orig.ident),
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
# STEP 9: SAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n================================================================\n")
cat("  STEP 9: SAVING\n")
cat("================================================================\n\n")

qs_save(cp003_clean, file.path(saved_dir, "cp003_fig3_annotated.qs2"))
cat("Saved:", file.path(saved_dir, "cp003_fig3_annotated.qs2"), "\n")
cat("  Cells:", ncol(cp003_clean), "\n")
cat("  Annotations:", paste(names(table(cp003_clean$Fig3_Annotation)), collapse = ", "), "\n\n")

qs_save(tara_cp003, file.path(saved_dir, "tara_cp003_fig3.qs2"))
cat("Saved:", file.path(saved_dir, "tara_cp003_fig3.qs2"), "\n")
cat("  Cells:", ncol(tara_cp003), "\n")
cat("  HIV-validated:", sum(tara_cp003$HIV_validated_status == "HIV-Specific (validated)"), "\n\n")

# Save annotation scores for reference (CSV → Analysis dir)
write.csv(scores, file.path(out_dir, "Cluster_Annotation_Scores.csv"), row.names = FALSE)

cat("================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("  qs2 objects: ", saved_dir, "\n")
cat("  CSVs/tables: ", out_dir, "\n")
cat("================================================================\n")