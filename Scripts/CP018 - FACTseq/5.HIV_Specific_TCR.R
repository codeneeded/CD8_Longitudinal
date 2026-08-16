############################################################
# CP018 FACT-Seq — define HIV-specific clonotypes
#
# Adapted from "CP003 - Longitudinal HIV Stim/6.HIV_Specific_TCR.R".
#
# FACT-Seq definition (as published for CP003): a cell is HIV-specific if it is
#   (i)  clonally expanded (clonalFrequency > 1), AND
#   (ii) NOT in the Naive/Bystander compartment
# The reasoning: clonal expansion coupled with active proliferation under HIV
# peptide stimulation reflects antigen-driven selection.
#
# IMPORTANT — MANUAL STEP: for CP003 the bystander clusters were MNN clusters
# 0 and 1, identified from marker expression. Cluster numbering is not portable
# between datasets - which is why the HIV call is now expansion-based only.
# annotation plots (script 3) before the HIV-specific call is meaningful.
# Bystander = high TCF7/LEF1/SELL/IL7R, no clonal expansion, low cytotoxicity.
#
# Input : QS_WITH_TCR
# Output: QS_ANNOTATED + HIV_ROOT tables/plots
############################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(qs2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratExtend)
library(scCustomize)

mkdirs(HIV_ROOT)

require_file(QS_WITH_TCR, "integrated object with TCR")
seu <- qs_read(QS_WITH_TCR)

# Fig_Annotation is created in script 3. Plots use LABELS, not cluster numbers.
if (!"Fig_Annotation" %in% colnames(seu@meta.data)) {
  stop("Fig_Annotation missing - re-run 3.Integration_Anno_Plots.R first ",
       "(it applies CP018_CLUSTER_ANNOTATION).", call. = FALSE)
}
seu$Fig_Annotation <- factor(as.character(seu$Fig_Annotation),
                             levels = intersect(ANNO_ORDER,
                                                unique(as.character(seu$Fig_Annotation))))

DefaultAssay(seu) <- "RNA"
Idents(seu) <- "mnn_clusters_rna"

# ══════════════════════════════════════════════════════════════════════════════
# SET THIS after reviewing script 3's annotation plots.
# NULL = not yet set; the script will refuse to write the HIV-specific call.
# ── HIV-specific definition ──────────────────────────────────────────────────
# Set in config.R: HIV_MIN_CLONE_SIZE, HIV_PER_TIMEPOINT, HIV_EXCLUDE_MAIT.
# NO CLUSTER NUMBERS ARE USED. The previous cluster-number approach broke when
# the object's mnn_clusters_rna came from a different resolution than the
# numbers were picked at; MAIT exclusion now goes by annotation LABEL.
# ══════════════════════════════════════════════════════════════════════════════

# -----------------------------
# 1) Bystander diagnostics (run FIRST, to choose the clusters above)
# -----------------------------
# Per-cluster mean module scores + expansion rate. The bystander cluster is the
# one with high stemness, near-zero expansion, and low cytotoxicity.
stem_genes  <- intersect(c("TCF7","LEF1","SELL","CCR7","IL7R","BACH2","BCL2"), rownames(seu))
cyto_genes  <- intersect(c("GZMB","GNLY","PRF1","NKG7","GZMA","GZMH"), rownames(seu))
exh_genes   <- intersect(c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4"), rownames(seu))
cycle_genes <- intersect(c("MKI67","TOP2A","PCNA"), rownames(seu))

seu <- AddModuleScore(seu, features = list(stem_genes),  name = "Stemness_score")
seu <- AddModuleScore(seu, features = list(cyto_genes),  name = "Cytotoxicity_score")
seu <- AddModuleScore(seu, features = list(exh_genes),   name = "Exhaustion_score")
seu <- AddModuleScore(seu, features = list(cycle_genes), name = "Cycling_score")

freq_num <- suppressWarnings(as.numeric(as.character(seu$clonalFrequency)))
seu$clonalFrequency_num <- freq_num
seu$is_expanded <- !is.na(freq_num) & freq_num > 1

cluster_profile <- seu@meta.data %>%
  group_by(mnn_clusters_rna) %>%
  summarise(
    n_cells        = n(),
    pct_expanded   = round(100 * mean(is_expanded), 1),
    pct_with_TCR   = round(100 * mean(!is.na(CTstrict)), 1),
    mean_stemness  = round(mean(Stemness_score1), 3),
    mean_cytotox   = round(mean(Cytotoxicity_score1), 3),
    mean_exhaust   = round(mean(Exhaustion_score1), 3),
    mean_cycling   = round(mean(Cycling_score1), 3),
    pct_S_G2M      = if ("Phase" %in% colnames(seu@meta.data))
                       round(100 * mean(Phase %in% c("S", "G2M")), 1) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_stemness))

write.csv(cluster_profile, file.path(HIV_ROOT, "CP018_cluster_profile_for_bystander_call.csv"),
          row.names = FALSE)

cat("\n=== Per-cluster profile (reference; the HIV call no longer uses cluster numbers) ===\n")
cat("Bystander = HIGH mean_stemness, LOW pct_expanded, LOW mean_cytotox\n\n")
print(as.data.frame(cluster_profile))

for (m in c("Stemness_score1", "Cytotoxicity_score1", "Exhaustion_score1", "Cycling_score1")) {
  ggsave(file.path(HIV_ROOT, paste0("Module_", sub("1$", "", m), "_by_cluster.png")),
         VlnPlot2(seu, features = m, group.by = "Fig_Annotation", show.mean = TRUE) +
           labs(x = NULL) + rotate_x(35),
         width = 11, height = 6, dpi = 400, bg = "white")
}

ggsave(file.path(HIV_ROOT, "CP018_UMAP_expanded_cells.png"),
       DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "is_expanded"),
       width = 8, height = 6, dpi = 400, bg = "white")

# -----------------------------
# 2) HIV-specific call
# -----------------------------
  # Expansion is counted WITHIN a timepoint: a clonotype with one cell at 2m and
  # one at 90m has a pooled frequency of 2 but never expanded at either visit.
  md_tmp <- seu@meta.data
  md_tmp$.ct <- as.character(md_tmp$CTstrict)
  md_tmp$.tp <- as.character(md_tmp$Timepoint)

  if (isTRUE(HIV_PER_TIMEPOINT)) {
    key <- paste(md_tmp$.ct, md_tmp$.tp, sep = "||")
    size_within <- table(key[!is.na(md_tmp$.ct)])
    clone_size  <- as.integer(size_within[key])
  } else {
    size_pool  <- table(md_tmp$.ct[!is.na(md_tmp$.ct)])
    clone_size <- as.integer(size_pool[md_tmp$.ct])
  }
  clone_size[is.na(md_tmp$.ct)] <- NA_integer_
  seu$clone_size_used <- clone_size

  # MAIT exclusion by LABEL, not cluster number.
  is_mait <- rep(FALSE, ncol(seu))
  if (isTRUE(HIV_EXCLUDE_MAIT)) {
    if ("Fig_Annotation" %in% colnames(seu@meta.data)) {
      is_mait <- grepl("MAIT", as.character(seu$Fig_Annotation), ignore.case = TRUE)
    } else if (exists("CP018_CLUSTER_ANNOTATION") && !is.null(CP018_CLUSTER_ANNOTATION)) {
      lab <- CP018_CLUSTER_ANNOTATION[as.character(seu$mnn_clusters_rna)]
      is_mait <- !is.na(lab) & grepl("MAIT", lab, ignore.case = TRUE)
    } else {
      warning("No annotation available - MAIT cells NOT excluded.")
    }
  }
  seu$is_MAIT_cell <- is_mait

  seu$HIV_Specific_TCR <- factor(
    ifelse(!is.na(clone_size) & clone_size > HIV_MIN_CLONE_SIZE & !is_mait,
           "HIV-Specific TCR", "Other"),
    levels = c("Other", "HIV-Specific TCR")
  )

  cat("\n=== HIV-specific call ===\n")
  cat("  definition: clone size >", HIV_MIN_CLONE_SIZE,
      if (isTRUE(HIV_PER_TIMEPOINT)) "within a timepoint" else "pooled across timepoints",
      if (isTRUE(HIV_EXCLUDE_MAIT)) ", MAIT excluded" else "", "\n")
  cat("  MAIT cells excluded:", sum(is_mait), "\n")
  print(table(seu$HIV_Specific_TCR, seu$Timepoint))

  cat("\nHIV-specific cells:\n"); print(table(seu$HIV_Specific_TCR))
  cat("\nBy timepoint:\n"); print(table(seu$Timepoint, seu$HIV_Specific_TCR))
  cat("\nBy library:\n");   print(table(seu$Sample_Subfolder, seu$HIV_Specific_TCR))

  # ── Per-cell table ──
  hiv_cells <- data.frame(
    Cell             = colnames(seu),
    Sample           = seu$Sample_Subfolder,
    Timepoint        = seu$Timepoint,
    Replicate        = seu$Replicate,
    Cluster          = as.character(seu$mnn_clusters_rna),
    CTstrict         = seu$CTstrict,
    clonalFrequency  = freq_num,
    HIV_Specific_TCR = as.character(seu$HIV_Specific_TCR),
    Phase            = if ("Phase" %in% colnames(seu@meta.data)) seu$Phase else NA,
    stringsAsFactors = FALSE
  ) %>%
    filter(HIV_Specific_TCR == "HIV-Specific TCR",
           !is.na(CTstrict), CTstrict != "",
           !is.na(clonalFrequency), clonalFrequency > 1)

  write.csv(hiv_cells, file.path(HIV_ROOT, "CP018_HIV_Specific_TCR_cells.csv"),
            row.names = FALSE)

  # ── Clone-level summary ──
  hiv_clone_summary <- hiv_cells %>%
    group_by(CTstrict) %>%
    summarise(
      n_cells             = n(),
      n_timepoints        = n_distinct(Timepoint),
      timepoints          = paste(sort(unique(Timepoint)), collapse = ","),
      n_replicates        = n_distinct(paste(Timepoint, Replicate)),
      max_clonalFrequency = max(clonalFrequency, na.rm = TRUE),
      clusters            = paste(sort(unique(Cluster)), collapse = ","),
      .groups = "drop"
    ) %>%
    arrange(desc(max_clonalFrequency), desc(n_cells))

  write.csv(hiv_clone_summary, file.path(HIV_ROOT, "CP018_HIV_Specific_TCR_clone_summary.csv"),
            row.names = FALSE)

  cat("\nUnique HIV-specific clonotypes:", nrow(hiv_clone_summary), "\n")
  cat("Detected at BOTH timepoints:", sum(hiv_clone_summary$n_timepoints > 1), "\n")
  cat("Found in BOTH replicates of a timepoint (higher confidence):",
      sum(hiv_clone_summary$n_replicates > 1), "\n")
  print(head(as.data.frame(hiv_clone_summary), 20))

  # A clonotype seen in both technical replicates is not a library artefact —
  # a confidence tier CP003 could not assign at its 101m timepoint.
  write.csv(
    hiv_clone_summary %>% filter(n_replicates > 1),
    file.path(HIV_ROOT, "CP018_HIV_Specific_clonotypes_replicate_confirmed.csv"),
    row.names = FALSE
  )

  # ── Plots ──
  ggsave(file.path(HIV_ROOT, "CP018_UMAP_HIV_Specific.png"),
         DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "HIV_Specific_TCR",
                  cols = c("Other" = "grey85", "HIV-Specific TCR" = "#D62728")),
         width = 8, height = 6, dpi = 400, bg = "white")

  ggsave(file.path(HIV_ROOT, "CP018_UMAP_HIV_Specific_by_Timepoint.png"),
         DimPlot2(seu, reduction = "umap.mnn.rna", group.by = "HIV_Specific_TCR",
                  split.by = "Timepoint",
                  cols = c("Other" = "grey85", "HIV-Specific TCR" = "#D62728")),
         width = 13, height = 6, dpi = 400, bg = "white")

  # Module scores: HIV-specific vs bystander (the CP003 Supplementary S3B panel)
  for (m in c("Stemness_score1", "Cytotoxicity_score1", "Exhaustion_score1")) {
    ggsave(file.path(HIV_ROOT, paste0("HIVspecific_vs_Other_", sub("1$", "", m), ".png")),
           VlnPlot2(seu, features = m, group.by = "HIV_Specific_TCR",
                    split.by = "Timepoint", show.mean = TRUE),
           width = 7, height = 5, dpi = 400, bg = "white")
  }

  # -----------------------------
  # CSV exports — HIV-specific call, fully auditable
  # -----------------------------
  csv_dir <- file.path(HIV_ROOT, "tables"); mkdirs(csv_dir)

  md <- seu@meta.data
  md$cell_barcode <- rownames(md)

  # 1. the call itself, per cell, with every input to the decision
  save_csv(md %>%
    transmute(cell_barcode, Sample_Subfolder, Timepoint, Replicate,
              cluster = as.character(mnn_clusters_rna),
              Fig_Annotation = if ("Fig_Annotation" %in% names(md)) as.character(Fig_Annotation) else NA,
              CTstrict, CTaa, clonalFrequency,
              is_expanded       = !is.na(clonalFrequency) & clonalFrequency > 1,
              annotation        = if ("Fig_Annotation" %in% names(md)) as.character(Fig_Annotation) else NA_character_,
              clone_size_used   = clone_size_used,
              is_MAIT           = is_MAIT_cell,
              HIV_Specific_TCR),
    file.path(csv_dir, "CP018_HIV_specific_call_per_cell.csv"))

  # 2. one row per HIV-specific clonotype
  hiv_clono <- md %>%
    filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict)) %>%
    group_by(CTstrict) %>%
    summarise(CTaa = dplyr::first(CTaa), CTgene = dplyr::first(CTgene),
              n_cells = n(),
              n_2m = sum(Timepoint == "2m"), n_90m = sum(Timepoint == "90m"),
              replicate_confirmed = dplyr::n_distinct(Replicate) > 1,
              clusters = paste(sort(unique(as.character(mnn_clusters_rna))), collapse = ";"),
              .groups = "drop") %>%
    arrange(desc(n_cells))
  save_csv(hiv_clono, file.path(csv_dir, "CP018_HIV_specific_clonotypes.csv"))

  # 3. summary counts — the numbers that go in the text
  save_csv(data.frame(
    metric = c("total_cells","cells_with_TCR","expanded_cells",
               "MAIT_cells_excluded",
               "HIV_specific_cells","HIV_specific_clonotypes",
               "HIV_specific_replicate_confirmed",
               "HIV_specific_cells_2m","HIV_specific_cells_90m"),
    value = c(nrow(md),
              sum(!is.na(md$CTstrict)),
              sum(md$clone_size_used > HIV_MIN_CLONE_SIZE, na.rm = TRUE),
              sum(md$is_MAIT_cell, na.rm = TRUE),
              sum(md$HIV_Specific_TCR == "HIV-Specific TCR", na.rm = TRUE),
              nrow(hiv_clono),
              sum(hiv_clono$replicate_confirmed),
              sum(md$HIV_Specific_TCR == "HIV-Specific TCR" & md$Timepoint == "2m", na.rm = TRUE),
              sum(md$HIV_Specific_TCR == "HIV-Specific TCR" & md$Timepoint == "90m", na.rm = TRUE))),
    file.path(csv_dir, "CP018_HIV_specific_summary.csv"))

  # 4. module scores per cell and per group
  score_cols <- grep("_score1$", colnames(md), value = TRUE)
  if (length(score_cols)) {
    save_csv(md %>% select(cell_barcode, Timepoint, Replicate,
                           cluster = mnn_clusters_rna, HIV_Specific_TCR,
                           all_of(score_cols)),
             file.path(csv_dir, "CP018_module_scores_per_cell.csv"))
    save_csv(md %>% group_by(HIV_Specific_TCR, Timepoint) %>%
               summarise(n = n(), across(all_of(score_cols),
                                         list(mean = ~round(mean(.x), 4),
                                              sd   = ~round(sd(.x), 4))),
                         .groups = "drop"),
             file.path(csv_dir, "CP018_module_scores_summary.csv"))
  }
  export_metadata(seu, file.path(csv_dir, "CP018_cell_metadata_annotated.csv"))

  qs_save(seu, file = QS_ANNOTATED)
  message("Saved: ", QS_ANNOTATED)
