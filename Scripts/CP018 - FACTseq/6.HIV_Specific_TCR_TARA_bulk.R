################################################################################
# CP018 FACT-Seq — final step: cross-reference HIV-specific clonotypes into the
# unstimulated TARA longitudinal data.
#
# Adapted from "CP003 - Longitudinal HIV Stim/7.HIV_Specific_TCR_TARA_bulk.R".
#
# This is what completes FACT-Seq: clonotypes validated by antigen-driven
# proliferation are looked up in the unstimulated dataset to establish
# (a) whether they exist at rest, (b) which CD8 subset they occupy, and
# (c) their resting functional state.
#
# WHY CP018 MATTERS: CP003 never suppressed within its sampled scRNA-seq
# timepoints (VL >500,000 at 12m and 24m). CP018 suppressed at 65 days
# (pre-ART VL 176,970; timepoints 1m/25m/42m). So this is the first look at
# HIV-specific clonotypes in a SUPPRESSED participant — directly testable
# against the Figure 4/5 finding that suppression installs type I IFN memory
# (IFIT1/IFIT3/ISG15/MX1) and preserves stemness (TCF7/SELL/BACH2).
#
# INPUT : QS_ANNOTATED (CP018 stim), TARA_ALL_annotated_final.qs2
# OUTPUT: HIV_ROOT/TARA_crossref/ tables + plots
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)
library(SeuratExtend)
library(patchwork)

out_dir <- file.path(HIV_ROOT, "TARA_crossref")
mkdirs(out_dir)

QS_TARA_ALL <- file.path(SAVED_DIR, "TARA_ALL_annotated_final.qs2")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load
# ══════════════════════════════════════════════════════════════════════════════
require_file(QS_ANNOTATED, "CP018 annotated stim object (run script 5 first)")
require_file(QS_TARA_ALL,  "TARA_ALL_annotated_final.qs2")

cp018 <- qs_read(QS_ANNOTATED)

if (!"HIV_Specific_TCR" %in% colnames(cp018@meta.data)) {
  stop("HIV_Specific_TCR not found — set BYSTANDER_CLUSTERS in script 5 and re-run it.")
}

TARA_ALL <- qs_read(QS_TARA_ALL)
TARA_ALL$PID <- sub("_.*$", "", TARA_ALL$orig.ident)

if (!"CP018" %in% unique(TARA_ALL$PID)) {
  stop("No CP018 cells in TARA_ALL. Available PIDs: ",
       paste(sort(unique(TARA_ALL$PID)), collapse = ", "))
}

tara_cp018 <- subset(TARA_ALL, subset = PID == "CP018")
rm(TARA_ALL); gc()

cat("CP018 stim cells:", ncol(cp018), "\n")
cat("TARA CP018 (unstimulated) cells:", ncol(tara_cp018), "\n")
cat("TARA CP018 timepoints:\n"); print(table(tara_cp018$orig.ident))

# The TARA annotation column changed names across manuscript versions
anno_col <- intersect(c("Annotation", "CD8_Annotation", "Manual_Annotation_refined"),
                      colnames(tara_cp018@meta.data))[1]
if (is.na(anno_col)) stop("No annotation column found in TARA object.")
cat("Using TARA annotation column:", anno_col, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Clone lists
# ══════════════════════════════════════════════════════════════════════════════
cp018_meta <- cp018@meta.data %>%
  transmute(
    Cell              = rownames(.),
    CTstrict          = CTstrict,
    HIV_Specific_TCR  = as.character(HIV_Specific_TCR),
    CP018_Sample      = Sample_Subfolder,
    CP018_Timepoint   = Timepoint,
    CP018_Replicate   = Replicate,
    CP018_Cluster     = as.character(mnn_clusters_rna),
    clonalFrequency   = suppressWarnings(as.numeric(as.character(clonalFrequency))),
    CP018_Phase       = if ("Phase" %in% colnames(.)) Phase else NA_character_
  )

hiv_clones <- cp018_meta %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict), CTstrict != "") %>%
  distinct(CTstrict) %>% pull(CTstrict)

all_stim_clones <- cp018_meta %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  distinct(CTstrict) %>% pull(CTstrict)

cat("\nCP018 HIV-specific clonotypes:", length(hiv_clones), "\n")
cat("CP018 all stim clonotypes:", length(all_stim_clones), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Match into unstimulated TARA
# ══════════════════════════════════════════════════════════════════════════════
if (!"CTstrict" %in% colnames(tara_cp018@meta.data)) {
  stop("TARA object has no CTstrict column — TCR data must be attached first.")
}

tara_meta <- tara_cp018@meta.data %>%
  transmute(
    Cell            = rownames(.),
    CTstrict        = CTstrict,
    TARA_Sample     = orig.ident,
    TARA_Cluster    = .data[[anno_col]],
    TARA_clonalFreq = suppressWarnings(as.numeric(as.character(clonalFrequency)))
  ) %>%
  filter(!is.na(CTstrict), CTstrict != "")

tara_meta$clone_status <- case_when(
  tara_meta$CTstrict %in% hiv_clones      ~ "HIV-Specific Match",
  tara_meta$CTstrict %in% all_stim_clones ~ "Stim Match (not HIV-specific)",
  TRUE                                    ~ "No Match"
)

matched <- tara_meta %>% filter(clone_status == "HIV-Specific Match")

cat("\n=== Cross-reference result ===\n")
cat("Unstimulated TARA CP018 cells carrying an HIV-specific clonotype:", nrow(matched), "\n")
cat("Unique HIV-specific clonotypes recovered:", n_distinct(matched$CTstrict),
    "of", length(hiv_clones), "\n")
if (nrow(matched) > 0) {
  cat("\nBy TARA timepoint:\n");  print(table(matched$TARA_Sample))
  cat("\nBy TARA CD8 subset:\n"); print(sort(table(matched$TARA_Cluster), decreasing = TRUE))
}

write.csv(tara_meta, file.path(out_dir, "CP018_TARA_all_cells_clone_status.csv"), row.names = FALSE)
write.csv(matched,   file.path(out_dir, "CP018_TARA_HIVspecific_matched_cells.csv"), row.names = FALSE)

# ── Clone-level master table ──
master <- matched %>%
  group_by(CTstrict) %>%
  summarise(
    n_cells_TARA      = n(),
    TARA_timepoints   = paste(sort(unique(TARA_Sample)), collapse = ","),
    n_TARA_timepoints = n_distinct(TARA_Sample),
    TARA_clusters     = paste(sort(unique(TARA_Cluster)), collapse = ","),
    max_TARA_freq     = suppressWarnings(max(TARA_clonalFreq, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(
    cp018_meta %>%
      filter(CTstrict %in% hiv_clones) %>%
      group_by(CTstrict) %>%
      summarise(
        n_cells_stim      = n(),
        stim_timepoints   = paste(sort(unique(CP018_Timepoint)), collapse = ","),
        n_stim_replicates = n_distinct(paste(CP018_Timepoint, CP018_Replicate)),
        max_stim_freq     = max(clonalFrequency, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "CTstrict"
  ) %>%
  arrange(desc(n_TARA_timepoints), desc(n_cells_TARA))

write.csv(master, file.path(out_dir, "CP018_Master_Clone_Match.csv"), row.names = FALSE)
if (nrow(master) > 0) print(head(as.data.frame(master), 20))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Resting functional state of matched cells
# ══════════════════════════════════════════════════════════════════════════════
# The CP003 Figure 3D analogue: % of matched cells expressing each key gene.
# Here it carries extra weight — CP018 is suppressed, so the IFN-memory module
# (IFIT1/IFIT3/ISG15/MX1) is the specific prediction to test.
if (nrow(matched) > 0) {

  tara_cp018$clone_match <- ifelse(
    colnames(tara_cp018) %in% matched$Cell, "HIV-Specific Match", "Other")

  gene_groups <- list(
    Cytotoxicity          = c("GZMB","GNLY","PRF1","NKG7","GZMA","GZMH"),
    `Cytokines/chemokines`= c("CCL3","CCL4","CCL5","IFNG","TNF"),
    `Effector TFs`        = c("TBX21","ZEB2","PRDM1","ID2"),
    `Effector markers`    = c("KLRG1","CX3CR1","FCGR3A","S1PR5"),
    Exhaustion            = c("PDCD1","LAG3","TIGIT","TOX","HAVCR2","CTLA4"),
    `Stemness/Naive`      = c("TCF7","SELL","BACH2","LEF1","IL7R","CCR7","BCL2"),
    `Type I IFN memory`   = c("IFIT1","IFIT3","ISG15","MX1")
  )

  expr <- LayerData(tara_cp018, assay = "RNA", layer = "data")
  match_idx <- which(tara_cp018$clone_match == "HIV-Specific Match")

  pct_rows <- lapply(names(gene_groups), function(grp) {
    genes <- intersect(gene_groups[[grp]], rownames(expr))
    if (length(genes) == 0) return(NULL)
    data.frame(
      Group = grp,
      Gene  = genes,
      pct_expressing = round(100 * Matrix::rowMeans(
        expr[genes, match_idx, drop = FALSE] > 0), 1),
      mean_expression = round(Matrix::rowMeans(
        expr[genes, match_idx, drop = FALSE]), 3),
      stringsAsFactors = FALSE
    )
  })
  pct_tab <- do.call(rbind, pct_rows)
  rownames(pct_tab) <- NULL

  write.csv(pct_tab, file.path(out_dir, "CP018_matched_cells_functional_genes.csv"),
            row.names = FALSE)
  cat("\n=== Resting functional state of matched HIV-specific cells (n=",
      length(match_idx), " cells) ===\n", sep = "")
  print(pct_tab)

  pct_tab$Gene <- factor(pct_tab$Gene, levels = pct_tab$Gene[order(pct_tab$Group, -pct_tab$pct_expressing)])

  p_func <- ggplot(pct_tab, aes(x = Gene, y = pct_expressing, fill = Group)) +
    geom_col(colour = "black", linewidth = 0.25) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    labs(title = "CP018 HIV-specific clonotypes — resting functional state",
         subtitle = paste0(length(match_idx),
                           " matched cells in unstimulated TARA data"),
         y = "% of cells expressing", x = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none",
          strip.text = element_text(size = 8))

  ggsave(file.path(out_dir, "CP018_matched_functional_genes.png"), p_func,
         width = 12, height = 5.5, dpi = 400, bg = "white")

  # Subset localisation (CP003 Figure 3F analogue)
  clust_tab <- as.data.frame(table(matched$TARA_Cluster)) %>%
    setNames(c("TARA_Cluster", "n_cells")) %>%
    mutate(pct = round(100 * n_cells / sum(n_cells), 1)) %>%
    arrange(desc(n_cells))
  write.csv(clust_tab, file.path(out_dir, "CP018_matched_cluster_distribution.csv"),
            row.names = FALSE)

  p_clust <- ggplot(clust_tab, aes(x = reorder(TARA_Cluster, n_cells), y = n_cells)) +
    geom_col(fill = "#D62728", colour = "black", linewidth = 0.25) +
    geom_text(aes(label = paste0(n_cells, " (", pct, "%)")), hjust = -0.1, size = 3) +
    coord_flip(clip = "off") +
    labs(title = "Where CP018 HIV-specific clonotypes sit at rest",
         x = NULL, y = "Cells") +
    theme_classic(base_size = 11) +
    theme(plot.margin = margin(5, 40, 5, 5))

  ggsave(file.path(out_dir, "CP018_matched_cluster_distribution.png"), p_clust,
         width = 7, height = 5, dpi = 400, bg = "white")

  # Module comparison: matched vs other, in the unstimulated data
  for (grp in c("Type I IFN memory", "Stemness/Naive", "Exhaustion", "Cytotoxicity")) {
    genes <- intersect(gene_groups[[grp]], rownames(tara_cp018))
    if (length(genes) < 2) next
    nm <- paste0("mod_", gsub("[^A-Za-z]", "", grp))
    tara_cp018 <- AddModuleScore(tara_cp018, features = list(genes), name = nm)
    ggsave(file.path(out_dir, paste0("CP018_TARA_", nm, "_matched_vs_other.png")),
           VlnPlot2(tara_cp018, features = paste0(nm, "1"),
                    group.by = "clone_match", show.mean = TRUE) +
             ggtitle(paste0(grp, " — unstimulated CP018")),
           width = 6, height = 5, dpi = 400, bg = "white")
  }

  qs_save(tara_cp018, file = file.path(SAVED_DIR, "tara_cp018_crossref.qs2"))
  message("Saved: ", file.path(SAVED_DIR, "tara_cp018_crossref.qs2"))

} else {
  cat("\nNo HIV-specific clonotypes recovered in the unstimulated data.\n",
      "For reference, CP003 recovered only 7 of 165 — a low hit rate is expected,\n",
      "since resting frequencies of antigen-specific clones are low.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — CP003 vs CP018 comparison scaffold
# ══════════════════════════════════════════════════════════════════════════════
# Shared clonotypes between two unrelated participants would be public TCRs
# (or a barcode-collision artefact) — worth checking either way.
if (file.exists(QS_CP003_ANNOT)) {
  cp003 <- qs_read(QS_CP003_ANNOT)
  cp003_hiv <- cp003@meta.data %>%
    filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict), CTstrict != "") %>%
    distinct(CTstrict) %>% pull(CTstrict)

  cross <- data.frame(
    participant = c("CP003", "CP018"),
    n_HIV_specific_clonotypes = c(length(cp003_hiv), length(hiv_clones)),
    ART_status_at_stim = c("Unsuppressed (VL >500,000)", "Suppressed (from day 65)"),
    stim_timepoints = c("2m, 101m", "2m, 90m")
  )
  write.csv(cross, file.path(out_dir, "CP003_vs_CP018_summary.csv"), row.names = FALSE)
  print(cross)

  shared_pub <- intersect(cp003_hiv, hiv_clones)
  cat("\nClonotypes shared between CP003 and CP018:", length(shared_pub), "\n")
  if (length(shared_pub) > 0) {
    write.csv(data.frame(CTstrict = shared_pub),
              file.path(out_dir, "CP003_CP018_shared_clonotypes.csv"), row.names = FALSE)
    cat("Inspect these — public TCR specificities, or a barcode artefact.\n")
  }
  rm(cp003); gc()
} else {
  message("CP003 annotated object not found; skipping cross-participant comparison.")
}

cat("\nDone. Outputs in:\n  ", out_dir, "\n")
