############################################################
# Figure 3 — Analysis Script
# HIV-specific clonal expansion validation & longitudinal tracking
#
# INPUTS (already saved from prior scripts):
#   - cp003_with_cc_epitope.qs2       (CP003 stimulation Seurat, with TCR + HIV labels + Trex)
#   - tara_cdnk_with_cp003_matches.qs2 (TARA cohort Seurat, with CP003 clone matching + Trex)
#   - Master_Clone_Match_CP003_vs_Tara.csv (clone-level match table)
#
# PURPOSE:
#   1. Consolidate + QC the HIV-specific clone set from CP003 stimulation
#   2. Track validated clones across TARA longitudinal timepoints
#   3. Reclassify Fig 2C "unknown" clones using functional validation
#   4. Characterize phenotype of matched clones IN THE TARA (unstimulated) data
#   5. Export analysis tables for Fig 3 panels
#
# OUTPUT directory: Fig3/analysis/
############################################################

library(Seurat)
library(SeuratExtend)
library(scRepertoire)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(qs2)
library(patchwork)

# =============================================================
# 0) PATHS — UPDATE THESE TO YOUR LOCAL PATHS
# =============================================================
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
fig3_dir  <- file.path(base_dir, "Manuscript/Fig3")
ana_dir   <- file.path(fig3_dir, "analysis")
panel_dir <- file.path(fig3_dir, "panels")

dir.create(ana_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# Input objects
qs_cp003 <- file.path(saved_dir, "cp003_with_cc_epitope.qs2")
qs_tara  <- file.path(saved_dir, "tara_cdnk_with_cp003_matches.qs2")
master_csv <- file.path(base_dir, "CP003/Shared_clone_Tables_Bulk/Master_Clone_Match_CP003_vs_Tara.csv")

# =============================================================
# 1) LOAD OBJECTS
# =============================================================
cat("Loading CP003 object...\n")
cp003 <- qs_read(qs_cp003)

cat("Loading TARA object...\n")
tara <- qs_read(qs_tara)

cat("Loading master match table...\n")
master_match <- read.csv(master_csv, stringsAsFactors = FALSE)

# =============================================================
# 2) SUMMARIZE CP003 STIMULATION EXPERIMENT
#    What came out of the HIV stimulation + sort?
# =============================================================

# --- 2A: Cell counts per library/timepoint ---
cp003_counts <- cp003@meta.data %>%
  group_by(Sample, Timepoint, CellType_Sort) %>%
  summarise(n_cells = n(), .groups = "drop")

cat("\n=== CP003 cell counts by library ===\n")
print(cp003_counts)

# --- 2B: Clone landscape from stimulation ---
# Cells with valid TCR
cp003_tcr <- cp003@meta.data %>%
  filter(!is.na(CTstrict), CTstrict != "")

cat("\nCells with TCR:", nrow(cp003_tcr), "of", ncol(cp003), "total\n")

# Clonally expanded (freq > 1) = our HIV-specific candidates
cp003_tcr$freq_num <- suppressWarnings(as.numeric(as.character(cp003_tcr$clonalFrequency)))

cp003_expanded <- cp003_tcr %>%
  filter(!is.na(freq_num), freq_num > 1)

cat("Expanded cells (freq > 1):", nrow(cp003_expanded), "\n")
cat("Unique expanded clonotypes:", n_distinct(cp003_expanded$CTstrict), "\n")

# HIV-specific (your definition: expanded + not cluster 0/1)
cp003_hiv <- cp003_tcr %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR")

cat("HIV-specific cells:", nrow(cp003_hiv), "\n")
cat("Unique HIV-specific clonotypes:", n_distinct(cp003_hiv$CTstrict), "\n")

# --- 2C: HIV-specific clones by timepoint ---
hiv_by_tp <- cp003_hiv %>%
  group_by(Timepoint) %>%
  summarise(
    n_cells = n(),
    n_unique_clones = n_distinct(CTstrict),
    .groups = "drop"
  )

cat("\n=== HIV-specific clones by timepoint ===\n")
print(hiv_by_tp)

# Clone-level summary: which timepoints each clone appears in
hiv_clone_timepoints <- cp003_hiv %>%
  group_by(CTstrict) %>%
  summarise(
    timepoints = paste(sort(unique(Timepoint)), collapse = ";"),
    n_timepoints = n_distinct(Timepoint),
    total_cells = n(),
    max_freq = max(freq_num, na.rm = TRUE),
    clusters = paste(sort(unique(as.character(mnn_snn_res.0.6))), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(total_cells))

write.csv(hiv_clone_timepoints,
          file.path(ana_dir, "CP003_HIV_clones_by_timepoint.csv"),
          row.names = FALSE)

# =============================================================
# 3) TRACK HIV-SPECIFIC CLONES INTO TARA LONGITUDINAL DATA
#    This is the key validation step
# =============================================================

# Get the HIV-specific clone IDs from CP003
hiv_clone_ids <- unique(cp003_hiv$CTstrict)
all_cp003_clone_ids <- unique(cp003_tcr$CTstrict)

# Check which are found in TARA
tara_meta <- tara@meta.data %>%
  filter(!is.na(CTstrict), CTstrict != "")

# Matches
hiv_in_tara <- tara_meta %>%
  filter(CTstrict %in% hiv_clone_ids)

all_in_tara <- tara_meta %>%
  filter(CTstrict %in% all_cp003_clone_ids)

cat("\n=== Clone matching: CP003 → TARA ===\n")
cat("HIV-specific clones found in TARA:", n_distinct(hiv_in_tara$CTstrict),
    "of", length(hiv_clone_ids), "\n")
cat("All CP003 clones found in TARA:", n_distinct(all_in_tara$CTstrict),
    "of", length(all_cp003_clone_ids), "\n")
cat("HIV-specific cells in TARA:", nrow(hiv_in_tara), "\n")

# --- 3A: Longitudinal tracking of matched HIV-specific clones ---
# Parse TARA sample names to extract PID and timepoint from orig.ident
tara_hiv_tracking <- hiv_in_tara %>%
  mutate(
    PID = sub("_.*", "", orig.ident),
    TARA_Timepoint = case_when(
      grepl("_entry$",  orig.ident, ignore.case = TRUE) ~ "entry",
      grepl("_2m$",     orig.ident, ignore.case = TRUE) ~ "2m",
      grepl("_12m$",    orig.ident, ignore.case = TRUE) ~ "12m",
      grepl("_V12$",    orig.ident, ignore.case = TRUE) ~ "12m",
      grepl("_24m$",    orig.ident, ignore.case = TRUE) ~ "24m",
      grepl("_V24$",    orig.ident, ignore.case = TRUE) ~ "24m",
      grepl("_60m$",    orig.ident, ignore.case = TRUE) ~ "60m",
      grepl("_101m$",   orig.ident, ignore.case = TRUE) ~ "101m",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(CTstrict, PID, TARA_Timepoint) %>%
  summarise(
    n_cells = n(),
    clusters = paste(sort(unique(as.character(CD8_Annotation))), collapse = ";"),
    .groups = "drop"
  )

# Which CP003 timepoints each clone was originally found in (stimulation)
tara_hiv_tracking <- tara_hiv_tracking %>%
  left_join(
    hiv_clone_timepoints %>% select(CTstrict, cp003_timepoints = timepoints),
    by = "CTstrict"
  )

write.csv(tara_hiv_tracking,
          file.path(ana_dir, "HIV_clones_tracked_in_TARA.csv"),
          row.names = FALSE)

# --- 3B: Persistence summary ---
# For each matched clone: across how many TARA timepoints does it appear?
clone_persistence <- tara_hiv_tracking %>%
  filter(PID == "CP003") %>%  # CP003's own longitudinal samples in TARA
  
  group_by(CTstrict) %>%
  summarise(
    tara_timepoints = paste(sort(unique(TARA_Timepoint)), collapse = ";"),
    n_tara_timepoints = n_distinct(TARA_Timepoint),
    total_tara_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  left_join(
    hiv_clone_timepoints %>% select(CTstrict, cp003_timepoints = timepoints,
                                    cp003_cells = total_cells, max_freq),
    by = "CTstrict"
  ) %>%
  arrange(desc(n_tara_timepoints), desc(total_tara_cells))

write.csv(clone_persistence,
          file.path(ana_dir, "HIV_clone_persistence_CP003_longitudinal.csv"),
          row.names = FALSE)

cat("\n=== Clone persistence in TARA CP003 samples ===\n")
cat("Clones found at 1 TARA timepoint:", sum(clone_persistence$n_tara_timepoints == 1), "\n")
cat("Clones found at 2 TARA timepoints:", sum(clone_persistence$n_tara_timepoints == 2), "\n")
cat("Clones found at 3+ TARA timepoints:", sum(clone_persistence$n_tara_timepoints >= 3), "\n")

# =============================================================
# 4) RECLASSIFICATION: Fig 2C "unknown" vs functionally validated
#    Compare Trex epitope predictions with stimulation-based labels
# =============================================================

# From master match table: what was the Trex prediction for each matched clone?
reclassification <- master_match %>%
  filter(CP003_clone_status == "HIV-Specific Match") %>%
  mutate(
    # Was this clone predicted as HIV by Trex at any edit distance?
    trex_hiv = grepl("HIV", paste(CP003_Species_ED0, CP003_Species_ED1, CP003_Species_ED2,
                                  Tara_Species_ED0, Tara_Species_ED1, Tara_Species_ED2),
                     ignore.case = TRUE),
    # Was it predicted as something else?
    trex_other = !trex_hiv & nchar(paste0(
      CP003_Species_ED0, CP003_Species_ED1, CP003_Species_ED2,
      Tara_Species_ED0, Tara_Species_ED1, Tara_Species_ED2
    )) > 6,  # crude check for non-empty species
    # Classification for panel
    trex_classification = case_when(
      trex_hiv ~ "Trex: HIV predicted",
      trex_other ~ "Trex: Non-HIV predicted",
      TRUE ~ "Trex: Unknown/No match"
    )
  )

cat("\n=== Reclassification of HIV-specific clones ===\n")
cat("Functionally validated HIV-specific clones:\n")
print(table(reclassification$trex_classification))

write.csv(reclassification,
          file.path(ana_dir, "HIV_clone_reclassification.csv"),
          row.names = FALSE)

# =============================================================
# 5) PHENOTYPE OF MATCHED CLONES IN TARA (UNSTIMULATED)
#    Since CP003 stim data can't be compared on the same UMAP,
#    we characterize the matched clones in the TARA object
# =============================================================

# Label cells in TARA by their clone match status (already done, but refine)
# Focus on CP003 samples within TARA
tara$PID <- sub("_.*", "", tara$orig.ident)

# Subset to CP003 samples in TARA for phenotype analysis
tara_cp003 <- subset(tara, PID == "CP003")

cat("\n=== TARA CP003 subset ===\n")
cat("Total cells:", ncol(tara_cp003), "\n")
cat("With TCR:", sum(!is.na(tara_cp003$CTstrict) & tara_cp003$CTstrict != ""), "\n")

# --- 5A: Cluster distribution of HIV-matched vs non-matched clones ---
tara_cp003_tcr <- as.data.frame(tara_cp003@meta.data)

# Seurat v5 can store some metadata columns as list-columns; flatten them
list_cols <- sapply(tara_cp003_tcr, is.list)
if (any(list_cols)) {
  cat("Flattening list-columns:", paste(names(list_cols)[list_cols], collapse = ", "), "\n")
  for (col in names(list_cols)[list_cols]) {
    tara_cp003_tcr[[col]] <- unlist(tara_cp003_tcr[[col]])
  }
}

tara_cp003_tcr <- tara_cp003_tcr %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  mutate(
    match_status = case_when(
      CTstrict %in% hiv_clone_ids ~ "HIV-Specific (validated)",
      CTstrict %in% all_cp003_clone_ids ~ "CP003 match (non-HIV)",
      TRUE ~ "No CP003 match"
    ),
    match_status = factor(match_status,
                          levels = c("No CP003 match", "CP003 match (non-HIV)",
                                     "HIV-Specific (validated)"))
  )

cluster_by_status <- tara_cp003_tcr %>%
  dplyr::count(match_status, CD8_Annotation) %>%
  group_by(match_status) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write.csv(cluster_by_status,
          file.path(ana_dir, "TARA_CP003_cluster_by_match_status.csv"),
          row.names = FALSE)

# --- 5C: Differential gene expression prep ---
# Add refined match status to TARA CP003 object for downstream DE
tara_cp003$HIV_validated_status <- "No CP003 match"
matched_cells <- rownames(tara_cp003@meta.data)[
  !is.na(tara_cp003$CTstrict) & tara_cp003$CTstrict %in% hiv_clone_ids
]
tara_cp003$HIV_validated_status[matched_cells] <- "HIV-Specific (validated)"

nonhiv_cells <- rownames(tara_cp003@meta.data)[
  !is.na(tara_cp003$CTstrict) &
    tara_cp003$CTstrict %in% all_cp003_clone_ids &
    !(tara_cp003$CTstrict %in% hiv_clone_ids)
]
tara_cp003$HIV_validated_status[nonhiv_cells] <- "CP003 match (non-HIV)"

# Only run DE if enough cells in each group
de_counts <- table(tara_cp003$HIV_validated_status)
cat("\n=== Cells per group for DE ===\n")
print(de_counts)

if (de_counts["HIV-Specific (validated)"] >= 10) {
  Idents(tara_cp003) <- "HIV_validated_status"
  
  # HIV-validated vs all other TCR+ cells
  de_hiv_vs_other <- FindMarkers(
    tara_cp003,
    ident.1 = "HIV-Specific (validated)",
    ident.2 = c("CP003 match (non-HIV)", "No CP003 match"),
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  de_hiv_vs_other$gene <- rownames(de_hiv_vs_other)
  
  write.csv(de_hiv_vs_other,
            file.path(ana_dir, "DE_HIV_validated_vs_other_in_TARA_CP003.csv"),
            row.names = FALSE)
  
  cat("\nTop DE genes (HIV-validated vs other):\n")
  print(head(de_hiv_vs_other, 20))
} else {
  cat("\nToo few HIV-validated cells for DE analysis. Skipping.\n")
}

# =============================================================
# 6) CELL CYCLE COMPARISON — matched clones in TARA
# =============================================================
cc_by_status <- tara_cp003_tcr %>%
  filter(!is.na(Phase)) %>%
  dplyr::count(match_status, Phase) %>%
  group_by(match_status) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write.csv(cc_by_status,
          file.path(ana_dir, "CellCycle_by_match_status_TARA_CP003.csv"),
          row.names = FALSE)

# =============================================================
# 7) ALLUVIAL DATA PREP
#    Build the data structure for tracking HIV-specific clones
#    across ALL CP003 timepoints (TARA entry/v12/v24 + stim 2m/101m)
# =============================================================

# Combine clone presence from both datasets
# TARA side: CP003 samples
tara_cp003_clones <- tara_cp003@meta.data %>%
  filter(!is.na(CTstrict), CTstrict != "",
         CTstrict %in% hiv_clone_ids) %>%
  mutate(
    source = "TARA_unstimulated",
    timepoint_label = case_when(
      grepl("entry", orig.ident, ignore.case = TRUE) ~ "Entry (~1-2m)",
      grepl("V12|12m", orig.ident, ignore.case = TRUE) ~ "12m",
      grepl("V24|24m", orig.ident, ignore.case = TRUE) ~ "24m",
      grepl("60m", orig.ident, ignore.case = TRUE) ~ "60m",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(CTstrict, timepoint_label, source) %>%
  summarise(n_cells = n(), .groups = "drop")

# CP003 stimulation side
cp003_stim_clones <- cp003_hiv %>%
  mutate(
    source = "CP003_stimulated",
    timepoint_label = case_when(
      Timepoint == "2m" ~ "2m (stim)",
      Timepoint == "101m" ~ "101m (stim)",
      TRUE ~ Timepoint
    )
  ) %>%
  group_by(CTstrict, timepoint_label, source) %>%
  summarise(n_cells = n(), .groups = "drop")

# Combined
alluvial_data <- bind_rows(tara_cp003_clones, cp003_stim_clones)

# Set timepoint order
alluvial_data$timepoint_label <- factor(
  alluvial_data$timepoint_label,
  levels = c("Entry (~1-2m)", "2m (stim)", "12m", "24m", "101m (stim)")
)

write.csv(alluvial_data,
          file.path(ana_dir, "Alluvial_HIV_clones_combined_timepoints.csv"),
          row.names = FALSE)

# Summary: which clones persist across how many timepoints total?
total_persistence <- alluvial_data %>%
  group_by(CTstrict) %>%
  summarise(
    timepoints = paste(sort(unique(as.character(timepoint_label))), collapse = " → "),
    n_timepoints = n_distinct(timepoint_label),
    total_cells = sum(n_cells),
    sources = paste(sort(unique(source)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_timepoints), desc(total_cells))

write.csv(total_persistence,
          file.path(ana_dir, "HIV_clone_total_persistence_summary.csv"),
          row.names = FALSE)

cat("\n=== Total persistence across all sources ===\n")
cat("Clones at 1 timepoint:", sum(total_persistence$n_timepoints == 1), "\n")
cat("Clones at 2 timepoints:", sum(total_persistence$n_timepoints == 2), "\n")
cat("Clones at 3 timepoints:", sum(total_persistence$n_timepoints == 3), "\n")
cat("Clones at 4+ timepoints:", sum(total_persistence$n_timepoints >= 4), "\n")

# =============================================================
# 8) SAVE UPDATED OBJECTS
# =============================================================
cat("\nSaving updated TARA CP003 subset...\n")
qs_save(tara_cp003, file = file.path(saved_dir, "tara_cp003_fig3_annotated.qs2"))

cat("\n=== Analysis complete. Outputs in:", ana_dir, "===\n")

# =============================================================
# SESSION INFO
# =============================================================
cat("\n")
sessionInfo()

