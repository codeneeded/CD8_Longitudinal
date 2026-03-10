############################################################
# CP003 HIV-Specific Clone Matching + Cell Cycle + Trex Epitope
# Across CP003 (seu) and tara_cdnk
############################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)
library(Trex)

# ----------------------------- #
# 0) Paths
# ----------------------------- #
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
out_dir   <- file.path(base_dir, "CP003/Shared_clone_Tables_Bulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qs_cp003_in  <- file.path(saved_dir, "seu_CP003_HIVSpecificTCR_annotated.qs")
rds_tara_in  <- file.path(saved_dir, "tara_cdnk_annotated.rds")
rds_cp003_out <- file.path(saved_dir, "cp003_with_cc_epitope.qs")
rds_tara_out  <- file.path(saved_dir, "tara_cdnk_with_cp003_matches.rds")

# ----------------------------- #
# 1) Load objects
# ----------------------------- #
cp003     <- qs2::qs_read(qs_cp003_in)
tara_cdnk <- readRDS(rds_tara_in)

# ----------------------------- #
# 2) Cell cycle scoring — CP003
# ----------------------------- #
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cp003 <- JoinLayers(cp003)
cp003 <- CellCycleScoring(
  cp003,
  s.features   = s.genes,
  g2m.features = g2m.genes,
  set.ident    = FALSE
)


# ----------------------------- #
# 4) Trex epitope annotation — CP003
# ----------------------------- #
cp003_trex_0 <- annotateDB(cp003, chains = "TRB")
cp003_trex_1 <- annotateDB(cp003, chains = "TRB", edit.distance = 1)
cp003_trex_2 <- annotateDB(cp003, chains = "TRB", edit.distance = 2)

# Pull TRB columns per edit distance and label them
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

cp003_trex_df <- pull_trex_meta(cp003_trex_0, "ED0") %>%
  left_join(pull_trex_meta(cp003_trex_1, "ED1"), by = c("Cell", "CTstrict")) %>%
  left_join(pull_trex_meta(cp003_trex_2, "ED2"), by = c("Cell", "CTstrict"))

# Add Trex results back to CP003 metadata
cp003 <- AddMetaData(cp003, cp003_trex_df %>% tibble::column_to_rownames("Cell"))

# ----------------------------- #
# 5) Trex epitope annotation — tara_cdnk
# ----------------------------- #
tara_trex_0 <- annotateDB(tara_cdnk, chains = "TRB")
tara_trex_1 <- annotateDB(tara_cdnk, chains = "TRB", edit.distance = 1)
tara_trex_2 <- annotateDB(tara_cdnk, chains = "TRB", edit.distance = 2)

tara_trex_df <- pull_trex_meta(tara_trex_0, "ED0") %>%
  left_join(pull_trex_meta(tara_trex_1, "ED1"), by = c("Cell", "CTstrict")) %>%
  left_join(pull_trex_meta(tara_trex_2, "ED2"), by = c("Cell", "CTstrict"))

tara_cdnk <- AddMetaData(tara_cdnk, tara_trex_df %>% tibble::column_to_rownames("Cell"))

# ----------------------------- #
# 6) Extract HIV-specific clone list from CP003
# ----------------------------- #
cp003_meta <- cp003@meta.data %>%
  transmute(
    Cell                 = rownames(.),
    CTstrict             = CTstrict,
    HIV_Specific_TCR     = HIV_Specific_TCR,
    Sample               = Sample,
    CP003_Cluster        = mnn_snn_res.0.6,
    clonalFrequency      = clonalFrequency,
    CP003_Phase          = Phase,
    CP003_S_Score        = S.Score,
    CP003_G2M_Score      = G2M.Score,
    CP003_Epitope_ED0    = Epitope.target_ED0,
    CP003_Epitope_ED1    = Epitope.target_ED1,
    CP003_Epitope_ED2    = Epitope.target_ED2,
    CP003_Species_ED0    = Epitope.species_ED0,
    CP003_Species_ED1    = Epitope.species_ED1,
    CP003_Species_ED2    = Epitope.species_ED2
  )

# Unique HIV-specific clones
cp003_hiv_clones <- cp003_meta %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR", !is.na(CTstrict), CTstrict != "") %>%
  distinct(CTstrict) %>%
  pull(CTstrict)

cat("Unique HIV-specific CP003 clones:", length(cp003_hiv_clones), "\n")

# All CP003 clones
cp003_all_clones <- cp003_meta %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  distinct(CTstrict) %>%
  pull(CTstrict)

cat("Total CP003 clones:", length(cp003_all_clones), "\n")

# Clone-level summary from CP003 (for joining onto match table later)
cp003_clone_origin <- cp003_meta %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict) %>%
  summarise(
    CP003_HIV_Specific      = paste(sort(unique(HIV_Specific_TCR)), collapse = ";"),
    CP003_Samples           = paste(sort(unique(Sample)),           collapse = ";"),
    CP003_Clusters          = paste(sort(unique(CP003_Cluster)),    collapse = ";"),
    CP003_n_cells           = n(),
    CP003_max_clonalFreq    = max(suppressWarnings(as.numeric(as.character(clonalFrequency))), na.rm = TRUE),
    CP003_Phase_breakdown   = paste(sort(unique(CP003_Phase)),      collapse = ";"),
    CP003_pct_S             = round(mean(CP003_Phase == "S",   na.rm = TRUE) * 100, 1),
    CP003_pct_G2M           = round(mean(CP003_Phase == "G2M", na.rm = TRUE) * 100, 1),
    CP003_pct_G1            = round(mean(CP003_Phase == "G1",  na.rm = TRUE) * 100, 1),
    CP003_mean_S_Score      = round(mean(CP003_S_Score,   na.rm = TRUE), 4),
    CP003_mean_G2M_Score    = round(mean(CP003_G2M_Score, na.rm = TRUE), 4),
    CP003_Epitope_ED0       = paste(sort(unique(na.omit(CP003_Epitope_ED0))), collapse = ";"),
    CP003_Epitope_ED1       = paste(sort(unique(na.omit(CP003_Epitope_ED1))), collapse = ";"),
    CP003_Epitope_ED2       = paste(sort(unique(na.omit(CP003_Epitope_ED2))), collapse = ";"),
    CP003_Species_ED0       = paste(sort(unique(na.omit(CP003_Species_ED0))), collapse = ";"),
    CP003_Species_ED1       = paste(sort(unique(na.omit(CP003_Species_ED1))), collapse = ";"),
    CP003_Species_ED2       = paste(sort(unique(na.omit(CP003_Species_ED2))), collapse = ";"),
    .groups = "drop"
  )
# ----------------------------- #
# 7) Label matching clones in tara_cdnk
# ----------------------------- #
tara_meta <- tara_cdnk@meta.data %>%
  transmute(
    Cell              = rownames(.),
    CTstrict          = CTstrict,
    Sample            = orig.ident,
    Tara_Cluster      = T_Cell_A,
    clonalFrequency   = clonalFrequency,
    Tara_Phase        = Phase,
    Tara_S_Score      = S.Score,
    Tara_G2M_Score    = G2M.Score,
    Tara_Epitope_ED0  = Epitope.target_ED0,
    Tara_Epitope_ED1  = Epitope.target_ED1,
    Tara_Epitope_ED2  = Epitope.target_ED2,
    Tara_Species_ED0  = Epitope.species_ED0,
    Tara_Species_ED1  = Epitope.species_ED1,
    Tara_Species_ED2  = Epitope.species_ED2
  ) %>%
  mutate(
    CP003_clone_status = case_when(
      !is.na(CTstrict) & CTstrict %in% cp003_hiv_clones ~ "HIV-Specific Match",
      !is.na(CTstrict) & CTstrict %in% cp003_all_clones ~ "Non-HIV Match",
      TRUE ~ "Other"
    ),
    CP003_clone_status = factor(
      CP003_clone_status,
      levels = c("Other", "Non-HIV Match", "HIV-Specific Match")
    )
  )

# Write back to tara_cdnk object
tara_cdnk$CP003_clone_status <- tara_meta$CP003_clone_status
tara_cdnk$Phase              <- tara_meta$Tara_Phase
tara_cdnk$S.Score            <- tara_meta$Tara_S_Score
tara_cdnk$G2M.Score          <- tara_meta$Tara_G2M_Score

# ----------------------------- #
# 8) Build master match table
# ----------------------------- #

# Tara-side clone summary (only matched clones)
tara_clone_summary <- tara_meta %>%
  filter(CP003_clone_status != "Other", !is.na(CTstrict), CTstrict != "") %>%
  group_by(CTstrict, CP003_clone_status) %>%
  summarise(
    Tara_Samples        = paste(sort(unique(trimws(as.character(Sample)))), collapse = ";"),
    Tara_n_samples      = n_distinct(Sample),
    Tara_Clusters       = paste(sort(unique(Tara_Cluster)), collapse = ";"),
    Tara_n_cells        = n(),
    Tara_max_clonalFreq = max(suppressWarnings(as.numeric(as.character(clonalFrequency))), na.rm = TRUE),
    Tara_pct_S          = round(mean(Tara_Phase == "S",   na.rm = TRUE) * 100, 1),
    Tara_pct_G2M        = round(mean(Tara_Phase == "G2M", na.rm = TRUE) * 100, 1),
    Tara_pct_G1         = round(mean(Tara_Phase == "G1",  na.rm = TRUE) * 100, 1),
    Tara_mean_S_Score   = round(mean(Tara_S_Score,   na.rm = TRUE), 4),
    Tara_mean_G2M_Score = round(mean(Tara_G2M_Score, na.rm = TRUE), 4),
    Tara_Epitope_ED0    = paste(sort(unique(na.omit(Tara_Epitope_ED0))), collapse = ";"),
    Tara_Epitope_ED1    = paste(sort(unique(na.omit(Tara_Epitope_ED1))), collapse = ";"),
    Tara_Epitope_ED2    = paste(sort(unique(na.omit(Tara_Epitope_ED2))), collapse = ";"),
    Tara_Species_ED0    = paste(sort(unique(na.omit(Tara_Species_ED0))), collapse = ";"),
    Tara_Species_ED1    = paste(sort(unique(na.omit(Tara_Species_ED1))), collapse = ";"),
    Tara_Species_ED2    = paste(sort(unique(na.omit(Tara_Species_ED2))), collapse = ";"),
    .groups = "drop"
  )

# Join CP003 origin info and merge sample columns
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

# ----------------------------- #
# Species-based clone status classification
# ----------------------------- #

species_cols <- c(
  "Tara_Species_ED0", "Tara_Species_ED1", "Tara_Species_ED2",
  "CP003_Species_ED0", "CP003_Species_ED1", "CP003_Species_ED2"
)

master_match_table <- master_match_table %>%
  mutate(
    # Check species predictions
    HIV_Predicted = apply(
      select(., all_of(species_cols)), 1,
      function(x) any(grepl("HIV", x, ignore.case = TRUE))
    ),
    NonHIV_Predicted = apply(
      select(., all_of(species_cols)), 1,
      function(x) {
        vals <- trimws(unlist(strsplit(paste(na.omit(x), collapse = ";"), ";")))
        vals <- vals[vals != ""]
        any(!grepl("HIV|Unknown", vals, ignore.case = TRUE))
      }
    ),
    NonHIV_Species = apply(
      select(., all_of(species_cols)), 1,
      function(x) {
        vals <- trimws(unlist(strsplit(paste(na.omit(x), collapse = ";"), ";")))
        vals <- vals[vals != ""]
        non_hiv <- unique(vals[!grepl("HIV|Unknown", vals, ignore.case = TRUE)])
        if (length(non_hiv) == 0) NA_character_ else paste(sort(non_hiv), collapse = ";")
      }
    ),
    # Assign clone status in priority order:
    # 1) HIV predicted in species           -> HIV-Specific Match
    # 2) Non-HIV predicted in species       -> Non-HIV Match
    # 3) No species prediction but CP003 expanded (n_cells > 1) and Tara species blank -> HIV-Specific Match
    # 4) Everything else                    -> keep existing status
    CP003_clone_status = case_when(
      HIV_Predicted                                                       ~ "HIV-Specific Match",
      NonHIV_Predicted                                                    ~ "Non-HIV Match",
      CP003_n_cells > 1 &
        (is.na(Tara_Species_ED0) | Tara_Species_ED0 == "") &
        (is.na(Tara_Species_ED1) | Tara_Species_ED1 == "") &
        (is.na(Tara_Species_ED2) | Tara_Species_ED2 == "")               ~ "HIV-Specific Match",
      TRUE                                                                ~ as.character(CP003_clone_status)
    ),
    CP003_clone_status = factor(
      CP003_clone_status,
      levels = c("Non-HIV Match", "HIV-Specific Match")
    )
  ) %>%
  arrange(desc(CP003_clone_status), desc(Tara_n_cells))
# ----------------------------- #
# 9) Save output tables
# ----------------------------- #

# Master table — all clones
write.csv(
  master_match_table,
  file = file.path(out_dir, "Master_Clone_Match_CP003_vs_Tara.csv"),
  row.names = FALSE
)

# HIV-specific only (CMV-flagged clones excluded)
write.csv(
  master_match_table %>% filter(CP003_clone_status == "HIV-Specific Match"),
  file = file.path(out_dir, "HIV_Specific_Clones_CP003_in_Tara.csv"),
  row.names = FALSE
)


# Per-cell table for matched clones in tara
write.csv(
  tara_meta %>% filter(CP003_clone_status != "Other"),
  file = file.path(out_dir, "Tara_Matched_Cells_per_cell.csv"),
  row.names = FALSE
)

# CP003 clone origin summary
write.csv(
  cp003_clone_origin,
  file = file.path(out_dir, "CP003_Clone_Origin_Summary.csv"),
  row.names = FALSE
)


# ----------------------------- #
# 10) Save updated Seurat objects
# ----------------------------- #
qs2::qs_save(cp003,     file = rds_cp003_out)
saveRDS(tara_cdnk,      file = rds_tara_out)

cat("Done. Outputs saved to:", out_dir, "\n")