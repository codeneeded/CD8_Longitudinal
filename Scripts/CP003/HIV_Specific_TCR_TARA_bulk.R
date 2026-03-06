############################################################
# Find CP003 HIV-expanding clones inside tara_cdnk
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs2)
})

# ----------------------------- #
# 0) Paths
# ----------------------------- #
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
out_dir <- file.path(base_dir, "CP003/Shared_clone_Tables_Bulk")
rds_in     <- file.path(saved_dir, "tara_cdnk_annotated.rds")

# EDIT THIS IF YOUR CP003 qs filename is different
qs_cp003_in <- file.path(saved_dir, "seu_CP003_HIVSpecificTCR_annotated.qs")

rds_out    <- file.path(saved_dir, "tara_cdnk_annotated_with_CP003_HIV_clones.rds")

# ----------------------------- #
# 1) Read objects
# ----------------------------- #
tara_cdnk <- readRDS(rds_in)
cp003     <- qs2::qs_read(qs_cp003_in)

# ----------------------------- #
# 2) Pull HIV-specific clones from CP003
# ----------------------------- #
# Assumes these columns already exist in CP003:
# cp003$HIV_Specific_TCR
# cp003$CTstrict

cp003_meta <- data.frame(
  Cell = colnames(cp003),
  CTstrict = cp003$CTstrict,
  HIV_Specific_TCR = cp003$HIV_Specific_TCR,
  Sample = cp003$Sample,
  Cluster = cp003$mnn_snn_res.0.6,
  clonalFrequency = cp003$clonalFrequency,
  stringsAsFactors = FALSE
)

cp003_hiv_clones <- cp003_meta %>%
  filter(
    HIV_Specific_TCR == "HIV-Specific TCR",
    !is.na(CTstrict),
    CTstrict != ""
  ) %>%
  distinct(CTstrict) %>%
  pull(CTstrict)

cat("Number of unique HIV-specific CP003 clones:", length(cp003_hiv_clones), "\n")

# All CP003 clones (regardless of HIV specificity)
cp003_all_clones <- cp003_meta %>%
  filter(
    !is.na(CTstrict),
    CTstrict != ""
  ) %>%
  distinct(CTstrict) %>%
  pull(CTstrict)

cat("Number of total CP003 clones:", length(cp003_all_clones), "\n")
# ----------------------------- #
# 3) Label matching cells in tara_cdnk
# ----------------------------- #

tara_meta <- data.frame(
  Cell = colnames(tara_cdnk),
  CTstrict = tara_cdnk$CTstrict,
  Sample = tara_cdnk$orig.ident,
  Cluster = tara_cdnk$T_Cell_A,
  clonalFrequency = tara_cdnk$clonalFrequency,
  stringsAsFactors = FALSE
)

tara_meta$CP003_clone_status <- ifelse(
  !is.na(tara_meta$CTstrict) & tara_meta$CTstrict %in% cp003_hiv_clones,
  "HIV specific Match",
  ifelse(
    !is.na(tara_meta$CTstrict) & tara_meta$CTstrict %in% cp003_all_clones,
    "Match",
    "Other"
  )
)

tara_meta$CP003_clone_status <- factor(
  tara_meta$CP003_clone_status,
  levels = c("Other", "Match", "HIV specific Match")
)

tara_cdnk$CP003_clone_status <- tara_meta$CP003_clone_status
# -----------------------------
# Dataframe of CP003 clone matches by cluster
# -----------------------------

match_df <- data.frame(
  Clone_Status = tara_cdnk$CP003_clone_status,
  Cluster = tara_cdnk$T_Cell_A,
  CTstrict = tara_cdnk$CTstrict,
  Sample = tara_cdnk$orig.ident,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(
    Clone_Status %in% c("Match", "HIV specific Match"),
    !is.na(CTstrict),
    CTstrict != ""
  )

# collapse to clone-level counts similar to table()
match_summary <- match_df |>
  dplyr::group_by(Clone_Status, CTstrict, Cluster) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    Samples = paste(sort(unique(Sample)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::arrange(Clone_Status, CTstrict, Cluster)

# save as csv
write.csv(
  match_summary,
  file = file.path(out_dir, "tara_cdnk_clone_match_summary.csv"),
  row.names = FALSE
)
