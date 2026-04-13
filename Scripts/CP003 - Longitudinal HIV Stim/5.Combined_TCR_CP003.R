############################################################
# TARA TCR only:
# - Read all TARA TCR contigs
# - KEEP original CP003 TARA samples
# - ADD merged CP003 VDJ samples as separate unique samples
# - Add metadata labels
# - Make alluvial plot for CP003
############################################################

suppressPackageStartupMessages({
  library(scRepertoire)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------- #
# 0) Output directory
# ----------------------------- #
out.dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CP003/VDJ/TCR/Shared_Clones/Combined"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- #
# 1) Read all TARA TCR contigs
# ----------------------------- #
base.path <- path.expand("~/Documents/10x_Genomics/")

get_folder_names <- function(in.path) {
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  folder_names[folder_names != ""]
}

f_names <- get_folder_names(base.path)

# Only TARA samples
TARA_names <- f_names[grepl("^C", f_names)]

read_tcr_contig <- function(sample_name, base.path) {
  t_file <- file.path(base.path, sample_name, "per_sample_outs", "TCR", "filtered_contig_annotations.csv")
  if (!file.exists(t_file)) {
    warning("Missing TCR file for: ", sample_name)
    return(NULL)
  }
  read.csv(t_file, stringsAsFactors = FALSE)
}

tara_tcr_list <- lapply(TARA_names, read_tcr_contig, base.path = base.path)
names(tara_tcr_list) <- TARA_names
tara_tcr_list <- tara_tcr_list[!vapply(tara_tcr_list, is.null, logical(1))]

# ----------------------------- #
# 2) Read and merge CP003 technical replicates
# ----------------------------- #
cp003.path <- "/home/akshay-iyer/CP003_multi_out/outs/per_sample_outs"

t_2m_001A  <- read.csv(file.path(cp003.path, "CP003_2m_001A",  "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)
t_2m_001B  <- read.csv(file.path(cp003.path, "CP003_2m_001B",  "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)
t_101m_003 <- read.csv(file.path(cp003.path, "CP003_101m_003", "vdj_t", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)

t_2m_cd8p <- bind_rows(t_2m_001A, t_2m_001B)

cp003_merged_list <- list(
  CP003_2m_CD8plus    = t_2m_cd8p,
  CP003_101m_CD8plus  = t_101m_003
)

# ----------------------------- #
# 3) Final TARA list:
#    all original TARA + merged CP003
# ----------------------------- #
tara_tcr_final_list <- c(tara_tcr_list, cp003_merged_list)

# ----------------------------- #
# 4) Combine TCR
# ----------------------------- #
combined.TCR.TARA <- combineTCR(
  tara_tcr_final_list,
  samples = names(tara_tcr_final_list),
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

# ----------------------------- #CP003_RNA_integrated_CCA_MNN_withTCR.qs2
# 5) Add metadata to combined object
# ----------------------------- #
samples_all <- names(combined.TCR.TARA)

PID_vec <- sub("^([^_]+).*", "\\1", samples_all)

# Timepoint parsing
Timepoint_vec <- dplyr::case_when(
  grepl("^CP003_2m_CD8plus$", samples_all) ~ "2m",
  grepl("^CP003_101m_CD8plus$", samples_all) ~ "101m",
  grepl("_entry$", samples_all, ignore.case = TRUE) ~ "entry",
  grepl("_2m$", samples_all, ignore.case = TRUE) ~ "2m",
  grepl("_12m$", samples_all, ignore.case = TRUE) ~ "12m",
  grepl("_24m$", samples_all, ignore.case = TRUE) ~ "24m",
  grepl("_60m$", samples_all, ignore.case = TRUE) ~ "60m",
  grepl("_101m$", samples_all, ignore.case = TRUE) ~ "101m",
  grepl("_v12$", samples_all, ignore.case = TRUE) ~ "12m",
  grepl("_v24$", samples_all, ignore.case = TRUE) ~ "24m",
  TRUE ~ NA_character_
)


Sample_Type_vec <- dplyr::case_when(
  grepl("^CP003_(2m|101m)_CD8(plus|minus)$", samples_all) ~ "CP003_merged_VDJ",
  TRUE ~ "TARA_original"
)

combined.TCR.TARA <- addVariable(combined.TCR.TARA, "PID", PID_vec)
combined.TCR.TARA <- addVariable(combined.TCR.TARA, "Timepoint", Timepoint_vec)
combined.TCR.TARA <- addVariable(combined.TCR.TARA, "Sample_Type", Sample_Type_vec)
combined.TCR.TARA <- addVariable(combined.TCR.TARA, "Sample_Label", samples_all)

# ----------------------------- #
# 6) Clonal overlap on ALL samples
# ----------------------------- #
p_overlap <- clonalOverlap(
  combined.TCR.TARA,
  cloneCall = "strict",
  method = "raw",
  palette = 'viridis'
)

ggsave(
  filename = file.path(out.dir, "TARA_Clonal_Overlap_Strict_Raw.png"),
  plot = p_overlap,
  width = 35,
  height = 11,
  dpi = 300,
  bg = "white"
)

# ----------------------------- #
# 7) Clonal compare for CP003 only
#    Keep CP003 entry/v12/v24 + added CD8+ samples
#    Exclude CD8- sample
# ----------------------------- #
cp003_samples <- c(
  "CP003_2m_CD8plus",
  "CP003_entry",
  "CP003_V12",
  "CP003_V24",
  "CP003_101m_CD8plus"
)

cp003_samples <- cp003_samples[cp003_samples %in% names(combined.TCR.TARA)]

p_compare <- clonalCompare(
  combined.TCR.TARA,
  top.clones = 20,
  samples = cp003_samples,
  order.by = cp003_samples,
  cloneCall = "strict",
  relabel.clones = TRUE,
  proportion = FALSE,
  graph = "alluvial"
)

ggsave(
  filename = file.path(out.dir, "CP003_Clonal_Comparison_Alluvial.png"),
  plot = p_compare,
  width = 15,
  height = 11,
  dpi = 300,
  bg = "white"
)
p_compare_2 <- clonalCompare(
  combined.TCR.TARA,
  top.clones = 20,
  samples = cp003_samples,
  order.by = cp003_samples,
  cloneCall = "strict",
  relabel.clones = TRUE,
  proportion = FALSE,
  graph = "alluvial",
  exportTable = T
)

write.csv(
  p_compare_2,
  file.path(out.dir, "CP003_Clonal_Comparison_Alluvial.csv"),
  row.names = FALSE
)
