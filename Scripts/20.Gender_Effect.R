library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(SeuratExtend)
library(scCustomize)
library(tidyr)
library(stringr)

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "TARA_ALL_post_annotation.rds")

TARA_ALL <- readRDS(rds_in)
DefaultAssay(TARA_ALL) <- 'RNA'

############################### Add Gender to metadata ########################


# 1) Build the lookup table from what you pasted
sex_key <- tibble::tribble(
  ~PID,    ~Sex,
  "CE001", "F",
  "CE002", "F",
  "CE003", "M",
  "CE004", "F",
  "CE005", "M",
  "CE006", "M",
  "CE007", "M",
  "CE008", "F",
  "CE009", "F",
  "CE010", "F",
  "CE011", "F",
  "CE012", "F",
  "CE013", "M",
  "CE014", "F",
  "CE015", "M",
  "CE016", "M",
  "CE017", "F",
  "CE018", "M",
  "CE019", "F",
  "CE020", "M",
  "CE021", "F",
  "CE022", "M",
  "CE023", "F",
  "CE024", "M",
  "CE025", "F",
  "CE026", "M",
  "CE027", "M",
  "CE028", "F",
  "CE030", "F",
  "CE031", "M",
  "CE032", "M",
  "CE033", "M",
  "CE034", "M",
  "CE035", "M",
  "CE036", "F",
  "CE037", "M",
  "CE038", "F",
  "CE039", "M",
  "CE040", "M",
  "CE041", "M",
  "CE042", "F",
  "CE043", "F",
  "CE044", "F",
  "CE045", "F",
  "CE046", "F",
  "CE048", "F",
  "CP001", "F",
  "CP002", "F",
  "CP003", "M",
  "CP004", "M",
  "CP005", "F",
  "CP006", "F",
  "CP007", "F",
  "CP008", "M",
  "CP009", "M",
  "CP010", "F",
  "CP011", "F",
  "CP012", "F",
  "CP013", "M",
  "CP014", "F",
  "CP015", "M",
  "CP016", "F",
  "CP017", "F",
  "CP018", "F",
  "CP019", "F",
  "CP020", "F",
  "CP021", "M",
  "CP022", "F",
  "CP023", "F",
  "CP024", "M",
  "CP025", "F",
  "CP026", "M",
  "CP027", "F",
  "CP028", "M",
  "CP029", "F",
  "CP030", "F",
  "CP031", "M",
  "CP032", "F",
  "CP033", "M",
  "CP034", "M",
  "CP035", "M",
  "CP036", "F",
  "CP037", "F",
  "CP038", "M",
  "CP039", "F",
  "CP040", "M",
  "CP041", "F",
  "CP042", "M",
  "CP043", "F",
  "CS005", "F",
  "CS015", "F"
)

# 2) Extract PID from orig.ident (everything before the first underscore)
TARA_ALL$PID_from_orig <- sub("_.*$", "", TARA_ALL$orig.ident)

# 3) Keep only the PIDs that are actually in your experiment/object
pids_in_obj <- sort(unique(TARA_ALL$PID_from_orig))
sex_key_rel <- sex_key %>% filter(PID %in% pids_in_obj)

# (Optional sanity check) see if any PIDs in object are missing from the sex table
missing_pids <- setdiff(pids_in_obj, sex_key$PID)
if (length(missing_pids) > 0) {
  warning("These PIDs are in TARA_ALL but missing from sex_key: ",
          paste(missing_pids, collapse = ", "))
}

# 4) Join Sex onto meta.data (by PID)
TARA_ALL@meta.data <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(sex_key_rel, by = c("PID_from_orig" = "PID")) %>%
  dplyr::rename(Gender = Sex) %>%   # or keep as Sex if you prefer
  tibble::column_to_rownames("cell")

library(dplyr)

# ---------------------------- #
# 0) Sanity check columns exist
# ---------------------------- #
stopifnot(all(c("Condition","Age","Viral_Load") %in% colnames(TARA_ALL@meta.data)))

# ---------------------------- #
# 1) Build unified Timepoint_Group
# ---------------------------- #
TARA_ALL$Timepoint_Group <- NA_character_

# HEU/HUU labeled directly
TARA_ALL$Timepoint_Group[TARA_ALL$Condition == "HEU"] <- "HEU"
TARA_ALL$Timepoint_Group[TARA_ALL$Condition == "HUU"] <- "HUU"

# HEI stratified by Age then Viral_Load
is_HEI <- TARA_ALL$Condition == "HEI"

TARA_ALL$Timepoint_Group[is_HEI] <- dplyr::case_when(
  # Pre-ART entry: <2 months (your rule)
  !is.na(TARA_ALL$Age[is_HEI]) & TARA_ALL$Age[is_HEI] < 2 ~ "PreART_Entry",
  
  # Post-ART suppressed/unsuppressed only when Age >= 2
  !is.na(TARA_ALL$Age[is_HEI]) & TARA_ALL$Age[is_HEI] >= 2 &
    !is.na(TARA_ALL$Viral_Load[is_HEI]) & TARA_ALL$Viral_Load[is_HEI] < 200 ~ "PostART_Suppressed",
  
  !is.na(TARA_ALL$Age[is_HEI]) & TARA_ALL$Age[is_HEI] >= 2 &
    !is.na(TARA_ALL$Viral_Load[is_HEI]) & TARA_ALL$Viral_Load[is_HEI] >= 200 ~ "PostART_Unsuppressed",
  
  TRUE ~ NA_character_
)

# ---------------------------- #
# 2) Factor order
# ---------------------------- #
TARA_ALL$Timepoint_Group <- factor(
  TARA_ALL$Timepoint_Group,
  levels = c("HEU","HUU","PreART_Entry","PostART_Suppressed","PostART_Unsuppressed")
)

# ---------------------------- #
# 3) QC tables
# ---------------------------- #
message("\nCounts by Condition x Timepoint_Group:")
print(table(TARA_ALL$Condition, TARA_ALL$Timepoint_Group, useNA = "ifany"))

message("\nCounts by Timepoint_Group:")
print(table(TARA_ALL$Timepoint_Group, useNA = "ifany"))


############## Gender effect on TCR ################
############################################################
# Gender effect on clonal expansion (clusters 8,9,27)
# - Whole UMAP with DimPlot2 overlays
# - Sample-level summaries (correct testing unit)
# - Wilcoxon + BH correction
# - Save plots + CSVs
############################################################

# ---------------------------- #
# 0) Output dirs
# ---------------------------- #
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Gender_Effect"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# 1) Define clusters of interest
# ---------------------------- #
clusters_of_interest <- c("8: CTL-like", "9: TRDV1+ CTL-like", "27: GZMK+ CD8 T cell")

# ---------------------------- #
# 2) Build expansion labels (no subsetting of object)
# ---------------------------- #
TARA_ALL$has_TCR  <- !is.na(TARA_ALL$CTstrict) & TARA_ALL$CTstrict != ""
TARA_ALL$expanded <- ifelse(
  TARA_ALL$has_TCR & !is.na(TARA_ALL$clonalFrequency) & TARA_ALL$clonalFrequency >= 2,
  "Expanded", "Not expanded"
)
TARA_ALL$expanded <- factor(TARA_ALL$expanded, levels = c("Not expanded","Expanded"))

TARA_ALL$expanded_COI <- "Other"
TARA_ALL$expanded_COI[TARA_ALL$Manual_Annotation %in% clusters_of_interest] <-
  as.character(TARA_ALL$expanded[TARA_ALL$Manual_Annotation %in% clusters_of_interest])
TARA_ALL$expanded_COI <- factor(TARA_ALL$expanded_COI, levels = c("Other","Not expanded","Expanded"))

# ---------------------------- #
# 3) Whole UMAP plots (DimPlot2)
# ---------------------------- #
# A) Highlight clusters 8/9/27 vs Other
TARA_ALL$COI_flag <- ifelse(TARA_ALL$Manual_Annotation %in% clusters_of_interest, "8/9/27", "Other")
TARA_ALL$COI_flag <- factor(TARA_ALL$COI_flag, levels = c("Other","8/9/27"))

p_umap_coi <- DimPlot2(
  TARA_ALL,
  cols = "default",
  group.by = "COI_flag",
  reduction = "wnn.umap"
) + ggtitle("Whole UMAP: clusters 8/9/27 highlighted")

ggsave(
  filename = file.path(out_dir, "plots", "UMAP_whole_COI_highlight.pdf"),
  plot = p_umap_coi, width = 8, height = 6
)
ggsave(
  filename = file.path(out_dir, "plots", "UMAP_whole_COI_highlight.png"),
  plot = p_umap_coi, width = 8, height = 6, dpi = 300
)

# B) Expansion overlay (8/9/27 only colored; everything else "Other")
p_umap_exp <- DimPlot2(
  TARA_ALL,
  cols = "default",
  group.by = "expanded_COI",
  reduction = "wnn.umap"
) + ggtitle("Whole UMAP: clonal expansion overlay (clusters 8/9/27 only)")

ggsave(
  filename = file.path(out_dir, "plots", "UMAP_whole_expansion_overlay_COI.pdf"),
  plot = p_umap_exp, width = 8, height = 6
)
ggsave(
  filename = file.path(out_dir, "plots", "UMAP_whole_expansion_overlay_COI.png"),
  plot = p_umap_exp, width = 8, height = 6, dpi = 300
)

# C) Same overlay split by Gender (if you want it)
p_umap_exp_split <- DimPlot2(
  TARA_ALL,
  cols = "default",
  group.by = "expanded_COI",
  reduction = "wnn.umap",
  split.by = "Gender"
) + ggtitle("Whole UMAP: expansion overlay split by Gender")

ggsave(
  filename = file.path(out_dir, "plots", "UMAP_whole_expansion_overlay_COI_splitByGender.pdf"),
  plot = p_umap_exp_split, width = 12, height = 6
)
ggsave(
  filename = file.path(out_dir, "plots", "UMAP_whole_expansion_overlay_COI_splitByGender.png"),
  plot = p_umap_exp_split, width = 12, height = 6, dpi = 300
)

# ---------------------------- #
# 4) Sample-level summaries (orig.ident is the unit)
# ---------------------------- #
# NOTE: This avoids pseudoreplication.
# Metrics computed within each orig.ident x cluster:
# - n_cells
# - n_tcr, pct_tcr
# - pct_expanded_among_tcr (clone size >=2)
# - mean/median clone size among TCR+ cells
df_sum <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(Manual_Annotation %in% clusters_of_interest) %>%
  group_by(orig.ident, Gender, Manual_Annotation) %>%
  summarise(
    n_cells = n(),
    n_tcr   = sum(has_TCR, na.rm = TRUE),
    pct_tcr = n_tcr / n_cells,
    pct_expanded_among_tcr = ifelse(
      n_tcr > 0,
      sum(expanded == "Expanded", na.rm = TRUE) / n_tcr,
      NA_real_
    ),
    mean_clone_size_tcr   = ifelse(
      n_tcr > 0,
      mean(clonalFrequency[has_TCR], na.rm = TRUE),
      NA_real_
    ),
    median_clone_size_tcr = ifelse(
      n_tcr > 0,
      median(clonalFrequency[has_TCR], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

write.csv(df_sum, file.path(out_dir, "tables", "SampleLevel_ClonalExpansion_Summary_COI.csv"), row.names = FALSE)

# ---------------------------- #
# 5) Gender tests (Wilcoxon) per cluster + BH correction
# ---------------------------- #
tests_pct <- df_sum %>%
  filter(!is.na(pct_expanded_among_tcr)) %>%
  group_by(Manual_Annotation) %>%
  summarise(
    nF = sum(Gender == "F"),
    nM = sum(Gender == "M"),
    p  = ifelse(nF >= 2 & nM >= 2,
                wilcox.test(pct_expanded_among_tcr ~ Gender)$p.value,
                NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

tests_mean <- df_sum %>%
  filter(!is.na(mean_clone_size_tcr)) %>%
  group_by(Manual_Annotation) %>%
  summarise(
    nF = sum(Gender == "F"),
    nM = sum(Gender == "M"),
    p  = ifelse(nF >= 2 & nM >= 2,
                wilcox.test(mean_clone_size_tcr ~ Gender)$p.value,
                NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

write.csv(tests_pct,  file.path(out_dir, "tables", "Wilcoxon_GenderEffect_pctExpandedAmongTCR.csv"), row.names = FALSE)
write.csv(tests_mean, file.path(out_dir, "tables", "Wilcoxon_GenderEffect_meanCloneSizeAmongTCR.csv"), row.names = FALSE)

# ---------------------------- #
# 6) Gender effect plots (sample-level)
# ---------------------------- #
# A) % Expanded among TCR+ cells
p_pct <- ggplot(df_sum, aes(x = Gender, y = pct_expanded_among_tcr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
  facet_wrap(~ Manual_Annotation, scales = "free_y") +
  ylab("% expanded among TCR+ cells (clone size ≥ 2)") +
  xlab("Gender") +
  ggtitle("Gender effect on clonal expansion (sample-level summaries)")

ggsave(file.path(out_dir, "plots", "GenderEffect_pctExpandedAmongTCR_COI.png"),
       p_pct, width = 10, height = 5, dpi = 300)

# B) Mean clone size among TCR+ cells
p_mean <- ggplot(df_sum, aes(x = Gender, y = mean_clone_size_tcr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
  facet_wrap(~ Manual_Annotation, scales = "free_y") +
  ylab("Mean clonalFrequency among TCR+ cells") +
  xlab("Gender") +
  ggtitle("Gender effect on mean clone size (sample-level summaries)")


ggsave(file.path(out_dir, "plots", "GenderEffect_meanCloneSizeAmongTCR_COI.png"),
       p_mean, width = 10, height = 5, dpi = 300)

# ---------------------------- #
# 7) Quick console summary
# ---------------------------- #
message("\nDONE.\nOutputs saved to: ", out_dir,
        "\n- plots/ : UMAP overlays + gender effect plots",
        "\n- tables/: sample-level summaries + Wilcoxon results\n")
print(tests_pct)
print(tests_mean)
TARA_ALL$Age
levels(as.factor(TARA_ALL$Condition))

############################# 
############################################################
# Gender effect on clonal expansion
# Stratified by Timepoint_Group
############################################################

# ---------------------------- #
# 0) Output directory
# ---------------------------- #
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Gender_Effect/Stratified_By_Timepoint"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

clusters_of_interest <- c("8: CTL-like", "9: TRDV1+ CTL-like", "27: GZMK+ CD8 T cell")

# ---------------------------- #
# 1) Sample-level summaries
# ---------------------------- #
df_sum <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(Manual_Annotation %in% clusters_of_interest) %>%
  group_by(orig.ident, Gender, Timepoint_Group, Manual_Annotation) %>%
  summarise(
    n_cells = n(),
    n_tcr   = sum(has_TCR, na.rm = TRUE),
    pct_tcr = n_tcr / n_cells,
    pct_expanded_among_tcr = ifelse(
      n_tcr > 0,
      sum(expanded == "Expanded", na.rm = TRUE) / n_tcr,
      NA_real_
    ),
    mean_clone_size_tcr = ifelse(
      n_tcr > 0,
      mean(clonalFrequency[has_TCR], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

write.csv(df_sum,
          file.path(out_dir, "tables", "SampleLevel_ClonalExpansion_StratifiedByTimepoint.csv"),
          row.names = FALSE)

# ---------------------------- #
# 2) Wilcoxon tests per cluster x timepoint
# ---------------------------- #
tests_pct <- df_sum %>%
  filter(!is.na(pct_expanded_among_tcr)) %>%
  group_by(Timepoint_Group, Manual_Annotation) %>%
  summarise(
    nF = sum(Gender == "F"),
    nM = sum(Gender == "M"),
    p  = ifelse(nF >= 2 & nM >= 2,
                wilcox.test(pct_expanded_among_tcr ~ Gender)$p.value,
                NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

write.csv(tests_pct,
          file.path(out_dir, "tables", "Wilcoxon_GenderEffect_pctExpanded_StratifedByTimepoint.csv"),
          row.names = FALSE)

# ---------------------------- #
# 3) Plot: % Expanded among TCR+
# ---------------------------- #
p_pct <- ggplot(df_sum,
                aes(x = Gender, y = pct_expanded_among_tcr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  facet_grid(Manual_Annotation ~ Timepoint_Group, scales = "free_y") +
  ylab("% Expanded among TCR+ cells") +
  xlab("Gender") +
  ggtitle("Gender effect on clonal expansion stratified by Timepoint_Group")


ggsave(file.path(out_dir, "plots",
                 "GenderEffect_pctExpanded_StratifedByTimepoint.png"),
       p_pct, width = 14, height = 8, dpi = 300)

# ---------------------------- #
# 4) Plot: Mean clone size
# ---------------------------- #
p_mean <- ggplot(df_sum,
                 aes(x = Gender, y = mean_clone_size_tcr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  facet_grid(Manual_Annotation ~ Timepoint_Group, scales = "free_y") +
  ylab("Mean clonalFrequency among TCR+ cells") +
  xlab("Gender") +
  ggtitle("Gender effect on mean clone size stratified by Timepoint_Group")


ggsave(file.path(out_dir, "plots",
                 "GenderEffect_meanCloneSize_StratifedByTimepoint.png"),
       p_mean, width = 14, height = 8, dpi = 300)

message("\nDONE. Outputs saved to:\n", out_dir)
print(tests_pct)
##% Clonal Expansion (among TCR⁺ cells) -  How much of the population is participating in any expansion at all?
## Mean Clonal Frequency (clone size) - How large the expansions are when they occur

############################################################
# Gender effect on Viral Load
# Stratified by Timepoint_Group
############################################################


out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Gender_Effect/ViralLoad_By_Timepoint"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

# ---------------------------------------------------------
# 1) Extract unique sample-level viral load
# ---------------------------------------------------------
df_vl <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(orig.ident, Gender, Condition, Timepoint_Group, Viral_Load) %>%
  distinct(orig.ident, .keep_all = TRUE) %>%     # one row per sample
  filter(Condition == "HEI") %>%                 # viral load meaningful for HEI
  filter(!is.na(Viral_Load))

write.csv(df_vl,
          file.path(out_dir, "tables", "SampleLevel_ViralLoad_HEI.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# 2) Log10 transform (recommended)
# ---------------------------------------------------------
df_vl$Viral_Load <- as.numeric(df_vl$Viral_Load)
df_vl$log10_VL <- log10(df_vl$Viral_Load + 1)

# ---------------------------------------------------------
# 3) Wilcoxon test per Timepoint_Group
# ---------------------------------------------------------
tests_vl <- df_vl %>%
  group_by(Timepoint_Group) %>%
  summarise(
    nF = sum(Gender == "F"),
    nM = sum(Gender == "M"),
    p  = ifelse(nF >= 2 & nM >= 2,
                wilcox.test(log10_VL ~ Gender)$p.value,
                NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

write.csv(tests_vl,
          file.path(out_dir, "tables", "Wilcoxon_GenderEffect_ViralLoad_byTimepoint.csv"),
          row.names = FALSE)

# ---------------------------------------------------------
# 4) Plot Viral Load by Gender within Timepoint_Group
# ---------------------------------------------------------
p_vl <- ggplot(df_vl,
               aes(x = Gender, y = log10_VL)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  facet_wrap(~ Timepoint_Group, scales = "free_y") +
  ylab("log10 Viral Load") +
  xlab("Gender") +
  ggtitle("Gender effect on Viral Load (HEI only)")

ggsave(file.path(out_dir, "plots",
                 "GenderEffect_ViralLoad_HEI_byTimepoint.png"),
       p_vl, width = 10, height = 5, dpi = 300)

message("\nDONE. Viral load gender effect analysis saved to:\n", out_dir)
print(tests_vl)

##############
############################################################
# Longitudinal CD57+ CD28- CD8s (ADT-based)
# - Gates using ADT (not RNA)
# - Summarize per sample (orig.ident) across timepoints
# - Longitudinal plots vs Age (months) grouped by PID
# - Saves outputs
############################################################

# ---------------------------- #
# 0) Output directory
# ---------------------------- #
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CD57_CD28_Longitudinal"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# 1) Identify ADT feature names (robust)
# ---------------------------- #
DefaultAssay(TARA_ALL) <- "ADT"
adt_feats <- rownames(TARA_ALL[["ADT"]])

# Try to find CD57 and CD28 features (handles TotalSeq / weird naming)
cd57_candidates <- grep("CD57|B3GAT1", adt_feats, value = TRUE, ignore.case = TRUE)
cd28_candidates <- grep("CD28",        adt_feats, value = TRUE, ignore.case = TRUE)

message("CD57 candidates in ADT:\n", paste(cd57_candidates, collapse = ", "))
message("CD28 candidates in ADT:\n", paste(cd28_candidates, collapse = ", "))

# Pick the first match by default (edit if you want a specific one)
if (length(cd57_candidates) == 0) stop("No ADT feature matching CD57/B3GAT1 found.")
if (length(cd28_candidates) == 0) stop("No ADT feature matching CD28 found.")

feat_cd57 <- cd57_candidates[1]
feat_cd28 <- cd28_candidates[1]

message("\nUsing ADT features:\n  CD57 = ", feat_cd57, "\n  CD28 = ", feat_cd28, "\n")

# ---------------------------- #
# 2) Define CD8 compartment (by your Manual_Annotation)
#    Edit this list if you want a different CD8 definition
# ---------------------------- #
cd8_clusters <- c(
  "1: Memory CD8 T cell",
  "6: Naïve CD8 T cell",
  "8: CTL-like",
  "9: TRDV1+ CTL-like",
  "27: GZMK+ CD8 T cell"
)

# ---------------------------- #
# 3) Helper: bimodal cutoff finder (dependency-light)
#    Uses density peaks and the trough between the top 2 peaks.
#    Falls back to median + 1*MAD if it can't find a clean trough.
# ---------------------------- #
get_bimodal_cutoff <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 200) return(stats::median(x) + stats::mad(x))
  
  d <- stats::density(x, na.rm = TRUE)
  y <- d$y
  # local maxima
  peaks <- which(diff(sign(diff(y))) == -2) + 1
  if (length(peaks) < 2) return(stats::median(x) + stats::mad(x))
  
  # take two highest peaks
  peaks <- peaks[order(y[peaks], decreasing = TRUE)][1:2]
  peaks <- sort(peaks)
  
  # trough between peaks
  region <- peaks[1]:peaks[2]
  trough_idx <- region[which.min(y[region])]
  cutoff <- d$x[trough_idx]
  
  # sanity fallback if cutoff is extreme
  if (!is.finite(cutoff)) cutoff <- stats::median(x) + stats::mad(x)
  cutoff
}

# ---------------------------- #
# 4) Pull ADT expression + metadata (cell-level),
#    then gate CD57+ CD28- within CD8 clusters
# ---------------------------- #
meta <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    PID = sub("_.*$", "", orig.ident)  # prefix before "_"
  )

adt_df <- FetchData(
  TARA_ALL,
  vars = c(feat_cd57, feat_cd28),
  slot = "data"              # ADT normalized data
) %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta, by = "cell")

# Keep only CD8 compartment for gating summaries
adt_cd8 <- adt_df %>%
  filter(Manual_Annotation %in% cd8_clusters)

# Compute global cutoffs (you can change to per-sample cutoffs if desired)
cut_cd57 <- get_bimodal_cutoff(adt_cd8[[feat_cd57]])
cut_cd28 <- get_bimodal_cutoff(adt_cd8[[feat_cd28]])

message("Computed cutoffs (global, CD8 cells only):\n",
        "  ", feat_cd57, " cutoff = ", signif(cut_cd57, 4), "\n",
        "  ", feat_cd28, " cutoff = ", signif(cut_cd28, 4), "\n")

adt_cd8 <- adt_cd8 %>%
  mutate(
    CD57_pos = .data[[feat_cd57]] > cut_cd57,
    CD28_pos = .data[[feat_cd28]] > cut_cd28,
    CD57pos_CD28neg = CD57_pos & !CD28_pos
  )

# ---------------------------- #
# 5) Sample-level longitudinal summary (orig.ident is the unit)
# ---------------------------- #
df_long <- adt_cd8 %>%
  group_by(orig.ident, PID, Condition, Timepoint_Group, Gender) %>%
  summarise(
    Age_months = median(Age, na.rm = TRUE),  # Age in months
    n_cd8 = n(),
    n_gate = sum(CD57pos_CD28neg, na.rm = TRUE),
    freq_gate = n_gate / n_cd8,
    median_CD57 = median(.data[[feat_cd57]], na.rm = TRUE),
    median_CD28 = median(.data[[feat_cd28]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(PID, Age_months)

write.csv(df_long, file.path(out_dir, "tables", "CD57pos_CD28neg_CD8_longitudinal_bySample.csv"), row.names = FALSE)

# ---------------------------- #
# 6) Plots: longitudinal frequency vs Age (months)
# ---------------------------- #

# A) Spaghetti plot by PID, colored by Gender, faceted by Timepoint_Group
p_freq_tp <- ggplot(df_long, aes(x = Age_months, y = freq_gate, group = PID, color = Gender)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.6, alpha = 0.85) +
  facet_wrap(~ Timepoint_Group, scales = "free_x") +
  ylab("Frequency of CD57+ CD28- among CD8 (ADT gate)") +
  xlab("Age (months)") +
  ggtitle("Longitudinal CD57+CD28- CD8s (ADT) by Timepoint_Group") +
  theme_bw()


ggsave(file.path(out_dir, "plots", "Longitudinal_CD57pos_CD28neg_CD8_byTimepointGroup.png"),
       p_freq_tp, width = 12, height = 6, dpi = 300)

# B) Same, but facet by Condition (HEI/HEU/HUU) to visually sanity-check
p_freq_cond <- ggplot(df_long, aes(x = Age_months, y = freq_gate, group = PID, color = Gender)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.6, alpha = 0.85) +
  facet_wrap(~ Condition, scales = "free_x") +
  ylab("Frequency of CD57+ CD28- among CD8 (ADT gate)") +
  xlab("Age (months)") +
  ggtitle("Longitudinal CD57+CD28- CD8s (ADT) by Condition") +
  theme_bw()


ggsave(file.path(out_dir, "plots", "Longitudinal_CD57pos_CD28neg_CD8_byCondition.png"),
       p_freq_cond, width = 10, height = 5, dpi = 300)

# ---------------------------- #
# 7) Whole UMAP overlay (DimPlot2): highlight gated cells, everything else gray
# ---------------------------- #
# Add gate label back into Seurat metadata
gate_vec <- rep(FALSE, nrow(TARA_ALL@meta.data))
names(gate_vec) <- rownames(TARA_ALL@meta.data)
gate_vec[adt_cd8$cell] <- adt_cd8$CD57pos_CD28neg

TARA_ALL$CD57pos_CD28neg_CD8Gate <- ifelse(gate_vec, "CD57+CD28- (CD8)", "Other")
TARA_ALL$CD57pos_CD28neg_CD8Gate <- factor(TARA_ALL$CD57pos_CD28neg_CD8Gate,
                                           levels = c("Other","CD57+CD28- (CD8)"))

p_umap_gate <- DimPlot2(
  TARA_ALL,
  cols = "default",
  group.by = "CD57pos_CD28neg_CD8Gate",
  reduction = "wnn.umap"
) + ggtitle("Whole UMAP: CD57+CD28- CD8 gate (ADT) highlighted")


ggsave(file.path(out_dir, "plots", "UMAP_whole_CD57posCD28neg_CD8Gate.png"),
       p_umap_gate, width = 8, height = 6, dpi = 300)

###############
############################################################
############################################################
# Per-PID longitudinal UMAP panels (side-by-side) + per-PID
# longitudinal bar/line summaries for CD57+CD28- (CD8 gate)
#
# - Uses your existing metadata column: CD57pos_CD28neg_CD8Gate
# - Colors by CD57+CD28- status (easy + clear)
# - Orders timepoints by sample median Age (months)
# - Dynamically scales figure width by #timepoints
# - Saves PNG only
############################################################


library(patchwork)

# ---------------------------- #
# Output directory
# ---------------------------- #
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/CD57_CD28_Longitudinal"
plot_dir <- file.path(out_dir, "plots", "CD57posCD28neg_byPID")
tab_dir  <- file.path(out_dir, "tables", "CD57posCD28neg_byPID")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir,  recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# Sanity checks
# ---------------------------- #
stopifnot(all(c("orig.ident","Age","Manual_Annotation","CD57pos_CD28neg_CD8Gate") %in% colnames(TARA_ALL@meta.data)))

red_use <- if ("wnn.umap" %in% names(TARA_ALL@reductions)) "wnn.umap" else "umap"
message("Using reduction: ", red_use)

safe_name <- function(x) {
  x %>%
    gsub("[^A-Za-z0-9_\\-]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}

# ---------------------------- #
# CD8 compartment definition used for frequency summaries
# (edit if you want different CD8 definition)
# ---------------------------- #
cd8_clusters <- c(
  "1: Memory CD8 T cell",
  "6: Naïve CD8 T cell",
  "8: CTL-like",
  "9: TRDV1+ CTL-like",
  "27: GZMK+ CD8 T cell"
)

# ---------------------------- #
# 1) Sample-level Age ordering table
# ---------------------------- #
sample_info <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(orig.ident) %>%
  summarise(
    PID = sub("_.*$", "", unique(orig.ident)[1]),
    Age_median = median(Age, na.rm = TRUE),
    Gender = {
      g <- unique(na.omit(Gender))
      if (length(g) == 1) g else paste(g, collapse = ",")
    },
    Condition = {
      cnd <- unique(na.omit(Condition))
      if (length(cnd) == 1) cnd else paste(cnd, collapse = ",")
    },
    Timepoint_Group = {
      tg <- unique(na.omit(Timepoint_Group))
      if (length(tg) == 1) tg else paste(tg, collapse = ",")
    },
    .groups = "drop"
  )

write.csv(sample_info, file.path(tab_dir, "SampleInfo_byOrigIdent.csv"), row.names = FALSE)

# ---------------------------- #
# 2) Cell-level -> sample-level frequency table (CD8 compartment)
# ---------------------------- #
df_long <- TARA_ALL@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(PID = sub("_.*$", "", orig.ident)) %>%
  filter(Manual_Annotation %in% cd8_clusters) %>%
  group_by(orig.ident, PID) %>%
  summarise(
    Age_months = median(Age, na.rm = TRUE),
    Gender = {
      g <- unique(na.omit(Gender))
      if (length(g) == 1) g else paste(g, collapse = ",")
    },
    Condition = {
      cnd <- unique(na.omit(Condition))
      if (length(cnd) == 1) cnd else paste(cnd, collapse = ",")
    },
    Timepoint_Group = {
      tg <- unique(na.omit(Timepoint_Group))
      if (length(tg) == 1) tg else paste(tg, collapse = ",")
    },
    n_cd8 = n(),
    n_CD57posCD28neg = sum(CD57pos_CD28neg_CD8Gate == "CD57+CD28- (CD8)", na.rm = TRUE),
    freq_CD57posCD28neg = ifelse(n_cd8 > 0, n_CD57posCD28neg / n_cd8, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(PID, Age_months)

write.csv(df_long, file.path(tab_dir, "CD57posCD28neg_CD8_Frequency_bySample.csv"), row.names = FALSE)

# ---------------------------- #
# 3) Per-PID plotting function  (FIXED)
# ---------------------------- #
plot_pid_panels <- function(pid, obj, df_long, sample_info, red_use, plot_dir, tab_dir) {
  
  samp_pid <- sample_info %>%
    dplyr::filter(PID == pid) %>%
    dplyr::arrange(Age_median)
  
  origs <- samp_pid$orig.ident
  if (length(origs) == 0) return(invisible(NULL))
  
  # --- A) UMAP panels (side-by-side) colored by gate status
  umap_plots <- list()
  for (s in origs) {
    seu_s <- subset(obj, subset = orig.ident == s)
    
    age_s <- samp_pid$Age_median[samp_pid$orig.ident == s][1]
    title_stub <- paste0(pid, " | ", s, " | Age=", round(age_s, 2), " mo")
    
    p <- DimPlot2(
      seu_s,
      cols = "default",
      group.by = "CD57pos_CD28neg_CD8Gate",
      reduction = red_use
    ) +
      ggtitle(title_stub) +
      theme(
        plot.title = element_text(size = 10),
        legend.title = element_blank()
      )
    
    umap_plots[[s]] <- p
  }
  
  umap_combined <- patchwork::wrap_plots(umap_plots, nrow = 1) +
    patchwork::plot_annotation(
      title = paste0("CD57+CD28- longitudinal UMAPs — ", pid),
      subtitle = "Colored by CD57pos_CD28neg_CD8Gate"
    )
  
  umap_w <- max(8, 4 * length(origs))
  umap_h <- 5
  
  ggsave(
    filename = file.path(plot_dir, paste0(safe_name(pid), "__UMAPs_CD57posCD28neg.png")),
    plot = umap_combined,
    width = umap_w, height = umap_h, dpi = 300
  )
  
  # --- B) Longitudinal summary per PID
  df_pid <- df_long %>%
    dplyr::filter(PID == pid) %>%
    dplyr::arrange(Age_months) %>%
    dplyr::mutate(
      orig_label = paste0(orig.ident, "\n(", round(Age_months, 1), "m)"),
      freq_delta_from_first = freq_CD57posCD28neg - dplyr::first(freq_CD57posCD28neg)
    )
  
  p_bar <- ggplot(df_pid, aes(x = reorder(orig_label, Age_months), y = freq_CD57posCD28neg)) +
    geom_col() +
    labs(
      title = paste0("CD57+CD28- among CD8 over time — ", pid),
      subtitle = paste0(
        "Gender=", unique(df_pid$Gender), " | Condition=", unique(df_pid$Condition),
        ifelse(any(!is.na(df_pid$Timepoint_Group)),
               paste0(" | Timepoint_Group=", paste(unique(df_pid$Timepoint_Group), collapse = ",")),
               "")
      ),
      x = "Sample (Age months)",
      y = "Frequency (CD57+CD28- among CD8)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.8))
  
  p_line <- ggplot(df_pid, aes(x = Age_months, y = freq_CD57posCD28neg)) +
    geom_line() +
    geom_point(size = 2) +
    labs(
      title = paste0("Trajectory — ", pid),
      x = "Age (months)",
      y = "Frequency (CD57+CD28- among CD8)"
    ) +
    theme_bw()
  
  p_delta <- ggplot(df_pid, aes(x = Age_months, y = freq_delta_from_first)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line() +
    geom_point(size = 2) +
    labs(
      title = paste0("Change from first timepoint — ", pid),
      x = "Age (months)",
      y = "Δ Frequency vs first"
    ) +
    theme_bw()
  
  summary_combined <- (p_bar / p_line / p_delta) + patchwork::plot_layout(heights = c(1.2, 1, 1))
  
  sum_w <- max(8, 2.2 + 1.2 * length(origs))
  sum_h <- 12
  
  ggsave(
    filename = file.path(plot_dir, paste0(safe_name(pid), "__Summary_CD57posCD28neg.png")),
    plot = summary_combined,
    width = sum_w, height = sum_h, dpi = 300
  )
  
  # Save per-PID table
  write.csv(df_pid, file.path(tab_dir, paste0(safe_name(pid), "__CD57posCD28neg_bySample.csv")), row.names = FALSE)
  
  invisible(TRUE)
}

# ---------------------------- #
# 4) Run for all PIDs present (FIXED)
# ---------------------------- #
pids <- sort(unique(sample_info$PID))
for (pid in pids) {
  plot_pid_panels(
    pid         = pid,
    obj         = TARA_ALL,
    df_long     = df_long,
    sample_info = sample_info,
    red_use     = red_use,
    plot_dir    = plot_dir,
    tab_dir     = tab_dir
  )
  message("DONE PID: ", pid)
}
