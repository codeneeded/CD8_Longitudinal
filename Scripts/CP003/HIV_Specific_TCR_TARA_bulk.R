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
library(ggpubr)
library(rstatix)
library(SeuratExtend)
# ----------------------------- #
# 0) Paths
# ----------------------------- #
base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
out_dir   <- file.path(base_dir, "CP003/Shared_clone_Tables_Bulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qs_cp003_in  <- file.path(saved_dir, "seu_CP003_HIVSpecificTCR_annotated.qs2")
rds_tara_in  <- file.path(saved_dir, "tara_cdnk_annotated.rds")
rds_cp003_out <- file.path(saved_dir, "cp003_with_cc_epitope.qs2")
rds_tara_out  <- file.path(saved_dir, "tara_cdnk_with_cp003_matches.qs2")

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
qs_save(cp003,file = rds_cp003_out)

############################
#Identification of Tscm
#############################
# --- Settings ---
plot_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Tscm"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

adt_tscm <- c("CD45RA", "IL7R", "CD27", "FAS", "SELL")

# ----------------------------- #
# Feature plots — ADT
# ----------------------------- #
DefaultAssay(tara_cdnk) <- 'ADT'
p_adt <- FeaturePlot(
  tara_cdnk,
  features  = adt_tscm,
  reduction = "wnn.umap",
  order     = TRUE,
  ncol      = 3
) & theme(legend.position = "right")

ggsave(
  file.path(plot_dir, "Tscm_ADT_FeaturePlot.png"),
  p_adt,
  width  = 15,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

# ----------------------------- #
# Define Tscm: CD45RA+ FAS+ CD27+ SELL+ IL7R+
# ADT is CLR-normalised — check distributions first
# ----------------------------- #
adt_expr <- GetAssayData(tara_cdnk, assay = "ADT", layer = "data")

# Check distributions to set thresholds
par(mfrow = c(2, 3))
for (m in adt_tscm) {
  hist(adt_expr[m, ], breaks = 50, main = m, xlab = "CLR expression")
}

# ----------------------------- #
# Tscm definition
# CD45RA+ FAS+  -- critical (separates from naive which is CD45RA+ FAS-)
# CD27+         -- high on Tscm
# SELL+         -- CD62L, lymph node homing
# IL7R+         -- survival signal
# ----------------------------- #
cd45ra_pos <- adt_expr["CD45RA", ] > 0.5
fas_pos    <- adt_expr["FAS",    ] > 0.5
cd27_pos   <- adt_expr["CD27",   ] > 0.5
sell_pos   <- adt_expr["SELL",   ] > 0.5
il7r_pos   <- adt_expr["IL7R",   ] > 0.5

is_tscm <- cd45ra_pos & fas_pos & cd27_pos & sell_pos & il7r_pos

cat("Cells meeting Tscm criteria:", sum(is_tscm), "\n")
cat("As % of total:", round(mean(is_tscm) * 100, 2), "%\n")

# ----------------------------- #
# Add to metadata
# ----------------------------- #
tara_cdnk$Tscm_Candidate <- ifelse(is_tscm, "Tscm", "Other")
tara_cdnk$Tscm_Candidate <- factor(tara_cdnk$Tscm_Candidate, levels = c("Other", "Tscm"))

# ----------------------------- #
# UMAP
# ----------------------------- #
p_tscm_umap <- DimPlot2(
  tara_cdnk,
  features  = "Tscm_Candidate",
  reduction = "wnn.umap",
  cols      = c("Other" = "grey85", "Tscm" = "#2196F3"),
  theme     = NoAxes()
)
p_tscm_umap
ggsave(
  file.path(plot_dir, "Tscm_Candidates_UMAP.png"),
  p_tscm_umap,
  width  = 8,
  height = 7,
  dpi    = 300,
  bg     = "white"
)

# ----------------------------- #
# Violin plots to validate
# ----------------------------- #
p_vln_adt <- VlnPlot2(
  tara_cdnk,
  features = adt_tscm,
  group.by = "Tscm_Candidate",
  assay    = "ADT",
  pt.size  = 0,
  show.mean = T,
  ncol     = 3
)


ggsave(
  file.path(plot_dir, "Tscm_ADT_VlnPlot.png"),
  p_vln_adt,
  width  = 15,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

# ----------------------------- #
# Summary table
# ----------------------------- #
tscm_summary <- tara_cdnk@meta.data %>%
  group_by(Tscm_Candidate, T_Cell_A) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(Tscm_Candidate), desc(n_cells))

tscm_summary
write.csv(
  tscm_summary,
  file = file.path(plot_dir, "Tscm_by_Cluster.csv"),
  row.names = FALSE
)

# Check FAS distribution specifically within CD45RA+ cells only
cd45ra_cells <- colnames(tara_cdnk)[adt_expr["CD45RA", ] > 0.5]

hist(
  adt_expr["FAS", cd45ra_cells],
  breaks = 100,
  main   = "FAS expression within CD45RA+ cells",
  xlab   = "CLR expression"
)

tscm_summary_2 <- tara_cdnk@meta.data %>%
  group_by(Tscm_Candidate, Timepoint_Group) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(Tscm_Candidate), desc(n_cells))

qs_save(tara_cdnk, file = rds_tara_out)


# ----------------------------- #
# Tscm proportion by timepoint group
# ----------------------------- #

library(ggpubr)
library(rstatix)

# ----------------------------- #
# Tscm proportion by timepoint group
# ----------------------------- #
tara_cdnk$PID <- sub("_.*", "", tara_cdnk$orig.ident, perl = TRUE)

# Calculate proportion of Tscm per sample per timepoint group
tscm_prop <- tara_cdnk@meta.data %>%
  group_by(PID, Timepoint_Group) %>%
  summarise(
    total_cells = n(),
    tscm_cells  = sum(Tscm_Candidate == "Tscm", na.rm = TRUE),
    pct_Tscm    = (tscm_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  filter(!is.na(Timepoint_Group)) %>%
  mutate(Timepoint_Group = factor(
    Timepoint_Group,
    levels = c("HUU", "HEU", "PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed")
  ))

# ----------------------------- #
# Remove outliers per group
# ----------------------------- #
tscm_prop_clean <- tscm_prop %>%
  group_by(Timepoint_Group) %>%
  filter(
    pct_Tscm >= quantile(pct_Tscm, 0.25) - 1.5 * IQR(pct_Tscm) &
      pct_Tscm <= quantile(pct_Tscm, 0.75) + 1.5 * IQR(pct_Tscm)
  ) %>%
  ungroup()

# ----------------------------- #
# Check pairing structure
# ----------------------------- #
cat("Samples per timepoint group:\n")
print(tscm_prop_clean %>% count(Timepoint_Group))


cat("\nSample to timepoint mapping:\n")
print(
  tara_cdnk@meta.data %>%
    select(PID, Timepoint_Group) %>%
    distinct() %>%
    arrange(PID)
)
# ----------------------------- #
# Significance testing
# All unpaired — sample sizes too small for reliable paired test
# ----------------------------- #
kw_test <- kruskal_test(pct_Tscm ~ Timepoint_Group, data = tscm_prop_clean)
cat("Kruskal-Wallis p =", kw_test$p, "\n")

pwc <- tscm_prop_clean %>%
  wilcox_test(pct_Tscm ~ Timepoint_Group, paired = FALSE, p.adjust.method = "none") %>%
  add_significance("p") %>%
  add_xy_position(x = "Timepoint_Group")

cat("\nPairwise comparisons:\n")
print(pwc %>% select(group1, group2, p, p.signif))

# ----------------------------- #
# Plot
# ----------------------------- #
p_tscm_prop <- ggplot(tscm_prop_clean,
                      aes(x = Timepoint_Group, y = pct_Tscm, fill = Timepoint_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    pwc %>% filter(p < 0.05),
    label      = "p.signif",
    tip.length = 0.01
  ) +
  scale_fill_manual(values = c(
    "HUU"                  = "#AEC6CF",
    "HEU"                  = "#FFD700",
    "PreART_Entry"         = "#FF8C69",
    "PostART_Suppressed"   = "#90EE90",
    "PostART_Unsuppressed" = "#E63946"
  )) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey30"),
    plot.caption    = element_text(hjust = 0.5, size = 9,  color = "grey50")
  ) +
  labs(
    x        = NULL,
    y        = "% Tscm cells",
    title    = "Tscm Proportion by Timepoint Group",
    subtitle = paste("Kruskal-Wallis p =", signif(kw_test$p, 3)),
    caption  = "Unpaired Wilcoxon, nominal p-values (unadjusted). Outliers removed per group (1.5x IQR)."
  )

# ----------------------------- #
# Save
# ----------------------------- #
ggsave(
  file.path(plot_dir, "Tscm_Proportion_by_TimepointGroup.png"),
  p_tscm_prop,
  width  = 10,
  height = 8,
  dpi    = 300,
  bg     = "white"
)

write.csv(
  tscm_prop_clean,
  file      = file.path(plot_dir, "Tscm_Proportion_by_TimepointGroup.csv"),
  row.names = FALSE
)

write.csv(
  pwc %>% mutate(across(where(is.list), ~ sapply(., paste, collapse = ";"))),
  file      = file.path(plot_dir, "Tscm_Proportion_Pairwise_Stats.csv"),
  row.names = FALSE
)
