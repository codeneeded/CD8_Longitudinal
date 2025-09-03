library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(forcats)
# -----------------------------
# Paths & object
# -----------------------------
# Path to your .RData file
tcr_file <- "/home/akshay-iyer/Documents/CD8_Longitudinal/saved_R_data/TARA_TCR_Combined.RData"

# Load objects
load(tcr_file)

##########

process_trb <- function(seu_obj, edit_distance_label) {
  md <- seu_obj@meta.data %>%
    filter(!is.na(orig.ident), !is.na(TRB_Epitope.species))
  
  # 1) HIV-only (split combos, keep only exact "HIV")
  hiv_counts <- md %>%
    filter(grepl("HIV", TRB_Epitope.species, ignore.case = TRUE)) %>%
    separate_rows(TRB_Epitope.species, sep = ";") %>%
    mutate(TRB_Epitope.species = str_trim(TRB_Epitope.species)) %>%
    filter(toupper(TRB_Epitope.species) == "HIV") %>%
    count(orig.ident, TRB_Epitope.species, name = "Count") %>%
    mutate(EditDistance = edit_distance_label, .before = 1)
  
  # 2) HIV-containing combos (do NOT split; keep raw categories)
  hiv_combo_counts <- md %>%
    filter(grepl("HIV", TRB_Epitope.species, ignore.case = TRUE)) %>%
    count(orig.ident, TRB_Epitope.species, name = "Count") %>%
    mutate(EditDistance = edit_distance_label, .before = 1)
  
  # 3) Restrict to patients with multiple timepoints (using prefix before first "_")
  hiv_combo_counts_multiTP <- hiv_combo_counts %>%
    mutate(PatientID = str_extract(orig.ident, "^[^_]+")) %>%
    group_by(PatientID) %>%
    filter(n_distinct(orig.ident) > 1) %>%
    ungroup()
  
  list(
    hiv_counts = hiv_counts,
    hiv_combo_counts = hiv_combo_counts,
    hiv_combo_counts_multiTP = hiv_combo_counts_multiTP
  )
}

# Run for TRB edit distances 0, 1, 2
trb_objs <- list(`0` = TARA_ALL_TRB_0, `1` = TARA_ALL_TRB_1, `2` = TARA_ALL_TRB_2)
out_list <- imap(trb_objs, ~ process_trb(.x, .y))

# Bind results across edit distances
hiv_counts_all <- bind_rows(map(out_list, "hiv_counts"))
hiv_combo_counts_all <- bind_rows(map(out_list, "hiv_combo_counts"))
hiv_combo_counts_multiTP_all <- bind_rows(map(out_list, "hiv_combo_counts_multiTP"))

# Inspect
head(hiv_counts_all)
head(hiv_combo_counts_all)
head(hiv_combo_counts_multiTP_all)

fix_timepoint <- function(df) {
  df %>%
    mutate(orig.ident = str_replace(orig.ident, "_V(\\d+)", "_\\1m"))
}

# Apply to all your outputs
hiv_counts_all <- fix_timepoint(hiv_counts_all)
hiv_combo_counts_all <- fix_timepoint(hiv_combo_counts_all)
hiv_combo_counts_multiTP_all <- fix_timepoint(hiv_combo_counts_multiTP_all)

################ Plotting ######################


# -Save CSV in the working directory you specified ---
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/VDJ/TCR/Trex"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

csv_path <- file.path(out_dir, "HIV_combo_counts_multiTP_by_origident_TRB0_1_2.csv")
readr::write_csv(hiv_combo_counts_multiTP_all, csv_path)

# --- 3) Plot: stacked counts per orig.ident, faceted by EditDistance ---
# Order orig.ident by overall total counts for readability
totals <- hiv_combo_counts_multiTP_all %>%
  group_by(orig.ident) %>%
  summarise(Total = sum(Count), .groups = "drop")

plot_df <- hiv_combo_counts_multiTP_all %>%
  left_join(totals, by = "orig.ident") %>%
  mutate(orig.ident = fct_reorder(orig.ident, Total))  # global order by total

plot_df <- hiv_combo_counts_multiTP_all %>%
  left_join(totals, by = "orig.ident") %>%
  mutate(
    orig.ident = fct_reorder(orig.ident, Total),
    EditDistance = factor(EditDistance,
                          levels = c("0","1","2"),
                          labels = c("Edit distance 0",
                                     "Edit distance 1",
                                     "Edit distance 2"))
  )

p <- ggplot(plot_df,
            aes(x = orig.ident, y = Count, fill = TRB_Epitope.species)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ EditDistance, scales = "free_x") +
  labs(
    title = "HIV-containing epitope combinations per sample",
    x = "Sample",
    y = "Count",
    fill = "Epitope combo"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

p
png_path <- file.path(out_dir, "HIV_combo_counts_multiTP_by_origident_TRB_facets.png")
ggsave(png_path, p, width = 12, height = 8, dpi = 300, bg = "white")

