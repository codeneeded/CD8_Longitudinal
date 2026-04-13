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

# =========================
# Trex per-child outputs
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(officer)
  library(flextable)
})

# -------- Paths --------
out_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/VDJ/TCR/Trex"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
plot_dir <- file.path(out_dir, "PerChild_Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
docx_path <- file.path(out_dir, "Trex_CloneCounts_and_ViralLoad.docx")

# -------- Column names (edit if needed) --------
orig_ident_col <- "orig.ident"
epitope_col    <- "TRB_Epitope.species"   # includes values like "HIV" or "HIV;CMV"
clone_col      <- "clonotype_id"          # <-- change if your Trex clone ID column uses a different name
viral_load_col <- "Viral_Load"            # <-- change to your VL col name if different

# -------- Helper: normalize timepoints CP018_V24 -> CP018_24m --------
normalize_time <- function(x) {
  x %>%
    str_replace("_V(\\d+)", "_\\1m") %>%
    str_replace("_entry$", "_0m")
}

# -------- Helper: split PatientID / Timepoint from orig.ident --------
split_pid_tp <- function(df) {
  df %>%
    mutate(
      !!orig_ident_col := normalize_time(.data[[orig_ident_col]]),
      PatientID = str_extract(.data[[orig_ident_col]], "^[^_]+"),
      Timepoint = str_replace(.data[[orig_ident_col]], "^[^_]+_", "")
    )
}

# -------- Core summarizer for a single TRB set (one edit distance) --------
summarize_trb <- function(seu_obj, edit_dist_label, clone_col = "CTstrict") {
  md <- seu_obj@meta.data
  
  # --- keep only cells with a clonotype (i.e., T cells with TCR) ---
  if (!clone_col %in% names(md)) stop("Clone ID column not found: ", clone_col)
  md[[clone_col]] <- as.character(md[[clone_col]])
  md <- md |> dplyr::filter(!is.na(.data[[clone_col]]), .data[[clone_col]] != "")
  
  # (optional) if you have a celltype column and want only T cells, uncomment & adapt:
  # md <- md |> dplyr::filter(grepl("T", CellType, ignore.case = TRUE))
  
  # --- normalize orig.ident to PatientID + Timepoint ---
  md$orig.ident <- stringr::str_replace(md$orig.ident, "_V(\\d+)", "_\\1m")
  md$orig.ident <- stringr::str_replace(md$orig.ident, "_entry$", "_0m")
  md$PatientID  <- stringr::str_extract(md$orig.ident, "^[^_]+")
  md$Timepoint  <- stringr::str_replace(md$orig.ident, "^[^_]+_", "")
  md$EditDistance <- as.character(edit_dist_label)
  
  # mark HIV anywhere in epitope combo
  epitope_col <- "TRB_Epitope.species"
  if (!epitope_col %in% names(md)) md[[epitope_col]] <- NA_character_
  md$HIV_any <- grepl("HIV", md[[epitope_col]], ignore.case = TRUE)
  
  # use your metadata for singlets
  if (!"clonalFrequency" %in% names(md)) stop("Metadata column 'clonalFrequency' not found.")
  md$clonalFrequency <- suppressWarnings(as.numeric(md$clonalFrequency))
  md$clonalFrequency[is.na(md$clonalFrequency)] <- 0
  
  # A) all clones (unique clonotypes per PatientID x Timepoint)
  all_clones <- md |>
    dplyr::distinct(PatientID, Timepoint, !!rlang::sym(clone_col)) |>
    dplyr::count(PatientID, Timepoint, name = "CloneCount") |>
    dplyr::mutate(Category = "All clones", EditDistance = as.character(edit_dist_label))
  
  # B) HIV-specific clones (any HIV in combo)
  hiv_clones <- md |>
    dplyr::filter(HIV_any) |>
    dplyr::distinct(PatientID, Timepoint, !!rlang::sym(clone_col)) |>
    dplyr::count(PatientID, Timepoint, name = "CloneCount") |>
    dplyr::mutate(Category = "HIV-specific clones", EditDistance = as.character(edit_dist_label))
  
  # C) exclude singlets via clonalFrequency > 1
  nonsinglet_clones <- md |>
    dplyr::filter(clonalFrequency > 1) |>
    dplyr::distinct(PatientID, Timepoint, !!rlang::sym(clone_col)) |>
    dplyr::count(PatientID, Timepoint, name = "CloneCount") |>
    dplyr::mutate(Category = "Clones excl. singlets", EditDistance = as.character(edit_dist_label))
  
  # epitope combos for plotting (keep raw combos, HIV-containing only)
  epitope_plot_df <- md |>
    dplyr::mutate(TRB_Epitope.species = tidyr::replace_na(.data[[epitope_col]], "Unknown")) |>
    dplyr::filter(grepl("HIV", .data[[epitope_col]], ignore.case = TRUE)) |>
    dplyr::count(PatientID, Timepoint, .data[[epitope_col]], name = "Count") |>
    dplyr::mutate(EditDistance = as.character(edit_dist_label))
  
  # viral load snapshot if present
  vl_df <- NULL
  if ("Viral_Load" %in% names(md)) {
    vl_df <- md |>
      dplyr::distinct(PatientID, Timepoint, Viral_Load) |>
      dplyr::rename(ViralLoad = Viral_Load) |>
      dplyr::mutate(EditDistance = as.character(edit_dist_label))
  }
  
  list(counts = dplyr::bind_rows(all_clones, hiv_clones, nonsinglet_clones),
       plot_df = epitope_plot_df,
       vl = vl_df)
}


# -------- Run for TRB 0/1/2 --------
trb_list <- list(`0` = TARA_ALL_TRB_0, `1` = TARA_ALL_TRB_1, `2` = TARA_ALL_TRB_2)
summ_list <- purrr::imap(trb_list, ~ summarize_trb(.x, .y, clone_col = "CTstrict"))

counts_all <- bind_rows(map(summ_list, "counts"))
plot_all   <- bind_rows(map(summ_list, "plot_df"))
vl_all     <- bind_rows(keep(map(summ_list, "vl"), ~ !is.null(.x)))

# =========================
# Per-child plots (faceted by edit distance)
# =========================
# For each child, make a stacked bar (Timepoint on x, fill by epitope combo), facet=EditDistance
unique_children <- sort(unique(plot_all$PatientID))

for (cid in unique_children) {
  dfc <- plot_all %>% filter(PatientID == cid)
  if (nrow(dfc) == 0) next
  
  # Order timepoints numerically when possible (0m, 1m, 12m, 24m, ...)
  num_tp <- str_extract(dfc$Timepoint, "\\d+(?=m$)") %>% as.integer()
  tp_levels <- dfc %>%
    mutate(num = num_tp) %>%
    arrange(num) %>%
    pull(Timepoint) %>%
    unique()
  dfc$Timepoint <- factor(dfc$Timepoint, levels = tp_levels)
  
  p <- ggplot(dfc, aes(x = Timepoint, y = Count, fill = .data[[epitope_col]])) +
    geom_col() +
    facet_wrap(~ EditDistance, ncol = 1,
               labeller = as_labeller(c(`0` = "Edit distance 0",
                                        `1` = "Edit distance 1",
                                        `2` = "Edit distance 2"))) +
    labs(title = paste0("HIV-containing epitope combos by timepoint — ", cid),
         x = "Timepoint", y = "Clone count", fill = "Epitope combo") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  ggsave(filename = file.path(plot_dir, paste0("Trex_HIVcombos_", cid, ".png")),
         plot = p, width = 8, height = 7, dpi = 300, bg = "white")
}

# =========================
# Word tables (horizontal by timepoint)
# =========================


# Ensure types
counts_all <- counts_all %>% mutate(EditDistance = as.character(EditDistance))
if (exists("vl_all")) vl_all <- vl_all %>% mutate(EditDistance = as.character(EditDistance))
# --- keep only CP* patients ---
counts_all <- counts_all %>% filter(grepl("^CP", PatientID))
if (exists("vl_all")) {
  vl_all <- vl_all %>% filter(grepl("^CP", PatientID))
}

# Timepoint order like 0m, 1m, 12m, 24m
tp_levels_all <- counts_all %>%
  distinct(Timepoint) %>%
  mutate(num = suppressWarnings(as.integer(str_extract(Timepoint, "\\d+(?=m$)")))) %>%
  arrange(num) %>% pull(Timepoint)

doc <- read_docx()

add_table_section <- function(doc, heading, df_wide) {
  doc <- body_add_par(doc, value = heading, style = "heading 2")
  ft  <- flextable(df_wide) |> autofit()
  doc <- body_add_flextable(doc, ft)
  return(doc)
}

make_block <- function(doc, ed, cat_label, counts_all) {
  df <- counts_all %>%
    filter(EditDistance == ed, Category == cat_label) %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_levels_all)) %>%
    arrange(PatientID, Timepoint)
  
  wide <- df %>%
    select(PatientID, Timepoint, CloneCount) %>%
    tidyr::pivot_wider(names_from = Timepoint, values_from = CloneCount) %>%
    arrange(PatientID)
  
  doc <- add_table_section(doc, cat_label, wide)
  return(doc)
}

# Build sections in the exact order you want
for ( ed in c("0","1","2") ) {
  doc <- body_add_par(doc, value = paste("Edit distance", ed), style = "heading 1")
  
  doc <- make_block(doc, ed, "All clones", counts_all)
  doc <- make_block(doc, ed, "HIV-specific clones", counts_all)
  doc <- make_block(doc, ed, "Clones excl. singlets", counts_all)
}

# Single viral-load table at the end (not faceted)
if (exists("vl_all") && nrow(vl_all) > 0) {
  vl_final <- vl_all %>%
    group_by(PatientID, Timepoint) %>%
    summarize(ViralLoad = dplyr::first(na.omit(ViralLoad)), .groups = "drop") %>%
    mutate(Timepoint = factor(Timepoint, levels = tp_levels_all)) %>%
    arrange(PatientID, Timepoint)
  
  vl_wide <- vl_final %>%
    tidyr::pivot_wider(names_from = Timepoint, values_from = ViralLoad) %>%
    arrange(PatientID)
  
  doc <- body_add_par(doc, value = "Viral load", style = "heading 1")
  doc <- add_table_section(doc, "Viral load by timepoint", vl_wide)
}

print(doc, target = docx_path)
