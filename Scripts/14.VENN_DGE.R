#### VENN DIAGRAM DGE ############
library(tidyverse)
library(ggVennDiagram)
library(fs)
library(stringr)
library(eulerr)
library(ggplot2)
library(readr)
################## DGE VENN DIAGRAMS #################################


# ---- Config ----
padj_thr <- 0.05
lfc_thr  <- 0.58
out_root <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/DGE_Summaries"  # <-- change if you like

dir_create(out_root)

input_dirs <- list( HEIvsHEU_PreART = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/HEIvsHEU_PreART", 
                    HEUvsHUU = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/HEUvsHUU", 
                    HighvsLowVL_PreART = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/HighvsLowVL_PreART", 
                    HighvsLowVL_ALL = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/HighvsLowVL_ALL", 
                    PostART_Suppressed_vs_PreART = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/PostART_Suppressed_vs_PreART", 
                    PostART_Unsuppressed_vs_PreART = "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression/PostART_Unsuppressed_vs_PreART" )

# ---- Helpers ----

# Read one markers CSV and standardize columns: gene, avg_log2FC, padj
read_markers_csv <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # Try to find gene column; prefer 'gene', else take the first column if it's character
  if (!"gene" %in% names(df)) {
    # Heuristic: use the first column as gene if it's character-like
    first_char_col <- names(df)[map_lgl(df, ~is.character(.x) || is.factor(.x))]
    if (length(first_char_col) >= 1) {
      df <- df %>% rename(gene = !!first_char_col[1])
    } else {
      # Last resort: if Seurat rownames were saved into an 'X1' column
      if ("X1" %in% names(df)) df <- df %>% rename(gene = X1)
    }
  }
  
  # Standardize log2FC column
  if ("avg_log2FC" %in% names(df)) {
    # ok
  } else if ("avg_logFC" %in% names(df)) {
    df <- df %>% rename(avg_log2FC = avg_logFC)
  } else {
    # Try any column containing 'log2FC'
    cand <- names(df)[str_detect(names(df), regex("log2?fc", ignore_case = TRUE))]
    if (length(cand) >= 1) df <- df %>% rename(avg_log2FC = !!cand[1])
  }
  
  # Standardize padj column
  if ("p_val_adj" %in% names(df)) {
    df <- df %>% rename(padj = p_val_adj)
  } else if ("padj" %in% names(df)) {
    df <- df %>% rename(padj = padj)
  } else {
    cand <- names(df)[str_detect(names(df), regex("adj", ignore_case = TRUE)) &
                        str_detect(names(df), regex("p", ignore_case = TRUE))]
    if (length(cand) >= 1) df <- df %>% rename(padj = !!cand[1])
  }
  
  # Return only if essentials exist
  stopifnot("gene" %in% names(df), "avg_log2FC" %in% names(df), "padj" %in% names(df))
  df
}

# Read all CSVs in a contrast directory; return a named list of data frames by "cluster key"
# We use the file stem (basename without .csv) as the cluster label, so your filenames define cluster names.
read_contrast_dir <- function(dir_path) {
  files <- dir_ls(dir_path, regexp = "\\.csv$", type = "file")
  if (length(files) == 0) return(list())
  set_names(
    map(files, read_markers_csv),
    nm = files %>% path_file() %>% path_ext_remove()
  )
}

# Significant gene vector
sig_genes <- function(df, padj_thr, lfc_thr) {
  df %>%
    filter(!is.na(padj), padj < padj_thr, abs(avg_log2FC) >= lfc_thr) %>%
    arrange(padj, desc(abs(avg_log2FC))) %>%
    pull(gene) %>% unique()
}

# Up/Down splits
split_up_down <- function(df, padj_thr, lfc_thr) {
  up <- df %>%
    filter(!is.na(padj), padj < padj_thr, avg_log2FC >= lfc_thr) %>%
    arrange(padj, desc(avg_log2FC))
  down <- df %>%
    filter(!is.na(padj), padj < padj_thr, avg_log2FC <= -lfc_thr) %>%
    arrange(padj, avg_log2FC)
  list(up = up, down = down)
}
# ---- 1) Venn: Post-ART Unsuppressed vs Pre-ART  vs  Post-ART Suppressed vs Pre-ART ----
library(eulerr)
library(ggplotify)

unsup_dir <- input_dirs$PostART_Unsuppressed_vs_PreART
supp_dir  <- input_dirs$PostART_Suppressed_vs_PreART

# Read all files
unsup_files <- fs::dir_ls(unsup_dir, regexp = "\\.csv$", type = "file")
supp_files  <- fs::dir_ls(supp_dir,  regexp = "\\.csv$", type = "file")

# Map: file stem -> base cluster key (strip contrast suffix)
base_key_unsup <- unsup_files %>%
  path_file() %>% path_ext_remove() %>% 
  str_remove("_PostART_Unsuppressed_vs_PreART$")
base_key_supp  <- supp_files %>%
  path_file() %>% path_ext_remove() %>% 
  str_remove("_PostART_Suppressed_vs_PreART$")

# Build named lists of data frames keyed by base cluster
unsup_list <- purrr::map(unsup_files, read_markers_csv) %>% rlang::set_names(base_key_unsup)
supp_list  <- purrr::map(supp_files,  read_markers_csv) %>% rlang::set_names(base_key_supp)

# Intersection of cluster keys present in BOTH contrasts
clusters <- intersect(names(unsup_list), names(supp_list))

venn_root <- fs::path(out_root, "Venn_PostART_vs_PreART")
fs::dir_create(fs::path(venn_root, "_summary"))

venn_counts <- tibble()         # summary counts (sig-based, for diagram)
all_clusters_table <- tibble()  # big tidy table (presence-based membership)

for (cl in clusters) {
  df_u <- unsup_list[[cl]]
  df_s <- supp_list[[cl]]
  
  # ---------- Diagram sets (SIG genes only) ----------
  genes_u_sig <- sig_genes(df_u, padj_thr, lfc_thr)
  genes_s_sig <- sig_genes(df_s, padj_thr, lfc_thr)
  
  venn.list <- list(
    "Post-ART Unsuppressed" = genes_u_sig,
    "Post-ART Suppressed"   = genes_s_sig
  )
  
  # counts for summary (diagram logic)
  overlap_sig <- intersect(genes_u_sig, genes_s_sig)
  u_only_sig  <- setdiff(genes_u_sig, genes_s_sig)
  s_only_sig  <- setdiff(genes_s_sig, genes_u_sig)
  
  venn_counts <- bind_rows(
    venn_counts,
    tibble(
      cluster           = cl,
      n_unsuppressed    = length(genes_u_sig),
      n_suppressed      = length(genes_s_sig),
      n_overlap         = length(overlap_sig),
      n_unsup_only      = length(u_only_sig),
      n_supp_only       = length(s_only_sig)
    )
  )
  
  # ---------- Plot with eulerr ----------
  fit <- euler(venn.list, shape = "circle")
  # Pretty cluster name & full title with both thresholds
  pretty_cl <- gsub("_", " ", cl)  # "0_CD4_T_cell" -> "0 CD4 T cell"
  title_txt <- sprintf("Cluster: %s  (padj < %s, |log2FC| \u2265 %s)",
                       pretty_cl,
                       format(padj_thr, trim = TRUE, scientific = FALSE),
                       format(lfc_thr, trim = TRUE, scientific = FALSE))
  
  gp <- plot(
    fit,
    fills = list(
      fill  = c("Post-ART Unsuppressed" = "#54928D",
                "Post-ART Suppressed"   = "#941C50"),
      alpha = 0.90
    ),
    edges = list(col = "#0F172A", lwd = 1),
    labels = FALSE,                                               # legend only
    quantities = list(col = "black", font = 2, cex = 1.15),       # region counts
    legend = TRUE
  )
  
  p <- as_ggplot(gp) +
    ggtitle(title_txt) +
    theme_void(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.margin     = margin(14, 22, 14, 22),
      legend.position = "bottom"
    )
  
  
  # Per-cluster folder
  cl_dir <- fs::path(venn_root, cl)
  fs::dir_create(cl_dir)
  
  ggsave(fs::path(cl_dir, paste0("Venn_", cl, "_eulerr.png")),
         p, width = 9, height = 5, dpi = 300)
  
  # ---------- Single tidy CSV per cluster (ALL genes; membership by PRESENCE) ----------
  comb <- full_join(
    df_u %>% select(gene, avg_log2FC_unsup = avg_log2FC, padj_unsup = padj),
    df_s %>% select(gene, avg_log2FC_supp = avg_log2FC, padj_supp  = padj),
    by = "gene"
  ) %>%
    mutate(
      in_unsup = !is.na(padj_unsup) | !is.na(avg_log2FC_unsup),
      in_supp  = !is.na(padj_supp)  | !is.na(avg_log2FC_supp),
      set_membership = dplyr::case_when(
        in_unsup & in_supp ~ "Overlap",
        in_unsup & !in_supp ~ "Unsuppressed_only",
        !in_unsup & in_supp ~ "Suppressed_only",
        TRUE ~ NA_character_
      ),
      # keep optional significance flags (NOT used for membership)
      sig_unsup = !is.na(padj_unsup) & padj_unsup < padj_thr &
        !is.na(avg_log2FC_unsup) & abs(avg_log2FC_unsup) >= lfc_thr,
      sig_supp  = !is.na(padj_supp)  & padj_supp  < padj_thr &
        !is.na(avg_log2FC_supp)  & abs(avg_log2FC_supp)  >= lfc_thr
    ) %>%
    select(gene, set_membership, padj_unsup, padj_supp, avg_log2FC_unsup, avg_log2FC_supp,
           sig_unsup, sig_supp)
  
  # write per-cluster CSV
  readr::write_csv(comb, fs::path(cl_dir, paste0("Venn_Table_", cl, ".csv")))
  
  # also collect into one big table with a cluster column at the front
  all_clusters_table <- bind_rows(
    all_clusters_table,
    comb %>% mutate(cluster = cl) %>% relocate(cluster, .before = 1)
  )
  
}

# ---------- Write summaries (robust) ----------
summary_dir <- fs::path(venn_root, "_summary")

# If a *file* named "_summary" exists, fail with a clear message
if (fs::file_exists(summary_dir) && !fs::is_dir(summary_dir)) {
  stop("A file named '_summary' exists at: ", summary_dir,
       "\nPlease remove/rename it so I can create the summary folder.")
}

# Create the directory (recursively) if missing
fs::dir_create(summary_dir, recurse = TRUE)

# (Optional) sanity check + diagnostics if something goes wrong
stopifnot(fs::dir_exists(summary_dir))

# Now write safely
readr::write_csv(venn_counts,      fs::path(summary_dir, "Venn_Counts.csv"))
readr::write_csv(all_clusters_table, fs::path(summary_dir, "All_Genes_Table.csv"))

# ---- 2) Up/Down counts for High VL vs Low VL (PreART and ALL) + gene tables ----

hl_sets <- c("HighvsLowVL_PreART", "HighvsLowVL_ALL")
for (key in hl_sets) {
  hl_dir <- input_dirs[[key]]
  if (is.null(hl_dir)) next
  hl_list <- read_contrast_dir(hl_dir)
  
  if (length(hl_list) == 0) next
  
  out_dir <- path(out_root, paste0("HighVL_vs_LowVL_", key))
  dir_create(c(out_dir, path(out_dir, "tables")))
  
  counts <- tibble()
  genes_long <- tibble()
  
  for (cl in names(hl_list)) {
    df <- hl_list[[cl]]
    splits <- split_up_down(df, padj_thr, lfc_thr)
    up   <- splits$up
    down <- splits$down
    
    counts <- bind_rows(counts,
                        tibble(
                          cluster = cl,
                          n_up   = nrow(up),
                          n_down = nrow(down),
                          n_sig  = nrow(up) + nrow(down)
                        )
    )
    
    # Gene-level table (keep stats)
    if (nrow(up) > 0) {
      genes_long <- bind_rows(
        genes_long,
        up %>% transmute(cluster = cl, direction = "Up", gene, avg_log2FC, padj)
      )
    }
    if (nrow(down) > 0) {
      genes_long <- bind_rows(
        genes_long,
        down %>% transmute(cluster = cl, direction = "Down", gene, avg_log2FC, padj)
      )
    }
  }
  
  write_csv(counts,     path(out_dir, "tables", paste0("HighLow_Counts_", key, ".csv")))
  write_csv(genes_long, path(out_dir, "tables", paste0("HighLow_GeneLists_", key, ".csv")))
}

message("Done. Outputs written to: ", out_root)
