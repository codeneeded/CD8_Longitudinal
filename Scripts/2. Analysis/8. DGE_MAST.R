# ---- Setup ----
library(Seurat)
library(dplyr)
library(ggplot2)


# Load annotated Seurat object
TARA_ALL <- readRDS("~/Documents/CD8_Longitudinal/saved_R_data/TARA_ALL_post_annotation.rds")

# Root output directory (everything goes here)
base_out <- "/home/akshay-iyer/Documents/CD8_Longitudinal/Differential_Gene_Expression"
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

# Ensure required metadata
DefaultAssay(TARA_ALL) <- "RNA"
Idents(TARA_ALL) <- "Manual_Annotation"
TARA_ALL$Comparison_Group <- TARA_ALL$Condition # HEI, HEU, HUU
TARA_ALL$PID <- sub("_.*", "", TARA_ALL$orig.ident, perl = TRUE)
TARA_ALL$Viral_Load <- suppressWarnings(as.numeric(as.character(TARA_ALL$Viral_Load)))

# Helpers
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
count_sig <- function(df, pcol = "p_val_adj", thr = 0.05) sum(df[[pcol]] < thr, na.rm = TRUE)

run_clusterwise_de <- function(obj, group_col, ident1, ident2, out_dir,
                               latent = c("nCount_RNA"),
                               min_cells_per_grp = 10,
                               title_prefix = "", file_stub = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  Idents(obj) <- "Manual_Annotation"
  DefaultAssay(obj) <- "RNA"
  
  de_list <- list()
  count_tbl <- data.frame()
  
  cl_levels <- levels(obj)
  for (cl in cl_levels) {
    message("Processing cluster: ", cl)
    cl_cells <- WhichCells(obj, idents = cl)
    sc <- subset(obj, cells = cl_cells)
    sc$.__grp <- sc[[group_col]]
    
    tab <- table(sc$.__grp)
    if (all(c(ident1, ident2) %in% names(tab)) && all(tab[c(ident1, ident2)] >= min_cells_per_grp)) {
      de <- FindMarkers(
        sc,
        ident.1 = ident1, ident.2 = ident2,
        group.by = ".__grp",
        test.use = "MAST",
        latent.vars = latent
      )
      comp_name <- if (is.null(file_stub)) paste0(ident1, "_vs_", ident2) else file_stub
      out_csv <- file.path(out_dir, paste0(safe_name(cl), "_", comp_name, ".csv"))
      write.csv(de, out_csv)
      de_list[[cl]] <- de
      count_tbl <- rbind(count_tbl, data.frame(Cluster = cl, DE_Genes = count_sig(de)))
    }
  }
  
  # Plot
  if (nrow(count_tbl) > 0) {
    p <- ggplot(count_tbl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      xlab("Cluster") + ylab("# of DE Genes (adj. p < 0.05)") +
      ggtitle(paste0(title_prefix, ": Clusters Ranked by DE Genes (", ident1, " vs ", ident2, ")"))
    ggsave(filename = file.path(out_dir, paste0("Cluster_DE_Gene_Counts_", if (is.null(file_stub)) paste0(ident1, "vs", ident2) else file_stub, ".png")),
           plot = p, width = 10, height = 6)
  }
  
  invisible(list(de = de_list, counts = count_tbl))
}

# =====================================================================
# 1) HEI vs HEU (ONLY Entry: Age <= 2) -> folder: HEIvsHEU_PreART
# No repeat measures at this restricted timepoint; do NOT control for PID.
# =====================================================================
hei_heu_pre_dir <- file.path(base_out, "HEIvsHEU_PreART")
TARA_pre <- subset(TARA_ALL, subset = Age <= 2 & Comparison_Group %in% c("HEI", "HEU"))

run_clusterwise_de(
  obj       = TARA_pre,
  group_col = "Comparison_Group",
  ident1    = "HEI",
  ident2    = "HEU",
  out_dir   = hei_heu_pre_dir,
  latent    = c("nCount_RNA"),
  title_prefix = "TARA (Pre-ART Entry ≤ 2m)"
)

# =====================================================================
# 2) HEU vs HUU (all are already at entry age) -> folder: HEUvsHUU
# No repeat measures; do NOT control for PID.
# =====================================================================
heu_huu_dir <- file.path(base_out, "HEUvsHUU")
TARA_heu_huu <- subset(TARA_ALL, subset = Comparison_Group %in% c("HEU", "HUU"))

run_clusterwise_de(
  obj       = TARA_heu_huu,
  group_col = "Comparison_Group",
  ident1    = "HEU",
  ident2    = "HUU",
  out_dir   = heu_huu_dir,
  latent    = c("nCount_RNA"),
  title_prefix = "TARA"
)

# =====================================================================
# 3) High vs Low VL at Entry ONLY (Age <= 2) -> folder: HighvsLowVL_PreART
# No repeat measures (entry only); do NOT control for PID.
# =====================================================================
highlow_pre_dir <- file.path(base_out, "HighvsLowVL_PreART")
TARA_entry <- subset(TARA_ALL, subset = Age %in% c(1, 2) & Condition == "HEI")  # VL is defined only for HEI
TARA_entry$Viral_Load_Category <- ifelse(TARA_entry$Viral_Load >= 100000, "High", "Low")
TARA_entry$Comparison_Group <- TARA_entry$Viral_Load_Category

run_clusterwise_de(
  obj       = TARA_entry,
  group_col = "Comparison_Group",
  ident1    = "High",
  ident2    = "Low",
  out_dir   = highlow_pre_dir,
  latent    = c("nCount_RNA"),
  title_prefix = "TARA HEI (Pre-ART Entry ≤ 2m)",
  file_stub = "HighvsLowVL_PreART"
)

# =====================================================================
# 4) High vs Low VL across ALL timepoints -> folder: HighvsLowVL_ALL
# Repeat measures exist; CONTROL for PID.
# =====================================================================
highlow_all_dir <- file.path(base_out, "HighvsLowVL_ALL")
TARA_hei_all <- subset(TARA_ALL, subset = Condition == "HEI")
TARA_hei_all$Viral_Load_Category <- ifelse(TARA_hei_all$Viral_Load >= 100000, "High", "Low")
TARA_hei_all$Comparison_Group <- TARA_hei_all$Viral_Load_Category

run_clusterwise_de(
  obj       = TARA_hei_all,
  group_col = "Comparison_Group",
  ident1    = "High",
  ident2    = "Low",
  out_dir   = highlow_all_dir,
  latent    = c("nCount_RNA", "PID"),
  title_prefix = "TARA HEI (All Timepoints)",
  file_stub = "HighvsLowVL_ALL"
)

# =====================================================================
# 5) Post-ART Suppressed vs Pre-ART Entry (HEI only)
#    -> folder: PostART_Suppressed_vs_PreART
# CONTROL for PID (longitudinal).
# =====================================================================
post_supp_dir <- file.path(base_out, "PostART_Suppressed_vs_PreART")
TARA_HEI <- subset(TARA_ALL, subset = Condition == "HEI")

TARA_HEI$Timepoint_Group <- ifelse(
  TARA_HEI$Age <= 2, "PreART_Entry",
  ifelse(TARA_HEI$Viral_Load < 200, "PostART_Suppressed",
         ifelse(!is.na(TARA_HEI$Viral_Load) & TARA_HEI$Viral_Load >= 200, "PostART_Unsuppressed", NA))
)

TARA_HEI_PP <- subset(TARA_HEI, subset = Timepoint_Group %in% c("PreART_Entry", "PostART_Suppressed"))
TARA_HEI_PP$Comparison_Group <- TARA_HEI_PP$Timepoint_Group

run_clusterwise_de(
  obj       = TARA_HEI_PP,
  group_col = "Comparison_Group",
  ident1    = "PostART_Suppressed",
  ident2    = "PreART_Entry",
  out_dir   = post_supp_dir,
  latent    = c("nCount_RNA", "PID"),
  title_prefix = "TARA HEI",
  file_stub = "PostART_Suppressed_vs_PreART"
)

# =====================================================================
# 6) Post-ART Unsuppressed vs Pre-ART Entry (HEI only)
#    -> folder: PostART_Unsuppressed_vs_PreART
# CONTROL for PID (longitudinal).
# =====================================================================
post_unsupp_dir <- file.path(base_out, "PostART_Unsuppressed_vs_PreART")
TARA_HEI_PU <- subset(TARA_HEI, subset = Timepoint_Group %in% c("PreART_Entry", "PostART_Unsuppressed"))
TARA_HEI_PU$Comparison_Group <- TARA_HEI_PU$Timepoint_Group

run_clusterwise_de(
  obj       = TARA_HEI_PU,
  group_col = "Comparison_Group",
  ident1    = "PostART_Unsuppressed",
  ident2    = "PreART_Entry",
  out_dir   = post_unsupp_dir,
  latent    = c("nCount_RNA", "PID"),
  title_prefix = "TARA HEI",
  file_stub = "PostART_Unsuppressed_vs_PreART"
)

