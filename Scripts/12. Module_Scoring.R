suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(readr)
  # library(lme4)   # <- optional if you decide to use mixed models
})

# -----------------------------
# Paths & object
# -----------------------------
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "TARA_ALL_post_annotation.rds")
base_out   <- file.path(base_dir, "Module_Score")
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

TARA_ALL <- readRDS(rds_in)
DefaultAssay(TARA_ALL) <- 'RNA'
assay_use   <- DefaultAssay(TARA_ALL)

# Edit if your cluster/group columns differ
cluster_col <- "Manual_Annotation" 
group_col   <- "Comparison_Group"
pid_col     <- "PID"
# Add Comparison_Group based on Condition
Idents(TARA_ALL) <- "Manual_Annotation"
TARA_ALL$Comparison_Group <- TARA_ALL$Condition # HEI, HEU, HUU
TARA_ALL$PID <- sub("_.*", "", TARA_ALL$orig.ident, perl = TRUE)
TARA_ALL$Viral_Load <- suppressWarnings(as.numeric(as.character(TARA_ALL$Viral_Load)))
TARA_ALL$PID
# Optional relabel on plots only (does not change data)
label_map <- c(HEI = "pHIV", HEU = "pHEU", HUU = "pHUU",
               PreART_Entry = "Pre-ART Entry",
               PostART_Suppressed = "Post-ART Suppressed",
               PostART_Unsuppressed = "Post-ART Unsuppressed",
               High = "High VL", Low = "Low VL")

relabel <- function(x) ifelse(x %in% names(label_map), label_map[x], x)

# -----------------------------
# Module definitions
# -----------------------------
modules <- list(
  # ---- TGF-beta ----
  TGFb = list(
    ReceptorsLigands = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3","ENG","ACVRL1"),
    SMAD_Core        = c("SMAD2","SMAD3","SMAD4","SMAD7","SKI","SKIL","PMEPA1"),
    Downstream       = c("FOXP3","IKZF2","CTLA4","LRRC32","ITGAE","CCR7","SELL","SERPINE1","SERPINE1","PTPN14","COL1A1","COL3A1"),
    Noncanonical     = c("MAPK14","MAPK8","MAPK9","MAP3K7","FOS","JUN","PIK3CD","AKT1","MTOR","RPTOR","NFKBIA","RELB")
  ),
  # ---- IL-10 / STAT3 ----
  IL10_STAT3 = list(
    ReceptorsTransducers = c("IL10RA","IL10RB","STAT3"),
    TFH_Adjacent_TFs     = c("BCL6","MAF","BATF","PIM1"),
    Treg_Overlap         = c("FOXP3","CTLA4","IKZF2"),
    Myeloid_M2_Hue       = c("MRC1","MSR1","CD163")
  ),
  # ---- Reservoir CD4 ----
  Reservoir_CD4 = list(
    Tfh_Core          = c("CXCR5","PDCD1","ICOS","BCL6","MAF","BATF","SH2D1A","IL21"),
    Tph_Discriminator = c("PRDM1","CXCL13","CCR2","PDCD1","ICOS"),
    Reservoir_Memory  = c("TCF7","CCR7","SELL","BCL2","IL7R")
  ),
  # ---- ISG ----
  ISG = list(
    TypeI  = c("ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","MX2","OAS1","OAS2","OAS3","RSAD2",
               "IFITM1","IFITM2","IFITM3","BST2","TRIM5","APOBEC3G","SAMHD1"),
    TypeII = c("STAT1","IRF1","GBP1","GBP5","CXCL9","CXCL10")
  ),
  # ---- Myeloid ----
  Myeloid = list(
    Inflammatory  = c("IL1B","TNF","NFKBIA","NFKBIZ","CCL2","CXCL8","SERPINA1","S100A8","S100A9"),
    Regulatory    = c("IL10RA","IL10RB","STAT3","SOCS3","TGFB1","TGFBR1","TGFBR2","MRC1","MSR1","CD163"),
    DC_Maturation = c("CCR7","LAMP3","FSCN1","CCL19","CD40","RELB")
  ),
  # ---- NK ----
  NK = list(
    Inhibition  = c("KLRC1","KLRD1"),
    Effector    = c("PRF1","GNLY","NKG7","GZMK"),
    Activation  = c("IFNG","XCL1","XCL2","CCL5")
  ),
  # ---- CD8 ----
  CD8 = list(
    Effector    = c("PRF1","GNLY","NKG7","CCL5","XCL1","XCL2","IFNG"),
    Memory      = c("TCF7","CCR7","SELL","BCL2","IL7R"),
    Exhaustion  = c("PDCD1","TIGIT","LAG3","HAVCR2","TOX","EOMES")
  ),
  # ---- B cells ----
  Bcell = list(
    GC_like      = c("BCL6","CXCR5","AICDA","S1PR2","RGS13","MEF2B","BACH2"),
    Atypical_ABC = c("TBX21","ITGAX","FCRL5","ZEB2"),
    Plasmablast  = c("PRDM1","XBP1","JCHAIN","MZB1","SDC1","TNFRSF17")
  ),
  # ---- Autophagy ----
  Autophagy = list(
    Initiation   = c("ULK1","RB1CC1","ATG13","ATG101","BECN1","PIK3C3","PIK3R4","ATG14","WIPI1","WIPI2","ATG2A","ATG2B"),
    Conjugation  = c("ATG5","ATG12","ATG16L1","ATG7","ATG3","MAP1LC3A","MAP1LC3B","MAP1LC3B2","GABARAP","GABARAPL1","GABARAPL2"),
    Fusion       = c("RAB7A","STX17","SNAP29","VAMP8","EPG5"),
    Lysosome     = c("LAMP1","LAMP2","CTSD","CTSB","ATP6V1A","ATP6V1B2","ATP6V0D1","PSAP","MCOLN1"),
    CLEAR_TFs    = c("TFEB","TFE3","MITF"),
    Cargo        = c("SQSTM1","NBR1","OPTN","CALCOCO2","TAX1BP1","TOLLIP"),
    Mitophagy    = c("PINK1","PRKN","BNIP3","BNIP3L","FUNDC1","PHB2","FKBP8","BCL2L13"),
    Xenophagy    = c("NOD2","RIPK2","TBK1","OPTN","CALCOCO2","SQSTM1","TAX1BP1","LRSAM1","ATG16L1","ATG5","ATG7"),
    cGAS_STING   = c("STING1","CGAS"),
    mTORC        = c("MTOR","RPTOR","MLST8","RHEB","TSC1","TSC2","DEPDC5"),
    AMPK         = c("PRKAA1","PRKAA2","PRKAB1","PRKAG1","STK11")
  ),
  # ---- New families you requested ----
  Cytokines = list(
    Core = c("IL6","CXCL8","IL10","TNF","IL18","IFNA1","IFNG","CXCL10","CXCL11",
             "IL1B","IL12A","IL12B","IL17A","IL23A","IL4","IL5","IL13",
             "CCL2","CCL3","CCL4","CCL5","CCL7","CCL20","CXCL9","CXCL12")
  ),
  ImmuneActivation = list(
    Core = c("TNFRSF1A","TNFRSF1B","IL2RA","CRP","GCH1","FGA","FGB","FGG",
             "SAA1","SAA2","HP","SERPINA1","IL1RN","IDO1","CD274","PDCD1LG2",
             "CXCL13","B2M","LAG3","HAVCR2")
  ),
  AdhesionMolecules = list(
    Core = c("VCAM1","ICAM1","SELL","SELE","ITGAL","ITGAM","ITGB2","PECAM1","CDH5","SPN","CD44","JAM3")
  ),
  GutTranslocation = list(
    Core = c("LBP","CD14","CD163","CCL2","REG3A","DEFA5","DEFB1","TJP1","OCLN","CLDN2","CLDN3","MUC2",
             "S100A8","S100A9","TLR2","TLR4","NOD2","MYD88")
  ),
  SenescenceMarkers = list(
    Core = c("CDKN2A","CDKN1A","TP53","GLB1","SERPINE1","IL6","MMP3","MMP9","CXCL8","CCL2")
  ),
  ExhaustionMarkers = list(
    Core = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","ENTPD1","TOX","EOMES","BATF")
  )
)

# -----------------------------
# Utilities
# -----------------------------
score_module <- function(obj, genes, name) {
  present <- intersect(genes, rownames(obj))
  if (length(present) == 0) {
    obj[[paste0("MS_", name)]] <- NA_real_
    return(list(obj=obj, present=character(0), missing=genes))
  }
  obj <- AddModuleScore(obj, features=list(present), name=paste0("MS_", name), assay=assay_use)
  col <- paste0("MS_", name, "1")
  obj[[paste0("MS_", name)]] <- obj[[col]]
  obj[[col]] <- NULL
  list(obj=obj, present=present, missing=setdiff(genes, present))
}

# Plot per-cell (unpaired) OR per-PID (paired) with optional paired lines
plot_module <- function(df, module_field, folder, title, paired_by=NULL) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(paired_by)) {
    # summarize per PID x cluster x group
    sumdf <- df %>%
      group_by(.data[[pid_col]], .data[[cluster_col]], .data[[group_col]]) %>%
      summarize(Score = mean(.data[[module_field]], na.rm=TRUE), .groups="drop")
    # keep only PIDs present in both groups within each facet
    keep <- sumdf %>%
      distinct(.data[[pid_col]], .data[[cluster_col]], .data[[group_col]]) %>%
      add_count(.data[[pid_col]], .data[[cluster_col]]) %>%
      filter(n >= 2) %>%    # has both groups
      distinct(.data[[pid_col]], .data[[cluster_col]])
    sumdf <- sumdf %>%
      inner_join(keep, by=c(pid_col, cluster_col))
    sumdf$GroupLabel <- relabel(sumdf[[group_col]])
    ","
    gg <- ggplot(sumdf, aes(x=GroupLabel, y=Score)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.9) +
      geom_point(aes(group=.data[[pid_col]]), position=position_jitter(width=0.15), size=0.7, alpha=0.5) +
      geom_line(aes(group=.data[[pid_col]]), alpha=0.2) +
      stat_compare_means(method="wilcox.test", label="p.format",
                         method.args=list(paired=TRUE)) +
      facet_wrap(as.formula(paste0("~", cluster_col)), scales="free_y") +
      theme_minimal(base_size=14) + xlab("") + ylab("Module Score") +
      ggtitle(title)
    ggsave(file.path(folder, paste0(module_field, ".png")), gg,
           width=20, height=18, dpi=300, bg="white")
  } else {
    df$GroupLabel <- relabel(df[[group_col]])
    gg <- ggplot(df, aes(x=GroupLabel, y=.data[[module_field]], fill=GroupLabel)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.85) +
      geom_jitter(width = 0.2, size = 0.4, alpha = 0.3) +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      facet_wrap(as.formula(paste0("~", cluster_col)), scales = "free_y") +
      theme_minimal(base_size = 14) +
      xlab("") + ylab("Module Score") + ggtitle(title)
    ggsave(file.path(folder, paste0(module_field, ".png")), gg,
           width=20, height=18, dpi=300, bg="white")
  }
}

# Run scoring for a given object+comparison name
run_scoring <- function(obj, out_dir, comparison_name, paired_by=NULL) {
  message("Running: ", comparison_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  coverage <- dplyr::tibble()
  # score all submodules
  for (family in names(modules)) {
    family_dir <- file.path(out_dir, family)
    dir.create(family_dir, showWarnings = FALSE)
    for (sub in names(modules[[family]])) {
      mod_name <- paste(family, sub, sep="_")
      res <- score_module(obj, modules[[family]][[sub]], mod_name)
      obj <- res$obj
      present <- res$present; missing <- res$missing
      coverage <- bind_rows(coverage,
                            tibble(Module=mod_name,
                                   Present=length(present), Missing=length(missing),
                                   Present_Genes=paste(present, collapse=";"),
                                   Missing_Genes=paste(missing, collapse=";"))
      )
      plot_module(
        df = obj@meta.data %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[cluster_col]])),
        module_field = paste0("MS_", mod_name),
        folder = family_dir,
        title = paste(comparison_name, "-", family, "/", sub),
        paired_by = paired_by
      )
    }
  }
  # Composite autophagy index (AMPK - mTORC)
  if (all(c("MS_Autophagy_AMPK", "MS_Autophagy_mTORC") %in% colnames(obj@meta.data))) {
    obj$MS_Autophagy_RegIndex <- obj$MS_Autophagy_AMPK - obj$MS_Autophagy_mTORC
    plot_module(
      df = obj@meta.data %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[cluster_col]])),
      module_field = "MS_Autophagy_RegIndex",
      folder = file.path(out_dir, "Autophagy"),
      title = paste(comparison_name, "- Autophagy Regulation Index (AMPK - mTORC)"),
      paired_by = paired_by
    )
  }
  write_tsv(coverage, file.path(out_dir, "Module_Gene_Coverage.tsv"))
#  saveRDS(obj, file.path(out_dir, paste0("seu_with_all_submodules_", gsub("[^A-Za-z0-9]+","_",comparison_name), ".rds")))
  invisible(TRUE)
}
# -----------------------------
# Define subsets for 6 comparisons
# -----------------------------
# 1) HEI vs HEU (Entry only; unpaired)
hei_heu_pre_dir <- file.path(base_out, "HEIvsHEU_PreART")
TARA_pre <- subset(TARA_ALL, subset = Age <= 2 & Comparison_Group %in% c("HEI","HEU"))

# 2) HEU vs HUU (entry; unpaired)
heu_huu_dir <- file.path(base_out, "HEUvsHUU")
TARA_heu_huu <- subset(TARA_ALL, subset = Comparison_Group %in% c("HEU","HUU"))

# 3) High vs Low VL at Entry only (HEI; unpaired)
highlow_pre_dir <- file.path(base_out, "HighvsLowVL_PreART")
TARA_entry <- subset(TARA_ALL, subset = Age %in% c(1,2) & Condition == "HEI")
TARA_entry$Viral_Load_Category <- ifelse(TARA_entry$Viral_Load >= 100000, "High", "Low")
TARA_entry$Comparison_Group <- TARA_entry$Viral_Load_Category

# 4) High vs Low VL across ALL timepoints (HEI; paired by PID)
highlow_all_dir <- file.path(base_out, "HighvsLowVL_ALL")
TARA_hei_all <- subset(TARA_ALL, subset = Condition == "HEI")
TARA_hei_all$Viral_Load_Category <- ifelse(TARA_hei_all$Viral_Load >= 100000, "High", "Low")
TARA_hei_all$Comparison_Group <- TARA_hei_all$Viral_Load_Category

# 5) Post-ART Suppressed vs Pre-ART Entry (HEI; paired by PID)
post_supp_dir <- file.path(base_out, "PostART_Suppressed_vs_PreART")
TARA_HEI <- subset(TARA_ALL, subset = Condition == "HEI")
TARA_HEI$Timepoint_Group <- ifelse(
  TARA_HEI$Age <= 2, "PreART_Entry",
  ifelse(TARA_HEI$Viral_Load < 200, "PostART_Suppressed",
         ifelse(!is.na(TARA_HEI$Viral_Load) & TARA_HEI$Viral_Load >= 200, "PostART_Unsuppressed", NA))
)
TARA_HEI_PP <- subset(TARA_HEI, subset = Timepoint_Group %in% c("PreART_Entry","PostART_Suppressed"))
TARA_HEI_PP$Comparison_Group <- TARA_HEI_PP$Timepoint_Group

# 6) Post-ART Unsuppressed vs Pre-ART Entry (HEI; paired by PID)
post_unsupp_dir <- file.path(base_out, "PostART_Unsuppressed_vs_PreART")
TARA_HEI_PU <- subset(TARA_HEI, subset = Timepoint_Group %in% c("PreART_Entry","PostART_Unsuppressed"))
TARA_HEI_PU$Comparison_Group <- TARA_HEI_PU$Timepoint_Group

# 7) Post-ART Suppressed vs Post-ART Unsuppressed (HEI only)
#    -> folder: PostART_Suppressed_vs_Unsuppressed
# CONTROL for PID (longitudinal).
# =====================================================================
post_sup_unsup_dir <- file.path(base_out, "PostART_Suppressed_vs_Unsuppressed")

TARA_HEI_SU <- subset(TARA_HEI, subset = Timepoint_Group %in% c("PostART_Suppressed","PostART_Unsuppressed"))
TARA_HEI_SU$Comparison_Group <- TARA_HEI_SU$Timepoint_Group
colnames(TARA_ALL@meta.data)

# -----------------------------
# Run all comparisons
# -----------------------------
run_scoring(TARA_pre,      hei_heu_pre_dir,        "HEI vs HEU (Entry)",             paired_by=NULL)
run_scoring(TARA_heu_huu,  heu_huu_dir,            "HEU vs HUU (Entry)",             paired_by=NULL)
run_scoring(TARA_entry,    highlow_pre_dir,        "High vs Low VL (Entry, HEI)",    paired_by=NULL)











############## Paired Comparisons ##################
run_scoring(TARA_hei_all,  highlow_all_dir,        "High vs Low VL (All timepoints, HEI)", paired_by=pid_col)
run_scoring(TARA_HEI_PP,   post_supp_dir,          "Post-ART Suppressed vs Pre-ART Entry (HEI)", paired_by=pid_col)
run_scoring(TARA_HEI_PU,   post_unsupp_dir,        "Post-ART Unsuppressed vs Pre-ART Entry (HEI)", paired_by=pid_col)
run_scoring(TARA_HEI_SU,
            post_sup_unsup_dir,
            "Post-ART Suppressed vs Post-ART Unsuppressed (HEI)",
            paired_by = NULL)

