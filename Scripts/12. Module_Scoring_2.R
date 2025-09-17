
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)


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

###################################### FUNCTIONS FOR MODULE SCORING ##############################################
# ===========================================================
# 0) MODULE SCORING (run once on TARA_ALL to create MS_* cols)
#    - Creates one column per submodule: MS_<Family>_<Sub>
#    - Avoids per-cell pseudoreplication later by always aggregating to PID×cluster×group
# ===========================================================
score_one_module <- function(obj, genes, name, assay_use = DefaultAssay(obj)) {
  present <- intersect(genes, rownames(obj))
  col_out <- paste0("MS_", name)
  if (length(present) == 0) {
    obj[[col_out]] <- NA_real_
    return(list(obj = obj, present = character(0), missing = genes))
  }
  obj <- Seurat::AddModuleScore(obj, features = list(present), name = paste0(col_out, "_tmp"), assay = assay_use)
  tmp <- paste0(col_out, "_tmp1")
  obj[[col_out]] <- obj[[tmp]]
  obj[[tmp]] <- NULL
  list(obj = obj, present = present, missing = setdiff(genes, present))
}

add_all_module_scores <- function(obj, modules_list, assay_use = DefaultAssay(obj), out_dir = NULL, comparison_name = "ALL") {
  coverage <- tibble()
  for (family in names(modules_list)) {
    for (sub in names(modules_list[[family]])) {
      name <- paste(family, sub, sep = "_")
      res <- score_one_module(obj, modules_list[[family]][[sub]], name, assay_use = assay_use)
      obj <- res$obj
      coverage <- bind_rows(
        coverage,
        tibble(
          Comparison = comparison_name,
          Module = paste0("MS_", name),
          Present = length(res$present),
          Missing = length(res$missing),
          Present_Genes = paste(res$present, collapse = ";"),
          Missing_Genes = paste(res$missing, collapse = ";")
        )
      )
    }
  }
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(coverage, file.path(out_dir, paste0("Module_Gene_Coverage_", gsub("[^A-Za-z0-9]+","_", comparison_name), ".tsv")))
  }
  obj
}

# ===========================================================
# 1) COMMON UTILITIES
# ===========================================================
compute_pid_means <- function(obj, module_field, cluster_col, group_col, pid_col) {
  md <- obj@meta.data
  stopifnot(module_field %in% colnames(md))
  md %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[cluster_col]]), !is.na(.data[[module_field]])) %>%
    group_by(.data[[pid_col]], .data[[cluster_col]], .data[[group_col]]) %>%
    summarize(Score = mean(.data[[module_field]], na.rm = TRUE), n_cells = dplyr::n(), .groups = "drop") %>%
    rename(PID = !!pid_col, Cluster = !!cluster_col, Group = !!group_col)
}

# Decide mode per cluster (explicit justification string returned in note)
choose_mode_for_comparison <- function(obj,
                                       group_col = "Comparison_Group",
                                       pid_col   = "PID",
                                       min_paired = 10L,
                                       max_unpaired_prop = 0.10) {
  tab <- obj@meta.data %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[pid_col]])) %>%
    dplyr::distinct(.data[[pid_col]], .data[[group_col]]) %>%
    dplyr::add_count(.data[[pid_col]], name = "ng") %>%
    dplyr::summarize(
      n_paired   = sum(ng == 2),
      n_unpaired = sum(ng == 1),
      .groups = "drop"
    )
  
  n_paired <- tab$n_paired %||% 0L
  n_unpaired <- tab$n_unpaired %||% 0L
  n_total <- n_paired + n_unpaired
  unpaired_prop <- ifelse(n_total > 0, n_unpaired / n_total, 1)
  
  if (n_paired >= min_paired && unpaired_prop <= max_unpaired_prop)
    return(list(mode = "paired", note = paste0("Comparison-level: paired≥", min_paired, " & unpaired_prop≤", max_unpaired_prop)))
  if (n_paired == 0L)
    return(list(mode = "unpaired", note = "Comparison-level: no paired PIDs"))
  return(list(mode = "mixed", note = "Comparison-level: mixed paired+unpaired; using mixed-effects"))
}

# ===========================================================
# 2) TEST-SPECIFIC ANALYZERS (kept separate for easy debugging)
# ===========================================================
analyze_unpaired <- function(pid_means, groupA, groupB) {
  dfA <- pid_means %>% dplyr::filter(Group == groupA & !is.na(Score))
  dfB <- pid_means %>% dplyr::filter(Group == groupB & !is.na(Score))
  
  # need at least one observation per group to run unpaired Wilcoxon
  if (nrow(dfA) < 1L || nrow(dfB) < 1L) {
    return(tibble::tibble(
      test = "wilcoxon_unpaired",
      n_total = nrow(dfA) + nrow(dfB),
      n_groupA = nrow(dfA), n_groupB = nrow(dfB), n_paired = 0L,
      estimate = NA_real_, estimate_type = "HL_shift",
      ci_lower = NA_real_, ci_upper = NA_real_,
      statistic = NA_real_, p_value = NA_real_,
      notes = "Insufficient observations in one or both groups"
    ))
  }
  
  wt <- suppressWarnings(wilcox.test(dfB$Score, dfA$Score,
                                     paired = FALSE, conf.int = TRUE, exact = FALSE))
  tibble::tibble(
    test = "wilcoxon_unpaired",
    n_total = nrow(dfA) + nrow(dfB),
    n_groupA = nrow(dfA), n_groupB = nrow(dfB), n_paired = 0L,
    estimate = unname(wt$estimate), estimate_type = "HL_shift",
    ci_lower = if (!is.null(wt$conf.int)) wt$conf.int[1] else NA_real_,
    ci_upper = if (!is.null(wt$conf.int)) wt$conf.int[2] else NA_real_,
    statistic = unname(wt$statistic), p_value = unname(wt$p.value),
    notes = NA_character_
  )
}
# Paired Wilcoxon on deltas — magnitude = median delta (GroupB − GroupA within PID)
analyze_paired <- function(pid_means, groupA, groupB) {
  wide <- pid_means %>%
    select(PID, Group, Score) %>%
    tidyr::pivot_wider(names_from = Group, values_from = Score) %>%
    filter(!is.na(.data[[groupA]]), !is.na(.data[[groupB]]))
  if (nrow(wide) < 2) {
    return(tibble(
      test = "wilcoxon_paired", n_total = nrow(wide) * 2L, n_groupA = nrow(wide), n_groupB = nrow(wide),
      n_paired = nrow(wide), estimate = NA_real_, estimate_type = "median_delta",
      ci_lower = NA_real_, ci_upper = NA_real_, statistic = NA_real_, p_value = NA_real_,
      notes = "Insufficient paired PIDs (<2)"
    ))
  }
  wt <- suppressWarnings(wilcox.test(wide[[groupB]], wide[[groupA]], paired = TRUE, conf.int = TRUE, exact = FALSE))
  tibble(
    test          = "wilcoxon_paired",
    n_total       = nrow(wide) * 2L,
    n_groupA      = nrow(wide),
    n_groupB      = nrow(wide),
    n_paired      = nrow(wide),
    estimate      = unname(wt$estimate),    # location shift of paired deltas
    estimate_type = "median_delta",
    ci_lower      = if (!is.null(wt$conf.int)) wt$conf.int[1] else NA_real_,
    ci_upper      = if (!is.null(wt$conf.int)) wt$conf.int[2] else NA_real_,
    statistic     = unname(wt$statistic),
    p_value       = unname(wt$p.value),
    notes         = NA_character_
  )
}

# Mixed effects — magnitude = β(GroupB) on score scale; controls for repeated PID
# Assumes: library(lmerTest) is loaded
# Tests: Score ~ Group (+ optional covars) + (1|PID)
# Returns Satterthwaite p-values for Group (level = groupB vs groupA)

analyze_mixed <- function(pid_means, groupA, groupB, covars = NULL) {
  # prep data
  df <- pid_means %>%
    dplyr::filter(!is.na(Score), !is.na(Group), !is.na(PID)) %>%
    dplyr::mutate(Group = factor(Group, levels = c(groupA, groupB)))
  
  # quick guards
  if (dplyr::n_distinct(df$Group) < 2L) {
    return(tibble::tibble(
      test="mixed_effects(Satt)", n_total=nrow(df),
      n_groupA=sum(df$Group==groupA), n_groupB=sum(df$Group==groupB),
      n_paired=NA_integer_, estimate=NA_real_, estimate_type="beta_GroupB",
      ci_lower=NA_real_, ci_upper=NA_real_, statistic=NA_real_, p_value=NA_real_,
      notes="Only one group present"
    ))
  }
  if (dplyr::n_distinct(df$PID) < 2L) {
    return(tibble::tibble(
      test="mixed_effects(Satt)", n_total=nrow(df),
      n_groupA=sum(df$Group==groupA), n_groupB=sum(df$Group==groupB),
      n_paired=NA_integer_, estimate=NA_real_, estimate_type="beta_GroupB",
      ci_lower=NA_real_, ci_upper=NA_real_, statistic=NA_real_, p_value=NA_real_,
      notes="Fewer than 2 PIDs"
    ))
  }
  
  # count truly paired PIDs (present in both groups in this cluster)
  n_paired <- pid_means %>%
    dplyr::distinct(PID, Group) %>%
    dplyr::add_count(PID, name = "ng") %>%
    dplyr::summarise(n_paired = sum(ng == 2L), .groups = "drop") %>%
    dplyr::pull(n_paired)
  if (length(n_paired) == 0) n_paired <- NA_integer_
  
  # formula
  rhs <- c("Group", covars)
  rhs <- rhs[!is.na(rhs) & nzchar(rhs)]
  fml <- stats::as.formula(paste("Score ~", paste(rhs, collapse = " + "), "+ (1|PID)"))
  
  # fit (use bobyqa for stability)
  fit <- try(
    lmerTest::lmer(fml, data = df, REML = FALSE,
                   control = lme4::lmerControl(optimizer = "bobyqa",
                                               optCtrl   = list(maxfun = 1e5))),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(tibble::tibble(
      test="mixed_effects(Satt)", n_total=nrow(df),
      n_groupA=sum(df$Group==groupA), n_groupB=sum(df$Group==groupB),
      n_paired=as.integer(n_paired), estimate=NA_real_, estimate_type="beta_GroupB",
      ci_lower=NA_real_, ci_upper=NA_real_, statistic=NA_real_, p_value=NA_real_,
      notes=paste("lmer error:", conditionMessage(attr(fit, "condition")))
    ))
  }
  
  # extract coefficient for GroupB vs GroupA
  coef_name <- paste0("Group", groupB)
  coefs <- as.data.frame(summary(fit)$coefficients)
  
  if (!(coef_name %in% rownames(coefs))) {
    return(tibble::tibble(
      test="mixed_effects(Satt)", n_total=nrow(df),
      n_groupA=sum(df$Group==groupA), n_groupB=sum(df$Group==groupB),
      n_paired=as.integer(n_paired), estimate=NA_real_, estimate_type="beta_GroupB",
      ci_lower=NA_real_, ci_upper=NA_real_, statistic=NA_real_, p_value=NA_real_,
      notes="Contrast not found (check factor levels)"
    ))
  }
  
  est  <- coefs[coef_name, "Estimate"]
  se   <- coefs[coef_name, "Std. Error"]
  stat <- coefs[coef_name, "t value"]
  pval <- coefs[coef_name, "Pr(>|t|)"]
  ci   <- est + c(-1, 1) * 1.96 * se
  
  # mark singular fits (variance ~ 0) but still return p-value
  sing <- try(lme4::isSingular(fit, tol = 1e-4), silent = TRUE)
  note <- if (isTRUE(sing)) "singular_fit" else NA_character_
  
  tibble::tibble(
    test          = "mixed_effects(Satt)",
    n_total       = nrow(df),
    n_groupA      = sum(df$Group == groupA),
    n_groupB      = sum(df$Group == groupB),
    n_paired      = as.integer(n_paired),
    estimate      = est,
    estimate_type = "beta_GroupB",
    ci_lower      = ci[1],
    ci_upper      = ci[2],
    statistic     = stat,
    p_value       = pval,
    notes         = note
  )
}

# ===========================================================
# 3) OPTIONAL PLOTTING (PID means; lines only for truly paired)
# ===========================================================
plot_module_faceted <- function(pid_means_all, groupA, groupB,
                                module_field, title, out_png,
                                colors = c("#1f77b4", "#d62728"),
                                test_mode = c("mixed", "paired", "unpaired"),
                                label_digits = 3) {
  test_mode <- match.arg(test_mode)
  
  # Build per-facet labels with p-values using the chosen test across ALL facets
  labels_list <- list()
  clusters <- sort(unique(pid_means_all$Cluster))
  for (cl in clusters) {
    pid_means_cl <- pid_means_all %>% dplyr::filter(Cluster == cl)
    stats_row <- switch(
      test_mode,
      "mixed"    = analyze_mixed (pid_means_cl, groupA, groupB, covars = NULL),
      "paired"   = analyze_paired(pid_means_cl, groupA, groupB),
      "unpaired" = analyze_unpaired(pid_means_cl, groupA, groupB)
    )
    lbl <- if (is.na(stats_row$p_value)) {
      paste0("p = NA (", stats_row$test, ")")
    } else {
      paste0("p = ", formatC(stats_row$p_value, digits = label_digits, format = "f"),
             " (", stats_row$test, ")")
    }
    labels_list[[length(labels_list) + 1]] <- tibble::tibble(Cluster = cl, label = lbl)
  }
  labels_df <- dplyr::bind_rows(labels_list)
  
  # Mark paired PIDs (present in both groups) within each cluster for lines
  paired_tbl <- pid_means_all %>%
    dplyr::distinct(Cluster, PID, Group) %>%
    dplyr::add_count(Cluster, PID, name = "ng") %>%
    dplyr::filter(ng == 2) %>%
    dplyr::select(Cluster, PID) %>%
    dplyr::mutate(is_paired = TRUE)
  
  df <- pid_means_all %>%
    dplyr::left_join(paired_tbl, by = c("Cluster", "PID")) %>%
    dplyr::mutate(
      is_paired = tidyr::replace_na(is_paired, FALSE),
      Group = factor(Group, levels = c(groupA, groupB))
    )
  
  p <- ggplot(df, aes(x = Group, y = Score, color = Group, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(aes(group = PID),
               position = position_jitter(width = 0.15),
               size = 0.9, alpha = 0.6) +
    # lines only for PIDs present in both groups within that cluster
    geom_line(data = df %>% dplyr::filter(is_paired),
              aes(group = PID), alpha = 0.25, color = "grey40") +
    facet_wrap(~ Cluster, scales = "free_y") +
    scale_color_manual(values = colors, guide = "none") +
    scale_fill_manual(values = colors, guide = "none") +
    theme_minimal(base_size = 13) +
    labs(x = NULL, y = "Module score (PID mean)",
         title = title,
         subtitle = paste0(gsub("^MS_", "", module_field), "  •  test = ", test_mode))
  
  # Add per-facet p-value text (same test used everywhere in this plot)
  p <- p + geom_text(data = labels_df,
                     aes(x = 1.5, y = Inf, label = label),
                     inherit.aes = FALSE, vjust = 1.15, size = 3.6)
  
  ggsave(out_png, p, width = 20, height = 16, dpi = 300, bg = "white")
}


# ===========================================================
# 4) PER-COMPARISON DRIVER
#    - Loops modules × clusters
#    - Decides mixed vs paired vs unpaired (with justification)
#    - Saves one CSV per comparison (and optional plots)
# ===========================================================
# Helper to parse family/submodule from "MS_Family_Submodule"
.parse_family_sub <- function(module_field) {
  s <- sub("^MS_", "", module_field)
  family <- sub("_.*$", "", s)
  submod <- sub("^[^_]*_", "", s)
  list(family = family, sub = submod)
}

# Vectorized relabel using your label_map
.relabel_vec <- function(x, label_map) {
  x <- as.character(x)
  hits <- x %in% names(label_map)
  x[hits] <- unname(label_map[x[hits]])
  x
}

# Faceted plot per module (all clusters), with manual group renaming and colors.
# - Uses ONE test mode for the whole comparison (pass test_mode = "mixed" | "paired" | "unpaired")
# - Colors are fully manual: pass a 2-color vector, optionally NAMED by the (renamed) group labels.
# - You can optionally rename groups for plotting via `label_map` (e.g., HEI -> pHIV).
# - You can optionally control which group appears first on the x-axis (and thus which color is first)
#   via `group_levels`. If omitted, it uses the (renamed) order of groupA, groupB.
# - P-value labels are placed ABOVE the points with dynamic headroom; tweak `top_padding_mult` if needed.
plot_module_faceted <- function(pid_means_all, groupA, groupB,
                                module_field, title, out_png,
                                colors = c("#1f77b4", "#d62728"),
                                test_mode = c("mixed", "paired", "unpaired"),
                                label_digits = 3,
                                label_map = c(HEI = "pHIV", HEU = "pHEU", HUU = "pHUU"),
                                group_levels = NULL,
                                top_padding_mult = 0.25,
                                family_folders = TRUE,
                                module_families = NULL) {
  test_mode <- match.arg(test_mode)
  
  # ----- (0) Optional: infer family and rewrite out path into subfolder -----
  infer_family <- function(module_field, module_families) {
    base <- sub("^MS_", "", module_field)
    if (is.null(module_families) || length(module_families) == 0) {
      # Fallback: first token before "_" (may be imperfect for names with underscores)
      return(strsplit(base, "_", fixed = TRUE)[[1]][1])
    }
    # Longest family name that prefixes "<family>_" in the base
    hits <- module_families[vapply(
      module_families,
      function(f) startsWith(base, paste0(f, "_")),
      logical(1)
    )]
    if (length(hits) == 0) "Misc" else hits[which.max(nchar(hits))]
  }
  
  final_out_png <- out_png
  if (isTRUE(family_folders)) {
    fam <- infer_family(module_field, module_families)
    base_dir <- dirname(out_png)
    file_nm  <- basename(out_png)
    fam_dir  <- file.path(base_dir, fam)
    dir.create(fam_dir, showWarnings = FALSE, recursive = TRUE)
    final_out_png <- file.path(fam_dir, file_nm)
  } else {
    dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  }
  
  # ----- (1) Per-facet p-values using ONE chosen test for the whole comparison -----
  clusters <- sort(unique(pid_means_all$Cluster))
  labels_list <- vector("list", length(clusters))
  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    pid_means_cl <- pid_means_all %>% dplyr::filter(Cluster == cl)
    stats_row <- switch(
      test_mode,
      "mixed"    = analyze_mixed (pid_means_cl, groupA, groupB, covars = NULL),
      "paired"   = analyze_paired(pid_means_cl, groupA, groupB),
      "unpaired" = analyze_unpaired(pid_means_cl, groupA, groupB)
    )
    lbl <- if (is.na(stats_row$p_value)) {
      paste0("p = NA (", stats_row$test, ")")
    } else {
      paste0("p = ", formatC(stats_row$p_value, digits = label_digits, format = "f"),
             " (", stats_row$test, ")")
    }
    labels_list[[i]] <- tibble::tibble(Cluster = cl, label = lbl)
  }
  labels_df <- dplyr::bind_rows(labels_list)
  
  # ----- (2) Prepare plotting data (rename groups; mark paired PIDs for lines) -----
  paired_tbl <- pid_means_all %>%
    dplyr::distinct(Cluster, PID, Group) %>%
    dplyr::add_count(Cluster, PID, name = "ng") %>%
    dplyr::filter(ng == 2) %>%
    dplyr::select(Cluster, PID) %>%
    dplyr::mutate(is_paired = TRUE)
  
  df <- pid_means_all %>%
    dplyr::left_join(paired_tbl, by = c("Cluster", "PID")) %>%
    dplyr::mutate(
      is_paired  = tidyr::replace_na(is_paired, FALSE),
      GroupLabel = dplyr::recode(Group, !!!label_map, .default = Group)
    )
  
  # x-axis order (controls color order too)
  if (is.null(group_levels)) {
    group_levels <- c(
      dplyr::recode(groupA, !!!label_map, .default = groupA),
      dplyr::recode(groupB, !!!label_map, .default = groupB)
    )
  }
  group_levels <- unique(group_levels)
  df$GroupLabel <- factor(df$GroupLabel, levels = group_levels)
  
  # Named color mapping aligned to renamed levels
  if (is.null(names(colors))) {
    colors_named <- stats::setNames(colors[seq_along(group_levels)], group_levels)
  } else {
    colors_named <- colors[group_levels]
  }
  
  # ----- (3) Compute y positions for labels with extra headroom -----
  yrng <- df %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarize(y_min = min(Score, na.rm = TRUE),
                     y_max = max(Score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(y_pad = (y_max - y_min + 1e-8) * top_padding_mult,
                  y_lab = y_max + y_pad)
  
  labels_df <- labels_df %>% dplyr::left_join(yrng, by = "Cluster")
  
  # ----- (4) Plot -----
  p <- ggplot(df, aes(x = GroupLabel, y = Score, color = GroupLabel, fill = GroupLabel)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(aes(group = PID),
               position = position_jitter(width = 0.15),
               size = 0.9, alpha = 0.6) +
    geom_line(data = df %>% dplyr::filter(is_paired),
              aes(group = PID), alpha = 0.25, color = "grey40") +
    facet_wrap(~ Cluster, scales = "free_y") +
    scale_color_manual(values = colors_named, guide = "none", breaks = group_levels) +
    scale_fill_manual(values = colors_named, guide = "none", breaks = group_levels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, top_padding_mult + 0.05))) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(plot.margin = unit(c(10, 10, 10, 10), "pt")) +
    labs(x = NULL, y = "Module score (PID mean)",
         title = title,
         subtitle = paste0(gsub("^MS_", "", module_field), "  •  test = ", test_mode))
  
  # place p-value text above panels center
  p <- p + geom_text(data = labels_df,
                     aes(x = mean(c(1, length(group_levels))), y = y_lab, label = label),
                     inherit.aes = FALSE, size = 3.6, vjust = 0)
  
  ggsave(final_out_png, p, width = 20, height = 16, dpi = 300, bg = "white")
}
run_stats_for_comparison <- function(obj,
                                     comparison_name,
                                     out_dir,
                                     cluster_col   = "Manual_Annotation",
                                     group_col     = "Comparison_Group",
                                     pid_col       = "PID",
                                     covars        = NULL,
                                     thresholds    = list(min_paired = 10L, max_unpaired_prop = 0.10),
                                     do_plots      = TRUE,
                                     plot_colors   = c("#1f77b4", "#d62728"),
                                     force_plot_test_mode = c("auto", "mixed", "paired", "unpaired")) {
  
  force_plot_test_mode <- match.arg(force_plot_test_mode)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # auto-detect module columns that start with MS_
  module_fields <- grep("^MS_", colnames(obj@meta.data), value = TRUE)
  if (!length(module_fields)) stop("No MS_* module columns found. Run add_all_module_scores() first.")
  
  # require exactly two groups in this subset
  grps <- obj@meta.data %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::distinct(.data[[group_col]]) %>% dplyr::pull(1) %>% sort()
  if (length(grps) != 2) stop("Comparison must contain exactly two groups; found: ", paste(grps, collapse = ", "))
  groupA <- grps[1]; groupB <- grps[2]
  
  # Decide ONCE per comparison which test to use for plots
  if (force_plot_test_mode == "auto") {
    cm_cmp <- choose_mode_for_comparison(
      obj,
      group_col = group_col,
      pid_col = pid_col,
      min_paired = thresholds$min_paired,
      max_unpaired_prop = thresholds$max_unpaired_prop
    )
    plot_test_mode <- cm_cmp$mode
  } else {
    plot_test_mode <- force_plot_test_mode
  }
  
  all_stats <- list()
  
  for (module_field in module_fields) {
    # PID×cluster×group means for this module (used both for stats & plot)
    pid_means_all <- compute_pid_means(obj, module_field, cluster_col, group_col, pid_col)
    
    # ---------- STATS (per-cluster, hybrid-aware; unchanged) ----------
    for (cl in sort(unique(pid_means_all$Cluster))) {
      pid_means_cl <- pid_means_all %>% dplyr::filter(Cluster == cl)
      
      # Decide mode per cluster for the CSV
      cm <- choose_mode_for_cluster(
        pid_means_cluster = pid_means_cl,
        min_paired = thresholds$min_paired,
        max_unpaired_prop = thresholds$max_unpaired_prop
      )
      chosen_mode_csv <- cm$mode; note <- cm$note
      
      stats_row <- switch(
        chosen_mode_csv,
        "mixed"    = analyze_mixed (pid_means_cl, groupA, groupB, covars = covars),
        "paired"   = analyze_paired(pid_means_cl, groupA, groupB),
        "unpaired" = analyze_unpaired(pid_means_cl, groupA, groupB),
        stop("Unknown mode: ", chosen_mode_csv)
      )
      
      stats_row <- stats_row %>%
        dplyr::mutate(
          comparison = comparison_name,
          cluster    = cl,
          module     = module_field,
          mode_note  = note
        ) %>%
        dplyr::select(comparison, cluster, module, test, mode_note,
                      n_total, n_groupA, n_groupB, n_paired,
                      estimate, estimate_type, ci_lower, ci_upper,
                      statistic, p_value)
      
      all_stats[[length(all_stats) + 1]] <- stats_row
    }
    
    # ---------- PLOT (one faceted plot per module; one test for whole comparison) ----------
    if (do_plots) {
      pl_dir <- file.path(out_dir, "Plots")
      dir.create(pl_dir, recursive = TRUE, showWarnings = FALSE)
      
      plot_module_faceted(
        pid_means_all = pid_means_all,
        groupA = groupA, groupB = groupB,
        module_field = module_field,
        title = paste0(comparison_name, "  (faceted by cluster)"),
        out_png = file.path(pl_dir, paste0(gsub("[^A-Za-z0-9]+", "_", module_field), ".png")),
        colors = plot_colors,
        test_mode = plot_test_mode,     # <<< one test per comparison
        label_digits = 3
      )
    }
  }
  
  stats_tbl <- dplyr::bind_rows(all_stats) %>%
    dplyr::mutate(p_adj = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH"))) %>%
    dplyr::arrange(cluster, module)
  
  csv_path <- file.path(out_dir, paste0("Module_Stats_", gsub("[^A-Za-z0-9]+", "_", comparison_name), ".csv"))
  readr::write_csv(stats_tbl, csv_path)
  message("Saved stats CSV: ", csv_path)
  
  invisible(stats_tbl)
}

# ===========================================================
# 5) RUN: (A) Score once, then (B) analyze each of your subsets
#     — Pairing decisions are per cluster; justification in CSV
# ===========================================================

# (A) Score all modules once on the full object (writes coverage TSV in base_out)
TARA_ALL <- add_all_module_scores(
  obj = TARA_ALL,
  modules_list = modules,
  assay_use = DefaultAssay(TARA_ALL),
  out_dir = base_out,
  comparison_name = "ALL"
)

# -----------------------------
# Define subsets for 7 comparisons
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



# (B) ANALYZE EACH OF YOUR 7 COMPARISONS
# Heuristic thresholds for pairing (tune if needed)
thr <- list(min_paired = 10L, max_unpaired_prop = 0.10)

# 1) HEI vs HEU (Entry; cross-sectional → mostly unpaired)
run_stats_for_comparison(
  obj = TARA_pre,
  comparison_name = "HEI vs HEU (Entry)",
  out_dir = hei_heu_pre_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  thresholds = thr,
  do_plots = TRUE,
  # names MUST match the post-rename labels used in the plot
  plot_colors = c("pHIV" = "#d62728",  # red
                  "pHEU" = "#1f77b4"), # blue
  force_plot_test_mode = "auto"
)

# 2) HEU vs HUU (Entry; cross-sectional → unpaired)
run_stats_for_comparison(
  obj = TARA_heu_huu,
  comparison_name = "HEU vs HUU (Entry)",
  out_dir = heu_huu_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  thresholds = thr,
  do_plots = TRUE,
  plot_colors = c("pHEU" = "#1f77b4",
                  "pHUU" = "#2ca02c"),
  force_plot_test_mode = "auto"
)


# 3) High vs Low VL (Entry, HEI; cross-sectional → unpaired)
run_stats_for_comparison(
  obj = TARA_entry,
  comparison_name = "High vs Low VL (Entry, HEI)",
  out_dir = highlow_pre_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  thresholds = thr,
  do_plots = TRUE,
  plot_colors = c("High" = "#ff7f0e",
                  "Low"  = "#9467bd")
)

# 4) High vs Low VL (All timepoints, HEI; longitudinal → prefer mixed)
# Names must match the plot labels ("High", "Low")
run_stats_for_comparison(
  obj = TARA_hei_all,
  comparison_name = "High vs Low VL (All timepoints, HEI)",
  out_dir = highlow_all_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  covars = c("Age"),           # optional covariate
  thresholds = thr,
  do_plots = TRUE,
  plot_colors = c("High" = "#ff7f0e",   # orange
                  "Low"  = "#9467bd"),  # blue
  force_plot_test_mode = "auto"
)

# 5) Post-ART Suppressed vs Pre-ART Entry (HEI; mixed)
# Names must match the RENAMED labels: "Post-ART Suppressed", "Pre-ART Entry"
run_stats_for_comparison(
  obj = TARA_HEI_PP,
  comparison_name = "Post-ART Suppressed vs Pre-ART Entry (HEI)",
  out_dir = post_supp_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  covars = c("Age"),
  thresholds = thr,
  do_plots = TRUE,
  plot_colors = c("Post-ART Suppressed" = "#d62728",  # red
                  "Pre-ART Entry"       = "#1f77b4"), # blue
  force_plot_test_mode = "auto"
)

# 6) Post-ART Unsuppressed vs Pre-ART Entry (HEI; mixed)
# Names must match the RENAMED labels: "Post-ART Unsuppressed", "Pre-ART Entry"
run_stats_for_comparison(
  obj = TARA_HEI_PU,
  comparison_name = "Post-ART Unsuppressed vs Pre-ART Entry (HEI)",
  out_dir = post_unsupp_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  covars = c("Age"),
  thresholds = thr,
  do_plots = TRUE,
  plot_colors = c("Post-ART Unsuppressed" = "#d62728", # red
                  "Pre-ART Entry"         = "#1f77b4"),# blue
  force_plot_test_mode = "auto"
)

# 7) Post-ART Suppressed vs Post-ART Unsuppressed (HEI; mixed)
# Names must match the RENAMED labels: "Post-ART Suppressed", "Post-ART Unsuppressed"
run_stats_for_comparison(
  obj = TARA_HEI_SU,
  comparison_name = "Post-ART Suppressed vs Post-ART Unsuppressed (HEI)",
  out_dir = post_sup_unsup_dir,
  cluster_col = cluster_col, group_col = group_col, pid_col = pid_col,
  covars = c("Age"),
  thresholds = thr,
  do_plots = TRUE,
  plot_colors = c("Post-ART Unsuppressed" = "#d62728", # red
                  "Post-ART Suppressed"   = "#1f77b4"),# blue
  force_plot_test_mode = "auto"
)

)

############### DEBUGGING ##########
# Quick p-value sanity check (no plotting)
# - obj: a comparison object (e.g., TARA_pre, TARA_hei_all, etc.)
# - module_field: e.g., "MS_CD8_Effector"; if NULL, it auto-picks the first MS_ column
# - test: "mixed", "paired", or "unpaired"
# - covars: only used for mixed (e.g., covars = c("Age"))
quick_pval_check <- function(obj,
                             module_field = NULL,
                             test = c("mixed","paired","unpaired"),
                             cluster_col = "Manual_Annotation",
                             group_col   = "Comparison_Group",
                             pid_col     = "PID",
                             covars      = NULL,
                             top_k_clusters = 5) {
  test <- match.arg(test)
  
  # pick a module if not provided
  if (is.null(module_field)) {
    ms <- grep("^MS_", colnames(obj@meta.data), value = TRUE)
    if (length(ms) == 0) stop("No MS_* columns found. Run add_all_module_scores() first.")
    module_field <- ms[1]
  }
  
  # confirm two groups
  grps <- sort(unique(obj@meta.data[[group_col]][!is.na(obj@meta.data[[group_col]])]))
  if (length(grps) != 2) stop("Need exactly two groups in this subset; found: ", paste(grps, collapse=", "))
  groupA <- grps[1]; groupB <- grps[2]
  
  # aggregate to PID x cluster x group
  pm <- compute_pid_means(obj, module_field, cluster_col, group_col, pid_col)
  
  # pick top K clusters by N
  top_clusters <- pm %>% dplyr::count(Cluster, sort = TRUE) %>% dplyr::slice_head(n = top_k_clusters) %>% dplyr::pull(Cluster)
  
  message("Module: ", module_field, " | Groups: ", groupA, " vs ", groupB,
          " | Test: ", test, if (exists("has_lmerTest")) paste0(" | has_lmerTest=", has_lmerTest) else "")
  
  out <- lapply(top_clusters, function(cl) {
    pmc <- pm %>% dplyr::filter(Cluster == cl)
    
    # Run selected test
    res <- switch(test,
                  mixed    = analyze_mixed (pmc, groupA, groupB, covars = covars),
                  paired   = analyze_paired(pmc, groupA, groupB),
                  unpaired = analyze_unpaired(pmc, groupA, groupB)
    )
    
    tibble::tibble(
      cluster   = cl,
      test      = res$test[1],
      n_groupA  = res$n_groupA[1],
      n_groupB  = res$n_groupB[1],
      p_value   = res$p_value[1],
      estimate  = res$estimate[1],
      note      = if ("notes" %in% names(res)) res$notes[1] else NA_character_
    )
  })
  
  dplyr::bind_rows(out)
}
# Optional but recommended:
# install.packages("lmerTest"); library(lmerTest)

debug_mixed_one <- function(obj,
                            module_field = "MS_CD8_Effector",
                            cluster_name = NULL,
                            cluster_col = "Manual_Annotation",
                            group_col   = "Comparison_Group",
                            pid_col     = "PID") {
  
  # 1) Figure out the two groups and compute PID×cluster×group means
  grps <- obj@meta.data[[group_col]]
  grps <- sort(unique(grps[!is.na(grps)]))
  if (length(grps) != 2) stop("Need exactly two groups in this subset; found: ", paste(grps, collapse=", "))
  groupA <- grps[1]; groupB <- grps[2]
  
  pm <- compute_pid_means(obj, module_field, cluster_col, group_col, pid_col)  # your existing helper
  
  # 2) Pick cluster
  if (is.null(cluster_name)) {
    cluster_name <- pm %>% dplyr::count(Cluster, sort = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(Cluster)
  }
  df <- pm %>% dplyr::filter(Cluster == !!cluster_name)
  
  # 3) Basic diagnostics
  message("=== DEBUG MIXED ===")
  message("Module: ", module_field)
  message("Comparison groups: ", groupA, " vs ", groupB)
  message("Cluster: ", cluster_name)
  message("Rows: ", nrow(df), " | PIDs: ", length(unique(df$PID)))
  print(df %>% dplyr::count(Group, name="n_by_group"))
  
  # Factor with explicit order
  df$Group <- factor(df$Group, levels = c(groupA, groupB))
  
  # 4) Try to fit the mixed model; show real error if it fails
  fml <- as.formula(Score ~ Group + (1|PID))
  message("Formula: ", deparse(fml))
  
  fit_full <- tryCatch({
    if (exists("has_lmerTest") && isTRUE(has_lmerTest)) {
      lmerTest::lmer(fml, data = df, REML = FALSE)
    } else {
      lme4::lmer(fml, data = df, REML = FALSE,
                 control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    }
  }, error = function(e) e)
  
  if (inherits(fit_full, "error")) {
    message("lmer ERROR: ", conditionMessage(fit_full))
    return(invisible(
      tibble::tibble(cluster = cluster_name,
                     test = "mixed_effects",
                     n_groupA = sum(df$Group == groupA),
                     n_groupB = sum(df$Group == groupB),
                     p_value = NA_real_,
                     estimate = NA_real_,
                     note = paste("lmer error:", conditionMessage(fit_full)))
    ))
  }
  
  # 5) Extract estimate + p-value (Satterthwaite if available, else LRT)
  coef_name <- paste0("Group", groupB)
  est <- NA_real_; se <- NA_real_; pval <- NA_real_; stat <- NA_real_; label <- "mixed_effects"
  
  if (exists("has_lmerTest") && isTRUE(has_lmerTest)) {
    sm <- summary(fit_full); coefs <- as.data.frame(sm$coefficients)
    if (coef_name %in% rownames(coefs)) {
      est  <- coefs[coef_name, "Estimate"]; se <- coefs[coef_name, "Std. Error"]
      stat <- coefs[coef_name, "t value"];   pval <- coefs[coef_name, "Pr(>|t|)"]
      label <- "mixed_effects(Satt)"
    } else {
      message("Coefficient not found in summary; falling back to LRT.")
    }
  }
  
  if (is.na(pval)) {
    # LRT fallback
    fit_red <- tryCatch(lme4::lmer(update(fml, . ~ . - Group), data = df, REML = FALSE), error = function(e) e)
    if (!inherits(fit_red, "error")) {
      aa <- tryCatch(anova(fit_red, fit_full, test = "Chisq"), error = function(e) e)
      if (!inherits(aa, "error") && nrow(aa) >= 2) {
        pval <- suppressWarnings(aa$`Pr(>Chisq)`[2])
        stat <- suppressWarnings(aa$Chisq[2])
        label <- "mixed_effects(LRT)"
      } else {
        message("LRT failed: ", if (inherits(aa, "error")) conditionMessage(aa) else "unknown")
      }
    } else {
      message("Reduced model failed: ", conditionMessage(fit_red))
    }
  }
  
  # 6) Return a small tibble
  tibble::tibble(
    cluster  = cluster_name,
    test     = label,
    n_groupA = sum(df$Group == groupA),
    n_groupB = sum(df$Group == groupB),
    p_value  = pval,
    estimate = est
  )
}
# Largest cluster in the comparison:
debug_mixed_one(TARA_hei_all, module_field = "MS_CD8_Effector")

# Or specify a particular cluster name exactly as it appears in your meta:
library(lme4)
# If you want Satterthwaite p-values (instead of LRT), also load:
# library(lmerTest)

# Pick one module and one cluster from your object
module_field <- "MS_CD8_Effector"
cluster_name <- "0: CD4 T cell"

# Aggregate PID × cluster × group means
pm <- compute_pid_means(TARA_hei_all,
                        module_field = module_field,
                        cluster_col  = "Manual_Annotation",
                        group_col    = "Comparison_Group",
                        pid_col      = "PID")

# Subset to one cluster
df <- pm %>% dplyr::filter(Cluster == cluster_name)

# Make sure Group is a factor with consistent order
grps <- sort(unique(df$Group))
groupA <- grps[1]; groupB <- grps[2]
df <- df %>% dplyr::mutate(Group = factor(Group, levels = c(groupA, groupB)))

# Build formula
fml <- Score ~ Group + (1 | PID)

# ---- Run mixed model ----
fit_full <- lmer(fml, data = df, REML = FALSE)

# ---- Likelihood ratio test against reduced model ----
fit_null <- lmer(Score ~ 1 + (1 | PID), data = df, REML = FALSE)
lrt <- anova(fit_null, fit_full)

# ---- Results ----
print(summary(fit_full))   # fixed effect estimates
print(lrt)                 # LRT p-value
lrt$`Pr(>Chisq)`
lrt$`Pr(>Chisq)`[2]
# Extract p-value
lrt_p <- lrt$`Pr(>Chisq)`[2]
cat("Cluster:", cluster_name, "| LRT p-value:", lrt_p, "\n")

gnitude = β(GroupB) on score scale; controls for repeated PID
analyze_mixed <- function(pid_means, groupA, groupB, covars = NULL) {
  df <- pid_means %>% dplyr::mutate(Group = factor(Group, levels = c(groupA, groupB)))
  rhs <- c("Group", covars)
  fml <- as.formula(paste("Score ~", paste(rhs, collapse = " + "), "+ (1|PID)"))
  
  # Fit full model
  fit_full <- try({
    if (has_lmerTest) lmerTest::lmer(fml, data = df, REML = FALSE) else lme4::lmer(fml, data = df, REML = FALSE)
  }, silent = TRUE)
  
  if (inherits(fit_full, "try-error")) {
    return(tibble::tibble(
      test="mixed_effects", n_total=nrow(df),
      n_groupA=sum(df$Group==groupA), n_groupB=sum(df$Group==groupB),
      n_paired=NA_integer_, estimate=NA_real_, estimate_type="beta_GroupB",
      ci_lower=NA_real_, ci_upper=NA_real_, statistic=NA_real_, p_value=NA_real_,
      notes="Mixed model failed to fit"
    ))
  }
  
  coef_name <- paste0("Group", groupB)
  
  # Extract estimate + Wald CI
  if (has_lmerTest) {
    sm <- summary(fit_full); coefs <- as.data.frame(sm$coefficients)
    if (!(coef_name %in% rownames(coefs))) {
      est <- NA_real_; se <- NA_real_; stat <- NA_real_; p_satt <- NA_real_
    } else {
      est <- coefs[coef_name, "Estimate"]
      se  <- coefs[coef_name, "Std. Error"]
      stat <- coefs[coef_name, "t value"]
      p_satt <- coefs[coef_name, "Pr(>|t|)"]
    }
  } else {
    cf <- lme4::fixef(fit_full)
    est <- if (coef_name %in% names(cf)) unname(cf[coef_name]) else NA_real_
    V <- try(as.matrix(vcov(fit_full)), silent = TRUE)
    se <- if (!inherits(V, "try-error") && !is.null(colnames(V)) && coef_name %in% colnames(V)) sqrt(V[coef_name, coef_name]) else NA_real_
    stat <- NA_real_; p_satt <- NA_real_
  }
  ci <- if (!is.na(se)) est + c(-1, 1) * 1.96 * se else c(NA_real_, NA_real_)
  
  # P-value: Satterthwaite if available; otherwise LRT (full vs reduced)
  if (!is.na(p_satt)) {
    pval <- p_satt
    test_label <- "mixed_effects(Satt)"
    stat_out <- stat
  } else {
    fml_reduced <- update(fml, . ~ . - Group)
    fit_reduced <- try(lme4::lmer(fml_reduced, data = df, REML = FALSE), silent = TRUE)
    if (inherits(fit_reduced, "try-error")) {
      pval <- NA_real_; test_label <- "mixed_effects"; stat_out <- NA_real_
    } else {
      a <- try(anova(fit_reduced, fit_full, test = "Chisq"), silent = TRUE)
      if (inherits(a, "try-error") || nrow(a) < 2) {
        pval <- NA_real_; test_label <- "mixed_effects"; stat_out <- NA_real_
      } else {
        pval <- suppressWarnings(a$`Pr(>Chisq)`[2])
        stat_out <- suppressWarnings(a$Chisq[2])
        test_label <- "mixed_effects(LRT)"
      }
    }
  }
  
  tibble::tibble(
    test          = test_label,
    n_total       = nrow(df),
    n_groupA      = sum(df$Group == groupA),
    n_groupB      = sum(df$Group == groupB),
    n_paired      = NA_integer_,just run
    estimate      = est,
    estimate_type = "beta_GroupB",
    ci_lower      = ci[1],
    ci_upper      = ci[2],
    statistic     = stat_out,
    p_value       = pval,
    notes         = ifelse(has_lmerTest, NA_character_, "lmerTest not loaded; LRT p-value used")
  )
}
