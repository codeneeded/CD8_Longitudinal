################################################################################
# TARA — Module Scoring Analysis
#
# Runs AFTER the annotation pipeline. Uses:
#   - TARA_ALL_annotated_final.qs2 (or checkpoint.qs2)
#   - DGE CSVs from Step 12 of annotation pipeline
#
# STRUCTURE:
#   PHASE 1: Module scoring + CSV checkpoint (review before continuing)
#   PHASE 2: Visualizations (volcanos, violins, heatmaps, correlations)
#   PHASE 3: Ratio analyses (IFN Memory/Stress, OXPHOS/Glycolysis)
#
# OUTPUT: ~/Documents/CD8_Longitudinal/Analysis/Module_Scoring_Analysis/
################################################################################


# =============================================================================
# STEP 0: LIBRARIES + PATHS
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(cowplot)
library(patchwork)
library(scales)
library(viridis)
library(grid)
library(rstatix)
library(qs2)

base_dir  <- "~/Documents/CD8_Longitudinal"
saved_dir <- file.path(base_dir, "saved_R_data")
dge_dir   <- file.path(base_dir, "Analysis", "Differential_Expression")

# Output tree
out_base      <- file.path(base_dir, "Analysis", "Module_Scoring_Analysis")
out_csv       <- file.path(out_base, "CSV_Checkpoints")
out_volcano   <- file.path(out_base, "Volcanos")
out_violin    <- file.path(out_base, "Violins")
out_heatmap   <- file.path(out_base, "Heatmaps")
out_corr      <- file.path(out_base, "Correlations")
out_ratio     <- file.path(out_base, "Ratio_Analysis")
out_receptor  <- file.path(out_base, "IFNAR_Receptor")

for (d in c(out_csv, out_volcano, out_violin, out_heatmap,
            out_corr, out_ratio, out_receptor)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("══ STEP 1: LOADING ══\n\n")

TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_annotated_final.qs2"))
DefaultAssay(TARA_ALL) <- "RNA"
Idents(TARA_ALL) <- "Annotation"

cat("Cells:", ncol(TARA_ALL), "\n")
cat("Clusters:", length(unique(TARA_ALL$Annotation)), "\n\n")


# =============================================================================
# STEP 2: MODULE DEFINITIONS (17 modules)
# =============================================================================

cat("══ STEP 2: MODULE DEFINITIONS ══\n\n")

modules <- list(
  
  # ── IFN Programmes (manuscript Table 1) ──
  
  IFN_Memory = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2",
                 "OAS1", "OAS2", "OAS3", "IFI6"),
  
  Stress_IFN_Acute = c("HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1",
                       "IFI27", "IFI44L"),
  
  IFN_Signalling = c("STAT1", "STAT2", "IRF7", "IRF9", "IRF1", "IRF3"),
  
  IFN_Negative_Reg = c("USP18", "SOCS1", "SOCS3"),
  
  # ── T cell functional programmes (manuscript Table 1) ──
  
  Exhaustion = c("TOX", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
                 "CTLA4", "ENTPD1"),
  
  Stemness_Naive = c("TCF7", "LEF1", "SELL", "CCR7", "BACH2",
                     "IL7R", "BCL2", "S1PR1", "KLF2"),
  
  Cytotoxicity = c("GZMB", "GZMA", "GZMH", "GNLY", "PRF1",
                   "NKG7", "FGFBP2"),
  
  Inflammatory_Chemokines = c("CCL3", "CCL4", "CCL4L2", "CCL5",
                              "XCL1", "XCL2"),
  
  Terminal_Differentiation = c("ZEB2", "PRDM1", "TBX21", "CX3CR1",
                               "S1PR5", "ID2", "EOMES"),
  
  # ── Activation & Cytokines ──
  
  Activation = c("HLA-DRA", "HLA-DRB1", "CD38", "FAS", "CD69",
                 "ICOS", "TNFRSF9"),
  
  Inflammatory_Cytokines = c("IL1B", "TNF", "IL6", "CXCL8", "IL18",
                             "CXCL10", "CXCL9", "CCL2", "CSF2", "IFNG"),
  
  # ── Metabolism ──
  
  Glycolysis_Hypoxia = c("HIF1A", "LDHA", "PKM", "GAPDH", "ENO1",
                         "PFKFB3", "PFKFB4", "SLC2A1", "SLC2A3", "ALDOA",
                         "HILPDA", "EGLN3", "IGFBP2", "VEGFA", "AK4"),
  
  OXPHOS = c("ATP5F1A", "ATP5F1B", "NDUFA4", "NDUFB8", "SDHB",
             "UQCRB", "COX7A2", "COX5A", "CS", "IDH2",
             "OGDH", "MDH2", "VDAC1", "SOD2"),
  
  Fatty_Acid_Oxidation = c("CPT1A", "ACADM", "HADHA", "ACADVL",
                           "PPARA", "PPARGC1A", "ACSL1"),
  
  # ── Cell fate & survival ──
  
  Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM7", "MCM3",
                    "TYMS", "CDK1", "CCNB1", "CCNA2", "BIRC5"),
  
  Apoptosis_Survival = c("BCL2", "BCL2L1", "MCL1", "BCL2A1",
                         "BAX", "BAK1", "BID", "BBC3",
                         "FAS", "FASLG", "CASP3", "CASP8"),
  
  Senescence = c("KLRG1", "CDKN1A", "CDKN2A", "B3GAT1",
                 "TP53", "RB1", "LMNB1", "HMGA1"),
  
  # ── Tissue homing ──
  
  Tissue_Residency = c("ITGAE", "ITGA1", "CD69", "ITGB7",
                       "ZNF683", "PRDM1", "CXCR6", "RUNX3")
)


# =============================================================================
# STEP 2b: MODULE DOCUMENTATION
# =============================================================================

module_docs <- c(
  "=" %>% strrep(79),
  "MODULE DEFINITIONS AND INTERPRETATION GUIDE",
  "TARA Cohort — Module Scoring Analysis",
  "=" %>% strrep(79),
  "",
  "This file documents all 17 functional gene modules scored in the analysis.",
  "For each module: gene list, biological meaning, and interpretation guide.",
  "",
  paste0("─" %>% strrep(79)),
  "",
  
  # IFN Memory
  "MODULE: IFN_Memory",
  "GENES: IFIT1, IFIT2, IFIT3, ISG15, MX1, MX2, OAS1, OAS2, OAS3, IFI6",
  "",
  "The 'good' IFN signature. Cells that were exposed to type I interferon",
  "during initial antigen encounter and retain a poised programme for rapid",
  "recall. These genes encode antiviral effector proteins: OAS enzymes activate",
  "RNase L to degrade viral RNA, MX GTPases restrict viral replication, IFIT",
  "proteins inhibit viral translation. Suppressed expanding CD8 clones showed",
  "significantly higher IFN Memory scores than unsuppressed (d = 0.25-0.43).",
  "Type I IFN is the critical 'third signal' for CD8 priming — without it,",
  "expansion and memory generation drop by 99%.",
  "",
  "INTERPRET: High = IFN-primed readiness for rapid recall upon re-encounter.",
  "",
  
  # Stress/IFN Acute
  "MODULE: Stress_IFN_Acute",
  "GENES: HSPA1A, HSPA1B, HSP90AA1, DNAJB1, IFI27, IFI44L",
  "",
  "The 'bad' IFN signature. Heat shock proteins (HSPA1A/B, HSP90AA1, DNAJB1)",
  "indicate protein misfolding stress from active viral replication. IFI27 and",
  "IFI44L are acute-phase ISGs marking cells under current IFN bombardment.",
  "This was the DOMINANT signature distinguishing suppressed from unsuppressed",
  "clones (d = -0.61 to -0.68 — larger than exhaustion).",
  "",
  "INTERPRET: High = active cellular damage from ongoing viral replication.",
  "KEY: High Memory + Low Stress = optimally primed. High Stress + Low Memory",
  "= viral damage without functional priming.",
  "",
  
  # IFN Signalling
  "MODULE: IFN_Signalling",
  "GENES: STAT1, STAT2, IRF7, IRF9, IRF1, IRF3",
  "",
  "The signal transduction machinery downstream of IFNAR1/IFNAR2 receptor.",
  "These TFs activate ISG expression. Distinguishes cells that are CAPABLE of",
  "IFN signalling (machinery present) from those actively RESPONDING (ISGs up).",
  "IRF7 is the master regulator of type I IFN production.",
  "",
  "INTERPRET: High = IFN signalling machinery active / IFN-responsive cell.",
  "",
  
  # IFN Negative Reg
  "MODULE: IFN_Negative_Reg",
  "GENES: USP18, SOCS1, SOCS3",
  "",
  "Feedback inhibitors that dampen IFN signalling. USP18 blocks IFNAR2;",
  "SOCS1/3 inhibit JAK-STAT. Could indicate appropriate feedback control or",
  "IFN resistance. Context with IFN Memory and Stress scores determines",
  "interpretation.",
  "",
  "INTERPRET: High + High Stress = trying to shut down damaging IFN.",
  "High + Low Stress = effective negative regulation (healthy feedback).",
  "",
  
  # Exhaustion
  "MODULE: Exhaustion",
  "GENES: TOX, PDCD1, HAVCR2, TIGIT, LAG3, CTLA4, ENTPD1",
  "",
  "Progressive loss of effector function through chronic antigen stimulation.",
  "TOX is the master regulator. Exhaustion correlated directly with viral load",
  "(rho = 0.70, p = 5.9e-4). ART actively reduced exhaustion below pre-ART.",
  "",
  "INTERPRET: High = chronic antigen-driven dysfunction. Track with viral load.",
  "",
  
  # Stemness
  "MODULE: Stemness_Naive",
  "GENES: TCF7, LEF1, SELL, CCR7, BACH2, IL7R, BCL2, S1PR1, KLF2",
  "",
  "Self-renewal capacity and quiescent memory potential. TCF-1 (TCF7) marks",
  "stem-like CD8 T cells that sustain long-term antiviral immunity. Its",
  "frequency predicts post-intervention viral control. ART preserved stemness",
  "in expanding clones (d = 0.35-0.46 vs unsuppressed).",
  "",
  "INTERPRET: High = self-renewal capacity preserved, good for cure strategies.",
  "",
  
  # Cytotoxicity
  "MODULE: Cytotoxicity",
  "GENES: GZMB, GZMA, GZMH, GNLY, PRF1, NKG7, FGFBP2",
  "",
  "Direct target cell killing through granzyme/perforin. Highest in pre-ART",
  "cells (maximally activated), partially resolved by ART.",
  "",
  "INTERPRET: High = active killing capacity. Context matters — high in a",
  "naive cell is abnormal, high in TEMRA/CTL is expected.",
  "",
  
  # Inflammatory Chemokines
  "MODULE: Inflammatory_Chemokines",
  "GENES: CCL3, CCL4, CCL4L2, CCL5, XCL1, XCL2",
  "",
  "Capacity to recruit other immune cells via MIP-1alpha/beta and lymphotactin.",
  "CCL3/CCL4 are also natural CCR5 ligands (the HIV co-receptor), so their",
  "expression has implications for HIV susceptibility of neighbouring cells.",
  "",
  "INTERPRET: High = actively recruiting immune cells / blocking CCR5.",
  "",
  
  # Terminal Differentiation
  "MODULE: Terminal_Differentiation",
  "GENES: ZEB2, PRDM1, TBX21, CX3CR1, S1PR5, ID2, EOMES",
  "",
  "Commitment to terminal effector fate with loss of plasticity. End-stage of",
  "the exhaustion trajectory where cells lose ability to respond to therapy.",
  "",
  "INTERPRET: High = terminally committed, limited therapeutic potential.",
  "",
  
  # Activation
  "MODULE: Activation",
  "GENES: HLA-DRA, HLA-DRB1, CD38, FAS, CD69, ICOS, TNFRSF9",
  "",
  "General immune activation. HLA-DR/CD38 co-expression is the classic measure",
  "of T cell activation in HIV. Persistent activation despite ART reflects",
  "ongoing viral replication or microbial translocation.",
  "",
  "INTERPRET: High = recently activated. Persistent high = pathological.",
  "",
  
  # Inflammatory Cytokines
  "MODULE: Inflammatory_Cytokines",
  "GENES: IL1B, TNF, IL6, CXCL8, IL18, CXCL10, CXCL9, CCL2, CSF2, IFNG",
  "",
  "Broad inflammatory cytokine production. Includes innate (IL-1beta, TNF,",
  "IL-6) and adaptive (IFN-gamma) cytokines. CXCL10 (IP-10) is used",
  "clinically as a biomarker of IFN-driven inflammation in HIV.",
  "",
  "INTERPRET: High = systemic inflammation. Expected in monocytes/pDCs.",
  "",
  
  # Glycolysis
  "MODULE: Glycolysis_Hypoxia",
  "GENES: HIF1A, LDHA, PKM, GAPDH, ENO1, PFKFB3, PFKFB4, SLC2A1, SLC2A3,",
  "       ALDOA, HILPDA, EGLN3, IGFBP2, VEGFA, AK4",
  "",
  "Whether cells have switched to anaerobic glucose metabolism. T cells do",
  "this normally upon activation, but exhausted cells get stuck in glycolysis",
  "and cannot switch back to mitochondrial metabolism — this prevents memory",
  "formation. A pre-ART naive CD8 cluster already showed this signature",
  "(PFKFB4, HILPDA, EGLN3, IGFBP2 elevated).",
  "",
  "INTERPRET: High in naive cells = abnormal (viremia-driven metabolic stress).",
  "High in effectors = expected. High in exhausted cells = metabolic trap.",
  "",
  
  # OXPHOS
  "MODULE: OXPHOS",
  "GENES: ATP5F1A, ATP5F1B, NDUFA4, NDUFB8, SDHB, UQCRB, COX7A2, COX5A,",
  "       CS, IDH2, OGDH, MDH2, VDAC1, SOD2",
  "",
  "Mitochondrial fitness — cells burning fuel efficiently through the electron",
  "transport chain and TCA cycle. This is the healthy metabolic state for",
  "stem-like and memory T cells. If ART preserves stemness, OXPHOS should",
  "track with it. Cells with high OXPHOS + high stemness have the metabolic",
  "capacity to survive long-term.",
  "",
  "INTERPRET: High = metabolically fit. Predict: suppressed > unsuppressed.",
  "The OXPHOS/Glycolysis ratio is as informative as the IFN Memory/Stress ratio.",
  "",
  
  # FAO
  "MODULE: Fatty_Acid_Oxidation",
  "GENES: CPT1A, ACADM, HADHA, ACADVL, PPARA, PPARGC1A, ACSL1",
  "",
  "The specific fuel source for memory T cell survival — burning fat to power",
  "mitochondria. CPT1A is the gatekeeper enzyme and a single-gene proxy for",
  "memory metabolic fitness. PGC-1alpha (PPARGC1A) is the master regulator",
  "of mitochondrial biogenesis.",
  "",
  "INTERPRET: High = using the memory fuel programme. Should track with TCF7.",
  "",
  
  # Proliferation
  "MODULE: Proliferation",
  "GENES: MKI67, TOP2A, PCNA, MCM7, MCM3, TYMS, CDK1, CCNB1, CCNA2, BIRC5",
  "",
  "Active cell division — S phase and G2/M markers. MKI67 is the classic",
  "proliferation marker. Distinguishes actively expanding populations from",
  "quiescent ones under different ART conditions.",
  "",
  "INTERPRET: High = actively dividing. In the context of ART, sustained",
  "proliferation in effector populations may indicate ongoing antigenic drive.",
  "",
  
  # Apoptosis/Survival
  "MODULE: Apoptosis_Survival",
  "GENES: BCL2, BCL2L1, MCL1, BCL2A1 (pro-survival)",
  "       BAX, BAK1, BID, BBC3 (pro-apoptosis)",
  "       FAS, FASLG, CASP3, CASP8 (death pathway)",
  "",
  "Balance between cell survival and programmed cell death. Pro-survival genes",
  "(BCL2 family) vs pro-apoptosis effectors and death receptors. Stem-like",
  "cells typically have high BCL2 (survival advantage). Net interpretation",
  "depends on the balance — a high overall score could reflect either active",
  "apoptosis or strong survival signalling.",
  "",
  "INTERPRET: Check individual gene contributions. High BCL2 + low BAX = survival.",
  "High FAS + CASP = apoptosis-prone. Will export individual gene scores for this.",
  "",
  
  # Senescence
  "MODULE: Senescence",
  "GENES: KLRG1, CDKN1A (p21), CDKN2A (p16), B3GAT1 (CD57),",
  "       TP53, RB1, LMNB1, HMGA1",
  "",
  "Replicative senescence — cells that have permanently exited the cell cycle",
  "but remain metabolically active. KLRG1 and CD57 (B3GAT1) are surface markers.",
  "p21 (CDKN1A) and p16 (CDKN2A) are the canonical cell cycle inhibitors.",
  "Senescence is distinct from exhaustion: senescent cells cannot proliferate",
  "but may retain effector function, while exhausted cells lose function but",
  "may still divide.",
  "",
  "INTERPRET: High = terminally arrested, no proliferative potential. Important",
  "for cure strategies — senescent HIV-specific cells cannot expand on demand.",
  "",
  
  # Tissue Residency
  "MODULE: Tissue_Residency",
  "GENES: ITGAE (CD103), ITGA1 (CD49a), CD69, ITGB7, ZNF683 (Hobit),",
  "       PRDM1 (Blimp-1), CXCR6, RUNX3",
  "",
  "Tissue-resident memory programme — cells anchored in tissues rather than",
  "circulating. CD103 and CD49a are adhesion molecules for tissue retention.",
  "CD69 blocks S1PR1-mediated tissue egress. Relevant because mucosal",
  "tissue-resident memory is key for HIV control at sites of viral replication.",
  "",
  "INTERPRET: High in PBMC = cells with tissue-homing programme (will migrate",
  "to tissues). The MAIT cells showed this signature (CD103+CD49a+).",
  "",
  
  paste0("─" %>% strrep(79)),
  "",
  "RATIO METRICS:",
  "",
  "IFN_Ratio = IFN_Memory - Stress_IFN_Acute",
  "  Positive = IFN-primed readiness (good). Negative = viral damage (bad).",
  "",
  "Metabolic_Ratio = OXPHOS - Glycolysis_Hypoxia",
  "  Positive = mitochondrial fitness (memory-like). Negative = glycolytic trap.",
  "",
  "Together these two ratios capture whether ART restores both immunological",
  "AND metabolic fitness — the two pillars of long-lived functional immunity.",
  ""
)

writeLines(module_docs,
           file.path(out_base, "Module_Definitions_and_Interpretation_Guide.txt"))
cat("Module documentation written\n\n")


# =============================================================================
# STEP 3: SCORE ALL MODULES
# =============================================================================

cat("══ STEP 3: MODULE SCORING ══\n\n")

DefaultAssay(TARA_ALL) <- "RNA"

for (m in names(modules)) {
  genes_present <- intersect(modules[[m]], rownames(TARA_ALL[["RNA"]]))
  cat(sprintf("  %s: %d/%d genes found\n", m, length(genes_present), length(modules[[m]])))
  
  if (length(genes_present) < 2) {
    TARA_ALL[[paste0("MS_", m)]] <- NA_real_
    next
  }
  
  TARA_ALL <- AddModuleScore(
    TARA_ALL, features = list(genes_present),
    name = paste0("MS_", m), assay = "RNA"
  )
  col_from <- paste0("MS_", m, "1")
  TARA_ALL[[paste0("MS_", m)]] <- TARA_ALL[[col_from]]
  TARA_ALL[[col_from]] <- NULL
}

# Ratio metrics
TARA_ALL$IFN_Ratio       <- TARA_ALL$MS_IFN_Memory - TARA_ALL$MS_Stress_IFN_Acute
TARA_ALL$Metabolic_Ratio <- TARA_ALL$MS_OXPHOS - TARA_ALL$MS_Glycolysis_Hypoxia

cat("\n✓ All 17 RNA modules scored + 2 ratios computed\n\n")


# =============================================================================
# STEP 3b: PROTEIN MODULE SCORING (DSB-normalized ADT)
# =============================================================================

cat("══ STEP 3b: PROTEIN MODULE SCORING ══\n\n")

# Protein modules: only include ADT antibodies that exist in the panel
# Use ADT feature names (which may differ slightly from gene symbols)
DefaultAssay(TARA_ALL) <- "ADT"
adt_features <- rownames(TARA_ALL[["ADT"]])

# Map: module name → ADT antibody names
# (ADT names must match exactly what's in the assay)
protein_modules <- list(
  
  Exhaustion_Protein = c("PDCD1", "HAVCR2", "TIGIT", "LAG3", "CTLA4", "ENTPD1"),
  # Missing: TOX (intracellular TF, no surface antibody)
  # Coverage: 6/7 genes — excellent
  
  Activation_Protein = c("HLA-DRA", "CD38", "FAS", "CD69", "ICOS", "TNFRSF9"),
  # Missing: HLA-DRB1 (covered by HLA-DRA)
  # Coverage: 6/7 — excellent
  
  Stemness_Protein = c("SELL", "CCR7", "IL7R"),
  # Missing: TCF7, LEF1, BACH2, BCL2, S1PR1, KLF2 (all intracellular)
  # Coverage: 3/9 — but these ARE the canonical surface markers
  
  Tissue_Residency_Protein = c("ITGAE", "ITGA1", "CD69", "ITGB7"),
  # Missing: ZNF683, PRDM1, CXCR6, RUNX3 (intracellular TFs)
  # Coverage: 4/8 — all adhesion/retention markers present
  
  Senescence_Protein = c("KLRG1", "B3GAT1"),
  # Missing: CDKN1A, CDKN2A, TP53, RB1, LMNB1, HMGA1 (all intracellular)
  # Coverage: 2/8 — but these are THE flow cytometry senescence markers
  
  Terminal_Diff_Protein = c("CX3CR1")
  # Coverage: 1/7 — minimal, but CX3CR1 is the canonical TEMRA marker
)

# Check which antibodies actually exist and score
for (m in names(protein_modules)) {
  ab_present <- intersect(protein_modules[[m]], adt_features)
  cat(sprintf("  %s: %d/%d antibodies found", m, length(ab_present),
              length(protein_modules[[m]])))
  
  if (length(ab_present) == 0) {
    TARA_ALL[[paste0("PS_", m)]] <- NA_real_
    cat(" — SKIPPED\n")
    next
  }
  
  # Print which are found/missing
  ab_missing <- setdiff(protein_modules[[m]], adt_features)
  if (length(ab_missing) > 0) cat(sprintf(" (missing: %s)", paste(ab_missing, collapse = ", ")))
  cat("\n")
  
  # Protein score = mean of DSB-normalized values per cell
  # (NOT AddModuleScore which assumes log-normalized RNA)
  adt_mat <- GetAssayData(TARA_ALL, assay = "ADT", slot = "data")[ab_present, , drop = FALSE]
  TARA_ALL[[paste0("PS_", m)]] <- colMeans(as.matrix(adt_mat))
}

DefaultAssay(TARA_ALL) <- "RNA"

cat("\n✓ Protein module scores computed (PS_ prefix)\n\n")

# =============================================================================
# STEP 3c: CREATE COMPOSITE COMPARISON COLUMN
# =============================================================================

# Timepoint_Group only covers HEI (PreART_Entry, PostART_Suppressed, PostART_Unsuppressed).
# HEU and HUU have NA. For cross-group comparisons at matched timepoints
# we need a single column that spans all groups.
# NOTE: Using base R — case_when can have scoping issues with Seurat metadata.

TARA_ALL$Exposure_ART <- TARA_ALL$Timepoint_Group  # copies HEI groups, NA for rest
ea <- TARA_ALL$Exposure_ART
cond <- TARA_ALL$Condition
ea[is.na(ea) & cond == "HEU"] <- "HEU"
ea[is.na(ea) & cond == "HUU"] <- "HUU"
TARA_ALL$Exposure_ART <- ea

cat("Exposure_ART levels:\n")
print(table(TARA_ALL$Exposure_ART, useNA = "ifany"))
cat("\n")


# =============================================================================
# STEP 4: CHECKPOINT — EXPORT CSVs
# =============================================================================

cat("══ STEP 4: CSV CHECKPOINT ══\n\n")

ms_cols <- c(grep("^MS_", colnames(TARA_ALL@meta.data), value = TRUE),
             grep("^PS_", colnames(TARA_ALL@meta.data), value = TRUE),
             "IFN_Ratio", "Metabolic_Ratio")

meta_cols <- c("Annotation", "orig.ident", "Condition",
               "Timepoint_Group", "Exposure_ART", "Viral_Load_num", "Age")

cell_scores <- TARA_ALL@meta.data[, c(meta_cols, ms_cols)]

# 4a: Per-sample per-cluster means
sample_means <- cell_scores %>%
  group_by(orig.ident, Annotation, Condition, Timepoint_Group, Viral_Load_num) %>%
  summarise(across(all_of(ms_cols), ~ mean(.x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop")
write.csv(sample_means,
          file.path(out_csv, "Module_Scores_PerSample_PerCluster.csv"),
          row.names = FALSE)

# 4b: Per-cluster by condition
cluster_cond <- cell_scores %>%
  group_by(Annotation, Condition) %>%
  summarise(across(all_of(ms_cols), ~ mean(.x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop")
write.csv(cluster_cond,
          file.path(out_csv, "Module_Scores_PerCluster_ByCondition.csv"),
          row.names = FALSE)

# 4c: Per-cluster by timepoint
cluster_tp <- cell_scores %>%
  filter(!is.na(Timepoint_Group)) %>%
  group_by(Annotation, Timepoint_Group) %>%
  summarise(across(all_of(ms_cols), ~ mean(.x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop")
write.csv(cluster_tp,
          file.path(out_csv, "Module_Scores_PerCluster_ByTimepoint.csv"),
          row.names = FALSE)

# 4c2: Per-cluster by Exposure_ART (the composite column — most useful)
cluster_ea <- cell_scores %>%
  filter(!is.na(Exposure_ART)) %>%
  group_by(Annotation, Exposure_ART) %>%
  summarise(across(all_of(ms_cols), ~ mean(.x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop")
write.csv(cluster_ea,
          file.path(out_csv, "Module_Scores_PerCluster_ByExposureART.csv"),
          row.names = FALSE)

# 4d: Pairwise stats — all modules × clusters × comparisons
# Using Exposure_ART column for cross-group comparisons at matched timepoints
comparisons_list <- list(
  # Cross-group (matched early timepoint ~1-2 months)
  list(col = "Exposure_ART", g1 = "PreART_Entry", g2 = "HEU",
       label = "PreART_Entry_vs_HEU"),
  list(col = "Exposure_ART", g1 = "HEU", g2 = "HUU",
       label = "HEU_vs_HUU"),
  list(col = "Exposure_ART", g1 = "PreART_Entry", g2 = "HUU",
       label = "PreART_Entry_vs_HUU"),
  
  # Longitudinal within HEI
  list(col = "Exposure_ART", g1 = "PostART_Suppressed", g2 = "PreART_Entry",
       label = "Suppressed_vs_PreART"),
  list(col = "Exposure_ART", g1 = "PostART_Unsuppressed", g2 = "PreART_Entry",
       label = "Unsuppressed_vs_PreART"),
  list(col = "Exposure_ART", g1 = "PostART_Suppressed", g2 = "PostART_Unsuppressed",
       label = "Suppressed_vs_Unsuppressed")
)

stat_results <- list()
for (comp in comparisons_list) {
  for (mod in ms_cols) {
    for (cl in unique(TARA_ALL$Annotation)) {
      cells <- TARA_ALL@meta.data[TARA_ALL$Annotation == cl, ]
      a_vals <- cells[cells[[comp$col]] == comp$g1, mod]
      b_vals <- cells[cells[[comp$col]] == comp$g2, mod]
      a_vals <- a_vals[!is.na(a_vals)]
      b_vals <- b_vals[!is.na(b_vals)]
      if (length(a_vals) < 5 || length(b_vals) < 5) next
      
      wt <- suppressWarnings(wilcox.test(a_vals, b_vals))
      d_hat <- (mean(a_vals) - mean(b_vals)) / sd(c(a_vals, b_vals))
      
      stat_results[[length(stat_results) + 1]] <- data.frame(
        Comparison = comp$label, Module = mod, Cluster = cl,
        n_g1 = length(a_vals), n_g2 = length(b_vals),
        mean_g1 = round(mean(a_vals), 4), mean_g2 = round(mean(b_vals), 4),
        cohens_d = round(d_hat, 3), p_value = wt$p.value,
        stringsAsFactors = FALSE)
    }
  }
}

stats_df <- bind_rows(stat_results)
stats_df$p_adj <- p.adjust(stats_df$p_value, method = "BH")
stats_df$sig <- ifelse(stats_df$p_adj < 0.001, "***",
                       ifelse(stats_df$p_adj < 0.01, "**",
                              ifelse(stats_df$p_adj < 0.05, "*", "ns")))
write.csv(stats_df,
          file.path(out_csv, "Module_Score_Statistics_AllComparisons.csv"),
          row.names = FALSE)

# 4e: Viral load correlations
vl_corr <- list()
hei_cells <- TARA_ALL@meta.data[!is.na(TARA_ALL$Viral_Load_num) &
                                  TARA_ALL$Viral_Load_num > 0, ]
for (mod in ms_cols) {
  for (cl in unique(hei_cells$Annotation)) {
    vals <- hei_cells[hei_cells$Annotation == cl, ]
    if (nrow(vals) < 10) next
    ct <- suppressWarnings(cor.test(log10(vals$Viral_Load_num), vals[[mod]],
                                    method = "spearman"))
    vl_corr[[length(vl_corr) + 1]] <- data.frame(
      Module = mod, Cluster = cl,
      spearman_rho = round(ct$estimate, 3),
      p_value = ct$p.value, n = nrow(vals),
      stringsAsFactors = FALSE)
  }
}
vl_corr_df <- bind_rows(vl_corr)
vl_corr_df$p_adj <- p.adjust(vl_corr_df$p_value, method = "BH")
write.csv(vl_corr_df,
          file.path(out_csv, "Module_VL_Correlations.csv"),
          row.names = FALSE)

# 4f: Individual gene scores for Apoptosis module (pro-survival vs pro-death)
apo_genes <- intersect(modules$Apoptosis_Survival, rownames(TARA_ALL[["RNA"]]))
if (length(apo_genes) > 0) {
  apo_expr <- AverageExpression(TARA_ALL, features = apo_genes,
                                group.by = "Annotation", assay = "RNA")$RNA
  write.csv(as.data.frame(apo_expr),
            file.path(out_csv, "Apoptosis_Individual_Genes_PerCluster.csv"))
}

# 4g: RNA vs Protein module score correlation
# For modules with both RNA (MS_) and Protein (PS_) scores, compute
# per-cluster Spearman correlation to test if RNA programme translates
# to surface protein expression
rna_protein_pairs <- list(
  Exhaustion         = c("MS_Exhaustion",         "PS_Exhaustion_Protein"),
  Activation         = c("MS_Activation",         "PS_Activation_Protein"),
  Stemness           = c("MS_Stemness_Naive",     "PS_Stemness_Protein"),
  Tissue_Residency   = c("MS_Tissue_Residency",   "PS_Tissue_Residency_Protein"),
  Senescence         = c("MS_Senescence",         "PS_Senescence_Protein")
)

rna_prot_corr <- list()
for (pair_name in names(rna_protein_pairs)) {
  rna_col <- rna_protein_pairs[[pair_name]][1]
  prot_col <- rna_protein_pairs[[pair_name]][2]
  
  if (!all(c(rna_col, prot_col) %in% colnames(TARA_ALL@meta.data))) next
  
  for (cl in unique(TARA_ALL$Annotation)) {
    cells <- TARA_ALL@meta.data[TARA_ALL$Annotation == cl, ]
    r_vals <- cells[[rna_col]]
    p_vals <- cells[[prot_col]]
    valid <- !is.na(r_vals) & !is.na(p_vals)
    if (sum(valid) < 20) next
    
    ct <- suppressWarnings(cor.test(r_vals[valid], p_vals[valid], method = "spearman"))
    
    rna_prot_corr[[length(rna_prot_corr) + 1]] <- data.frame(
      Module = pair_name, Cluster = cl,
      n_cells = sum(valid),
      spearman_rho = round(ct$estimate, 3),
      p_value = ct$p.value,
      mean_RNA = round(mean(r_vals[valid]), 4),
      mean_Protein = round(mean(p_vals[valid]), 4),
      stringsAsFactors = FALSE)
  }
}

if (length(rna_prot_corr) > 0) {
  rna_prot_df <- bind_rows(rna_prot_corr)
  rna_prot_df$p_adj <- p.adjust(rna_prot_df$p_value, method = "BH")
  rna_prot_df$concordant <- ifelse(
    sign(rna_prot_df$mean_RNA) == sign(rna_prot_df$mean_Protein) &
      rna_prot_df$spearman_rho > 0.1, "Yes",
    ifelse(rna_prot_df$spearman_rho < -0.1, "Discordant", "Weak"))
  write.csv(rna_prot_df,
            file.path(out_csv, "RNA_vs_Protein_Correlation.csv"),
            row.names = FALSE)
  cat("✓ RNA-Protein correlation exported\n")
}

# 4h: Protein module scores by cluster × condition/timepoint
ps_cols <- grep("^PS_", colnames(TARA_ALL@meta.data), value = TRUE)
if (length(ps_cols) > 0) {
  prot_by_cond <- TARA_ALL@meta.data %>%
    select(Annotation, Condition, all_of(ps_cols)) %>%
    group_by(Annotation, Condition) %>%
    summarise(across(all_of(ps_cols), ~ mean(.x, na.rm = TRUE)),
              n_cells = n(), .groups = "drop")
  write.csv(prot_by_cond,
            file.path(out_csv, "Protein_Scores_PerCluster_ByCondition.csv"),
            row.names = FALSE)
  
  prot_by_tp <- TARA_ALL@meta.data %>%
    filter(!is.na(Timepoint_Group)) %>%
    select(Annotation, Timepoint_Group, all_of(ps_cols)) %>%
    group_by(Annotation, Timepoint_Group) %>%
    summarise(across(all_of(ps_cols), ~ mean(.x, na.rm = TRUE)),
              n_cells = n(), .groups = "drop")
  write.csv(prot_by_tp,
            file.path(out_csv, "Protein_Scores_PerCluster_ByTimepoint.csv"),
            row.names = FALSE)
  cat("✓ Protein score summaries exported\n")
}

# Save scored object
qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_ModuleScored.qs2"))

message("\n",
        "══════════════════════════════════════════════════════════════════\n",
        " PHASE 1 CHECKPOINT: CSVs exported for review.\n",
        " Directory: ", out_csv, "\n",
        "\n",
        " Key files:\n",
        "   - Module_Scores_PerCluster_ByCondition.csv\n",
        "   - Module_Scores_PerCluster_ByTimepoint.csv\n",
        "   - Module_Score_Statistics_AllComparisons.csv  (Cohen's d + BH p-adj)\n",
        "   - Module_VL_Correlations.csv                  (Spearman vs viral load)\n",
        "   - Apoptosis_Individual_Genes_PerCluster.csv   (pro- vs anti-apoptosis)\n",
        "   - RNA_vs_Protein_Correlation.csv              (RNA-protein concordance)\n",
        "   - Protein_Scores_PerCluster_ByCondition.csv   (protein surface levels)\n",
        "   - Protein_Scores_PerCluster_ByTimepoint.csv   (protein by ART status)\n",
        "\n",
        " Share these for interpretation, then continue to PHASE 2.\n",
        "══════════════════════════════════════════════════════════════════\n"
)

# ▼▼▼ UNCOMMENT TO STOP HERE ▼▼▼
# stop("PHASE 1 CHECKPOINT — review CSVs before continuing")

# To resume:
# TARA_ALL <- qs_read(file.path(saved_dir, "TARA_ALL_ModuleScored.qs2"))


################################################################################
#
#  PHASE 2: VISUALIZATIONS
#
################################################################################


# =============================================================================
# STEP 5: VOLCANO PLOTS
# =============================================================================

cat("\n══ STEP 5: VOLCANO PLOTS ══\n\n")

highlight_cols <- c(
  "Exhaustion"        = "#A32D2D",  "Stemness"          = "#185FA5",
  "Cytotoxicity"      = "#3B6D11",  "Stress/IFN acute"  = "#BA7517",
  "IFN Memory"        = "#0F6E56",  "Chemokines"        = "#993556",
  "Terminal diff."    = "#534AB7",  "Activation"        = "#D85A30",
  "Glycolysis"        = "#CC6677",  "OXPHOS"            = "#44AA99",
  "Other"             = "grey80"
)

gene_cats <- list(
  "Exhaustion"       = modules$Exhaustion,
  "Stemness"         = modules$Stemness_Naive,
  "Cytotoxicity"     = modules$Cytotoxicity,
  "Stress/IFN acute" = modules$Stress_IFN_Acute,
  "IFN Memory"       = modules$IFN_Memory,
  "Chemokines"       = modules$Inflammatory_Chemokines,
  "Terminal diff."   = modules$Terminal_Differentiation,
  "Activation"       = modules$Activation,
  "Glycolysis"       = modules$Glycolysis_Hypoxia,
  "OXPHOS"           = modules$OXPHOS
)

plot_volcano <- function(de_df, title_text,
                         direction_labels = c("Higher in group 1", "Higher in group 2"),
                         out_path) {
  
  de_df$gene <- rownames(de_df)
  de_df$neg_log10_padj <- -log10(de_df$p_val_adj + 1e-300)
  
  de_df$highlight <- "Other"
  for (cn in names(gene_cats)) {
    de_df$highlight[de_df$gene %in% gene_cats[[cn]]] <- cn
  }
  de_df$highlight <- factor(de_df$highlight, levels = names(highlight_cols))
  
  cat_genes <- unlist(gene_cats)
  de_df$label <- ""
  de_df$label[de_df$gene %in% cat_genes & de_df$p_val_adj < 0.05] <-
    de_df$gene[de_df$gene %in% cat_genes & de_df$p_val_adj < 0.05]
  
  x_vals <- abs(de_df$avg_log2FC[is.finite(de_df$avg_log2FC)])
  x_lim  <- quantile(x_vals, 0.995, na.rm = TRUE) * 1.15
  y_vals <- de_df$neg_log10_padj[is.finite(de_df$neg_log10_padj)]
  y_cap  <- quantile(y_vals, 0.995, na.rm = TRUE) * 1.10
  
  de_df$x_plot <- pmax(pmin(de_df$avg_log2FC, x_lim * 0.98), -x_lim * 0.98)
  de_df$y_plot <- pmin(de_df$neg_log10_padj, y_cap)
  
  p <- ggplot(de_df, aes(x = x_plot, y = y_plot, color = highlight)) +
    geom_point(data = de_df %>% filter(highlight == "Other"),
               size = 1.5, alpha = 0.2) +
    geom_point(data = de_df %>% filter(highlight != "Other"),
               size = 5, alpha = 0.85) +
    geom_label_repel(
      data = de_df %>% filter(label != ""),
      aes(label = label, fill = highlight), color = "white",
      size = 5, fontface = "bold.italic",
      max.overlaps = 40, segment.size = 0.3, segment.color = "grey40",
      min.segment.length = 0.1, box.padding = 0.4,
      label.size = 0.15, show.legend = FALSE,
      force = 2, force_pull = 0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    scale_color_manual(values = highlight_cols, name = "Category", drop = FALSE) +
    scale_fill_manual(values = highlight_cols, guide = "none") +
    scale_x_continuous(limits = c(-x_lim, x_lim), oob = squish) +
    scale_y_continuous(limits = c(0, y_cap), oob = squish) +
    labs(x = expression(log[2]~fold~change),
         y = expression(-log[10]~adjusted~italic(p)),
         title = title_text) +
    theme_cowplot(font_size = 16) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      legend.key.size = unit(1.1, "cm"),
      axis.text = element_text(size = 14),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))
  
  arrow_r <- grobTree(
    linesGrob(x = unit(c(0.52, 0.92), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
              gp = gpar(col = "#52B788", lwd = 7)),
    textGrob(direction_labels[1], x = unit(0.72, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#52B788", fontsize = 18, fontface = "bold")))
  arrow_l <- grobTree(
    linesGrob(x = unit(c(0.48, 0.08), "npc"), y = unit(c(0.5, 0.5), "npc"),
              arrow = arrow(length = unit(0.8, "cm"), type = "closed"),
              gp = gpar(col = "#E76F51", lwd = 7)),
    textGrob(direction_labels[2], x = unit(0.28, "npc"), y = unit(0.0, "npc"),
             gp = gpar(col = "#E76F51", fontsize = 18, fontface = "bold")))
  
  p_final <- p +
    theme(plot.margin = margin(10, 10, 65, 10)) +
    annotation_custom(grob = arrow_r, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.22, ymax = -y_cap * 0.10) +
    annotation_custom(grob = arrow_l, xmin = -Inf, xmax = Inf,
                      ymin = -y_cap * 0.22, ymax = -y_cap * 0.10) +
    coord_cartesian(ylim = c(0, y_cap), xlim = c(-x_lim, x_lim), clip = "off")
  
  ggsave(out_path, p_final, width = 14, height = 11, dpi = 300, bg = "white")
}

dge_comps <- list(
  list(dir = "HEIvsHEU_PreART", title = "HEI vs HEU (Pre-ART)",
       arrows = c("Higher in HEI", "Higher in HEU")),
  list(dir = "HEUvsHUU", title = "HEU vs HUU",
       arrows = c("Higher in HEU", "Higher in HUU")),
  list(dir = "Suppressed_vs_PreART", title = "Suppressed vs Pre-ART",
       arrows = c("Higher in Suppressed", "Higher in Pre-ART")),
  list(dir = "Unsuppressed_vs_PreART", title = "Unsuppressed vs Pre-ART",
       arrows = c("Higher in Unsuppressed", "Higher in Pre-ART")),
  list(dir = "Suppressed_vs_Unsuppressed", title = "Suppressed vs Unsuppressed",
       arrows = c("Higher in Suppressed", "Higher in Unsuppressed"))
)

for (comp in dge_comps) {
  comp_path <- file.path(dge_dir, comp$dir)
  if (!dir.exists(comp_path)) { cat("  Skipping", comp$dir, "\n"); next }
  de_files <- list.files(comp_path, pattern = "\\.csv$", full.names = TRUE)
  vol_sub <- file.path(out_volcano, comp$dir)
  dir.create(vol_sub, recursive = TRUE, showWarnings = FALSE)
  
  for (f in de_files) {
    de <- tryCatch(read.csv(f, row.names = 1), error = function(e) NULL)
    if (is.null(de) || nrow(de) < 10) next
    cl_name <- gsub("\\.csv$", "", basename(f))
    plot_volcano(de, paste0(cl_name, " — ", comp$title),
                 comp$arrows, file.path(vol_sub, paste0(cl_name, "_volcano.png")))
  }
  cat("  ✓", comp$dir, "\n")
}


# =============================================================================
# STEP 6: HEATMAPS
# =============================================================================

cat("\n══ STEP 6: HEATMAPS ══\n\n")

plot_module_heatmap <- function(data, group_col, out_path, title_text) {
  mat <- data %>%
    select(Annotation, all_of(group_col), starts_with("MS_"),
           IFN_Ratio, Metabolic_Ratio) %>%
    pivot_longer(cols = c(starts_with("MS_"), "IFN_Ratio", "Metabolic_Ratio"),
                 names_to = "Module", values_to = "Score") %>%
    group_by(Annotation, .data[[group_col]], Module) %>%
    summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  mat$Module <- gsub("^MS_", "", mat$Module)
  
  p <- ggplot(mat, aes(x = .data[[group_col]], y = Annotation, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~ Module, scales = "free_x", nrow = 1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Mean\nScore") +
    labs(title = title_text, x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 7, face = "bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(out_path, p,
         width = max(24, 2.5 * length(unique(mat$Module))),
         height = max(8, 0.4 * length(unique(mat$Annotation))),
         dpi = 300, bg = "white")
}

ea_meta <- TARA_ALL@meta.data[!is.na(TARA_ALL$Exposure_ART), ]

plot_module_heatmap(ea_meta, "Exposure_ART",
                    file.path(out_heatmap, "Heatmap_AllModules_ByExposureART.png"),
                    "All Module Scores × Cluster × Exposure/ART Group")

cat("✓ Heatmaps done\n\n")


# =============================================================================
# STEP 7: VIOLIN PLOTS
# =============================================================================

cat("══ STEP 7: VIOLINS ══\n\n")

plot_module_violin <- function(obj, module_col, group_col, groups, group_colors,
                               title_text, out_path) {
  md <- obj@meta.data
  md <- md[md[[group_col]] %in% groups & !is.na(md[[module_col]]), ]
  md[[group_col]] <- factor(md[[group_col]], levels = groups)
  if (nrow(md) < 20) return(invisible(NULL))
  
  if (length(groups) == 2) {
    pw <- md %>%
      group_by(Annotation) %>%
      filter(n_distinct(.data[[group_col]]) == 2) %>%
      wilcox_test(as.formula(paste(module_col, "~", group_col))) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj") %>%
      add_xy_position(x = group_col)
  } else { pw <- NULL }
  
  p <- ggplot(md, aes(x = .data[[group_col]], y = .data[[module_col]],
                      fill = .data[[group_col]])) +
    geom_violin(alpha = 0.7, scale = "width", linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9,
                 fill = "white", linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "black") +
    facet_wrap(~ Annotation, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = group_colors, name = "") +
    labs(title = title_text, y = gsub("^MS_", "", module_col), x = "") +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
      legend.position = "top", legend.text = element_text(size = 12),
      panel.spacing = unit(0.8, "lines"),
      panel.grid.major.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA))
  
  if (!is.null(pw) && nrow(pw) > 0) {
    pw_sig <- pw[pw$p.adj.signif != "ns", ]
    if (nrow(pw_sig) > 0) {
      p <- p + stat_pvalue_manual(pw_sig, label = "p.adj.signif",
                                  tip.length = 0.02, size = 3.5, hide.ns = TRUE)
    }
  }
  
  n_cl <- length(unique(md$Annotation))
  ggsave(out_path, p,
         width = min(22, 4.5 * min(5, n_cl)),
         height = max(5, 3.5 * ceiling(n_cl / 5)),
         dpi = 300, bg = "white")
}

tp_cols  <- c("PreART_Entry" = "#E76F51", "PostART_Suppressed" = "#2A9D8F",
              "PostART_Unsuppressed" = "#E9C46A")
sup_cols <- c("PostART_Suppressed" = "#2A9D8F", "PostART_Unsuppressed" = "#E76F51")
entry_heu_cols <- c("PreART_Entry" = "#E76F51", "HEU" = "#264653")
heu_huu_cols   <- c("HEU" = "#264653", "HUU" = "#457B9D")

# Full 5-way palette
ea_cols <- c("HUU" = "#457B9D", "HEU" = "#264653", "PreART_Entry" = "#E76F51",
             "PostART_Suppressed" = "#2A9D8F", "PostART_Unsuppressed" = "#E9C46A")

all_modules <- c(ms_cols[grep("^MS_", ms_cols)], "IFN_Ratio", "Metabolic_Ratio")

for (mod in all_modules) {
  mod_safe <- gsub("^MS_", "", mod)
  vln_sub <- file.path(out_violin, mod_safe)
  dir.create(vln_sub, recursive = TRUE, showWarnings = FALSE)
  
  # 1. PreART_Entry vs HEU (age-matched, effect of HIV infection)
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("PreART_Entry", "HEU"), entry_heu_cols,
                     paste0(mod_safe, " — PreART Entry vs HEU (age-matched)"),
                     file.path(vln_sub, paste0(mod_safe, "_PreART_vs_HEU.png")))
  
  # 2. HEU vs HUU (effect of HIV exposure without infection)
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("HEU", "HUU"), heu_huu_cols,
                     paste0(mod_safe, " — HEU vs HUU"),
                     file.path(vln_sub, paste0(mod_safe, "_HEU_vs_HUU.png")))
  
  # 3. 3-way HEI longitudinal
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"), tp_cols,
                     paste0(mod_safe, " — Pre-ART / Suppressed / Unsuppressed"),
                     file.path(vln_sub, paste0(mod_safe, "_3way_HEI.png")))
  
  # 4. Suppressed vs Unsuppressed (binary)
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("PostART_Suppressed", "PostART_Unsuppressed"), sup_cols,
                     paste0(mod_safe, " — Suppressed vs Unsuppressed"),
                     file.path(vln_sub, paste0(mod_safe, "_Sup_vs_Unsup.png")))
  
  # 5. Full 5-way (all groups)
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("HUU", "HEU", "PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"),
                     ea_cols,
                     paste0(mod_safe, " — All Groups"),
                     file.path(vln_sub, paste0(mod_safe, "_5way_AllGroups.png")))
}

# --- Protein module violins ---
cat("  Protein module violins...\n")

ps_to_plot <- grep("^PS_", colnames(TARA_ALL@meta.data), value = TRUE)
ps_to_plot <- ps_to_plot[!sapply(ps_to_plot, function(x) all(is.na(TARA_ALL@meta.data[[x]])))]

for (mod in ps_to_plot) {
  mod_safe <- gsub("^PS_", "Protein_", mod)
  vln_sub <- file.path(out_violin, mod_safe)
  dir.create(vln_sub, recursive = TRUE, showWarnings = FALSE)
  
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("PreART_Entry", "HEU"), entry_heu_cols,
                     paste0(mod_safe, " (DSB) — PreART Entry vs HEU"),
                     file.path(vln_sub, paste0(mod_safe, "_PreART_vs_HEU.png")))
  
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"), tp_cols,
                     paste0(mod_safe, " (DSB) — Pre-ART / Suppressed / Unsuppressed"),
                     file.path(vln_sub, paste0(mod_safe, "_3way_HEI.png")))
  
  plot_module_violin(TARA_ALL, mod, "Exposure_ART",
                     c("PostART_Suppressed", "PostART_Unsuppressed"), sup_cols,
                     paste0(mod_safe, " (DSB) — Suppressed vs Unsuppressed"),
                     file.path(vln_sub, paste0(mod_safe, "_Sup_vs_Unsup.png")))
}

# --- RNA vs Protein scatter per module ---
cat("  RNA vs Protein scatter plots...\n")

rna_prot_scatter_dir <- file.path(out_violin, "RNA_vs_Protein")
dir.create(rna_prot_scatter_dir, recursive = TRUE, showWarnings = FALSE)

for (pair_name in names(rna_protein_pairs)) {
  rna_col <- rna_protein_pairs[[pair_name]][1]
  prot_col <- rna_protein_pairs[[pair_name]][2]
  if (!all(c(rna_col, prot_col) %in% colnames(TARA_ALL@meta.data))) next
  
  # Per-cluster mean scatter colored by Exposure_ART
  scatter_df <- TARA_ALL@meta.data %>%
    filter(!is.na(Exposure_ART)) %>%
    group_by(Annotation, Exposure_ART) %>%
    summarise(RNA = mean(.data[[rna_col]], na.rm = TRUE),
              Protein = mean(.data[[prot_col]], na.rm = TRUE),
              n = n(), .groups = "drop")
  
  p_rp <- ggplot(scatter_df, aes(x = RNA, y = Protein, color = Exposure_ART)) +
    geom_point(aes(size = n), alpha = 0.8) +
    geom_text_repel(aes(label = Annotation), size = 3, max.overlaps = 25,
                    show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, linetype = "dashed") +
    scale_color_manual(values = ea_cols, name = "Group") +
    scale_size_continuous(range = c(2, 8), name = "n cells") +
    labs(x = paste0(pair_name, " RNA Score"),
         y = paste0(pair_name, " Protein Score (DSB)"),
         title = paste0(pair_name, " — RNA vs Protein per Cluster")) +
    theme_cowplot(font_size = 13) +
    theme(legend.position = "top",
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(rna_prot_scatter_dir,
                   paste0(pair_name, "_RNA_vs_Protein.png")),
         p_rp, width = 12, height = 10, dpi = 300, bg = "white")
}

cat("✓ Violins done\n\n")


# =============================================================================
# STEP 8: RATIO ANALYSES (IFN + Metabolic)
# =============================================================================

cat("══ STEP 8: RATIO ANALYSIS ══\n\n")

plot_ratio_scatter <- function(data, x_col, y_col, group_col, group_colors,
                               x_lab, y_lab, title_text, out_path,
                               quadrant_labels = NULL) {
  summary_df <- data %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(Annotation, .data[[group_col]]) %>%
    summarise(x = mean(.data[[x_col]], na.rm = TRUE),
              y = mean(.data[[y_col]], na.rm = TRUE),
              n = n(), .groups = "drop")
  
  p <- ggplot(summary_df, aes(x = x, y = y, color = .data[[group_col]])) +
    geom_point(aes(size = n), alpha = 0.8) +
    geom_text_repel(aes(label = Annotation), size = 3, max.overlaps = 25,
                    show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    scale_color_manual(values = group_colors, name = "") +
    scale_size_continuous(range = c(2, 10), name = "n cells") +
    labs(x = x_lab, y = y_lab, title = title_text) +
    theme_cowplot(font_size = 14) +
    theme(legend.position = "top",
          plot.background = element_rect(fill = "white", color = NA))
  
  if (!is.null(quadrant_labels)) {
    p <- p +
      annotate("text", x = Inf, y = Inf, label = quadrant_labels[1],
               hjust = 1.1, vjust = 1.1, fontface = "italic",
               color = "#2A9D8F", size = 4.5) +
      annotate("text", x = Inf, y = -Inf, label = quadrant_labels[2],
               hjust = 1.1, vjust = -0.1, fontface = "italic",
               color = "#E76F51", size = 4.5)
  }
  
  ggsave(out_path, p, width = 13, height = 11, dpi = 300, bg = "white")
}

# IFN ratio scatter — by Exposure_ART (5 groups)
plot_ratio_scatter(
  TARA_ALL@meta.data, "MS_Stress_IFN_Acute", "MS_IFN_Memory", "Exposure_ART",
  ea_cols, "Mean Stress/IFN Acute", "Mean IFN Memory",
  "IFN Memory vs Stress Balance — All Groups",
  file.path(out_ratio, "IFN_Memory_vs_Stress_AllGroups.png"),
  c("IFN-primed (good)", "Stressed (bad)")
)

# Metabolic ratio scatter — by Exposure_ART
plot_ratio_scatter(
  TARA_ALL@meta.data, "MS_Glycolysis_Hypoxia", "MS_OXPHOS", "Exposure_ART",
  ea_cols, "Mean Glycolysis/Hypoxia", "Mean OXPHOS",
  "Metabolic Fitness — All Groups",
  file.path(out_ratio, "OXPHOS_vs_Glycolysis_AllGroups.png"),
  c("Metabolically fit", "Glycolytic trap")
)

# Export ratio CSV
ratio_sum <- TARA_ALL@meta.data %>%
  filter(!is.na(Exposure_ART)) %>%
  group_by(Annotation, Exposure_ART) %>%
  summarise(
    IFN_Memory = mean(MS_IFN_Memory, na.rm = TRUE),
    Stress = mean(MS_Stress_IFN_Acute, na.rm = TRUE),
    IFN_Ratio = mean(IFN_Ratio, na.rm = TRUE),
    OXPHOS = mean(MS_OXPHOS, na.rm = TRUE),
    Glycolysis = mean(MS_Glycolysis_Hypoxia, na.rm = TRUE),
    Metabolic_Ratio = mean(Metabolic_Ratio, na.rm = TRUE),
    n = n(), .groups = "drop")
write.csv(ratio_sum,
          file.path(out_csv, "Ratio_Summary_ByExposureART.csv"),
          row.names = FALSE)

cat("✓ Ratio analysis done\n\n")


# =============================================================================
# STEP 9: VL CORRELATION BUBBLE PLOT
# =============================================================================

cat("══ STEP 9: VL CORRELATIONS ══\n\n")

vl_corr_df <- read.csv(file.path(out_csv, "Module_VL_Correlations.csv"))
vl_corr_df$Module <- gsub("^MS_", "", vl_corr_df$Module)
vl_corr_df$sig <- vl_corr_df$p_adj < 0.05

p_bubble <- ggplot(vl_corr_df,
                   aes(x = Module, y = Cluster, fill = spearman_rho,
                       size = -log10(p_adj + 1e-10))) +
  geom_point(shape = 21, stroke = 0.3) +
  geom_point(data = vl_corr_df %>% filter(sig),
             shape = 21, stroke = 1.2, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  scale_size_continuous(range = c(1, 8), name = expression(-log[10]~p[adj])) +
  labs(title = "Module Score vs Viral Load", x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid = element_line(color = "grey95"),
    plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(out_corr, "VL_Correlation_Bubble_AllModules.png"),
       p_bubble,
       width = max(14, 1.2 * length(unique(vl_corr_df$Module))),
       height = max(8, 0.4 * length(unique(vl_corr_df$Cluster))),
       dpi = 300, bg = "white")

cat("✓ VL correlations done\n\n")


# =============================================================================
# STEP 10: IFNAR RECEPTOR PROFILING
# =============================================================================

cat("══ STEP 10: IFNAR ══\n\n")

DefaultAssay(TARA_ALL) <- "RNA"

for (gene in c("IFNAR1", "IFNAR2")) {
  if (!gene %in% rownames(TARA_ALL[["RNA"]])) next
  
  expr <- FetchData(TARA_ALL, vars = c(gene, "Annotation", "Exposure_ART"))
  expr <- expr[!is.na(expr$Exposure_ART), ]
  expr$Exposure_ART <- factor(expr$Exposure_ART,
                              levels = c("HUU", "HEU", "PreART_Entry", "PostART_Suppressed", "PostART_Unsuppressed"))
  
  p <- ggplot(expr, aes(x = Annotation, y = .data[[gene]],
                        fill = Exposure_ART)) +
    geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3, trim = TRUE,
                position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
                 fill = "white", linewidth = 0.3,
                 position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = ea_cols, name = "Group") +
    labs(title = paste0(gene, " Expression × Cluster × Group"),
         y = paste0(gene, " (RNA)"), x = "") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "top",
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(out_receptor, paste0(gene, "_Violin_AllGroups.png")),
         p, width = 22, height = 8, dpi = 300, bg = "white")
}

cat("✓ IFNAR done\n\n")


# =============================================================================
# STEP 11: SAVE
# =============================================================================

#qs_save(TARA_ALL, file.path(saved_dir, "TARA_ALL_ModuleScored.qs2"))

message("\n",
        "══════════════════════════════════════════════════════════════\n",
        " MODULE SCORING ANALYSIS COMPLETE\n",
        "══════════════════════════════════════════════════════════════\n",
        " Output: ", out_base, "\n",
        "\n",
        " CSV_Checkpoints/  — RNA + protein scores, stats, VL corr, RNA-protein concordance\n",
        " Volcanos/         — per-cluster per-comparison (10 gene categories)\n",
        " Violins/          — 17 RNA + 6 protein modules × 3 comparisons each\n",
        "                     + RNA vs Protein scatter per module\n",
        " Heatmaps/         — cluster × condition overview\n",
        " Ratio_Analysis/   — IFN + Metabolic balance scatters\n",
        " Correlations/     — VL correlation bubble\n",
        " IFNAR_Receptor/   — IFNAR1/2 profiling\n",
        "\n",
        " Module_Definitions_and_Interpretation_Guide.txt\n",
        "══════════════════════════════════════════════════════════════\n"
)