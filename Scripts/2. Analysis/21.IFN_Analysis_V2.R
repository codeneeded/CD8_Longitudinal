library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratExtend)   # DotPlot2 + VlnPlot2
library(stringr)
# --------------------------
# Paths + load object
# --------------------------
base_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal"
out_dir  <- file.path(base_dir, "IFN_Cytokine_Analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_in <- file.path(base_dir, "saved_R_data", "TARA_ALL_post_annotation.rds")
TARA_ALL <- readRDS(rds_in)
DefaultAssay(TARA_ALL) <- "RNA"

# --------------------------
# Make sure Viral_Load is numeric
# --------------------------
TARA_ALL$Viral_Load <- suppressWarnings(as.numeric(as.character(TARA_ALL$Viral_Load)))

# --------------------------
# Viral_Load_Category (High/Low for HEI; pass through HEU/HUU)
# --------------------------
TARA_ALL@meta.data <- TARA_ALL@meta.data %>%
  mutate(
    Viral_Load_Category = case_when(
      Condition == "HEI" & !is.na(Viral_Load) & Viral_Load >= 100000 ~ "High",
      Condition == "HEI" & !is.na(Viral_Load) & Viral_Load < 100000  ~ "Low",
      Condition == "HEU" ~ "HEU",
      Condition == "HUU" ~ "HUU",
      TRUE ~ NA_character_
    )
  )

# --------------------------
# Timepoint_Group (HEI PreART / PostART suppressed / PostART unsuppressed; HEU/HUU passthrough)
# --------------------------
TARA_ALL@meta.data <- TARA_ALL@meta.data %>%
  mutate(
    Timepoint_Group = case_when(
      Condition == "HEI" & !is.na(Age) & Age <= 2 ~ "PreART_Entry",
      Condition == "HEI" & !is.na(Viral_Load) & Viral_Load < 200 ~ "PostART_Suppressed",
      Condition == "HEI" & !is.na(Viral_Load) & Viral_Load >= 200 ~ "PostART_Unsuppressed",
      Condition %in% c("HEU", "HUU") ~ Condition,
      TRUE ~ NA_character_
    )
  )

# --------------------------
# (Optional) A combined HEI-only variable you asked for:
# HEI_PreART_Entry / HEI_PostART_Suppressed / HEI_PostART_Unsuppressed
# plus HEU/HUU passthrough
# --------------------------
TARA_ALL@meta.data <- TARA_ALL@meta.data %>%
  mutate(
    HEI_Timepoint_Status = case_when(
      Condition == "HEI" & Timepoint_Group == "PreART_Entry" ~ "HEI_PreART_Entry",
      Condition == "HEI" & Timepoint_Group == "PostART_Suppressed" ~ "HEI_PostART_Suppressed",
      Condition == "HEI" & Timepoint_Group == "PostART_Unsuppressed" ~ "HEI_PostART_Unsuppressed",
      Condition == "HEU" ~ "HEU",
      Condition == "HUU" ~ "HUU",
      TRUE ~ NA_character_
    )
  )

# --------------------------
# Set Idents to your cell types
# --------------------------
Idents(TARA_ALL) <- "Manual_Annotation"
# ======================================================
# STEP 2.5 — Subset desired immune lineages
# ======================================================

# View available clusters (optional sanity check)
levels(as.factor(TARA_ALL$Manual_Annotation))

# --------------------------
# Define keep / drop patterns
# --------------------------
keep_pat <- "(CD4|CD8|CTL|NK)"
drop_pat <- "(DN T|Gamma Delta|Monocyte|B cell|Plasmablast|pDC|APC)"

clusters_keep <- levels(factor(TARA_ALL$Manual_Annotation))

# keep CD4/CD8/NK/CTL
clusters_keep <- clusters_keep[
  grepl(keep_pat, clusters_keep, ignore.case = TRUE)
]

# remove unwanted groups
clusters_keep <- clusters_keep[
  !grepl(drop_pat, clusters_keep, ignore.case = TRUE)
]

# OPTIONAL: remove TRDV1+ CTL-like (gamma-delta–ish)
clusters_keep <- setdiff(clusters_keep, "9: TRDV1+ CTL-like")

# --------------------------
# Subset object
# --------------------------
TARA_ALL <- subset(
  TARA_ALL,
  subset = Manual_Annotation %in% clusters_keep
)

TARA_ALL <- droplevels(TARA_ALL)

# Reset identities
Idents(TARA_ALL) <- "Manual_Annotation"

# --------------------------
# Check what remains
# --------------------------
levels(as.factor(TARA_ALL$Manual_Annotation))
table(TARA_ALL$Manual_Annotation)

# ======================================================
# STEP 3 — IFNAR1 expression by cell type + contrasts
# ======================================================

# --------------------------
# Output folder for Step 3
# --------------------------
step3_dir <- file.path(out_dir, "Step3_IFNAR1")
dir.create(step3_dir, recursive = TRUE, showWarnings = FALSE)

# Make sure Idents are set
Idents(TARA_ALL) <- "Manual_Annotation"

# # --------------------------
# STEP 3A: DotPlot2 (NO split) to identify IFNAR1+ cell types
# --------------------------
step3_dir <- file.path(out_dir, "Step3_IFNAR1")
dir.create(step3_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(TARA_ALL) <- "RNA"
Idents(TARA_ALL) <- "Manual_Annotation"


p_dot_ifnar12 <- SeuratExtend::DotPlot2(
  TARA_ALL,
  features = c("IFNAR1", "IFNAR2"),
  group.by = "Manual_Annotation",
  flip = FALSE
) +
  RotatedAxis() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Cell Type",
    y = "Gene"
  ) +
  ggtitle("IFNAR1/IFNAR2 expression by cell type")

ggsave(
  file.path(step3_dir, "DotPlot2_IFNAR1_IFNAR2_byCellType_noSplit.png"),
  p_dot_ifnar12,
  dpi = 600,
  width = 12,
  height = 8,
  bg = "white"
)

# ======================================================
# STEP 3B–3D — IFNAR1 + IFNAR2 VlnPlot2 contrasts (minimal styling)
# ======================================================

DefaultAssay(TARA_ALL) <- "RNA"
Idents(TARA_ALL) <- "Manual_Annotation"

# ---- minimal helper: same plot style you tested ----
plot_ifnar_vln <- function(obj, gene, split_col, groups, out_png, title) {

  md <- obj@meta.data
  keep_cells <- rownames(md)[ md[[split_col]] %in% groups ]

  obj2 <- subset(obj, cells = keep_cells)
  obj2 <- droplevels(obj2)

  DefaultAssay(obj2) <- "RNA"
  Idents(obj2) <- "Manual_Annotation"

  # enforce order
  obj2[[split_col]] <- factor(obj2[[split_col]][,1], levels = groups)

  p <- SeuratExtend::VlnPlot2(
    obj2,
    features    = gene,
    split.by    = split_col,
    stat.method = "wilcox.test",
    show.mean   = TRUE,
    pt.alpha    = 0
  ) +
    ggtitle(title)

  ggsave(out_png, p, dpi = 600, width = 16, height = 7, bg = "white")
  invisible(p)
}

# genes to run
for (gene in c("IFNAR1","IFNAR2")) {

  # --------------------------
  # 3B) Condition contrasts
  # --------------------------
  plot_ifnar_vln(
    TARA_ALL, gene,
    "Condition", c("HEI","HEU"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_HEI_vs_HEU.png")),
    paste0(gene, " — HEI vs HEU")
  )

  plot_ifnar_vln(
    TARA_ALL, gene,
    "Condition", c("HEU","HUU"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_HEU_vs_HUU.png")),
    paste0(gene, " — HEU vs HUU")
  )

  plot_ifnar_vln(
    TARA_ALL, gene,
    "Condition", c("HEI","HUU"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_HEI_vs_HUU.png")),
    paste0(gene, " — HEI vs HUU")
  )

  # --------------------------
  # 3C) High vs Low VL (HEI only)
  # --------------------------
  hei_only <- subset(TARA_ALL, subset = Condition == "HEI")

  plot_ifnar_vln(
    hei_only, gene,
    "Viral_Load_Category", c("High","Low"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_High_vs_LowVL_HEIonly.png")),
    paste0(gene, " — High vs Low VL (HEI only)")
  )

  # --------------------------
  # 3D) Timepoint contrasts (HEI only)
  # --------------------------
  plot_ifnar_vln(
    hei_only, gene,
    "Timepoint_Group", c("PreART_Entry","PostART_Suppressed"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_PreART_vs_PostSupp.png")),
    paste0(gene, " — PreART vs PostART Suppressed")
  )

  plot_ifnar_vln(
    hei_only, gene,
    "Timepoint_Group", c("PreART_Entry","PostART_Unsuppressed"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_PreART_vs_PostUnsupp.png")),
    paste0(gene, " — PreART vs PostART Unsuppressed")
  )

  plot_ifnar_vln(
    hei_only, gene,
    "Timepoint_Group", c("PostART_Suppressed","PostART_Unsuppressed"),
    file.path(step3_dir, paste0("VlnPlot2_", gene, "_PostSupp_vs_PostUnsupp.png")),
    paste0(gene, " — Post Suppressed vs Unsuppressed")
  )
}

message("Step 3 IFNAR1 + IFNAR2 VlnPlot2 contrasts complete. Outputs in: ", step3_dir)

# ======================================================
# STEP 4 — Module scoring + VlnPlot2 for all contrasts
# (minimal styling; same look as your working VlnPlot2 call)
# Outputs -> /home/akshay-iyer/Documents/CD8_Longitudinal/IFN_Cytokine_Analysis/Step4_ModuleScores
# ======================================================

step4_dir <- file.path(out_dir, "Step4_ModuleScores")
dir.create(step4_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(TARA_ALL) <- "RNA"
Idents(TARA_ALL) <- "Manual_Annotation"

# --------------------------
# Modules
# --------------------------
modules <- list(
  Activation_HLADRA_CD38 = c("HLA-DRA", "CD38"),
  
  Primary_IFN_I_module = c("ISG15","IFIT1","IFIT2","IFIT3","IFI6","IFI27","MX1","MX2","OAS1","OAS2","OAS3","BST2"),
  IFN_Effector_module  = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","MX2","OAS1","OAS2","OAS3"),
  IFN_signalling_module = c("STAT1","STAT2","IRF7","IRF9","IRF1","IRF3"),
  IFN_negative_regulation_module = c("USP18","SOCS1","SOCS3"),
  CD8_NK_interface_genes = c("GZMB","PRF1","NKG7","GNLY","TYROBP"),
  Immune_stress_transcription_factors = c("JUN","FOS","FOSL2","ATF3","BATF"),
  IFN_linked_transcription_factors = c("STAT1","IRF1")
)

cytokine_module <- c(
  "IL1B","IL1A","TNF","IL6","CXCL8","IL18",
  "CXCL10","CXCL9",
  "CCL2","CCL3","CCL4","CCL5",
  "CSF2","IFNG",
  "IL12A","IL12B","IL23A",
  "IL1RN"
)

# Add cytokines module into same list
modules[["Inflammatory_Cytokines_Module"]] <- cytokine_module

# --------------------------
# Score function (stores as MS_<ModuleName>)
# --------------------------
score_module <- function(obj, genes, module_name) {
  genes_present <- intersect(genes, rownames(obj))
  if (length(genes_present) == 0) {
    obj[[paste0("MS_", module_name)]] <- NA_real_
    return(obj)
  }
  
  obj <- AddModuleScore(
    object   = obj,
    features = list(genes_present),
    name     = paste0("MS_", module_name),
    assay    = DefaultAssay(obj)
  )
  
  # AddModuleScore makes "MS_<name>1" — rename to "MS_<name>"
  tmp <- paste0("MS_", module_name, "1")
  obj[[paste0("MS_", module_name)]] <- obj[[tmp]]
  obj[[tmp]] <- NULL
  
  obj
}

# --------------------------
# Minimal VlnPlot2 helper (same look as your working call)
# --------------------------
plot_ms_vln <- function(obj, feature, split_col, groups, out_png, title) {
  md <- obj@meta.data
  keep_cells <- rownames(md)[ md[[split_col]] %in% groups ]
  
  obj2 <- subset(obj, cells = keep_cells)
  obj2 <- droplevels(obj2)
  
  DefaultAssay(obj2) <- "RNA"
  Idents(obj2) <- "Manual_Annotation"
  
  # enforce order
  obj2[[split_col]] <- factor(obj2[[split_col]][,1], levels = groups)
  
  p <- SeuratExtend::VlnPlot2(
    obj2,
    features    = feature,
    split.by    = split_col,
    stat.method = "wilcox.test",
    show.mean   = TRUE,
    pt.alpha    = 0
  ) +
    ggtitle(title)
  
  ggsave(out_png, p, dpi = 600, width = 16, height = 7, bg = "white")
  invisible(p)
}

# --------------------------
# Contrasts to run
# --------------------------
condition_contrasts <- list(
  c("HEI","HEU"),
  c("HEU","HUU"),
  c("HEI","HUU")
)

vl_contrast <- c("High","Low")

timepoint_contrasts <- list(
  c("PreART_Entry","PostART_Suppressed"),
  c("PreART_Entry","PostART_Unsuppressed"),
  c("PostART_Suppressed","PostART_Unsuppressed")
)

# --------------------------
# Score + plot loop
# --------------------------
hei_only <- subset(TARA_ALL, subset = Condition == "HEI")
hei_only <- droplevels(hei_only)

for (m in names(modules)) {
  
  message("Scoring module: ", m)
  
  # Score on full object
  TARA_ALL <- score_module(TARA_ALL, modules[[m]], m)
  
  ms_col <- paste0("MS_", m)
  
  # Make a module folder
  mod_dir <- file.path(step4_dir, m)
  dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 4A) Condition contrasts on full object
  for (cc in condition_contrasts) {
    a <- cc[1]; b <- cc[2]
    plot_ms_vln(
      obj = TARA_ALL,
      feature = ms_col,
      split_col = "Condition",
      groups = c(a,b),
      out_png = file.path(mod_dir, paste0(ms_col, "__", a, "_vs_", b, "__Condition.png")),
      title = paste0(ms_col, " — ", a, " vs ", b, " (Condition)")
    )
  }
  
  # Score HEI-only object too (so ms exists there)
  hei_only <- score_module(hei_only, modules[[m]], m)
  
  # 4B) High vs Low VL (HEI only)
  plot_ms_vln(
    obj = hei_only,
    feature = ms_col,
    split_col = "Viral_Load_Category",
    groups = vl_contrast,
    out_png = file.path(mod_dir, paste0(ms_col, "__High_vs_Low__Viral_Load_Category_HEIonly.png")),
    title = paste0(ms_col, " — High vs Low VL (HEI only)")
  )
  
  # 4C) Timepoint contrasts (HEI only)
  for (tc in timepoint_contrasts) {
    a <- tc[1]; b <- tc[2]
    plot_ms_vln(
      obj = hei_only,
      feature = ms_col,
      split_col = "Timepoint_Group",
      groups = c(a,b),
      out_png = file.path(mod_dir, paste0(ms_col, "__", a, "_vs_", b, "__Timepoint_Group_HEIonly.png")),
      title = paste0(ms_col, " — ", a, " vs ", b, " (HEI only, Timepoint_Group)")
    )
  }
}


message("Step 4 done. Outputs in: ", step4_dir)
