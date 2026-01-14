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
library(clustree)
library(scRepertoire)
library(monocle3)
library(SeuratWrappers)
library(ggrepel)

##############################3
base_dir   <- "~/Documents/CD8_Longitudinal"
saved_dir  <- file.path(base_dir, "saved_R_data")
rds_in     <- file.path(saved_dir, "tara_cdnk_annotated.rds")
tara_cdnk <- readRDS(rds_in)
levels(as.factor(tara_cdnk$Timepoint_Group))
levels(as.factor(tara_cdnk$orig.ident))
levels(as.factor(tara_cdnk$T_Cell_A))
levels(as.factor(tara_cdnk$clonalFrequency))

## ---------------------------
## Step 1 — Define clonal status + exclude no-TCR cells
## ---------------------------


# tara_cdnk already exists in your environment
# 1) Keep only cells with a valid clonalFrequency:
#    - exclude NA (no TCR)
#    - exclude anything non-numeric or <=0 if it exists
tara_cdnk$clonalFrequency_num <- suppressWarnings(as.numeric(tara_cdnk$clonalFrequency))

tara_cdnk_tcr <- subset(
  tara_cdnk,
  subset = !is.na(clonalFrequency_num) & clonalFrequency_num > 0
)
tara_cdnk_tcr$PID <- sub("_.*", "", tara_cdnk_tcr$orig.ident)

# 2) Label clonally expanding vs non-expanding
tara_cdnk_tcr$ClonalStatus <- ifelse(
  tara_cdnk_tcr$clonalFrequency_num > 1,
  "Clonally_Expanding",
  "Non_Expanding"  # == 1
)
tara_cdnk_tcr$ClonalStatus <- factor(
  tara_cdnk_tcr$ClonalStatus,
  levels = c("Non_Expanding", "Clonally_Expanding")
)

# 3) Quick sanity checks
table(tara_cdnk$clonalFrequency, useNA = "ifany")
table(tara_cdnk_tcr$ClonalStatus)
# ----------------------------
# Output directory (question-first)
# ----------------------------
base_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/VDJ"

# One folder for this question
q_dir <- file.path(base_dir, "Q1_Clonal_vs_NonClonal_Expansion")
dir.create(q_dir, recursive = TRUE, showWarnings = FALSE)

# 1) UMAP: overall
p_umap_overall <- DimPlot2(
  tara_cdnk_tcr,
  group.by  = "ClonalStatus",
  reduction = "wnn.umap",
  pt.size   = 0.1
)

ggsave(
  filename = file.path(q_dir, "UMAP__ClonalStatus__Overall.png"),
  plot     = p_umap_overall,
  width    = 7,
  height   = 6,
  dpi      = 300,
  bg='white'
)

# 2) UMAP: split by timepoint
p_umap_split_tp <- DimPlot2(
  tara_cdnk_tcr,
  group.by  = "ClonalStatus",
  split.by  = "Timepoint_Group",
  reduction = "wnn.umap",
  pt.size   = 0.1
)

ggsave(
  filename = file.path(q_dir, "UMAP__ClonalStatus__SplitBy_Timepoint_Group.png"),
  plot     = p_umap_split_tp,
  width    = 12,
  height   = 6,
  dpi      = 300,
  bg='white'
)

# 1) UMAP colored by cluster (whatever your current Idents are)
p_cluster <- DimPlot2(
  tara_cdnk_tcr,
  group.by  = "T_Cell_A",
  reduction = "wnn.umap",
  box=T,
  pt.size   = 0.1,
  label     = TRUE
)

ggsave(
  filename = file.path(q_dir, "UMAP__Clusters_TCR.png"),
  plot     = p_cluster,
  width    = 15,
  height   = 9,
  dpi      = 300,
  bg='white'
)

# 2) ClonalStatus split by cluster
p_clonal_split_cluster <- DimPlot2(
  tara_cdnk_tcr,
  group.by  = "ClonalStatus",
  split.by  = "T_Cell_A",
  reduction = "wnn.umap",
  pt.size   = 0.1
)

ggsave(
  filename = file.path(q_dir, "UMAP__ClonalStatus__SplitBy_Cluster.png"),
  plot     = p_clonal_split_cluster,
  width    = 14,
  height   = 10,
  dpi      = 300,
  bg='white'
)
DimPlot2(
  tara_cdnk_tcr,
  group.by  = "ClonalStatus",
  split.by  = "orig.ident",
  reduction = "wnn.umap",
  pt.size   = 0.1
)

ggsave(
  filename = file.path(q_dir, "UMAP__ClonalStatus__SplitBy_Sample.png"),
  width    = 14,
  height   = 10,
  dpi      = 300,
  bg='white'
)
tara_cdnk_tcr$T_Cell_A
DimPlot2(
  tara_cdnk_tcr,
  group.by  = "T_Cell_A",
  split.by  = "orig.ident",
  reduction = "wnn.umap",
  pt.size   = 0.1
)
ggsave(
  filename = file.path(q_dir, "UMAP__Annotated__SplitBy_Sample.png"),
  width    = 14,
  height   = 10,
  dpi      = 300,
  bg='white'
)
############################################################
# 0) Paths + assay
############################################################
base_dir <- "/home/akshay-iyer/Documents/CD8_Longitudinal/T_Cell_Subsets/VDJ"
q_dir    <- file.path(base_dir, "Q1_Clonal_vs_NonClonal_Expansion")
out_dir  <- file.path(q_dir, "DGE_MAST")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(tara_cdnk_tcr) <- "RNA"

############################################################
# 1) Minimal metadata prep
############################################################
tara_cdnk_tcr$PID <- sub("_.*", "", tara_cdnk_tcr$orig.ident)

tara_cdnk_tcr$ClonalStatus    <- factor(tara_cdnk_tcr$ClonalStatus)
tara_cdnk_tcr$Timepoint_Group <- factor(tara_cdnk_tcr$Timepoint_Group) # used as ART status grouping
tara_cdnk_tcr$T_Cell_A        <- factor(tara_cdnk_tcr$T_Cell_A)

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

make_block_dirs <- function(parent, block_name) {
  block_dir <- file.path(parent, block_name)
  dge_dir   <- file.path(block_dir, "DGE")
  plt_dir   <- file.path(block_dir, "Plots")
  dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plt_dir, recursive = TRUE, showWarnings = FALSE)
  list(block_dir = block_dir, dge_dir = dge_dir, plt_dir = plt_dir)
}

############################################################
# 2) Volcano plot (your style, better labels)
############################################################
plot_volcano <- function(de_table,
                         out_png,
                         padj_thr = 0.05,
                         lfc_thr  = 0.1,
                         top_n    = 10,
                         title    = NULL,
                         xlab     = "avg_log2FC") {
  
  df <- de_table
  if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
  
  p_col  <- c("p_val_adj","p_val","PValue","padj")[c("p_val_adj","p_val","PValue","padj") %in% colnames(df)][1]
  fc_col <- c("avg_log2FC","log2FC","avg_logFC")[c("avg_log2FC","log2FC","avg_logFC") %in% colnames(df)][1]
  if (is.na(p_col) || is.na(fc_col)) stop("Could not find p-value and/or logFC columns in de_table.")
  
  df$neglog10_padj <- -log10(df[[p_col]] + 1e-300)
  
  df$sig <- ifelse(df[[p_col]] < padj_thr & df[[fc_col]] >  lfc_thr, "Higher in Expanding",
                   ifelse(df[[p_col]] < padj_thr & df[[fc_col]] < -lfc_thr, "Higher in Non-expanding", "NS"))
  
  lab_up <- subset(df, df[[p_col]] < padj_thr & df[[fc_col]] >  lfc_thr)
  lab_dn <- subset(df, df[[p_col]] < padj_thr & df[[fc_col]] < -lfc_thr)
  
  lab_up <- lab_up[order(-lab_up$neglog10_padj), ][seq_len(min(nrow(lab_up), top_n)), , drop = FALSE]
  lab_dn <- lab_dn[order(-lab_dn$neglog10_padj), ][seq_len(min(nrow(lab_dn), top_n)), , drop = FALSE]
  lab_genes <- unique(c(lab_up$gene, lab_dn$gene))
  
  if (is.null(title)) title <- "Volcano"
  
  p <- ggplot(df, aes(x = .data[[fc_col]], y = neglog10_padj, color = sig)) +
    geom_point(size = 1.5) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
    geom_hline(yintercept = -log10(padj_thr), linetype = 2) +
    scale_color_manual(values = c(
      "Higher in Expanding"      = "#D55E00",
      "Higher in Non-expanding"  = "#0072B2",
      "NS"                       = "grey"
    )) +
    labs(title = title, x = xlab, y = "-log10(adj p)") +
    theme_bw() +
    geom_text_repel(
      data = subset(df, gene %in% lab_genes),
      aes(label = gene),
      max.overlaps = Inf,
      size = 3
    )
  
  ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  invisible(p)
}

############################################################
# 3) Safe FindMarkers wrapper (never crashes on missing groups)
############################################################
run_FindMarkers_safe <- function(seu, group_by, ident1, ident2,
                                 latent_vars,
                                 label,
                                 logfc.threshold = 0.1,
                                 min.pct = 0.05,
                                 min_cells_per_group = 30) {
  
  md <- seu[[]]
  if (!(group_by %in% colnames(md))) {
    return(list(ok = FALSE, label = label, reason = paste0("Missing group.by column: ", group_by)))
  }
  
  tab <- table(md[[group_by]])
  n1 <- if (ident1 %in% names(tab)) unname(tab[[ident1]]) else 0
  n2 <- if (ident2 %in% names(tab)) unname(tab[[ident2]]) else 0
  
  if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
    return(list(ok = FALSE, label = label,
                reason = paste0("Too few cells: ", ident1, "=", n1, ", ", ident2, "=", n2),
                n1 = n1, n2 = n2))
  }
  
  latent_vars <- latent_vars[latent_vars %in% colnames(md)]
  
  res <- FindMarkers(
    seu,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    group.by = group_by,
    latent.vars = latent_vars,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct
  )
  
  res$gene <- rownames(res)
  list(ok = TRUE, label = label, res = res, n1 = n1, n2 = n2, latent_used = latent_vars)
}

write_res_and_plot <- function(out, dirs, base_name, title, xlab) {
  if (!out$ok) return(invisible(NULL))
  
  csv_path <- file.path(dirs$dge_dir, paste0(base_name, ".csv"))
  write.csv(out$res, csv_path, row.names = FALSE)
  
  png_path <- file.path(dirs$plt_dir, paste0(base_name, "__Volcano.png"))
  plot_volcano(out$res, png_path, title = title, xlab = xlab)
}

############################################################
# A) OVERALL: Clonal vs Non-clonal
############################################################
dirs_A <- make_block_dirs(out_dir, "Overall__Clonal_vs_NonClonal")

manifest_A <- list()

label_A <- "Clonally_Expanding_vs_Non_Expanding"
out_A <- run_FindMarkers_safe(
  seu = tara_cdnk_tcr,
  group_by = "ClonalStatus",
  ident1 = "Clonally_Expanding",
  ident2 = "Non_Expanding",
  latent_vars = c("nCount_RNA", "Timepoint_Group", "PID"),
  label = label_A
)

manifest_A[[1]] <- data.frame(
  Analysis = "Overall",
  Comparison = label_A,
  OK = out_A$ok,
  Reason = ifelse(out_A$ok, "", out_A$reason),
  n_Ident1 = ifelse(is.null(out_A$n1), NA, out_A$n1),
  n_Ident2 = ifelse(is.null(out_A$n2), NA, out_A$n2),
  Latent = if (out_A$ok) paste(out_A$latent_used, collapse = ";") else "",
  stringsAsFactors = FALSE
)

write_res_and_plot(
  out_A, dirs_A,
  base_name = "Overall__Clonally_Expanding_vs_Non_Expanding",
  title = "Overall: Clonally Expanding vs Non-expanding",
  xlab = "avg_log2FC (Expanding vs Non-expanding)"
)

write.csv(bind_rows(manifest_A), file.path(dirs_A$block_dir, "MANIFEST.csv"), row.names = FALSE)

############################################################
# B) SPLIT BY CLUSTER (T_Cell_A): Clonal vs Non-clonal within each cluster
############################################################
dirs_B <- make_block_dirs(out_dir, "SplitBy_Cluster_T_Cell_A__Clonal_vs_NonClonal")
manifest_B <- list()

i <- 0
for (cl in levels(tara_cdnk_tcr$T_Cell_A)) {
  seu_sub <- subset(tara_cdnk_tcr, subset = T_Cell_A == cl)
  
  comp_name <- paste0("Cluster__", safe_name(cl), "__Clonally_Expanding_vs_Non_Expanding")
  
  out <- run_FindMarkers_safe(
    seu = seu_sub,
    group_by = "ClonalStatus",
    ident1 = "Clonally_Expanding",
    ident2 = "Non_Expanding",
    latent_vars = c("nCount_RNA", "Timepoint_Group", "PID"),
    label = comp_name
  )
  
  i <- i + 1
  manifest_B[[i]] <- data.frame(
    Analysis = "SplitBy_Cluster_T_Cell_A",
    Cluster = cl,
    Comparison = "Clonally_Expanding_vs_Non_Expanding",
    OK = out$ok,
    Reason = ifelse(out$ok, "", out$reason),
    n_Expanding = ifelse(is.null(out$n1), NA, out$n1),
    n_NonExpanding = ifelse(is.null(out$n2), NA, out$n2),
    Latent = if (out$ok) paste(out$latent_used, collapse = ";") else "",
    stringsAsFactors = FALSE
  )
  
  if (out$ok) {
    write_res_and_plot(
      out, dirs_B,
      base_name = comp_name,
      title = paste0("Cluster (T_Cell_A): ", cl, " | Expanding vs Non-expanding"),
      xlab = "avg_log2FC (Expanding vs Non-expanding)"
    )
  }
}
write.csv(bind_rows(manifest_B), file.path(dirs_B$block_dir, "MANIFEST.csv"), row.names = FALSE)

############################################################
# C) SPLIT BY ART STATUS (Timepoint_Group): Clonal vs Non-clonal within each ART status
############################################################
dirs_C <- make_block_dirs(out_dir, "SplitBy_ART_Status__Clonal_vs_NonClonal")
manifest_C <- list()

i <- 0
for (tp in levels(tara_cdnk_tcr$Timepoint_Group)) {
  seu_sub <- subset(tara_cdnk_tcr, subset = Timepoint_Group == tp)
  
  comp_name <- paste0("ART_Status__", safe_name(tp), "__Clonally_Expanding_vs_Non_Expanding")
  
  out <- run_FindMarkers_safe(
    seu = seu_sub,
    group_by = "ClonalStatus",
    ident1 = "Clonally_Expanding",
    ident2 = "Non_Expanding",
    latent_vars = c("nCount_RNA", "PID"),  # Timepoint_Group constant in subset
    label = comp_name
  )
  
  i <- i + 1
  manifest_C[[i]] <- data.frame(
    Analysis = "SplitBy_ART_Status",
    ART_Status = tp,
    Comparison = "Clonally_Expanding_vs_Non_Expanding",
    OK = out$ok,
    Reason = ifelse(out$ok, "", out$reason),
    n_Expanding = ifelse(is.null(out$n1), NA, out$n1),
    n_NonExpanding = ifelse(is.null(out$n2), NA, out$n2),
    Latent = if (out$ok) paste(out$latent_used, collapse = ";") else "",
    stringsAsFactors = FALSE
  )
  
  if (out$ok) {
    write_res_and_plot(
      out, dirs_C,
      base_name = comp_name,
      title = paste0("ART status: ", tp, " | Expanding vs Non-expanding"),
      xlab = "avg_log2FC (Expanding vs Non-expanding)"
    )
  }
}
write.csv(bind_rows(manifest_C), file.path(dirs_C$block_dir, "MANIFEST.csv"), row.names = FALSE)

############################################################
# D) EXPANDING ONLY: Pairwise ART status comparisons (overall)
############################################################
dirs_D <- make_block_dirs(out_dir, "ExpandingOnly__ART_Status_Comparisons")
manifest_D <- list()

seu_exp <- subset(tara_cdnk_tcr, subset = ClonalStatus == "Clonally_Expanding")
tp_levels <- levels(droplevels(seu_exp$Timepoint_Group))
tp_pairs  <- combn(tp_levels, 2, simplify = FALSE)

i <- 0
for (pair in tp_pairs) {
  tp1 <- pair[1]; tp2 <- pair[2]
  
  seu_pair <- subset(seu_exp, subset = Timepoint_Group %in% c(tp1, tp2))
  
  comp_name <- paste0("ExpandingOnly__ART_Status__", safe_name(tp1), "_vs_", safe_name(tp2))
  
  out <- run_FindMarkers_safe(
    seu = seu_pair,
    group_by = "Timepoint_Group",
    ident1 = tp1,
    ident2 = tp2,
    latent_vars = c("nCount_RNA", "PID"),
    label = comp_name
  )
  
  i <- i + 1
  manifest_D[[i]] <- data.frame(
    Analysis = "ExpandingOnly__ART_Status_Comparisons",
    Comparison = paste0(tp1, "_vs_", tp2),
    OK = out$ok,
    Reason = ifelse(out$ok, "", out$reason),
    n_tp1 = ifelse(is.null(out$n1), NA, out$n1),
    n_tp2 = ifelse(is.null(out$n2), NA, out$n2),
    Latent = if (out$ok) paste(out$latent_used, collapse = ";") else "",
    stringsAsFactors = FALSE
  )
  
  if (out$ok) {
    write_res_and_plot(
      out, dirs_D,
      base_name = comp_name,
      title = paste0("Expanding only: ", tp1, " vs ", tp2),
      xlab = paste0("avg_log2FC (", tp1, " vs ", tp2, ")")
    )
  }
}
write.csv(bind_rows(manifest_D), file.path(dirs_D$block_dir, "MANIFEST.csv"), row.names = FALSE)

############################################################
# E) EXPANDING ONLY: Suppressed vs Unsuppressed split by cluster (robust)
############################################################
dirs_E <- make_block_dirs(out_dir, "ExpandingOnly__Suppressed_vs_Unsuppressed__SplitBy_Cluster")
manifest_E <- list()

tp_all  <- levels(droplevels(seu_exp$Timepoint_Group))
tp_supp <- tp_all[grepl("supp", tp_all, ignore.case = TRUE)]
tp_uns  <- tp_all[grepl("unsupp|un_supp|unsup", tp_all, ignore.case = TRUE)]

tp_supp <- tp_supp[1]
tp_uns  <- tp_uns[1]

if (!is.na(tp_supp) && !is.na(tp_uns)) {
  
  i <- 0
  for (cl in levels(droplevels(seu_exp$T_Cell_A))) {
    
    # Find cells for this cluster + these two ART statuses
    cells_keep <- WhichCells(
      seu_exp,
      expression = T_Cell_A == cl & Timepoint_Group %in% c(tp_supp, tp_uns)
    )
    
    # If no cells, just record and continue (no crash)
    if (length(cells_keep) == 0) {
      i <- i + 1
      manifest_E[[i]] <- data.frame(
        Analysis = "ExpandingOnly__Suppressed_vs_Unsuppressed__SplitBy_Cluster",
        Cluster = cl,
        Comparison = paste0(tp_supp, "_vs_", tp_uns),
        OK = FALSE,
        Reason = "No cells found for this Cluster x (Suppressed/Unsuppressed)",
        n_supp = 0,
        n_uns = 0,
        Latent = "",
        stringsAsFactors = FALSE
      )
      next
    }
    
    seu_sub <- subset(seu_exp, cells = cells_keep)
    
    comp_name <- paste0("Cluster__", safe_name(cl), "__", safe_name(tp_supp), "_vs_", safe_name(tp_uns))
    
    out <- run_FindMarkers_safe(
      seu = seu_sub,
      group_by = "Timepoint_Group",
      ident1 = tp_supp,
      ident2 = tp_uns,
      latent_vars = c("nCount_RNA", "PID"),
      label = comp_name,
      min_cells_per_group = 30   # adjust if needed
    )
    
    i <- i + 1
    manifest_E[[i]] <- data.frame(
      Analysis = "ExpandingOnly__Suppressed_vs_Unsuppressed__SplitBy_Cluster",
      Cluster = cl,
      Comparison = paste0(tp_supp, "_vs_", tp_uns),
      OK = out$ok,
      Reason = ifelse(out$ok, "", out$reason),
      n_supp = ifelse(is.null(out$n1), NA, out$n1),
      n_uns  = ifelse(is.null(out$n2), NA, out$n2),
      Latent = if (out$ok) paste(out$latent_used, collapse = ";") else "",
      stringsAsFactors = FALSE
    )
    
    if (out$ok) {
      write_res_and_plot(
        out, dirs_E,
        base_name = comp_name,
        title = paste0("Expanding only | Cluster: ", cl, " | ", tp_supp, " vs ", tp_uns),
        xlab = paste0("avg_log2FC (", tp_supp, " vs ", tp_uns, ")")
      )
    }
  }
  
  write.csv(bind_rows(manifest_E),
            file.path(dirs_E$block_dir, "MANIFEST.csv"),
            row.names = FALSE)
  
} else {
  warning("Could not automatically identify suppressed/unsuppressed levels in Timepoint_Group. Check levels(seu_exp$Timepoint_Group).")
}


############################################################
# Global manifest (all blocks)
############################################################
global_manifest <- bind_rows(
  bind_rows(manifest_A),
  bind_rows(manifest_B),
  bind_rows(manifest_C),
  bind_rows(manifest_D),
  bind_rows(manifest_E)
)

write.csv(global_manifest, file.path(out_dir, "MANIFEST__ALL_BLOCKS.csv"), row.names = FALSE)

