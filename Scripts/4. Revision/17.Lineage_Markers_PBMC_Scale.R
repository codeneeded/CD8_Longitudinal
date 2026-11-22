################################################################################
# Figure 1C CD8 clusters: lineage-marker dot plots on a total-PBMC scale, and
# Naive 2 compared between the two clusterings (reviewer 2, major 1b / 1d)
#
#   A dot plot z-scores each feature ACROSS THE PLOTTED GROUPS ONLY. When every
#   group is a CD8 subcluster, z-scoring removes the signal that establishes
#   identity and whichever subcluster is lowest renders blue. Including the other
#   PBMC lineages as reference rows puts the scale back on a total-PBMC footing.
#
# CORRECTION (this version): the first attempt used CD8_Annotation from
# TARA_cd8_HEI_annotated_final.qs2 - the DEDICATED 10-cluster re-clustering, and
# HEI-only. Figure 1C uses neither. It uses the whole-PBMC `Annotation` column
# restricted to the six cd8_cluster_names, on the Figure 1 working subset
# (HEU + HUU + HEI at PreART_Entry). This version reproduces that cell set and
# that vocabulary, so the panel answers the comment about Figure 1C.
#
# Cluster names and both marker panels are READ OUT OF Figure_1.R rather than
# retyped, so this panel cannot drift from the published one.
#
# Naive 2: compares
#      the Figure 1 cluster "Naive 2 CD8" against the dedicated-clustering
#      "Naive CD8 2" by actual cell overlap, not by count alone.
#
# Out: Reviewer_Comments/Figures/R2M1b_05_bubble_RNA_PBMCscale.png
#      Reviewer_Comments/Figures/R2M1b_06_bubble_ADT_PBMCscale.png
#      Reviewer_Comments/Tables/R2M1d_05_Naive2_wholePBMC_vs_dedicated.csv
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr); library(ggplot2)})
# Repository root: CD8_BASE if set, otherwise the first folder above this script
# that holds saved_R_data/ (works under Rscript, source() and an IDE).
R <- local({
  b <- Sys.getenv("CD8_BASE"); if (nzchar(b)) return(normalizePath(b))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f <- if (length(f)) gsub("~\\+~", " ", f[1]) else tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  d <- if (length(f) && nzchar(f)) dirname(normalizePath(f)) else getwd()
  while (!dir.exists(file.path(d, "saved_R_data")) && dirname(d) != d) d <- dirname(d)
  if (!dir.exists(file.path(d, "saved_R_data"))) stop("set CD8_BASE to the repository root", call. = FALSE)
  d
})
FIG <- file.path(R, "Reviewer_Comments", "Figures")
TAB <- file.path(R, "Reviewer_Comments", "Tables")
for (d in c(FIG, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- definitions taken from Figure_1.R, never retyped --------------------------
f1 <- readLines(file.path(R, "Scripts", "3. Manuscript", "Figure_1.R"), warn = FALSE)
grab <- function(name) {
  i <- grep(paste0("^", name, "\\s*<-\\s*(setNames\\()?\\s*c\\("), f1)
  if (!length(i)) stop("cannot find ", name, " in Figure_1.R", call. = FALSE)
  i <- i[1]; depth <- 0L; out <- character()
  for (k in seq(i, length(f1))) {
    ln <- f1[k]; bare <- gsub('"[^"]*"', "", ln); out <- c(out, ln)
    depth <- depth + lengths(regmatches(bare, gregexpr("\\(", bare))) -
                     lengths(regmatches(bare, gregexpr("\\)", bare)))
    if (depth <= 0L && k > i - 1L) break
  }
  eval(parse(text = paste(out, collapse = "\n")))
}
cd8_cluster_names   <- grab("cd8_cluster_names")
rna_heatmap_markers <- grab("rna_heatmap_markers")
adt_heatmap_markers <- grab("adt_heatmap_markers")
adt_gene_to_protein <- grab("adt_gene_to_protein")
cat("Figure 1 CD8 clusters:", paste(cd8_cluster_names, collapse = " | "), "\n")

# Figure 1C's marker list carries CD8A but NOT CD4, so on its own it cannot
# answer the reviewer's literal question - "is Naive 2 really CD8+ T cells? What
# is the CD4 expression in this cluster?". CD4 and CD8B are appended as identity
# markers. The published Figure 1C set is untouched and comes first; these two
# are additions, and the panel title says so.
rna_feats <- c(names(rna_heatmap_markers), "CD8B", "CD4")
adt_feats <- c(names(adt_heatmap_markers), "CD4")
adt_gene_to_protein["CD4"] <- "CD4"

all <- qs_read(file.path(R, "saved_R_data", "TARA_ALL_annotated_final.qs2"))

# ---- the Figure 1 working subset, reproduced exactly ---------------------------
keep_cells <- colnames(all)[
  all$Condition == "HEU" | all$Condition == "HUU" |
  (all$Condition == "HEI" & all$Timepoint_Group == "PreART_Entry")]
sub <- subset(all, cells = keep_cells)
cat("Figure 1 working subset:", ncol(sub), "cells\n"); print(table(sub$Condition))

ann <- as.character(sub$Annotation)
is_cd8 <- ann %in% cd8_cluster_names
cat("\nCD8 cells in the Figure 1 subset:", sum(is_cd8), "\n")
print(table(ann[is_cd8]))

# non-CD8 lineages become reference rows so the colour scale spans total PBMC
broad <- as.character(sub$Broad_CellType)
grp <- ifelse(is_cd8, ann, paste0("[ref] ", broad))
grp[!is_cd8 & grepl("CD8", broad)] <- "[ref] CD8 T (other)"
drop <- grepl("^\\[ref\\] (Other|NA|pDC)$", grp) | is.na(grp)
sub$bubble_grp <- grp
sub <- subset(sub, cells = colnames(sub)[!drop])
sub$bubble_grp <- factor(sub$bubble_grp,
  levels = c(cd8_cluster_names[cd8_cluster_names %in% grp],
             sort(unique(grp[grepl("^\\[ref\\]", grp) & !drop]))))
cat("\nGroups plotted:\n"); print(table(sub$bubble_grp))

mk <- function(assay, feats, rename, title, fn, w = 16, h = 7) {
  DefaultAssay(sub) <- assay
  f <- feats[feats %in% rownames(sub)]
  p <- DotPlot(sub, features = f, group.by = "bubble_grp",
               cols = c("lightgrey", "#08519C"), dot.scale = 7) + RotatedAxis() +
    guides(colour = guide_colourbar(title = "Mean expression\n(z-score across\nthe groups shown)"),
           size = guide_legend(title = "Percent expressed")) +
    labs(title = title, x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 15),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 12, face = "bold"))
  if (!is.null(rename)) {
    lv <- levels(p$data$features.plot)
    p <- p + scale_x_discrete(labels = ifelse(lv %in% names(rename), rename[lv], lv))
  }
  ggsave(file.path(FIG, fn), p, width = w, height = h, dpi = 300, bg = "white")
  write.csv(p$data, file.path(FIG, sub("\\.png$", "_source_data.csv", fn)), row.names = FALSE)
  cat("  wrote", fn, "\n")
}
mk("RNA", rna_feats, NULL,
   "Figure 1C CD8 clusters on a total-PBMC expression scale (RNA; CD8B and CD4 added)",
   "R2M1b_05_bubble_RNA_PBMCscale.png")
mk("ADT", adt_feats, adt_gene_to_protein,
   "Figure 1C CD8 clusters on a total-PBMC expression scale (protein, DSB; CD4 added)",
   "R2M1b_06_bubble_ADT_PBMCscale.png")

# ---- C13: overlap between the two "Naive 2" clusters ---------------------------
cd8d <- qs_read(file.path(R, "saved_R_data", "TARA_cd8_HEI_annotated_final.qs2"))
bare <- function(x) sub("^.*_", "", x)
k_f1  <- paste(as.character(sub$orig.ident), bare(colnames(sub)))[as.character(sub$Annotation) == "Naive 2 CD8"]
k_all <- paste(as.character(all$orig.ident), bare(colnames(all)))[as.character(all$Annotation) == "Naive 2 CD8"]
k_ded <- paste(as.character(cd8d$orig.ident), bare(colnames(cd8d)))[as.character(cd8d$CD8_Annotation) == "Naïve CD8 2"]
cat("\n=== C13 ===\n")
cat(sprintf("  'Naive 2 CD8' in the Figure 1 subset : %d cells\n", length(k_f1)))
cat(sprintf("  'Naive 2 CD8' in the whole object    : %d cells\n", length(k_all)))
cat(sprintf("  dedicated 'Naive CD8 2'              : %d cells\n", length(k_ded)))
cat(sprintf("  shared with the Figure 1 cluster     : %d\n", length(intersect(k_f1, k_ded))))
cat(sprintf("  shared with the whole-object cluster : %d\n", length(intersect(k_all, k_ded))))
write.csv(data.frame(
  cluster = c("Naive 2 CD8 (Figure 1 subset)", "Naive 2 CD8 (whole object)",
              "Naive CD8 2 (dedicated)", "shared: Figure 1 vs dedicated",
              "shared: whole object vs dedicated"),
  n_cells = c(length(k_f1), length(k_all), length(k_ded),
              length(intersect(k_f1, k_ded)), length(intersect(k_all, k_ded)))),
  file.path(TAB, "R2M1d_05_Naive2_wholePBMC_vs_dedicated.csv"), row.names = FALSE)
