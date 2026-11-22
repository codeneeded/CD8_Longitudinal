################################################################################
# Annotation defence for the integrated CP003 + CP018 clusters
#
# Same dot-plot form as the other annotation panels (Scripts/4. Revision/17.Lineage_Markers_PBMC_Scale.R):
# Seurat DotPlot, lightgrey -> #08519C, z-score across the groups shown, percent
# expressed as size. Markers are grouped by the label they defend, in the same
# order as the clusters, so each label's evidence sits on its own diagonal block.
#
# The per-cell annotation is READ from the source data written by
# 14.Joint_UMAPs_Annotated.R, so this panel and the UMAP cannot disagree, and the
# cells removed there (12.abCD8_Filter.R) are absent here too.
#
# Out: Reviewer_Comments/Figures/Fig3B_annotation_dot_plot.png
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
FIG <- file.path(R,"Reviewer_Comments","Figures")
o <- qs_read(file.path(R,"saved_R_data","CP003_CP018_abCD8.qs2"))
src <- read.csv(file.path(FIG,"Fig3B_UMAP_source_data.csv"), check.names=FALSE)
o <- subset(o, cells = intersect(colnames(o), src$cell))
stopifnot(ncol(o) == nrow(src))
lab <- setNames(as.character(src$joint_anno), src$cell)
ORD <- c("Naive/Bystander","Tpex","Transitional Tem","TEMRA/Effector","Cycling Effector","Tc17")
stopifnot(setequal(unique(lab), ORD))
o$joint_anno <- factor(unname(lab[colnames(o)]), levels = rev(ORD))
cat("cells:", ncol(o), "| assays:", paste(Assays(o), collapse=", "), "\n")

# one block of defining genes per label, in cluster order
RNA_BLOCKS <- list(
  "Lineage"           = c("CD3E","CD8A","CD8B","CD4"),
  "Naive/Bystander"   = c("CCR7","SELL","LEF1","TCF7","IL7R"),
  "Tpex"              = c("TOX","IKZF2","ZNF683","KLRC2","XCL1"),
  "Transitional Tem"  = c("CD38","IL2RA","HAVCR2","LAG3","GZMK"),
  "TEMRA/Effector"    = c("GZMB","GZMH","PRF1","GNLY","CX3CR1","FCGR3A"),
  "Cycling Effector"  = c("MKI67","TOP2A","TYMS"),
  "Tc17"              = c("KLRB1","CCR6","IL23R","RORC"))
DefaultAssay(o) <- "RNA"
if (length(Layers(o[["RNA"]])) > 1) o <- JoinLayers(o, assay="RNA")
feats <- unique(unlist(RNA_BLOCKS)); feats <- feats[feats %in% rownames(o)]
missing <- setdiff(unique(unlist(RNA_BLOCKS)), feats)
if (length(missing)) cat("  not in object:", paste(missing, collapse=", "), "\n")

mk <- function(assay, f, title, fn, w = 19, h = 7.5) {
  DefaultAssay(o) <- assay
  p <- DotPlot(o, features = f, group.by = "joint_anno",
               cols = c("lightgrey", "#08519C"), dot.scale = 7) + RotatedAxis() +
    guides(colour = guide_colourbar(title = "Mean expression\n(z-score across\nthe groups shown)"),
           size = guide_legend(title = "Percent expressed")) +
    labs(title = title, x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 15),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 12, face = "bold"))
  ggsave(file.path(FIG, fn), p, width = w, height = h, dpi = 300, bg = "white")
  write.csv(p$data, file.path(FIG, sub("\\.png$", "_source_data.csv", fn)), row.names = FALSE)
  cat("  wrote", fn, "\n")
  invisible(p$data)
}
pd <- mk("RNA", feats, "Integrated CP003 and CP018 αβ CD8+ T cells: defining genes",
         "Fig3B_annotation_dot_plot.png")

# block boundaries as a source table, so the evidence per label is auditable
blk <- do.call(rbind, lapply(names(RNA_BLOCKS), function(b)
  data.frame(block=b, gene=intersect(RNA_BLOCKS[[b]], feats))))
write.csv(blk, file.path(FIG, "Fig3B_annotation_dot_plot_marker_blocks.csv"), row.names=FALSE)

# No protein panel: the FACT-seq libraries carry RNA only.

# the numbers that carry each label, from the plotted data itself
chk <- function(g, cl) { r <- pd[pd$features.plot==g & pd$id==cl, ]; if (nrow(r)) round(r$pct.exp,1) else NA }
cat("\n=== defining values read back from the plotted data ===\n")
for (x in list(c("CCR7","Naive/Bystander"), c("TOX","Tpex"), c("IKZF2","Tpex"), c("CD38","Transitional Tem"),
               c("GZMB","TEMRA/Effector"), c("GZMH","TEMRA/Effector"), c("MKI67","Cycling Effector"),
               c("KLRB1","Tc17"), c("CCR6","Tc17"), c("CD8B","Tc17")))
  cat(sprintf("  %-8s in %-20s %5.1f%%\n", x[1], x[2], chk(x[1], x[2])))
