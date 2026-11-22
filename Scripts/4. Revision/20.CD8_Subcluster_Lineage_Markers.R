################################################################################
# Supplementary Figure 5: lineage-marker dot plots for the ten dedicated CD8
# clusters, with the non-CD8 PBMC lineages as reference rows
#
# Same correction as for the Figure 1C panels (17.Lineage_Markers_PBMC_Scale.R): a dot plot z-scores each
# feature across the groups plotted, so a panel containing only CD8 subclusters
# cannot establish that they are CD8. Including the non-CD8 PBMC lineages as
# reference rows puts the scale on a total-PBMC footing.
#
# The 10 dedicated CD8 clusters live in TARA_cd8_HEI_annotated_final.qs2
# (CD8_Annotation); the lineages come from TARA_ALL_annotated_final.qs2
# (Broad_CellType), restricted to the same HEI cells so the comparison is within
# one cohort rather than across different participants.
#
# Out: Reviewer_Comments/Figures/FigS5_lineage_markers_RNA.png
#      Reviewer_Comments/Figures/FigS5_lineage_markers_ADT.png
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

cd8 <- qs_read(file.path(R,"saved_R_data","TARA_cd8_HEI_annotated_final.qs2"))
all <- qs_read(file.path(R,"saved_R_data","TARA_ALL_annotated_final.qs2"))
all <- subset(all, cells = colnames(all)[all$Condition == "HEI"])
cat("HEI PBMC:", ncol(all), " dedicated CD8:", ncol(cd8), "\n")

bc  <- as.character(all$Broad_CellType)
cd8_bc <- colnames(cd8)
grp_all <- ifelse(colnames(all) %in% cd8_bc, NA, paste0("[ref] ", bc))
keep <- !is.na(grp_all) & !grepl("^\\[ref\\] (Other|NA|pDC)$", grp_all)
ref <- subset(all, cells = colnames(all)[keep])
ref$bubble_grp <- grp_all[keep]
cd8$bubble_grp <- as.character(cd8$CD8_Annotation)
cat("reference lineage rows:", paste(sort(unique(ref$bubble_grp)), collapse=" | "), "\n")
cat("CD8 cluster rows:", length(unique(cd8$bubble_grp)), "\n")

mk <- function(assay, feats, title, fn, w=17, h=8) {
  DefaultAssay(ref) <- assay; DefaultAssay(cd8) <- assay
  f <- intersect(feats, intersect(rownames(ref), rownames(cd8)))
  mg <- merge(DietSeurat(cd8, assays=assay, features=f),
              DietSeurat(ref, assays=assay, features=f))
  if (length(Layers(mg[[assay]])) > 1) mg <- JoinLayers(mg, assay=assay)
  mg$bubble_grp <- factor(mg$bubble_grp,
    levels = c(sort(unique(as.character(cd8$bubble_grp))),
               sort(unique(as.character(ref$bubble_grp)))))
  p <- DotPlot(mg, features=f, group.by="bubble_grp",
               cols=c("lightgrey","#08519C"), dot.scale=7) + RotatedAxis() +
    guides(colour=guide_colourbar(title="Mean expression\n(z-score across\nthe groups shown)"),
           size=guide_legend(title="Percent expressed")) +
    labs(title=title, x=NULL, y=NULL) +
    theme(plot.title=element_text(face="bold", size=15),
          axis.text.x=element_text(size=11),
          axis.text.y=element_text(size=12, face="bold"))
  ggsave(file.path(FIG, fn), p, width=w, height=h, dpi=300, bg="white")
  write.csv(p$data, file.path(FIG, sub("\\.png$","_source_data.csv", fn)), row.names=FALSE)
  cat("  wrote", fn, "\n")
}
mk("RNA", c("CD3D","CD3E","CD8A","CD8B","CD4","TRAC","TRDC","KLRB1","SLC4A10",
            "NKG7","GNLY","KLRD1","NCAM1","FCGR3A","KIR2DL3",
            "MS4A1","CD79A","LYZ","CD14","FCN1","PPBP","HBB",
            "CCR7","SELL","TCF7","LEF1","IL7R","GZMB","PRF1","TOX","PDCD1","MKI67"),
   "Dedicated CD8 clusters and PBMC lineages on one expression scale (RNA)",
   "FigS5_lineage_markers_RNA.png")
mk("ADT", c("CD3","CD8A","CD4","CD19","CD14","CD16","CD56","CD45RA","CD45RO",
            "CD27","CD28","CD95","KLRG1","KIR3DL1","CX3CR1","PD-1","TIGIT","HLA-DRA"),
   "Dedicated CD8 clusters and PBMC lineages on one expression scale (protein, DSB)",
   "FigS5_lineage_markers_ADT.png")
