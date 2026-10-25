################################################################################
# Lesley comments 12 and 13
#
# C12: "the relative expression using only CD8 subclusters makes it seem like some
#       of the clusters are not even CD8 T cells due to low (blue) scaled mean
#       expression... he/she wants to see relative expression of these markers in
#       these clusters but relative to total PBMC."
#
#   The original S1A/Fig1C bubble plots z-score each feature ACROSS THE PLOTTED
#   GROUPS ONLY. When the only groups are CD8 subclusters, every one of them is
#   CD8-positive, so z-scoring removes the very signal that establishes identity
#   and the lowest subcluster renders blue. Adding the other PBMC lineages as
#   reference groups puts the scale back on a total-PBMC footing: a CD8
#   subcluster is then blue only if it is genuinely low FOR A PBMC.
#
# C13: "does naive 2 have the same number of cells and composition?" - compares
#      the whole-PBMC "Naive 2 CD8" against the dedicated-clustering "Naive CD8 2"
#      by actual cell overlap, not just by count.
#
# Out: Reviewer_Comments/Figures/R2M1b_05_bubble_RNA_PBMCscale.png
#      Reviewer_Comments/Figures/R2M1b_06_bubble_ADT_PBMCscale.png
#      Reviewer_Comments/Tables/R2M1d_05_Naive2_wholePBMC_vs_dedicated.csv
#      (these answer reviewer comments, so they live with the other
#       reviewer panels rather than in the exploratory tree)
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr); library(ggplot2)})
R   <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
FIG <- file.path(R, "Reviewer_Comments", "Figures")
TAB <- file.path(R, "Reviewer_Comments", "Tables")
for (d in c(FIG, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
OUT <- FIG

all <- qs_read(file.path(R,"saved_R_data","TARA_ALL_annotated_final.qs2"))
cd8 <- qs_read(file.path(R,"saved_R_data","TARA_cd8_HEI_annotated_final.qs2"))

# ---- C13: overlap, not just counts -------------------------------------------
bare <- function(x) sub("^.*_", "", x)
a <- data.frame(bc = bare(colnames(all)), samp = as.character(all$orig.ident),
                ann = as.character(all$Annotation))
d <- data.frame(bc = bare(colnames(cd8)), samp = as.character(cd8$orig.ident),
                ann = as.character(cd8$CD8_Annotation))
a$key <- paste(a$samp, a$bc); d$key <- paste(d$samp, d$bc)
wp <- a$key[a$ann == "Naive 2 CD8"]
dd <- d$key[d$ann == "Naïve CD8 2"]
cat("=== C13 ===\n")
cat(sprintf("  whole-PBMC 'Naive 2 CD8'      : %d cells\n", length(wp)))
cat(sprintf("  dedicated  'Naive CD8 2'      : %d cells\n", length(dd)))
cat(sprintf("  shared cells                  : %d\n", length(intersect(wp, dd))))
cat(sprintf("  %% of dedicated drawn from it   : %.1f%%\n", 100*length(intersect(wp,dd))/max(1,length(dd))))
cat(sprintf("  %% of whole-PBMC cluster kept   : %.1f%%\n", 100*length(intersect(wp,dd))/max(1,length(wp))))
cat("  where the rest of the whole-PBMC cluster went in the dedicated clustering:\n")
rest <- d$ann[d$key %in% setdiff(wp, dd)]
print(head(sort(table(rest), decreasing = TRUE), 6))
cat(sprintf("  (and %d of its cells are absent from the CD8 object entirely - the protein gate)\n",
            length(setdiff(wp, d$key))))
write.csv(data.frame(metric = c("whole_PBMC_Naive2CD8","dedicated_NaiveCD8_2","shared"),
                     n = c(length(wp), length(dd), length(intersect(wp,dd)))),
          file.path(TAB,"R2M1d_05_Naive2_wholePBMC_vs_dedicated.csv"), row.names = FALSE)

# ---- C12: bubble plot on a total-PBMC scale ----------------------------------
lab <- as.character(all$Annotation)
cd8_map <- setNames(as.character(cd8$CD8_Annotation), paste(as.character(cd8$orig.ident), bare(colnames(cd8))))
akey <- paste(as.character(all$orig.ident), bare(colnames(all)))
grp  <- ifelse(akey %in% names(cd8_map), cd8_map[akey], NA_character_)
# non-CD8 lineages collapsed to a single reference row each
broad <- as.character(all$Broad_CellType)
grp[is.na(grp)] <- paste0("[ref] ", broad[is.na(grp)])
all$bubble_grp <- grp
keep <- !grepl("^\\[ref\\] (Other|pDC)$", grp)
cat("\n=== C12 groups ===\n"); print(sort(table(grp[keep]), decreasing = TRUE))

ord_cd8 <- sort(unique(grp[grp %in% cd8_map]))
ord_ref <- sort(unique(grp[grepl("^\\[ref\\]", grp) & keep]))
sub <- subset(all, cells = colnames(all)[keep])
sub$bubble_grp <- factor(sub$bubble_grp, levels = c(ord_cd8, ord_ref))

mk <- function(assay, feats, title, fn, w = 15, h = 7.5) {
  DefaultAssay(sub) <- assay
  f <- feats[feats %in% rownames(sub)]
  p <- DotPlot(sub, features = f, group.by = "bubble_grp",
               cols = c("lightgrey","#08519C"), dot.scale = 7) + RotatedAxis() +
    guides(colour = guide_colourbar(title = "Mean expression\n(z-score across ALL\ngroups shown)"),
           size = guide_legend(title = "Percent expressed")) +
    labs(title = title, x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", size = 15),
          axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 10))
  ggsave(file.path(OUT, fn), p, width = w, height = h, dpi = 300, bg = "white")
  write.csv(p$data, file.path(FIG, sub("\\.png$","_source_data.csv", fn)), row.names = FALSE)
  cat("  wrote", fn, "\n")
}
mk("RNA", c("CD8A","CD8B","CD4","CD3D","CD3E","SELL","CCR7","TCF7","LEF1","IL7R",
            "GZMB","GNLY","NKG7","PRF1","KLRB1","FCGR3A","TOX","PDCD1"),
   "CD8 subclusters on a total-PBMC expression scale (RNA)",
   "R2M1b_05_bubble_RNA_PBMCscale.png")
mk("ADT", c("CD8A","CD4","CD3D","CD45RA","CD45RO","CD27","CD28","IL7R","FAS",
            "PDCD1","TIGIT","NCAM1","FCGR3A","KLRG1"),
   "CD8 subclusters on a total-PBMC expression scale (protein, DSB)",
   "R2M1b_06_bubble_ADT_PBMCscale.png")
