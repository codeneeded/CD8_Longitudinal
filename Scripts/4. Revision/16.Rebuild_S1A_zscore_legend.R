################################################################################
# Revision — rebuild Supplementary S1A with the colour scale named
#
# ANSWERS: Reviewer 2 minor 6, and Lesley's comment 11 ("where is the updated
#          figure S1A and legend?").
#
# THE CONFUSION. The reviewer read the -1 to 2 colour range as DSB units and
# noted that DSB-normalised protein should sit at 5-20 SD above background for
# a positive population. The range is low because DotPlot plots the z-score of
# the group mean ACROSS the plotted groups, not the DSB value. Nothing about
# the data changes here; the colour bar is simply named for what it shows.
#
# Written to Reviewer_Comments/Figures so the manuscript panel is untouched.
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(ggplot2)})
source("config_revision.R")

seu <- qs_read(QS_ALL)
DefaultAssay(seu) <- "ADT"

adt <- c("CD3D","CD4","CD8A","CD19","MS4A1","CD79B","IGHM","IGHD","NCAM1","KLRF1",
         "KLRD1","CD14","FCGR3A","ITGAM","FCER1A","CD1C","CLEC4C","HLA-DRA","CD38",
         "IL7R","TCR-AB","TCR-vD2","CD56","PDCD1","TIGIT","CD45RA","CD45RO","CD27",
         "CD28","ITGB7","FAS")
disp <- c(CD3D="CD3",CD4="CD4",CD8A="CD8a",CD19="CD19",MS4A1="CD20",CD79B="CD79b",
          IGHM="IgM",IGHD="IgD",NCAM1="CD56",KLRF1="NKp80",KLRD1="CD94",CD14="CD14",
          FCGR3A="CD16",ITGAM="CD11b",FCER1A="FcER1a",CD1C="CD1c",CLEC4C="BDCA-2",
          `HLA-DRA`="HLA-DR",CD38="CD38",IL7R="CD127",`TCR-AB`="TCRab",
          `TCR-vD2`="TCR Vd2",CD56="CD56",PDCD1="PD-1",TIGIT="TIGIT",CD45RA="CD45RA",
          CD45RO="CD45RO",CD27="CD27",CD28="CD28",ITGB7="Integrin b7",FAS="CD95")
avail <- adt[adt %in% rownames(seu)]
cat("ADT features plotted:", length(avail), "of", length(adt), "\n")

ord <- c("CD4 T cells","CD8 T cells","\u03b3\u03b4 / MAIT","NK cells","B cells",
         "Monocytes","pDC","Other")
seu$Broad_CellType <- factor(as.character(seu$Broad_CellType),
                             levels = intersect(ord, unique(as.character(seu$Broad_CellType))))

p <- DotPlot(seu, features = avail, group.by = "Broad_CellType",
             cols = c("lightgrey", "#08519C"), dot.scale = 8) + RotatedAxis() +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(disp), disp[x], x)) +
  guides(colour = guide_colourbar(title = "Mean expression\n(z-score across groups)"),
         size = guide_legend(title = "Percent expressed")) +
  labs(title = "S1A: Broad PBMC annotation \u2014 surface protein markers (ADT)",
       x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 11))

save_fig("R2m6_02_S1A_PBMC_Lineage_DotPlot_zscore", p, width = 18, height = 8)
save_src("R2m6_02_S1A_PBMC_Lineage_DotPlot_zscore", p$data)
cat("\ncolour bar now names the scale; DSB values themselves are unchanged\n")
