################################################################################
# Lesley items 1, 13 and 24 - proposed replacements, written to Reviewer_Comments
#
#  1: S1B sub-gating - legend, dots and axis text are illegible at manuscript
#     size. The text was already 22-28 pt, but the canvas was 18 x 9 in, so
#     everything shrinks when the panel is placed at column width. Fix is the
#     canvas, not the font: 13 x 7 in, larger dots, large legend keys.
# 13: Fig 3D and 3F for both donors.
# 24: S5C - say what "source cluster" means in the facet label and the caption.
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(patchwork); library(tidyr)
})
R   <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
FIG <- file.path(R,"Reviewer_Comments","Figures"); TAB <- file.path(R,"Reviewer_Comments","Tables")

cd8 <- qs_read(file.path(R,"saved_R_data","TARA_cd8_HEI_annotated_final.qs2"))
DefaultAssay(cd8) <- "ADT"; adt <- GetAssayData(cd8, layer="data")
ann <- as.character(cd8$CD8_Annotation)

# ---- item 1 -----------------------------------------------------------------
# Replotted from the PUBLISHED panel's own source data, so the cells, the
# grouping and the labels are identical to S1B by construction - no risk of
# rebuilding it on the wrong object or the wrong annotation vocabulary.
# The text was already 22-28 pt; it reads small because the canvas is 18 x 9 in,
# so the fix is the canvas and the mark size, not the font size.
S1 <- file.path(R,"Manuscript","Supplementary 1","source_data")
g  <- read.csv(file.path(S1,"S1B_CD8_SubGating_panel1_data_source_data.csv"), check.names=FALSE)
cat("item 1: replotted from published source data,", nrow(g), "cells\n")
print(table(g$Gate))
g$Gate <- factor(g$Gate, levels=c("Naive 1","Naive Intermediate","Tscm"))
# Naive 1 is a real population, so it gets a real colour. Grey stays for
# absence, reference rows and the unexpressed end of a scale.
COLS <- c("Naive 1"="#2CA02C", "Naive Intermediate"="#1F78B4", "Tscm"="#E31A1C")
# draw the bulk naive population first so the two gated populations stay visible
g <- g[order(g$Gate), ]
sg <- function(x, xlab, ttl) ggplot(g, aes(.data[[x]], CD45RO, colour=Gate)) +
  geom_point(size=2.4, alpha=.75) +
  scale_colour_manual(values=COLS, name=NULL) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40", linewidth=.7) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=.7) +
  labs(title=ttl, x=xlab, y="CD45RO (DSB)") +
  guides(colour=guide_legend(override.aes=list(size=7, alpha=1))) +
  theme_classic(base_size=19) +
  theme(plot.title=element_text(face="bold", size=20, hjust=.5),
        axis.title=element_text(face="bold", size=18),
        axis.text=element_text(size=16, colour="black"),
        legend.text=element_text(size=17), legend.position="bottom")
p1 <- (sg("FAS","CD95/FAS (DSB)","CD45RO vs CD95") | sg("CD45RA","CD45RA (DSB)","CD45RO vs CD45RA")) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave(file.path(FIG,"L3_item01_S1B_subgating_legible.png"), p1,
       width=13, height=7, dpi=300, bg="white")
cat("  wrote L3_item01_S1B_subgating_legible.png\n")

# ---- item 24 ----------------------------------------------------------------
if ("source_cluster" %in% names(cd8@meta.data)) {
  m <- cd8@meta.data
  sub_ <- ann %in% c("Tscm CD8","Naïve Intermediate CD8")
  d24 <- data.frame(CD45RO=as.numeric(adt["CD45RO",sub_]), FAS=as.numeric(adt["FAS",sub_]),
                    pop=ann[sub_], src=as.character(m$source_cluster[sub_]))
  nm <- tapply(d24$pop, d24$src, function(v) paste(sort(unique(v)), collapse=" / "))
  # Item 24: the facet said "(source cluster N)" with no definition anywhere.
  # Short enough to fit the strip, and no arrow glyph: U+2192 renders as a box
  # in the default device font.
  d24$facet <- paste0("WNN cluster ", d24$src, " before sub-gating\n(",
                      gsub(" CD8", "", nm[d24$src]), ")")
  p24 <- ggplot(d24, aes(FAS, CD45RO, colour=pop)) +
    geom_point(size=2.6, alpha=.7) +
    scale_colour_manual(values=c("Tscm CD8"="#E31A1C","Naïve Intermediate CD8"="#1F78B4"),
      name=NULL, labels=c("Tscm CD8"="Tscm (CD45RO\u2212 CD95+)",
                          "Naïve Intermediate CD8"="Na\u00efve Intermediate (CD45RO+)")) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey40", linewidth=.7) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=.7) +
    facet_wrap(~facet, nrow=1) +
    guides(colour=guide_legend(override.aes=list(size=7, alpha=1))) +
    labs(title="Surface sub-gating of the na\u00efve clusters",
         subtitle=NULL,
         x="CD95/FAS (DSB)", y="CD45RO (DSB)") +
    theme_classic(base_size=18) +
    theme(plot.title=element_text(face="bold", size=20),
          plot.subtitle=element_text(size=14, colour="grey25"),
          strip.background=element_rect(fill="grey95", colour=NA),
          strip.text=element_text(size=12.5, face="bold", lineheight=1.05),
          axis.title=element_text(face="bold", size=17),
          axis.text=element_text(size=15, colour="black"),
          legend.position="bottom", legend.text=element_text(size=16))
  ggsave(file.path(FIG,"L3_item24_S5C_sourcecluster_legend.png"), p24,
         width=17, height=7, dpi=300, bg="white")
  cat("  wrote L3_item24_S5C_sourcecluster_legend.png | facets:",
      length(unique(d24$facet)), "\n")
}

# Item 13 (Figure 3D and 3F, both donors) is built in 23.Figure3_FACTseq_both_donors.R.
