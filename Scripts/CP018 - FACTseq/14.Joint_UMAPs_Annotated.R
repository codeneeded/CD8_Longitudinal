################################################################################
# Integrated CP003 + CP018 UMAPs, annotated, one panel per file
#
# CELLS: alpha-beta CD8 T cells only (12.abCD8_Filter.R), re-clustered on
# those cells alone (13.Recluster_abCD8.R), as the main analysis removes
# gamma-delta and MAIT cells before clustering CD8 T cells.
#
# LABELS are the published Figure 3 labels. Each joint cluster takes the Figure 3
# label carried by most of its CP003 cells (CP003 is the Figure 3 participant); the
# agreement is written to Fig3B_label_agreement.csv. The one cluster with
# too few CP003 cells (68 cells, 7 from CP003; KLRB1 100%, CCR6 44%, IL23R 19%,
# SLC4A10 3%) is Tc17. Colours are those of Figure_3.R.
#
# Format follows Figure_1.R exactly: point size 0.8 / alpha 0.7, white centroid
# labels with leader lines, corner UMAP arrows, theme_void, coord_fixed, 12x11.
#
# Out: Reviewer_Comments/Figures/Fig3B_UMAP_{annotated,HIV_responsive,donor}.png
#      Reviewer_Comments/Tables/Fig3B_{excluded_cells,cluster_annotation,label_agreement}.csv
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(ggrepel)
})
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
FIG <- file.path(R,"Reviewer_Comments","Figures"); TAB <- file.path(R,"Reviewer_Comments","Tables")
o <- qs_read(file.path(R,"saved_R_data","CP003_CP018_integrated.qs2")); m <- o@meta.data

source(file.path(R, "Scripts", "CP018 - FACTseq", "12.abCD8_Filter.R"))
why <- joint_exclusion(o)
ex <- as.data.frame(table(reason = why[!is.na(why)], participant = m[names(why)[!is.na(why)], "participant"]))
ex <- tidyr::pivot_wider(ex, names_from = participant, values_from = Freq)
ex$Total <- ex$CP003 + ex$CP018
write.csv(ex, file.path(TAB, "Fig3B_excluded_cells.csv"), row.names = FALSE)
keep <- joint_keep(o)
o <- qs_read(file.path(R, "saved_R_data", "CP003_CP018_abCD8.qs2")); m <- o@meta.data
stopifnot("re-clustered object is out of date: re-run 21a" = setequal(colnames(o), keep))
cat("excluded:", sum(!is.na(why)), "| alpha-beta CD8 kept:", ncol(o), "\n")

LABS <- c("Naive/Bystander","Tpex","Transitional Tem","TEMRA/Effector","Cycling Effector","Tc17")
# Rule, not a hand-typed table, so a re-clustering cannot silently mislabel:
# each cluster takes the published Figure 3 label held by most of its CP003 cells.
# A cluster with fewer than 20 CP003 cells cannot be named that way; it is named
# Tc17 only if it is KLRB1 (CD161)+ in >= 90% of cells with CCR6 present, otherwise
# the script stops.
cl_ <- as.character(m$abcd8_clusters)
X_  <- GetAssayData(o, assay = "RNA", layer = "counts")
ANNO <- sapply(levels(m$abcd8_clusters), function(k) {
  i3 <- cl_ == k & m$participant == "CP003"
  if (sum(i3) >= 20) return(names(which.max(table(as.character(m$Fig3_Annotation[i3])))))
  kl <- mean(X_["KLRB1", cl_ == k] > 0); cc <- mean(X_["CCR6", cl_ == k] > 0)
  if (kl >= 0.9 && cc >= 0.2) return("Tc17")
  stop("cluster ", k, " has ", sum(i3), " CP003 cells and no rule-based label")
})
stopifnot(all(ANNO %in% LABS))
ag <- m[m$participant == "CP003", ]
ag$label <- unname(ANNO[as.character(ag$abcd8_clusters)])
ag <- as.data.frame(table(cluster = ag$abcd8_clusters, label = ag$label,
                          published_Fig3 = as.character(ag$Fig3_Annotation)))
ag <- ag[ag$Freq > 0, ]
ag <- do.call(rbind, lapply(split(ag, ag$cluster), function(x) { x$pct <- round(100*x$Freq/sum(x$Freq), 1); x[order(-x$Freq), ] }))
write.csv(ag, file.path(TAB, "Fig3B_label_agreement.csv"), row.names = FALSE)
cat("labels:", paste(names(ANNO), ANNO, sep="=", collapse="  "), "\n")

m$joint_anno <- factor(unname(ANNO[as.character(m$abcd8_clusters)]), levels=LABS)
# HIV-responsive, one rule for both donors (the published CP003 rule, as stated in
# the Methods): the cell's clonotype is seen in more than one cell of its library
# (clonalFrequency, computed per library) and the cell lies outside the
# Naive/Bystander cluster of THIS clustering. The per-donor calls carried in the
# object used different rules (CP018: > 2 cells per timepoint, no Naive/Bystander
# exclusion), so they are replaced here and every downstream panel reads this.
cf_ <- suppressWarnings(as.numeric(as.character(m$clonalFrequency)))
m$HIV_Specific_TCR <- ifelse(!is.na(cf_) & cf_ > 1 & m$joint_anno != "Naive/Bystander",
                             "HIV-Specific TCR", "Other")
cat("HIV-responsive (one rule):", paste(names(table(m$participant[m$HIV_Specific_TCR == "HIV-Specific TCR"])),
    table(m$participant[m$HIV_Specific_TCR == "HIV-Specific TCR"]), collapse = " | "), "\n")
stopifnot("unlabelled cluster" = !any(is.na(m$joint_anno)))

emb <- as.data.frame(Embeddings(o,"umap.abcd8")); names(emb) <- c("UMAP1","UMAP2")
d <- cbind(emb, m[, c("joint_anno","participant","HIV_Specific_TCR")])
d$resp <- d$HIV_Specific_TCR == "HIV-Specific TCR"
write.csv(cbind(cell=rownames(d), d), file.path(FIG,"Fig3B_UMAP_source_data.csv"), row.names=FALSE)
write.csv(as.data.frame(m %>% count(cluster = abcd8_clusters, joint_anno, participant) %>%
          tidyr::pivot_wider(names_from=participant, values_from=n, values_fill=0)),
          file.path(TAB,"Fig3B_cluster_annotation.csv"), row.names=FALSE)

# house arrow placement, from Figure_1.R
pad <- 0.08 * diff(range(d$UMAP1, na.rm=TRUE))
x_arrow <- min(d$UMAP1, na.rm=TRUE) - pad
y_arrow <- min(d$UMAP2, na.rm=TRUE) - pad
XL <- c(x_arrow - 2.6, max(d$UMAP1, na.rm=TRUE) + pad)
YL <- c(y_arrow - 2.6, max(d$UMAP2, na.rm=TRUE) + 1.6*pad)
W_IN <- 12; H_IN <- round(W_IN * diff(YL) / diff(XL) + 0.4, 1)
arrows_ <- function(p) p +
  annotate("segment", x=x_arrow, xend=x_arrow+3.5, y=y_arrow, yend=y_arrow,
           arrow=arrow(length=unit(0.4,"cm"), type="closed"), linewidth=1.2, color="black") +
  annotate("text", x=x_arrow+1.75, y=y_arrow-1.3, label="UMAP1", size=7, fontface="bold") +
  annotate("segment", x=x_arrow, xend=x_arrow, y=y_arrow, yend=y_arrow+3.5,
           arrow=arrow(length=unit(0.4,"cm"), type="closed"), linewidth=1.2, color="black") +
  annotate("text", x=x_arrow-1.3, y=y_arrow+1.75, label="UMAP2", size=7, fontface="bold", angle=90) +
  theme_void() + theme(legend.position="none", plot.margin=margin(10,10,10,10)) +
  coord_fixed(xlim=XL, ylim=YL)

PAL <- c("Naive/Bystander"="#2166AC", "Tpex"="#F4A261", "Transitional Tem"="#1A9850",
         "TEMRA/Effector"="#E63946", "Cycling Effector"="#9D4EDD", "Tc17"="#17BECF")   # Figure_3.R; Tc17 new
cent <- d %>% group_by(joint_anno) %>% summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2), .groups="drop")
p1 <- ggplot(d, aes(UMAP1, UMAP2, color=joint_anno)) +
  geom_point(size=0.8, alpha=0.7) + scale_color_manual(values=PAL) +
  geom_label_repel(data=cent, aes(UMAP1, UMAP2, label=joint_anno), fill="white", color="black",
                   size=7, fontface="bold", label.size=0.7, label.padding=unit(0.4,"lines"),
                   box.padding=2.2, point.padding=0.6, segment.color="grey40", segment.size=0.8,
                   min.segment.length=0, max.overlaps=Inf, force=14, seed=1,
                   show.legend=FALSE)
ggsave(file.path(FIG,"Fig3B_UMAP_annotated.png"), arrows_(p1),
       width=W_IN, height=H_IN, dpi=300, bg="white")

p2 <- ggplot(d[!d$resp, ], aes(UMAP1, UMAP2)) +
  geom_point(size=0.8, alpha=0.5, colour="grey82") +
  geom_point(data=d[d$resp, ], aes(UMAP1, UMAP2), size=0.9, alpha=0.85, colour="#C00000")
ggsave(file.path(FIG,"Fig3B_UMAP_HIV_responsive.png"), arrows_(p2),
       width=W_IN, height=H_IN, dpi=300, bg="white")

DON <- c(CP003="#1F77B4", CP018="#FF7F0E")
# Two donors: a legend, not centroid labels - a donor's median position can fall
# on the other donor's cells.
p3 <- ggplot(d, aes(UMAP1, UMAP2, color=participant)) +
  geom_point(size=0.8, alpha=0.7) + scale_color_manual(values=DON, name=NULL) +
  guides(color=guide_legend(override.aes=list(size=7, alpha=1)))
ggsave(file.path(FIG,"Fig3B_UMAP_donor.png"),
       arrows_(p3) + theme(legend.position=c(0.92, 0.9), legend.text=element_text(size=22, face="bold")),
       width=W_IN, height=H_IN, dpi=300, bg="white")

cat("cells:", nrow(d), " responsive:", sum(d$resp), "\n")
print(as.data.frame(m %>% count(joint_anno, participant) %>%
      tidyr::pivot_wider(names_from=participant, values_from=n, values_fill=0)), row.names=FALSE)
cat("\n  wrote three panels\n")
