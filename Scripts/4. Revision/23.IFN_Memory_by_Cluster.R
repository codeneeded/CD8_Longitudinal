################################################################################
# Type I IFN memory module score by cluster in the longitudinal CD8 analysis,
# for the Results statement on IFN memory at the naive terminus and in effectors
#
# The module-score panels in the current figures come
# from the CP003 stimulation data, whereas Naive CD8 3 is a cluster in the
# LONGITUDINAL dedicated analysis. This builds the panel from the longitudinal
# object, using the IFN_Memory gene set already defined by the project's own
# module-scoring script so the panel and the pipeline agree.
#
# Out: Reviewer_Comments/Figures/IFN_memory_by_cluster_longitudinal.png
#      Reviewer_Comments/Tables/IFN_memory_by_cluster.csv
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(rstatix); library(patchwork)
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

# the project's own definition, from Scripts/2. Analysis/3.Module_Scoring.R
IFN_MEM <- c("IFIT1","IFIT2","IFIT3","ISG15","MX1","MX2","OAS1","OAS2","OAS3","IFI6")
o <- qs_read(file.path(R,"saved_R_data","TARA_cd8_HEI_annotated_final.qs2"))
DefaultAssay(o) <- "RNA"
if (length(Layers(o[["RNA"]])) > 1) o <- JoinLayers(o, assay="RNA")
g <- intersect(IFN_MEM, rownames(o))
cat("IFN_Memory genes present:", length(g), "of", length(IFN_MEM), "-", paste(g, collapse=", "), "\n")
stopifnot(length(g) >= 5)
o <- AddModuleScore(o, features=list(g), name="IFNmem_", seed=42)

ORD <- c("Naïve CD8","Naïve CD8 2","Naïve CD8 3","Naïve CD8 4","Naïve Intermediate CD8",
         "Tscm CD8","TEM CD8","TEMRA CD8","CD16+ Effector CD8","KIR+ innate-like CD8")
d <- data.frame(cluster=factor(as.character(o$CD8_Annotation), levels=ORD),
                score=o$IFNmem_1)
d <- d[!is.na(d$cluster), ]
# naive compartment vs everything downstream, for the "established early" claim
NAIVE <- ORD[1:5]
d$stage <- factor(ifelse(as.character(d$cluster) %in% NAIVE, "Naïve compartment",
                         "Effector / memory"), levels=c("Naïve compartment","Effector / memory"))
tb <- d %>% group_by(cluster, stage) %>%
  summarise(n_cells=n(), median_score=round(median(score),4),
            mean_score=round(mean(score),4), .groups="drop")
# each cluster against the naive clusters OTHER than Naive CD8 3, the comparison
# the sentence rests on
ref <- d$score[as.character(d$cluster) %in% setdiff(NAIVE, "Naïve CD8 3")]
tb$p_vs_other_naive <- vapply(as.character(tb$cluster), function(cl) {
  x <- d$score[as.character(d$cluster)==cl]
  if (cl %in% setdiff(NAIVE,"Naïve CD8 3")) return(NA_real_)
  suppressWarnings(wilcox.test(x, ref)$p.value)
}, numeric(1))
tb$p_adj <- p.adjust(tb$p_vs_other_naive, method="BH")
tb$higher <- ifelse(is.na(tb$p_vs_other_naive), NA,
                    tb$median_score > median(ref))
cat("\n=== IFN memory by cluster (reference = naive clusters other than Naive CD8 3) ===\n")
print(as.data.frame(tb), row.names=FALSE)
write.csv(tb, file.path(TAB,"IFN_memory_by_cluster.csv"), row.names=FALSE)

# Median and mean disagree for Naive CD8 3 (median below the other naive
# clusters, mean above), which is the signature of a right-tailed distribution:
# a SUBSET of cells carries the programme rather than the whole cluster. The
# fraction above a fixed high threshold is the statistic that captures that, so
# it is reported alongside the distribution rather than instead of it.
thr <- quantile(ref, 0.90)
frac <- d %>% group_by(cluster) %>%
  summarise(pct_IFN_high = round(100*mean(score > thr), 1), .groups="drop")
tb <- left_join(tb, frac, by="cluster")
cat("\n=== fraction of cells above the 90th percentile of the naive reference ===\n")
print(as.data.frame(tb[, c("cluster","stage","n_cells","median_score","mean_score","pct_IFN_high","p_adj")]),
      row.names=FALSE)
write.csv(tb, file.path(TAB,"IFN_memory_by_cluster.csv"), row.names=FALSE)

lab <- tb %>% mutate(txt=ifelse(is.na(p_adj), "reference",
                     ifelse(p_adj < 0.001, "p < 0.001", paste0("p = ", signif(p_adj,2)))))
pA <- ggplot(d, aes(cluster, score, fill=stage)) +
  geom_violin(scale="width", linewidth=.35, colour="grey30") +
  geom_boxplot(width=.14, outlier.shape=NA, fill="white", linewidth=.4) +
  geom_text(data=lab, aes(cluster, max(d$score)*1.02, label=txt), inherit.aes=FALSE,
            size=3.8, angle=0, vjust=0) +
  scale_fill_manual(values=c("Naïve compartment"="#9ECAE1","Effector / memory"="#E6550D"), name=NULL) +
  labs(title="Type I interferon memory module across the CD8 compartment",
       x=NULL, y="IFN memory module score") +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold", size=17),
        axis.text.x=element_text(angle=35, hjust=1, size=12, colour="black"),
        legend.position="top")
pB <- ggplot(d, aes(stage, score, fill=stage)) +
  geom_violin(scale="width", linewidth=.35, colour="grey30", show.legend=FALSE) +
  geom_boxplot(width=.16, outlier.shape=NA, fill="white", linewidth=.4) +
  scale_fill_manual(values=c("Naïve compartment"="#9ECAE1","Effector / memory"="#E6550D")) +
  labs(title="Naïve and effector compartments", x=NULL, y=NULL) +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold", size=15), legend.position="none",
        axis.text.x=element_text(size=12, colour="black"))
pC <- ggplot(tb, aes(cluster, pct_IFN_high, fill=stage)) +
  geom_col(width=.7) +
  geom_text(aes(label=sprintf("%.1f%%", pct_IFN_high)), vjust=-0.35, size=4, fontface="bold") +
  scale_fill_manual(values=c("Naïve compartment"="#9ECAE1","Effector / memory"="#E6550D"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  labs(title="Proportion of cells above the reference threshold",
       subtitle="above the 90th percentile of the naïve reference",
       x=NULL, y="% of cells in the cluster") +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold", size=16),
        plot.subtitle=element_text(size=12, colour="grey25"),
        axis.text.x=element_text(angle=35, hjust=1, size=12, colour="black"),
        legend.position="none")
ggsave(file.path(FIG,"IFN_memory_by_cluster_longitudinal.png"),
       (pA / pC) + plot_layout(heights=c(1,.85)), width=15, height=13, dpi=300, bg="white")
write.csv(d, file.path(FIG,"IFN_memory_by_cluster_longitudinal_source_data.csv"), row.names=FALSE)
n3 <- tb[tb$cluster=="Naïve CD8 3", ]
cat("\nNaive CD8 3: median", n3$median_score, " p_adj", signif(n3$p_adj,3),
    " higher than other naive:", n3$higher, "\n")
eff <- tb[tb$stage=="Effector / memory", ]
cat("effector clusters higher than other naive:", sum(eff$higher, na.rm=TRUE), "of", nrow(eff), "\n")
cat("  wrote the panel and table\n")
