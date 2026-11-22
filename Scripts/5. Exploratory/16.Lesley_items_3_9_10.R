################################################################################
# Lesley items 3, 9 and 10 - proposed replacements, written to Reviewer_Comments
#
#  9: left panel - remove the "Cluster span" title, x-axis to "Number of
#     Clusters Occupied", restore the number on top of the bar for "1".
# 10: right panel - title to "Shared clonotypes between clusters".
#  3: the naive comparison panels become Tscm vs Naive Intermediate (left) and
#     Naive 1 vs Naive Intermediate (right), with larger legend dots and font.
#     This is NOT a retitle: the published panel compares Naive 1 vs Naive 2, so
#     the comparison is recomputed from the CD8 object.
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(patchwork)})
R   <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
FIG <- file.path(R,"Reviewer_Comments","Figures")
S2  <- file.path(R,"Manuscript","Supplementary 2","source_data")

p1 <- read.csv(file.path(S2,"S2C_CrossCluster_ClonalSharing_panel1_data_source_data.csv"), check.names=FALSE)
p2 <- read.csv(file.path(S2,"S2C_CrossCluster_ClonalSharing_panel2_data_source_data.csv"), check.names=FALSE)
cat("panel1 cols:", paste(names(p1), collapse=", "), "\n")
cat("panel2 cols:", paste(names(p2), collapse=", "), "\n")
# panel1 is PER CLONOTYPE (one row per CTstrict), not pre-aggregated, so it has
# to be counted into the histogram the panel actually shows.
spanc <- p1 %>% mutate(bucket = ifelse(n_clusters >= 3, "3+", as.character(n_clusters))) %>%
  count(bucket, name="n_clonotypes") %>%
  mutate(bucket = factor(bucket, levels = c("1","2","3+")))
cat("  cluster span:\n"); print(as.data.frame(spanc), row.names=FALSE)
gA <- ggplot(spanc, aes(bucket, n_clonotypes)) +
  geom_col(fill="#4C72B0", width=.7) +
  # item 9: the number is restored on EVERY bar, including "1"
  geom_text(aes(label=n_clonotypes), vjust=-0.35, size=5.6, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(0,.16))) +
  labs(title=NULL, x="Number of Clusters Occupied", y="Number of clonotypes") +
  theme_classic(base_size=17) +
  theme(axis.title=element_text(face="bold"), axis.text=element_text(colour="black", size=14))

# panel2 is a cluster-by-cluster sharing matrix; drop the diagonal and keep the
# upper triangle so each pair appears once.
sh <- p2 %>% filter(Cluster1 != Cluster2, n_shared > 0) %>%
  mutate(a = pmin(Cluster1, Cluster2), b = pmax(Cluster1, Cluster2)) %>%
  distinct(a, b, n_shared) %>% mutate(pair = paste(a, "&", b)) %>%
  arrange(desc(n_shared))
cat("  cluster pairs sharing clonotypes:", nrow(sh), "\n")
gB <- ggplot(sh, aes(reorder(pair, n_shared), n_shared)) +
  geom_col(fill="#C44E52", width=.7) + coord_flip() +
  geom_text(aes(label=n_shared), hjust=-0.25, size=4.8, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  labs(title="Shared clonotypes between clusters", x=NULL, y="Number of clonotypes") +
  theme_classic(base_size=16) +
  theme(plot.title=element_text(face="bold", size=18),
        axis.text=element_text(colour="black", size=13))
ggsave(file.path(FIG,"L3_item09_10_cluster_span_and_sharing.png"), gA | gB,
       width=17, height=6.4, dpi=300, bg="white")
cat("  wrote L3_item09_10_cluster_span_and_sharing.png\n")

# ---- item 3 --------------------------------------------------------------
cd8 <- qs_read(file.path(R,"saved_R_data","TARA_cd8_HEI_annotated_final.qs2"))
DefaultAssay(cd8) <- "ADT"
adt <- GetAssayData(cd8, layer="data")
ann <- as.character(cd8$CD8_Annotation)
cat("  populations:", paste(sort(unique(ann)), collapse=" | "), "\n")
pick <- function(pat) { h <- grep(pat, unique(ann), value=TRUE); if (!length(h)) NA else h[1] }
NI <- pick("Intermediate"); TS <- pick("Tscm"); N1 <- pick("^Na.*ve CD8$")
if (is.na(N1)) N1 <- pick("Na.*ve CD8 1")
cat("  comparing:", TS, "vs", NI, " and ", N1, "vs", NI, "\n")
mkpair <- function(a, b, ttl) {
  ca <- which(ann==a); cb <- which(ann==b)
  stopifnot(length(ca)>0, length(cb)>0)
  d <- data.frame(Marker=rownames(adt),
                  A=rowMeans(adt[, ca, drop=FALSE] > 0)*100,
                  B=rowMeans(adt[, cb, drop=FALSE] > 0)*100)
  d$Diff <- d$A - d$B
  d <- d[order(-abs(d$Diff)), ][1:20, ]
  d$Marker <- factor(d$Marker, levels=rev(d$Marker))
  dl <- tidyr::pivot_longer(d, c(A,B), names_to="grp", values_to="pct")
  dl$grp <- factor(ifelse(dl$grp=="A", a, b), levels=c(a,b))
  ggplot(dl, aes(pct, Marker, colour=grp)) +
    geom_line(aes(group=Marker), colour="grey70", linewidth=.7) +
    geom_point(size=4.6) +
    scale_colour_manual(values=setNames(c("#55A868","#C44E52"), c(a,b)), labels=SHORT, name=NULL) +
    labs(title=ttl, x="% of cells positive (DSB > 0)", y=NULL) +
    theme_classic(base_size=16) +
    theme(plot.title=element_text(face="bold", size=18),
          legend.position="bottom", legend.text=element_text(size=14),
          axis.text.y=element_text(size=13, face="bold", colour="black")) +
    guides(colour=guide_legend(override.aes=list(size=6)))
}
# Titles and legend labels are the ones Lesley asked for (item 3), not the
# CD8_Annotation strings, which carry a "CD8" suffix and call Naive 1 "Naive CD8".
SHORT <- setNames(c("Tscm", "Naïve Intermediate", "Naïve 1"), c(TS, NI, N1))
g3 <- mkpair(TS, NI, "Tscm vs Naïve Intermediate") | mkpair(N1, NI, "Naïve 1 vs Naïve Intermediate")
ggsave(file.path(FIG,"L3_item03_naive_comparisons.png"), g3,
       width=16, height=8, dpi=300, bg="white")
cat("  wrote L3_item03_naive_comparisons.png\n")
