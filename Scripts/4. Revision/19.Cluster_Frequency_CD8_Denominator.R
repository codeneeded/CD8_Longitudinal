################################################################################
# Figures 1F / S2D and 1G / S2E with total CD8+ T cells as the denominator,
# with p values on the S2D box plots
#
# Published panels divide by total_pbmc (every cell in the sample). These divide
# by total CD8, defined as the sum of the six Figure 1C clusters in that sample -
# so each panel reads as "share of the CD8 compartment" rather than "share of
# PBMC". Both are computed side by side and the significance is compared, so the
# question "should the interpretation change" is answered from the numbers.
#
# Cluster names come out of Figure_1.R so they cannot drift.
# Out: Reviewer_Comments/Figures/Fig1F_S2D_cluster_frequency_CD8_denominator.png
#      Reviewer_Comments/Figures/Fig1G_S2E_VL_correlation_CD8_denominator.png
#      Reviewer_Comments/Tables/Fig1F_1G_denominator_comparison.csv
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(rstatix); library(ggpubr)
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

f1 <- readLines(file.path(R,"Scripts","3. Manuscript","Figure_1.R"), warn=FALSE)
grab <- function(name) {
  # Anchor on the assignment only. cd8_short_labels is written as
  # `setNames(` then `c(` on the NEXT line, so requiring c( on the same
  # line silently returns NA and the walker then fails on seq(NA, n).
  i <- grep(paste0("^", name, "\\s*<-"), f1)[1]
  if (is.na(i)) stop("cannot find ", name, " in Figure_1.R", call. = FALSE)
  depth <- 0L; out <- character()
  for (k in seq(i, length(f1))) {
    ln <- f1[k]; bare <- gsub('"[^"]*"', "", ln); out <- c(out, ln)
    depth <- depth + lengths(regmatches(bare, gregexpr("\\(", bare))) -
                     lengths(regmatches(bare, gregexpr("\\)", bare))); if (depth <= 0L) break
  }
  eval(parse(text = paste(out, collapse="\n")))
}
cd8_cluster_names <- grab("cd8_cluster_names"); cd8_short_labels <- grab("cd8_short_labels")

all <- qs_read(file.path(R,"saved_R_data","TARA_ALL_annotated_final.qs2"))
keep <- colnames(all)[ all$Condition=="HEU" | all$Condition=="HUU" |
                      (all$Condition=="HEI" & all$Timepoint_Group=="PreART_Entry") ]
sub <- subset(all, cells = keep)
sub$Condition <- factor(sub$Condition, levels=c("HUU","HEU","HEI"))
m <- as.data.frame(sub@meta.data)
cat("Figure 1 subset:", nrow(m), "cells\n")

tot_pbmc <- m %>% group_by(orig.ident, Condition) %>% summarise(total_pbmc=n(), .groups="drop")
cf <- m %>% filter(Annotation %in% cd8_cluster_names) %>%
  group_by(orig.ident, Condition, Annotation) %>% summarise(cluster_cells=n(), .groups="drop") %>%
  group_by(orig.ident, Condition) %>%
  complete(Annotation=cd8_cluster_names, fill=list(cluster_cells=0)) %>% ungroup() %>%
  left_join(tot_pbmc, by=c("orig.ident","Condition")) %>%
  group_by(orig.ident) %>% mutate(total_cd8 = sum(cluster_cells)) %>% ungroup() %>%
  mutate(pct_pbmc = 100*cluster_cells/total_pbmc,
         pct_cd8  = 100*cluster_cells/total_cd8,
         Facet_Label = factor(cd8_short_labels[as.character(Annotation)], levels=cd8_short_labels))

# ---- Figure 1F / S2D: composition by condition, CD8 denominator, with p values -----
st <- cf %>% group_by(Facet_Label) %>% wilcox_test(pct_cd8 ~ Condition) %>%
      adjust_pvalue(method="none") %>% add_significance("p") %>% ungroup()
# scales::pvalue() renders p > 0.999 as ">0.999", which with the "p = " prefix
# reads "p = >0.999". Build the label explicitly instead.
st <- st %>% mutate(plab = ifelse(p < 0.001, "p < 0.001",
                           ifelse(p > 0.999, "p > 0.999", paste0("p = ", signif(p, 2)))))
st <- st %>% add_xy_position(x="Condition", scales="free")
p6 <- ggplot(cf, aes(Condition, pct_cd8, fill=Condition)) +
  geom_boxplot(outlier.shape=NA, width=.65, alpha=.85) +
  geom_jitter(width=.16, size=2.4, shape=21, fill="white", colour="black", stroke=.5) +
  stat_pvalue_manual(st, label="plab", tip.length=.008, size=4.3,
                     bracket.nudge.y=0, step.increase=.055) +
  scale_y_continuous(expand = expansion(mult = c(.05, .22))) +
  facet_wrap(~Facet_Label, scales="free_y", nrow=1) +
  scale_fill_manual(values=c(HUU="#4C72B0", HEU="#55A868", HEI="#C44E52")) +
  labs(title="CD8 cluster frequency as a share of total CD8+ T cells",
       x=NULL, y="% of total CD8+ T cells") +
  theme_classic(base_size=16) +
  theme(plot.title=element_text(face="bold", size=17),
        strip.text=element_text(face="bold", size=14), legend.position="none",
        axis.text.x=element_text(size=13, colour="black"))
ggsave(file.path(FIG,"Fig1F_S2D_cluster_frequency_CD8_denominator.png"), p6,
       width=20, height=6.4, dpi=300, bg="white")
write.csv(cf, file.path(FIG,"Fig1F_S2D_cluster_frequency_CD8_denominator_source_data.csv"), row.names=FALSE)

# ---- Figure 1G / S2E: VL correlation at pre-ART HEI, both denominators ----------------
hei <- m %>% filter(Condition=="HEI")
vl <- hei %>% group_by(orig.ident) %>% summarise(VL=first(Viral_Load_num), .groups="drop") %>%
      mutate(log10_VL = log10(VL))
cv <- cf %>% filter(Condition=="HEI") %>% left_join(vl, by="orig.ident") %>% filter(!is.na(log10_VL))
cor_cd8  <- cv %>% group_by(Facet_Label) %>% cor_test(log10_VL, pct_cd8,  method="spearman") %>% ungroup()
cor_pbmc <- cv %>% group_by(Facet_Label) %>% cor_test(log10_VL, pct_pbmc, method="spearman") %>% ungroup()
lab <- cor_cd8 %>% mutate(stat_label=paste0("rho = ", sprintf("%.2f",cor), "\np = ", signif(p,2)),
                          x=-Inf, y=Inf)
p7 <- ggplot(cv, aes(log10_VL, pct_cd8)) +
  geom_smooth(method="lm", se=TRUE, linewidth=.8, colour="grey40", fill="grey85", alpha=.3) +
  geom_point(shape=21, size=4.2, stroke=.5, fill="#C44E52", colour="white") +
  geom_text(data=lab, aes(x=x, y=y, label=stat_label), inherit.aes=FALSE,
            hjust=-0.06, vjust=1.2, size=5.6, lineheight=1.2, fontface="bold", colour="grey15") +
  facet_wrap(~Facet_Label, scales="free", nrow=1) +
  labs(title="CD8 cluster frequency against pre-ART viral load, CD8 denominator",
       x=expression(log[10]~"viral load"), y="% of total CD8+ T cells") +
  theme_classic(base_size=16) +
  theme(plot.title=element_text(face="bold", size=17),
        strip.text=element_text(face="bold", size=14))
ggsave(file.path(FIG,"Fig1G_S2E_VL_correlation_CD8_denominator.png"), p7,
       width=22, height=6.2, dpi=300, bg="white")
write.csv(cv, file.path(FIG,"Fig1G_S2E_VL_correlation_CD8_denominator_source_data.csv"), row.names=FALSE)

# ---- the comparison that answers "should the interpretation change" ----------
cmp <- cor_pbmc %>% select(Facet_Label, rho_PBMC=cor, p_PBMC=p) %>%
  left_join(cor_cd8 %>% select(Facet_Label, rho_CD8=cor, p_CD8=p), by="Facet_Label") %>%
  mutate(sig_PBMC = p_PBMC < .05, sig_CD8 = p_CD8 < .05,
         changes  = ifelse(sig_PBMC != sig_CD8, "YES - significance flips", "no"))
cat("\n=== ITEM 7: viral-load correlation, published PBMC denominator vs CD8 denominator ===\n")
print(as.data.frame(cmp %>% mutate(across(where(is.numeric), ~round(.,4)))), row.names=FALSE)
cat("\n=== ITEM 6: condition comparison under the CD8 denominator ===\n")
print(as.data.frame(st %>% select(Facet_Label, group1, group2, p) %>% mutate(p=signif(p,3))), row.names=FALSE)
st_pbmc <- cf %>% group_by(Facet_Label) %>% wilcox_test(pct_pbmc ~ Condition) %>% ungroup()
cat("\n=== ITEM 6: the same comparison under the PUBLISHED PBMC denominator ===\n")
print(as.data.frame(st_pbmc %>% select(Facet_Label, group1, group2, p) %>% mutate(p=signif(p,3))), row.names=FALSE)
out <- st %>% select(Facet_Label, group1, group2, p_CD8=p) %>%
  left_join(st_pbmc %>% select(Facet_Label, group1, group2, p_PBMC=p), by=c("Facet_Label","group1","group2")) %>%
  mutate(sig_CD8=p_CD8<.05, sig_PBMC=p_PBMC<.05,
         changes=ifelse(sig_CD8!=sig_PBMC,"YES - significance flips","no"))
write.csv(bind_rows(
  cmp %>% mutate(comparison="VL correlation (pre-ART HEI)") %>%
    select(comparison, cluster=Facet_Label, stat_PBMC=rho_PBMC, p_PBMC, stat_CD8=rho_CD8, p_CD8, changes),
  out %>% mutate(comparison=paste0("condition ", group1, " vs ", group2),
                 stat_PBMC=NA_real_, stat_CD8=NA_real_) %>%
    select(comparison, cluster=Facet_Label, stat_PBMC, p_PBMC, stat_CD8, p_CD8, changes)),
  file.path(TAB,"Fig1F_1G_denominator_comparison.csv"), row.names=FALSE)
cat("\n=== condition comparisons where significance flips ===\n")
print(as.data.frame(out %>% filter(changes!="no") %>% select(Facet_Label,group1,group2,p_PBMC,p_CD8)), row.names=FALSE)
cat("\nwrote both panels, their source data, and Fig1F_1G_denominator_comparison.csv\n")
