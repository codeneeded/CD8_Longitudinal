################################################################################
# FACT-seq with a second donor: what the figure needs
#
# CP018 is SUPPRESSED, CP003 UNSUPPRESSED, so adding it turns Figure 3 from a
# single-donor description into a two-donor contrast and the ART axis becomes
# the organising variable. This builds the panel that contrast needs, which the
# current figure has no equivalent of: the HIV-responsive fraction per donor per
# timepoint. Absolute counts hide it, because the two donors have different cell
# yields at different timepoints.
#
# Out: Reviewer_Comments/Figures/L3_factseq_responsive_fraction.png
#      Reviewer_Comments/Tables/L3_factseq_two_donor_summary.csv
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(patchwork)})
R   <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
FIG <- file.path(R,"Reviewer_Comments","Figures"); TAB <- file.path(R,"Reviewer_Comments","Tables")
o <- qs_read(file.path(R,"saved_R_data","CP003_CP018_integrated.qs2")); m <- o@meta.data
# alpha-beta CD8 cells only, defined once in 00.Joint_abCD8_filter.R (same cells as item 11).
source(file.path(R, "Scripts", "5. Exploratory", "00.Joint_abCD8_filter.R"))
m <- m[joint_keep(o), ]
src <- read.csv(file.path(FIG, "L3_jointUMAP_source_data.csv"), check.names = FALSE)
stopifnot(setequal(rownames(m), src$cell)); m <- m[src$cell, ]
m$HIV_Specific_TCR <- ifelse(src$resp, "HIV-Specific TCR", "Other")   # one rule for both donors, set in 21

s <- m %>% mutate(resp = HIV_Specific_TCR == "HIV-Specific TCR") %>%
  group_by(participant, ART_group, Timepoint, cellranger) %>%
  summarise(n_cells=n(), n_resp=sum(resp),
            n_clono=n_distinct(CTstrict[resp & !is.na(CTstrict)]), .groups="drop") %>%
  mutate(pct_resp = round(100*n_resp/n_cells, 1),
         visit = factor(ifelse(Timepoint == "2m", "2 months",
                      "Late\n(CP003 101m, CP018 90m)"),
                 levels = c("2 months", "Late\n(CP003 101m, CP018 90m)")),
         label = paste0(participant, "\n", ART_group))
cat("=== two-donor FACT-seq summary ===\n"); print(as.data.frame(s), row.names=FALSE)
write.csv(s, file.path(TAB,"L3_factseq_two_donor_summary.csv"), row.names=FALSE)

DON <- c(CP003="#1F77B4", CP018="#FF7F0E")
pA <- ggplot(s, aes(visit, pct_resp, group=participant, colour=participant)) +
  geom_line(linewidth=1.3) + geom_point(size=5) +
  geom_text(aes(label=sprintf("%.1f%%", pct_resp)), nudge_x=0.13, hjust=0, size=4.8, fontface="bold", show.legend=FALSE) +
  scale_colour_manual(values=DON, guide="none") +
  scale_y_continuous(limits=c(0, NA), expand=expansion(mult=c(0,.1))) +
  scale_x_discrete(expand=expansion(add=c(.4,.7))) +
  labs(title="HIV-responsive cells as a proportion of αβ CD8+ T cells",
       subtitle="CP003 unsuppressed, CP018 suppressed",
       x=NULL, y="% of αβ CD8+ T cells") +
  theme_classic(base_size=16) +
  theme(plot.title=element_text(face="bold", size=17),
        plot.subtitle=element_text(size=13, colour="grey25"),
        axis.text=element_text(colour="black", size=13))
pl <- s %>% select(participant, visit, n_resp, n_clono) %>%
  tidyr::pivot_longer(c(n_resp, n_clono), names_to="what", values_to="v") %>%
  mutate(what=recode(what, n_resp="HIV-responsive cells", n_clono="HIV-responsive clonotypes"))
pB <- ggplot(pl, aes(visit, v, fill=participant)) +
  geom_col(position=position_dodge(.8), width=.7) +
  geom_text(aes(label=v), position=position_dodge(.8), vjust=-0.35, size=4.2, fontface="bold") +
  facet_wrap(~what, scales="free_y") +
  scale_fill_manual(values=DON, name="Donor") +
  scale_y_continuous(expand=expansion(mult=c(0,.2))) +
  labs(title="Counts", x=NULL, y=NULL) +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold", size=16),
        strip.background=element_rect(fill="grey95", colour=NA),
        strip.text=element_text(face="bold", size=13),
        axis.text=element_text(colour="black", size=12))
ggsave(file.path(FIG,"L3_factseq_responsive_fraction.png"),
       (pA | pB) + plot_layout(widths=c(1,1.35), guides="collect"),
width=17, height=6.8, dpi=300, bg="white")
cat("  wrote L3_factseq_responsive_fraction.png\n")
