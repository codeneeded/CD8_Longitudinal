################################################################################
# Revision script 11 - the two composition panels that were still promised
#
#   Rev 2 major 1c : Figure 1E shows only TEMRA/CTL, the one population that
#                    reached significance. Show EVERY CD8 cluster's proportion
#                    by exposure group so the non-significant comparisons are
#                    visible and the panel cannot read as selective.
#   Rev 2 minor 15 : each CD8 subset's proportion of the compartment at every
#                    sampled timepoint, as an alluvial.
#
# Both are computable from the saved objects - no new data, no re-run of the
# manuscript figure scripts (which cannot run here: they hardcode
# ~/Documents/CD8_Longitudinal and need scRepertoire/Trex/monocle3).
#
# COHORT SCOPE MATTERS. Figure 1E is drawn on the ENTRY-visit subset, exactly
# as Figure_1.R does:
#     Condition == "HEU" | Condition == "HUU" |
#     (Condition == "HEI" & Timepoint_Group == "PreART_Entry")
# That is 79,705 cells from 16 participants. Summing the whole object instead
# gives 126,435 across 26 samples and silently mixes in follow-up visits.
################################################################################

source("/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/4. Revision/config_revision.R")
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(ggalluvial)
})

################################################################################
# Rev 2 major 1c - every CD8 cluster by exposure group
################################################################################
o <- qs_read(QS_ALL)
md <- o@meta.data
entry <- md[md$Condition %in% c("HEU", "HUU") |
            (md$Condition == "HEI" & md$Timepoint_Group == "PreART_Entry"), ]
stopifnot(nrow(entry) > 0)
message("entry subset: ", nrow(entry), " cells, ",
        length(unique(entry$PID)), " participants")

cd8_labs <- grep("CD8", unique(as.character(entry$Annotation)), value = TRUE)

# Proportion is computed PER PARTICIPANT and then compared between groups: a
# pooled proportion would be dominated by whichever participant contributed
# most cells.
# droplevels() is essential: Annotation is a factor over all 31 PBMC
# populations, and complete() re-expands to every LEVEL, not every observed
# value - which silently added 24 non-CD8 rows of zeros and made the BH
# correction span 31 tests instead of the 7 CD8 populations actually compared.
comp <- entry %>%
  filter(Annotation %in% cd8_labs) %>%
  mutate(Annotation = droplevels(factor(Annotation))) %>%
  count(PID, Condition, Annotation, name = "n") %>%
  group_by(PID) %>%
  mutate(pct_of_CD8 = 100 * n / sum(n)) %>%
  ungroup() %>%
  left_join(entry %>% count(PID, name = "n_pbmc"), by = "PID") %>%
  mutate(pct_of_PBMC = 100 * n / n_pbmc) %>%
  complete(nesting(PID, Condition), Annotation,
           fill = list(n = 0, pct_of_CD8 = 0, pct_of_PBMC = 0))
stopifnot(nlevels(comp$Annotation) == length(cd8_labs))

# The published Figure 1E reports PAIRWISE Wilcoxon on percentage of total
# PBMC (HEI vs HEU p = 0.022, HEI vs HUU p = 0.026), unadjusted. Use the same
# statistic so the reviewer can place our numbers against the published ones,
# and report the BH-adjusted value alongside rather than instead - with 7
# populations tested the adjustment matters and hiding it would be worse than
# showing that the effect does not survive it.
pw <- function(df, a, b) {
  x <- df$pct_of_PBMC[df$Condition == a]; y <- df$pct_of_PBMC[df$Condition == b]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y)$p.value)
}
stats <- comp %>%
  group_by(Annotation) %>%
  group_modify(~ tibble(
    n_participants = n_distinct(.x$PID),
    p_HEI_vs_HEU   = pw(.x, "HEI", "HEU"),
    p_HEI_vs_HUU   = pw(.x, "HEI", "HUU"),
    p_kruskal_3grp = tryCatch(
      kruskal.test(pct_of_PBMC ~ factor(Condition), data = .x)$p.value,
      error = function(e) NA_real_))) %>%
  ungroup() %>%
  mutate(q_HEI_vs_HEU = p.adjust(p_HEI_vs_HEU, method = "BH"),
         q_HEI_vs_HUU = p.adjust(p_HEI_vs_HUU, method = "BH"))

lab <- stats %>%
  transmute(Annotation,
            tag = ifelse(is.na(p_HEI_vs_HEU), "n.t.",
                  sprintf("HEU p=%.3f\nHUU p=%.3f", p_HEI_vs_HEU, p_HEI_vs_HUU)))

comp$Condition <- factor(comp$Condition, levels = c("HUU", "HEU", "HEI"))

p1c <- ggplot(comp, aes(x = Condition, y = pct_of_PBMC, fill = Condition)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.30, linewidth = 0.6) +
  geom_jitter(width = 0.14, size = 1.9, alpha = 0.85, colour = "grey20") +
  geom_text(data = lab, aes(x = 1.5, y = Inf, label = tag),
            inherit.aes = FALSE, vjust = 1.15, size = 3.0, lineheight = 0.95) +
  facet_wrap(~ Annotation, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(HUU = "#7FB069", HEU = "#F4A259", HEI = "#C1503F")) +
  labs(title = "CD8 population frequencies by HIV exposure group",
       subtitle = paste(strwrap(paste0(
         "Each point is one participant's percentage of their own total PBMC, ",
         "the same denominator as the published Figure 1E (entry visit; ",
         nrow(entry), " PBMCs from ", length(unique(entry$PID)),
         " participants). p values are pairwise Wilcoxon against HEI, again ",
         "matching the published panel; Benjamini-Hochberg adjusted values ",
         "across all ", nrow(stats), " populations are in the source table, ",
         "where no population reaches q < 0.05."),
         width = 104), collapse = "\n"),
       x = NULL, y = "% of that participant's total PBMC") +
  theme_rev() +
  theme(legend.position = "none")

save_fig("R2M1c_01_all_CD8_clusters_by_group", p1c, width = 13, height = 8)
save_tab("R2M1c_01_all_CD8_clusters_by_group", comp)
save_tab("RevFig_1E_all_CD8_group_stats", stats)
cat("\n=== Rev 2 major 1c: per-cluster group comparison ===\n")
print(as.data.frame(stats), row.names = FALSE)

rm(o); invisible(gc(verbose = FALSE))

################################################################################
# Rev 2 minor 15 - CD8 subset proportion at each timepoint, alluvial
################################################################################
c8 <- qs_read(QS_CD8)
mc <- c8@meta.data
stopifnot(all(c("CD8_Annotation", "Timepoint_Group", "PID") %in% names(mc)))

TP_ORDER <- c("PreART_Entry", "PostART_Unsuppressed", "PostART_Suppressed")
TP_LABEL <- c(PreART_Entry = "Pre-ART\n(entry)",
              PostART_Unsuppressed = "On ART\n(unsuppressed)",
              PostART_Suppressed = "On ART\n(suppressed)")

# Proportion within each timepoint, averaged over participants so a single
# deeply-sequenced sample cannot set the shape of the ribbon.
prop <- mc %>%
  count(PID, Timepoint_Group, CD8_Annotation, name = "n") %>%
  group_by(PID, Timepoint_Group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  complete(nesting(PID, Timepoint_Group), CD8_Annotation,
           fill = list(n = 0, pct = 0)) %>%
  group_by(Timepoint_Group, CD8_Annotation) %>%
  summarise(pct = mean(pct), n_participants = n_distinct(PID), .groups = "drop") %>%
  filter(Timepoint_Group %in% TP_ORDER) %>%
  mutate(tp = factor(TP_LABEL[Timepoint_Group],
                     levels = unname(TP_LABEL[TP_ORDER])))

# Renormalise so each timepoint sums to 100 after averaging.
prop <- prop %>% group_by(tp) %>% mutate(pct = 100 * pct / sum(pct)) %>% ungroup()

cols <- if (exists("CD8_COLORS")) CD8_COLORS[levels(factor(prop$CD8_Annotation))] else NULL
# The strata must be NARROW and the flows drawn UNDER them, or the ribbons are
# hidden: with the same population order at every stage the flows run
# horizontally and a wide stratum covers them completely. Ordering strata by
# size within each stage (decreasing) also makes the crossings visible, which
# is the whole point of an alluvial rather than three stacked bars.
p15 <- ggplot(prop, aes(x = tp, y = pct, stratum = CD8_Annotation,
                        alluvium = CD8_Annotation, fill = CD8_Annotation)) +
  # colour = NA on geom_flow makes the ribbons vanish entirely - the layer draws
  # but with no visible ink, leaving what looks like three stacked bars. Omit it
  # and let the fill scale carry the ribbon; use a faint border instead.
  geom_flow(alpha = 0.45, width = 0.16, decreasing = TRUE,
            colour = "grey85", linewidth = 0.15) +
  geom_stratum(width = 0.16, colour = "white", linewidth = 0.3,
               decreasing = TRUE) +
  { if (!is.null(cols) && !any(is.na(cols)))
      scale_fill_manual(values = cols, name = NULL) else
      scale_fill_brewer(palette = "Set3", name = NULL) } +
  labs(title = "CD8 subset composition across the sampled course",
       subtitle = paste(strwrap(paste0(
         "Each band is one CD8 population's share of the compartment at that ",
         "timepoint, averaged across participants so a single deeply sequenced ",
         "sample cannot set the ribbon shape. Timepoint groups are treatment ",
         "stages, not equally spaced times, and participants contribute to ",
         "different stages - so read the panel as composition per stage rather ",
         "than as a rate of change."), width = 104), collapse = "\n"),
       x = NULL, y = "% of the CD8 compartment") +
  theme_rev() +
  theme(legend.position = "right")

save_fig("R2m15_01_S_CD8_composition_alluvial", p15, width = 12, height = 8)
save_tab("R2m15_01_S_CD8_composition_alluvial", prop)
cat("\n=== Rev 2 minor 15: composition per treatment stage ===\n")
print(as.data.frame(prop %>% arrange(tp, desc(pct))), row.names = FALSE)

cat("\nFigures: ", REV_FIG, "\nTables:  ", REV_TAB, "\n")
