################################################################################
# Revision — disaggregate the main results by inferred sex   (Editor point E-5)
#
# WHAT SAGER ASKS, AND WHAT IT DOES NOT
#   SAGER asks that sex be REPORTED and that results be disaggregated where the
#   data allow. It does not ask for an underpowered comparison. Those are
#   different obligations and this cohort can meet the first but not the second.
#
# WHY NOT THE SECOND — the cross-tabulation settles it, not a judgement call:
#     HEI  8 female / 3 male
#     HEU  2 female / 1 male
#     HUU  2 female / 0 male
#   One cell is EMPTY and another is n = 1. A sex x exposure-group comparison
#   over 16 participants cannot be estimated, and reporting a p value from it
#   would be worse than reporting nothing. So:
#     - the sex distribution is reported in full (the SAGER disclosure)
#     - a DESCRIPTIVE within-HEI split is given, with group sizes on the figure
#     - no test is run across exposure groups, and the response says why
#
# OUT
#   Tables/E5_06_composition_by_sex_descriptive.csv
#   Figures/E5_02_composition_by_sex.png
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(scales)
})
source("config_revision.R")

sx <- read.csv(file.path(REV_TAB, "E5_03_sex_inference_per_participant.csv"))
o  <- load_qs("cd8_annotated_scored")
ANN <- if ("CD8_Annotation_Rev" %in% names(o@meta.data)) "CD8_Annotation_Rev" else "CD8_Annotation"

md <- o@meta.data
md$participant <- sub("_.*$", "", md$orig.ident)
md <- left_join(md, sx[, c("participant", "sex_inferred")], by = "participant")
stopifnot(!any(is.na(md$sex_inferred)))

# per-participant composition: the unit of analysis is the participant, not the
# cell, so a large library cannot dominate the summary.
comp <- md %>%
  group_by(participant, Condition, sex_inferred, pop = .data[[ANN]]) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  group_by(participant, Condition, sex_inferred) %>%
  mutate(pct_of_CD8 = 100 * n / sum(n)) %>%
  ungroup()

cat("=== group sizes (participants) ===\n")
print(comp %>% distinct(participant, Condition, sex_inferred) %>%
        count(Condition, sex_inferred) %>% as.data.frame(), row.names = FALSE)

# ---- descriptive summary, HEI only ------------------------------------------
hei <- comp %>% filter(Condition == "HEI")
n_f <- n_distinct(hei$participant[hei$sex_inferred == "female"])
n_m <- n_distinct(hei$participant[hei$sex_inferred == "male"])
cat(sprintf("\nHEI only: %d female, %d male participants\n", n_f, n_m))

summ <- hei %>% group_by(pop, sex_inferred) %>%
  summarise(n_participants = n_distinct(participant),
            median_pct = round(median(pct_of_CD8), 2),
            min_pct = round(min(pct_of_CD8), 2),
            max_pct = round(max(pct_of_CD8), 2), .groups = "drop")

wide <- summ %>%
  select(pop, sex_inferred, median_pct) %>%
  pivot_wider(names_from = sex_inferred, values_from = median_pct,
              names_prefix = "median_pct_") %>%
  mutate(difference = round(median_pct_female - median_pct_male, 2))

cat("\n=== HEI composition by inferred sex — DESCRIPTIVE, no test ===\n")
print(as.data.frame(wide), row.names = FALSE)

out <- summ %>% left_join(wide %>% select(pop, difference), by = "pop") %>%
  mutate(comparison = "descriptive only",
         why_no_test = paste0("n = ", n_f, " female vs ", n_m,
                              " male HEI participants; sex x exposure-group ",
                              "cross-tabulation has an empty cell (HUU male = 0) ",
                              "and an n = 1 cell (HEU male), so no test is ",
                              "estimable across groups"))

mkdirs(REV_TAB); mkdirs(REV_FIG)
write.csv(out, file.path(REV_TAB, "E5_06_composition_by_sex_descriptive.csv"), row.names = FALSE)

# ---- figure -----------------------------------------------------------------
# Show every participant as a point. With n = 3 male the reader must see the
# individual values; a bar of medians would imply precision that is not there.
pl <- hei %>% mutate(pop = factor(pop, levels = rev(sort(unique(pop)))))
p <- ggplot(pl, aes(x = pct_of_CD8, y = pop, colour = sex_inferred)) +
  geom_point(size = 4, alpha = 0.85,
             position = position_dodge(width = 0.6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.4,
               position = position_dodge(width = 0.6), show.legend = FALSE,
               aes(colour = sex_inferred)) +
  scale_colour_manual(values = c(female = "#B2182B", male = "#2166AC"),
                      name = "Inferred sex",
                      labels = c(female = paste0("female (n=", n_f, ")"),
                                 male = paste0("male (n=", n_m, ")"))) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "CD8 composition by inferred genetic sex, HEI participants",
       subtitle = paste0("each point is one participant; bar = median. DESCRIPTIVE ONLY \u2014 ",
                         "no test is reported.\nThe sex \u00d7 exposure-group table has an empty ",
                         "cell (HUU male = 0) and an n = 1 cell (HEU male)."),
       x = "percentage of that participant's CD8 cells", y = NULL) +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 14, colour = "grey25"),
        legend.position = "top")

ggsave(file.path(REV_FIG, "E5_02_composition_by_sex.png"), plot = p,
       width = 14, height = 9.5, dpi = 300, bg = "white")

cat("\nwrote E5_06_composition_by_sex_descriptive.csv and E5_02_composition_by_sex.png\n")
cat("\nFOR THE RESPONSE: report the distribution and this descriptive split.\n",
    "Do NOT report a sex effect. State that the cohort cannot support a\n",
    "sex-stratified comparison and give the cross-tabulation as the reason.\n", sep = "")
