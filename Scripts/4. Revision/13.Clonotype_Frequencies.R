################################################################################
# Revision — clonotype dynamics as FREQUENCIES rather than counts
#                                                   (Reviewer 2, major 2b)
#
# THE COMMENT
#   "Clonotype dynamics should be shown as frequencies rather than counts."
#
# WHY THE REVIEWER IS RIGHT
#   A count of cells in expanded clonotypes rises and falls with how many cells
#   were recovered from that sample, and recovery varies several-fold across
#   these libraries. Two timepoints with identical repertoire structure can
#   therefore show different counts, and a real change can be hidden by a
#   compensating change in depth. The frequency - cells in expanded clonotypes
#   as a proportion of cells with a recovered TCR - removes that dependence.
#
#   The DENOMINATOR matters and is easy to get wrong. It must be cells with a
#   recovered TCR, not all cells: dividing by all cells scales each timepoint by
#   its TCR recovery rate and reintroduces exactly the artefact being removed.
#
# OUT
#   Tables/R2M2b_01_clonotype_frequency_by_timepoint.csv
#   Figures/R2M2b_02_clonotype_frequency.png
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(scales); library(ggrepel)
})
source("config_revision.R")

o  <- load_qs("cd8_annotated_scored")
md <- o@meta.data
md$.cell <- rownames(md)

ctcol <- grep("^CTstrict$|^CTaa$|^CTnt$", names(md), value = TRUE)[1]
stopifnot(!is.na(ctcol))
cat("clonotype column:", ctcol, "\n")

md$.ct <- as.character(md[[ctcol]])
md$.hastcr <- !is.na(md$.ct) & nzchar(md$.ct)
cat("TCR-recovered cells:", sum(md$.hastcr), "of", nrow(md), "\n")

md$participant <- sub("_.*$", "", md$orig.ident)
md$timepoint   <- sub("^[^_]*_", "", md$orig.ident)

# Expanded = clonotype observed in more than one cell WITHIN A SAMPLE, the
# definition stated in Methods (see R2-M2c/2d). Pooling across samples would
# call a clonotype expanded on the strength of two different draws.
tc <- md %>% filter(.hastcr) %>% count(orig.ident, .ct, name = "n_cells")
exp_ct <- tc %>% filter(n_cells > 1)

per_sample <- md %>% filter(.hastcr) %>%
  group_by(participant, timepoint, orig.ident) %>%
  summarise(n_tcr = n(), .groups = "drop") %>%
  left_join(exp_ct %>% group_by(orig.ident) %>%
              summarise(n_expanded_clonotypes = n(),
                        n_cells_in_expanded = sum(n_cells), .groups = "drop"),
            by = "orig.ident") %>%
  mutate(across(c(n_expanded_clonotypes, n_cells_in_expanded), ~tidyr::replace_na(., 0)),
         pct_cells_in_expanded = round(100 * n_cells_in_expanded / n_tcr, 2),
         pct_clonotypes_expanded = round(
           100 * n_expanded_clonotypes /
             sapply(orig.ident, function(s) n_distinct(tc$.ct[tc$orig.ident == s])), 2))

# order timepoints numerically, not alphabetically (1m, 2m, 12m, 24m, 42m, 44m)
tp_num <- function(x) as.numeric(gsub("[^0-9.]", "", x))
per_sample <- per_sample %>% mutate(tp_months = tp_num(timepoint)) %>%
  arrange(participant, tp_months)

cat("\n=== frequency of cells in expanded clonotypes, per sample ===\n")
print(as.data.frame(per_sample %>%
  select(participant, timepoint, n_tcr, n_cells_in_expanded,
         pct_cells_in_expanded, pct_clonotypes_expanded)), row.names = FALSE)

# The point of the comment: does the count mislead where the frequency does not?
cat("\n=== count vs frequency: do they rank samples differently? ===\n")
cmp <- per_sample %>%
  mutate(rank_count = rank(-n_cells_in_expanded),
         rank_freq  = rank(-pct_cells_in_expanded),
         rank_shift = rank_count - rank_freq) %>%
  arrange(desc(abs(rank_shift)))
print(as.data.frame(cmp %>% select(participant, timepoint, n_tcr,
        n_cells_in_expanded, pct_cells_in_expanded, rank_count, rank_freq,
        rank_shift) %>% head(8)), row.names = FALSE)
cat(sprintf("\nTCR-recovered cells range %d-%d across samples (%.1f-fold), which is why\n",
            min(per_sample$n_tcr), max(per_sample$n_tcr),
            max(per_sample$n_tcr) / min(per_sample$n_tcr)))
cat("the count and the frequency do not rank samples the same way.\n")

mkdirs(REV_TAB); mkdirs(REV_FIG)
write.csv(per_sample, file.path(REV_TAB, "R2M2b_01_clonotype_frequency_by_timepoint.csv"),
          row.names = FALSE)

# ---- figure: longitudinal participants only, frequency on the y axis -------
lon <- per_sample %>% group_by(participant) %>% filter(n() > 1) %>% ungroup()
p <- ggplot(lon, aes(x = tp_months, y = pct_cells_in_expanded,
                     colour = participant, group = participant)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 4) +
  # Repelled, not offset by a fixed vjust: several participants share a
  # timepoint at a similar frequency and fixed offsets overprinted each other
  # (29.2% over 29.5%, 1.0% over 0.0%).
  ggrepel::geom_text_repel(aes(label = sprintf("%.1f%%", pct_cells_in_expanded)),
            size = 4.3, show.legend = FALSE, min.segment.length = 0,
            box.padding = 0.7, point.padding = 0.45, seed = 1,
            max.overlaps = Inf, force = 5, segment.size = 0.3,
            segment.colour = "grey55") +
  # Sampling ages cluster in pairs (24/25 and 42/44 months), so labelling every
  # distinct value overprinted the pairs. Ticks name the visit; the exact age
  # stays in the source table.
  scale_x_continuous(breaks = c(1, 12, 24.5, 43),
                     labels = c("1m", "12m", "24-25m", "42-44m"),
                     minor_breaks = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, max(lon$pct_cells_in_expanded) * 1.15)) +
  labs(title = "Cells in expanded clonotypes, as a percentage of TCR-recovered cells",
       subtitle = paste0("cells in expanded clonotypes as a percentage of cells with a ",
                         "recovered TCR\ndenominator is TCR-recovered cells per sample, ",
                         "which range ", min(per_sample$n_tcr), "-",
                         max(per_sample$n_tcr), " (",
                         sprintf("%.1f", max(per_sample$n_tcr)/min(per_sample$n_tcr)),
                         "-fold)"),
       x = "age at sampling", y = "cells in expanded clonotypes",
       colour = "Participant") +
  theme_classic(base_size = 19) +
  theme(plot.title = element_text(face = "bold", size = 21),
        plot.subtitle = element_text(size = 14.5, colour = "grey25"))

ggsave(file.path(REV_FIG, "R2M2b_02_clonotype_frequency.png"), plot = p,
       width = 14, height = 9, dpi = 300, bg = "white")
cat("\nwrote R2M2b_01_clonotype_frequency_by_timepoint.csv and",
    "R2M2b_02_clonotype_frequency.png\n")
