################################################################################
# Revision — pathway enrichment across all lineages    (Reviewer 1, major 2)
#
# Reads the table written by 14a.Lineage_Pathway_Aggregate.py, which collects
# the per-lineage enrichment that was already run. No new statistics here.
#
# OUT  Reviewer_Comments/Figures/R1M2_07_lineage_pathway_heatmap.png
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})
source("config_revision.R")

f <- file.path(REV_TAB, "R1M2_07_lineage_pathway_summary.csv")
stopifnot(file.exists(f))
d <- read.csv(f, check.names = FALSE)

# Two contrasts, because they separate two different questions and give
# different answers: HEI vs HEU before treatment is the effect of infection,
# suppressed vs unsuppressed on treatment is the effect of ongoing viraemia.
KEYS <- c(HEIvsHEU_PreART          = "HEI vs HEU, pre-ART\n(effect of infection)",
          Suppressed_vs_Unsuppressed = "suppressed vs unsuppressed, on ART\n(effect of viraemia)")
d <- d %>% filter(contrast %in% names(KEYS),
                  collection == "MSigDB_Hallmark_2020")
stopifnot(nrow(d) > 0, n_distinct(d$contrast) == 2)

# One value per family x term: the strongest enrichment in that family.
m <- d %>% group_by(contrast, family, term) %>%
  summarise(padj = min(padj), n_pops = n_distinct(lineage), .groups = "drop") %>%
  mutate(neglog10 = -log10(padj))

# Keep terms seen in at least two families, so the panel shows shared
# programmes rather than one-off hits. Single-family terms stay in the CSV.
# Keep a term if it spans at least three families in either contrast, so the
# panel shows programmes that are genuinely system-wide.
keep <- m %>% group_by(contrast, term) %>% summarise(nf = n(), .groups = "drop") %>%
  group_by(term) %>% summarise(mx = max(nf), .groups = "drop") %>%
  filter(mx >= 3) %>% pull(term)
m <- m %>% filter(term %in% keep)

ford <- c("CD4", "CD8", "gd T", "NK", "B cell", "Monocyte", "MAIT", "pDC")
m$family <- factor(m$family, levels = intersect(ford, unique(m$family)))
tord <- m %>% group_by(term) %>% summarise(s = max(neglog10)) %>%
  arrange(s) %>% pull(term)
m$term <- factor(m$term, levels = tord)
m$contrast_lab <- factor(KEYS[m$contrast], levels = KEYS)

p <- ggplot(m, aes(x = family, y = term, fill = neglog10)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = n_pops), size = 4.6, colour = "grey15") +
  scale_fill_gradient(low = "#FDE0D9", high = "#A11C1C",
                      name = expression(-log[10]~italic(p)[adj])) +
  facet_wrap(~ contrast_lab, nrow = 1) +
  labs(title = "Infection and ongoing viraemia drive different system-wide programmes",
       subtitle = paste0("MSigDB Hallmark, Benjamini-Hochberg p < 0.05, across all 30 ",
                         "annotated populations\n",
                         "number in each cell = populations in that family reaching ",
                         "significance; blank = not significant"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 17) +
  theme(plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 13.5, colour = "grey25"),
        panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold"),
        legend.position = "right")

mkdirs(REV_FIG)
ggsave(file.path(REV_FIG, "R1M2_07_lineage_pathway_heatmap.png"), plot = p,
       width = 19, height = 10, dpi = 300, bg = "white")

cat("contrasts:", paste(levels(m$contrast_lab), collapse = " | "), "\n")
cat("families:", paste(levels(m$family), collapse = ", "), "\n")
cat("shared terms plotted:", length(keep), "\n")
cat("wrote R1M2_07_lineage_pathway_heatmap.png\n")
