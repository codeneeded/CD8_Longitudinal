################################################################################
# Revision — what the NON-CD8 lineages show        (Reviewer 1, major point 2)
#
# The reviewer asks what is in the CD4, NK, B, monocyte and other compartments.
# This panel is those lineages only; CD8 is excluded so the comparison is not
# dominated by the compartment the paper already covers.
#
# Reads R1M2_10_nonCD8_lineage_brief.csv (written by the aggregation step) and
# the per-term table. Nothing new is computed.
#
# OUT  Reviewer_Comments/Figures/R1M2_11_nonCD8_lineage_brief.png
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})
source("config_revision.R")

b <- read.csv(file.path(REV_TAB, "R1M2_10_nonCD8_lineage_brief.csv"),
              check.names = FALSE)
stopifnot(nrow(b) > 0, !"CD8" %in% b$lineage)

long <- b %>%
  select(lineage, verdict, n_populations,
         infection = DE_infection_HEIvsHEU,
         viraemia  = DE_viraemia_suppVsUnsupp,
         exposure  = DE_exposure_HEUvsHUU) %>%
  pivot_longer(c(infection, viraemia, exposure),
               names_to = "contrast", values_to = "n_DE")

# pDC has no differential expression at all - it failed the cells-per-group
# threshold. Keep the row visible and label it, rather than dropping it and
# leaving the reader to wonder where pDC went.
long$contrast <- factor(long$contrast,
  levels = c("infection", "viraemia", "exposure"),
  labels = c("HEI vs HEU\n(infection)", "suppressed vs\nunsuppressed (viraemia)",
             "HEU vs HUU\n(exposure)"))
ord <- b %>% arrange(desc(DE_infection_HEIvsHEU)) %>% pull(lineage)
long$lineage <- factor(long$lineage, levels = rev(ord))

VC <- c("Of interest" = "#1B7837", "Interesting, does not fit" = "#7FBC41",
        "Little of interest" = "#BDBDBD", "Nothing of interest" = "#8C8C8C",
        "Not assessed" = "#E08214")
lab <- b %>% transmute(lineage = factor(lineage, levels = rev(ord)), verdict,
                       n_populations)

p <- ggplot(long, aes(x = n_DE, y = lineage, fill = contrast)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(data = lab, inherit.aes = FALSE,
            aes(x = -450, y = lineage, label = paste0(n_populations, " pop")),
            hjust = 1, size = 4.4, colour = "grey30") +
  geom_point(data = lab, inherit.aes = FALSE,
             aes(x = -1500, y = lineage, colour = verdict), size = 6) +
  scale_fill_manual(values = c("#4575B4", "#C0392B", "#F4A582"), name = NULL) +
  scale_colour_manual(values = VC, name = "assessment") +
  # Annotations sit at negative x so they fall outside the plotting area, but
  # the axis must not advertise negative counts - break only at real values and
  # let the annotations draw beyond the panel.
  scale_x_continuous(labels = comma, breaks = pretty(c(0, max(long$n_DE)), n = 4),
                     limits = c(-2100, max(long$n_DE) * 1.04),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(title = "Differential expression outside the CD8 compartment",
       subtitle = paste0("significant gene-by-population calls (MAST, ",
                         "Benjamini-Hochberg p < 0.05), summed over the ",
                         "populations in each lineage\n",
                         "pDC did not meet the minimum cells-per-group ",
                         "threshold in any contrast and was not assessed"),
       x = "significant gene-by-population calls", y = NULL) +
  theme_classic(base_size = 17) +
  theme(plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 13, colour = "grey25"),
        legend.position = "right",
        axis.text.y = element_text(face = "bold"),
        axis.line.x = element_line(),
        plot.margin = margin(12, 14, 12, 12))

mkdirs(REV_FIG)
ggsave(file.path(REV_FIG, "R1M2_11_nonCD8_lineage_brief.png"), plot = p,
       width = 15.5, height = 8.5, dpi = 300, bg = "white")

cat("lineages:", paste(ord, collapse = ", "), "\n")
cat("wrote R1M2_11_nonCD8_lineage_brief.png\n")
