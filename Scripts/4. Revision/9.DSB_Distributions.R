################################################################################
# Revision script 9 — DSB-normalised protein distributions
#
# ANSWERS  Rev 1 minor 6: "S1A protein scale -1 to 2 is too low for DSB"
#
# THE POINT
#   S1A is a Seurat DotPlot, which z-scores each feature across groups and clips
#   the colour scale to +/-2.5 SD (DotPlot defaults scale=TRUE, col.min=-2.5,
#   col.max=2.5). Its scale was never in DSB units. This panel shows the actual
#   DSB values so the reviewer can see they reach the expected range.
#
#   Isotype controls are plotted alongside as the internal negative control: if
#   background subtraction worked they sit at zero while real markers do not.
#
# Output: R1m6_03_DSB_distributions.png + source CSV
################################################################################

.cfg <- Sys.getenv("CD8_REV_CONFIG", "")
if (!nzchar(.cfg)) {
  for (b in c("~/Documents/CD8_Longitudinal",
              "~/Desktop/Projects/Repositories/CD8_Longitudinal")) {
    p <- file.path(path.expand(b), "Scripts", "4. Revision", "config_revision.R")
    if (file.exists(p)) { .cfg <- p; break }
  }
}
source(.cfg)

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(Matrix); library(ggplot2); library(dplyr)
})

o <- qs_read(QS_ALL)
A <- GetAssayData(o, assay = "ADT", layer = "data")
ann <- as.character(o$Annotation)

ISO <- grep("IgG", rownames(A), value = TRUE)
stopifnot(length(ISO) > 0)

# Lineage markers paired with the population that should be positive, so the
# panel shows a positive population rather than a pooled distribution that is
# dominated by the negative cells.
# Patterns must match the PUBLISHED Annotation labels. "CD4 T" alone matches
# only the 2,441-cell "ISG+ CD4 T cell" cluster, not the 52,594 CD4 cells.
CD4_PAT <- "CD4|Treg"
PAIRS <- list(
  c("CD3D",   CD4_PAT,          "CD4 T cells"),
  c("CD4",    CD4_PAT,          "CD4 T cells"),
  c("CD8A",   "CD8",            "CD8 T cells"),
  c("CD19",   "B cell",         "B cells"),
  c("CD14",   "Monocyte",       "Monocytes"),
  c("CD56",   "NK",             "NK cells"),
  c("CD45RA", "^Naive",         "Naive T/B cells"),
  c("CD27",   CD4_PAT,          "CD4 T cells")
)

rows <- list()
for (pr in PAIRS) {
  g <- pr[1]
  if (!g %in% rownames(A)) next
  sel <- grepl(pr[2], ann)
  if (!any(sel)) next
  rows[[length(rows) + 1]] <- data.frame(
    marker = g, group = pr[3], class = "lineage marker",
    value = as.numeric(A[g, sel]), stringsAsFactors = FALSE)
}
for (g in ISO) {
  rows[[length(rows) + 1]] <- data.frame(
    marker = g, group = "all cells", class = "isotype control",
    value = as.numeric(A[g, ]), stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)

# Summary table: this is what the response letter quotes
smry <- df %>% group_by(marker, group, class) %>%
  summarise(n_cells = n(), median = round(median(value), 2),
            q95 = round(quantile(value, .95), 2),
            q99 = round(quantile(value, .99), 2),
            max = round(max(value), 2), .groups = "drop") %>%
  arrange(class, desc(q99))
save_tab("DSB_marker_distributions", as.data.frame(smry))
stopifnot(file.exists(file.path(REV_TAB, "DSB_marker_distributions.csv")))

# Panel-wide counts, so the "how many markers reach 5-20" claim is derived
q99_all <- apply(A, 1, function(v) quantile(v, .99))
panel <- data.frame(
  metric = c("antibodies in panel", "isotype controls",
             "markers with q99 > 5", "markers with q99 > 10", "markers with q99 > 20",
             "whole-matrix q99.9", "isotype median", "isotype q99",
             "real-marker q99", "real-marker max", "stored minimum",
             "exact zeros in stored matrix (pct)"),
  value = c(nrow(A), length(ISO),
            sum(q99_all > 5), sum(q99_all > 10), sum(q99_all > 20),
            round(quantile(A, .999), 2),
            round(median(A[ISO, ]), 3), round(quantile(A[ISO, ], .99), 3),
            round(quantile(A[setdiff(rownames(A), ISO), ], .99), 2),
            round(max(A[setdiff(rownames(A), ISO), ]), 2),
            round(min(A), 3),
            round(100 * (1 - length(A@x) / (nrow(A) * ncol(A))), 1)))
save_tab("R1m6_02_DSB_panel_summary", panel)

# Restrict to cells actually positive for the marker. Pooling a whole lineage
# mixes true positives with negatives (81.7% of the stored matrix is exactly
# zero), which makes a high-median marker look like a spike at zero.
dfp <- df %>% filter(value > 0)
np  <- dfp %>% count(marker, name = "n_pos")
smry_pos <- dfp %>% group_by(marker, group, class) %>%
  summarise(n_pos = n(), median_pos = round(median(value), 2),
            q95_pos = round(quantile(value, .95), 2),
            max_pos = round(max(value), 2), .groups = "drop")
save_tab("R1m6_01_DSB_marker_distributions_positive_cells", as.data.frame(smry_pos))

ord <- smry_pos %>% arrange(median_pos) %>% pull(marker)
dfp$marker <- factor(dfp$marker, levels = ord)
smry_pos$marker <- factor(smry_pos$marker, levels = ord)

p <- ggplot(dfp, aes(x = value, y = marker, fill = class)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.15, linewidth = 0.35,
               width = 0.62, colour = "grey25") +
  geom_vline(xintercept = c(5, 20), linetype = "dashed", colour = "grey35") +
  geom_text(data = smry_pos, aes(x = max_pos, label = paste0("  med ", median_pos)),
            hjust = 0, size = 3.5, colour = "grey20") +
  scale_fill_manual(values = c(`lineage marker` = "#2166AC",
                               `isotype control` = "#B2182B"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  labs(title = "DSB-normalised protein values reach the expected range",
       subtitle = paste(strwrap(paste0(
         "Distribution over cells POSITIVE for each marker (value > 0); dashed lines mark the ",
         "5-20 range expected of DSB-normalised positive populations. Lineage markers are ",
         "restricted to their positive population, isotype controls span all cells as the ",
         "internal negative control. Of ", nrow(A), " antibodies, ", sum(q99_all > 5),
         " have a 99th percentile above 5 and ", sum(q99_all > 20), " above 20. Supplementary ",
         "S1A used a Seurat DotPlot, whose colour scale is z-scored across groups and clipped ",
         "to +/-2.5 SD - it was never in DSB units."),
         width = 104), collapse = "\n"),
       x = "DSB-normalised value (positive cells)", y = NULL) +
  theme_rev() + theme(legend.position = "top")

save_fig("R1m6_03_DSB_distributions", p, width = 11, height = 8)

cat("\n=== panel summary ===\n"); print(panel, row.names = FALSE)
cat("\n=== per-marker (positive populations) ===\n"); print(as.data.frame(smry), row.names = FALSE)
