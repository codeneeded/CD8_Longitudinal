################################################################################
# Revision — infer participant sex from the transcriptome  (Editor point E-5)
#
# THE ASK
#   "Follow SAGER guidance and disaggregate results by sex and gender where
#    these data were collected."
#
# WHY THIS SCRIPT EXISTS
#   Sex was never brought into the analysis data: no sex or gender field exists
#   in any saved object or table. The honest first answer was that the field is
#   absent. But sex-chromosome gene expression is IN the data, and single-cell
#   RNA-seq calls genetic sex essentially unambiguously, so the editor's ask is
#   answerable without the clinical record.
#
#   Note the distinction the SAGER guidance itself draws: this infers SEX from
#   sex-chromosome expression. It says nothing about gender, which cannot be
#   inferred from a transcriptome and would have to come from the participants.
#   Both the response and the manuscript must say "genetic sex, inferred" and
#   not "sex and gender".
#
# THE INTERNAL CONTROL THAT MAKES THIS TRUSTWORTHY
#   Participants were sampled repeatedly. Each participant's libraries are
#   independent measurements of one person, so every timepoint MUST return the
#   same call. Disagreement within a participant would mean the method is
#   unreliable here (or that samples are mislabelled) - either way it must stop
#   the analysis rather than be averaged away. That check is an assertion, not
#   a diagnostic printout.
#
# METHOD
#   Per library, the fraction of cells with any detected counts for
#     XIST                                    (X-inactivation, female-specific)
#     RPS4Y1 DDX3Y UTY EIF1AY KDM5D USP9Y NLGN4Y ZFY   (Y-linked)
#   A fraction, not a mean, because these are sparse counts and a fraction is
#   robust to depth differences between libraries; both series are included and
#   they differ ~5x in depth.
#
# OUT
#   Reviewer_Comments/Tables/E5_03_sex_inference_per_participant.csv
#   Reviewer_Comments/Figures/E5_01_sex_inference.png
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(scales)
})
source("config_revision.R")

XIST <- "XIST"
YGEN <- c("RPS4Y1", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y", "NLGN4Y", "ZFY")

o <- qs_read(QS_ALL)
DefaultAssay(o) <- "RNA"
if (length(Layers(o[["RNA"]])) > 3) o <- JoinLayers(o, assay = "RNA")

need <- c(XIST, YGEN)
have <- intersect(need, rownames(o[["RNA"]]))
stopifnot(XIST %in% have, length(intersect(YGEN, have)) >= 4)
cat("genes used:", paste(have, collapse = ", "), "\n")

cnt <- GetAssayData(o, assay = "RNA", layer = "counts")[have, , drop = FALSE]
md  <- o@meta.data
md$.cell <- colnames(o)

# per-library detection fractions
lib <- unique(md$orig.ident)
rows <- lapply(lib, function(s) {
  cells <- md$.cell[md$orig.ident == s]
  sub <- cnt[, cells, drop = FALSE]
  yg <- intersect(YGEN, rownames(sub))
  data.frame(
    library   = s,
    n_cells   = length(cells),
    xist_frac = Matrix::rowSums(sub[XIST, , drop = FALSE] > 0) / length(cells),
    y_frac    = sum(Matrix::colSums(sub[yg, , drop = FALSE]) > 0) / length(cells),
    row.names = NULL
  )
})
tb <- bind_rows(rows)

# Per-gene Y detection as well, so a single dropout gene cannot drive a call.
ydet <- sapply(intersect(YGEN, rownames(cnt)), function(g)
  sapply(lib, function(s) {
    cells <- md$.cell[md$orig.ident == s]
    sum(cnt[g, cells] > 0) / length(cells)
  }))
ydet <- as.data.frame(ydet); ydet$library <- rownames(ydet)
tb <- left_join(tb, ydet, by = "library")

# ---- the call ----------------------------------------------------------------
# Female: XIST widely detected, Y essentially absent. Male: the reverse.
# The margin is the separation between the two fractions; a small margin means
# the library cannot be called and must be reported as such rather than forced.
tb <- tb %>%
  mutate(margin = abs(xist_frac - y_frac),
         sex_call = case_when(
           xist_frac > y_frac & margin >= 0.10 ~ "female",
           y_frac > xist_frac & margin >= 0.10 ~ "male",
           TRUE                                ~ "indeterminate"))

# participant id: strip the timepoint suffix
tb$participant <- sub("_.*$", "", tb$library)

cat("\n=== per-library calls ===\n")
print(tb %>% select(library, n_cells, xist_frac, y_frac, margin, sex_call) %>%
        mutate(across(c(xist_frac, y_frac, margin), ~round(., 3))) %>%
        arrange(participant <- sub("_.*$", "", library)), row.names = FALSE)

# ---- the internal control ---------------------------------------------------
per_part <- tb %>% group_by(participant) %>%
  summarise(n_lib = n(),
            calls = paste(sort(unique(sex_call)), collapse = "/"),
            consistent = n_distinct(sex_call[sex_call != "indeterminate"]) <= 1,
            .groups = "drop")

cat("\n=== per-participant consistency across timepoints ===\n")
print(as.data.frame(per_part), row.names = FALSE)

bad <- per_part %>% filter(!consistent)
if (nrow(bad)) {
  cat("\n** INCONSISTENT participants:\n"); print(as.data.frame(bad), row.names = FALSE)
  stop("A participant returned different sex calls at different timepoints. ",
       "Do not report these calls: either the method is unreliable in this ",
       "dataset or sample labels are crossed. Resolve before using E-5.")
}
cat("\nall participants internally consistent across timepoints\n")

final <- tb %>% group_by(participant) %>%
  summarise(n_libraries = n(),
            sex_inferred = {
              v <- sex_call[sex_call != "indeterminate"]
              if (!length(v)) "indeterminate" else unique(v)
            },
            min_margin = round(min(margin), 3),
            mean_xist_frac = round(mean(xist_frac), 3),
            mean_y_frac = round(mean(y_frac), 3),
            .groups = "drop")

cat("\n=== participant-level result ===\n")
print(as.data.frame(final), row.names = FALSE)
cat("\ndistribution:\n"); print(table(final$sex_inferred))

mkdirs(REV_TAB); mkdirs(REV_FIG)
write.csv(tb,    file.path(REV_TAB, "E5_04_sex_inference_per_library.csv"), row.names = FALSE)
write.csv(final, file.path(REV_TAB, "E5_03_sex_inference_per_participant.csv"), row.names = FALSE)

# ---- figure -----------------------------------------------------------------
# The point of the panel is separation: two clouds with nothing between them is
# what makes the call defensible, so plot the two fractions against each other
# rather than a bar of the result.
p <- ggplot(tb, aes(x = xist_frac, y = y_frac, colour = sex_call)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
  geom_point(size = 4.6, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = library), size = 4.6, max.overlaps = 40,
                           show.legend = FALSE, min.segment.length = 0.2) +
  scale_colour_manual(values = c(female = "#B2182B", male = "#2166AC",
                                 indeterminate = "grey50"), name = "Inferred sex") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(tb$xist_frac, tb$y_frac) * 1.12)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(tb$xist_frac, tb$y_frac) * 1.12)) +
  labs(title = "Genetic sex inferred from sex-chromosome expression",
       subtitle = paste0("each point is one library; ", nrow(final),
                         " participants, ", nrow(tb), " libraries\n",
                         "dashed line = equal detection; distance from it is the call margin"),
       x = "cells detecting XIST",
       y = "cells detecting any Y-linked gene") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 22),
        plot.subtitle = element_text(size = 16, colour = "grey25"),
        legend.position = "right")

ggsave(file.path(REV_FIG, "E5_01_sex_inference.png"), plot = p,
       width = 13, height = 9.5, dpi = 300, bg = "white")

cat("\nwrote:\n  E5_03_sex_inference_per_participant.csv\n",
    "  E5_04_sex_inference_per_library.csv\n  E5_01_sex_inference.png\n", sep = "")
cat("\nNOTE for the response: this is GENETIC SEX from sex-chromosome\n",
    "expression. Gender cannot be inferred from a transcriptome and is not\n",
    "reported. Say 'inferred genetic sex', never 'sex and gender'.\n", sep = "")
