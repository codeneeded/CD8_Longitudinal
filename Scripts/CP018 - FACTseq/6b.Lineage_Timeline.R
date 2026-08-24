################################################################################
# CP018 — full lineage timeline: pre-ART entry through 90 months
#
# Joins the UNSTIMULATED eCITE-seq timepoints (CP018_1m = pre-ART entry,
# CP018_25m and CP018_42m = suppressed) to the two FACT-Seq stimulation
# timepoints (2m, 90m), for the antigen-validated clonotypes recovered by
# script 6.
#
# MATCHING: normalised CTstrict (see normalise_ct() in config.R). Raw CTstrict
# cannot match across these experiments because Cell Ranger 10.1.0 writes an
# empty D-gene field where the older alignments write ".NA.".
#
# Input : tara_cp018_crossref.qs2 (script 6), QS_ANNOTATED (script 5)
# Output: TL_ROOT — timeline table, alluvial, resting-phenotype test
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(ggalluvial); library(scales)
})

TL_ROOT <- file.path(HIV_ROOT, "Lineage_Timeline"); mkdirs(TL_ROOT)
csv_tl  <- file.path(TL_ROOT, "tables");            mkdirs(csv_tl)

CROSSREF <- file.path(SAVED_DIR, "tara_cp018_crossref.qs2")
require_file(CROSSREF,     "script 6 cross-reference object")
require_file(QS_ANNOTATED, "CP018 annotated stim object")

# Timepoint ordering across BOTH experiments, in real time since ART start.
TP_LEVELS <- c("1m (pre-ART)", "2m (stim)", "25m", "42m", "90m (stim)")
TP_MAP <- c(CP018_1m = "1m (pre-ART)", CP018_25m = "25m", CP018_42m = "42m",
            `2m` = "2m (stim)", `90m` = "90m (stim)")

################################################################################
# 1. Assemble per-clonotype, per-timepoint detections
################################################################################
o <- qs_read(CROSSREF); u <- o@meta.data
require_ct_canon(o, "unstimulated crossref object")
u$ctn <- u[[CT_CANON]]
CT_COL <- "predicted.celltype.l1"
matched_ct <- unique(u$ctn[u$clone_match == "HIV-Specific Match" & !is.na(u$ctn)])

stim <- qs_read(QS_ANNOTATED); sm <- stim@meta.data
require_ct_canon(stim, "CP018 stim object")
sm$ctn <- sm[[CT_CANON]]

cat("Antigen-validated clonotypes recovered in unstimulated data:",
    length(matched_ct), "\n")

rows <- list()
for (cc in matched_ct) {
  su <- u[!is.na(u$ctn) & u$ctn == cc, ]
  ss <- sm[!is.na(sm$ctn) & sm$ctn == cc, ]
  cdr3 <- if (nrow(ss)) ss$CTaa[1] else su$CTaa[1]

  for (tp in unique(su$orig.ident)) {
    k <- su$orig.ident == tp
    rows[[length(rows) + 1]] <- data.frame(
      clone_cdr3 = cdr3, ctn = cc,
      experiment = "Unstimulated (eCITE-seq)",
      timepoint  = unname(TP_MAP[as.character(tp)]),
      ART_status = as.character(su$Timepoint_Group[k][1]),
      cells      = sum(k),
      state      = paste(sort(unique(as.character(su$Annotation[k]))), collapse = "; "),
      stringsAsFactors = FALSE)
  }
  for (tp in unique(ss$Timepoint)) {
    k <- ss$Timepoint == tp
    rows[[length(rows) + 1]] <- data.frame(
      clone_cdr3 = cdr3, ctn = cc,
      experiment = "Stimulated (FACT-Seq)",
      timepoint  = unname(TP_MAP[as.character(tp)]),
      ART_status = "PostART_Suppressed",
      cells      = sum(k),
      state      = paste(sort(unique(as.character(ss$Fig_Annotation[k]))), collapse = "; "),
      stringsAsFactors = FALSE)
  }
}
TL <- do.call(rbind, rows) %>%
  mutate(timepoint = factor(timepoint, levels = TP_LEVELS)) %>%
  arrange(clone_cdr3, timepoint)

save_csv(TL, file.path(csv_tl, "CP018_lineage_timeline.csv"))
cat("\n=== Lineage timeline ===\n")
print(as.data.frame(TL %>% mutate(
        clone = ifelse(nchar(clone_cdr3) > 30,
                       paste0(substr(clone_cdr3, 1, 29), "..."), clone_cdr3)) %>%
        select(clone, experiment, timepoint, ART_status, cells, state)),
      row.names = FALSE)

# min()/max() are undefined on an unordered factor - use the level index.
span <- TL %>% mutate(.i = as.integer(timepoint)) %>% group_by(clone_cdr3) %>%
  summarise(n_timepoints = n_distinct(timepoint),
            first = TP_LEVELS[min(.i)], last = TP_LEVELS[max(.i)],
            experiments = n_distinct(experiment), total_cells = sum(cells),
            .groups = "drop") %>% arrange(desc(n_timepoints), desc(total_cells))
save_csv(span, file.path(csv_tl, "CP018_lineage_span.csv"))
cat("\n=== Span per clonotype ===\n"); print(as.data.frame(span), row.names = FALSE)

################################################################################
# 2. Timeline figure — cells detected per timepoint, per clonotype
################################################################################
# Every clonotype gets a row at every timepoint, so absence is drawn as zero
# rather than as a gap the eye reads as missing data.
grid <- expand.grid(clone_cdr3 = unique(TL$clone_cdr3), timepoint = TP_LEVELS,
                    stringsAsFactors = FALSE) %>%
  left_join(TL %>% select(clone_cdr3, timepoint, cells, experiment, state),
            by = c("clone_cdr3", "timepoint")) %>%
  mutate(cells = ifelse(is.na(cells), 0, cells),
         timepoint = factor(timepoint, levels = TP_LEVELS),
         short = ifelse(nchar(clone_cdr3) > 30,
                        paste0(substr(clone_cdr3, 1, 29), "\u2026"), clone_cdr3))
save_csv(grid, file.path(csv_tl, "CP018_timeline_source_data.csv"))

p_tl <- ggplot(grid, aes(x = as.integer(timepoint), y = cells,
                         group = short, colour = short)) +
  # Shade the pre-ART timepoint. On a DISCRETE x-axis, annotate() must use the
  # level's integer position, not a continuous coordinate.
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = -Inf, ymax = Inf,
           fill = "grey88", alpha = 0.5) +
  annotate("label", x = 1, y = Inf, label = "pre-ART", vjust = 1.1,
           size = 3.3, colour = "grey25", label.size = 0, fill = "white") +
  geom_line(linewidth = 1.05, alpha = 0.9) +
  geom_point(aes(shape = experiment), size = 3.4) +
  geom_text(data = subset(grid, cells > 0), aes(label = cells),
            vjust = -0.9, size = 3.5, show.legend = FALSE) +
  scale_shape_manual(values = c("Unstimulated (eCITE-seq)" = 16,
                                "Stimulated (FACT-Seq)" = 17),
                     na.translate = FALSE, name = NULL) +
  scale_colour_brewer(palette = "Dark2", name = "Clonotype (CDR3\u03b1_CDR3\u03b2)") +
  scale_x_continuous(breaks = seq_along(TP_LEVELS), labels = TP_LEVELS,
                     limits = c(0.5, length(TP_LEVELS) + 0.5)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Cells detected",
       title = "Antigen-validated CD8 lineages, pre-ART entry to 90 months",
       subtitle = paste0(nrow(span), " clonotypes tracked across two experiments. ",
                         "Zero = not detected, NOT proven absent\n",
                         "(see detection-power controls).")) +
  theme_cp018(base_size = 13) +
  theme(plot.title.position = "plot",
        plot.title    = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        legend.position = "right", legend.box = "vertical")
save_png(file.path(TL_ROOT, "CP018_lineage_timeline.png"), p_tl,
         width = 11, height = 6.5)

################################################################################
# 3. Alluvial: cell state at each detection
################################################################################
allu <- TL %>%
  mutate(state_simple = ifelse(state == "" | is.na(state), "not annotated", state),
         short = ifelse(nchar(clone_cdr3) > 24,
                        paste0(substr(clone_cdr3, 1, 23), "\u2026"), clone_cdr3))
if (n_distinct(allu$timepoint) > 1) {
  p_al <- ggplot(allu, aes(x = timepoint, y = cells, stratum = state_simple,
                           alluvium = short, fill = short)) +
    geom_flow(alpha = 0.65, width = 0.32) +
    geom_stratum(width = 0.32, colour = "black", linewidth = 0.3, fill = "grey96") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.1) +
    scale_fill_brewer(palette = "Dark2", name = "Clonotype") +
    labs(x = NULL, y = "Cells",
         title = "Cell state of antigen-validated lineages over time") +
    theme_cp018(base_size = 13) + theme(plot.title.position = "plot")
  save_png(file.path(TL_ROOT, "CP018_lineage_state_alluvial.png"), p_al,
           width = 11, height = 7)
}

################################################################################
# 4. Resting phenotype: are these clones Tscm at rest?
################################################################################
# The stimulated data could not separate genuine Tscm from stimulation-lowered
# IL7R. The unstimulated data can: compare matched CD8 cells against all other
# CD8 cells (NOT all lineages - pooling CD4/B/NK/monocytes would inflate every
# T-cell signal).
cd8    <- !is.na(u[[CT_COL]]) & u[[CT_COL]] == "CD8 T"
is_mat <- u$clone_match == "HIV-Specific Match"
cat("\nCD8 T cells in unstimulated CP018:", sum(cd8),
    "| antigen-validated among them:", sum(cd8 & is_mat), "\n")

DefaultAssay(o) <- "RNA"
PANELS <- list(
  `Stemness / Tscm` = c("IL7R","SELL","CCR7","TCF7","LEF1","BACH2","CD27","CD28","FAS"),
  `Exhaustion`      = c("PDCD1","TOX","LAG3","TIGIT","HAVCR2","CTLA4","ENTPD1"),
  `Effector`        = c("GZMB","GNLY","NKG7","KLRG1","PRF1","GZMH"))
X <- GetAssayData(o, layer = "data")

rest <- lapply(names(PANELS), function(pn) {
  g <- intersect(PANELS[[pn]], rownames(X))
  if (!length(g)) return(NULL)
  data.frame(panel = pn, gene = g,
             matched_CD8 = round(Matrix::rowMeans(X[g, cd8 &  is_mat, drop = FALSE]), 3),
             other_CD8   = round(Matrix::rowMeans(X[g, cd8 & !is_mat, drop = FALSE]), 3),
             row.names = NULL, stringsAsFactors = FALSE)
})
REST <- do.call(rbind, rest)
REST$fold_vs_other <- round(REST$matched_CD8 / (REST$other_CD8 + 0.01), 2)
save_csv(REST, file.path(csv_tl, "CP018_resting_phenotype_matched_CD8.csv"))
cat("\n=== Resting phenotype: antigen-validated vs other CD8 (unstimulated) ===\n")
print(REST, row.names = FALSE)

# Is there a Tscm cluster in the unstimulated annotation, and are any matched
# cells in it? A named Tscm population with zero matched cells is a direct
# negative answer, not an absence of evidence.
if ("Annotation" %in% names(u)) {
  tscm_lv <- grep("Tscm", unique(as.character(u$Annotation)), value = TRUE)
  if (length(tscm_lv)) {
    n_tscm     <- sum(u$Annotation %in% tscm_lv, na.rm = TRUE)
    n_tscm_mat <- sum(u$Annotation %in% tscm_lv & is_mat, na.rm = TRUE)
    cat("\nTscm cluster(s) in unstimulated data:", paste(tscm_lv, collapse = ", "),
        "\n  cells:", n_tscm, "| antigen-validated among them:", n_tscm_mat, "\n")
    save_csv(data.frame(tscm_cluster = paste(tscm_lv, collapse = "; "),
                        n_cells = n_tscm, n_antigen_validated = n_tscm_mat),
             file.path(csv_tl, "CP018_Tscm_cluster_occupancy.csv"))
  }
}

state_tab <- as.data.frame(table(Timepoint = u$orig.ident[is_mat],
                                 State     = as.character(u$Annotation[is_mat])))
state_tab <- state_tab[state_tab$Freq > 0, ]
save_csv(state_tab, file.path(csv_tl, "CP018_matched_state_by_timepoint.csv"))
cat("\n=== Matched cells: state by timepoint ===\n")
print(state_tab, row.names = FALSE)

cat("\nOutputs in:\n  ", TL_ROOT, "\n")
cat("\nINTERPRETATION GUARD: n = 3 clonotypes / <10 cells. Directions are\n",
    "descriptive; a zero at any timepoint means NOT DETECTED, not absent.\n", sep = "")
