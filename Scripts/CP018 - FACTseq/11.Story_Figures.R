################################################################################
# CP018 story figures — designed around where the data is actually strong
#
# WHY THIS REPLACES THE ATLAS-MAPPING FRAMING
#   The earlier panels leaned on assigning CP018's antigen-specific cells to an
#   external PBMC atlas and reading the cluster NAME ("Tex CD8") as a phenotype.
#   Three things are wrong with that:
#     1. n = 7 cells. Percentages on 7 cells are not a finding.
#     2. The atlas cluster does not meet our own exhaustion criterion
#        (exhaustion:stemness 0.71, threshold 1.5) - the NAME is not the biology.
#     3. It throws away the well-powered data. Antigen-specific cells available:
#          CP003  stim 1515  vs  unstim CD8-restricted 121   -> 12.5x more
#          CP018  stim  311  vs  unstim CD8-restricted   7   -> 44.4x more
#        (ratios computed from the objects; see Story_power_ratio.csv, which the
#         script writes so this comment cannot drift from the data again.)
#
#   So: the STIMULATED data carries every phenotype claim; the unstimulated data
#   is used only for what it is uniquely good at - anchoring the TIMELINE.
#
# THE DESIGN
#   Each participant is compared to ITS OWN non-specific cells, and the two
#   within-participant effect sizes are then compared. This is the only contrast
#   that survives n=1-per-group: a global batch or aligner shift moves a
#   participant's specific and non-specific cells together and cancels out.
#   (It matters here - CP018's background IFN score is genuinely higher, so
#   ABSOLUTE levels between participants are not interpretable.)
#
# THE RESULT
#   Of five programmes, four are CONSERVED (same direction, similar magnitude):
#   cytotoxicity, effector, exhaustion, stemness. Exactly one INVERTS:
#   Type I IFN memory (Cohen's d +0.50 in CP003 vs -0.53 in CP018).
#   That is the paper's Figure 4/5 axis, measured in antigen-validated cells.
#
# Output: STORY root - panels + source CSV each
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(RColorBrewer); library(scales)
})

STORY <- file.path(CP018_ANALYSIS, "Story_Figures"); mkdirs(STORY)
csv_s <- file.path(STORY, "tables");                mkdirs(csv_s)

# ── Readable typography ──────────────────────────────────────────────────────
# The manuscript panels use base_size 36 because they are rendered at 16x12in.
# Reusing that size on a 7in panel gives the unreadable "giant text, tiny data"
# result. Type is sized for the canvas each panel is actually saved at.
theme_story <- function(base = 13) {
  theme_classic(base_size = base) +
    theme(plot.title       = element_text(size = base + 4, face = "bold"),
          plot.subtitle    = element_text(size = base, colour = "grey30"),
          plot.caption     = element_text(size = base - 3, colour = "grey45",
                                          hjust = 0, lineheight = 1.15),
          axis.title       = element_text(size = base + 1, face = "bold"),
          axis.text        = element_text(size = base),
          legend.title     = element_text(size = base, face = "bold"),
          legend.text      = element_text(size = base - 1),
          strip.text       = element_text(size = base + 1, face = "bold"),
          strip.background = element_blank(),
          plot.title.position = "plot")
}
# Captions must be hard-wrapped: ggplot does not wrap them, so a long single
# line runs off the canvas and is silently clipped.
wrap_cap <- function(x, width = 95) paste(strwrap(x, width = width), collapse = "\n")

P_COLS <- c("CP003 (never suppressed)" = "#B2182B",
            "CP018 (suppressed)"       = "#2166AC")
P_SHORT <- c(CP003 = "#B2182B", CP018 = "#2166AC")

################################################################################
# 1. Score both participants' stimulated data on the same programmes
################################################################################
MODULES <- list(
  `Type I IFN memory` = c("ISG15","IFI6","MX1","OAS1","STAT1","IRF7","IFI44L"),
  Stemness            = c("TCF7","LEF1","BACH2","ID3","SELL","CCR7","IL7R","FOXP1","MYB"),
  Exhaustion          = c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1","CD160"),
  Cytotoxicity        = c("GZMB","GZMA","GZMH","PRF1","GNLY","NKG7","KLRD1"),
  Effector            = c("TBX21","KLRG1","CX3CR1","ZEB2","PRDM1"))

collect_scores <- function(obj, label) {
  DefaultAssay(obj) <- "RNA"
  if (length(Layers(obj[["RNA"]])) > 3) obj <- JoinLayers(obj, assay = "RNA")
  X <- GetAssayData(obj, layer = "data")
  out <- data.frame(
    participant = label,
    hiv = obj$HIV_Specific_TCR == "HIV-Specific TCR",
    stringsAsFactors = FALSE)
  for (m in names(MODULES)) {
    g <- intersect(MODULES[[m]], rownames(X))
    if (!length(g)) { out[[m]] <- NA_real_; next }
    out[[m]] <- as.numeric(Matrix::colMeans(X[g, , drop = FALSE]))
  }
  attr(out, "gene_coverage") <- vapply(MODULES, function(v)
    sum(v %in% rownames(X)) / length(v), numeric(1))
  out
}

cp03 <- qs_read(QS_CP003_ANNOT_FIG3)
s03  <- collect_scores(cp03, "CP003"); rm(cp03); gc(verbose = FALSE)
cp18 <- qs_read(QS_ANNOTATED)
s18  <- collect_scores(cp18, "CP018")

scores <- bind_rows(s03, s18)
scores$participant <- factor(scores$participant, levels = c("CP003", "CP018"))
save_csv(scores %>% group_by(participant, hiv) %>%
           summarise(across(all_of(names(MODULES)), ~ round(mean(.x), 4)),
                     n = n(), .groups = "drop"),
         file.path(csv_s, "Story_module_means.csv"))

# Record the power ratio that motivates this design, computed rather than
# asserted - an earlier version of this header carried a typed range that did
# not match the data.
un18_pow <- qs_read(file.path(SAVED_DIR, "tara_cp018_crossref.qs2"))
un03_pow <- qs_read(file.path(SAVED_DIR, "tara_cp003_fig3.qs2"))
CD8_CL <- c("Naive 1 CD8","Naive 2 CD8","Naive Intermediate CD8",
            "Tscm CD8","TEMRA/CTL CD8","Tex CD8")
pow <- data.frame(
  participant = c("CP003", "CP018"),
  stim_antigen_specific = c(sum(scores$participant == "CP003" & scores$hiv),
                            sum(scores$participant == "CP018" & scores$hiv)),
  unstim_matched_CD8 = c(
    sum(un03_pow$HIV_validated_status == "HIV-Specific (validated)" &
        as.character(un03_pow$Annotation) %in% CD8_CL),
    sum(un18_pow$clone_match == "HIV-Specific Match" &
        as.character(un18_pow$Annotation) %in% CD8_CL)),
  stringsAsFactors = FALSE)
pow$fold_more_in_stim <- round(pow$stim_antigen_specific / pow$unstim_matched_CD8, 1)
save_csv(pow, file.path(csv_s, "Story_power_ratio.csv"))
cat("\n=== Why the stimulated data carries the phenotype claims ===\n")
print(pow, row.names = FALSE)
rm(un03_pow); gc(verbose = FALSE)

################################################################################
# 2. PANEL A — the headline: within-participant effect size, side by side
################################################################################
# Cohen's d of antigen-specific vs that participant's OWN non-specific cells.
# Plotted as a paired slope so a sign flip is the visual event.
cohen_d <- function(x, y) {
  sp <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
             (length(x) + length(y) - 2))
  (mean(x) - mean(y)) / sp
}

eff <- do.call(rbind, lapply(names(MODULES), function(m) {
  ds <- vapply(c("CP003", "CP018"), function(p) {
    s <- scores[scores$participant == p, ]
    cohen_d(s[[m]][s$hiv], s[[m]][!s$hiv])
  }, numeric(1))
  fit <- lm(as.formula(paste0("`", m, "` ~ hiv * participant")), data = scores)
  cf  <- summary(fit)$coefficients
  data.frame(module = m, CP003 = ds[["CP003"]], CP018 = ds[["CP018"]],
             delta_d = ds[["CP018"]] - ds[["CP003"]],
             sign_flip = sign(ds[["CP003"]]) != sign(ds[["CP018"]]),
             interaction_p = cf["hivTRUE:participantCP018", "Pr(>|t|)"],
             stringsAsFactors = FALSE)
}))
eff$interaction_p_adj <- p.adjust(eff$interaction_p, "BH")
eff$n_CP003 <- sum(scores$participant == "CP003" & scores$hiv)
eff$n_CP018 <- sum(scores$participant == "CP018" & scores$hiv)
save_csv(eff, file.path(csv_s, "PanelA_within_participant_effect_sizes.csv"))
cat("\n=== Within-participant effect sizes (Cohen's d) ===\n")
print(eff %>% mutate(across(c(CP003, CP018, delta_d), ~ round(.x, 3)),
                     interaction_p_adj = signif(interaction_p_adj, 3)) %>%
        arrange(desc(abs(delta_d))) %>%
        select(module, CP003, CP018, delta_d, sign_flip, interaction_p_adj),
      row.names = FALSE)

eff_long <- eff %>%
  select(module, CP003, CP018, sign_flip, interaction_p_adj) %>%
  pivot_longer(c(CP003, CP018), names_to = "participant", values_to = "d") %>%
  mutate(module = factor(module, levels = eff$module[order(eff$delta_d)]))

star <- function(p) ifelse(p < 1e-10, "****", ifelse(p < 1e-4, "***",
                    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
lab_a <- eff %>% mutate(module = factor(module, levels = levels(eff_long$module)),
                        lab = paste0(star(interaction_p_adj)),
                        xpos = pmax(CP003, CP018))

panelA <- ggplot(eff_long, aes(x = d, y = module)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5) +
  geom_line(aes(group = module,
                colour = ifelse(sign_flip, "flip", "conserved")),
            linewidth = 1.5, alpha = 0.75) +
  geom_point(aes(fill = participant), shape = 21, size = 5.4,
             colour = "white", stroke = 1.2) +
  geom_text(data = lab_a, aes(x = xpos, label = lab), hjust = -0.45,
            size = 4.2, fontface = "bold", colour = "grey25") +
  scale_colour_manual(values = c(flip = "#D95F02", conserved = "grey78"),
                      guide = "none") +
  scale_fill_manual(values = P_SHORT, name = NULL,
                    labels = c("CP003 (never suppressed)", "CP018 (suppressed)")) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.16))) +
  theme_story(13) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey93")) +
  labs(x = "Cohen's d  (antigen-specific vs that participant's own other cells)",
       y = NULL,
       title = wrap_cap(paste0("Type I IFN memory is reduced in the SUPPRESSED ",
                               "participant's antigen-specific cells; four other ",
                               "programmes are unchanged"), 78),
       subtitle = "Antigen-specific vs non-specific cells, measured within each participant",
       caption = wrap_cap(paste0("Orange link = opposite direction. Stars: ",
                        "BH-adjusted interaction. n = ", eff$n_CP003[1],
                        " (CP003) and ", eff$n_CP018[1],
                        " (CP018) antigen-specific cells."), 105))
save_png(file.path(STORY, "A_effect_size_inversion.png"), panelA,
         width = 10.4, height = 6.1, dpi = 400)

################################################################################
# 3. PANEL B — the inverting programme, shown as distributions
################################################################################
ifn <- scores %>%
  transmute(participant, hiv,
            score = .data[["Type I IFN memory"]],
            grp = ifelse(hiv, "Antigen-specific", "Other cells"),
            pl  = ifelse(participant == "CP003", "CP003 (never suppressed)",
                                                 "CP018 (suppressed)"))
save_csv(ifn %>% group_by(pl, grp) %>%
           summarise(n = n(), mean = round(mean(score), 4),
                     median = round(median(score), 4), .groups = "drop"),
         file.path(csv_s, "PanelB_IFN_memory_distributions.csv"))

# Within-participant test; between-participant absolute levels are NOT compared
# because the backgrounds differ (see the caption).
bt <- ifn %>% group_by(pl) %>%
  summarise(p = wilcox.test(score ~ grp)$p.value,
            hi = max(score), .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, "BH"), lab = paste0("p = ", signif(p_adj, 2)))

panelB <- ggplot(ifn, aes(x = grp, y = score, fill = grp)) +
  geom_violin(scale = "width", width = 0.72, alpha = 0.62,
              colour = "grey30", linewidth = 0.4) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white",
               colour = "grey20", linewidth = 0.45) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8,
               fill = "black", colour = "white") +
  geom_text(data = bt, aes(x = 1.5, y = hi * 1.02, label = lab),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~ pl, nrow = 1) +
  scale_fill_manual(values = c(`Antigen-specific` = "#E6A23C",
                               `Other cells` = "grey80"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  theme_story(13) +
  labs(x = NULL, y = "Type I IFN memory score",
       title = wrap_cap(paste0("Interferon-stimulated genes are raised in the ",
                               "UNSUPPRESSED participant's antigen-specific cells ",
                               "and lowered in the SUPPRESSED participant's"), 82),
       subtitle = "Type I IFN memory score per cell",
       caption = wrap_cap(paste0("Diamond = mean, box = median and quartiles. ",
                        "Wilcoxon within each participant."), 100))
save_png(file.path(STORY, "B_IFN_memory_inversion.png"), panelB,
         width = 9.2, height = 5, dpi = 400)

################################################################################
# 4. PANEL C — lineage traces, colour NOT carried across a detection gap
################################################################################
# A continuous coloured line across a timepoint where the lineage was not
# detected asserts something the data does not support. Coloured segments are
# drawn ONLY between consecutive DETECTED timepoints; a gap is a faint grey
# dashed connector, which reads as "same lineage, not observed in between".
CP018_TP <- c("1 mo\n(pre-ART)", "2 mo\n(HIV-stim)", "25 mo", "42 mo",
              "90 mo\n(HIV-stim)")
CP003_TP <- c("2 mo\n(HIV-stim)", "1 mo", "12 mo", "24 mo", "101 mo\n(HIV-stim)")

trace_panel <- function(dat, tp_levels, title, subtitle) {
  if (!nrow(dat)) return(ggplot() + theme_void())
  ord <- dat %>% group_by(clone) %>% summarise(tot = sum(n_cells), .groups = "drop") %>%
    arrange(desc(tot))
  d <- dat %>% mutate(tp = factor(as.character(tp), levels = tp_levels),
                      x  = as.integer(tp),
                      lab = factor(substr(clone, 1, 24),
                                   levels = substr(ord$clone, 1, 24))) %>%
    filter(!is.na(x)) %>% arrange(lab, x)

  # Consecutive detected pairs: coloured where adjacent, grey dashed where the
  # lineage skips one or more timepoints.
  segs <- d %>% group_by(lab) %>%
    mutate(x2 = lead(x), n2 = lead(n_cells)) %>%
    filter(!is.na(x2)) %>%
    mutate(gap = (x2 - x) > 1) %>% ungroup()

  pal <- colorRampPalette(brewer.pal(max(3, min(nrow(ord), 8)), "Dark2"))(nrow(ord))
  names(pal) <- levels(d$lab)

  ggplot(d, aes(x = x, y = lab)) +
    # dashed bridge FIRST so coloured segments draw over it
    geom_segment(data = segs %>% filter(gap),
                 aes(x = x, xend = x2, y = lab, yend = lab),
                 inherit.aes = FALSE, linetype = "22",
                 linewidth = 0.7, colour = "grey72") +
    geom_segment(data = segs %>% filter(!gap),
                 aes(x = x, xend = x2, y = lab, yend = lab, colour = lab),
                 inherit.aes = FALSE, linewidth = 2.1, alpha = 0.85) +
    # sqrt area scaling with a FLOOR: scale_size_area maps a 1-cell point to a
    # radius thinner than the connecting line, so the count becomes unreadable
    # exactly where it matters most (the single-cell late detections).
    geom_point(aes(size = n_cells, colour = lab)) +
    geom_text(aes(label = n_cells), colour = "white",
              size = 3.1, fontface = "bold") +
    scale_size(range = c(4.6, 11), trans = "sqrt", guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = seq_along(tp_levels), labels = tp_levels,
                       limits = c(0.6, length(tp_levels) + 0.4),
                       expand = c(0, 0)) +
    theme_story(12) +
    theme(axis.text.y = element_text(family = "mono", size = 9),
          axis.text.x = element_text(size = 10, lineheight = 0.92),
          panel.grid.major.y = element_line(colour = "grey94")) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle)
}

# Build the lineage tables at CLONOTYPE level (specificity is a receptor
# property; see script 10 for why the per-cell filter truncates lineages).
lineage_tbl <- function(stim_md, ct, tpmap, un_md, un_ct, un_tpmap,
                        stim_hiv, un_hiv) {
  s <- stim_md; s$.ct <- s[[ct]]; s <- s[!is.na(s$.ct) & s$.ct != "NA", ]
  u <- un_md;   u$.ct <- u[[un_ct]]; u <- u[!is.na(u$.ct) & u$.ct != "NA", ]
  keep <- union(unique(s$.ct[stim_hiv(s)]), unique(u$.ct[un_hiv(u)]))
  bind_rows(
    s %>% filter(.ct %in% keep) %>% mutate(tp = tpmap(.)) %>%
      count(clone = .ct, tp, name = "n_cells"),
    u %>% filter(.ct %in% keep) %>% mutate(tp = un_tpmap(.)) %>%
      count(clone = .ct, tp, name = "n_cells")) %>%
    group_by(clone) %>% filter(n_distinct(tp) >= 2) %>% ungroup()
}

un18 <- qs_read(file.path(SAVED_DIR, "tara_cp018_crossref.qs2"))
L18 <- lineage_tbl(
  cp18@meta.data, CT_CANON,
  function(d) ifelse(d$Timepoint == "2m", "2 mo\n(HIV-stim)", "90 mo\n(HIV-stim)"),
  un18@meta.data, CT_CANON,
  function(d) unname(c(CP018_1m = "1 mo\n(pre-ART)", CP018_25m = "25 mo",
                       CP018_42m = "42 mo")[as.character(d$orig.ident)]),
  function(d) d$HIV_Specific_TCR == "HIV-Specific TCR",
  function(d) d$clone_match == "HIV-Specific Match")
save_csv(L18, file.path(csv_s, "PanelC_CP018_lineages.csv"))

TRACE_SUB <- wrap_cap(paste0("Coloured link = consecutive detections | ",
  "grey dashed = same lineage, not detected in between"), 88)
trace18 <- trace_panel(L18, CP018_TP, "CP018 (suppressed)", TRACE_SUB)
save_png(file.path(STORY, "C_CP018_lineage_traces.png"), trace18,
         width = 9.5, height = 3.6, dpi = 400)

cp03b <- qs_read(QS_CP003_ANNOT_FIG3); cp03b <- add_ct_canon(cp03b, verbose = FALSE)
un03  <- qs_read(file.path(SAVED_DIR, "tara_cp003_fig3.qs2"))
un03  <- add_ct_canon(un03, verbose = FALSE)
L03 <- lineage_tbl(
  cp03b@meta.data, CT_CANON,
  function(d) ifelse(d$Months == 2, "2 mo\n(HIV-stim)", "101 mo\n(HIV-stim)"),
  un03@meta.data, CT_CANON,
  function(d) paste0(d$Age, " mo"),
  function(d) d$HIV_Specific_TCR == "HIV-Specific TCR",
  function(d) d$HIV_validated_status == "HIV-Specific (validated)")
save_csv(L03, file.path(csv_s, "PanelC_CP003_lineages.csv"))

trace03 <- trace_panel(L03, CP003_TP, "CP003 (never suppressed)", TRACE_SUB)

save_png(file.path(STORY, "C_lineage_traces_both.png"),
         (trace03 / trace18) +
           plot_layout(heights = c(nrow(distinct(L03, clone)),
                                   nrow(distinct(L18, clone))) + 1.4) +
           plot_annotation(
             title = wrap_cap(paste0("HIV-specific clonotypes are detected from ",
                               "before ART initiation through 90 months of ",
                               "suppression"), 88),
             subtitle = paste0("Point area = cells recovered; a gap means the ",
                               "clonotype was not sampled at that draw"),
             theme = theme(plot.title = element_text(size = 16, face = "bold"),
                           plot.subtitle = element_text(size = 11.5, colour = "grey30"))),
         width = 10, height = 8.4, dpi = 400)

################################################################################
# 5. PANEL D — persistence rate, the one thing the timeline data supports
################################################################################
rate <- data.frame(
  pl = c("CP003 (never suppressed)", "CP018 (suppressed)"),
  hiv_clonotypes = c(
    n_distinct(cp03b@meta.data[[CT_CANON]][
      cp03b$HIV_Specific_TCR == "HIV-Specific TCR" &
      !is.na(cp03b@meta.data[[CT_CANON]])]),
    n_distinct(cp18@meta.data[[CT_CANON]][
      cp18$HIV_Specific_TCR == "HIV-Specific TCR" &
      !is.na(cp18@meta.data[[CT_CANON]])])),
  persisting = c(n_distinct(L03$clone), n_distinct(L18$clone)),
  stringsAsFactors = FALSE) %>%
  mutate(pct = round(100 * persisting / hiv_clonotypes, 1))
save_csv(rate, file.path(csv_s, "PanelD_persistence_rate.csv"))
cat("\n=== Persistence rate ===\n"); print(rate, row.names = FALSE)

ft <- fisher.test(matrix(c(rate$persisting,
                           rate$hiv_clonotypes - rate$persisting), nrow = 2))
panelD <- ggplot(rate, aes(x = pl, y = pct, fill = pl)) +
  geom_col(width = 0.5, colour = "grey25", linewidth = 0.4) +
  geom_text(aes(label = paste0(pct, "%\n", persisting, "/", hiv_clonotypes)),
            vjust = -0.3, size = 4.1, fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = P_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
  theme_story(13) +
  theme(axis.text.x = element_text(size = 11.5)) +
  labs(x = NULL, y = "% of antigen-specific clonotypes\ndetected at \u22652 timepoints",
       title = wrap_cap(paste0("HIV-specific clonotypes persist equally well with ",
                               "and without viral suppression"), 52),
       subtitle = paste0("Fisher exact p = ", signif(ft$p.value, 2)),
       caption = wrap_cap(paste0("Clonotypes detected at two or more of that ",
                        "participant's timepoints."), 58))
save_png(file.path(STORY, "D_persistence_rate.png"), panelD,
         width = 6.4, height = 5, dpi = 400)

################################################################################
# 6. Assembled story figure
################################################################################
save_png(file.path(STORY, "STORY_assembled.png"),
         ((panelA | panelB) / (panelD | trace18)) +
           plot_layout(heights = c(1, 0.92)) +
           plot_annotation(
             tag_levels = "A",
             title = wrap_cap(paste0("Retention of HIV-specific CD8 clonotypes is ",
                               "unchanged under suppression; their interferon ",
                               "programme is not"), 92),
             subtitle = wrap_cap(paste0("A: effect sizes per programme. B: the ",
                               "interferon score itself. C: persistence rate. ",
                               "D: individual clonotypes over time."), 110),
             # patchwork draws its annotation outside the child plots' margins,
             # so padding must be set on the composite theme itself.
             theme = theme(plot.title = element_text(size = 19, face = "bold"),
                           plot.subtitle = element_text(size = 12.5, colour = "grey25",
                                                        lineheight = 1.2),
                           plot.margin = margin(14, 80, 14, 14))) &
           theme(plot.tag = element_text(size = 17, face = "bold")),
         width = 19, height = 11.5, dpi = 350)

cat("\nOutputs in:\n  ", STORY, "\n")
cat("\nFRAMING NOTE: no claim here rests on assigning cells to an external\n",
    "atlas or on reading a cluster NAME as a phenotype. Every phenotype\n",
    "statement is measured in the stimulated data, within participant.\n", sep = "")
