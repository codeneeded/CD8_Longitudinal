################################################################################
# Figure 3 mirror for CP018 — panel-for-panel, plus CP003+CP018 combined
#
# WHY THIS SCRIPT EXISTS
#   The earlier CP018 plots were built from first principles and did not match
#   the manuscript's Figure 3. They are therefore not directly comparable to the
#   published panels. This script reproduces Figure 3's actual structure,
#   palette and styling for CP018, reading the construction from
#   "Scripts/3. Manuscript/Figure_3.R" rather than reinventing it.
#
# THE KEY STRUCTURAL POINT I HAD MISSED
#   Figure 3C does NOT plot the stimulation timepoints alone. It interleaves the
#   STIMULATED and UNSTIMULATED timepoints on a single real-time axis:
#
#     CP003:  2 months (HIV-stim) | 1 months | 12 months | 24 months | 101 months (HIV-stim)
#     CP018:  1 months (pre-ART)  | 2 months (HIV-stim) | 25 months | 42 months | 90 months (HIV-stim)
#
#   That is why the published alluvial shows persistence where mine showed
#   almost none: a clonotype absent from the other stim timepoint can still be
#   caught at an unstimulated one. CP018 additionally has a PRE-ART timepoint,
#   which CP003 does not - so its axis starts before treatment.
#
# PANELS (mirroring Fig3B1-F)
#   B1  UMAP, annotated clusters              D  functional gene profile
#   B2  UMAP, HIV-specific overlay            E  persistence histogram
#   C   alluvial, interleaved timepoints      F  unstim cluster composition
#   Combined: C, E and F with CP003 and CP018 side by side.
#
# Input : QS_ANNOTATED, tara_cp018_crossref.qs2, QS_CP003_ANNOT_FIG3,
#         tara_cp003_fig3.qs2
# Output: FIG3M root - panels, combined panels, source-data CSV per panel
################################################################################

source("~/Desktop/Projects/Repositories/CD8_Longitudinal/Scripts/CP018 - FACTseq/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(ggalluvial); library(patchwork)
  library(RColorBrewer); library(scales)
})

FIG3M   <- file.path(CP018_ANALYSIS, "Figure3_Mirror"); mkdirs(FIG3M)
csv_f3  <- file.path(FIG3M, "tables");                  mkdirs(csv_f3)
CMB     <- file.path(FIG3M, "Combined");                mkdirs(CMB)

CROSSREF   <- file.path(SAVED_DIR, "tara_cp018_crossref.qs2")
TARA_CP003 <- file.path(SAVED_DIR, "tara_cp003_fig3.qs2")
require_file(QS_ANNOTATED, "CP018 annotated stim object (script 5)")
require_file(CROSSREF,     "CP018 unstim crossref (script 6)")

# ── Manuscript styling, copied from Figure_3.R so panels are interchangeable ──
# base_size 36 with 32-48pt text: the manuscript panels are rendered large and
# downscaled, so type stays legible in a multi-panel composite.
FIG3_BASE <- 36
HIV_RED   <- "#E63946"

theme_fig3 <- function(base_size = FIG3_BASE) {
  theme_classic(base_size = base_size) +
    theme(axis.text      = element_text(size = base_size - 4, face = "bold"),
          axis.title     = element_text(size = base_size, face = "bold"),
          plot.title     = element_text(size = base_size + 8, face = "bold"),
          plot.subtitle  = element_text(size = base_size - 4, colour = HIV_RED,
                                        face = "bold"),
          plot.title.position = "plot")
}

hiv_cols <- c("Other" = "#E0E0E0", "HIV-Specific TCR" = HIV_RED)

# CP003's Figure 3 palette. CP018 shares the labels that genuinely correspond
# (established by pseudobulk correlation in script 3b) and adds its own for the
# populations with no CP003 counterpart, so a shared label carries a shared
# colour across the two participants' panels.
FIG3_COLS <- c(
  "Naive/Bystander"     = "#2166AC",
  "Transitional Tem"    = "#1A9850",
  "Tpex"                = "#F4A261",
  "TEMRA/Effector"      = HIV_RED,
  "Cycling Effector"    = "#9D4EDD",
  "Tex"                 = "#6D4C41",
  # CP018-only labels
  "Naive (hypoxic)"     = "#7FB3D5",
  "Activated Stem-like" = "#F6C445",
  "MAIT"                = "#00A9A5"
)

################################################################################
# 1. Load, and build the interleaved real-time axis
################################################################################
cp18  <- qs_read(QS_ANNOTATED)
un18  <- qs_read(CROSSREF)
require_ct_canon(cp18, "CP018 stim object")
require_ct_canon(un18, "CP018 unstim crossref")

# Real months since ART initiation, so stim and unstim points share one axis.
CP018_TP <- c("1 months\n(pre-ART)", "2 months\n(HIV-stim)",
              "25 months", "42 months", "90 months\n(HIV-stim)")

# ── SPECIFICITY IS A PROPERTY OF THE RECEPTOR, NOT OF THE CELL ───────────────
# Figure_3.R filters CELLS carrying the HIV label. That silently truncates a
# lineage: CP018's largest persister (TRAV9-2, 15 cells at 2m) has exactly ONE
# cell at 90m, and a 1-cell clone fails the ">2 cells within a timepoint"
# expansion rule there - so that cell is labelled "Other" and the 90m end of a
# genuinely persistent clonotype DISAPPEARS from the alluvial.
#
# Measured: CP018 has 1 such split-label clonotype, CP003 has 30. Once a TCR is
# called antigen-specific anywhere, every cell carrying that receptor belongs to
# the same lineage. So promote the call to the clonotype and track all its cells.
# Both constructions are computed; the paper's is kept for comparability.
md18 <- cp18@meta.data
md18$.ct <- md18[[CT_CANON]]
md18 <- md18[!is.na(md18$.ct) & md18$.ct != "NA", ]
hiv_ct18 <- unique(md18$.ct[md18$HIV_Specific_TCR == "HIV-Specific TCR"])

stim18_paper <- md18 %>%
  filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
  mutate(tp = ifelse(Timepoint == "2m", "2 months\n(HIV-stim)",
                                        "90 months\n(HIV-stim)")) %>%
  group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

stim18 <- md18 %>%
  filter(.ct %in% hiv_ct18) %>%
  mutate(tp = ifelse(Timepoint == "2m", "2 months\n(HIV-stim)",
                                        "90 months\n(HIV-stim)")) %>%
  group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

n_recovered_by_promotion <- nrow(anti_join(stim18, stim18_paper,
                                           by = c("clone", "tp")))
cat("cell-filter rows:", nrow(stim18_paper),
    "| clonotype-level rows:", nrow(stim18),
    "| detections recovered:", n_recovered_by_promotion, "\n")

un18_lab <- c(CP018_1m = "1 months\n(pre-ART)", CP018_25m = "25 months",
              CP018_42m = "42 months")
mu18 <- un18@meta.data
mu18$.ct <- mu18[[CT_CANON]]
mu18 <- mu18[!is.na(mu18$.ct) & mu18$.ct != "NA", ]
hiv_ct18 <- union(hiv_ct18, unique(mu18$.ct[mu18$clone_match == "HIV-Specific Match"]))

unstim18 <- mu18 %>%
  filter(.ct %in% hiv_ct18) %>%
  mutate(tp = unname(un18_lab[as.character(orig.ident)])) %>%
  group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

unstim18_paper <- mu18 %>%
  filter(clone_match == "HIV-Specific Match") %>%
  mutate(tp = unname(un18_lab[as.character(orig.ident)])) %>%
  group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

# Re-derive the stim side against the FULL clonotype set (the union above can
# add a clonotype first called in the unstimulated data).
stim18 <- md18 %>%
  filter(.ct %in% hiv_ct18) %>%
  mutate(tp = ifelse(Timepoint == "2m", "2 months\n(HIV-stim)",
                                        "90 months\n(HIV-stim)")) %>%
  group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

allu18_paper <- bind_rows(stim18_paper, unstim18_paper) %>%
  filter(!is.na(clone), !is.na(tp)) %>%
  mutate(tp = factor(tp, levels = CP018_TP))

allu18 <- bind_rows(stim18, unstim18) %>%
  filter(!is.na(clone), !is.na(tp)) %>%
  mutate(tp = factor(tp, levels = CP018_TP))

cat("CP018 HIV-specific clonotypes:", n_distinct(stim18$clone),
    "| recovered in unstim:", n_distinct(unstim18$clone), "\n")

################################################################################
# 2. Panel C — alluvial, interleaved timepoints (mirrors Fig3C)
################################################################################
multi18 <- allu18 %>% group_by(clone) %>%
  filter(n_distinct(tp) >= 2) %>% ungroup()
save_csv(multi18, file.path(csv_f3, "Fig3C_CP018_alluvial_source.csv"))

# ── geom_flow ONLY CONNECTS ADJACENT AXES ────────────────────────────────────
# A clonotype absent at an intermediate timepoint has no stratum there, so its
# ribbon BREAKS and the later timepoint renders as a disconnected island. That
# is a plotting artefact, not absence of persistence: CP018's TRAV9-2 lineage
# occupies 1m, 2m and 90m but not 25m or 42m, so the 2m->90m link was invisible
# even though the data contains it.
#
# Fix: complete the clonotype x timepoint grid with zeros and draw the ribbon
# with geom_alluvium, which renders ONE continuous band per lineage across all
# axes. A zero-height band is invisible, so non-detected timepoints are drawn at
# a thin FLOOR width - explicitly marked, since a visible ribbon must not be
# read as cells having been detected there.
# Ribbon thickness where a lineage is NOT detected. Kept small because this
# value now also sets the (invisible) stratum height at those timepoints, and a
# large floor would visibly inflate the stack.
ALLU_FLOOR <- 0.18

make_alluvial <- function(dat, tp_levels, title, participant,
                          bridge_gaps = TRUE) {
  if (!nrow(dat)) {
    return(ggplot() + annotate("text", 1, 1, size = 12,
             label = paste0(participant, ": no clonotype at \u22652 timepoints")) +
             theme_void())
  }
  ord <- dat %>% group_by(clone) %>% summarise(tot = sum(n_cells), .groups = "drop") %>%
    arrange(desc(tot))
  n_cl <- nrow(ord)
  pal  <- colorRampPalette(brewer.pal(max(3, min(n_cl, 12)), "Set3"))(n_cl)
  names(pal) <- ord$clone

  if (!bridge_gaps) {
    # Paper-faithful rendering, kept for comparison.
    return(
      ggplot(dat, aes(x = tp, y = n_cells,
                      stratum = clone, alluvium = clone, fill = clone)) +
        geom_stratum(width = 0.3, colour = "grey30", linewidth = 0.3) +
        geom_flow(alpha = 0.5, width = 0.3) +
        scale_fill_manual(values = pal, guide = "none") +
        scale_x_discrete(expand = c(0.08, 0.08), drop = FALSE, limits = tp_levels) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_fig3() +
        theme(axis.text.x = element_text(size = FIG3_BASE - 4, hjust = 0.5,
                                         lineheight = 0.9, face = "bold")) +
        labs(x = NULL, y = "Number of cells", title = title,
             subtitle = paste0(n_cl, " clonotypes at \u22652 timepoints",
                               "  (ribbons break at non-detected timepoints)")))
  }

  # Complete grid: every lineage gets a row at every timepoint.
  full <- expand.grid(clone = ord$clone, tp = tp_levels,
                      stringsAsFactors = FALSE) %>%
    left_join(dat %>% mutate(tp = as.character(tp)) %>%
                select(clone, tp, n_cells),
              by = c("clone", "tp")) %>%
    mutate(n_cells   = ifelse(is.na(n_cells), 0, n_cells),
           detected  = n_cells > 0,
           y_ribbon  = pmax(n_cells, ALLU_FLOOR),
           tp        = factor(tp, levels = tp_levels),
           clone     = factor(clone, levels = ord$clone))

  # CRITICAL: both layers must stack from the SAME data, or their y-orders
  # diverge. Drawing strata from a filtered copy (detected rows only) while the
  # alluvium stacks the complete grid puts a lineage's ribbon and its bar at
  # different heights - the colour appears to jump between timepoints and no
  # longer identifies the lineage. So geom_stratum reads `full` too, and the
  # undetected rows are made invisible by mapping them to zero linewidth/alpha
  # rather than by removing them.
  full <- full %>% mutate(
    bar_alpha = ifelse(detected, 1, 0),
    bar_line  = ifelse(detected, 0.35, 0))

  ggplot(full, aes(x = tp, alluvium = clone, fill = clone,
                   stratum = clone, y = y_ribbon)) +
    # `decreasing` MUST be left at its default on both layers. Setting it to
    # FALSE here re-sorted the ribbons by cell count at each timepoint, while
    # the bars kept the fixed lineage order - so a ribbon arrived at a different
    # height from its own bar (measured: 733 px off at 90 months) and the
    # colours appeared to swap. Default ordering is by stratum factor level,
    # which both layers share.
    geom_alluvium(alpha = 0.5, width = 0.3, knot.pos = 0.28) +
    geom_stratum(aes(alpha = bar_alpha, linewidth = bar_line),
                 width = 0.3, colour = "grey25") +
    scale_alpha_identity() +
    scale_linewidth_identity() +
    # Open ticks mark timepoints a lineage passes THROUGH without detection.
    # No per-timepoint "not detected" marker here: several lineages share one
    # timepoint, so their markers collide and cannot be attributed. The lineage
    # trace panel (Fig3C2) carries that information unambiguously, one row per
    # lineage. Here the thin floor ribbon is the only gap indicator.
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(expand = c(0.08, 0.08), drop = FALSE, limits = tp_levels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_fig3() +
    theme(axis.text.x = element_text(size = FIG3_BASE - 4, hjust = 0.5,
                                     lineheight = 0.9, face = "bold"),
          plot.subtitle = element_text(size = FIG3_BASE - 12, colour = HIV_RED,
                                       face = "bold", lineheight = 1.05),
          plot.margin = margin(12, 26, 12, 12)) +
    labs(x = NULL, y = "Number of cells", title = title,
         # Short enough not to clip; the mechanism is documented in the
         # script header and the trace panel shows it per lineage.
         subtitle = paste0(n_cl, " clonotypes detected at \u22652 timepoints\n",
                           "thin ribbon = lineage not detected there"))
}

panelC18 <- make_alluvial(multi18, CP018_TP,
                          paste0("Three HIV-specific clonotypes carry\n",
                                 "through 90 months of suppression"), "CP018")
save_png(file.path(FIG3M, "Fig3C_CP018_Alluvial.png"), panelC18,
         width = 17, height = 13, dpi = 400)

################################################################################
# 2b. LINEAGE TRACE — the gap-proof view of the same data
################################################################################
# An alluvial has to decide what to draw where a lineage is not detected, and any
# choice can mislead. A trace plot has no such problem: one line per lineage
# across the full time axis, a filled point where the clonotype WAS detected
# (sized by cell count), an open point where it was not. Nothing is implied
# about the gap beyond "not detected in this sample".
make_trace <- function(dat, tp_levels, title, subtitle) {
  if (!nrow(dat)) return(ggplot() + theme_void())
  ord <- dat %>% group_by(clone) %>% summarise(tot = sum(n_cells), .groups = "drop") %>%
    arrange(desc(tot))
  full <- expand.grid(clone = ord$clone, tp = tp_levels, stringsAsFactors = FALSE) %>%
    left_join(dat %>% mutate(tp = as.character(tp)) %>% select(clone, tp, n_cells),
              by = c("clone", "tp")) %>%
    mutate(n_cells  = ifelse(is.na(n_cells), 0, n_cells),
           detected = n_cells > 0,
           tp       = factor(tp, levels = tp_levels),
           # Short CDR3 label; keep BOTH chains when a clonotype is dual-chain.
           lab      = factor(substr(clone, 1, 26), levels = substr(ord$clone, 1, 26)))
  pal <- colorRampPalette(brewer.pal(max(3, min(nrow(ord), 8)), "Dark2"))(nrow(ord))
  names(pal) <- levels(full$lab)

  ggplot(full, aes(x = tp, y = lab, group = lab, colour = lab)) +
    geom_line(linewidth = 1.4, alpha = 0.55, orientation = "y") +
    geom_point(data = full %>% filter(!detected), shape = 1, size = 5,
               stroke = 1.2, colour = "grey55") +
    geom_point(data = full %>% filter(detected), aes(size = n_cells)) +
    geom_text(data = full %>% filter(detected), aes(label = n_cells),
              colour = "white", size = 4.6, fontface = "bold") +
    scale_size_continuous(range = c(7, 20), guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    theme_fig3(base_size = 24) +
    theme(axis.text.x = element_text(size = 20, hjust = 0.5, lineheight = 0.9,
                                     face = "bold"),
          axis.text.y = element_text(family = "mono", size = 17),
          panel.grid.major.y = element_line(colour = "grey92")) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle)
}

trace18 <- make_trace(multi18, CP018_TP, "CP018 lineage traces",
  "Filled point = detected (number = cells) | open circle = not detected in this sample")
save_png(file.path(FIG3M, "Fig3C2_CP018_LineageTrace.png"), trace18,
         width = 15, height = 7, dpi = 400)
save_csv(multi18 %>% arrange(clone, tp),
         file.path(csv_f3, "Fig3C2_CP018_lineage_trace_source.csv"))

################################################################################
# 3. Panel E — persistence histogram over BOTH experiments (mirrors Fig3E)
################################################################################
# Counted over ALL HIV-specific clonotypes, not only persistent ones - so the
# height of the "1 timepoint" bar is part of the result rather than omitted.
all18 <- unique(stim18$clone)
persist18 <- data.frame(clone = all18) %>%
  left_join(allu18 %>% group_by(clone) %>%
              summarise(n_timepoints = n_distinct(tp), .groups = "drop"),
            by = "clone") %>%
  mutate(n_timepoints = ifelse(is.na(n_timepoints), 0L, n_timepoints))
psum18 <- persist18 %>% count(n_timepoints) %>%
  mutate(n_timepoints = factor(n_timepoints))
save_csv(psum18, file.path(csv_f3, "Fig3E_CP018_persistence_source.csv"))

make_persist <- function(psum, n_total, title) {
  ggplot(psum, aes(x = n_timepoints, y = n)) +
    geom_col(fill = HIV_RED, width = 0.6) +
    geom_text(aes(label = n), vjust = -0.4, size = 12, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    theme_fig3() +
    labs(x = "Timepoints Detected", y = "N Clonotypes", title = title,
         subtitle = paste0(n_total, " HIV-specific clonotypes"))
}
panelE18 <- make_persist(psum18, nrow(persist18), "Clonotype Persistence")
save_png(file.path(FIG3M, "Fig3E_CP018_Persistence.png"), panelE18,
         width = 12, height = 10, dpi = 400)

################################################################################
# 4. Panel F — unstim cluster composition, CD8 clusters only (mirrors Fig3F)
################################################################################
# Figure_3.R restricts to CD8 clusters: a few validated cells land in MAIT/NK
# clusters at the PBMC annotation level. Excluded here too, and the count of
# exclusions is reported rather than dropped silently.
CD8_CLUSTERS <- c("Naive 1 CD8", "Naive 2 CD8", "Naive Intermediate CD8",
                  "Tscm CD8", "TEMRA/CTL CD8", "Tex CD8")

val18_all <- un18@meta.data %>% filter(clone_match == "HIV-Specific Match")
val18_cd8 <- val18_all %>% filter(as.character(Annotation) %in% CD8_CLUSTERS)
cat("CP018 validated cells:", nrow(val18_all),
    "| in CD8 clusters:", nrow(val18_cd8),
    "| excluded (non-CD8):", nrow(val18_all) - nrow(val18_cd8), "\n")

# FIXED COLOUR PER CLUSTER, not per rank. Assigning fill by rank order means
# the same cluster gets a different colour in each panel - which reverses the
# visual reading of the comparison. Also allows a shared y limit so bar LENGTH
# is comparable between panels rather than each being auto-scaled.
TARA_CLUSTER_COLS <- c(
  "TEMRA/CTL CD8"          = "#00A9A5",
  "Tex CD8"                = "#E15759",
  "Tscm CD8"               = "#F6C445",
  "Naive 1 CD8"            = "#2166AC",
  "Naive 2 CD8"            = "#7FB3D5",
  "Naive Intermediate CD8" = "#A6CEE3"
)

make_clusterbar <- function(val_cd8, title, all_n, y_max = NA,
                            use_pct = FALSE) {
  if (!nrow(val_cd8)) {
    return(ggplot() + annotate("text", 1, 1, size = 12,
             label = "no validated cells in CD8 clusters") + theme_void())
  }
  d <- val_cd8 %>% count(Annotation = as.character(Annotation), name = "n") %>%
    mutate(pct = round(100 * n / sum(n), 1))
  # Order by the fixed vocabulary, so the same cluster sits in the same row
  # position in both panels.
  lev <- rev(intersect(names(TARA_CLUSTER_COLS), d$Annotation))
  d$Annotation <- factor(d$Annotation, levels = lev)
  yv <- if (use_pct) "pct" else "n"
  ylab <- if (use_pct) "% of validated CD8 cells" else "N Cells"

  # A bar near the axis maximum pushes its outside label past the panel edge and
  # the text gets clipped mid-string. Place the label INSIDE for long bars.
  d$inside <- d[[yv]] > 0.62 * (if (!is.na(y_max)) y_max else max(d[[yv]]))
  p <- ggplot(d, aes(x = Annotation, y = .data[[yv]], fill = Annotation)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = paste0(n, " (", pct, "%)"),
                  hjust = ifelse(inside, 1.06, -0.08),
                  colour = inside),
              size = 10, fontface = "bold", show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_manual(values = TARA_CLUSTER_COLS, drop = FALSE) +
    coord_flip() + theme_fig3() +
    labs(x = NULL, y = ylab, title = title,
         subtitle = paste0(nrow(val_cd8), " CD8+ cells from ",
                           n_distinct(val_cd8[[CT_CANON]]), " validated clonotypes"))
  if (!is.na(y_max)) {
    p + scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02)))
  } else {
    p + scale_y_continuous(expand = expansion(mult = c(0, 0.42)))
  }
}
panelF18 <- make_clusterbar(val18_cd8, "Cluster Distribution in TARA", nrow(val18_all))
save_png(file.path(FIG3M, "Fig3F_CP018_TARA_Clusters.png"), panelF18,
         width = 14, height = 10, dpi = 400)
save_csv(val18_cd8 %>% count(Annotation = as.character(Annotation), name = "n_cells"),
         file.path(csv_f3, "Fig3F_CP018_cluster_source.csv"))

################################################################################
# 4b. RECONCILING THE TWO "Tex" LABELS — they are not the same measurement
################################################################################
# Panel F says CP018's validated cells are mostly "Tex CD8". Elsewhere this
# project reports that CP018's own FACT-Seq clustering has NO exhausted
# population. Both are true, because the two labels come from DIFFERENT
# clusterings with different namespaces:
#
#   "Tex CD8"  - a cluster in the PBMC-atlas annotation of the UNSTIMULATED
#                eCITE-seq data, defined across ALL TARA participants.
#   "Tex"      - a cluster in CP018's OWN FACT-Seq clustering (res 0.4).
#                CP018 has no such cluster; CP003 does.
#
# So a cell can sit in the atlas's "Tex CD8" cluster while its participant has
# no Tex cluster of its own. The atlas label is an assignment to a shared
# reference, not a claim about CP018's internal structure.
#
# The stronger check: is the atlas "Tex CD8" cluster exhausted by the SAME
# criterion used on the FACT-Seq clusters (exhaustion:stemness > 1.5)? Measure
# it rather than trusting the name.
DefaultAssay(un18) <- "RNA"
if (length(Layers(un18[["RNA"]])) > 3) un18 <- JoinLayers(un18, assay = "RNA")
Xu <- GetAssayData(un18, layer = "data")
EXH  <- intersect(c("TOX","PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1"), rownames(Xu))
STEM <- intersect(c("TCF7","IL7R","SELL","CCR7","ID3","LEF1","BACH2"),        rownames(Xu))

atlas_scores <- do.call(rbind, lapply(CD8_CLUSTERS, function(k) {
  m <- !is.na(un18$Annotation) & as.character(un18$Annotation) == k
  if (sum(m) < 5) return(NULL)
  ex <- mean(Matrix::colMeans(Xu[EXH,  m, drop = FALSE]))
  st <- mean(Matrix::colMeans(Xu[STEM, m, drop = FALSE]))
  data.frame(atlas_cluster = k, n_cells = sum(m),
             exhaustion = round(ex, 3), stemness = round(st, 3),
             exh_stem_ratio = round(ex / (st + 0.01), 2),
             meets_Tex_criterion = ex / (st + 0.01) > 1.5,
             n_validated_here = sum(m & un18$clone_match == "HIV-Specific Match"),
             stringsAsFactors = FALSE)
}))
atlas_scores <- atlas_scores[order(-atlas_scores$exh_stem_ratio), ]
save_csv(atlas_scores, file.path(csv_f3, "Atlas_Tex_cluster_exhaustion_check.csv"))
cat("\n=== Is the unstim atlas 'Tex CD8' cluster exhausted by our own criterion? ===\n")
print(atlas_scores, row.names = FALSE)

tex_row <- atlas_scores[atlas_scores$atlas_cluster == "Tex CD8", ]
if (nrow(tex_row) && !tex_row$meets_Tex_criterion) {
  cat("\nIMPORTANT: the atlas 'Tex CD8' cluster does NOT meet the exhaustion\n",
      "criterion applied to the FACT-Seq clusters (ratio ", tex_row$exh_stem_ratio,
      " vs threshold 1.5). It is the most exhausted-LEANING CD8 cluster in the\n",
      "atlas, but it retains stemness. Do NOT read panel F as 'CP018's\n",
      "antigen-specific cells are terminally exhausted'. The defensible\n",
      "statement is that they occupy the atlas's most differentiated,\n",
      "highest-inhibitory-receptor CD8 compartment.\n", sep = "")
}

p_atlas <- ggplot(atlas_scores,
                  aes(x = reorder(atlas_cluster, exh_stem_ratio), y = exh_stem_ratio,
                      fill = meets_Tex_criterion)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.68) +
  geom_hline(yintercept = 1.5, linetype = "dashed", linewidth = 0.9, colour = "grey25") +
  # Label placed LEFT of the line: to its right there is no room before the
  # panel edge and the text was being clipped mid-word.
  annotate("text", x = 0.7, y = 1.45, hjust = 1, vjust = 0.4, size = 6.5,
           colour = "grey25", fontface = "italic",
           label = "terminal-exhaustion\nthreshold (1.5)") +
  geom_text(aes(label = paste0(exh_stem_ratio, "  (n=", n_cells,
                               ", validated=", n_validated_here, ")")),
            hjust = -0.06, size = 7, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72"), guide = "none") +
  scale_y_continuous(limits = c(0, 2.2), expand = expansion(mult = c(0, 0.05))) +
  coord_flip(clip = "off") + theme_fig3(base_size = 24) +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  labs(x = NULL, y = "Exhaustion : stemness ratio",
       title = "Atlas CD8 clusters in the unstimulated data",
       subtitle = paste0("No atlas CD8 cluster crosses the terminal-exhaustion ",
                         "threshold, including the one named 'Tex CD8'"))
save_png(file.path(FIG3M, "Atlas_Tex_exhaustion_check.png"), p_atlas,
         width = 15, height = 8, dpi = 350)

################################################################################
# 5. Panel D — functional gene profile of validated cells (mirrors Fig3D)
################################################################################
FUNC_GENES <- c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7",
                "IFNG", "CCL3", "CCL4", "CCL5",
                "TBX21", "EOMES", "PRDM1",
                "CX3CR1", "KLRG1", "FAS",
                "PDCD1", "LAG3", "HAVCR2", "TOX")
CAT_MAP <- c(rep("Cytotoxicity", 5), rep("Cytokines/Chemokines", 4),
             rep("Effector TFs", 3), rep("Effector Markers", 3),
             rep("Exhaustion", 4))
names(CAT_MAP) <- FUNC_GENES
CAT_COLS <- c("Cytotoxicity" = "#E15759", "Cytokines/Chemokines" = "#FF9DA7",
              "Effector TFs" = "#F28E2B", "Effector Markers" = "#EDC948",
              "Exhaustion" = "#9C755F")

gene_profile <- function(obj, cells, label) {
  DefaultAssay(obj) <- "RNA"
  if (length(Layers(obj[["RNA"]])) > 3) obj <- JoinLayers(obj, assay = "RNA")
  X <- GetAssayData(obj, layer = "data")
  g <- intersect(FUNC_GENES, rownames(X))
  do.call(rbind, lapply(g, function(gg) {
    e <- X[gg, cells]
    data.frame(gene = gg, pct_expressing = round(100 * mean(e > 0), 1),
               mean_expr = round(mean(e), 3), participant = label,
               stringsAsFactors = FALSE)
  }))
}

val18_bc <- rownames(val18_all)
gs18 <- gene_profile(un18, val18_bc, "CP018")
gs18$category <- factor(CAT_MAP[gs18$gene], levels = names(CAT_COLS))
gs18$gene <- factor(gs18$gene, levels = rev(FUNC_GENES))
save_csv(gs18, file.path(csv_f3, "Fig3D_CP018_functional_source.csv"))

make_funcplot <- function(gs, n_cells, title) {
  ggplot(gs, aes(x = gene, y = pct_expressing, fill = category)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = CAT_COLS, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, 100)) +
    coord_flip() + theme_fig3(base_size = 30) +
    theme(legend.position = "right", legend.text = element_text(size = 22)) +
    labs(x = NULL, y = "% of cells expressing", title = title,
         subtitle = paste0(n_cells, " validated cells"))
}
panelD18 <- make_funcplot(gs18, length(val18_bc), "Functional Profile at Rest")
save_png(file.path(FIG3M, "Fig3D_CP018_Functional_Profile.png"), panelD18,
         width = 13, height = 14, dpi = 400)

################################################################################
# 6. Panels B1 / B2 — UMAPs in the manuscript's style
################################################################################
red18 <- if ("umap.mnn.rna" %in% Reductions(cp18)) "umap.mnn.rna" else
         grep("umap", Reductions(cp18), value = TRUE)[1]
cp18$Fig_Annotation <- factor(as.character(cp18$Fig_Annotation),
                              levels = intersect(names(FIG3_COLS),
                                                 unique(as.character(cp18$Fig_Annotation))))

pB1 <- DimPlot(cp18, reduction = red18, group.by = "Fig_Annotation",
               cols = FIG3_COLS, pt.size = 0.5, label = FALSE) +
  labs(title = "CP018 annotated clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_fig3(base_size = 30) +
  theme(legend.text = element_text(size = 22), legend.position = "right",
        axis.text = element_blank(), axis.ticks = element_blank())
save_png(file.path(FIG3M, "Fig3B1_CP018_UMAP_Annotated.png"), pB1,
         width = 13, height = 10, dpi = 400)

pB2 <- DimPlot(cp18, reduction = red18, group.by = "HIV_Specific_TCR",
               cols = hiv_cols, pt.size = 0.5, order = "HIV-Specific TCR") +
  labs(title = "HIV-specific TCR", x = "UMAP 1", y = "UMAP 2") +
  theme_fig3(base_size = 30) +
  theme(legend.text = element_text(size = 24), legend.position = "bottom",
        axis.text = element_blank(), axis.ticks = element_blank())
save_png(file.path(FIG3M, "Fig3B2_CP018_UMAP_HIV_Specific.png"), pB2,
         width = 11, height = 11, dpi = 400)

################################################################################
# 7. COMBINED: CP003 beside CP018, same construction on both
################################################################################
# This is the panel the comparison actually needs - the two participants drawn
# with identical code, so a difference in the picture is a difference in the
# data rather than in how it was plotted.
if (file.exists(QS_CP003_ANNOT_FIG3) && file.exists(TARA_CP003)) {
  cp03 <- qs_read(QS_CP003_ANNOT_FIG3)
  un03 <- qs_read(TARA_CP003)
  cp03 <- add_ct_canon(cp03, verbose = FALSE)
  un03 <- add_ct_canon(un03, verbose = FALSE)

  CP003_TP <- c("2 months\n(HIV-stim)", "1 months", "12 months",
                "24 months", "101 months\n(HIV-stim)")

  # Same clonotype-level promotion as CP018, so any difference between the two
  # panels is a difference in the data, not in the rule applied.
  md03 <- cp03@meta.data; md03$.ct <- md03[[CT_CANON]]
  md03 <- md03[!is.na(md03$.ct) & md03$.ct != "NA", ]
  mu03 <- un03@meta.data; mu03$.ct <- mu03[[CT_CANON]]
  mu03 <- mu03[!is.na(mu03$.ct) & mu03$.ct != "NA", ]
  hiv_ct03 <- union(
    unique(md03$.ct[md03$HIV_Specific_TCR == "HIV-Specific TCR"]),
    unique(mu03$.ct[mu03$HIV_validated_status == "HIV-Specific (validated)"]))

  tp03_stim <- function(d) ifelse(d$Months == 2, "2 months\n(HIV-stim)",
                                                 "101 months\n(HIV-stim)")
  stim03_paper <- md03 %>% filter(HIV_Specific_TCR == "HIV-Specific TCR") %>%
    mutate(tp = tp03_stim(.)) %>%
    group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")
  unstim03_paper <- mu03 %>%
    filter(HIV_validated_status == "HIV-Specific (validated)") %>%
    mutate(tp = paste0(Age, " months")) %>%
    group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

  stim03 <- md03 %>% filter(.ct %in% hiv_ct03) %>%
    mutate(tp = tp03_stim(.)) %>%
    group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")
  unstim03 <- mu03 %>% filter(.ct %in% hiv_ct03) %>%
    mutate(tp = paste0(Age, " months")) %>%
    group_by(clone = .ct, tp) %>% summarise(n_cells = n(), .groups = "drop")

  allu03_paper <- bind_rows(stim03_paper, unstim03_paper) %>%
    filter(!is.na(clone), !is.na(tp)) %>%
    mutate(tp = factor(tp, levels = CP003_TP))
  allu03 <- bind_rows(stim03, unstim03) %>%
    filter(!is.na(clone), !is.na(tp)) %>%
    mutate(tp = factor(tp, levels = CP003_TP))
  multi03 <- allu03 %>% group_by(clone) %>% filter(n_distinct(tp) >= 2) %>% ungroup()
  multi03_paper <- allu03_paper %>% group_by(clone) %>%
    filter(n_distinct(tp) >= 2) %>% ungroup()

  # ── Paper-faithful reproduction, kept so the published number is verifiable ──
  multi18_paper <- allu18_paper %>% group_by(clone) %>%
    filter(n_distinct(tp) >= 2) %>% ungroup()
  rule_cmp <- data.frame(
    participant = c("CP003", "CP018"),
    persisting_cell_filter = c(n_distinct(multi03_paper$clone),
                               n_distinct(multi18_paper$clone)),
    persisting_clonotype_level = c(n_distinct(multi03$clone),
                                   n_distinct(multi18$clone)),
    stringsAsFactors = FALSE)
  save_csv(rule_cmp, file.path(csv_f3, "Persistence_rule_comparison.csv"))
  cat("\n=== Effect of the specificity rule on persistence ===\n")
  print(rule_cmp, row.names = FALSE)
  cat("The cell-filter column reproduces the published Fig3C subtitle (CP003 = 9),\n",
      "confirming the construction was read correctly; the clonotype-level column\n",
      "is the corrected count.\n", sep = "")

  save_png(file.path(CMB, "Combined_C2_LineageTrace_CP003_vs_CP018.png"),
           (make_trace(multi03, CP003_TP, "CP003 (never suppressed)",
                       "filled = detected (n cells) | open = not detected") /
            make_trace(multi18, CP018_TP, "CP018 (suppressed)",
                       "filled = detected (n cells) | open = not detected")) +
             plot_annotation(
               title = "Antigen-specific lineage traces across all timepoints",
               subtitle = paste0("Gap-proof view: a lineage absent from an ",
                                 "intermediate sample keeps its line, so a ",
                                 "first-to-last span is readable directly."),
               theme = theme(plot.title = element_text(size = 34, face = "bold"),
                             plot.subtitle = element_text(size = 22))),
           width = 19, height = 15, dpi = 320)

  save_png(file.path(CMB, "Combined_C_Alluvial_paperRule_CP003_vs_CP018.png"),
           (make_alluvial(multi03_paper, CP003_TP, "CP003 (paper's cell filter)",
                          "CP003", bridge_gaps = FALSE) |
            make_alluvial(multi18_paper, CP018_TP, "CP018 (paper's cell filter)",
                          "CP018", bridge_gaps = FALSE)) +
             plot_annotation(
               title = "Paper-faithful reproduction: specificity filtered per CELL",
               subtitle = paste0("Kept for comparability. A clonotype's low-count ",
                                 "timepoint drops out here - see the clonotype-level ",
                                 "version for the corrected view."),
               theme = theme(plot.title = element_text(size = 40, face = "bold"),
                             plot.subtitle = element_text(size = 24))),
           width = 34, height = 13, dpi = 300)
  save_csv(multi03, file.path(csv_f3, "Fig3C_CP003_alluvial_source.csv"))

  # ── Combined C: alluvials side by side ──
  pC03 <- make_alluvial(multi03, CP003_TP, "CP003 (never suppressed)", "CP003")
  pC18 <- make_alluvial(multi18, CP018_TP, "CP018 (suppressed)",      "CP018")
  save_png(file.path(CMB, "Combined_C_Alluvial_CP003_vs_CP018.png"),
           (pC03 | pC18) +
             plot_annotation(
               title = "HIV-specific clonotype persistence across stimulated and unstimulated timepoints",
               subtitle = paste0("Identical construction on both participants. ",
                                 "CP018's axis begins PRE-ART, which CP003's does not."),
               theme = theme(plot.title = element_text(size = 42, face = "bold"),
                             plot.subtitle = element_text(size = 28))),
           width = 32, height = 13, dpi = 300)

  # ── Combined E: persistence histograms ──
  all03 <- unique(stim03$clone)
  persist03 <- data.frame(clone = all03) %>%
    left_join(allu03 %>% group_by(clone) %>%
                summarise(n_timepoints = n_distinct(tp), .groups = "drop"),
              by = "clone") %>%
    mutate(n_timepoints = ifelse(is.na(n_timepoints), 0L, n_timepoints))
  psum03 <- persist03 %>% count(n_timepoints) %>% mutate(n_timepoints = factor(n_timepoints))

  pE03 <- make_persist(psum03, nrow(persist03), "CP003 (never suppressed)")
  pE18 <- make_persist(psum18, nrow(persist18), "CP018 (suppressed)")
  save_png(file.path(CMB, "Combined_E_Persistence_CP003_vs_CP018.png"),
           (pE03 | pE18) +
             plot_annotation(title = "Clonotype persistence across both experiments",
               theme = theme(plot.title = element_text(size = 42, face = "bold"))),
           width = 26, height = 11, dpi = 300)

  # Persistence as PROPORTION, since the two differ in clonotype count
  # (CP003 179 vs CP018 72) - raw bars would mislead.
  prop_cmp <- bind_rows(
    psum03 %>% mutate(participant = "CP003 (never suppressed)",
                      pct = 100 * n / sum(n)),
    psum18 %>% mutate(participant = "CP018 (suppressed)",
                      pct = 100 * n / sum(n)))
  save_csv(prop_cmp, file.path(csv_f3, "Combined_E_persistence_proportion.csv"))

  pEp <- ggplot(prop_cmp, aes(x = n_timepoints, y = pct, fill = participant)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "black", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(pct, 1), "%\n(n=", n, ")")),
              position = position_dodge(width = 0.8), vjust = -0.25, size = 7,
              fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = c("CP003 (never suppressed)" = "#B2182B",
                                 "CP018 (suppressed)" = "#2166AC"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    theme_fig3(base_size = 26) +
    theme(legend.position = "top", legend.text = element_text(size = 22)) +
    labs(x = "Timepoints detected", y = "% of participant's HIV-specific clonotypes",
         title = "Persistence, normalised within participant",
         # Counts computed here, never hardcoded: the denominators depend on the
         # HIV-specific rule in force, and a stale literal would misstate them.
         subtitle = paste0("Proportions: CP003 has ", nrow(persist03),
                           " HIV-specific clonotypes vs CP018's ", nrow(persist18)))
  save_png(file.path(CMB, "Combined_E_Persistence_proportion.png"), pEp,
           width = 15, height = 10, dpi = 350)

  # ── Combined F: unstim cluster composition ──
  val03_all <- un03@meta.data %>% filter(HIV_validated_status == "HIV-Specific (validated)")
  val03_cd8 <- val03_all %>% filter(as.character(Annotation) %in% CD8_CLUSTERS)
  # Plotted as PERCENTAGE on a shared 0-100 axis. CP003 has 121 validated CD8
  # cells vs CP018's 7, so a raw count axis makes CP018's bars invisible and an
  # auto-scaled count axis makes the two panels silently non-comparable.
  pF03 <- make_clusterbar(val03_cd8, "CP003 (never suppressed)", nrow(val03_all),
                          y_max = 100, use_pct = TRUE)
  pF18 <- make_clusterbar(val18_cd8, "CP018 (suppressed)",      nrow(val18_all),
                          y_max = 100, use_pct = TRUE)
  save_png(file.path(CMB, "Combined_F_TARA_Clusters_CP003_vs_CP018.png"),
           (pF03 | pF18) +
             plot_annotation(
               title = "Resting-state cluster distribution of antigen-validated clones",
               subtitle = paste0("Percentage on a shared axis, same CD8 restriction ",
                                 "and same colour per cluster in both panels.\n",
                                 "CAUTION: CP018 n=", nrow(val18_cd8),
                                 " cells - percentages are unstable at this depth."),
               theme = theme(plot.title = element_text(size = 40, face = "bold"),
                             plot.subtitle = element_text(size = 24))),
           width = 30, height = 11, dpi = 300)

  comp_cmp <- bind_rows(
    val03_cd8 %>% count(Annotation = as.character(Annotation), name = "n") %>%
      mutate(participant = "CP003", pct = round(100 * n / sum(n), 1)),
    val18_cd8 %>% count(Annotation = as.character(Annotation), name = "n") %>%
      mutate(participant = "CP018", pct = round(100 * n / sum(n), 1)))
  save_csv(comp_cmp, file.path(csv_f3, "Combined_F_cluster_composition.csv"))
  cat("\n=== Resting cluster distribution, both participants ===\n")
  print(as.data.frame(comp_cmp), row.names = FALSE)

  # ── Combined D: functional profile, matched genes ──
  gs03 <- gene_profile(un03, rownames(val03_all), "CP003")
  gs03$category <- factor(CAT_MAP[gs03$gene], levels = names(CAT_COLS))
  gs_both <- bind_rows(gs03, gs18)
  gs_both$gene <- factor(gs_both$gene, levels = rev(FUNC_GENES))
  save_csv(gs_both, file.path(csv_f3, "Combined_D_functional_both.csv"))

  pDb <- ggplot(gs_both, aes(x = gene, y = pct_expressing, fill = participant)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    scale_fill_manual(values = c(CP003 = "#B2182B", CP018 = "#2166AC"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 100)) +
    coord_flip() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_fig3(base_size = 24) +
    theme(legend.position = "top", strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 18, face = "bold"),
          strip.background = element_blank()) +
    labs(x = NULL, y = "% of validated cells expressing",
         title = "Functional profile of antigen-validated clones at rest",
         subtitle = paste0("CP003 n=", nrow(val03_all),
                           " cells | CP018 n=", nrow(val18_all), " cells"))
  save_png(file.path(CMB, "Combined_D_Functional_CP003_vs_CP018.png"), pDb,
           width = 15, height = 18, dpi = 350)

  # ── Headline numbers, both participants, one table ──
  summ <- data.frame(
    participant = c("CP003", "CP018"),
    ART_status  = c("Never suppressed", "Suppressed"),
    stim_timepoints = c("2m, 101m", "2m, 90m"),
    unstim_timepoints = c("1m, 12m, 24m", "1m (pre-ART), 25m, 42m"),
    hiv_specific_clonotypes = c(length(all03), length(all18)),
    recovered_in_unstim = c(n_distinct(unstim03$clone), n_distinct(unstim18$clone)),
    pct_recovered = c(round(100 * n_distinct(unstim03$clone) / length(all03), 1),
                      round(100 * n_distinct(unstim18$clone) / length(all18), 1)),
    persisting_2plus_tp = c(n_distinct(multi03$clone), n_distinct(multi18$clone)),
    pct_persisting = c(round(100 * n_distinct(multi03$clone) / length(all03), 1),
                       round(100 * n_distinct(multi18$clone) / length(all18), 1)),
    # TWO different denominators, both reported so neither can be quoted as the
    # other: all-lineage validated cells, and the CD8-restricted subset that
    # panel F actually plots (Figure_3.R excludes MAIT/NK at the PBMC level).
    validated_cells_all_lineages = c(nrow(val03_all), nrow(val18_all)),
    validated_cells_CD8_only     = c(nrow(val03_cd8), nrow(val18_cd8)),
    excluded_non_CD8 = c(nrow(val03_all) - nrow(val03_cd8),
                         nrow(val18_all) - nrow(val18_cd8)),
    stringsAsFactors = FALSE)
  save_csv(summ, file.path(csv_f3, "Combined_headline_numbers.csv"))
  cat("\n=== Headline numbers, matched construction ===\n")
  print(as.data.frame(summ), row.names = FALSE)

  rm(cp03, un03); gc(verbose = FALSE)
} else {
  cat("\nCP003 objects not found - CP018-only panels written, combined skipped.\n")
}

cat("\nOutputs:\n  ", FIG3M, "\n  ", CMB, "\n")
cat("\nNOTE: the alluvial axis interleaves stimulated and unstimulated\n",
    "timepoints in real time, matching Figure 3C. Persistence measured on\n",
    "stim timepoints ALONE is far lower and is not what the manuscript shows.\n",
    sep = "")
