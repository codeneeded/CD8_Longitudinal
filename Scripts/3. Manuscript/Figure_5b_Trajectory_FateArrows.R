################################################################################
# Figure 5, ALTERNATIVE trajectory panel — fate arrows on the pseudotime UMAP
#
# WHAT THIS IS
#   A second, cleaner rendering of the SAME monocle3 trajectory already shown in
#   Fig5A, restyled after the miR-29a PDAC panel: cells coloured by pseudotime,
#   an arrow per branch, the root marked in red, boxed state labels, and corner
#   axis arrows instead of an axis frame.
#
#   Fig5A_Monocle3_pseudotime_UMAP.png is NOT touched. This is an extra file.
#
# WHAT IS NOT RE-COMPUTED  (this is the important part)
#   Nothing. No trajectory is re-inferred, no clustering or annotation is
#   redone, slingshot is not run. The pseudotime plotted here is read straight
#   from the existing
#     Analysis/CD8_subset/07_trajectory/monocle3/Monocle3_pseudotime_per_cell.csv
#   and the coordinates are the object's own wnn.umap. The panel therefore shows
#   the same numbers as Fig5A.
#
# WHAT THE ARROWS ARE, AND WHAT THEY ARE NOT
#   They are a drawing of the ORDERING the panel is coloured by: annotated
#   states are placed at their median UMAP position, joined by a minimum
#   spanning tree in UMAP space, rooted at the lowest-pseudotime state, and each
#   edge is oriented from lower to higher median pseudotime.
#
#   They are NOT monocle3's principal graph. That graph exists only inside the
#   cds object, which the pipeline does not save (4.CD8_Subseting.R builds it at
#   line 1316 and keeps only the per-cell pseudotime), so drawing it would mean
#   re-running learn_graph()/order_cells(). Say "arrows summarise the monocle3
#   pseudotime ordering across states", not "the monocle3 graph".
#
# TWO MISTAKES THIS FORMAT MADE IN THE miR-29a PAPER, AVOIDED HERE
#   1. Arrow direction was taken from the trajectory tool's per-curve order
#      rather than from the pseudotime the panel is COLOURED by. Those can
#      disagree, and an arrow once rendered running backwards against its own
#      colour scale. Here every edge is oriented by the displayed pseudotime and
#      a guard stops the script if any branch would run backwards.
#   2. Waypoints were re-sorted by median pseudotime WITHIN a branch, which
#      reorders nodes because a per-state median can fall between two other
#      states' medians - the arrow then zig-zags instead of following its
#      branch. Here the branch order is the tree's own root-to-leaf path.
#
# OUT  Manuscript/Figure 5/Fig5A_alt_FateArrows_pseudotime_UMAP.png
#      Manuscript/Figure 5/Fig5A_alt_FateArrows_waypoints.csv
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2)
  library(ggplot2); library(dplyr)
  library(ggnewscale); library(ggrepel); library(viridis)
})

# Export the data behind every panel (Editor point E-3).
.sd_shim_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) dirname(normalizePath(f[1])) else {
    d <- NULL
    for (i in rev(seq_len(sys.nframe()))) {
      o <- sys.frame(i)$ofile
      if (!is.null(o)) { d <- dirname(normalizePath(o)); break }
    }
    if (is.null(d)) getwd() else d
  }
})
.sd_shim <- file.path(.sd_shim_dir, "00.Source_Data_Shim.R")
if (file.exists(.sd_shim)) source(.sd_shim) else
  warning("source-data shim not found at ", .sd_shim,
          "; panels will not export", call. = FALSE, immediate. = TRUE)
CD8_ROOT <- if (exists("sd_repo_root")) sd_repo_root(.sd_shim_dir) else
  Sys.getenv("CD8_BASE", unset = normalizePath(
    file.path(.sd_shim_dir, "..", ".."), mustWork = FALSE))

base_dir     <- Sys.getenv("CD8_BASE", unset = CD8_ROOT)
saved_dir    <- file.path(base_dir, "saved_R_data")
fig5_dir     <- file.path(base_dir, "Manuscript", "Figure 5")
analysis_dir <- file.path(base_dir, "Analysis", "CD8_subset")
dir.create(fig5_dir, recursive = TRUE, showWarnings = FALSE)

SEED <- 42
set.seed(SEED)

# Same palette as Figure_4.R, so a state keeps one colour across the figure set.
cluster_colors <- c(
  "Naïve CD8"              = "#2D5F8A",
  "Naïve CD8 2"            = "#1B9E77",
  "Naïve CD8 3"            = "#7570B3",
  "Naïve CD8 4"            = "#6BAED6",
  "Tscm CD8"               = "#2CA02C",
  "Naïve Intermediate CD8" = "#98D8C8",
  "TEM CD8"                = "#D95F02",
  "TEMRA CD8"              = "#E6AB02",
  "CD16+ Effector CD8"     = "#D62728",
  "KIR+ innate-like CD8"   = "#8C564B"
)
ROOT_STATE <- "Naïve CD8"

# ---------------------------------------------------------------- 1) inputs ---
TARA_cd8 <- qs_read(file.path(saved_dir, "TARA_cd8_HEI_annotated_final.qs2"))
cat("loaded:", ncol(TARA_cd8), "cells\n")

umap_coords <- Embeddings(TARA_cd8, reduction = "wnn.umap")

# Pseudotime comes from the EXISTING monocle3 run - nothing is recomputed.
pt_csv <- file.path(analysis_dir, "07_trajectory", "monocle3",
                    "Monocle3_pseudotime_per_cell.csv")
if (!file.exists(pt_csv))
  stop("monocle3 pseudotime not found at\n  ", pt_csv,
       "\nRun the trajectory step in \"2. Analysis/4.CD8_Subseting.R\" first.",
       call. = FALSE)
m3 <- read.csv(pt_csv, check.names = FALSE)
cat("monocle3 pseudotime read:", nrow(m3), "cells\n")

emb <- as.data.frame(umap_coords[, 1:2])
colnames(emb) <- c("UMAP1", "UMAP2")
emb$cell  <- rownames(emb)
emb$state <- as.character(TARA_cd8$CD8_Annotation)
emb$pt    <- m3$monocle3_pseudotime[match(emb$cell, m3$cell)]

n_pt <- sum(is.finite(emb$pt))
cat("cells with a pseudotime:", n_pt, "of", nrow(emb), "\n")
if (n_pt < 0.5 * nrow(emb))
  stop("fewer than half the cells carry a monocle3 pseudotime - the cell ",
       "barcodes in the CSV do not match this object", call. = FALSE)
stopifnot(ROOT_STATE %in% emb$state)

# --------------------------------------------- 2) states, ordered by the pt ---
ctr <- emb %>%
  filter(is.finite(pt)) %>%
  group_by(state) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
            .pt = median(pt), n = n(), .groups = "drop") %>%
  as.data.frame()
cat("\nstates by median monocle3 pseudotime:\n")
print(ctr %>% arrange(.pt) %>%
        transmute(state, n, median_pt = round(.pt, 3)), row.names = FALSE)

# ------------------------------------------------------- 3) the branch tree ---
# Each state attaches to its nearest neighbour in UMAP space among states that
# are SIGNIFICANTLY EARLIER in pseudotime (Wilcoxon on the per-cell values,
# p < ALPHA). Two reasons for the significance test rather than a plain
# nearest-earlier rule:
#
#   - It stops the panel asserting an order the data does not support. In this
#     dataset CD16+ Effector (median 24.27) and KIR+ innate-like (24.35) differ
#     by 0.08 with p = 0.93, so a plain rule chains KIR+ downstream of CD16+ and
#     the figure claims a differentiation step that is not there. With the test
#     they attach to TEMRA instead and read as alternative endpoints, which is
#     what the pseudotime actually shows. Every other consecutive pair here
#     separates at p < 1e-7.
#   - It makes each branch monotone in pseudotime by construction, so the
#     backwards-arrow guard below can only fire on a real problem.
#
# Set LINEAGES to override the structure entirely, e.g.
#   LINEAGES <- list(c("Naïve CD8", "TEM CD8", "TEMRA CD8", "CD16+ Effector CD8"))
ALPHA <- 0.05

# THE BRANCHES ARE THE MANUSCRIPT'S, NOT DERIVED HERE.
#
# Results paragraph on Figure 5A states them explicitly: "Three branches
# diverged from a shared root: an effector trajectory running from Naïve CD8
# through Naïve CD8 4 and TEM CD8 to TEMRA CD8 and terminating at the KIR3DL1+
# innate like and CD16+ Effector populations (pseudotime 3.1 -> 24.4), and two
# naïve branches ... The first terminated at Naïve CD8 2 ... The second
# terminated at Naïve CD8 3".
#
# WHY NOT DERIVE THEM. The auto-derived tree (set LINEAGES <- NULL to see it)
# routes the two naïve branches through Naïve Intermediate CD8 and Tscm CD8,
# because monocle3 puts Tscm at pseudotime 9.27 and Naïve CD8 2 / 3 at 9.50 and
# 12.68. Drawn as arrows that reads as stem-cell memory giving rise to naïve
# cells, which is backwards and is not what the manuscript claims. The marker
# programmes show why the ordering there cannot bear an arrow: across all six
# naïve/stem states the differentiated-minus-naïve score spans only 0.35
# (-1.586 to -1.236) and Tscm sits INSIDE that range, between Naïve CD8 2
# (-1.466) and Naïve CD8 3 (-1.236), whereas the one real transition, into TEM
# CD8, is a jump of +1.92. Within the naïve compartment monocle3 pseudotime is
# not tracking differentiation, so the manuscript's "shared root" is the
# defensible topology: both naïve branches leave Naïve CD8 directly.
#
# Naïve Intermediate CD8 and Tscm CD8 are on no branch in the manuscript's
# description, so they are drawn as labelled nodes without arrows.
#
# Every edge below is checked for monotone pseudotime by the guard further down.
LINEAGES <- list(
  c("Naïve CD8", "Naïve CD8 4", "TEM CD8", "TEMRA CD8", "CD16+ Effector CD8"),
  c("Naïve CD8", "Naïve CD8 4", "TEM CD8", "TEMRA CD8", "KIR+ innate-like CD8"),
  c("Naïve CD8", "Naïve CD8 2"),
  c("Naïve CD8", "Naïve CD8 3")
)

K <- nrow(ctr)
D <- as.matrix(dist(ctr[, c("UMAP1", "UMAP2")]))
dimnames(D) <- list(ctr$state, ctr$state)
pt_by_state <- split(emb$pt[is.finite(emb$pt)], emb$state[is.finite(emb$pt)])

root_i <- which(ctr$state == ROOT_STATE)
if (ctr$.pt[root_i] != min(ctr$.pt))
  warning("the root state ", ROOT_STATE, " is not the lowest-pseudotime state (",
          ctr$state[which.min(ctr$.pt)], " is). Arrows still start at the root.",
          call. = FALSE, immediate. = TRUE)

earlier_p <- function(a, b) {   # is a significantly earlier than b?
  suppressWarnings(wilcox.test(pt_by_state[[a]], pt_by_state[[b]],
                               alternative = "less")$p.value)
}
par <- rep(NA_integer_, K)
edge_log <- list()
for (j in order(ctr$.pt)) {
  if (j == root_i) next
  cand <- setdiff(which(ctr$.pt < ctr$.pt[j]), j)
  if (!length(cand)) { par[j] <- root_i; next }
  cand <- cand[order(D[j, cand])]          # nearest first
  chosen <- NA_integer_; pv <- NA_real_
  for (cc in cand) {
    p <- earlier_p(ctr$state[cc], ctr$state[j])
    if (is.finite(p) && p < ALPHA) { chosen <- cc; pv <- p; break }
  }
  if (is.na(chosen)) { chosen <- root_i; pv <- NA_real_ }
  par[j] <- chosen
  edge_log[[length(edge_log) + 1L]] <- data.frame(
    parent = ctr$state[chosen], child = ctr$state[j],
    d_umap = round(D[j, chosen], 3),
    d_median_pt = round(ctr$.pt[j] - ctr$.pt[chosen], 3),
    p_parent_earlier = signif(pv, 3), stringsAsFactors = FALSE)
}
cat("\ntree edges (parent must be significantly earlier, p <", ALPHA, "):\n")
print(bind_rows(edge_log), row.names = FALSE)

use_declared <- !is.null(LINEAGES)
if (use_declared) {
  miss <- setdiff(unlist(LINEAGES), ctr$state)
  if (length(miss))
    stop("LINEAGES names states that are not in this object: ",
         paste(miss, collapse = ", "), call. = FALSE)
  cat("\nusing the manuscript's declared branches (", length(LINEAGES), ")\n", sep = "")
  off <- setdiff(ctr$state, unlist(LINEAGES))
  if (length(off))
    cat("off-branch states, drawn as labelled nodes only: ",
        paste(off, collapse = ", "), "\n", sep = "")
}

# One arrow per branch, root to terminus.
if (use_declared) {
  poly <- bind_rows(lapply(LINEAGES, function(L) {
    q <- ctr[match(L, ctr$state), ]
    q$fate <- tail(L, 1)
    q$k <- seq_len(nrow(q))
    q
  }))
  FATE <- vapply(LINEAGES, function(L) tail(L, 1), "")
} else {
  leaves <- setdiff(seq_len(K), unique(par[!is.na(par)]))
  path_of <- function(i) { p <- i; while (!is.na(par[p[1]])) p <- c(par[p[1]], p); p }
  poly <- bind_rows(lapply(leaves, function(lf) {
    idx <- path_of(lf)
    q <- ctr[idx, ]
    q$fate <- ctr$state[lf]
    q$k <- seq_len(nrow(q))
    q
  }))
  FATE <- ctr$state[leaves]
}
poly$fate <- factor(poly$fate, levels = FATE)

# Guard: no branch may run backwards against the colour scale (mistake 1).
bad <- poly %>% group_by(fate) %>%
  filter(n() >= 2, any(diff(.pt) <= 0)) %>% pull(fate) %>% unique()
if (length(bad)) {
  cat("\nbranches that would run backwards in pseudotime:\n")
  print(poly %>% filter(fate %in% bad) %>%
          transmute(fate = as.character(fate), k, state,
                    median_pt = round(.pt, 3)), row.names = FALSE)
  stop("arrow would run backwards in pseudotime on: ",
       paste(bad, collapse = ", "),
       "\nThe nearest-neighbour tree disagrees with the pseudotime ordering ",
       "here; inspect the table above before forcing a drawing.", call. = FALSE)
}
cat("\nbranches drawn:", length(FATE), "\n")

# Arrowhead on the final segment only, so the shaft is not over-drawn.
tipp <- poly %>% group_by(fate) %>% slice_max(k, n = 2) %>% arrange(fate, k)

# Fate colour = the colour of the state the lineage ENDS in, so the arrow and
# its destination label share a colour and the reader decodes one key, not two.
LC <- setNames(unname(cluster_colors[FATE]), FATE)

# ------------------------------------------------------------------ 4) plot ---
p <- ggplot(emb[order(emb$pt), ], aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = pt), size = 0.55, alpha = 0.7, stroke = 0) +
  scale_colour_viridis_c(
    option = "inferno", name = "Pseudotime", na.value = "grey88",
    breaks = range(emb$pt, na.rm = TRUE),
    labels = c("naive", "terminal"),
    # The guide must be declared ON the scale: a guides() call placed before
    # new_scale_colour() targets an aesthetic the rename moves, and the
    # colourbar silently disappears.
    guide = guide_colourbar(barwidth = 0.8, barheight = 6, order = 1)) +
  ggnewscale::new_scale_colour() +
  geom_path(data = poly, aes(UMAP1, UMAP2, colour = fate, group = fate),
            inherit.aes = FALSE, linewidth = 1.25) +
  geom_path(data = tipp, aes(UMAP1, UMAP2, colour = fate, group = fate),
            inherit.aes = FALSE, linewidth = 1.25,
            arrow = arrow(length = unit(0.18, "inches"), type = "closed")) +
  geom_point(data = ctr,
             aes(UMAP1, UMAP2), inherit.aes = FALSE,
             size = 2.4, colour = "grey25") +
  # The root in red: every branch starts here, and the arrows only read as a
  # tree once that single origin is findable.
  #
  # WHITE HALO, NOT A BARE RED DOT. The house palette gives CD16+ Effector
  # "#D62728" and the miR-29a root marker is "#D7191C" - close enough that the
  # root read as another arrow head rather than as the start of the trajectory.
  # A white ring under a slightly smaller red disc separates the two without
  # changing the state palette, which has to stay consistent with Figure 4.
  geom_point(data = ctr[ctr$state == ROOT_STATE, ],
             aes(UMAP1, UMAP2), inherit.aes = FALSE, shape = 21,
             fill = "white", colour = "white", size = 8.5) +
  geom_point(data = ctr[ctr$state == ROOT_STATE, ],
             aes(UMAP1, UMAP2), inherit.aes = FALSE, shape = 21,
             fill = "#D7191C", colour = "grey15", size = 5.2, stroke = 1.1) +
  scale_colour_manual(values = LC, name = "Branch endpoint",
                      guide = guide_legend(order = 2,
                                           override.aes = list(linewidth = 1.8))) +
  geom_label_repel(data = ctr, aes(UMAP1, UMAP2, label = state),
                   inherit.aes = FALSE, colour = "grey15", fill = "white",
                   alpha = 0.9, size = 5.2, fontface = "bold",
                   label.size = 0.45, label.padding = unit(0.3, "lines"),
                   box.padding = 1.2, point.padding = 0.5,
                   segment.colour = "grey40", segment.size = 0.5,
                   min.segment.length = 0, max.overlaps = 30,
                   force = 5, force_pull = 0.2, seed = SEED) +
  theme_void(base_size = 13) +
  theme(legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(6, 6, 18, 18))

# Room beside the cloud for the repel solver.
xr <- diff(range(emb$UMAP1)); yr <- diff(range(emb$UMAP2))
p <- p + expand_limits(x = range(emb$UMAP1) + c(-1, 1) * 0.10 * xr,
                       y = range(emb$UMAP2) + c(-1, 1) * 0.10 * yr)

# Corner axis arrows, verbatim from Figure_4.R. Each arm is 0.12 of its OWN axis
# range, so they render equal only on a SQUARE canvas - hence 10 x 10 below.
af <- 0.12
x0 <- min(emb$UMAP1) - xr * 0.05
y0 <- min(emb$UMAP2) - yr * 0.08
p <- p +
  annotate("segment", x = x0, xend = x0 + xr * af, y = y0, yend = y0,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, colour = "black") +
  annotate("text", x = x0 + xr * af / 2, y = y0 - yr * 0.04,
           label = "UMAP1", size = 5, fontface = "bold") +
  annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + yr * af,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 0.9, colour = "black") +
  annotate("text", x = x0 - xr * 0.04, y = y0 + yr * af / 2,
           label = "UMAP2", size = 5, fontface = "bold", angle = 90) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig5_dir, "Fig5A_alt_FateArrows_pseudotime_UMAP.png"),
       plot = p, width = 10, height = 10, dpi = 320, bg = "white")

# The waypoint table IS the figure's own annotation layer, so it ships with it.
wp <- poly %>%
  transmute(fate = as.character(fate), step = k, state,
            UMAP1 = round(UMAP1, 4), UMAP2 = round(UMAP2, 4),
            median_pseudotime = round(.pt, 4))
write.csv(wp, file.path(fig5_dir, "Fig5A_alt_FateArrows_waypoints.csv"),
          row.names = FALSE)

cat("\nwrote Fig5A_alt_FateArrows_pseudotime_UMAP.png (10 x 10 in)\n")
cat("the original Fig5A_Monocle3_pseudotime_UMAP.png is untouched\n")
cat("\nlineages drawn:\n")
print(wp, row.names = FALSE)
if (exists("sd_report")) sd_report()
