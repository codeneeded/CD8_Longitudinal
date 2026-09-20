################################################################################
# Revision — the manuscript panels that reviewer asks would change
#
# WHY THIS EXISTS
#   Four responses said the change "requires re-rendering that panel" and left
#   it among the items still to complete. That pushes the work onto whoever
#   reads the response. These panels are rendered here instead, so each of
#   those four asks is answered with a figure rather than an intention.
#
# WHAT IT DOES NOT DO
#   Nothing here is applied to the manuscript. Each panel is an OPTION: it goes
#   into the response document, and the response states what the manuscript
#   change would be if the panel is adopted. The published panels are untouched.
#
# ASKS ANSWERED
#   R2-M1b  Figure 1C as a bubble plot: fraction expressing and level
#   R2-M1e  cluster frequency as a percentage of that sample's CD8 compartment
#   R2-M2a  Figure 2A clonal expansion on a CD8-only embedding
#   R2-m2   cell number annotated inside the UMAP panel, not only in the legend
#
# MARKER LISTS AND CLUSTER NAMES ARE READ OUT OF Figure_1.R
#   Re-typing them here would let the option panel drift from the published one,
#   and a bubble plot of different markers is not an answer to 1b. They are
#   parsed from the manuscript script so the two cannot disagree.
#
# OUT  Reviewer_Comments/Figures/ + Tables/  (never Manuscript/)
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2)
  library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork); library(ggrepel)
})

source("config_revision.R")

MS_DIR <- file.path(BASE, "Scripts", "3. Manuscript")
stopifnot(dir.exists(MS_DIR))
f1 <- readLines(file.path(MS_DIR, "Figure_1.R"), warn = FALSE)

# Pull a top-level `name <- c(...)` block out of Figure_1.R and evaluate it.
grab <- function(name) {
  i <- grep(paste0("^", name, "\\s*<-\\s*(setNames\\()?\\s*c\\("), f1)
  if (!length(i)) stop("cannot find ", name, " in Figure_1.R", call. = FALSE)
  i <- i[1]; depth <- 0L; out <- character()
  for (k in seq(i, length(f1))) {
    ln <- f1[k]
    bare <- gsub('"[^"]*"', "", ln)          # ignore parens inside strings
    out <- c(out, ln)
    depth <- depth + lengths(regmatches(bare, gregexpr("\\(", bare))) -
                     lengths(regmatches(bare, gregexpr("\\)", bare)))
    if (depth <= 0L && k > i - 1L) break
  }
  eval(parse(text = paste(out, collapse = "\n")))
}
cd8_cluster_names   <- grab("cd8_cluster_names")
adt_heatmap_markers <- grab("adt_heatmap_markers")
rna_heatmap_markers <- grab("rna_heatmap_markers")
cd8_short_labels    <- setNames(c("Naive 1", "Naive 2", "Naive Intermediate",
                                  "Tscm", "TEMRA/CTL", "Tex"), cd8_cluster_names)
cat("read from Figure_1.R:", length(cd8_cluster_names), "clusters,",
    length(adt_heatmap_markers), "ADT and", length(rna_heatmap_markers), "RNA markers\n")

theme_opt <- function(base = 20) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(face = "bold", size = base + 4),
          plot.subtitle = element_text(size = base - 2, colour = "grey25"),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = base - 2),
          legend.position = "right")
}

# ---------------------------------------------------------------- inputs ------
ALL <- qs_read(QS_ALL)
cat("whole object:", ncol(ALL), "cells\n")

# Figure 1 uses HEU and HUU entire, and HEI restricted to pre-ART entry.
keep <- (ALL$Condition %in% c("HEU", "HUU")) |
        (ALL$Condition == "HEI" & ALL$Timepoint_Group == "PreART_Entry")
SUB <- ALL[, keep]
SUB$Condition <- factor(SUB$Condition, levels = c("HUU", "HEU", "HEI"))
cat("Figure 1 subset:", ncol(SUB), "cells (published panel A states 79,705)\n")

CD8 <- SUB[, SUB$Annotation %in% cd8_cluster_names]
CD8$Annotation <- droplevels(factor(CD8$Annotation, levels = cd8_cluster_names))
cat("CD8-associated clusters:", ncol(CD8), "cells\n")

n_lab <- function(n) format(n, big.mark = ",")

# ================================================================ R2-M1b ======
# Figure 1C shows scaled group means as two heatmaps. The ask is a bubble plot
# carrying BOTH the fraction of cells expressing and the level. Same markers,
# same clusters, same order as the published panel.
bubble <- function(obj, assay, markers, title, subtitle) {
  DefaultAssay(obj) <- assay
  feats <- intersect(names(markers), rownames(obj))
  # unname() is load-bearing: Idents<- reads a named vector's names as cell
  # barcodes, warns "Cannot find cells provided" and keeps the old identities.
  Idents(obj) <- factor(unname(cd8_short_labels[as.character(obj$Annotation)]),
                        levels = unname(cd8_short_labels))
  stopifnot(identical(levels(Idents(obj)), unname(cd8_short_labels)))
  d <- DotPlot(obj, features = feats, assay = assay)$data
  d$module <- factor(markers[as.character(d$features.plot)],
                     levels = unique(markers[feats]))
  d$features.plot <- factor(as.character(d$features.plot), levels = rev(feats))
  ggplot(d, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    scale_size_area(name = "% of cells\nexpressing", max_size = 11,
                    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_colour_gradient2(name = "Scaled mean\nexpression",
                           low = "#2166AC", mid = "grey92", high = "#B2182B",
                           midpoint = 0) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_opt(19) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid.major = element_line(colour = "grey93"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1))
}

p_adt <- bubble(CD8, "ADT", adt_heatmap_markers,
  "Figure 1C as a bubble plot: surface protein",
  paste0("dot area = percentage of cells above zero after DSB normalisation\n",
         "colour = mean expression, scaled across clusters\n",
         n_lab(ncol(CD8)), " CD8-associated cells, ",
         length(cd8_cluster_names), " whole-PBMC clusters"))
ggsave(file.path(REV_FIG, "R2M1b_03_Fig1C_bubble_ADT.png"), plot = p_adt,
       width = 15, height = 13.5, dpi = 300, bg = "white")

p_rna <- bubble(CD8, "RNA", rna_heatmap_markers,
  "Figure 1C as a bubble plot: transcript",
  paste0("dot area = percentage of cells with a detected transcript\n",
         "colour = mean expression, scaled across clusters\n",
         n_lab(ncol(CD8)), " CD8-associated cells, ",
         length(cd8_cluster_names), " whole-PBMC clusters"))
ggsave(file.path(REV_FIG, "R2M1b_04_Fig1C_bubble_RNA.png"), plot = p_rna,
       width = 15, height = 13.5, dpi = 300, bg = "white")

bub_src <- bind_rows(
  transform(p_adt$data, modality = "surface protein (ADT)"),
  transform(p_rna$data, modality = "transcript (RNA)"))[
  , c("modality", "features.plot", "id", "pct.exp", "avg.exp", "avg.exp.scaled", "module")]
names(bub_src)[2:3] <- c("feature", "cluster")
write.csv(bub_src, file.path(REV_TAB, "R2M1b_03_Fig1C_bubble_source.csv"),
          row.names = FALSE)
cat("R2-M1b: two bubble panels,", nrow(bub_src), "feature-by-cluster rows\n")

# ================================================================ R2-M1e ======
# The published frequency panels divide by total PBMC. The ask is a percentage
# of the CD8 compartment, which removes the dependence on how much of each
# sample was CD8 in the first place. Both denominators are computed so the
# difference is visible rather than asserted.
md <- SUB@meta.data
md$is_cd8 <- md$Annotation %in% cd8_cluster_names
tot <- md %>% group_by(orig.ident, Condition) %>%
  summarise(n_pbmc = n(), n_cd8 = sum(is_cd8), .groups = "drop")

freq <- md %>% filter(is_cd8) %>%
  group_by(orig.ident, Condition, Annotation) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  complete(nesting(orig.ident, Condition), Annotation = cd8_cluster_names,
           fill = list(n_cells = 0)) %>%
  left_join(tot, by = c("orig.ident", "Condition")) %>%
  mutate(pct_of_cd8  = 100 * n_cells / n_cd8,
         pct_of_pbmc = 100 * n_cells / n_pbmc,
         cluster = factor(cd8_short_labels[as.character(Annotation)],
                          levels = cd8_short_labels))

# Wilcoxon across exposure groups, per cluster, BH-corrected across clusters.
pairs <- list(c("HUU", "HEU"), c("HUU", "HEI"), c("HEU", "HEI"))
st <- bind_rows(lapply(levels(freq$cluster), function(cl) {
  d <- freq[freq$cluster == cl, ]
  bind_rows(lapply(pairs, function(pr) {
    a <- d$pct_of_cd8[d$Condition == pr[1]]; b <- d$pct_of_cd8[d$Condition == pr[2]]
    if (length(a) < 3 || length(b) < 3) return(NULL)
    data.frame(cluster = cl, group1 = pr[1], group2 = pr[2],
               p = suppressWarnings(wilcox.test(a, b)$p.value))
  }))
}))
st$p_adj <- p.adjust(st$p, method = "BH")

p_pct <- ggplot(freq, aes(x = Condition, y = pct_of_cd8)) +
  geom_boxplot(aes(fill = Condition), outlier.shape = NA, alpha = 0.35,
               width = 0.62, linewidth = 0.5) +
  geom_jitter(width = 0.16, height = 0, size = 2.3, alpha = 0.85,
              colour = "grey15") +
  facet_wrap(~ cluster, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(HUU = "#4DAF4A", HEU = "#377EB8", HEI = "#E41A1C"),
                    guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.16))) +
  labs(x = NULL, y = "Percentage of that sample's CD8 compartment",
       title = "CD8 cluster frequency as a percentage of total CD8",
       subtitle = paste0("one point per sample (n = ", nrow(tot),
                         "); box = median and IQR. Denominator is each sample's ",
                         "own CD8 count\n",
                         "no pairwise comparison reaches BH-adjusted p < 0.05 ",
                         "(minimum p_adj = ", sprintf("%.2f", min(st$p_adj)), ")")) +
  theme_opt(19)
ggsave(file.path(REV_FIG, "R2M1e_02_Fig1E_pct_of_CD8.png"), plot = p_pct,
       width = 16, height = 10, dpi = 300, bg = "white")

write.csv(freq %>% select(orig.ident, Condition, cluster, n_cells, n_cd8, n_pbmc,
                          pct_of_cd8, pct_of_pbmc),
          file.path(REV_TAB, "R2M1e_02_cluster_frequency_both_denominators.csv"),
          row.names = FALSE)
write.csv(st, file.path(REV_TAB, "R2M1e_03_frequency_statistics.csv"),
          row.names = FALSE)
cat("R2-M1e: percentage-of-CD8 panel; min BH p_adj",
    sprintf("%.3f", min(st$p_adj)), "\n")

# ================================================================ R2-M2a ======
# Figure 2A plots clone size on the whole-PBMC embedding. The ask is a CD8-only
# UMAP.
#
# WHAT THE AVAILABLE DATA ALLOWS, STATED PLAINLY
#   The only CD8-only embedding in this study is the dedicated CD8 re-clustering,
#   which was run on HEI participants alone (18,843 cells, 11 participants). So a
#   CD8-only version of Figure 2A cannot carry the HEU and HUU facets the
#   published panel has; building one would mean re-integrating CD8 cells across
#   all three groups, which is a new analysis rather than a re-render. This panel
#   therefore shows the HEI CD8 compartment on its own embedding, faceted by
#   treatment stage, which is the comparison those cells support.
DED <- qs_read(QS_CD8)
cat("dedicated CD8 object:", ncol(DED), "cells,",
    length(unique(DED$orig.ident)), "samples\n")

cs_levels <- c("Single (0 < X <= 1)", "Small (1 < X <= 5)", "Medium (5 < X <= 20)",
               "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 500)")
cs_short  <- c("Single", "Small (2-5)", "Medium (6-20)", "Large (21-100)",
               "Hyperexpanded (>100)")
cs_cols   <- setNames(c("#FEE0B6", "#F1A340", "#E08214", "#B35806", "#7F3B08"),
                      cs_short)

emb <- as.data.frame(Embeddings(DED, reduction = "wnn.umap"))
names(emb) <- c("UMAP1", "UMAP2")
emb$cloneSize <- as.character(DED$cloneSize)
emb$stage <- recode(as.character(DED$Timepoint_Group),
                    PreART_Entry = "Pre-ART",
                    PostART_Suppressed = "Suppressed",
                    PostART_Unsuppressed = "Unsuppressed")
emb$stage <- factor(emb$stage, levels = c("Pre-ART", "Suppressed", "Unsuppressed"))
emb$cs <- ifelse(is.na(emb$cloneSize), NA_character_,
                 cs_short[match(emb$cloneSize, cs_levels)])
emb$cs <- factor(emb$cs, levels = cs_short)

n_by_stage <- emb %>% group_by(stage) %>%
  summarise(n = n(), n_tcr = sum(!is.na(cs)), .groups = "drop") %>%
  mutate(lab = paste0(stage, "\n", n_lab(n), " cells | ", n_lab(n_tcr), " with a TCR"))
emb <- left_join(emb, n_by_stage[, c("stage", "lab")], by = "stage")
emb$lab <- factor(emb$lab, levels = n_by_stage$lab)

p_cd8umap <- ggplot() +
  geom_point(data = emb[is.na(emb$cs), ], aes(UMAP1, UMAP2),
             colour = "grey86", size = 0.45) +
  geom_point(data = emb[!is.na(emb$cs), ], aes(UMAP1, UMAP2, colour = cs),
             size = 0.6) +
  facet_wrap(~ lab, nrow = 1) +
  scale_colour_manual(values = cs_cols, name = "Clonotype size",
                      drop = FALSE,
                      guide = guide_legend(override.aes = list(size = 5))) +
  coord_equal() +
  labs(title = "Clonotype expansion on the CD8-only embedding",
       subtitle = paste0("dedicated CD8 re-clustering, HEI participants only (",
                         n_lab(ncol(DED)), " cells, ",
                         length(unique(DED$orig.ident)), " samples, ",
                         length(unique(DED$orig.ident %>% substr(1, 5))),
                         " participants); grey = no TCR recovered\n",
                         "faceted by treatment stage, the comparison these cells ",
                         "support")) +
  theme_opt(19) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave(file.path(REV_FIG, "R2M2a_02_CD8only_clonotype_UMAP.png"), plot = p_cd8umap,
       width = 19, height = 8, dpi = 300, bg = "white")
write.csv(emb[, c("UMAP1", "UMAP2", "stage", "cs")],
          file.path(REV_TAB, "R2M2a_02_CD8only_clonotype_UMAP_source.csv"),
          row.names = FALSE)
cat("R2-M2a: CD8-only clonotype UMAP,", nrow(emb), "cells\n")

# ================================================================= R2-m2 ======
# The counts are in the legends. The ask is the count on the panel itself.
count_umap <- function(obj, reduction, group, title, sub) {
  e <- as.data.frame(Embeddings(obj, reduction = reduction))[, 1:2]
  names(e) <- c("UMAP1", "UMAP2")
  e$grp <- as.character(obj[[group]][, 1])
  ctr <- e %>% group_by(grp) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), n = n(),
              .groups = "drop")
  ggplot(e, aes(UMAP1, UMAP2)) +
    geom_point(aes(colour = grp), size = 0.4, alpha = 0.75, show.legend = FALSE) +
    # The count sits on the panel, beside each population. REPELLED, not placed
    # at the centroid: with 31 populations the boxes overlapped and partly hid
    # each other's digits, so a count could be misread (860 appeared as 8).
    ggrepel::geom_label_repel(
      data = ctr, aes(label = paste0(grp, "\nn = ", n_lab(n))),
      size = 4.3, fontface = "bold", label.size = 0.28,
      fill = alpha("white", 0.9), lineheight = 0.95,
      box.padding = 0.45, point.padding = 0.2, label.padding = 0.22,
      min.segment.length = 0, segment.size = 0.35, segment.colour = "grey35",
      max.overlaps = Inf, force = 6, seed = 1) +
    annotate("label", x = min(e$UMAP1), y = max(e$UMAP2), hjust = 0, vjust = 1,
             label = paste0("total n = ", n_lab(nrow(e))), size = 6.4,
             fontface = "bold", fill = alpha("white", 0.9), linewidth = 0.4) +
    coord_equal() +
    labs(title = title, subtitle = sub, x = NULL, y = NULL) +
    theme_opt(19) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank())
}

p_c1 <- count_umap(SUB, "wnn.umap", "Annotation",
  "Whole-PBMC UMAP with cell numbers on the panel",
  "every population labelled with its own count; total on the panel")
p_c2 <- count_umap(DED, "wnn.umap", "CD8_Annotation",
  "CD8 UMAP with cell numbers on the panel",
  "dedicated CD8 clustering, HEI participants")
ggsave(file.path(REV_FIG, "R2m2_02_UMAP_counts_on_panel.png"),
       plot = (p_c1 | p_c2), width = 30, height = 15, dpi = 300, bg = "white")
cat("R2-m2: counts annotated on both first-appearance UMAPs\n")

cat("\nAll panels written to", REV_FIG, "\n")
cat("NOTHING in Manuscript/ was written. Each panel is an option for the\n",
    "response document; the manuscript change is stated but not applied.\n", sep = "")
if (exists("sd_report")) sd_report()
