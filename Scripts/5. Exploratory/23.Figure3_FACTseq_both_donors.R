################################################################################
# Figure 3 and Supplementary Figure 3, FACT-Seq panels, rebuilt on the
# integrated CP003 + CP018 alpha-beta CD8+ T cells
#
# Panels B1/B2 are 21.Joint_UMAPs_annotated.R; this script rebuilds the rest
# with the construction and styling of Scripts/3. Manuscript/Figure_3.R.
#
# CELLS   alpha-beta CD8 T cells (00.Joint_abCD8_filter.R); cluster labels are
#         read from the source data of 21, so every panel uses the same cells
#         and the same labels as the UMAPs.
# HIV-RESPONSIVE  the published rule applied to both donors in 21 (clonotype in
#         >1 cell of its library, outside Naive/Bystander); read from 21. In the
#         unstimulated data a cell is matched when its clonotype is one of its
#         donor's HIV-responsive clonotypes.
# CLONOTYPE MATCHING  on CTstrict normalised for the Cell Ranger convention
#         (normalise_ct, CP018 config.R). CP003 was processed with Cell Ranger
#         9.0.1 and CP018 with 10.1.0, which write an empty D gene differently;
#         unnormalised strings match none of CP018's clonotypes.
# UNSTIMULATED DATA  CP003: tara_cp003_fig3.qs2 (1, 12, 24 months).
#                    CP018: tara_cp018_crossref.qs2 (1 month pre-ART, 25, 42 months).
#
# Out: Reviewer_Comments/Figures/L3_Fig3{C,D,E,F}_*_both_donors.png
#      Reviewer_Comments/Figures/L3_Fig3supp_{A..F}_*.png  (+ _source_data.csv)
#      Reviewer_Comments/Tables/L3_Fig3supp_module_score_statistics_both_donors.csv
# 'Fig3supp' is the supplementary figure of Figure 3 (repository folder
# Manuscript/Supplementary 3; Supplementary Figure 4 in Lesley's comments).
################################################################################
suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr); library(ggplot2)
  library(ggalluvial); library(patchwork); library(RColorBrewer)
})
R   <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
sd  <- file.path(R, "saved_R_data")
FIG <- file.path(R, "Reviewer_Comments", "Figures")
suppressMessages(source(file.path(R, "Scripts", "CP018 - FACTseq", "config.R")))   # normalise_ct()

save_panel <- function(p, fn, w, h, data, dpi = 300) {
  ggsave(file.path(FIG, fn), p, width = w, height = h, dpi = dpi, bg = "white")
  write.csv(data, file.path(FIG, sub("\\.png$", "_source_data.csv", fn)), row.names = FALSE)
  cat("  wrote", fn, "\n")
}

# ---- cells, labels, palette (Figure_3.R order and colours; Tc17 added) -------
src <- read.csv(file.path(FIG, "L3_jointUMAP_source_data.csv"), check.names = FALSE)
o <- qs_read(file.path(sd, "CP003_CP018_abCD8.qs2"))
stopifnot("21 and the object disagree on cells: re-run 21" = setequal(colnames(o), src$cell))
o <- o[, src$cell]
o$HIV_Specific_TCR <- ifelse(src$resp, "HIV-Specific TCR", "Other")   # one rule, set in 21
ANNO_ORDER <- c("Naive/Bystander", "Transitional Tem", "Tpex", "TEMRA/Effector",
                "Cycling Effector", "Tc17")
ANNO_COLS  <- c("Naive/Bystander" = "#2166AC", "Transitional Tem" = "#1A9850", "Tpex" = "#F4A261",
                "TEMRA/Effector" = "#E63946", "Cycling Effector" = "#9D4EDD", "Tc17" = "#17BECF")
o$Annotation <- factor(unname(src$joint_anno), levels = ANNO_ORDER)
stopifnot(!anyNA(o$Annotation))
o$ct <- unname(normalise_ct(o$CTstrict))
o$Visit <- factor(ifelse(o$Months == 2, "2 months", "Late"), levels = c("2 months", "Late"))
DONOR_LAB <- c(CP003 = "CP003 (unsuppressed)", CP018 = "CP018 (suppressed)")
o$Donor <- factor(unname(DONOR_LAB[o$participant]), levels = unname(DONOR_LAB))
m <- o@meta.data
resp <- m$HIV_Specific_TCR == "HIV-Specific TCR"
hiv_ct <- lapply(split(m$ct[resp], m$participant[resp]), function(x) unique(x[!is.na(x)]))
cat("HIV-responsive clonotypes: CP003", length(hiv_ct$CP003), "| CP018", length(hiv_ct$CP018), "\n")

hiv_cols   <- c("Other" = "grey70", "HIV-Specific TCR" = "#E63946")
hiv_labs   <- c("Other" = "Bystander", "HIV-Specific TCR" = "HIV-responsive")
phase_cols <- c("G1" = "#AEC6CF", "S" = "#FFD700", "G2M" = "#E63946")
tp_cols    <- c("2 months" = "#E63946", "Late" = "#2166AC")

# ---- unstimulated cross-reference, both donors -----------------------------
t3 <- qs_read(file.path(sd, "tara_cp003_fig3.qs2"))
u8 <- qs_read(file.path(sd, "tara_cp018_crossref.qs2"))
unstim_meta <- function(obj, donor, tp_fun) {
  x <- obj@meta.data
  x$ct <- normalise_ct(x$CTstrict)
  x$cell <- rownames(x)
  x$participant <- donor
  x$tp <- tp_fun(x)
  x$matched <- !is.na(x$ct) & x$ct %in% hiv_ct[[donor]]
  x
}
un3 <- unstim_meta(t3, "CP003", function(x) ifelse(x$Age == 1, "1 month\n(pre-ART)", paste0(x$Age, " months")))
un8 <- unstim_meta(u8, "CP018", function(x)
  unname(c(CP018_1m = "1 month\n(pre-ART)", CP018_25m = "25 months", CP018_42m = "42 months")[as.character(x$orig.ident)]))
stopifnot(!anyNA(un8$tp))
cat("matched unstimulated cells: CP003", sum(un3$matched), "| CP018", sum(un8$matched), "\n")
cat("CP003 matched cells (published 124):", sum(un3$matched), "\n")

# ---- 3C alluvial -------------------------------------------------------------
TP <- list(CP003 = c("1 month\n(pre-ART)", "2 months\n(HIV-stim)", "12 months", "24 months", "101 months\n(HIV-stim)"),
           CP018 = c("1 month\n(pre-ART)", "2 months\n(HIV-stim)", "25 months", "42 months", "90 months\n(HIV-stim)"))
stim_tp <- function(p, mo) paste0(mo, " months\n(HIV-stim)")
stopifnot(identical(names(hiv_ct), c("CP003", "CP018")))
allu <- bind_rows(
  # A clonotype HIV-responsive at either stimulation timepoint is tracked at EVERY
  # timepoint where any of its cells is detected - stimulated or unstimulated -
  # so a clonotype seen once at the other stimulation timepoint is not dropped.
  m[!is.na(m$ct) & m$ct %in% hiv_ct[[1]] & m$participant == names(hiv_ct)[1] |
    !is.na(m$ct) & m$ct %in% hiv_ct[[2]] & m$participant == names(hiv_ct)[2], ] %>%
    transmute(participant, CTstrict = ct, tp = stim_tp(participant, Months)),
  un3[un3$matched, c("participant", "ct", "tp")] %>% rename(CTstrict = ct),
  un8[un8$matched, c("participant", "ct", "tp")] %>% rename(CTstrict = ct)) %>%
  count(participant, CTstrict, tp, name = "n_cells")
multi <- allu %>% group_by(participant, CTstrict) %>% filter(n_distinct(tp) >= 2) %>% ungroup()
# geom_flow only joins ADJACENT axes, so a clonotype missing at an intermediate
# timepoint lost its ribbon and its later detection appeared unconnected (e.g. a
# CP018 clonotype at 1m and 2m seen again at 90m). As in 10.Figure3_Mirror_CP018.R,
# the clonotype x timepoint grid is completed and drawn with geom_alluvium, one
# continuous band per clonotype; where it is not detected the band runs at a thin
# floor width and no bar is drawn. Strata and ribbons stack from the same data.
alluvial_one <- function(p) {
  d <- multi %>% filter(participant == p)
  ord <- d %>% group_by(CTstrict) %>% summarise(t = sum(n_cells), .groups = "drop") %>% arrange(desc(t))
  pal <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(max(nrow(ord), 3))[seq_len(nrow(ord))], ord$CTstrict)
  floor_h <- 0.008 * max(tapply(d$n_cells, d$tp, sum))
  full <- expand.grid(CTstrict = ord$CTstrict, tp = TP[[p]], stringsAsFactors = FALSE) %>%
    left_join(d %>% select(CTstrict, tp, n_cells), by = c("CTstrict", "tp")) %>%
    mutate(n_cells = ifelse(is.na(n_cells), 0, n_cells), detected = n_cells > 0,
           y = pmax(n_cells, floor_h), tp = factor(tp, levels = TP[[p]]),
           CTstrict = factor(CTstrict, levels = ord$CTstrict),
           bar_alpha = ifelse(detected, 1, 0), bar_line = ifelse(detected, 0.3, 0))
  # only bridge between a clonotype's first and last detection
  span <- full %>% group_by(CTstrict) %>%
    summarise(lo = min(as.integer(tp)[detected]), hi = max(as.integer(tp)[detected]), .groups = "drop")
  full <- full %>% left_join(span, by = "CTstrict") %>%
    mutate(y = ifelse(as.integer(tp) < lo | as.integer(tp) > hi, 0, y))
  ggplot(full, aes(x = tp, y = y, alluvium = CTstrict, stratum = CTstrict, fill = CTstrict)) +
    geom_alluvium(alpha = 0.5, width = 0.3, knot.pos = 0.28) +
    geom_stratum(aes(alpha = bar_alpha, linewidth = bar_line), width = 0.3, colour = "grey30") +
    scale_alpha_identity() + scale_linewidth_identity() +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(expand = c(0.08, 0.08), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    theme_classic(base_size = 30) +
    theme(axis.text.x = element_text(size = 22, lineheight = 0.9, face = "bold"),
          axis.text.y = element_text(size = 22), axis.title = element_text(size = 26, face = "bold"),
          plot.title = element_text(size = 34, face = "bold"),
          plot.subtitle = element_text(size = 24, color = "grey25", face = "bold")) +
    labs(x = NULL, y = "Number of cells", title = DONOR_LAB[[p]],
         subtitle = paste0(nrow(ord), " clonotypes detected at \u22652 timepoints"))
}
pC <- (alluvial_one("CP003") | alluvial_one("CP018")) +
  plot_annotation(title = "HIV-responsive clonotype persistence",
                  theme = theme(plot.title = element_text(size = 40, face = "bold")))
save_panel(pC, "L3_Fig3C_alluvial_both_donors.png", 30, 12, multi)

# ---- 3D functional profile, matched cells in unstimulated data --------------
func_genes <- c("GZMB","GZMA","PRF1","GNLY","NKG7","IFNG","CCL3","CCL4","CCL5",
                "TBX21","EOMES","PRDM1","CX3CR1","KLRG1","FAS","PDCD1","LAG3","HAVCR2","TOX")
category <- c(rep("Cytotoxicity",5), rep("Cytokines/Chemokines",4), rep("Effector TFs",3),
              rep("Effector Markers",3), rep("Exhaustion",4))
cat_cols <- c("Cytotoxicity"="#E15759","Cytokines/Chemokines"="#FF9DA7","Effector TFs"="#F28E2B",
              "Effector Markers"="#EDC948","Exhaustion"="#9C755F")
prof <- bind_rows(lapply(list(list(t3, un3, "CP003"), list(u8, un8, "CP018")), function(a) {
  obj <- a[[1]]; DefaultAssay(obj) <- "RNA"
  if (length(Layers(obj[["RNA"]])) > 1) obj <- JoinLayers(obj, assay = "RNA")
  cells <- a[[2]]$cell[a[[2]]$matched]
  X <- GetAssayData(obj, assay = "RNA", layer = "data")[func_genes, cells, drop = FALSE]
  data.frame(participant = a[[3]], gene = func_genes, category = category, n_cells = length(cells),
             pct_expressing = round(100 * Matrix::rowMeans(X > 0), 1))
}))
prof$Donor <- factor(paste0(DONOR_LAB[prof$participant], "\nn = ", prof$n_cells, " cells"))
prof$gene <- factor(prof$gene, levels = rev(func_genes))
prof$category <- factor(prof$category, levels = names(cat_cols))
pD <- ggplot(prof, aes(gene, pct_expressing, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0(pct_expressing, "%")), hjust = -0.08, size = 6, fontface = "bold") +
  facet_wrap(~Donor, nrow = 1) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_flip() + theme_classic(base_size = 26) +
  theme(axis.text.y = element_text(size = 20, face = "italic"), axis.text.x = element_text(size = 20),
        strip.background = element_blank(), strip.text = element_text(size = 22, face = "bold"),
        legend.position = "bottom", legend.text = element_text(size = 20),
        plot.title = element_text(size = 32, face = "bold")) +
  labs(x = NULL, y = "% cells expressing",
       title = "Functional profile of HIV-responsive clonotypes in unstimulated samples") +
  guides(fill = guide_legend(nrow = 1))
save_panel(pD, "L3_Fig3D_functional_profile_both_donors.png", 22, 13, prof)

# ---- 3E persistence ---------------------------------------------------------
pers <- allu %>% group_by(participant, CTstrict) %>% summarise(n_timepoints = n_distinct(tp), .groups = "drop") %>%
  count(participant, n_timepoints, name = "n") %>%
  complete(participant, n_timepoints = 1:max(n_timepoints), fill = list(n = 0))
ntot <- pers %>% group_by(participant) %>% summarise(t = sum(n))
pers$Donor <- factor(paste0(DONOR_LAB[pers$participant], "\n", ntot$t[match(pers$participant, ntot$participant)],
                            " HIV-responsive clonotypes"))
pE <- ggplot(pers, aes(factor(n_timepoints), n)) +
  geom_col(fill = "#E63946", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 10, fontface = "bold") +
  facet_wrap(~Donor, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 28) +
  theme(axis.text = element_text(size = 24, face = "bold"), strip.background = element_blank(),
        strip.text = element_text(size = 24, face = "bold"), plot.title = element_text(size = 34, face = "bold")) +
  labs(x = "Timepoints detected", y = "N clonotypes", title = "Clonotype persistence")
save_panel(pE, "L3_Fig3E_persistence_both_donors.png", 18, 10, pers)

# ---- 3F cluster distribution in unstimulated data ----------------------------
CL6 <- c("Naive 1 CD8","Naive 2 CD8","Naive Intermediate CD8","Tscm CD8","TEMRA/CTL CD8","Tex CD8")
CL6_COLS <- c("Naive 1 CD8"="#2166AC","Naive 2 CD8"="#67A9CF","Naive Intermediate CD8"="#D6604D",
              "Tscm CD8"="#1A9850","TEMRA/CTL CD8"="#B2182B","Tex CD8"="#762A83")   # Figure_1.R
f3 <- bind_rows(un3[un3$matched, c("participant","Annotation","ct")], un8[un8$matched, c("participant","Annotation","ct")]) %>%
  mutate(Annotation = as.character(Annotation))
f3_excl <- f3 %>% filter(!Annotation %in% CL6) %>% count(participant, Annotation, name = "n_excluded_nonCD8")
f3d <- f3 %>% filter(Annotation %in% CL6) %>% count(participant, Annotation, name = "n") %>%
  group_by(participant) %>% mutate(pct = round(100 * n / sum(n), 1), n_cells = sum(n),
                                   n_clono = n_distinct(f3$ct[f3$participant == first(participant) & f3$Annotation %in% CL6])) %>% ungroup()
cat("CP003 TEMRA/CTL matched cells (published 117):", f3d$n[f3d$participant == "CP003" & f3d$Annotation == "TEMRA/CTL CD8"], "\n")
f3d$Donor <- factor(paste0(DONOR_LAB[f3d$participant], "\n", f3d$n_cells, " CD8+ cells from ", f3d$n_clono, " clonotypes"))
f3d$Annotation <- factor(f3d$Annotation, levels = rev(CL6))
pF <- ggplot(f3d, aes(Annotation, n, fill = Annotation)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n, " (", pct, "%)")), hjust = -0.08, size = 8, fontface = "bold") +
  facet_wrap(~Donor, nrow = 1, scales = "free_x") +
  scale_x_discrete(drop = FALSE) + scale_fill_manual(values = CL6_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.55))) +
  coord_flip() + theme_classic(base_size = 26) +
  theme(axis.text.y = element_text(size = 22, face = "bold"), strip.background = element_blank(),
        strip.text = element_text(size = 22, face = "bold"), plot.title = element_text(size = 32, face = "bold")) +
  labs(x = NULL, y = "N cells", title = "Cluster distribution in unstimulated samples")
save_panel(pF, "L3_Fig3F_TARA_clusters_both_donors.png", 20, 8, bind_rows(f3d, f3_excl))

# ---- S3A cell cycle ---------------------------------------------------------
cc <- m %>% filter(!is.na(Phase)) %>% count(Donor, Annotation, Phase) %>%
  group_by(Donor, Annotation) %>% mutate(proportion = n / sum(n)) %>% ungroup()
pS3A <- ggplot(cc, aes(Annotation, proportion, fill = Phase)) +
  geom_col(width = 0.7) + facet_wrap(~Donor, nrow = 1) +
  scale_fill_manual(values = phase_cols) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
        strip.background = element_blank(), strip.text = element_text(size = 20, face = "bold"),
        legend.position = "top", plot.title = element_text(size = 28, face = "bold")) +
  labs(x = NULL, y = "Proportion", fill = "Cell cycle", title = "Cell cycle by cluster annotation")
save_panel(pS3A, "L3_Fig3supp_A_cell_cycle_both_donors.png", 16, 8, cc)

# ---- module scores (S3B, S3F) -----------------------------------------------
DefaultAssay(o) <- "RNA"; if (length(Layers(o[["RNA"]])) > 1) o <- JoinLayers(o, assay = "RNA")
o <- NormalizeData(o, verbose = FALSE)
modules <- list(Cytotoxicity = c("GZMB","GZMA","PRF1","GNLY","NKG7","FASLG"),
                Exhaustion   = c("TOX","PDCD1","HAVCR2","LAG3","CTLA4","TIGIT","CD244","ENTPD1"),
                Stemness     = c("TCF7","LEF1","SELL","CCR7","IL7R","BACH2","BCL2","S1PR1"),
                `IFN Memory` = c("IFIT1","IFIT3","ISG15","MX1","OAS1","STAT1"))
for (k in names(modules))
  o <- AddModuleScore(o, features = list(intersect(modules[[k]], rownames(o))),
                      name = paste0(gsub(" ", "_", k), "_mod"), seed = 42)
mc <- paste0(gsub(" ", "_", names(modules)), "_mod1")
long <- o@meta.data %>% select(Donor, Visit, HIV_Specific_TCR, all_of(mc)) %>%
  pivot_longer(all_of(mc), names_to = "Module", values_to = "Score") %>%
  mutate(Module = factor(gsub("_", " ", sub("_mod1$", "", Module)), levels = names(modules)),
         HIV_Specific_TCR = factor(HIV_Specific_TCR, levels = c("Other", "HIV-Specific TCR")))
vio <- function(d, fillvar, cols, title, labs_ = waiver()) {
  ggplot(d, aes(Module, Score, fill = .data[[fillvar]])) +
    geom_violin(scale = "width", width = 0.8, alpha = 0.8, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.6, position = position_dodge(0.8)) +
    facet_wrap(~Donor, nrow = 1) + scale_fill_manual(values = cols, labels = labs_, name = NULL) +
    theme_classic(base_size = 24) +
    theme(axis.text.x = element_text(size = 18, face = "bold"), strip.background = element_blank(),
          strip.text = element_text(size = 20, face = "bold"), legend.position = "top",
          plot.title = element_text(size = 28, face = "bold")) +
    labs(x = NULL, y = "Module score", title = title)
}
fmt_p <- function(x) ifelse(x < 2.2e-16, "p < 2.2e-16", paste0("p = ", formatC(x, format = "g", digits = 2)))
mod_stats <- function(d, grp, comparison) {
  d %>% group_by(Donor, Module) %>%
    summarise(n1 = sum(.data[[grp]] == levels(.data[[grp]])[1]), n2 = sum(.data[[grp]] == levels(.data[[grp]])[2]),
              p = tryCatch(wilcox.test(Score ~ .data[[grp]])$p.value, error = function(e) NA_real_),
              y = max(Score) + 0.25, .groups = "drop") %>%
    mutate(comparison = comparison, test = "Wilcoxon rank-sum, two-sided", label = fmt_p(p))
}
add_p <- function(g, st) g + geom_text(data = st, aes(Module, y, label = label), inherit.aes = FALSE,
                                       size = 5.5, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)))
st_b <- mod_stats(long, "HIV_Specific_TCR", "HIV-responsive vs bystander")
save_panel(add_p(vio(long, "HIV_Specific_TCR", hiv_cols, "HIV-responsive vs bystander (stimulation data)", hiv_labs), st_b),
           "L3_Fig3supp_B_module_scores_HIVvsOther_both_donors.png", 18, 9,
           long %>% group_by(Donor, Module, HIV_Specific_TCR) %>%
             summarise(n = n(), median = round(median(Score), 3), .groups = "drop"))
lr <- long %>% filter(HIV_Specific_TCR == "HIV-Specific TCR")
nv <- o@meta.data %>% filter(HIV_Specific_TCR == "HIV-Specific TCR") %>% count(Donor, Visit)
st_f <- mod_stats(lr, "Visit", "HIV-responsive cells, 2 months vs late")
write.csv(bind_rows(st_b, st_f) %>% select(comparison, Donor, Module, n1, n2, test, p),
          file.path(R, "Reviewer_Comments", "Tables", "L3_Fig3supp_module_score_statistics_both_donors.csv"), row.names = FALSE)
print(as.data.frame(bind_rows(st_b, st_f) %>% transmute(comparison, Donor = sub(" .*", "", Donor), Module, n1, n2, p = signif(p, 2))), row.names = FALSE)
save_panel(add_p(vio(lr, "Visit", tp_cols, "HIV-responsive cells: 2 months vs late (CP003 101 months, CP018 90 months)"), st_f) +
             labs(subtitle = paste(sprintf("%s %s: n = %d", sub(" .*", "", nv$Donor), nv$Visit, nv$n), collapse = "   ")) +
             theme(plot.subtitle = element_text(size = 18, colour = "#E63946", face = "bold")),
           "L3_Fig3supp_F_module_scores_2m_vs_late_both_donors.png", 18, 9,
           lr %>% group_by(Donor, Module, Visit) %>% summarise(n = n(), median = round(median(Score), 3), .groups = "drop"))

# ---- S3C dot plot: published gene set first, Tc17 genes appended -------------
dot_genes <- c("TCF7","LEF1","SELL","IL7R","BACH2","GZMK","EOMES","GZMB","GNLY","PRF1","NKG7","CX3CR1","KLRG1",
               "TOX","PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","MKI67","TOP2A","CCL3","CCL4","CCL5","IFNG","CD69","HLA-DRA",
               "KLRB1","CCR6","RORC","IL23R")
dot_genes <- intersect(dot_genes, rownames(o))
pS3C <- DotPlot(o, features = dot_genes, group.by = "Annotation", cols = c("grey90", "#E15759"), dot.scale = 10) +
  coord_flip() + theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(size = 18, angle = 40, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 16, face = "italic"), plot.title = element_text(size = 24, face = "bold")) +
  guides(colour = guide_colourbar(title = "Mean expression\n(z-score across\nthe groups shown)"),
         size = guide_legend(title = "Percent expressed")) +
  labs(x = NULL, y = NULL, title = "Key marker expression by annotation",
       subtitle = "Submitted gene set; KLRB1, CCR6, RORC and IL23R added") +
  theme(plot.subtitle = element_text(size = 17))
save_panel(pS3C, "L3_Fig3supp_C_dot_plot_both_donors.png", 14, 14, pS3C$data)

# ---- S3D HIV-responsive % per cluster, S3E timepoint composition -------------
hp <- m %>% group_by(Donor, Annotation) %>%
  summarise(n_cells = n(), n_resp = sum(HIV_Specific_TCR == "HIV-Specific TCR"), .groups = "drop") %>%
  mutate(pct_hiv = round(100 * n_resp / n_cells, 1))
pS3D <- ggplot(hp, aes(Annotation, pct_hiv, fill = Annotation)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct_hiv, "%")), vjust = -0.5, size = 7, fontface = "bold") +
  facet_wrap(~Donor, nrow = 1) + scale_fill_manual(values = ANNO_COLS) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
        strip.background = element_blank(), strip.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 28, face = "bold")) +
  labs(x = NULL, y = "% of cells in the cluster\nclassified HIV-responsive", title = "HIV-responsive cell enrichment by cluster")
save_panel(pS3D, "L3_Fig3supp_D_HIV_responsive_by_cluster_both_donors.png", 16, 8, hp)

tpc <- m %>% count(Donor, Annotation, Visit) %>% group_by(Donor, Annotation) %>%
  mutate(pct = n / sum(n)) %>% ungroup()
pS3E <- ggplot(tpc, aes(Annotation, pct, fill = Visit)) +
  geom_col(width = 0.7) + facet_wrap(~Donor, nrow = 1) +
  scale_fill_manual(values = tp_cols, labels = c("2 months", "Late (CP003 101m, CP018 90m)")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
        strip.background = element_blank(), strip.text = element_text(size = 20, face = "bold"),
        legend.position = "top", plot.title = element_text(size = 28, face = "bold")) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Timepoint composition by cluster")
save_panel(pS3E, "L3_Fig3supp_E_timepoint_composition_both_donors.png", 16, 8, tpc)

# ---- numbers the text needs --------------------------------------------------
cat("\n=== per donor ===\n")
print(as.data.frame(m %>% group_by(participant) %>%
  summarise(cells = n(), resp_cells = sum(HIV_Specific_TCR == "HIV-Specific TCR"),
            resp_2m = sum(HIV_Specific_TCR == "HIV-Specific TCR" & Months == 2),
            resp_late = sum(HIV_Specific_TCR == "HIV-Specific TCR" & Months != 2), .groups = "drop")), row.names = FALSE)
cat("\n=== per cluster: cells, % HIV-responsive, % from 2 months ===\n")
print(as.data.frame(m %>% group_by(Annotation) %>%
  summarise(cells = n(), pct_hiv = round(100*mean(HIV_Specific_TCR == "HIV-Specific TCR"),1),
            pct_2m = round(100*mean(Months == 2),1), .groups = "drop")), row.names = FALSE)
cat("\n=== persistence ===\n"); print(as.data.frame(pers[, c("participant","n_timepoints","n")]), row.names = FALSE)
cat("\n=== 3F ===\n"); print(as.data.frame(f3d[, c("participant","Annotation","n","pct")]), row.names = FALSE)
print(as.data.frame(f3_excl), row.names = FALSE)
