################################################################################
# Revision step 3 — score the pilot: are there really HIV+ cells?
#
# ANSWERS: Reviewer 1, major point 3.
#
# A raw HIV count per cell is NOT an infected-cell call. Ambient viral RNA from
# lysed cells is present throughout the suspension and gets into every droplet,
# so some barcodes carry viral reads with no infection. This script runs the
# three controls that separate signal from that background, and writes a verdict
# that says explicitly whether the full 26-library realignment is worth doing.
#
#   CONTROL 1  empty droplets      ambient background rate
#   CONTROL 2  lineage restriction HIV should be in CD4 T cells, not B/Mono
#   CONTROL 3  host transcriptome  adding viral contigs must not perturb GRCh38
#
# INPUT   cellranger_out/pilot_hiv_<sample>/outs/{raw,filtered}_feature_bc_matrix
# OUTPUT  Reviewer_Comments/Tables/Pilot_HIV_*.csv
#         Reviewer_Comments/Figures/RevFig_HIV_pilot_*.png
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

PILOT_SAMPLE <- Sys.getenv("PILOT_SAMPLE", "CP006_1m")
CR_OUT <- Sys.getenv("OUT_PARENT", path.expand("~/Desktop/Projects/cellranger_out"))
PILOT   <- file.path(CR_OUT, paste0("pilot_hiv_", PILOT_SAMPLE), "outs")
HIV_FEATURES <- c("HIV1_HXB2_subB", "HIV1_subC")

require_obj(file.path(PILOT, "filtered_feature_bc_matrix"))

cat("\n================ PILOT:", PILOT_SAMPLE, "================\n")

################################################################################
# Load both matrices. The RAW one is essential - it holds the empty droplets
# that measure ambient background.
################################################################################
filt <- Read10X(file.path(PILOT, "filtered_feature_bc_matrix"))
raw  <- Read10X(file.path(PILOT, "raw_feature_bc_matrix"))
if (is.list(filt)) filt <- filt[["Gene Expression"]]
if (is.list(raw))  raw  <- raw[["Gene Expression"]]

present <- intersect(HIV_FEATURES, rownames(filt))
if (!length(present)) {
  stop("No HIV features in the matrix. The reference did not include the viral\n",
       "contigs, or mkref dropped them. Re-run 1.Build_HIV_Reference.sh and\n",
       "check its VERIFYING section.", call. = FALSE)
}
cat("HIV features found:", paste(present, collapse = ", "), "\n")
cat("cells (filtered):", ncol(filt), " barcodes (raw):", ncol(raw), "\n")

hiv_filt <- Matrix::colSums(filt[present, , drop = FALSE])
hiv_raw  <- Matrix::colSums(raw[present, , drop = FALSE])
umi_raw  <- Matrix::colSums(raw)

################################################################################
# CONTROL 1 — ambient background from empty droplets
################################################################################
cat("\n--- CONTROL 1: ambient viral RNA in empty droplets ---\n")

is_cell  <- colnames(raw) %in% colnames(filt)
# Empty droplets: low total UMI and not called as a cell. The 100-UMI ceiling
# keeps genuine low-quality cells out of the "empty" set.
is_empty <- !is_cell & umi_raw > 0 & umi_raw < 100

amb <- data.frame(
  droplet_class = c("called cells", "empty droplets"),
  n_droplets    = c(sum(is_cell), sum(is_empty)),
  n_HIV_pos     = c(sum(hiv_raw[is_cell] > 0), sum(hiv_raw[is_empty] > 0)),
  total_HIV_UMI = c(sum(hiv_raw[is_cell]), sum(hiv_raw[is_empty])),
  median_UMI    = c(median(umi_raw[is_cell]), median(umi_raw[is_empty]))
)
amb$pct_HIV_pos    <- round(100 * amb$n_HIV_pos / amb$n_droplets, 3)
# HIV UMI per 10,000 total UMI - the rate that is comparable between the two
# classes despite very different sequencing depth.
amb$HIV_per_10k_UMI <- round(1e4 * amb$total_HIV_UMI /
                             c(sum(umi_raw[is_cell]), sum(umi_raw[is_empty])), 4)
print(amb, row.names = FALSE)
save_tab("Pilot_HIV_control1_ambient", amb)

ambient_rate <- amb$HIV_per_10k_UMI[2]
cell_rate    <- amb$HIV_per_10k_UMI[1]
enrichment   <- if (ambient_rate > 0) round(cell_rate / ambient_rate, 2) else NA
cat("\n  depth-normalised enrichment in cells vs empty droplets:",
    ifelse(is.na(enrichment), "empty droplets had zero HIV UMI", enrichment), "\n")
cat("  Interpretation: a ratio near 1 means the viral reads in cells are\n")
cat("  ambient contamination. Genuine infection gives a ratio well above 1.\n")

################################################################################
# CONTROL 2 — is HIV where it biologically should be?
################################################################################
cat("\n--- CONTROL 2: lineage of HIV+ cells ---\n")

seu <- CreateSeuratObject(filt, project = PILOT_SAMPLE, min.cells = 0, min.features = 200)
seu$HIV_UMI <- hiv_filt[colnames(seu)]
seu$HIV_pos <- seu$HIV_UMI > 0

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

# Score the canonical lineages by marker expression rather than importing an
# annotation, so this pilot stands alone.
lin_markers <- list(
  CD4_T = c("CD3E","CD3D","CD4","IL7R","CCR7"),
  CD8_T = c("CD3E","CD3D","CD8A","CD8B","GZMK"),
  B     = c("MS4A1","CD79A","CD79B","IGHM"),
  Mono  = c("LYZ","CD14","S100A8","FCN1"),
  NK    = c("NKG7","GNLY","KLRD1","NCAM1"))
for (nm in names(lin_markers)) {
  g <- intersect(lin_markers[[nm]], rownames(seu))
  seu <- AddModuleScore(seu, features = list(g), name = paste0("lin_", nm), verbose = FALSE)
}
score_cols <- grep("^lin_", colnames(seu@meta.data), value = TRUE)
cl_scores <- seu@meta.data %>%
  group_by(cluster = seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop")
lin_call <- setNames(
  sub("^lin_", "", sub("1$", "", score_cols[apply(cl_scores[score_cols], 1, which.max)])),
  cl_scores$cluster)
seu$lineage <- lin_call[as.character(seu$seurat_clusters)]

lin_tab <- seu@meta.data %>%
  group_by(lineage) %>%
  summarise(n_cells = n(), n_HIV_pos = sum(HIV_pos),
            total_HIV_UMI = sum(HIV_UMI), .groups = "drop") %>%
  mutate(pct_HIV_pos = round(100 * n_HIV_pos / n_cells, 3)) %>%
  arrange(desc(n_HIV_pos))
print(as.data.frame(lin_tab), row.names = FALSE)
save_tab("Pilot_HIV_control2_lineage", lin_tab)

# Fisher test: are HIV+ cells enriched in CD4 T cells specifically?
cd4_test <- NA_real_
if (any(lin_tab$lineage == "CD4_T") && sum(lin_tab$n_HIV_pos) >= 2) {
  m <- matrix(c(
    lin_tab$n_HIV_pos[lin_tab$lineage == "CD4_T"],
    sum(lin_tab$n_HIV_pos) - lin_tab$n_HIV_pos[lin_tab$lineage == "CD4_T"],
    lin_tab$n_cells[lin_tab$lineage == "CD4_T"] - lin_tab$n_HIV_pos[lin_tab$lineage == "CD4_T"],
    sum(lin_tab$n_cells) - lin_tab$n_cells[lin_tab$lineage == "CD4_T"] -
      (sum(lin_tab$n_HIV_pos) - lin_tab$n_HIV_pos[lin_tab$lineage == "CD4_T"])),
    nrow = 2)
  ft <- fisher.test(m)
  cd4_test <- ft$p.value
  cat(sprintf("\n  CD4-T enrichment of HIV+ cells: Fisher p = %.3g, OR = %.2f\n",
              ft$p.value, ft$estimate))
} else {
  cat("\n  Too few HIV+ cells for a lineage enrichment test.\n")
}

################################################################################
# CONTROL 3 — did adding viral contigs perturb the human alignment?
################################################################################
cat("\n--- CONTROL 3: host transcriptome unchanged ---\n")

host_qc <- data.frame(
  metric = c("cells called", "median genes/cell", "median UMI/cell",
             "total human genes detected"),
  pilot_GRCh38_HIV = c(ncol(seu), median(seu$nFeature_RNA), median(seu$nCount_RNA),
                       sum(Matrix::rowSums(filt) > 0)))

# Compare against the published alignment for the same sample if available.
if (file.exists(QS_ALL)) {
  pub <- qs_read(QS_ALL)
  pm <- pub@meta.data[pub@meta.data$orig.ident == PILOT_SAMPLE, ]
  if (nrow(pm)) {
    host_qc$published_GRCh38only <- c(nrow(pm), median(pm$nFeature_RNA),
                                      median(pm$nCount_RNA), NA)
    host_qc$pct_change <- round(100 * (host_qc$pilot_GRCh38_HIV -
                                host_qc$published_GRCh38only) /
                                host_qc$published_GRCh38only, 1)
  }
  rm(pub); gc(verbose = FALSE)
}
print(host_qc, row.names = FALSE)
save_tab("Pilot_HIV_control3_host_qc", host_qc)
cat("\n  Note: the published object is post-QC and post-doublet-removal, so a\n")
cat("  lower published cell count is expected. Genes and UMI per cell are the\n")
cat("  comparable metrics; those should be within a few percent.\n")

################################################################################
# Figures
################################################################################
p1 <- ggplot(seu@meta.data, aes(x = lineage, fill = HIV_pos)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(`FALSE` = "grey85", `TRUE` = "#D55E00"),
                    labels = c("HIV-", "HIV+"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  labs(title = paste0("HIV-mapping reads by lineage - ", PILOT_SAMPLE),
       subtitle = paste0(sum(seu$HIV_pos), " of ", ncol(seu),
                         " cells carry >=1 HIV UMI; ", n_cells_label(seu)),
       x = NULL, y = "Fraction of cells",
       caption = paste0("Ambient control: ", cell_rate, " vs ", ambient_rate,
                        " HIV UMI per 10,000 total UMI in cells vs empty droplets.")) +
  theme_rev(11) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_fig("RevFig_HIV_pilot_lineage", p1, width = 8, height = 5.5)
save_src("RevFig_HIV_pilot_lineage", lin_tab)

um <- as.data.frame(Embeddings(seu, "umap"))
names(um) <- c("UMAP_1","UMAP_2")
um$HIV_pos <- seu$HIV_pos; um$lineage <- seu$lineage; um$HIV_UMI <- seu$HIV_UMI
p2 <- ggplot(um, aes(UMAP_1, UMAP_2)) +
  geom_point(data = subset(um, !HIV_pos), colour = "grey88", size = 0.35) +
  geom_point(data = subset(um, HIV_pos), aes(size = HIV_UMI),
             colour = "#D55E00", alpha = 0.9) +
  scale_size_continuous(range = c(1.8, 5), name = "HIV UMI") +
  labs(title = paste0("HIV+ cells on the UMAP - ", PILOT_SAMPLE),
       subtitle = paste0(n_cells_label(seu), "; grey = no HIV reads"),
       caption = "Unannotated pilot clustering; lineages called by marker module score.") +
  theme_rev(11)
save_fig("RevFig_HIV_pilot_umap", p2, width = 7.5, height = 6.5)
save_src("RevFig_HIV_pilot_umap", um)

################################################################################
# VERDICT — say plainly whether the full realignment is worth doing
################################################################################
n_pos     <- sum(seu$HIV_pos)
n_cd4_pos <- if (any(lin_tab$lineage == "CD4_T"))
               lin_tab$n_HIV_pos[lin_tab$lineage == "CD4_T"] else 0
amb_clean <- is.na(enrichment) || enrichment >= 3

verdict <- if (n_pos >= 5 && amb_clean && n_cd4_pos >= 1) {
  "PROCEED - HIV is detectable above ambient background"
} else if (n_pos >= 1 && amb_clean) {
  "BORDERLINE - detectable but too few cells to phenotype from one sample"
} else if (n_pos >= 1) {
  "AMBIENT - viral reads present but not enriched over empty droplets"
} else {
  "NEGATIVE - no HIV-mapping reads in called cells"
}

res <- data.frame(
  sample = PILOT_SAMPLE, viral_load = 1e7,
  cells = ncol(seu), HIV_pos_cells = n_pos,
  pct_HIV_pos = round(100 * n_pos / ncol(seu), 3),
  HIV_pos_CD4_T = n_cd4_pos,
  ambient_enrichment = enrichment,
  cd4_enrichment_p = cd4_test,
  verdict = verdict)
save_tab("Pilot_HIV_VERDICT", res)

cat("\n================ VERDICT ================\n")
print(t(res))
cat("\n", verdict, "\n\n", sep = "")

if (startsWith(verdict, "NEGATIVE") || startsWith(verdict, "AMBIENT")) {
  cat("Do NOT realign the remaining 25 libraries on this basis.\n")
  cat("This sample has the highest viral load in the cohort (10^7 copies/mL,\n")
  cat("pre-ART), so a negative here is close to decisive for the rest.\n\n")
  cat("For the response letter, the honest answer is a bounded negative:\n")
  cat("  - standard 10x poly-dT chemistry under-captures HIV RNA, which is\n")
  cat("    frequently not polyadenylated; published work needing reliable\n")
  cat("    detection added HIV-specific capture primers at the RT step\n")
  cat("  - our libraries were not made that way, so this is a chemistry\n")
  cat("    limitation, not evidence that infected cells are absent\n")
  cat("  - state the detection ceiling and the cells examined, and offer the\n")
  cat("    targeted-capture experiment as future work\n")
} else if (startsWith(verdict, "BORDERLINE")) {
  cat("Detectable but sparse. Realigning the 15 viremic samples is justified;\n")
  cat("the 6 suppressed ones are unlikely to add cells. Expect a descriptive\n")
  cat("result - the pooled yield will not support differential expression.\n")
} else {
  cat("Proceed to the 15 viremic samples. Keep the ambient control per sample;\n")
  cat("it is the analysis the reviewer will scrutinise, not the raw count.\n")
}

cat("\nTables:  ", REV_TAB, "\n")
cat("Figures: ", REV_FIG, "\n")
