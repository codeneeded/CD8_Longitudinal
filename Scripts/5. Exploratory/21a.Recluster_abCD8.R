################################################################################
# Re-cluster the integrated CP003 + CP018 data on alpha-beta CD8 cells only
#
# The main analysis removes gamma-delta and MAIT cells BEFORE clustering CD8
# T cells (Scripts/2. Analysis/4.CD8_Subseting.R). The joint clusters in
# CP003_CP018_integrated.qs2 were computed with those cells present, and they
# shaped the clusters: one cluster was 55% gamma-delta, another 41% MAIT. So the
# cells are filtered first (00.Joint_abCD8_filter.R) and neighbours, clusters and
# UMAP are recomputed with the parameters used to build the integrated object
# (Scripts/CP018 - FACTseq/9.Integrate_CP003_CP018.R): MNN-integrated reduction,
# dims 1:30, k = 20, Louvain resolution 0.4; UMAP 30 neighbours, cosine, min.dist 0.3.
#
# Out: saved_R_data/CP003_CP018_abCD8.qs2
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr)})
R <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
source(file.path(R, "Scripts", "5. Exploratory", "00.Joint_abCD8_filter.R"))
o <- qs_read(file.path(R, "saved_R_data", "CP003_CP018_integrated.qs2"))
o <- subset(o, cells = joint_keep(o))
o <- FindNeighbors(o, reduction = "integrated", dims = 1:30, k.param = 20, verbose = FALSE)
o <- FindClusters(o, resolution = 0.4, algorithm = 1, cluster.name = "abcd8_clusters", verbose = FALSE)
o <- RunUMAP(o, reduction = "integrated", dims = 1:30, n.neighbors = 30, metric = "cosine",
             min.dist = 0.3, reduction.name = "umap.abcd8", verbose = FALSE)
qs_save(o, file.path(R, "saved_R_data", "CP003_CP018_abCD8.qs2"))
cat("alpha-beta CD8 cells:", ncol(o), "| clusters:", nlevels(o$abcd8_clusters), "\n")
