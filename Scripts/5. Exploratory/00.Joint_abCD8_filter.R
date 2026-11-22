################################################################################
# Single definition of the alpha-beta CD8 cells used by every two-donor panel
# (CP003 + CP018, CP003_CP018_integrated.qs2). Sourced by scripts 14, 17, 20-23.
#
# Removed, in this order of precedence:
#   1. joint cluster 10 - heat-shock genes (HSPA1A/HSPA1B), a handling-stress state
#   2. joint cluster 8  - CD4-expressing activated cells and low-quality cells
#   3. MAIT             - TRAV1-2 joined to TRAJ33, TRAJ12 or TRAJ20. This is the
#                         receptor the Va7.2 antibody detects in the main analysis;
#                         FACT-seq has no antibody data, so the TCR is used.
#   4. gamma-delta     - TRDC detected at >= 2 UMI in a cell with no TCR alpha
#                         chain recovered. Cells lacking a recovered alpha chain
#                         but without TRDC are kept: most carry a beta chain or
#                         TRAC transcripts, i.e. alpha-beta cells whose alpha chain
#                         was not captured.
#
# WHY NOT THE MAIN ANALYSIS' GENE RULE (TRDV1/TRDV2/TRGV9/TRDC > 0): TRGV9 is
# transcribed by alpha-beta T cells, because the gamma locus rearranges before
# alpha-beta commitment. Here it flags 858 cells carrying a full alpha-beta TCR.
################################################################################
joint_exclusion <- function(o) {
  m <- o@meta.data
  aa <- as.character(m$CTaa); aa[is.na(aa)] <- "NA"
  a_part <- sub("_.*$", "", aa)
  has_alpha <- a_part != "NA" & nzchar(a_part)
  gen <- as.character(m$CTgene); gen[is.na(gen)] <- ""
  mait <- grepl("(^|;)TRAV1-2\\.TRAJ(33|12|20)", sub("_.*$", "", gen))
  trdc <- as.numeric(GetAssayData(o, assay = "RNA", layer = "counts")["TRDC", rownames(m)])
  gd <- !has_alpha & trdc >= 2
  cl <- as.character(m$joint_clusters)
  why <- ifelse(cl == "10", "Heat-shock cluster",
         ifelse(cl == "8",  "CD4-expressing and low-quality cluster",
         ifelse(mait,       "MAIT (TRAV1-2 with TRAJ33/12/20)",
         ifelse(gd,         "Gamma-delta (TRDC, no TCR alpha chain)", NA))))
  setNames(why, rownames(m))
}
joint_keep <- function(o) { w <- joint_exclusion(o); names(w)[is.na(w)] }
