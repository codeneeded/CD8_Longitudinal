################################################################################
# Supplementary Tables S3 and S4: per-cluster tables for the manually gated
# central-memory (Tcm) and stem-cell memory (Tscm) populations
#
# The existing R1m1_01_CD8_central_memory_check.csv is a 7-row SUMMARY. A
# supplementary table needs the per-cluster breakdown, so both are rebuilt in
# that form here, with the summary rows kept.
#
# The Tscm gate is the one actually used for annotation in
# Scripts/2. Analysis/4.CD8_Subseting.R: within the naive clusters,
# CD45RO <= 0 AND FAS > 0 in DSB units.
#
# Out: Reviewer_Comments/Tables/TableS3_Tcm_profile_by_cluster.csv
#      Reviewer_Comments/Tables/TableS4_Tscm_gate_by_source_cluster.csv
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr)})
# Repository root: CD8_BASE if set, otherwise the first folder above this script
# that holds saved_R_data/ (works under Rscript, source() and an IDE).
R <- local({
  b <- Sys.getenv("CD8_BASE"); if (nzchar(b)) return(normalizePath(b))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f <- if (length(f)) gsub("~\\+~", " ", f[1]) else tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  d <- if (length(f) && nzchar(f)) dirname(normalizePath(f)) else getwd()
  while (!dir.exists(file.path(d, "saved_R_data")) && dirname(d) != d) d <- dirname(d)
  if (!dir.exists(file.path(d, "saved_R_data"))) stop("set CD8_BASE to the repository root", call. = FALSE)
  d
})
TAB <- file.path(R, "Reviewer_Comments", "Tables")

cd8 <- qs_read(file.path(R,"saved_R_data","TARA_cd8_HEI_annotated_final.qs2"))
DefaultAssay(cd8) <- "RNA"; if (length(Layers(cd8[["RNA"]]))>1) cd8 <- JoinLayers(cd8, assay="RNA")
rna <- FetchData(cd8, vars=c("CCR7","SELL","CD27"), layer="counts")
DefaultAssay(cd8) <- "ADT"
adt <- FetchData(cd8, vars=intersect(c("CD45RA","CD45RO","FAS"), rownames(cd8)), layer="data")
ann <- as.character(cd8$CD8_Annotation)
N   <- ncol(cd8)

# ---- Tcm, exactly as the existing summary table defines it --------------------
tcm_tbl <- rna$CCR7>0 & rna$SELL>0 & rna$CD27>0 & adt$CD45RA < 1
# ---- Tcm as the MANUSCRIPT sentence describes it (adds surface CD45RO+) -------
tcm_ms  <- tcm_tbl & adt$CD45RO > 0
cat("=== Tcm gate: table definition vs manuscript sentence ===\n")
cat(sprintf("  table  (CCR7+ SELL+ CD27+ RNA, CD45RA DSB < 1)          : %d cells (%.2f%%)\n",
            sum(tcm_tbl), 100*sum(tcm_tbl)/N))
cat(sprintf("  manuscript (same, PLUS surface CD45RO > 0)              : %d cells (%.2f%%)\n",
            sum(tcm_ms), 100*sum(tcm_ms)/N))

mk <- function(mask, label, fn, gate_text) {
  per <- tibble(cluster = ann, hit = mask) %>%
    group_by(cluster) %>%
    summarise(n_cells_in_cluster = n(), n_gated = sum(hit),
              pct_of_cluster = round(100*sum(hit)/n(), 2), .groups="drop") %>%
    arrange(desc(n_gated))
  per$pct_of_all_gated <- round(100*per$n_gated/max(1,sum(mask)), 2)
  out <- bind_rows(
    tibble(cluster="ALL CD8", n_cells_in_cluster=N, n_gated=sum(mask),
           pct_of_cluster=round(100*sum(mask)/N,2), pct_of_all_gated=100),
    per)
  write.csv(out, file.path(TAB, fn), row.names=FALSE)
  cat("\n===", label, "===\n"); cat("  gate:", gate_text, "\n")
  cat(sprintf("  %d of %d CD8 cells (%.2f%%), spanning %d of %d clusters\n",
      sum(mask), N, 100*sum(mask)/N, sum(per$n_gated>0), nrow(per)))
  print(as.data.frame(per), row.names=FALSE)
  cat("  wrote", fn, "\n")
}
mk(tcm_tbl, "Tcm transcriptional profile", "TableS3_Tcm_profile_by_cluster.csv",
   "CCR7+ SELL+ CD27+ by RNA detection, surface CD45RA < 1 (DSB)")

# ---- Tscm, the gate actually used for annotation ------------------------------
# A per-cluster table of the Tscm gate against the FINAL annotation is circular:
# the gate IS the annotation, so it returns 438 of 438 in one cluster. The
# informative table is what the gate did to the naive clusters it was applied to
# - i.e. the cluster each Tscm cell belonged to BEFORE gating. That column is
# source_cluster, the facet shown in Supplementary Figure S5C.
sc <- if ("source_cluster" %in% names(cd8@meta.data)) as.character(cd8$source_cluster) else NA
cat("\n=== Tscm surface gate, by the cluster each cell came FROM ===\n")
if (all(is.na(sc))) {
  cat("  source_cluster not present in the object; cannot build the pre-gate table\n")
} else {
  tscm_gate <- adt$CD45RO <= 0 & adt$FAS > 0
  gated <- grepl("Tscm", ann)
  d <- tibble(source = sc, final = ann,
              ro = adt$CD45RO, fas = adt$FAS) %>% filter(!is.na(source))
  per <- d %>% group_by(source) %>%
    summarise(n_cells_in_source_cluster = n(),
              n_to_Tscm          = sum(grepl("Tscm", final)),
              n_to_NaiveInt      = sum(grepl("Intermediate", final)),
              n_kept_original    = sum(!grepl("Tscm|Intermediate", final)),
              pct_to_Tscm        = round(100*sum(grepl("Tscm", final))/n(), 2),
              .groups="drop") %>% arrange(desc(n_to_Tscm))
  write.csv(per, file.path(TAB, "TableS4_Tscm_gate_by_source_cluster.csv"), row.names=FALSE)
  print(as.data.frame(per), row.names=FALSE)
  cat("  wrote TableS4_Tscm_gate_by_source_cluster.csv\n")
  cat("\n  ITEM 24 - what source_cluster means: the WNN cluster a cell belonged to BEFORE\n")
  cat("  the surface sub-gating reassigned it. Legend wording: 'Cluster of origin prior to\n")
  cat("  CD45RO/CD95 surface sub-gating'.\n")
}
cat("\n=== MS COMMENT 142: why Tscm is gated and Tcm is not ===\n")
cat("  CCR7 in the ADT panel:", "CCR7" %in% rownames(cd8[["ADT"]]), "\n")
for (mk_ in c("CD45RO","FAS","CD45RA")) cat(sprintf("  %-8s in ADT panel: %s\n", mk_, mk_ %in% rownames(cd8[["ADT"]])))
cat("  -> Tscm's defining surface markers (CD45RO, CD95/FAS) are on the panel, so it can be\n")
cat("     gated on protein. CCR7, which defines Tcm, is NOT on the panel, so Tcm can only be\n")
cat("     scored from RNA - and that profile does not resolve a separable population.\n")
