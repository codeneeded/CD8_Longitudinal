################################################################################
# EXPLORATORY — write the CD8 GEX matrix for CoNGA as Matrix Market
#
# CoNGA reads 10x-style mtx or h5ad. Writing mtx from R avoids a SeuratDisk
# dependency; the companion python step converts it to h5ad.
# Only cells present in the clones file are written, so the two inputs agree.
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(Matrix)})
# Repository root: CD8_BASE if set, otherwise the first folder above this script
# that holds saved_R_data/ (works under Rscript, source() and an IDE).
REPO <- local({
  b <- Sys.getenv("CD8_BASE"); if (nzchar(b)) return(normalizePath(b))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f <- if (length(f)) gsub("~\\+~", " ", f[1]) else tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  d <- if (length(f) && nzchar(f)) dirname(normalizePath(f)) else getwd()
  while (!dir.exists(file.path(d, "saved_R_data")) && dirname(d) != d) d <- dirname(d)
  if (!dir.exists(file.path(d, "saved_R_data"))) stop("set CD8_BASE to the repository root", call. = FALSE)
  d
})
PROJ <- dirname(dirname(REPO))   # the projects folder that holds Repositories/
OUT  <- file.path(PROJ, "Analysis_Exploratory/CoNGA")
MTX  <- file.path(OUT, "gex"); dir.create(MTX, recursive = TRUE, showWarnings = FALSE)

seu <- qs_read(file.path(REPO, "saved_R_data", "TARA_cd8_HEI_annotated_final.qs2"))
bc  <- read.delim(file.path(OUT, "tara_cd8_clones.tsv.barcode_mapping.tsv"),
                  stringsAsFactors = FALSE)
keep <- unique(unlist(strsplit(bc$barcodes, ",", fixed = TRUE)))
keep <- intersect(keep, colnames(seu))
cat("cells in clones file:", length(keep), "of", ncol(seu), "\n")

DefaultAssay(seu) <- "RNA"
if (length(Layers(seu[["RNA"]])) > 1) seu <- JoinLayers(seu, assay = "RNA")
cnt <- GetAssayData(seu, assay = "RNA", layer = "counts")[, keep, drop = FALSE]
cat("matrix:", nrow(cnt), "genes x", ncol(cnt), "cells\n")

writeMM(cnt, file.path(MTX, "matrix.mtx"))
write.table(data.frame(rownames(cnt), rownames(cnt), "Gene Expression"),
            file.path(MTX, "features.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(cnt), file.path(MTX, "barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("wrote matrix.mtx, features.tsv, barcodes.tsv\n")
