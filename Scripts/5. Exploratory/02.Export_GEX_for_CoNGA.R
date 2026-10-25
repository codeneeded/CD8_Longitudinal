################################################################################
# EXPLORATORY — write the CD8 GEX matrix for CoNGA as Matrix Market
#
# CoNGA reads 10x-style mtx or h5ad. Writing mtx from R avoids a SeuratDisk
# dependency; the companion python step converts it to h5ad.
# Only cells present in the clones file are written, so the two inputs agree.
################################################################################
suppressPackageStartupMessages({library(Seurat); library(qs2); library(Matrix)})
REPO <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
OUT  <- "/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/CoNGA"
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
