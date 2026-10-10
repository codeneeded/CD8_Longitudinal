################################################################################
# Revision — the 32 vRNA-positive cells on the published whole-PBMC UMAP
#
# ANSWERS: Lesley's comment 6 on the response document — "could you show the
#          distribution in a UMAP overlay? This would emphasize the point that
#          there are very few cells to analyze."
#
# An earlier version of this panel existed for the pilot participant only
# (7 cells, R1M3_04). This one carries every vRNA-positive cell found in the
# cohort realignment.
#
# JOIN KEY. The object's rownames are "<original library>_<barcode>" and the
# library prefix is NOT orig.ident (CP022_1m's cells are rownamed
# CP022_entry_<barcode>), so cells are matched on (orig.ident, bare barcode).
# Pasting the sample name onto the barcode matches nothing while looking clean.
#
# Out: Reviewer_Comments/Figures/R1M3_05_HIV_cohort_UMAP_overlay.png
#      + its source CSV
################################################################################

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(ggrepel)
})
source("config_revision.R")

pos <- read.csv(file.path(REV_TAB, "HIV_cohort_positive_cells.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
cat("vRNA-positive cells to place:", nrow(pos), "\n")

seu <- qs_read(QS_ALL)
md  <- seu@meta.data
md$barcode <- sub("^.*_", "", rownames(md))
stopifnot(all(grepl("^[ACGT]+-[0-9]+$", md$barcode)))

emb <- Embeddings(seu, reduction = "wnn.umap")
md$UMAP1 <- emb[, 1]; md$UMAP2 <- emb[, 2]

key_md  <- paste(as.character(md$orig.ident), md$barcode, sep = "|")
key_pos <- paste(pos$sample, pos$barcode, sep = "|")
ix <- match(key_pos, key_md)
cat("matched into the published object:", sum(!is.na(ix)), "of", nrow(pos), "\n")
if (!any(!is.na(ix)))
  stop("no vRNA-positive cell joined the object - check the barcode key", call. = FALSE)

hit <- pos[!is.na(ix), ]
hit$UMAP1 <- md$UMAP1[ix[!is.na(ix)]]
hit$UMAP2 <- md$UMAP2[ix[!is.na(ix)]]
hit$lineage <- as.character(md$Annotation[ix[!is.na(ix)]])

bg <- data.frame(UMAP1 = md$UMAP1, UMAP2 = md$UMAP2)

# One mark per cell, plain red, as the house rule requires: no size or shape
# encoding of burden, which needs extra legends and is better read off the table.
p <- ggplot() +
  geom_point(data = bg, aes(UMAP1, UMAP2), colour = "grey86", size = 0.25, stroke = 0) +
  geom_point(data = hit, aes(UMAP1, UMAP2), colour = "#C00000", size = 3.4, stroke = 0) +
  ggrepel::geom_text_repel(
    data = hit[hit$viral_UMI >= 3, ],
    aes(UMAP1, UMAP2, label = paste0(sample, " (", viral_UMI, ")")),
    size = 4.4, colour = "grey15", box.padding = 0.7, min.segment.length = 0,
    segment.colour = "grey55", max.overlaps = Inf, seed = 1) +
  coord_equal() +
  # The count must distinguish "found" from "placeable": 32 cells carry viral
  # reads, 17 of them survived the quality-control filtering that built the
  # annotated object and therefore have UMAP coordinates. Titling this "n = 17"
  # alone would read as though only 17 were detected.
  labs(title = paste0("HIV-1 transcript-positive cells on the published UMAP"),
       subtitle = paste0(nrow(hit), " of ", nrow(pos), " vRNA-positive cells retained in the ",
                         "annotated object (", format(nrow(md), big.mark = ","),
                         " cells); the other ", nrow(pos) - nrow(hit),
                         " were removed by quality control"),
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(size = 13, colour = "grey25"),
        axis.text = element_blank(), axis.ticks = element_blank())

save_fig("R1M3_05_HIV_cohort_UMAP_overlay", p, width = 10, height = 8.4)
save_src("R1M3_05_HIV_cohort_UMAP_overlay",
         hit[, c("sample", "barcode", "viral_UMI", "subB", "subC",
                 "lineage", "pct_of_transcriptome", "UMAP1", "UMAP2")])

cat("\nplaced", nrow(hit), "cells;",
    nrow(pos) - nrow(hit), "were removed by quality control before annotation\n")
cat("lineages carrying a vRNA-positive cell:\n")
print(sort(table(hit$lineage), decreasing = TRUE))
