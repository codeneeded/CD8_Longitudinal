################################################################################
# EXPLORATORY — export the CD8 object in the form CoNGA expects
#
# NOT FOR THE PAPER unless the result turns out to be interesting. CoNGA asks
# whether transcriptome neighbourhoods and TCR neighbourhoods are associated
# with each other more than chance allows; that is a different question from
# anything the reviewers raised, so this lives outside Scripts/2-4 and writes
# outside Analysis/ and Reviewer_Comments/.
#
# WHY WE CANNOT USE setup_10x_for_conga.py
#   That script reads per-sample filtered_contig_annotations.csv. Those files
#   are not on this machine (the same gap that stopped Figure_2.R). Everything
#   it would extract is already in the Seurat object's scRepertoire columns, so
#   the clones file is rebuilt from those instead.
#
#   CTgene  TRAV38-2/DV8.TRAJ47.TRAC_TRBV7-9.NA.TRBJ2-5.TRBC2
#   CTaa    CAVSVLSGGGADGLTF_CASSLSMNTEAFF
#   CTnt    TGTGCT...CTCACCTTT_TGTGCCAGCAGCTTGAGTATG...
#
#   Alpha and beta are separated by "_", segments within a chain by ".".
#   The beta slot is V . D . J . C with D left as NA, so J is the THIRD field.
#
# ALLELE SUFFIX. CoNGA's gene database is keyed on IMGT ids carrying an allele
# ("TRAV1*01"). scRepertoire reports the gene without one, so "*01" is appended
# and every gene is then checked against CoNGA's own db file; any cell with a
# gene CoNGA does not know is dropped rather than silently mis-scored.
#
# Out: Projects/Analysis_Exploratory/CoNGA/
################################################################################

suppressPackageStartupMessages({library(Seurat); library(qs2); library(dplyr)})

REPO <- "/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal"
PKGS <- "/home/ammasakshay/Desktop/Projects/Packages"
OUT  <- "/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/CoNGA"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

seu <- qs_read(file.path(REPO, "saved_R_data", "TARA_cd8_HEI_annotated_final.qs2"))
m   <- seu@meta.data
m$barcode <- rownames(m)
cat("CD8 cells:", nrow(m), "\n")

has <- function(x) !is.na(x) & grepl("_", x) & !grepl("^NA_|_NA$|^_|_$", x)
m <- m[has(m$CTgene) & has(m$CTaa) & has(m$CTnt), ]
cat("cells with a complete paired receptor:", nrow(m), "\n")

split2 <- function(x, i) vapply(strsplit(as.character(x), "_", fixed = TRUE),
                               function(p) if (length(p) >= i) p[[i]] else NA_character_, "")
seg    <- function(x, i) vapply(strsplit(as.character(x), ".", fixed = TRUE),
                               function(p) if (length(p) >= i) p[[i]] else NA_character_, "")

a_gene <- split2(m$CTgene, 1); b_gene <- split2(m$CTgene, 2)
m$va_gene <- seg(a_gene, 1); m$ja_gene <- seg(a_gene, 2)
m$vb_gene <- seg(b_gene, 1); m$jb_gene <- seg(b_gene, 3)
m$cdr3a <- split2(m$CTaa, 1); m$cdr3b <- split2(m$CTaa, 2)
m$cdr3a_nucseq <- split2(m$CTnt, 1); m$cdr3b_nucseq <- split2(m$CTnt, 2)

# MULTI-CHAIN CELLS. scRepertoire joins a second alpha or beta with ";" in
# CTaa/CTnt, but CTgene reports only the FIRST V/J pair for that chain - so the
# fields are not parallel and the second CDR3 cannot be matched to a V gene.
# Taking the first element of every ";" list keeps V, J and CDR3 internally
# consistent (827 clones carry a second alpha, 241 a second beta). CoNGA's
# cdr3a2 slot is left empty rather than filled with an unassignable sequence.
first_of <- function(x) vapply(strsplit(as.character(x), ";", fixed = TRUE),
                              function(p) p[[1]], "")
n_dual_a <- sum(grepl(";", m$cdr3a)); n_dual_b <- sum(grepl(";", m$cdr3b))
for (cc in c("cdr3a","cdr3b","cdr3a_nucseq","cdr3b_nucseq")) m[[cc]] <- first_of(m[[cc]])
cat("cells with a second alpha:", n_dual_a, " second beta:", n_dual_b,
    "- primary chain kept\n")

need <- c("va_gene","ja_gene","vb_gene","jb_gene",
          "cdr3a","cdr3b","cdr3a_nucseq","cdr3b_nucseq")
m <- m[complete.cases(m[, need]), ]
cat("with all eight fields parsed:", nrow(m), "\n")

# The object carries BOTH spellings of the dual TRAV/TRDV genes - "TRAV14/DV4"
# and "TRAV14DV4" - because scRepertoire passes through whatever Cell Ranger
# wrote. CoNGA's db uses the slashed form only, so the slashless one is
# normalised rather than discarded; it cost 841 cells on the first run.
fix_dv <- function(g) sub("^(TRAV[0-9]+(?:-[0-9]+)?)(DV[0-9]+)$", "\\1/\\2", g, perl = TRUE)
add_allele <- function(g) ifelse(grepl("\\*", g), g, paste0(g, "*01"))
for (cc in c("va_gene","ja_gene","vb_gene","jb_gene"))
  m[[cc]] <- add_allele(fix_dv(m[[cc]]))

db_path <- file.path(PKGS, "conga-master", "conga", "tcrdist", "db", "combo_xcr.tsv")
stopifnot(file.exists(db_path))
db    <- read.delim(db_path, stringsAsFactors = FALSE)
known <- db$id[db$organism == "human"]
cat("human genes in CoNGA db:", length(known), "\n")

gene_ok <- m$va_gene %in% known & m$ja_gene %in% known &
           m$vb_gene %in% known & m$jb_gene %in% known
bad <- setdiff(unique(unlist(m[!gene_ok, c("va_gene","ja_gene","vb_gene","jb_gene")])), known)
cat("cells dropped for an unrecognised gene:", sum(!gene_ok), "\n")
if (length(bad)) cat("  unrecognised:", paste(head(sort(bad), 12), collapse = ", "), "\n")
m <- m[gene_ok, ]
cat("cells into the clones file:", nrow(m), "\n")

m$clone_key <- paste(m$va_gene, m$ja_gene, m$cdr3a,
                     m$vb_gene, m$jb_gene, m$cdr3b, sep = "|")
cl <- m %>% group_by(clone_key) %>%
  summarise(subject      = names(sort(table(orig.ident), decreasing = TRUE))[1],
            clone_size   = n(),
            va_gene = first(va_gene), ja_gene = first(ja_gene),
            vb_gene = first(vb_gene), jb_gene = first(jb_gene),
            cdr3a   = first(cdr3a),   cdr3b   = first(cdr3b),
            cdr3a_nucseq = first(cdr3a_nucseq),
            cdr3b_nucseq = first(cdr3b_nucseq), .groups = "drop") %>%
  arrange(desc(clone_size))
cl$clone_id <- paste0("clone_", seq_len(nrow(cl)))
cl$va2_gene <- "-"; cl$ja2_gene <- "-"; cl$cdr3a2 <- "-"; cl$cdr3a2_nucseq <- "-"

cols <- c("clone_id","subject","clone_size","va_gene","ja_gene","va2_gene",
          "ja2_gene","vb_gene","jb_gene","cdr3a","cdr3a_nucseq","cdr3a2",
          "cdr3a2_nucseq","cdr3b","cdr3b_nucseq")
clones_file <- file.path(OUT, "tara_cd8_clones.tsv")
write.table(cl[, cols], clones_file, sep = "\t", quote = FALSE, row.names = FALSE)

bc <- m %>% left_join(cl[, c("clone_key","clone_id")], by = "clone_key") %>%
  group_by(clone_id) %>%
  summarise(barcodes = paste(barcode, collapse = ","), .groups = "drop")
write.table(bc, paste0(clones_file, ".barcode_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

meta_cols <- intersect(c("barcode","orig.ident","PID","Condition","Cohort",
                         "Timepoint_Group","CD8_Annotation","Annotation",
                         "cloneSize","CTstrict"), names(m))
write.csv(m[, meta_cols], file.path(OUT, "tara_cd8_cellmeta.csv"), row.names = FALSE)

cat("\nclones:", nrow(cl), " cells:", nrow(m),
    " largest clone:", max(cl$clone_size),
    " expanded (>1 cell):", sum(cl$clone_size > 1), "\n")
