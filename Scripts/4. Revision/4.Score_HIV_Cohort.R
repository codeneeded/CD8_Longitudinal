################################################################################
# Revision step 4 — score HIV detection across the entry-visit realignment
#
# ANSWERS: Reviewer 1, major point 3 (cohort-wide, after 2b.Align_Entry_HIV.sh).
#
# WHY THIS IS SEPARATE FROM 3.Score_HIV_Detection.R
#   That script scores ONE pilot directory (pilot_hiv_<sample>) and is left
#   untouched so the pilot result stays reproducible. This one reads every
#   hiv_<sample>/ produced by 2b and scores them together. Single purpose,
#   separately runnable: it reads alignment output and writes tables, and does
#   no alignment of its own.
#
# THE THREE CONTROLS, unchanged in meaning from the pilot
#   CONTROL 1  ambient      viral UMIs inside called cells vs in empty droplets.
#                           A sample whose viral UMIs are mostly OUTSIDE cells
#                           is ambient background, not infection.
#   CONTROL 2  lineage      HIV should sit in CD4 T cells. Lineage comes from
#                           the project's own annotation (TARA_ALL), joined by
#                           barcode - cellranger 10.1.0 `count` writes no
#                           cell_types.csv, and re-clustering these libraries
#                           would invent a second, inconsistent annotation.
#   CONTROL 3  dose         cells carrying more viral UMIs should be MORE
#                           lineage-restricted. If 1-UMI and 10-UMI cells look
#                           alike, the signal is noise.
#
# READING THE OUTPUT
#   A sample with few vRNA UMIs and a high ambient fraction is negative at
#   this detection ceiling, which is not absence of infected cells.
#   Detection is bounded by what 3' short-read GEX captures; the pilot's
#   own numbers (CP006_1m: 761 UMIs, 90% cell-associated, 19 cells) are the
#   positive reference to compare against.
#
# MEMORY
#   The raw matrix is ~38,608 x ~473,000 and would be several GB if loaded.
#   Only two rows are ever needed, so the .mtx is streamed through awk and only
#   the viral entries reach R. Nothing large is held.
#
# INPUT   cellranger_out/hiv_<sample>/outs/{raw,filtered}_feature_bc_matrix
#         saved_R_data/TARA_ALL_annotated_final.qs2   (annotation only)
# OUTPUT  Reviewer_Comments/Tables/HIV_cohort_per_sample.csv
#         Reviewer_Comments/Tables/HIV_cohort_positive_cells.csv
#         Reviewer_Comments/Tables/HIV_cohort_controls.csv
#         Reviewer_Comments/Figures/R1M3_01_HIV_cohort_detection.png
#
# RUN     Rscript "Scripts/4. Revision/4.Score_HIV_Cohort.R"
#         OUT_PARENT=/path/to/cellranger_out Rscript ...   # if elsewhere
################################################################################

.cfg <- Sys.getenv("CD8_REV_CONFIG", "")
if (!nzchar(.cfg)) {
  for (b in c("~/Documents/CD8_Longitudinal",
              "~/Desktop/Projects/Repositories/CD8_Longitudinal")) {
    p <- file.path(path.expand(b), "Scripts", "4. Revision", "config_revision.R")
    if (file.exists(p)) { .cfg <- p; break }
  }
}
if (!nzchar(.cfg)) stop("cannot find config_revision.R; set CD8_REV_CONFIG", call. = FALSE)
source(.cfg)

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(ggplot2); library(tidyr)
})

CR_OUT <- Sys.getenv("OUT_PARENT", path.expand("~/Desktop/Projects/cellranger_out"))
MIN_UMI_HI <- as.integer(Sys.getenv("HIV_MIN_UMI", "3"))   # CONTROL 3 threshold

dirs <- sort(Sys.glob(file.path(CR_OUT, "hiv_*", "outs")))
if (!length(dirs))
  stop("no hiv_*/outs directories under ", CR_OUT,
       "\n  run 2b.Align_Entry_HIV.sh first, or set OUT_PARENT.", call. = FALSE)
samples <- sub("^hiv_", "", basename(dirname(dirs)))
cat("samples found:", length(samples), "-", paste(samples, collapse = ", "), "\n\n")

# ---------------------------------------------------------------- readers ----
# Return the viral UMI count per barcode for one matrix directory, WITHOUT
# loading the matrix. Streams the .mtx and keeps only the viral feature rows.
viral_per_barcode <- function(mdir) {
  fpath <- file.path(mdir, "features.tsv.gz")
  bpath <- file.path(mdir, "barcodes.tsv.gz")
  mpath <- file.path(mdir, "matrix.mtx.gz")
  stopifnot(file.exists(fpath), file.exists(bpath), file.exists(mpath))

  ft  <- read.delim(gzfile(fpath), header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
  idx <- which(grepl("^HIV1_", ft[[1]]))
  names(idx) <- ft[[1]][idx]
  # A missing contig yields a convincing zero, so refuse rather than report it.
  if (length(idx) != 2)
    stop(mdir, ": expected 2 viral features, found ", length(idx),
         ". This matrix was not built against GRCh38_HIV1.", call. = FALSE)

  bc <- readLines(gzfile(bpath))

  # awk: drop comment lines, drop the dimension line, keep only viral rows.
  cmd <- sprintf(
    "zcat %s | awk 'BEGIN{d=0} /^%%/ {next} d==0 {d=1; next} ($1==%d || $1==%d) {print $1\" \"$2\" \"$3}'",
    shQuote(mpath), idx[1], idx[2])
  con <- pipe(cmd, open = "r")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  df <- tryCatch(
    utils::read.table(con, col.names = c("row", "col", "val")),
    error = function(e) data.frame(row = integer(), col = integer(), val = integer()))

  list(n_barcodes = length(bc), barcodes = bc, idx = idx, entries = df)
}

# ------------------------------------------------------- per-sample pass ----
per_sample <- list(); pos_cells <- list()

for (i in seq_along(dirs)) {
  s <- samples[i]
  fdir <- file.path(dirs[i], "filtered_feature_bc_matrix")
  rdir <- file.path(dirs[i], "raw_feature_bc_matrix")
  cat("[", s, "] reading\n", sep = "")

  f <- viral_per_barcode(fdir)
  r <- viral_per_barcode(rdir)

  contig <- names(f$idx)
  f_by_contig <- tapply(f$entries$val, factor(f$entries$row, levels = f$idx), sum)
  r_by_contig <- tapply(r$entries$val, factor(r$entries$row, levels = r$idx), sum)
  f_by_contig[is.na(f_by_contig)] <- 0
  r_by_contig[is.na(r_by_contig)] <- 0

  in_cell   <- sum(f$entries$val)
  all_bc    <- sum(r$entries$val)
  ambient   <- all_bc - in_cell
  pos_bc    <- length(unique(f$entries$col))

  per_sample[[s]] <- data.frame(
    sample            = s,
    n_cells           = f$n_barcodes,
    n_barcodes_raw    = r$n_barcodes,
    viral_UMI_cell_assoc= in_cell,
    viral_UMI_all     = all_bc,
    viral_UMI_ambient = ambient,
    pct_ambient       = if (all_bc > 0) round(100 * ambient / all_bc, 1) else NA_real_,
    n_viral_pos_cells = pos_bc,
    pct_cells_viral   = round(100 * pos_bc / f$n_barcodes, 4),
    subB_UMI          = as.integer(f_by_contig[[1]]),
    subC_UMI          = as.integer(f_by_contig[[2]]),
    stringsAsFactors  = FALSE)

  if (nrow(f$entries)) {
    pc <- f$entries %>%
      group_by(col) %>%
      summarise(viral_UMI = sum(val),
                subB = sum(val[row == f$idx[1]]),
                subC = sum(val[row == f$idx[2]]), .groups = "drop") %>%
      # No composite cell_id here: the object's rowname prefix is the original
      # library name, not the sample, so a pasted key would not join.
      mutate(sample = s, barcode = f$barcodes[col]) %>%
      select(sample, barcode, viral_UMI, subB, subC)
    pos_cells[[s]] <- as.data.frame(pc)
  }
}

samp_tab <- bind_rows(per_sample)
pos_tab  <- if (length(pos_cells)) bind_rows(pos_cells) else
            data.frame(sample = character(), barcode = character(),
                       viral_UMI = integer(), subB = integer(), subC = integer())

# ------------------------------------------- CONTROL 2: lineage, by join ----
# The realigned libraries are the SAME libraries already in TARA_ALL, so the
# published annotation applies directly. Joining is strictly better than
# re-clustering these nine alone, which would produce labels inconsistent with
# every other panel. Cells absent from the object were removed by QC there and
# are reported as such rather than silently dropped.
cat("\njoining the project annotation for lineage\n")
seu <- qs_read(QS_ALL)
md  <- seu@meta.data

# JOIN KEY: (orig.ident, bare barcode) — NOT paste0(sample, "_", barcode).
#
# The object's rownames are "<library>_<barcode>", and the library prefix is
# the ORIGINAL library name, not orig.ident: CP022_1m's cells are rownamed
# CP022_entry_<barcode>, CP017_2m's are CP017_entry_<barcode>. Only 38,478 of
# 126,435 rownames happen to start with their own orig.ident, so building the
# key by pasting the sample name matched ZERO of 32 viral-positive cells and
# would have silently reported every one as absent from the object — killing
# CONTROL 2 while looking like a clean run. Strip to the bare barcode instead
# and join within sample; the barcode is unique inside each orig.ident.
md$barcode <- sub("^.*_", "", rownames(md))
stopifnot(all(grepl("^[ACGT]+-[0-9]+$", md$barcode)))

lin_col <- if ("Annotation" %in% names(md)) "Annotation" else "predicted.celltype.l1"
keep <- c("orig.ident", "barcode", lin_col, "Broad_CellType", "nCount_RNA")
keep <- intersect(keep, names(md))
ann <- md[, keep]
names(ann)[names(ann) == lin_col] <- "lineage"
names(ann)[names(ann) == "orig.ident"] <- "sample"
ann$sample <- as.character(ann$sample)
if (anyDuplicated(paste(ann$sample, ann$barcode)))
  stop("barcodes are not unique within a sample; the join key is unsafe.", call. = FALSE)
rm(seu, md); invisible(gc(verbose = FALSE))

pos_tab <- left_join(pos_tab, ann, by = c("sample", "barcode"))
pos_tab$in_published_object <- !is.na(pos_tab$lineage)
# A join that silently matches nothing is the failure mode this guards against.
if (nrow(pos_tab) && !any(pos_tab$in_published_object))
  stop("the annotation join matched 0 of ", nrow(pos_tab), " viral-positive cells. ",
       "Check the barcode format in both sources before trusting CONTROL 2.", call. = FALSE)
cat("  annotation joined for", sum(pos_tab$in_published_object), "of", nrow(pos_tab),
    "viral-positive cells;", sum(!pos_tab$in_published_object),
    "were removed by QC in the published object\n")
pos_tab$pct_of_transcriptome <- ifelse(
  is.na(pos_tab$nCount_RNA), NA_real_,
  round(100 * pos_tab$viral_UMI / pos_tab$nCount_RNA, 3))

# ------------------------------------------------- CONTROL 3: dose series ----
dose <- pos_tab %>%
  filter(in_published_object) %>%
  mutate(bin = ifelse(viral_UMI >= MIN_UMI_HI,
                      paste0(">=", MIN_UMI_HI, " UMI"), "1-2 UMI")) %>%
  group_by(bin, lineage) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(bin) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
  arrange(bin, desc(n))

controls <- bind_rows(
  samp_tab %>% transmute(
    control = "1 ambient", sample,
    value = paste0(viral_UMI_cell_assoc, " of ", viral_UMI_all,
                   " vRNA UMIs cell-associated (", 100 - pct_ambient, "%)"),
    verdict = ifelse(is.na(pct_ambient), "no vRNA detected",
              ifelse(viral_UMI_all < 10, "below evaluable threshold",
              ifelse(pct_ambient <= 25, "cell-associated", "ambient-dominated")))),
  pos_tab %>% filter(in_published_object) %>%
    group_by(sample) %>%
    summarise(n_pos = n(),
              n_T = sum(grepl("CD4|CD8|T cell|Treg|MAIT|gd", lineage, ignore.case = TRUE)),
              .groups = "drop") %>%
    transmute(control = "2 lineage", sample,
              value = paste0(n_T, " of ", n_pos, " viral+ cells are T cells"),
              verdict = ifelse(n_pos < 3, "below evaluable threshold",
                        ifelse(n_T / n_pos >= 0.5, "T-cell restricted",
                               "not lineage restricted"))))

# -------------------------------------------------------------- outputs ----
mkdirs(REV_TAB, REV_FIG)
write.csv(samp_tab, file.path(REV_TAB, "HIV_cohort_per_sample.csv"), row.names = FALSE)
write.csv(pos_tab,  file.path(REV_TAB, "HIV_cohort_positive_cells.csv"), row.names = FALSE)
write.csv(controls, file.path(REV_TAB, "HIV_cohort_controls.csv"), row.names = FALSE)
write.csv(dose,     file.path(REV_TAB, "HIV_cohort_dose_response.csv"), row.names = FALSE)

pl <- samp_tab %>%
  mutate(sample = factor(sample, levels = sample[order(-viral_UMI_cell_assoc)])) %>%
  select(sample, `Cell-associated` = viral_UMI_cell_assoc,
         `Ambient` = viral_UMI_ambient) %>%
  pivot_longer(-sample, names_to = "where", values_to = "UMI")

p <- ggplot(pl, aes(sample, UMI, fill = where)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Cell-associated" = "#C00000",
                               "Ambient" = "grey70")) +
  labs(x = NULL, y = "HIV-1 UMIs", fill = NULL,
       title = "HIV-1 transcript detection at the entry visit") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top")
save_fig("R1M3_01_HIV_cohort_detection", p, width = 9, height = 6)
save_src("R1M3_01_HIV_cohort_detection", pl)

# --------------------------------------------------------------- report ----
cat("\n================ PER SAMPLE ================\n")
print(samp_tab[, c("sample","n_cells","viral_UMI_cell_assoc","viral_UMI_ambient",
                   "pct_ambient","n_viral_pos_cells","subB_UMI","subC_UMI")],
      row.names = FALSE)
cat("\n================ CONTROLS ==================\n")
print(as.data.frame(controls), row.names = FALSE)
cat("\n================ CONTROL 3: dose response ==\n")
print(as.data.frame(dose), row.names = FALSE)
cat("\nviral+ cells total:", nrow(pos_tab),
    "| with published annotation:", sum(pos_tab$in_published_object), "\n")
cat("wrote 4 tables to", REV_TAB, "\n")
cat("\nSamples below the evaluable threshold are negative at this detection\n",
    "ceiling; that is not equivalent to absence of infected cells.\n", sep = "")
