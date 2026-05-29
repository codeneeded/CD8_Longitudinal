# ═══════════════════════════════════════════════════════════════════════════
# 01_threshold_inspection.R
#
# Determines the BCR clonal clustering threshold by FIRST inspecting the shape
# of the heavy-chain distance-to-nearest distribution, then SELECTING the
# method appropriate to that shape:
#
#   bimodal  (a low-distance clonal mode exists)  -> find the valley
#            via density; fall back to gmm; take the value in a sane range.
#   unimodal (no clonal mode — only the "unrelated" hump) -> do NOT run gmm
#            (it fails to converge on one-component data and hangs); report
#            minimal clonal expansion and use a conservative fixed threshold.
#
# WHY SHAPE-FIRST: on this cohort the distribution proved UNIMODAL (single hump
# centred ~0.31, nothing at low distance). GMM hung trying to fit two components
# to one. Shared V/J/junction-length groups existed but their CDR3s differ by
# ~30% — convergent recombination, not clonal expansion. Choosing the method
# from the shape avoids the hang and picks a defensible threshold automatically.
#
# This script writes thresholds.tsv (consumed by the clone-calling script) and
# prints a methods-ready verdict. It does not call clones.
# ═══════════════════════════════════════════════════════════════════════════

library(dplyr)
library(ggplot2)
library(patchwork)
library(alakazam)
library(shazam)

# ── Paths ──────────────────────────────────────────────────────────────────
data_path <- "~/Documents/CD8_Longitudinal/VDJ/BCR"                 # Stage 1 output (read)
base_path <- "~/Documents/CD8_Longitudinal/Analysis/VDJ/BCR"        # analysis home (write)
meta_file <- file.path(base_path, "sample_metadata.csv")
fig_path  <- file.path(base_path, "threshold")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

# ── Decision parameters (edit if you disagree) ───────────────────────────────
# Region searched for a clonal mode / valley:
CLONAL_MAX     <- 0.15   # a clonal mode, if present, should peak below this
SANE_MIN       <- 0.05   # accept an empirical valley only within
SANE_MAX       <- 0.25   #   [SANE_MIN, SANE_MAX]
# Minimum fraction of cells below CLONAL_MAX to call the distribution bimodal:
BIMODAL_MIN_FRAC <- 0.05 # >=5% of neighboured cells in the low region
# Fixed fallback when unimodal (conservative; won't over-merge):
FIXED_THRESHOLD <- 0.10

CANDIDATES <- c(0.10, 0.12, 0.15, 0.18, 0.20)


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 — LOAD + QC  (identical to the main clone-calling script)
# ═══════════════════════════════════════════════════════════════════════════

.probe <- list.files(data_path, pattern = "_db-pass\\.tsv$",
                     recursive = TRUE, full.names = TRUE)
cat("Stage 1 output under:", path.expand(data_path), "\n")
cat("Found", length(.probe), "_db-pass.tsv files\n")
if (length(.probe) == 0) stop("No *_db-pass.tsv found; check data_path.")

meta <- read.csv(meta_file, stringsAsFactors = FALSE)

load_dbpass <- function(sample_id) {
  f <- file.path(data_path, sample_id, paste0(sample_id, "_db-pass.tsv"))
  if (!file.exists(f)) { cat("  MISSING:", f, "\n"); return(NULL) }
  db <- read.delim(f, stringsAsFactors = FALSE)
  if (nrow(db) == 0) return(NULL)
  db$sample_id <- sample_id
  db
}

cat("Loading and QC-filtering...\n")
bcr <- bind_rows(lapply(meta$sample, load_dbpass)) %>%
  left_join(meta, by = c("sample_id" = "sample"))
bcr$cell_id_unique <- paste(bcr$sample_id, bcr$cell_id, sep = "_")

bcr <- bcr %>% filter(
  (grepl("^IGHV", v_call) & grepl("^IGHJ", j_call) & (is.na(c_call) | grepl("^IGH", c_call))) |
    (grepl("^IGKV", v_call) & grepl("^IGKJ", j_call) & (is.na(c_call) | grepl("^IGK", c_call))) |
    (grepl("^IGLV", v_call) & grepl("^IGLJ", j_call) & (is.na(c_call) | grepl("^IGL", c_call)))
) %>% filter(productive)
multi_heavy <- table(filter(bcr, locus == "IGH")$cell_id_unique)
bcr <- bcr %>% filter(!cell_id_unique %in% names(multi_heavy)[multi_heavy > 1])
heavy_cells <- filter(bcr, locus == "IGH")$cell_id_unique
light_cells <- filter(bcr, locus %in% c("IGK", "IGL"))$cell_id_unique
bcr <- bcr %>% filter(!cell_id_unique %in% setdiff(light_cells, heavy_cells))

bcr_heavy <- bcr %>% filter(locus == "IGH")
cat("Clean heavy chains:", nrow(bcr_heavy),
    "across", length(unique(bcr_heavy$participant)), "participants\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2 — DISTANCE TO NEAREST
# ═══════════════════════════════════════════════════════════════════════════

cat("Computing distance-to-nearest (pooled, single-cell mode)...\n")
dtn <- distToNearest(
  bcr_heavy,
  sequenceColumn = "junction", vCallColumn = "v_call", jCallColumn = "j_call",
  model = "ham", normalize = "len",
  cellIdColumn = "cell_id_unique", locusColumn = "locus",
  onlyHeavy = TRUE, nproc = 4
)
d <- dtn$dist_nearest[!is.na(dtn$dist_nearest)]
cat("Cells with a nearest-neighbour distance:", length(d),
    "of", nrow(dtn), "\n")
cat("Distance summary:\n"); print(summary(d)); cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 — INSPECT SHAPE, THEN CHOOSE METHOD
# ───────────────────────────────────────────────────────────────────────────
# NOTE: We decide bimodal vs unimodal from the data, then route to the right
# method. The test: is there appreciable mass below CLONAL_MAX, AND a density
# valley separating a low mode from the main hump?
# ═══════════════════════════════════════════════════════════════════════════

frac_low <- mean(d < CLONAL_MAX)
cat(sprintf("Fraction of neighboured cells with distance < %.2f: %.3f\n",
            CLONAL_MAX, frac_low))

# Smooth density, look for a local minimum (valley) in the clonal search window
dens <- density(d, n = 512)
# candidate valley: local minima of the density within [SANE_MIN, SANE_MAX]
in_win <- dens$x >= SANE_MIN & dens$x <= SANE_MAX
valley_x <- NA_real_
if (sum(in_win) > 3) {
  y <- dens$y[in_win]; x <- dens$x[in_win]
  # local minima: lower than both neighbours
  is_min <- c(FALSE, y[-c(1, length(y))] < y[-c(1,2)] &
                y[-c(1, length(y))] < y[-c(length(y)-1, length(y))], FALSE)
  if (any(is_min)) valley_x <- x[which(is_min)[1]]
}

# Decide shape
is_bimodal <- (frac_low >= BIMODAL_MIN_FRAC) && !is.na(valley_x)

cat("\n──────────────────────────────────────────────\n")
cat("SHAPE DIAGNOSIS:\n")
cat(sprintf("  low-distance fraction (< %.2f): %.3f  (need >= %.2f)\n",
            CLONAL_MAX, frac_low, BIMODAL_MIN_FRAC))
cat(sprintf("  density valley in [%.2f, %.2f]: %s\n",
            SANE_MIN, SANE_MAX, ifelse(is.na(valley_x), "none found",
                                       round(valley_x, 4))))
cat(sprintf("  => distribution looks: %s\n", ifelse(is_bimodal, "BIMODAL", "UNIMODAL")))
cat("──────────────────────────────────────────────\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 — SELECT THRESHOLD BY THE DIAGNOSED SHAPE
# ═══════════════════════════════════════════════════════════════════════════

chosen     <- NA_real_
method_used <- NA_character_
verdict     <- ""

if (is_bimodal) {
  # ── BIMODAL branch: find the valley. Density first, gmm as backup. ──
  cat("Bimodal -> determining valley (density, then gmm backup)...\n")
  
  thr_den <- tryCatch(findThreshold(d, method = "density"),
                      error = function(e) { cat("  density error:", e$message, "\n"); NULL })
  den_val <- if (!is.null(thr_den)) thr_den@threshold else NA_real_
  
  if (!is.na(den_val) && den_val >= SANE_MIN && den_val <= SANE_MAX) {
    chosen <- den_val; method_used <- "density"
  } else {
    # gmm is only attempted in the bimodal branch (safe: two components exist)
    cat("  density value out of range; trying gmm...\n")
    thr_gmm <- tryCatch(
      findThreshold(d, method = "gmm", model = "gamma-norm"),
      error = function(e) { cat("  gmm error:", e$message, "\n"); NULL })
    gmm_val <- if (!is.null(thr_gmm)) thr_gmm@threshold else NA_real_
    if (!is.na(gmm_val) && gmm_val >= SANE_MIN && gmm_val <= SANE_MAX) {
      chosen <- gmm_val; method_used <- "gmm"
    } else if (!is.na(valley_x)) {
      chosen <- valley_x; method_used <- "density_valley_manual"
    } else {
      chosen <- FIXED_THRESHOLD; method_used <- "fixed_fallback"
    }
  }
  verdict <- sprintf(
    "Bimodal nearest-neighbour distribution; clonal threshold set to %.3f (%s).",
    chosen, method_used)
  
} else {
  # ── UNIMODAL branch: NO clonal mode. Do NOT run gmm. ──
  cat("Unimodal -> no clonal mode detected. Skipping gmm (would not converge).\n")
  chosen <- FIXED_THRESHOLD; method_used <- "fixed_unimodal"
  verdict <- paste0(
    "The heavy-chain nearest-neighbour distance distribution was unimodal ",
    "with no low-distance clonal mode (median ", round(median(d), 2),
    "; ", round(100 * frac_low, 1), "% of cells below ", CLONAL_MAX, "), ",
    "indicating minimal clonal expansion. A conservative fixed clonal ",
    "distance threshold of ", FIXED_THRESHOLD, " (normalized Hamming) was ",
    "applied; shared V/J/junction-length groupings reflect convergent ",
    "recombination rather than clonal relatedness.")
}

cat("\n════════════════════════════════════════════════════════\n")
cat(sprintf("CHOSEN THRESHOLD: %.4f   (method: %s)\n", chosen, method_used))
cat("════════════════════════════════════════════════════════\n")
cat("VERDICT (methods-ready):\n", verdict, "\n\n", sep = "")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 5 — PER-THRESHOLD COUNTS (context for the supplement)
# ═══════════════════════════════════════════════════════════════════════════

counts <- sapply(CANDIDATES, function(t) sum(d <= t))
count_tab <- data.frame(threshold = CANDIDATES,
                        cells_with_neighbour = counts,
                        pct_of_neighboured = round(100 * counts / length(d), 2))
cat("Cells gaining a clonal relative at each candidate threshold:\n")
print(count_tab); cat("\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 6 — PLOTS (PNG + PDF)
# ═══════════════════════════════════════════════════════════════════════════

plot_df <- data.frame(dist_nearest = d)

p_full <- ggplot(plot_df, aes(x = dist_nearest)) +
  geom_histogram(binwidth = 0.01, fill = "grey75", colour = "white") +
  geom_vline(xintercept = chosen, colour = "#E41A1C", linewidth = 1) +
  annotate("text", x = chosen + 0.01, y = Inf,
           label = sprintf("threshold = %.3f (%s)", chosen, method_used),
           hjust = 0, vjust = 1.5, colour = "#E41A1C", size = 3.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Distance-to-nearest (full range)",
       subtitle = ifelse(is_bimodal, "Diagnosed BIMODAL", "Diagnosed UNIMODAL — no clonal mode"),
       x = "Normalized Hamming distance", y = "Count") +
  theme_bw(base_size = 11)

p_zoom <- ggplot(subset(plot_df, dist_nearest <= 0.4), aes(x = dist_nearest)) +
  geom_histogram(binwidth = 0.005, fill = "#4472C4", colour = "white") +
  geom_vline(xintercept = CANDIDATES, colour = "grey60", linetype = "dotted") +
  geom_vline(xintercept = chosen, colour = "#E41A1C", linewidth = 1) +
  scale_x_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.05)) +
  labs(title = "Zoomed 0-0.4 (clonal search region)",
       subtitle = "dotted grey = candidate thresholds; red = chosen",
       x = "Normalized Hamming distance", y = "Count") +
  theme_bw(base_size = 11)

combined <- p_full / p_zoom
ggsave(file.path(fig_path, "threshold_inspection.png"), combined, width = 9, height = 9, dpi = 150)
ggsave(file.path(fig_path, "threshold_inspection.pdf"), combined, width = 9, height = 9)
cat("Saved: threshold_inspection.png + .pdf\n")


# ═══════════════════════════════════════════════════════════════════════════
# STEP 7 — WRITE thresholds.tsv + verdict
# ═══════════════════════════════════════════════════════════════════════════

all_samples <- unique(bcr_heavy$sample_id)
thresholds_out <- data.frame(sample = all_samples, threshold = round(chosen, 4))
write.table(thresholds_out, file = file.path(base_path, "thresholds.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

writeLines(c(
  paste0("chosen_threshold\t", round(chosen, 4)),
  paste0("method\t", method_used),
  paste0("shape\t", ifelse(is_bimodal, "bimodal", "unimodal")),
  paste0("median_distance\t", round(median(d), 4)),
  paste0("frac_below_", CLONAL_MAX, "\t", round(frac_low, 4)),
  paste0("verdict\t", verdict)
), file.path(base_path, "threshold_decision.txt"))

cat("\nWrote:\n")
cat("  ", file.path(base_path, "thresholds.tsv"), "(", round(chosen, 4),
    "for all", length(all_samples), "samples)\n")
cat("  ", file.path(base_path, "threshold_decision.txt"), "(decision record + verdict)\n")
cat("\nNext: set CLONE_THRESHOLD =", round(chosen, 4),
    "in the clone-calling script and run it.\n")
