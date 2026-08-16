################################################################################
# Install the R stack for the CD8_Longitudinal pipeline (new machine)
#
# Run in a terminal, NOT RStudio (compilation output is easier to debug):
#     Rscript "Scripts/CP018 - FACTseq/install_R_stack.R"
#
# Takes 30-90 min — most packages compile from source. Safe to re-run: it skips
# anything already installed.
#
# SYSTEM LIBRARIES FIRST. Several packages need these headers or they fail with
# "cannot find -lxml2" / "configuration failed". Install before running this:
#
#   sudo apt update && sudo apt install -y \
#     libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
#     libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
#     libtiff5-dev libjpeg-dev libhdf5-dev libglpk-dev libgsl-dev \
#     libgeos-dev libmagick++-dev cmake
################################################################################

options(Ncpus = max(1, parallel::detectCores() - 1))   # parallel compile
options(repos = c(CRAN = "https://cloud.r-project.org"))

message("R version: ", R.version.string)
message("Library:   ", .libPaths()[1])

need <- function(p) !requireNamespace(p, quietly = TRUE)

# ── 1. Infrastructure ────────────────────────────────────────────────────────
for (p in c("BiocManager", "remotes", "devtools")) {
  if (need(p)) install.packages(p)
}

# Bioconductor release must match the R version; let BiocManager decide.
BiocManager::install(version = BiocManager::version(), ask = FALSE, update = FALSE)

# ── 2. CRAN packages ─────────────────────────────────────────────────────────
cran_pkgs <- c(
  "Seurat", "SeuratObject",
  "qs2",                    # object serialisation used by every script
  "Matrix", "dplyr", "tidyr", "stringr", "readr",
  "ggplot2", "patchwork", "ggrepel", "ggalluvial", "ggpubr", "rstatix",
  "RColorBrewer", "viridis", "clustree", "hdf5r", "R.utils"
)
to_get <- cran_pkgs[vapply(cran_pkgs, need, logical(1))]
if (length(to_get)) {
  message("\n== CRAN: ", paste(to_get, collapse = ", "))
  install.packages(to_get)
}

# ── 3. Bioconductor packages ─────────────────────────────────────────────────
bioc_pkgs <- c(
  "scDblFinder",            # doublet detection (script 1)
  "batchelor",              # FastMNN integration (script 3)
  "scRepertoire",           # TCR analysis (script 4)
  "SingleCellExperiment", "SummarizedExperiment",
  "glmGamPoi",              # speeds up SCTransform if used
  "MAST",                   # differential expression (manuscript scripts)
  "ComplexHeatmap"
)
to_get <- bioc_pkgs[vapply(bioc_pkgs, need, logical(1))]
if (length(to_get)) {
  message("\n== Bioconductor: ", paste(to_get, collapse = ", "))
  BiocManager::install(to_get, ask = FALSE, update = FALSE)
}

# ── 4. GitHub-only packages ──────────────────────────────────────────────────
# These are the ones the sandbox could not install.
gh <- c(
  "SeuratExtend"     = "huayc09/SeuratExtend",       # DimPlot2 / VlnPlot2 / ClusterDistrBar
  "scCustomize"      = "samuel-marsh/scCustomize",   # QC_Plots_* / FeaturePlot_scCustom
  "SeuratWrappers"   = "satijalab/seurat-wrappers",  # RunFastMNN wrapper
  "Azimuth"          = "satijalab/azimuth",          # reference annotation
  "Trex"             = "ncborcherding/Trex"          # used by CP003 script 7
)
for (p in names(gh)) {
  if (need(p)) {
    message("\n== GitHub: ", p, " (", gh[[p]], ")")
    tryCatch(remotes::install_github(gh[[p]], upgrade = "never"),
             error = function(e) message("  FAILED: ", conditionMessage(e)))
  }
}

# Monocle3 — only needed for the TARA trajectory scripts, not the CP018 pipeline.
# Uncomment if you plan to re-run trajectory analysis:
BiocManager::install(c("BiocGenerics","DelayedArray","DelayedMatrixStats",
                        "limma","lme4","S4Vectors","SingleCellExperiment",
                        "SummarizedExperiment","batchelor","HDF5Array",
                        "terra","ggrastr"), ask = FALSE, update = FALSE)
 remotes::install_github("cole-trapnell-lab/monocle3", upgrade = "never")

# ── 5. Verify ────────────────────────────────────────────────────────────────
check <- c("Seurat","SeuratObject","qs2","scRepertoire","scDblFinder","batchelor",
           "clustree","SeuratExtend","scCustomize","SeuratWrappers","Azimuth",
           "ggalluvial","hdf5r","Matrix","dplyr","ggplot2","patchwork","ggrepel")

cat("\n", strrep("=", 60), "\nINSTALLATION SUMMARY\n", strrep("=", 60), "\n", sep = "")
res <- data.frame(package = check, version = NA_character_, ok = FALSE)
for (i in seq_along(check)) {
  v <- tryCatch(as.character(packageVersion(check[i])), error = function(e) NA_character_)
  res$version[i] <- ifelse(is.na(v), "MISSING", v)
  res$ok[i] <- !is.na(v)
}
print(res, row.names = FALSE)

missing <- res$package[!res$ok]
if (length(missing)) {
  cat("\nSTILL MISSING:", paste(missing, collapse = ", "), "\n")
  cat("Usually a system library. Check the compile error above and install the\n",
      "matching -dev package, then re-run this script.\n", sep = "")
} else {
  cat("\nAll packages installed. The CP018 pipeline can run in this R.\n")
}

# Azimuth additionally downloads a ~1 GB PBMC reference on FIRST USE
# (from zenodo/seurat servers). Pre-fetch it now so script 1 does not stall:
if (requireNamespace("Azimuth", quietly = TRUE)) {
  cat("\nAzimuth installed. Its PBMC reference downloads on first RunAzimuth()\n",
      "call (~1 GB) — expect a pause the first time script 1 runs.\n", sep = "")
}
