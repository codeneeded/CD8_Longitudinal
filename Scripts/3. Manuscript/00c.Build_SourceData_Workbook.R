################################################################################
# Collect every *_source_data.csv into one workbook           (Editor point E-3)
#
#   Rscript "00c.Build_SourceData_Workbook.R"
#
# Writes TARA_SourceData.xlsx with one sheet per panel plus a README sheet
# mapping sheet names to figure files. Sheet names are capped at 31 characters
# because Excel rejects longer ones, and collisions are numbered.
################################################################################

# Resolve the repository root rather than assume it. The default base_dir in
# the figure scripts is ~/Documents/CD8_Longitudinal, which is not where this
# checkout lives - assuming it made this script fail outright. Try, in order:
# CD8_BASE, the repo two levels above this script, then the old default.
# sys.frame(1) does not exist under `Rscript`, so resolve from --file= first
# and fall back to the frame stack for source().
script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  # Rscript replaces spaces in --file= with "~+~"; this repo has a space in
  # "Scripts/3. Manuscript", so the path is unusable until it is decoded.
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) dirname(normalizePath(f[1])) else {
    d <- NULL
    for (i in rev(seq_len(sys.nframe()))) {
      o <- sys.frame(i)$ofile
      if (!is.null(o)) { d <- dirname(normalizePath(o)); break }
    }
    if (is.null(d)) getwd() else d
  }
})
candidates <- unique(c(
  Sys.getenv("CD8_BASE", unset = NA),
  normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE),
  path.expand("~/Documents/CD8_Longitudinal")))
candidates <- candidates[!is.na(candidates)]

base_dir <- NA
for (cand in candidates) {
  if (any(dir.exists(file.path(cand, c("Manuscript", "Analysis"))))) {
    base_dir <- cand; break
  }
}
if (is.na(base_dir))
  stop("no Manuscript/ or Analysis/ found. Tried:\n  ",
       paste(candidates, collapse = "\n  "),
       "\nSet CD8_BASE to the repository root.")
cat("repository root:", base_dir, "\n")

# E-3 asks for the numbers behind every MAIN and SUPPLEMENTARY panel, so the
# workbook is scoped to the manuscript figures. Reviewer_Comments/Figures holds
# revision-response panels, which are separate deliverables already shipped
# with their own CSVs; include them only on request.
INCLUDE_REVISION <- nzchar(Sys.getenv("INCLUDE_REVISION_FIGS"))
want <- c("Manuscript", "Analysis")
if (INCLUDE_REVISION) want <- c(want, file.path("Reviewer_Comments", "Figures"))
roots <- file.path(base_dir, want)
roots <- roots[dir.exists(roots)]

files <- unlist(lapply(roots, function(r)
  list.files(r, pattern = "_source_data\\.csv$", recursive = TRUE,
             full.names = TRUE)))
# _superseded holds retired figures; their data must never enter the workbook.
files <- files[!grepl("(^|/)_superseded(/|$)", files)]
files <- sort(unique(files))

# Report per root, so "0 under Manuscript/" reads as "the figure scripts have
# not been run yet" rather than as a silent near-empty workbook.
for (rt_ in roots) {
  n <- sum(startsWith(files, rt_))
  cat(sprintf("  %-28s %3d source-data files\n",
              sub(paste0(base_dir, "/"), "", rt_, fixed = TRUE), n))
}
cat(sprintf("total: %d\n", length(files)))
if (!length(files))
  stop("no source-data files found. Run 00b.Add_Source_Data_Export.R, then each\n",
       "  figure script - the panels export as they are saved.")
if (!any(startsWith(files, file.path(base_dir, "Manuscript"))))
  warning("no manuscript panels have exported yet: run the Figure_*.R scripts.",
          call. = FALSE, immediate. = TRUE)

# Checked here, not at the top, so the root and the file count are reported
# even when the package is missing - that is the information needed to decide
# whether anything is wrong beyond a missing dependency.
if (!requireNamespace("openxlsx", quietly = TRUE))
  stop("openxlsx is required to write the workbook.\n",
       "  install.packages(\"openxlsx\")\n",
       "The per-panel CSVs above are already written and complete without it.")
library(openxlsx)

wbk <- createWorkbook()
used <- character(0)
idx <- data.frame(sheet = character(0), figure = character(0),
                  rows = integer(0), cols = integer(0), path = character(0),
                  stringsAsFactors = FALSE)

for (f in files) {
  d <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) { cat("  skip (empty):", basename(f), "\n"); next }
  nm <- sub("_source_data$", "", tools::file_path_sans_ext(basename(f)))
  nm <- gsub("[\\[\\]\\*/\\\\?:]", "_", nm)
  nm <- substr(nm, 1, 31)
  if (nm %in% used) {
    k <- 2
    while (paste0(substr(nm, 1, 28), "_", k) %in% used) k <- k + 1
    nm <- paste0(substr(nm, 1, 28), "_", k)
  }
  used <- c(used, nm)
  addWorksheet(wbk, nm)
  writeData(wbk, nm, d)
  idx <- rbind(idx, data.frame(sheet = nm, figure = basename(f),
    rows = nrow(d), cols = ncol(d),
    path = sub(path.expand(base_dir), "", f, fixed = TRUE),
    stringsAsFactors = FALSE))
}

addWorksheet(wbk, "README")
writeData(wbk, "README", data.frame(note = c(
  "Source data for the TARA CD8 longitudinal manuscript.",
  "One sheet per figure panel. The figure column gives the panel the sheet belongs to;",
  "path gives its location in the repository (github.com/codeneeded/CD8_Longitudinal).",
  "Generated by Scripts/3. Manuscript/00c.Build_SourceData_Workbook.R",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")))))
writeData(wbk, "README", idx, startRow = 8)
worksheetOrder(wbk) <- c(length(names(wbk)), seq_len(length(names(wbk)) - 1))

out <- file.path(path.expand(base_dir), "TARA_SourceData.xlsx")
saveWorkbook(wbk, out, overwrite = TRUE)
cat(sprintf("\nwrote %s\n  %d panel sheets + README\n", out, nrow(idx)))
print(utils::head(idx[, c("sheet", "rows", "cols")], 10))
