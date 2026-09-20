################################################################################
# Replace the defunct `slot =` argument with `layer =`
#
# WHAT THIS IS AND WHY IT IS SEPARATE
#   This is NOT part of the source-data export. It fixes a pre-existing
#   incompatibility between the figure scripts and the installed SeuratObject:
#
#     Error: The `slot` argument of `GetAssayData()` was deprecated in
#            SeuratObject 5.0.0 and is now defunct.
#
#   It is kept in its own script because it changes analysis code rather than
#   adding an export, so it should be reviewed and run deliberately.
#
# WHAT IT CHANGES
#   `slot = ` -> `layer = ` inside GetAssayData(), AverageExpression() and
#   AggregateExpression() calls only. Both replacements were verified against
#   this repository's own objects on SeuratObject 5.4.0:
#     GetAssayData(o, layer = "data")                      -> 137 x 18843
#     AverageExpression(o, ..., layer = "data")             -> 5 x 10
#   `slot` is defunct for both, so this is a required change, not a preference.
#
# WHAT IT DOES NOT TOUCH
#   Any other use of the word slot, and any `slot =` inside a comment.
#
#   Rscript "00d.Fix_Deprecated_Seurat_Args.R"          # apply
#   Rscript "00d.Fix_Deprecated_Seurat_Args.R" --check  # report only
################################################################################
self_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (i in rev(seq_len(sys.nframe()))) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(dirname(normalizePath(o)))
  }
  getwd()
}
here <- self_dir()
CHECK <- "--check" %in% commandArgs(trailingOnly = TRUE)

targets <- c("Figure_1.R", "Figure_2.R", "Figure_3.R", "Figure_4.R", "Figure_5.R",
             "Supliment_4_9.R")

total <- 0L
for (f in targets) {
  p <- file.path(here, f)
  if (!file.exists(p)) next
  txt <- readLines(p, warn = FALSE)
  is_comment <- grepl("^\\s*#", txt)

  # Only lines that are part of a call to one of these functions. The call can
  # span lines, so a line carrying `slot =` counts when the nearest preceding
  # opening call is one of them.
  FUNS <- "GetAssayData|AverageExpression|AggregateExpression"
  opens <- grepl(FUNS, txt)
  owner <- rep(NA_integer_, length(txt))
  cur <- NA_integer_
  depth <- 0L
  for (i in seq_along(txt)) {
    if (opens[i]) cur <- i
    owner[i] <- cur
    # a line ending a top-level call closes ownership
    o_ <- lengths(regmatches(txt[i], gregexpr("\\(", txt[i])))
    c_ <- lengths(regmatches(txt[i], gregexpr("\\)", txt[i])))
    depth <- depth + o_ - c_
    if (depth <= 0L) { depth <- 0L; if (!opens[i]) cur <- NA_integer_ }
  }

  hits <- which(grepl("\\bslot\\s*=", txt) & !is_comment & !is.na(owner))
  if (!length(hits)) { cat(sprintf("  %-18s nothing to change\n", f)); next }

  if (CHECK) {
    cat(sprintf("  %-18s %d call(s) would change:\n", f, length(hits)))
    for (h in hits) cat(sprintf("      %4d: %s\n", h, trimws(txt[h])))
    total <- total + length(hits)
    next
  }

  bak <- paste0(p, ".preSeurat5.bak")
  if (!file.exists(bak)) file.copy(p, bak)
  txt[hits] <- gsub("\\bslot(\\s*)=", "layer\\1=", txt[hits])
  writeLines(txt, p)

  left <- sum(grepl("\\bslot\\s*=", readLines(p, warn = FALSE)) &
              !grepl("^\\s*#", readLines(p, warn = FALSE)))
  cat(sprintf("  %-18s %d changed, %d slot= left in code\n", f, length(hits), left))
  total <- total + length(hits)
}

if (CHECK) {
  cat(sprintf("\n%d call(s) would change. Re-run without --check to apply.\n", total))
} else {
  cat(sprintf("\n%d call(s) changed. Backups written as *.preSeurat5.bak\n", total))
  cat("Parse check:\n")
  bad <- 0L
  for (f in targets) {
    p <- file.path(here, f)
    if (!file.exists(p)) next
    ok <- tryCatch({ parse(p); TRUE }, error = function(e) { cat("  FAIL", f, conditionMessage(e), "\n"); FALSE })
    if (ok) cat(sprintf("  ok    %s\n", f)) else bad <- bad + 1L
  }
  if (bad) quit(status = 1)
}
