################################################################################
# Source-data export shim                                    (Editor point E-3)
#
# WHAT THIS DOES
#   Overrides ggsave() in the calling script's environment so that every panel
#   saved also writes the data behind it to <figure_name>_source_data.csv, in a
#   source_data/ folder beside the figure. Nothing else about the figures
#   changes - the same file, at the same size, with the same content.
#
# WHY IT IS A WRAPPER AND NOT 71 EDITS
#   The figure scripts contain 71 ggsave() calls and 16 write.csv() calls, so
#   most panels had no source table. Editing each call site is 71 chances to
#   introduce a typo into working figure code. Wrapping the one function that
#   every panel already goes through gets complete coverage with one line added
#   per script.
#
# HOW TO USE
#   Add this immediately after the library() block:
#       source(file.path(dirname(sys.frame(1)$ofile), "00.Source_Data_Shim.R"))
#   or simply:
#       source("00.Source_Data_Shim.R")
#   The patcher (00b.Add_Source_Data_Export.R) inserts it for you.
#
# WHAT IT CANNOT CAPTURE
#   pheatmap() and grid-drawn output do not pass through ggsave and carry no
#   tidy data frame. Those panels are listed in the manifest with
#   captured = FALSE so the gap is visible rather than silent.
################################################################################

# Resolve the repository root from the script's own location.
#
# WHY THIS IS HERE
#   The figure scripts hardcoded "~/Documents/CD8_Longitudinal". That path does
#   not exist in this checkout, so every script died at its first qs_read. The
#   repo root is two levels above Scripts/3. Manuscript, so it can be derived
#   rather than assumed. CD8_BASE still wins if set.
sd_repo_root <- function(script_dir = getwd()) {
  cands <- c(Sys.getenv("CD8_BASE", unset = NA),
             normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE),
             path.expand("~/Documents/CD8_Longitudinal"))
  cands <- unique(cands[!is.na(cands) & nzchar(cands)])
  for (c_ in cands)
    if (dir.exists(file.path(c_, "saved_R_data")) ||
        dir.exists(file.path(c_, "Manuscript")))
      return(c_)
  stop("cannot locate the CD8_Longitudinal repository root. Tried:\n  ",
       paste(cands, collapse = "\n  "),
       "\nSet CD8_BASE to the repository root and re-run.", call. = FALSE)
}

.sd_dir <- function(fig_path) {
  d <- file.path(dirname(fig_path), "source_data")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.sd_manifest <- new.env(parent = emptyenv())
.sd_manifest$rows <- list()

# Pull every distinct data frame a plot draws from.
#
# PATCHWORK NESTS, AND THE FIRST VERSION ONLY WENT ONE LEVEL DEEP
#   Figure 1B is (a | b) / (c | d). Each patchwork node carries the LAST plot
#   added as its own $data and the earlier ones in $patches$plots, recursively;
#   the outermost node of a "/" composition is a plot_filler with no data at
#   all. A single-level walk therefore found nothing and reported "no tidy
#   data" for a panel whose four sub-plots were all reachable. This recurses,
#   visiting children before self so sheets come out in draw order.
.sd_max_depth <- 8L

.sd_extract <- function(p) {
  acc <- list()
  n <- 0L
  single <- !inherits(p, "patchwork")
  walk <- function(node, depth) {
    if (depth > .sd_max_depth) return(invisible(NULL))
    kids <- tryCatch(node$patches$plots, error = function(e) NULL)
    if (length(kids)) for (k in kids) walk(k, depth + 1L)
    got <- tryCatch(.sd_extract_one(node), error = function(e) list())
    if (length(got)) {
      n <<- n + 1L
      for (nm in names(got)) {
        key <- if (single) nm else paste0("panel", n, "_", nm)
        acc[[key]] <<- got[[nm]]
      }
    }
  }
  walk(p, 0L)
  acc
}

# Layers below this many rows are treated as annotation - cluster labels, scale
# bars, arrows - and not written. Figure 1A produced three such files per panel
# (7 label positions, a 1-row segment, a 1-row point), which would be 200+ junk
# sheets across 71 panels. Plot-level data is ALWAYS kept, however small.
.sd_min_layer_rows <- function() getOption("sd_min_layer_rows", 10L)

.sd_extract_one <- function(p) {
  out <- list()
  base <- tryCatch(p$data, error = function(e) NULL)
  if (is.data.frame(base) && nrow(base)) out[["data"]] <- base
  layers <- tryCatch(p$layers, error = function(e) list())
  for (i in seq_along(layers)) {
    ld <- tryCatch(layers[[i]]$data, error = function(e) NULL)
    # waiver() means "inherit the plot's data" - already captured above
    if (!is.data.frame(ld) || !nrow(ld)) next
    if (nrow(ld) < .sd_min_layer_rows()) next
    dup <- FALSE
    for (have in out) if (identical(dim(have), dim(ld))) dup <- TRUE
    if (!dup) out[[paste0("layer", i)]] <- ld
  }
  out
}

# Drop list-columns: a Seurat/ggplot data frame can carry them and write.csv
# fails on them, which would take the figure down with it.
.sd_flatten <- function(df) {
  keep <- vapply(df, function(col) is.atomic(col) || is.factor(col), logical(1))
  df <- df[, keep, drop = FALSE]
  if (ncol(df) == 0) return(NULL)
  # Column names inherit the figure's axis labels, which are deliberately split
  # over two lines ("Naive\nIntermediate", "TEMRA/\nCTL"). An embedded newline
  # in a spreadsheet header is unreadable, so flatten it to a space.
  names(df) <- gsub("[\r\n]+", " ", names(df))
  df
}

ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  res <- ggplot2::ggsave(filename = filename, plot = plot, ...)
  if (is.character(filename) && length(filename) == 1L)
    .sd_last_dir$path <- dirname(filename)
    .sd_flush_pending(.sd_last_dir$path)

  captured <- FALSE
  n_rows <- 0L

  # Clear this panel's previous exports BEFORE writing new ones. Without this a
  # rerun leaves orphans from older code: the first version of this shim wrote
  # three annotation-layer files per panel, and after those were suppressed the
  # stale files stayed on disk and were still collected into the workbook -
  # mixing data from two different versions of the same figure.
  local({
    stem <- tools::file_path_sans_ext(basename(filename))
    d <- file.path(dirname(filename), "source_data")
    if (dir.exists(d)) {
      old <- list.files(d, pattern = paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])",
                        "\\\\\\1", stem), "(_.*)?_source_data\\.csv$"),
                        full.names = TRUE)
      if (length(old)) file.remove(old)
    }
  })

  tabs <- tryCatch(.sd_extract(plot), error = function(e) list())
  if (length(tabs)) {
    stem <- tools::file_path_sans_ext(basename(filename))
    dir <- .sd_dir(filename)
    for (nm in names(tabs)) {
      df <- .sd_flatten(tabs[[nm]])
      if (is.null(df)) next
      suffix <- if (nm == "data") "" else paste0("_", nm)
      out <- file.path(dir, paste0(stem, suffix, "_source_data.csv"))
      tryCatch({
        utils::write.csv(df, out, row.names = FALSE)
        captured <- TRUE
        n_rows <- n_rows + nrow(df)
      }, error = function(e) message("  [source data] ", basename(out), ": ", e$message))
    }
  }

  .sd_manifest$rows[[length(.sd_manifest$rows) + 1L]] <- data.frame(
    figure = basename(filename), captured = captured, rows = n_rows,
    stringsAsFactors = FALSE)

  if (!captured) message("  [source data] no tidy data for ", basename(filename))
  invisible(res)
}

################################################################################
# pheatmap panels
#
# The heatmap panels never reach ggsave: they are drawn with
# png() -> grid.draw(p$gtable) -> dev.off(), so the ggsave wrapper reports
# nothing for them and their matrices were the one real gap in coverage.
#
# pheatmap() takes the matrix as its first argument, so wrapping it captures
# that matrix directly. Naming is the only difficulty: in Figure_1.R the
# pheatmap() calls run BEFORE the png() that draws them, so there is no open
# device to name the file after. The output directory of the most recent panel
# written in the same script is used instead, which is the figure's own
# directory because scripts write panels in order.
################################################################################
.sd_last_dir <- new.env(parent = emptyenv())
.sd_last_dir$path <- NULL

# Matrices seen before any output directory is known. S4A in Supliment_4_9.R is
# the FIRST panel in its script and is a pheatmap, so nothing had been saved yet
# and the write had nowhere to go - it was reported "not captured". Hold it here
# and flush when the next png()/ggsave() reveals the directory.
.sd_last_dir$pending <- list()

.sd_write_matrix <- function(mat, stem, dir) {
  df <- tryCatch(as.data.frame(as.matrix(mat)), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(list(ok = FALSE, rows = 0L))
  df <- cbind(feature = rownames(df), df)
  names(df) <- gsub("[\r\n]+", " ", names(df))
  d <- file.path(dir, "source_data")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(d, paste0(stem, "_source_data.csv"))
  tryCatch({ utils::write.csv(df, out, row.names = FALSE)
             list(ok = TRUE, rows = nrow(df)) },
           error = function(e) { message("  [source data] ", basename(out), ": ",
                                         e$message)
                                 list(ok = FALSE, rows = 0L) })
}

.sd_flush_pending <- function(dir) {
  if (!length(.sd_last_dir$pending)) return(invisible(NULL))
  if (is.null(dir) || !dir.exists(dir)) return(invisible(NULL))
  for (p in .sd_last_dir$pending) {
    r <- .sd_write_matrix(p$mat, p$stem, dir)
    .sd_manifest$rows[[length(.sd_manifest$rows) + 1L]] <- data.frame(
      figure = paste0(p$stem, " (pheatmap)"), captured = r$ok, rows = r$rows,
      stringsAsFactors = FALSE)
    if (r$ok) message("  [source data] ", p$stem, " -> ", basename(dir))
  }
  .sd_last_dir$pending <- list()
}

png <- function(filename = "Rplot%03d.png", ...) {
  if (is.character(filename) && length(filename) == 1L)
    .sd_last_dir$path <- dirname(filename)
    .sd_flush_pending(.sd_last_dir$path)
  grDevices::png(filename = filename, ...)
}

pheatmap <- function(mat, ...) {
  res <- pheatmap::pheatmap(mat, ...)

  # Supliment_4_9.R passes the matrix inline as as.matrix(comp_mat), so a raw
  # deparse yields "as.matrix(comp_mat)". Take the innermost identifier instead.
  nm <- tryCatch(deparse(substitute(mat))[1], error = function(e) "matrix")
  inner <- regmatches(nm, gregexpr("[A-Za-z][A-Za-z0-9._]*", nm))[[1]]
  inner <- setdiff(inner, c("as.matrix", "as", "matrix", "data.frame", "t", "log", "log2",
                            "log10", "scale", "apply", "unname", "as.data.frame"))
  nm <- if (length(inner)) tail(inner, 1) else "matrix"
  nm <- gsub("[^A-Za-z0-9._-]", "_", nm)
  dots <- list(...)
  # pheatmap's own filename= wins when the caller supplied one
  dir <- if (!is.null(dots$filename) && is.character(dots$filename))
           dirname(dots$filename) else .sd_last_dir$path
  stem <- if (!is.null(dots$filename) && is.character(dots$filename))
            tools::file_path_sans_ext(basename(dots$filename)) else
            paste0("heatmap_", nm)

  if (is.null(dir) || !dir.exists(dir)) {
    # No output directory yet - defer rather than lose the matrix.
    .sd_last_dir$pending[[length(.sd_last_dir$pending) + 1L]] <-
      list(mat = mat, stem = stem)
    message("  [source data] heatmap ", nm,
            " deferred until the first panel is saved")
    return(invisible(res))
  }
  r  <- .sd_write_matrix(mat, stem, dir)
  ok <- r$ok
  nr <- r$rows

  .sd_manifest$rows[[length(.sd_manifest$rows) + 1L]] <- data.frame(
    figure = paste0(stem, " (pheatmap)"), captured = ok, rows = nr,
    stringsAsFactors = FALSE)
  if (!ok) message("  [source data] heatmap matrix not written: ", nm)
  invisible(res)
}

# Call at the end of a figure script to see coverage for that script.
sd_report <- function() {
  if (!length(.sd_manifest$rows)) { cat("no panels saved\n"); return(invisible(NULL)) }
  m <- do.call(rbind, .sd_manifest$rows)
  cat(sprintf("\nsource data: %d of %d panels captured\n",
              sum(m$captured), nrow(m)))
  miss <- m$figure[!m$captured]
  if (length(miss)) cat("  not captured:", paste(miss, collapse = ", "), "\n")
  invisible(m)
}

message("  [source data] ggsave wrapped; panels export to source_data/ beside each figure")
