################################################################################
# Patch the figure scripts to export source data              (Editor point E-3)
#
# Inserts one source() line into each figure script, after its OPENING library
# block, and makes base_dir overridable so the scripts run from any checkout.
# Idempotent: running it twice changes nothing. Writes a .bak once per file.
#
#   Rscript "00b.Add_Source_Data_Export.R"
#
# PLACEMENT IS THE WHOLE PROBLEM
#   A first version inserted after the LAST library() call in the file. Figure_2
#   has library(tidyverse) at line 1074 and Figure_4 has library(viridis) at 575,
#   so the shim landed after most panels had already been written: 19 of 19
#   panels in Figure_2 and 7 of 8 in Figure_4 were saved unwrapped, silently.
#   The shim must precede the FIRST ggsave, so it goes after the first
#   contiguous library block and the result is verified, not assumed.
################################################################################
self_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  # Rscript replaces spaces in --file= with "~+~"; this repo has a space in
  # "Scripts/3. Manuscript", so the path is unusable until it is decoded.
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (i in rev(seq_len(sys.nframe()))) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(dirname(normalizePath(o)))
  }
  getwd()
}
here <- self_dir()
if (!file.exists(file.path(here, "00.Source_Data_Shim.R"))) here <- getwd()

targets <- c("Figure_1.R", "Figure_2.R", "Figure_3.R", "Figure_4.R", "Figure_5.R",
             "Supliment_4_9.R")
PREAMBLE <- c(
  "",
  "# Export the data behind every panel (Editor point E-3).",
  "# Resolves this script's own directory under BOTH `Rscript foo.R` and",
  "# `source(\"foo.R\")`. sys.frame(1) alone fails under Rscript: there is no",
  "# frame 1 at top level, which halted every figure script on the first try.",
  ".sd_shim_dir <- local({",
  "  a <- commandArgs(trailingOnly = FALSE)",
  "  f <- sub(\"^--file=\", \"\", a[grep(\"^--file=\", a)])",
  "  f <- gsub(\"~+~\", \" \", f, fixed = TRUE)   # Rscript encodes spaces",
  "  if (length(f)) dirname(normalizePath(f[1])) else {",
  "    d <- NULL",
  "    for (i in rev(seq_len(sys.nframe()))) {",
  "      o <- sys.frame(i)$ofile",
  "      if (!is.null(o)) { d <- dirname(normalizePath(o)); break }",
  "    }",
  "    if (is.null(d)) getwd() else d",
  "  }",
  "})",
  ".sd_shim <- file.path(.sd_shim_dir, \"00.Source_Data_Shim.R\")",
  "if (file.exists(.sd_shim)) source(.sd_shim) else",
  "  warning(\"source-data shim not found at \", .sd_shim,",
  "          \"; panels will not export\", call. = FALSE, immediate. = TRUE)",
  "# Repository root, derived from this script's location. Replaces the",
  "# hardcoded ~/Documents/CD8_Longitudinal, which does not exist here.",
  "CD8_ROOT <- if (exists(\"sd_repo_root\")) sd_repo_root(.sd_shim_dir) else",
  "  Sys.getenv(\"CD8_BASE\", unset = normalizePath(",
  "    file.path(.sd_shim_dir, \"..\", \"..\"), mustWork = FALSE))",
  "")

# End of the first contiguous library block: start at the first library() line
# and walk forward while lines are library/require/comment/blank.
opening_block_end <- function(txt) {
  is_lib <- grepl("^\\s*(library|require|suppressPackageStartupMessages)\\s*\\(", txt)
  first <- which(is_lib)
  if (!length(first)) return(0L)
  i <- first[1]
  last <- i
  while (i < length(txt)) {
    i <- i + 1L
    if (is_lib[i]) { last <- i; next }
    if (grepl("^\\s*(#.*)?$", txt[i])) next
    break
  }
  last
}

fail <- 0L
for (f in targets) {
  p <- file.path(here, f)
  if (!file.exists(p)) { cat(sprintf("  MISSING  %s\n", f)); next }

  # A previous run may have placed the shim wrongly. Restore from .bak so the
  # correct placement is computed against the original file.
  bak <- paste0(p, ".bak")
  if (file.exists(bak) && any(grepl("00.Source_Data_Shim.R", readLines(p, warn = FALSE),
                                    fixed = TRUE))) {
    file.copy(bak, p, overwrite = TRUE)
    cat(sprintf("  restored from .bak   %s\n", f))
  }
  txt <- readLines(p, warn = FALSE)
  if (!file.exists(bak)) file.copy(p, bak)

  at <- opening_block_end(txt)
  if (at == 0L) { cat(sprintf("  NO library() BLOCK  %s - skipped\n", f)); fail <- fail + 1L; next }

  txt <- append(txt, PREAMBLE, after = at)

  # Every occurrence of the hardcoded root becomes CD8_ROOT. Comments are left
  # alone; only code strings are rewritten.
  LIT <- '"~/Documents/CD8_Longitudinal"'
  is_comment <- grepl("^\\s*#", txt)
  nrep <- sum(grepl(LIT, txt, fixed = TRUE) & !is_comment)
  txt[!is_comment] <- gsub(LIT, "CD8_ROOT", txt[!is_comment], fixed = TRUE)

  txt <- c(txt, "", "if (exists(\"sd_report\")) sd_report()")
  writeLines(txt, p)

  # VERIFY: no ggsave may precede the shim, or that panel exports nothing.
  chk    <- readLines(p, warn = FALSE)
  shim_l <- grep("00.Source_Data_Shim.R", chk, fixed = TRUE)[1]
  gg     <- grep("ggsave\\(", chk)
  before <- sum(gg < shim_l)
  ok     <- before == 0L
  left <- sum(grepl('"~/Documents/CD8_Longitudinal"', chk, fixed = TRUE) &
              !grepl("^\\s*#", chk))
  if (left) { ok <- FALSE; cat(sprintf("    %d hardcoded path(s) remain in %s\n", left, f)) }
  if (!ok) fail <- fail + 1L
  cat(sprintf("  %s  %-18s shim@%-5d panels %2d  unwrapped %d%s\n",
              if (ok) "patched " else "**FAILED", f, shim_l, length(gg), before,
              sprintf("  %d path literal(s) -> CD8_ROOT", nrep)))
}

if (fail) {
  cat(sprintf("\n%d script(s) would export incomplete source data. Not safe to run.\n", fail))
  quit(status = 1)
}
cat("\nAll scripts patched with the shim ahead of every panel.\n")
cat("Now run each figure script as usual, then build the workbook:\n")
cat('  Rscript "00c.Build_SourceData_Workbook.R"\n')
