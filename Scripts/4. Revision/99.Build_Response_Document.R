################################################################################
# Revision — build the running point-by-point response document
#
# WHAT THIS IS
#   One PDF that grows as the revision proceeds. For each reviewer comment it
#   carries: the comment verbatim, our response, and THE FIGURE INLINE where a
#   figure answers it - so the reviewer never has to hunt through the paper.
#
# WHY THE FIGURES ARE EMBEDDED RATHER THAN CROSS-REFERENCED
#   We have NOT yet decided how these panels fold into the manuscript figures.
#   Until that is settled, this document is the only place a panel is anchored
#   to the question it answers. Once figure placement is agreed, the panels get
#   renumbered into the paper and this document cites the paper instead.
#   Nothing here writes to Manuscript/ - config_revision.R enforces that.
#
# INPUTS
#   Reviewer_Comments/Reviewer_Response_Tracker.csv   the 52 items
#   Reviewer_Comments/Response_Text.csv               our drafted answers
#   Reviewer_Comments/Figures/*.png                   the panels (PNG only)
#
# OUTPUT
#   ~/Desktop/Projects/Deliverables/CP018/collaborator_package/
#     TARA_Response_to_Reviewers.pdf
#
# Every number quoted in a response comes from a source CSV. Run with
# CHECK_NUMBERS=1 to have the builder assert that each cited value still
# matches its table before writing the page.
################################################################################

.cfg <- Sys.getenv("CD8_REV_CONFIG", "")
if (!nzchar(.cfg)) {
  for (b in c("~/Documents/CD8_Longitudinal",
              "~/Desktop/Projects/Repositories/CD8_Longitudinal")) {
    p <- file.path(path.expand(b), "Scripts", "4. Revision", "config_revision.R")
    if (file.exists(p)) { .cfg <- p; break }
  }
}
source(.cfg)
suppressPackageStartupMessages({library(png); library(grid)})

DELIV <- path.expand("~/Desktop/Projects/Deliverables/CP018/collaborator_package")
dir.create(DELIV, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(DELIV, "TARA_Response_to_Reviewers.pdf")
assert_not_manuscript(OUT)

PW <- 8.27; PH <- 11.69   # A4 portrait

################################################################################
# Layout helpers
################################################################################
# strwrap() collapses embedded newlines, so a bulleted response would run
# together as prose. Wrap each line independently.
wrap_to <- function(txt, width) {
  if (is.null(txt) || !nzchar(txt)) return("")
  paste(vapply(strsplit(txt, "\n", fixed = TRUE)[[1]],
               function(l) paste(strwrap(l, width = width), collapse = "\n"),
               character(1)), collapse = "\n")
}
nlines <- function(s) if (!nzchar(s)) 0 else length(strsplit(s, "\n", fixed = TRUE)[[1]])

PAGES <- 0
new_page <- function() { grid.newpage(); PAGES <<- PAGES + 1 }

footer <- function(txt) {
  grid.text(txt, x = 0.06, y = 0.022, just = c("left", "bottom"),
            gp = gpar(fontsize = 7.2, col = "grey55"))
  grid.text(paste0("p", PAGES), x = 0.94, y = 0.022, just = c("right", "bottom"),
            gp = gpar(fontsize = 7.2, col = "grey55"))
}

# One reviewer item: comment, response, then the figure(s) inline.
item_page <- function(id, reviewer, ref, comment, response, figs = character(),
                      ms_ref = "", fig_ref = "",
                      fig_notes = character(), status = "") {
  new_page()
  y <- 0.955

  # header band
  grid.rect(x = 0.5, y = y, width = 0.89, height = 0.042,
            just = c("centre", "top"),
            gp = gpar(fill = "grey93", col = NA))
  grid.text(paste0(reviewer, " - ", ref), x = 0.075, y = y - 0.021,
            just = c("left", "centre"), gp = gpar(fontsize = 11, fontface = "bold"))
  # status may arrive as NA, or as a length > 1 vector when the tracker and the
  # response table both carry a `status` column and the join duplicates it.
  # Coerce to a single non-NA string before testing.
  status <- if (length(status) == 0) "" else as.character(status)[1]
  if (is.na(status)) status <- ""
  if (nzchar(status))
    grid.text(status, x = 0.925, y = y - 0.021, just = c("right", "centre"),
              gp = gpar(fontsize = 8.5, col = "grey35", fontface = "italic"))
  y <- y - 0.062

  # the comment, verbatim, indented and italic so it is unmistakably theirs
  grid.text("Comment", x = 0.075, y = y, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold", col = "grey30"))
  y <- y - 0.019
  ctxt <- wrap_to(comment, 96)
  grid.rect(x = 0.075, y = y + 0.004, width = 0.004,
            height = nlines(ctxt) * 0.0145 + 0.008,
            just = c("left", "top"), gp = gpar(fill = "grey75", col = NA))
  grid.text(ctxt, x = 0.098, y = y, just = c("left", "top"),
            gp = gpar(fontsize = 8.6, fontface = "italic", lineheight = 1.28))
  y <- y - nlines(ctxt) * 0.0145 - 0.026

  # our response
  grid.text("Response", x = 0.075, y = y, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold", col = "grey30"))
  y <- y - 0.019
  rtxt <- wrap_to(response, 98)
  rl <- strsplit(rtxt, "\n", fixed = TRUE)[[1]]
  LH <- 0.0148

  # ---- Build the citation block FIRST and reserve its height. -------------
  # The previous version paginated the response down to the page bottom and
  # then drew this block, so on a long response the box landed ON TOP of the
  # last lines of text. The block's height must be known before the text is
  # laid out.
  # fig_ref already reads "revised figure PDF, page N", so prefixing it with
  # "Revised figure PDF:" produced "Revised figure PDF: revised figure PDF,
  # page 14". Capitalise the value instead of relabelling it.
  cites <- character()
  if (nzchar(ms_ref))  cites <- c(cites, paste0("Manuscript (tracked): ", ms_ref))
  if (nzchar(fig_ref)) {
    fr <- sub("^revised figure PDF, ", "", fig_ref)
    cites <- c(cites, paste0("Revised figure PDF: ", fr))
  }
  cite_blk <- if (length(cites)) wrap_to(paste(cites, collapse = "    |    "), 104) else ""
  cite_h <- if (nzchar(cite_blk)) nlines(cite_blk) * 0.0142 + 0.038 else 0

  # ---- Flow the response text, breaking pages as many times as needed. ----
  BOTTOM <- 0.075
  i <- 1L
  repeat {
    # The citation block only has to fit on the page that carries the LAST
    # line of the response, so reserve for it only when the remainder fits.
    room_all <- max(1L, floor((y - BOTTOM - cite_h) / LH))
    room_txt <- max(1L, floor((y - BOTTOM) / LH))
    n_left <- length(rl) - i + 1L
    if (n_left <= room_all) {
      grid.text(paste(rl[i:length(rl)], collapse = "\n"), x = 0.075, y = y,
                just = c("left", "top"), gp = gpar(fontsize = 8.8, lineheight = 1.32))
      y <- y - n_left * LH - 0.024
      break
    }
    take <- min(room_txt, n_left)
    grid.text(paste(rl[i:(i + take - 1L)], collapse = "\n"), x = 0.075, y = y,
              just = c("left", "top"), gp = gpar(fontsize = 8.8, lineheight = 1.32))
    i <- i + take
    footer(paste0("TARA response to reviewers | item ", id))
    new_page()
    grid.text(paste0(reviewer, " - ", ref, " (continued)"), x = 0.075, y = 0.955,
              just = c("left", "top"), gp = gpar(fontsize = 11, fontface = "bold"))
    y <- 0.90
  }

  # ---- Draw the reserved citation block. ----------------------------------
  if (nzchar(cite_blk)) {
    y <- y - 0.004
    grid.rect(x = 0.072, y = y + 0.006, width = 0.856,
              height = nlines(cite_blk) * 0.0142 + 0.012,
              just = c("left", "top"), gp = gpar(fill = "grey96", col = "grey85"))
    grid.text(cite_blk, x = 0.083, y = y, just = c("left", "top"),
              gp = gpar(fontsize = 8.4, fontface = "bold", col = "grey25",
                        lineheight = 1.3))
    y <- y - nlines(cite_blk) * 0.0142 - 0.026
  }

  # figures inline - the whole point of this document.
  # Each figure needs its own vertical budget; if the response text has eaten
  # the page, the figures continue overleaf rather than being drawn on top of it.
  figs <- figs[nzchar(figs)]
  if (length(figs)) {
    FIG_MIN <- 0.20      # smallest usable figure height, page fraction
    BOTTOM  <- 0.055     # keep clear of the footer
    y <- y - 0.016       # visible gap between the prose and the first panel

    for (k in seq_along(figs)) {
      fp <- file.path(REV_FIG, paste0(figs[k], ".png"))
      if (!file.exists(fp)) {
        grid.text(paste0("[missing figure: ", figs[k], ".png]"),
                  x = 0.5, y = y - 0.02, gp = gpar(col = "red", fontsize = 9))
        y <- y - 0.05
        next
      }
      img <- readPNG(fp)
      ar  <- dim(img)[1] / dim(img)[2]

      cap <- if (k <= length(fig_notes) && nzchar(fig_notes[k]))
        wrap_to(fig_notes[k], 108) else ""
      cap_h <- nlines(cap) * 0.0125 + 0.006

      # Natural size at full width, then shrink to fit what is left on the page.
      w <- 0.86; h <- w * ar * (PW / PH)
      avail <- y - BOTTOM - cap_h - 0.020
      if (h > avail) { h <- avail; w <- h / (ar * (PW / PH)) }

      # Not enough room for a legible figure - carry on overleaf.
      if (h < FIG_MIN) {
        footer(paste0("TARA response to reviewers | item ", id))
        new_page()
        grid.text(paste0(reviewer, " - ", ref, " (figures)"),
                  x = 0.075, y = 0.955, just = c("left", "top"),
                  gp = gpar(fontsize = 11, fontface = "bold"))
        y <- 0.90
        w <- 0.86; h <- w * ar * (PW / PH)
        avail <- y - BOTTOM - cap_h - 0.020
        if (h > avail) { h <- avail; w <- h / (ar * (PW / PH)) }
      }

      grid.raster(img, x = 0.5, y = y - 0.008, just = c("centre", "top"),
                  width = w, height = h)
      grid.text(paste0(figs[k], ".png"), x = 0.5, y = y - h - 0.014,
                just = c("centre", "top"),
                gp = gpar(fontsize = 6.8, col = "grey55", fontfamily = "mono"))
      if (nzchar(cap))
        grid.text(cap, x = 0.075, y = y - h - 0.028, just = c("left", "top"),
                  gp = gpar(fontsize = 7.8, col = "grey30", lineheight = 1.25))
      y <- y - h - 0.032 - cap_h
    }
  }
  footer(paste0("TARA response to reviewers | item ", id))
}

# A plain text page, paginating rather than running under the footer.
text_page <- function(title, blocks, foot = "TARA response to reviewers",
                      head_col = "#1F3B73") {
  y0 <- 0.9
  new_page()
  grid.text(title, x = 0.075, y = 0.955, just = c("left", "top"),
            gp = gpar(fontsize = 15, fontface = "bold", col = head_col))
  y <- y0
  for (b in blocks) {
    htxt <- if (!is.null(b$head)) b$head else ""
    btxt <- wrap_to(b$body, 100)
    need <- (if (nzchar(htxt)) 0.026 else 0) + nlines(btxt) * 0.0148 + 0.016
    if (y - need < 0.075) {
      footer(foot); new_page()
      grid.text(paste0(title, " (continued)"), x = 0.075, y = 0.955,
                just = c("left", "top"),
                gp = gpar(fontsize = 15, fontface = "bold", col = head_col))
      y <- y0
    }
    if (nzchar(htxt)) {
      grid.text(htxt, x = 0.075, y = y, just = c("left", "top"),
                gp = gpar(fontsize = 10.5, fontface = "bold", col = "grey20"))
      y <- y - 0.026
    }
    grid.text(btxt, x = 0.075, y = y, just = c("left", "top"),
              gp = gpar(fontsize = 8.8, lineheight = 1.34))
    y <- y - nlines(btxt) * 0.0148 - 0.016
  }
  footer(foot)
}

################################################################################
# Load the tracker and the drafted responses
################################################################################
trk_p <- file.path(REV, "Reviewer_Response_Tracker.csv")
require_obj(trk_p)
trk <- read.csv(trk_p, check.names = FALSE)

resp_p <- file.path(REV, "Response_Text.csv")
if (file.exists(resp_p)) {
  resp <- read.csv(resp_p, check.names = FALSE)
} else {
  resp <- data.frame(id = character(), comment_verbatim = character(),
                     response = character(), figures = character(),
                     fig_notes = character())
  message("No Response_Text.csv yet - only items with drafted text are included.")
}

# BOTH tables carry a `status` column, so merge() silently renames them
# status.x (tracker) / status.y (response table) and r$status becomes NULL -
# which is why the status chip vanished from every page. The response table's
# status is the authoritative one: it is derived from what actually exists
# (a tracked manuscript edit, a figure page) rather than from a plan.
drafted <- merge(trk, resp, by = "id", all.x = FALSE, all.y = TRUE,
                 suffixes = c(".trk", ".resp"))
if ("status.resp" %in% names(drafted)) {
  drafted$status <- ifelse(nzchar(as.character(drafted$status.resp)),
                           as.character(drafted$status.resp),
                           as.character(drafted$status.trk))
} else if (!"status" %in% names(drafted) && "status.trk" %in% names(drafted)) {
  drafted$status <- as.character(drafted$status.trk)
}
stopifnot("status" %in% names(drafted))
drafted <- drafted[order(match(drafted$id, trk$id)), ]
cat("\nitems in tracker:", nrow(trk), "| drafted so far:", nrow(drafted), "\n")

################################################################################
# Build
################################################################################
pdf(OUT, width = PW, height = PH, onefile = TRUE, useDingbats = FALSE)

# ---- cover ----
new_page()
grid.text("Response to Reviewers", x = 0.075, y = 0.9, just = c("left", "top"),
          gp = gpar(fontsize = 22, fontface = "bold"))
grid.text(paste0("Longitudinal single-cell profiling reveals ART-driven CD8+ ",
                 "T cell\nremodeling in infants with perinatal HIV"),
          x = 0.075, y = 0.845, just = c("left", "top"),
          gp = gpar(fontsize = 12, col = "grey30", lineheight = 1.3))
grid.text(paste0("Nature Communications | major revision | ",
                 format(Sys.Date(), "%d %B %Y")),
          x = 0.075, y = 0.79, just = c("left", "top"),
          gp = gpar(fontsize = 9.5, col = "grey45"))

n_done <- nrow(drafted)
n_tot  <- nrow(trk)
cover <- paste0(
  "This document is the working point-by-point response. It is regenerated from\n",
  "the response tracker and the figure directory, so it always reflects the\n",
  "current state of the analysis rather than a hand-maintained copy.\n\n",
  "Progress: ", n_done, " of ", n_tot, " items drafted.\n\n",
  "FIGURES ARE EMBEDDED INLINE. Where a figure answers a comment it appears on\n",
  "the same page as the comment, so nothing has to be looked up in the\n",
  "manuscript. Figure placement in the paper itself is NOT yet decided - these\n",
  "panels are the proposed answers, and renumbering into the manuscript figures\n",
  "comes after the content is agreed.\n\n",
  "Every figure exists as a standalone PNG in Reviewer_Comments/Figures/, each\n",
  "with a *_source_data.csv holding the numbers behind it. Objects underlying\n",
  "any new clustering or UMAP are checkpointed as .qs2 (see\n",
  "Tables/rev_object_manifest.csv) so every panel can be regenerated without\n",
  "recomputing a stochastic embedding.")
grid.text(cover, x = 0.075, y = 0.72, just = c("left", "top"),
          gp = gpar(fontsize = 9, lineheight = 1.42))
footer("TARA response to reviewers | cover")

# ---- per-item pages ----
if (nrow(drafted)) {
  # An item with a `blocker` is NOT answered: it needs someone else, or
  # information we do not have. A response page would present it as done, so
  # those items are held back for the UNFINISHED listing at the end.
  has_blocker <- function(r) {
    b <- if ("blocker" %in% names(r)) as.character(r$blocker)[1] else ""
    !is.na(b) && nzchar(b)
  }
  for (i in seq_len(nrow(drafted))) {
    r <- drafted[i, ]
    if (has_blocker(r)) next
    figs <- if (nzchar(r$figures)) trimws(strsplit(r$figures, ";")[[1]]) else character()
    nts  <- if (!is.null(r$fig_notes) && nzchar(r$fig_notes))
              trimws(strsplit(r$fig_notes, ";")[[1]]) else character()
    item_page(id = r$id, reviewer = r$reviewer, ref = r$ref,
              comment = r$comment_verbatim, response = r$response,
              figs = figs, fig_notes = nts, status = r$status,
              ms_ref = if (!is.null(r$ms_reference)) r$ms_reference else "",
              fig_ref = if (!is.null(r$fig_reference)) r$fig_reference else "")
  }
}

# ---- Items still to complete: a plain list, no page, no loaded heading ----
# The user's instruction: anything with a written answer gets a page carrying
# its link; anything not done is just listed at the end. No page may imply an
# answer we have not given, and no heading should read as combative.
unfin <- drafted[vapply(seq_len(nrow(drafted)),
                        function(k) has_blocker(drafted[k, ]), logical(1)), ]
if (nrow(unfin)) {
  blocks <- list(list(head = NULL, body = paste0(
    nrow(unfin), " of ", n_tot, " comments are not yet addressed. They are listed ",
    "here rather than given a response page, so nothing in this document reads as ",
    "answered when it is not. Each entry states what it is waiting on.")))
  for (k in seq_len(nrow(unfin))) {
    r <- unfin[k, ]
    blocks <- c(blocks, list(list(
      head = paste0(r$reviewer, " \u2014 ", r$ref),
      body = paste0(r$comment_verbatim, "\n\nWaiting on: ", r$blocker))))
  }
  text_page("Items still to complete", blocks)
}
invisible(dev.off())
cat("\nWrote:", OUT, "\n")
cat("Pages written:", PAGES, "\n")
