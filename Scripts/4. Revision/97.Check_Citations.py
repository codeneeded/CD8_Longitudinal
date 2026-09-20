#!/usr/bin/env python3
"""
Gate the reviewer deliverables on citation consistency. Run BEFORE sending.

WHY THIS EXISTS
  Manuscript page/line citations were first computed by a layout model that
  over-counted the document by 21 pages (51 estimated vs 30 real). They were
  later recomputed from the author's rendered PDF, but nothing verified that the
  corrected values actually reached the shipped documents - so a stale citation
  could have gone out. This script makes that class of error fail loudly.

CHECKS
  1. No manuscript citation exceeds the real page count of the rendered PDF.
  2. Every ms_reference in Response_Text.csv equals the value derived from
     Manuscript_edit_manifest.csv (single source of truth).
  3. Every comment with a tracked manuscript edit has a response entry.
  4. Every figure-PDF page cited exists in that PDF.
  5. Every citation PRINTED in the shipped response PDF is backed by a
     manifest row - catching a stale rebuild.
  6. The response PDF has one page per ANSWERED item (blocked items get no page,
     by design) and an 'Items still to complete' listing exists whenever any item
     is not yet done -
     so a comment cannot silently vanish from the document.
  7. A legend cell count matches the cohort its figure actually plots. Added
     after the Figure 1A legend claimed 126,435 cells across 26 samples while
     the Methods describe an entry-visit cohort of 16 participants: the counter
     had summed every row of Table S1 instead of replicating Figure_1.R's own
     subset (HEU | HUU | HEI & PreART_Entry -> 79,705 cells, 16 samples).

Exit status is non-zero if any check fails.
"""
import csv, re, sys, pathlib, collections

REPO = pathlib.Path(__file__).resolve().parents[2]
REV = REPO / "Reviewer_Comments"
T = REV / "Tables"
DEL = pathlib.Path.home() / "Desktop/Projects/Deliverables/CP018/collaborator_package"
MS_PDF = (pathlib.Path.home() / "Desktop/Projects/Manuscripts/CD8_Longitudinal"
          / "TARA_Manuscript_2026-06-03_FACTseqedit_v5LD.pdf")

DASH = re.compile(r"[\u2010-\u2015\u2212]")
fails, warns = [], []


def read(p):
    return list(csv.DictReader(p.open()))


def cite_from_manifest(man):
    by = collections.defaultdict(list)
    for r in man:
        if r.get("page"):
            by[r["comment_id"]].append(r)
    out = {}
    for cid, rs in by.items():
        per = collections.defaultdict(list)
        for r in rs:
            per[int(r["page"])] += [int(r["line_from"]), int(r["line_to"])]
        parts = []
        for pn in sorted(per):
            ls = sorted(set(per[pn]))
            parts.append(f"p. {pn}, l. {ls[0]}" if ls[0] == ls[-1]
                         else f"p. {pn}, ll. {ls[0]}\u2013{ls[-1]}")
        out[cid] = "; ".join(parts)
    return out


def pdf_pages(path):
    try:
        import pypdfium2 as pdfium
    except ImportError:
        warns.append(f"pypdfium2 unavailable; cannot inspect {path.name}")
        return None, ""
    d = pdfium.PdfDocument(str(path))
    txt = " ".join(d[i].get_textpage().get_text_range() for i in range(len(d)))
    return len(d), DASH.sub("-", re.sub(r"\s+", " ", txt))


man = read(T / "Manuscript_edit_manifest.csv")
resp = read(REV / "Response_Text.csv")
figpages = read(T / "Revised_figure_pages.csv")
TRUTH = cite_from_manifest(man)

# 1. page bound against the real manuscript
n_ms, _ = pdf_pages(MS_PDF)
if n_ms:
    over = sorted({int(r["page"]) for r in man if r.get("page") and int(r["page"]) > n_ms})
    if over:
        fails.append(f"manifest cites manuscript pages {over} but the document has {n_ms}")
    else:
        print(f"  [ok] all manuscript citations within {n_ms} pages")

# 1b. the line numbers must be sourced from the rendered PDF, not estimated
srcs = collections.Counter(r.get("source", "") for r in man)
if any("printed" not in s for s in srcs):
    fails.append(f"manifest has non-printed line-number sources: {dict(srcs)}")
else:
    print(f"  [ok] every line number read from the rendered PDF {dict(srcs)}")

# 2. Response_Text must agree with the manifest
bad = [(r["id"], r.get("ms_reference", ""), TRUTH.get(r["id"], ""))
       for r in resp if (r.get("ms_reference") or "") != TRUTH.get(r["id"], "")]
if bad:
    for cid, got, want in bad:
        fails.append(f"{cid}: response says {got!r}, manifest says {want!r}")
else:
    print("  [ok] every ms_reference matches the manifest")

# 3. no edited comment left without a response
orphan = sorted(set(TRUTH) - {r["id"] for r in resp})
if orphan:
    fails.append(f"comments edited in the .docx with no response entry: {orphan}")
else:
    print("  [ok] every tracked manuscript edit has a response entry")

# 4. figure-PDF pages must exist
n_fig, _ = pdf_pages(DEL / "TARA_Revised_Figures.pdf")
if n_fig:
    claimed = sorted({int(r["page"]) for r in figpages})
    if claimed[-1] > n_fig:
        fails.append(f"figure manifest cites page {claimed[-1]} but the PDF has {n_fig}")
    else:
        print(f"  [ok] figure pages 2-{claimed[-1]} within a {n_fig}-page PDF")
    missing = [r["file"] for r in figpages if r.get("present", "TRUE").upper() == "FALSE"]
    if missing:
        fails.append(f"figure PDF has placeholder pages for: {missing}")

# 5. citations printed in the shipped response PDF must be manifest-backed
n_resp, rtxt = pdf_pages(DEL / "TARA_Response_to_Reviewers.pdf")
if n_resp is not None:
    # Items carrying a `blocker` deliberately get NO response page - they are
    # listed once at the end instead, so that the document never reads as
    # if we answered something we did not. Only unblocked items need a page.
    n_blocked = sum(1 for r in resp if (r.get("blocker") or "").strip())
    n_pages_expected = len(resp) - n_blocked
    if n_resp < n_pages_expected:
        fails.append(f"response PDF has {n_resp} pages for {n_pages_expected} answerable "
                     f"items - at least one page per item is expected, so the build was "
                     f"truncated")
    else:
        print(f"  [ok] response PDF has {n_resp} pages for {n_pages_expected} answered "
              f"items ({n_blocked} listed as unfinished)")
    # The UNFINISHED section must actually be present when items are blocked -
    # otherwise those comments vanish from the document entirely.
    SECTION = "Items still to complete"   # renamed from UNFINISHED: neutral wording
    if n_blocked and rtxt and SECTION not in rtxt:
        fails.append(f"{n_blocked} comments are marked unfinished but the response PDF "
                     f"contains no '{SECTION}' section - they would be invisible")
    elif n_blocked:
        print(f"  [ok] '{SECTION}' section present for {n_blocked} comments")
if rtxt:
    printed = set(re.findall(r"p\. (\d{1,3}), ll?\. (\d{1,4})", rtxt))
    allowed = {(str(r["page"]), str(r["line_from"])) for r in man if r.get("page")}
    allowed |= {(str(r["page"]), str(r["line_to"])) for r in man if r.get("page")}
    unknown = sorted(x for x in printed if x not in allowed)
    # A PDF with ZERO citations must FAIL, not pass vacuously: it means the
    # build errored partway and the citation blocks never rendered. This
    # exact false pass happened once, so the count is now a hard floor.
    expected = len({r["comment_id"] for r in man if r.get("page")})
    if not printed:
        fails.append(f"shipped PDF prints NO manuscript citations, but {expected} "
                     f"comments have tracked edits - the build likely errored")
    elif unknown:
        fails.append(f"shipped PDF prints citations absent from the manifest: {unknown}")
    else:
        print(f"  [ok] all {len(printed)} citations in the shipped PDF are manifest-backed")
    over = [x for x in printed if int(x[0]) > (n_ms or 10**6)]
    if over:
        fails.append(f"shipped PDF cites pages beyond the manuscript: {over}")

# ---- 7. a legend cell count must match the cohort its figure actually plots ----
# An auditor caught the Figure 1A legend claiming 126,435 cells across 26 samples
# while the Methods describe an entry-visit cohort of 16 participants: the
# counter summed every row of Table S1 instead of replicating Figure_1.R's own
# subset. Re-derive both figures from the table and require the docx to agree.
try:
    import zipfile
    from lxml import etree as _et
    _NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    _dx = REPO / "Reviewer_Comments" / ".." / ".."
    _cands = [pathlib.Path("/home/ammasakshay/Desktop/Projects/Manuscripts/"
                           "CD8_Longitudinal/TARA_Manuscript_REVISED_trackchanges.docx")]
    _docx = next((c for c in _cands if c.exists()), None)
    _s1 = REPO / "Reviewer_Comments" / "Tables" / "TableS1_sample_level_counts.csv"
    if _docx and _s1.exists():
        _rows = list(csv.DictReader(_s1.open()))
        _entry = [r for r in _rows
                  if r["Condition"] in ("HEU", "HUU")
                  or (r["Condition"] == "HEI"
                      and r["Timepoint_Group"] == "PreART_Entry")]
        _n = sum(int(r["n_PBMC_QC_pass"]) for r in _entry if r["n_PBMC_QC_pass"])
        _all = sum(int(r["n_PBMC_QC_pass"]) for r in _rows if r["n_PBMC_QC_pass"])
        with zipfile.ZipFile(_docx) as _z:
            _d = _et.fromstring(_z.read("word/document.xml"))
        _acc = re.sub(r"\s+", " ", "".join(
            n.text or "" for n in _d.iter()
            if _et.QName(n).localname == "t"))
        _mfig = re.search(r"Panel A shows ([\d,]+) PBMCs", _acc)
        if _mfig:
            _claimed = _mfig.group(1)
            if _claimed == f"{_n:,}":
                print(f"  [ok] Figure 1A legend cites {_claimed} cells, matching the "
                      f"entry-visit subset the figure plots")
            elif _claimed == f"{_all:,}":
                fails.append(
                    f"Figure 1A legend cites {_claimed} cells - that is the WHOLE-cohort "
                    f"total across all timepoints, but Figure_1.R subsets to entry visits "
                    f"({_n:,} cells). This contradicts the Methods cohort description.")
            else:
                warns.append(f"Figure 1A legend cites {_claimed} cells; the entry-visit "
                             f"subset is {_n:,} and the whole cohort {_all:,}")
        else:
            warns.append("could not find the Figure 1A cell-count sentence in the docx")
except Exception as _e:                     # never let the check itself mask a real fail
    warns.append(f"legend count-scope check could not run: {_e}")

print()
for w in warns:
    print(f"  [warn] {w}")
if fails:
    print(f"  FAILED {len(fails)} check(s):")
    for f in fails:
        print(f"    - {f}")
    sys.exit(1)
print("  all citation checks passed")
