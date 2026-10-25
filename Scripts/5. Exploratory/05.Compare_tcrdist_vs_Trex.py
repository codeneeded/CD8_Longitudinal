#!/usr/bin/env python
"""
EXPLORATORY - compare tcrdist specificity scoring with the published Trex call.

NOT FOR THE PAPER, and NOT a replacement for the published method. Trex
(Borcherding's R package, annotateDB(); TRB CDR3 only, Levenshtein edit
distance <= 2, against VDJdb / McPAS-TCR / IEDB) remains the project's
specificity call. Trex is NOT TCRex.

WHY A TARGETED MINIMUM, NOT THE TOP-1 HIT
  Script 04 recorded each clonotype's single nearest database neighbour. That
  is the wrong instrument for the question "do the Trex CMV calls survive",
  because the reference is dominated by SARS-CoV-2 (107k of 310k human beta
  records, 34%). A clonotype can sit 12 tdus from a CMV record and still have
  a SARS-CoV-2 record at 11 tdus win the top-1 slot, which would read as "CMV
  not confirmed" when the CMV evidence is in fact unchanged.

  So this script computes, for every clonotype, the MINIMUM tcrdist to each
  organism group separately. The CMV question is then asked directly - how far
  is this clonotype from the nearest CMV record - instead of being inferred
  from which organism happened to win.

VOCABULARY
  The two methods do not share a label vocabulary. Trex reports
  TRB_Epitope.species, which for McPAS-derived records is a PATHOLOGY
  ("Immunization", "Multiple sclerosis"); the index reports antigen_organism.
  Both sides are mapped onto a coarse set of organism groups, and the
  agreement analysis is restricted to the pathogen groups where both
  vocabularies have an unambiguous equivalent. Trex's non-pathogen labels
  (Immunization, named autoimmune and cancer indications) have no organism
  analogue and are reported separately rather than forced into a match.

  The CMV regex is the one the project already uses in
  Scripts/2. Analysis/5.CMV_vs_Expanding_Clones.R, unchanged:
      CMV|cytomegalo|herpesvirus 5|HHV-?5

Out: Projects/Analysis_Exploratory/TCRdist/
       tcrdist_vs_trex_comparison.csv     per-clonotype, both calls
       tcrdist_group_min_distances.csv    per-clonotype min distance per group
       tables/*.csv                       the comparison tables
"""
from __future__ import annotations

import os
import re
from pathlib import Path

OUT = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist")
os.environ.setdefault("TCR_EXPLORER_DATA", str(OUT / "tcr_explorer_data"))
(OUT / "tables").mkdir(parents=True, exist_ok=True)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from tcr_explorer import similarity as S  # noqa: E402
from tcr_explorer import tcrdist_engine as te  # noqa: E402
from tcr_explorer.data_paths import records_index_path  # noqa: E402

N_PROC = 12

# Coarse organism groups. CMV is the project's own regex, unchanged.
GROUPS: dict[str, str] = {
    "CMV":              r"CMV|cytomegalo|herpesvirus 5|HHV-?5",
    "EBV":              r"EBV|Epstein|herpesvirus 4|HHV-?4",
    "Influenza":        r"influenza",
    "SARS-CoV-2":       r"SARS-CoV-?2|SARS coronavirus 2|severe acute respiratory syndrome.*2",
    "HIV":              r"\bHIV\b|immunodeficiency virus",
    "M.tuberculosis":   r"tuberculosis|M\.? ?tub",
    "HBV":              r"\bHBV\b|hepatitis B",
    "HCV":              r"\bHCV\b|hepatitis C",
    "YFV":              r"\bYFV\b|yellow fever",
    "Self/Human":       r"HomoSapiens|Homo sapiens|SelfAg",
}
# Groups whose label exists in BOTH vocabularies, so agreement is meaningful.
PATHOGEN_GROUPS = ["CMV", "EBV", "Influenza", "SARS-CoV-2", "HIV",
                   "M.tuberculosis", "HBV", "HCV", "YFV"]
CMV_RE = re.compile(GROUPS["CMV"], re.I)


def group_mask(series: pd.Series, pattern: str) -> np.ndarray:
    return series.fillna("").str.contains(pattern, case=False, regex=True).to_numpy()


def trex_groups(label: str) -> set[str]:
    """Map a (possibly multi-valued, ';'-separated) Trex species label to groups."""
    if not isinstance(label, str) or not label.strip():
        return set()
    out = set()
    for part in label.split(";"):
        for g, pat in GROUPS.items():
            if re.search(pat, part, re.I):
                out.add(g)
    return out


# ---------------------------------------------------------------- reference
def build_reference():
    df = pd.read_parquet(records_index_path())
    h = df[df.species.fillna("").str.lower().eq("human")]

    beta = (h[h.chain.eq("beta")]
            .dropna(subset=["cdr3_aa", "v_gene"])
            .loc[:, ["cdr3_aa", "v_gene", "antigen_organism"]])
    bu = beta.drop_duplicates(["cdr3_aa", "v_gene"]).reset_index(drop=True)
    key = beta.cdr3_aa + "\x00" + beta.v_gene
    ukey = bu.cdr3_aa + "\x00" + bu.v_gene
    pos = pd.Series(np.arange(len(bu)), index=ukey.to_numpy())
    idx = pos.reindex(key.to_numpy()).to_numpy()
    bmask = {g: np.zeros(len(bu), bool) for g in GROUPS}
    for g, pat in GROUPS.items():
        m = group_mask(beta.antigen_organism, pat)
        bmask[g][idx[m]] = True

    paired = S._build_paired_index(h, "human")
    paired = paired.dropna(subset=["cdr3_a_aa", "v_a_gene", "cdr3_b_aa", "v_b_gene"])
    pcols = ["cdr3_a_aa", "v_a_gene", "cdr3_b_aa", "v_b_gene"]
    pu = paired.drop_duplicates(pcols).reset_index(drop=True)
    pkey = paired[pcols].agg("\x00".join, axis=1)
    upkey = pu[pcols].agg("\x00".join, axis=1)
    ppos = pd.Series(np.arange(len(pu)), index=upkey.to_numpy())
    pidx = ppos.reindex(pkey.to_numpy()).to_numpy()
    pmask = {g: np.zeros(len(pu), bool) for g in GROUPS}
    for g, pat in GROUPS.items():
        m = group_mask(paired.antigen_organism, pat)
        pmask[g][pidx[m]] = True

    return bu, bmask, pu, pmask


_B = _BM = _P = _PM = None


def _init(bu, bmask, pu, pmask):
    global _B, _BM, _P, _PM
    _B, _BM, _P, _PM = bu, bmask, pu, pmask


def _group_min(dists, masks):
    d = np.array([np.nan if x is None else x for x in dists], dtype="float64")
    out = {}
    for g, m in masks.items():
        sub = d[m]
        sub = sub[~np.isnan(sub)]
        out[g] = float(sub.min()) if sub.size else np.nan
    return out


def score_row(rec):
    res = {"row_id": rec["row_id"]}
    if rec["cdr3b"] and rec["vb_gene"] and te.resolve_v_loops(str(rec["vb_gene"])) is not None:
        d, _ = te.tcrdist_one_to_many(str(rec["cdr3b"]), str(rec["vb_gene"]),
                                      _B.cdr3_aa.tolist(), _B.v_gene.tolist())
        for g, v in _group_min(d, _BM).items():
            res[f"beta_min_{g}"] = v
    if all(rec.get(k) for k in ("cdr3a", "va_gene", "cdr3b", "vb_gene")):
        q = te.PairedTCR(str(rec["cdr3a"]), str(rec["va_gene"]),
                         str(rec["cdr3b"]), str(rec["vb_gene"]))
        cands = [te.PairedTCR(a, b, c, dd) for a, b, c, dd in
                 zip(_P.cdr3_a_aa, _P.v_a_gene, _P.cdr3_b_aa, _P.v_b_gene)]
        d, _ = te.tcrdist_paired_one_to_many(q, cands)
        for g, v in _group_min(d, _PM).items():
            res[f"paired_min_{g}"] = v
    return res


def main() -> int:
    scored = pd.read_csv(OUT / "tcrdist_scored_expanded_clones.csv")
    print(f"clonotypes: {len(scored)}", flush=True)

    bu, bmask, pu, pmask = build_reference()
    print(f"reference: {len(bu)} unique beta, {len(pu)} unique paired", flush=True)
    for g in GROUPS:
        print(f"   {g:16s} beta={bmask[g].sum():6d}  paired={pmask[g].sum():5d}", flush=True)

    recs = scored[["row_id", "cdr3a", "va_gene", "cdr3b", "vb_gene"]].where(
        pd.notna(scored[["row_id", "cdr3a", "va_gene", "cdr3b", "vb_gene"]]), None
    ).to_dict("records")

    from multiprocessing import Pool
    with Pool(N_PROC, initializer=_init, initargs=(bu, bmask, pu, pmask)) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(score_row, recs, chunksize=4), 1):
            res.append(r)
            if i % 100 == 0:
                print(f"  {i}/{len(recs)}", flush=True)

    mins = pd.DataFrame(res)
    out = scored.merge(mins, on="row_id", how="left")
    out.to_csv(OUT / "tcrdist_group_min_distances.csv", index=False)
    print(f"wrote {OUT / 'tcrdist_group_min_distances.csv'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
