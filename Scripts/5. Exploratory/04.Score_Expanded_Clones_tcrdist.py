#!/usr/bin/env python
"""
EXPLORATORY - score expanded TARA clonotypes against the epitope databases
with tcrdist, for comparison with the published Trex call.

NOT FOR THE PAPER, and NOT a replacement for the published method. The
project's specificity call is Trex (Borcherding's R package, annotateDB()),
TRB CDR3 only, Levenshtein edit distance <= 2, against VDJdb / McPAS-TCR /
IEDB, stored as TRB_Epitope.species. Nothing here overwrites that; the Trex
call is carried through the table unchanged and only compared against.

NOTE: Trex is NOT TCRex. Different tools; not conflated anywhere here.

ENGINE
  TCR Explorer (local checkout, installed with the [tcrdist] extra) computes
  the authoritative tcrdist metric via pwseqdist - the same distance engine
  tcrdist3 delegates to - using the vendored tcrdist3 germline CDR table for
  the CDR1/CDR2/CDR2.5 loops. Without the extra it silently falls back to a
  BLOSUM CDR3 distance, so the engine actually used is recorded per row in
  `engine_beta` / `engine_paired` and must read "tcrdist". Any row that falls
  back is visible in the output table, not hidden.

TWO DISTANCES PER CLONOTYPE, AND WHY
  The task asks for paired scoring where both chains are present. We do that,
  but we ALSO score beta-only for every row, because the two numbers answer
  different questions:

    dist_beta    beta chain only. This is the like-for-like comparison with
                 Trex, which only ever looks at TRB CDR3.
    dist_paired  alpha + beta. TCR Explorer defines the paired distance as the
                 SUM of the two single-chain tcrdists, so it lives on roughly
                 twice the scale and is penalised whenever the alpha chain
                 disagrees - which is most of the time, because the reference
                 databases are beta-heavy.

  Reporting only the paired distance would make tcrdist look far more
  conservative than it is, for a reason that has nothing to do with the beta
  chain Trex actually used.

THRESHOLDS (no universal cutoff exists; these are stated, not discovered)
  beta    <= 36 tdus. In the tcrdist3 meta-clonotype framework
          (Mayer-Blackwell et al., eLife 2021) beta-chain radii were optimised
          per centroid and ranged 10-36 tdus, median 22; tcrdist3's own
          documentation uses radius = 36 as the standard bulk search radius.
          36 is therefore the LOOSEST radius that framework endorses - a
          deliberately generous choice, so that "tcrdist does not confirm it"
          cannot be an artefact of an overly strict cutoff.
  paired  <= 72 tdus, i.e. 2 x 36, because the paired distance is the sum of
          two single-chain distances on the same scale.

  Both are reported with a sensitivity sweep so the authors can see exactly
  how the counts move when the threshold is tightened or loosened.

Out: Projects/Analysis_Exploratory/TCRdist/
       tcrdist_scored_expanded_clones.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

OUT = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist")
# The index must come from the persistent copy built by tcr-explorer-refresh
# under OUT; the platformdirs default is not durable in this sandbox.
os.environ.setdefault("TCR_EXPLORER_DATA", str(OUT / "tcr_explorer_data"))

import pandas as pd  # noqa: E402
from tcr_explorer import similarity as S  # noqa: E402
from tcr_explorer.data_paths import records_index_path  # noqa: E402

TOP_K = 5
N_PROC = 12


def _na(x):
    return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() in ("", "NA", "nan")


def score_one(rec: dict) -> dict:
    """Score one clonotype beta-only and (when both chains are present) paired."""
    out = dict(
        row_id=rec["row_id"],
        engine_beta=None, dist_beta=None, n_cand_beta=0,
        beta_epitope=None, beta_antigen=None, beta_organism=None,
        engine_paired=None, dist_paired=None, n_cand_paired=0,
        paired_epitope=None, paired_antigen=None, paired_organism=None,
        warn_beta="", warn_paired="",
    )

    # ---- beta-only (directly comparable to Trex) ----
    if not _na(rec["cdr3b"]) and not _na(rec["vb_gene"]):
        jb = "" if _na(rec["jb_gene"]) else str(rec["jb_gene"])
        n, eng, tot, w = S.find_similar_tcrs(
            str(rec["cdr3b"]), str(rec["vb_gene"]), jb,
            species="human", chain="beta", top_k=TOP_K)
        out["engine_beta"] = eng
        out["n_cand_beta"] = tot
        out["warn_beta"] = ";".join(x.code for x in w)
        if n:
            b = n[0]
            out.update(dist_beta=b.distance, beta_epitope=b.epitope_aa,
                       beta_antigen=b.antigen, beta_organism=b.antigen_organism)

    # ---- paired ----
    if (not _na(rec["cdr3a"]) and not _na(rec["va_gene"])
            and not _na(rec["cdr3b"]) and not _na(rec["vb_gene"])):
        n, eng, tot, w = S.find_similar_paired_tcrs(
            str(rec["cdr3a"]), str(rec["va_gene"]),
            str(rec["cdr3b"]), str(rec["vb_gene"]),
            species="human", top_k=TOP_K)
        out["engine_paired"] = eng
        out["n_cand_paired"] = tot
        out["warn_paired"] = ";".join(x.code for x in w)
        if n:
            p = n[0]
            out.update(dist_paired=p.distance, paired_epitope=p.epitope_aa,
                       paired_antigen=p.antigen, paired_organism=p.antigen_organism)
    return out


def main() -> int:
    idx = records_index_path()
    if not idx.exists():
        sys.exit(f"records index missing at {idx}; run tcr-explorer-refresh first")
    print(f"index: {idx}", flush=True)

    clones = pd.read_csv(OUT / "expanded_clones.csv")
    print(f"expanded clonotypes: {len(clones)}", flush=True)
    clones["row_id"] = range(len(clones))

    recs = clones[["row_id", "cdr3a", "va_gene", "cdr3b", "vb_gene", "jb_gene"]].to_dict("records")

    from multiprocessing import Pool
    with Pool(N_PROC) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(score_one, recs, chunksize=4), 1):
            res.append(r)
            if i % 100 == 0:
                print(f"  scored {i}/{len(recs)}", flush=True)

    scored = clones.merge(pd.DataFrame(res), on="row_id", how="left")
    scored.to_csv(OUT / "tcrdist_scored_expanded_clones.csv", index=False)
    print(f"\nwrote {OUT / 'tcrdist_scored_expanded_clones.csv'}", flush=True)

    # Engine audit: any BLOSUM fallback must be visible, never hidden.
    print("\nengine_beta:\n", scored.engine_beta.value_counts(dropna=False), flush=True)
    print("\nengine_paired:\n", scored.engine_paired.value_counts(dropna=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
