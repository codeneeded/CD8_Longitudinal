#!/usr/bin/env python
"""
EXPLORATORY - CoNGA TCR clumping on the TARA CD8 repertoire.

NOT FOR THE PAPER. Trex remains the project's published specificity call and
is untouched. Trex is NOT TCRex.

WHY CLUMPING RATHER THAN A DISTANCE CUTOFF
  Script 06 built a tcrdist network at fixed radii (60 / 92 / 120) and found
  it almost entirely disconnected. A fixed radius is an arbitrary instrument:
  it cannot tell "these two TCRs are close" from "TCRs of this V/J composition
  are close to everything". CoNGA's clumping test is calibrated instead - for
  each clonotype it compares the OBSERVED number of neighbours within a radius
  against the number expected for a TCR of the same V/J usage and CDR3 length,
  from a background of `num_random_samples` rearrangements. A clonotype is a
  clumping hit when it has significantly MORE neighbours than chance allows.
  That is the better test of convergent selection, so a negative result from
  the raw network is not the last word.

INPUT - the EXISTING CD8 clones file, reused rather than re-exported
  Analysis_Exploratory/CoNGA/tara_cd8_clones.tsv (9,003 clonotypes), written
  by Scripts/5. Exploratory/01.Export_for_CoNGA.R with all three parsing traps
  already handled. The background model needs a repertoire, not just the
  expanded tail, so running clumping on all 9,003 CD8 clonotypes and then
  intersecting the hits with the 980 expanded clonotypes gives the test more
  power than running it on the expanded set alone.

  NOTE ON COMPARTMENT: this file is the CD8 compartment, while the expanded
  set in scripts 03-06 comes from the whole-PBMC object. The intersection is
  therefore over expanded clonotypes that are present in the CD8 object; the
  count is reported rather than assumed.

CoNGA is used from the checkout (sys.path, no pip install) and the
tcrdist_cpp binaries are already compiled.

Out: Projects/Analysis_Exploratory/TCRdist/
       tables/conga_tcr_clumping_hits.csv
       tables/conga_clumping_vs_annotation.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

OUT = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist")
CONGA_DIR = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/CoNGA")
(OUT / "tables").mkdir(parents=True, exist_ok=True)
sys.path.insert(0, "/home/ammasakshay/Desktop/Projects/Packages/conga-master")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

import conga  # noqa: E402
import conga.tcr_clumping  # noqa: E402
import conga.util  # noqa: E402

RADII = [24, 48, 72, 96]     # CoNGA's own defaults, on CoNGA's tcrdist scale
NUM_RANDOM = 50000
PVALUE = 1.0                 # keep everything; filter afterwards


def strip_allele(g: str) -> str:
    return str(g).split("*")[0]


def main() -> int:
    if not conga.util.tcrdist_cpp_available():
        sys.exit("tcrdist_cpp binaries not available")

    clones = pd.read_csv(CONGA_DIR / "tara_cd8_clones.tsv", sep="\t")
    print(f"CD8 clonotypes: {len(clones)}", flush=True)

    tcrs = [((r.va_gene, r.ja_gene, r.cdr3a, str(r.cdr3a_nucseq).lower()),
             (r.vb_gene, r.jb_gene, r.cdr3b, str(r.cdr3b_nucseq).lower()))
            for r in clones.itertuples(index=False)]

    res = conga.tcr_clumping.find_tcr_clumping(
        tcrs, "human",
        tmpfile_prefix=str(OUT / "tmp_clumping"),
        radii=RADII, num_random_samples=NUM_RANDOM, pvalue_threshold=PVALUE,
        verbose=False)
    print(f"\nclumping rows returned: {len(res)}", flush=True)
    if len(res) == 0:
        print("no clumping hits at any radius", flush=True)
        res.to_csv(OUT / "tables" / "conga_tcr_clumping_hits.csv", index=False)
        return 0

    res = res.merge(clones[["clone_id", "subject", "clone_size"]].reset_index()
                    .rename(columns={"index": "clone_index"}),
                    on="clone_index", how="left")
    res.to_csv(OUT / "tables" / "conga_tcr_clumping_hits.csv", index=False)

    sig = res[res.pvalue_adj < 0.05]
    print(f"significant (adj p < 0.05): {len(sig)} rows, "
          f"{sig.clone_index.nunique()} distinct clonotypes", flush=True)
    if len(sig):
        print("\nby radius:", flush=True)
        print(sig.groupby("nbr_radius").clone_index.nunique().to_string(), flush=True)
        if "clumping_group" in sig.columns:
            gs = sig.groupby("clumping_group").clone_index.nunique().sort_values(ascending=False)
            print(f"\nclumping groups: {len(gs)}; sizes: {gs.head(10).tolist()}", flush=True)

    # ---- intersect with the expanded set and its Trex annotation ----
    exp = pd.read_csv(OUT / "tcrdist_group_min_distances.csv")
    key = lambda d, a, b, c, e, f, g: (  # noqa: E731
        d[a].map(strip_allele) + "|" + d[b].map(strip_allele) + "|" + d[c].astype(str) + "|"
        + d[e].map(strip_allele) + "|" + d[f].map(strip_allele) + "|" + d[g].astype(str))
    clones["rkey"] = key(clones, "va_gene", "ja_gene", "cdr3a", "vb_gene", "jb_gene", "cdr3b")
    expp = exp[exp.chain_mode.eq("paired")].copy()
    expp["rkey"] = key(expp, "va_gene", "ja_gene", "cdr3a", "vb_gene", "jb_gene", "cdr3b")
    expp = expp.drop_duplicates("rkey")

    hit_keys = set(clones.loc[clones.index.isin(sig.clone_index.unique()), "rkey"]) if len(sig) else set()
    inter = expp[expp.rkey.isin(set(clones.rkey))].copy()
    inter["clumping_hit"] = inter.rkey.isin(hit_keys)
    print(f"\nexpanded paired clonotypes found in the CD8 clones file: "
          f"{len(inter)} of {len(expp)}", flush=True)
    print(f"  of these, CoNGA clumping hits: {int(inter.clumping_hit.sum())}", flush=True)

    rows = []
    if len(inter) and inter.clumping_hit.any():
        a = int((inter.clumping_hit & inter.trex_annotated).sum())
        b = int((inter.clumping_hit & ~inter.trex_annotated).sum())
        c = int((~inter.clumping_hit & inter.trex_annotated).sum())
        d = int((~inter.clumping_hit & ~inter.trex_annotated).sum())
        orr, p = fisher_exact([[a, b], [c, d]])
        print(f"\n  clumping hits: {a} Trex-annotated / {b} unannotated", flush=True)
        print(f"  non-hits:      {c} Trex-annotated / {d} unannotated", flush=True)
        print(f"  Fisher OR={orr:.2f}  p={p:.4f}", flush=True)
        rows.append(dict(comparison="clumping_hit_vs_trex_annotated",
                         hit_annotated=a, hit_unannotated=b,
                         nonhit_annotated=c, nonhit_unannotated=d,
                         odds_ratio=round(float(orr), 3), p_value=float(p)))
    pd.DataFrame(rows).to_csv(OUT / "tables" / "conga_clumping_vs_annotation.csv", index=False)
    inter.to_csv(OUT / "tables" / "expanded_clonotypes_clumping_status.csv", index=False)
    print(f"\nwrote {OUT / 'tables'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
