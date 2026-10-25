#!/usr/bin/env python
"""
EXPLORATORY - derive the tcrdist cutoff from database ground truth, instead of
importing a number from another study.

NOT FOR THE PAPER. Trex remains the published specificity call. Trex is NOT TCRex.

WHY DERIVE IT
  tcrdist is a weighted sum over CDR loops, so its scale depends on how many
  chains and loops enter the sum. Published neighbour radii run from ~50
  (CDR3-only) to ~120 (paired); a number lifted from one paper encodes that
  paper's chain count, not a universal boundary. The authors' suggested 92 is
  therefore swept and reported, but NOT treated as given.

METHOD
  The reference databases carry TCRs with a known epitope, which is labelled
  ground truth. For each chain configuration separately:
    1. keep epitopes with >= MIN_PER_EPITOPE TCRs,
    2. compute pairwise tcrdist within that set,
    3. label each pair positive if the two TCRs share an epitope,
    4. ROC over the negated distance; report AUC, the Youden-J cutoff, and the
       cutoffs at 1% and 5% false-positive rate.

  TWO CONFIGURATIONS, TWO CUTOFFS. A cutoff calibrated on paired distances is
  invalid for single-chain scoring, because the paired distance is the sum of
  two chains. Scripts 04/05 score beta-only; script 06 builds the paired
  network. So both are calibrated, each with the SAME engine it is applied
  with:
     beta   -> TCR Explorer / pwseqdist (the engine used in scripts 04-05)
     paired -> CoNGA TcrDistCalculator (the engine used in script 06, and the
               scale CoNGA's clumping radii live on)
  Those two engines differ systematically (script 06 measured r=0.97, CoNGA
  median +12 tdus), which is exactly why each cutoff is learned on its own
  engine rather than converted.

  SUBSAMPLING. Pairwise work is quadratic, so each epitope is capped at
  MAX_PER_EPITOPE TCRs drawn at random with a fixed seed, and the paired
  configuration is additionally capped at MAX_EPITOPES_PAIRED epitopes because
  the CoNGA calculator is a per-pair Python call. The caps are reported in the
  output table; they reduce precision, not validity, since the pairs remain a
  random sample of the labelled set.

  The cutoff is learned on REFERENCE data and then applied to the TARA
  clonotypes. It is not tuned on TARA.

Out: Projects/Analysis_Exploratory/TCRdist/tables/
       threshold_calibration_roc.csv      summary per configuration
       threshold_roc_curve_<cfg>.csv      the ROC curve points
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

OUT = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist")
os.environ.setdefault("TCR_EXPLORER_DATA", str(OUT / "tcr_explorer_data"))
(OUT / "tables").mkdir(parents=True, exist_ok=True)
sys.path.insert(0, "/home/ammasakshay/Desktop/Projects/Packages/tcr-explorer-main/src")
sys.path.insert(0, "/home/ammasakshay/Desktop/Projects/Packages/conga-master")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_curve, roc_auc_score  # noqa: E402
from tcr_explorer import similarity as S  # noqa: E402
from tcr_explorer import tcrdist_engine as te  # noqa: E402
from tcr_explorer.data_paths import records_index_path  # noqa: E402

MIN_PER_EPITOPE = 20
MAX_PER_EPITOPE = 60
MAX_EPITOPES_BETA = 60      # 484 epitopes clear the >=20 bar; all of them would
MAX_EPITOPES_PAIRED = 40    # give a 23.6k x 23.6k matrix (558M pairs)
SEED = 0
# Cutoffs at which positive predictive value is reported. The beta grid spans
# the tcrdist3 meta-clonotype range (10-36) and above; the paired grid spans
# the range the authors raised (60-120) plus the tighter end.
PRECISION_CUTOFFS = {
    "beta":   [0, 6, 12, 18, 24, 30, 36, 48, 60, 74],
    "paired": [0, 12, 24, 36, 48, 60, 72, 92, 120, 150, 200],
}


def resolvable(genes):
    """Keep only V genes the tcrdist germline table can resolve.

    The databases carry ambiguous multi-gene calls such as 'TRBV6-1/6-5/6-6',
    which have no single germline CDR set; tcrdist_one_to_many raises on them
    rather than guessing. They are dropped from the ground-truth set and the
    number dropped is reported, because silently keeping them would either
    crash or, worse, assign one arbitrary allele's loops.
    """
    from tcr_explorer import tcrdist_engine as _te
    return genes.map(lambda g: _te.resolve_v_loops(str(g)) is not None)


def summarise(name, dist, same, caps):
    ok = ~np.isnan(dist)
    dist, same = dist[ok], same[ok]
    y = same.astype(int)
    score = -dist                      # closer = more likely same epitope
    auc = roc_auc_score(y, score)
    fpr, tpr, thr = roc_curve(y, score)
    cut = -thr
    j = tpr - fpr
    youden = float(cut[int(np.argmax(j))])

    def at_fpr(target):
        i = int(np.searchsorted(fpr, target, side="right") - 1)
        i = max(i, 0)
        return float(cut[i]), float(tpr[i]), float(fpr[i])

    c1, t1, f1 = at_fpr(0.01)
    c5, t5, f5 = at_fpr(0.05)

    # PRECISION AT SMALL DISTANCES.
    # AUC answers "does distance rank same-epitope pairs above different-epitope
    # pairs across the whole range", which is the wrong question for a cutoff:
    # most same-epitope TCR pairs are genuinely dissimilar (they use unrelated
    # V/CDR3 solutions to the same peptide), which caps AUC no matter how
    # informative the metric is at the short end. What a threshold actually
    # needs is the POSITIVE PREDICTIVE VALUE at that threshold - if two TCRs
    # are this close, how often do they really share an epitope - and its
    # enrichment over the base rate.
    base = float(y.mean())
    prec_rows = []
    for t in PRECISION_CUTOFFS[name]:
        sel = dist <= t
        n_sel = int(sel.sum())
        ppv = float(y[sel].mean()) if n_sel else float("nan")
        prec_rows.append(dict(cutoff=t, n_pairs=n_sel, ppv=round(ppv, 4),
                              enrichment=round(ppv / base, 1) if n_sel else float("nan"),
                              recall=round(float(y[sel].sum() / y.sum()), 4) if n_sel else 0.0))
    pr = pd.DataFrame(prec_rows)
    print(f"  base rate (any pair shares an epitope): {100*base:.2f}%", flush=True)
    print("  precision at fixed cutoffs:", flush=True)
    print("    " + pr.to_string(index=False).replace("\n", "\n    "), flush=True)
    pr.insert(0, "config", name)
    pr.to_csv(OUT / "tables" / f"threshold_precision_{name}.csv", index=False)
    print(f"\n=== {name} ===", flush=True)
    print(f"  pairs={len(y)}  positives={y.sum()} ({100*y.mean():.2f}%)  AUC={auc:.4f}", flush=True)
    print(f"  Youden-J cutoff      : {youden:.0f}  (TPR={tpr[np.argmax(j)]:.3f}, FPR={fpr[np.argmax(j)]:.3f})",
          flush=True)
    print(f"  cutoff at 1% FPR     : {c1:.0f}  (TPR={t1:.3f})", flush=True)
    print(f"  cutoff at 5% FPR     : {c5:.0f}  (TPR={t5:.3f})", flush=True)
    pd.DataFrame({"cutoff": cut, "tpr": tpr, "fpr": fpr}).to_csv(
        OUT / "tables" / f"threshold_roc_curve_{name}.csv", index=False)
    return dict(config=name, n_pairs=int(len(y)), n_positive=int(y.sum()),
                pct_positive=round(100 * float(y.mean()), 3), auc=round(float(auc), 4),
                youden_cutoff=round(youden), youden_tpr=round(float(tpr[np.argmax(j)]), 3),
                youden_fpr=round(float(fpr[np.argmax(j)]), 3),
                cutoff_fpr1=round(c1), tpr_at_fpr1=round(t1, 3),
                cutoff_fpr5=round(c5), tpr_at_fpr5=round(t5, 3), **caps)


def main() -> int:
    rng = np.random.default_rng(SEED)
    df = pd.read_parquet(records_index_path())
    h = df[df.species.fillna("").str.lower().eq("human")]
    rows = []

    # ---------------- beta-only, TCR Explorer engine ----------------
    b = (h[h.chain.eq("beta")].dropna(subset=["cdr3_aa", "v_gene", "epitope_aa"])
         .drop_duplicates(["cdr3_aa", "v_gene"]))
    n0 = len(b)
    b = b[resolvable(b.v_gene)]
    print(f"beta: dropped {n0 - len(b)} records with an unresolvable V gene", flush=True)
    vc = b.epitope_aa.value_counts()
    keep = vc[vc >= MIN_PER_EPITOPE].index
    print(f"beta: {len(keep)} epitopes clear >= {MIN_PER_EPITOPE} TCRs; "
          f"keeping the {min(len(keep), MAX_EPITOPES_BETA)} most represented", flush=True)
    keep = keep[:MAX_EPITOPES_BETA]
    b = b[b.epitope_aa.isin(keep)]
    b = (b.groupby("epitope_aa", group_keys=False)
         .apply(lambda g: g.sample(min(len(g), MAX_PER_EPITOPE), random_state=SEED)))
    print(f"beta ground truth: {len(b)} TCRs over {b.epitope_aa.nunique()} epitopes", flush=True)

    cd, vg, ep = b.cdr3_aa.tolist(), b.v_gene.tolist(), b.epitope_aa.to_numpy()
    n = len(cd)
    D = np.full((n, n), np.nan)
    for i in range(n):
        d, _ = te.tcrdist_one_to_many(cd[i], vg[i], cd, vg)
        D[i] = [np.nan if x is None else x for x in d]
    iu = np.triu_indices(n, k=1)
    rows.append(summarise("beta", D[iu], ep[iu[0]] == ep[iu[1]],
                          dict(n_tcrs=n, n_epitopes=int(b.epitope_aa.nunique()),
                               min_per_epitope=MIN_PER_EPITOPE,
                               max_per_epitope=MAX_PER_EPITOPE, engine="tcr_explorer_pwseqdist")))

    # ---------------- paired, CoNGA engine ----------------
    from conga.tcrdist.tcr_distances import TcrDistCalculator
    p = S._build_paired_index(h, "human")
    p = p.dropna(subset=["cdr3_a_aa", "v_a_gene", "cdr3_b_aa", "v_b_gene", "epitope_aa"])
    p = p.drop_duplicates(["cdr3_a_aa", "v_a_gene", "cdr3_b_aa", "v_b_gene"])
    n0 = len(p)
    p = p[resolvable(p.v_a_gene) & resolvable(p.v_b_gene)]
    print(f"paired: dropped {n0 - len(p)} records with an unresolvable V gene", flush=True)
    vc = p.epitope_aa.value_counts()
    keep_all = vc[vc >= MIN_PER_EPITOPE].index
    print(f"paired: {len(keep_all)} epitopes clear >= {MIN_PER_EPITOPE} TCRs; "
          f"keeping the {min(len(keep_all), MAX_EPITOPES_PAIRED)} most represented", flush=True)
    keep = keep_all[:MAX_EPITOPES_PAIRED]
    p = p[p.epitope_aa.isin(keep)]
    p = (p.groupby("epitope_aa", group_keys=False)
         .apply(lambda g: g.sample(min(len(g), MAX_PER_EPITOPE), random_state=SEED)))
    print(f"paired ground truth: {len(p)} TCRs over {p.epitope_aa.nunique()} epitopes", flush=True)

    al = lambda g: g if "*" in str(g) else f"{g}*01"  # noqa: E731
    # CoNGA needs a J gene; the paired index does not carry one, and the J term
    # is not part of the tcrdist sum (CDR1/2/2.5 come from V, CDR3 from the
    # junction), so a fixed placeholder is used and recorded here explicitly.
    tcrs = [((al(r.v_a_gene), "TRAJ1*01", r.cdr3_a_aa),
             (al(r.v_b_gene), "TRBJ1-1*01", r.cdr3_b_aa))
            for r in p.itertuples(index=False)]
    tdc = TcrDistCalculator("human")
    ep = p.epitope_aa.to_numpy()
    n = len(tcrs)
    D = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            try:
                D[i, j] = D[j, i] = float(tdc(tcrs[i], tcrs[j]))
            except Exception:
                pass
    iu = np.triu_indices(n, k=1)
    rows.append(summarise("paired", D[iu], ep[iu[0]] == ep[iu[1]],
                          dict(n_tcrs=n, n_epitopes=int(p.epitope_aa.nunique()),
                               min_per_epitope=MIN_PER_EPITOPE,
                               max_per_epitope=MAX_PER_EPITOPE, engine="conga_tcrdist")))

    pd.DataFrame(rows).to_csv(OUT / "tables" / "threshold_calibration_roc.csv", index=False)
    print(f"\nwrote {OUT / 'tables' / 'threshold_calibration_roc.csv'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
