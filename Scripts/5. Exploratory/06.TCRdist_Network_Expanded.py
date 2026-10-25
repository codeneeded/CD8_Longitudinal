#!/usr/bin/env python
"""
EXPLORATORY - tcrdist similarity network over the expanded TARA clonotypes.

NOT FOR THE PAPER. Trex (annotateDB(), TRB CDR3 only, Levenshtein ED <= 2)
remains the project's published specificity call and is untouched here; it is
only used as a node colouring. Trex is NOT TCRex.

THE QUESTION THIS SERVES (Reviewer 2, major points 3e and 5b)
  The response currently CONCEDES that 61-70% of expanded clonotypes are
  unannotated because public clonotypes dominate the databases. If the
  unannotated clonotypes sit in structured tcrdist communities, and those
  communities also contain annotated members, then the unannotated fraction
  reflects database COVERAGE rather than absence of shared specificity - a
  tested statement rather than a concession. If they are scattered with no
  community structure, that is equally a real answer and is reported as such.

NODES
  Expanded clonotypes (>= 2 cells in a sample) WITH BOTH CHAINS, because the
  edge metric is the paired tcrdist. Of the 980 expanded sample x clonotype
  pairs, 647 are paired; 296 are beta-only and 37 carry no usable beta. The
  333 non-paired clonotypes are therefore absent from this graph, and every
  count below is out of the paired set, not 980. A beta-only network would not
  be comparable because the paired distance is the SUM of two chains and so
  sits on roughly twice the scale.

  ONE NODE PER UNIQUE CLONOTYPE, not per sample x clonotype. The 647 paired
  rows collapse to 559 distinct CTstrict: 83 clonotypes are expanded at more
  than one timepoint and so appear in several rows. Left uncollapsed they
  generate edges at distance 0 to themselves - on the first run 93 of the 242
  edges at R=92 (38%) were one clonotype meeting itself at another timepoint.
  Those are longitudinal persistence, not convergent specificity, and counting
  them as network structure would manufacture exactly the result this analysis
  is supposed to test. This also follows the project convention that
  specificity is a property of the receptor, not of the cell.

EDGES
  Paired tcrdist below a threshold. The authors suggest < 92. That is NOT
  hardcoded as the only value: the graph is rebuilt at 60, 92 and 120 and the
  edge count, component count, largest-component size and community count are
  reported at each, so the choice can be made on this dataset's own distance
  distribution rather than on a number carried in from another study.

COMMUNITIES
  Leiden (leidenalg, RBConfiguration). CoNGA's own calibrated clumping test is
  run separately in script 08 - that asks whether a clonotype has MORE
  neighbours than chance allows, which is a better instrument than any fixed
  radius. The two are complementary and both are reported.

Out: Projects/Analysis_Exploratory/TCRdist/
       network_nodes.csv, network_edges_<R>.csv
       tables/network_threshold_sweep.csv
       tables/network_community_composition.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

OUT = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist")
os.environ.setdefault("TCR_EXPLORER_DATA", str(OUT / "tcr_explorer_data"))
(OUT / "tables").mkdir(parents=True, exist_ok=True)
sys.path.insert(0, "/home/ammasakshay/Desktop/Projects/Packages/tcr-explorer-main/src")
# CoNGA is used from the checkout, not pip-installed (house convention).
sys.path.insert(0, "/home/ammasakshay/Desktop/Projects/Packages/conga-master")

import igraph as ig  # noqa: E402
import leidenalg  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from tcr_explorer import tcrdist_engine as te  # noqa: E402

THRESHOLDS = [60, 92, 120]
SEED = 0
MIN_COMMUNITY = 5   # below this we report the size, never name it as a finding


def _allele(g):
    g = str(g)
    return g if "*" in g else g + "*01"


def build_matrix_conga(nodes: pd.DataFrame) -> np.ndarray:
    """Paired tcrdist on CoNGA's scale.

    WHY CoNGA'S METRIC AND NOT TCR EXPLORER'S HERE. Both implement "tcrdist"
    but they do not return the same number: over all 155,961 clonotype pairs
    they correlate r=0.97 with CoNGA systematically HIGHER (median +12 tdus,
    max difference 79).

    The divergence is in CDR3 LENGTH handling, not in substitutions. Measured
    directly, holding V/J and the alpha chain fixed:
        CASSLSMNTEAFF vs CASSLSMNTQAFF  (one substitution)  CoNGA  6, Explorer  6
        CASSLSMNTEAFF vs CASSLSMNTEAFFF (one residue longer) CoNGA 72, Explorer 48
        CASSLSMNTEAFF vs CASSLSMNTEAFFFF (two longer)        CoNGA 84, Explorer 60
    Substitutions score identically; a CDR3 length difference costs about 24
    tdus more under CoNGA. (CoNGA's tcrdist.hh applies
    gap_penalty_cdr3_region_ = 12 per length unit and explicitly does NOT
    weight it further, while TCR Explorer passes gap_penalty = 4 into
    pwseqdist; the two therefore do not coincide, and the earlier assumption
    in this file that they did was wrong.)

    A threshold is only meaningful on the scale it was derived on. The
    authors' suggested 92, and CoNGA's own clumping radii (24/48/72/96) used
    in script 07, are both CoNGA-scale numbers - so the network is built on
    CoNGA's metric and the two analyses are directly comparable. TCR
    Explorer's metric is retained as a cross-check, not as the edge weight.
    """
    from conga.tcrdist.tcr_distances import TcrDistCalculator
    tdc = TcrDistCalculator("human")
    tcrs = [((_allele(r.va_gene), _allele(r.ja_gene), r.cdr3a),
             (_allele(r.vb_gene), _allele(r.jb_gene), r.cdr3b))
            for r in nodes.itertuples(index=False)]
    n = len(tcrs)
    D = np.full((n, n), np.nan)
    for i in range(n):
        D[i, i] = 0.0
        for j in range(i + 1, n):
            try:
                d = float(tdc(tcrs[i], tcrs[j]))
            except Exception:
                d = np.nan
            D[i, j] = D[j, i] = d
    return D


def build_matrix_explorer(nodes: pd.DataFrame) -> np.ndarray:
    """Same graph on TCR Explorer's (tcrdist3/pwseqdist) scale, for cross-check."""
    tcrs = [te.PairedTCR(a, b, c, d) for a, b, c, d in
            zip(nodes.cdr3a, nodes.va_gene, nodes.cdr3b, nodes.vb_gene)]
    n = len(tcrs)
    D = np.full((n, n), np.nan)
    for i, q in enumerate(tcrs):
        d, _ = te.tcrdist_paired_one_to_many(q, tcrs)
        D[i] = [np.nan if x is None else x for x in d]
    D = np.fmin(D, D.T)
    np.fill_diagonal(D, 0.0)
    return D


def graph_at(D: np.ndarray, R: float) -> ig.Graph:
    iu = np.triu_indices(D.shape[0], k=1)
    keep = (D[iu] < R) & ~np.isnan(D[iu])
    edges = list(zip(iu[0][keep], iu[1][keep]))
    g = ig.Graph(n=D.shape[0], edges=edges)
    g.es["weight"] = [float(D[i, j]) for i, j in edges]
    return g


def communities(g: ig.Graph):
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition, seed=SEED, n_iterations=10)
    return np.array(part.membership), part.modularity


def main() -> int:
    d = pd.read_csv(OUT / "tcrdist_group_min_distances.csv")
    paired = d[d.chain_mode.eq("paired")].reset_index(drop=True).copy()
    n_beta_only = int((d.chain_mode == "beta_only").sum())
    n_unscorable = int((d.chain_mode == "unscorable").sum())
    print(f"expanded sample x clonotype pairs: {len(d)}  paired rows: {len(paired)}", flush=True)
    print(f"  excluded: beta-only {n_beta_only}, unscorable {n_unscorable}", flush=True)

    # one node per unique receptor (see header)
    paired["participant"] = paired["sample"].str.split("_").str[0]
    first = lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan  # noqa: E731
    nodes = (paired.sort_values("n_cells", ascending=False)
             .groupby("CTstrict", as_index=False)
             .agg(n_cells=("n_cells", "sum"),
                  n_samples=("sample", "nunique"),
                  samples=("sample", lambda s: ";".join(sorted(set(s)))),
                  participant=("participant", lambda s: ";".join(sorted(set(s)))),
                  n_participants=("participant", "nunique"),
                  va_gene=("va_gene", first), ja_gene=("ja_gene", first), cdr3a=("cdr3a", first),
                  vb_gene=("vb_gene", first), jb_gene=("jb_gene", first), cdr3b=("cdr3b", first),
                  **{"TRB_Epitope.species": ("TRB_Epitope.species", first)},
                  Annotation=("Annotation", first),
                  CD8_Annotation=("CD8_Annotation", first),
                  trex_annotated=("trex_annotated", "max"),
                  trex_cmv=("trex_cmv", "max"),
                  beta_min_CMV=("beta_min_CMV", "min"),
                  dist_beta=("dist_beta", "min"),
                  dist_paired=("dist_paired", "min"))
             .reset_index(drop=True))
    print(f"  collapsed to unique clonotypes (network nodes): {len(nodes)} "
          f"({int((paired.CTstrict.value_counts() > 1).sum())} expanded at >1 timepoint)", flush=True)

    D = build_matrix_conga(nodes)
    np.save(OUT / "paired_tcrdist_matrix_conga.npy", D)
    Dx = build_matrix_explorer(nodes)
    np.save(OUT / "paired_tcrdist_matrix_explorer.npy", Dx)

    iu = np.triu_indices(len(nodes), k=1)
    off, offx = D[iu], Dx[iu]
    ok = ~np.isnan(off) & ~np.isnan(offx)
    print(f"\nmetric cross-check over {ok.sum()} clonotype pairs: "
          f"r={np.corrcoef(off[ok], offx[ok])[0,1]:.4f}, "
          f"CoNGA minus TCR Explorer: median {np.median(off[ok]-offx[ok]):+.0f}, "
          f"max |diff| {np.abs(off[ok]-offx[ok]).max():.0f}", flush=True)
    print("\npaired tcrdist (CoNGA scale) between expanded clonotypes:", flush=True)
    for q in [0.01, 0.1, 1, 5, 25, 50]:
        print(f"  {q:5.2f}th pct: {np.nanpercentile(off, q):7.1f}", flush=True)
    Dm = D.copy(); np.fill_diagonal(Dm, np.inf)
    nn = np.nanmin(Dm, axis=1)
    for R in THRESHOLDS:
        print(f"  clonotypes with >=1 neighbour < {R}: {(nn < R).sum()}/{len(nn)} "
              f"({100*(nn < R).mean():.1f}%)", flush=True)

    # ---- threshold sweep ----
    rows, memb = [], {}
    for R in THRESHOLDS:
        g = graph_at(D, R)
        m, mod = communities(g)
        memb[R] = m
        comps = g.connected_components()
        sizes = sorted((len(c) for c in comps), reverse=True)
        csizes = pd.Series(m).value_counts()
        rows.append(dict(
            threshold=R, n_nodes=len(nodes), n_edges=g.ecount(),
            density=round(2 * g.ecount() / (len(nodes) * (len(nodes) - 1)), 5),
            n_components=len(comps),
            largest_component=sizes[0],
            pct_in_largest=round(100 * sizes[0] / len(nodes), 1),
            n_singletons=sum(1 for s in sizes if s == 1),
            n_communities=len(csizes),
            n_communities_ge5=int((csizes >= MIN_COMMUNITY).sum()),
            modularity=round(mod, 4)))
        nodes[f"community_{R}"] = m
        print(f"\nR={R}: edges={g.ecount()} components={len(comps)} "
              f"largest={sizes[0]} ({100*sizes[0]/len(nodes):.1f}%) "
              f"communities={len(csizes)} (>={MIN_COMMUNITY}: {(csizes>=MIN_COMMUNITY).sum()}) "
              f"modularity={mod:.3f}", flush=True)

        ei, ej = (np.array([e.tuple for e in g.es]).T if g.ecount()
                  else (np.array([], int), np.array([], int)))
        pa = nodes.participant.to_numpy()
        cross = pa[ei] != pa[ej] if len(ei) else np.array([], bool)
        pd.DataFrame({"source": ei, "target": ej,
                      "dist": [float(D[i, j]) for i, j in zip(ei, ej)],
                      "source_participant": pa[ei] if len(ei) else [],
                      "target_participant": pa[ej] if len(ei) else [],
                      "cross_participant": cross}
                     ).to_csv(OUT / f"network_edges_{R}.csv", index=False)
        rows[-1]["n_edges_cross_participant"] = int(cross.sum())
        rows[-1]["n_edges_within_participant"] = int((~cross).sum())

    sweep = pd.DataFrame(rows)
    sweep.to_csv(OUT / "tables" / "network_threshold_sweep.csv", index=False)

    # ---- annotated vs unannotated across communities ----
    comp_rows = []
    for R in THRESHOLDS:
        nodes["_c"] = nodes[f"community_{R}"]
        big = nodes.groupby("_c").filter(lambda x: len(x) >= MIN_COMMUNITY)
        for c, grp in big.groupby("_c"):
            sp = grp.loc[grp.trex_annotated, "TRB_Epitope.species"]
            # UNION across the community, not the per-clonotype value. Each node
            # carries its own n_participants (almost always 1, since one
            # clonotype usually sits in one infant); taking .max() of those would
            # report 1 for a community whose members come from two different
            # participants - which is precisely the cross-participant
            # convergence this table exists to show, and would contradict the
            # `participants` field built on the same row.
            parts = sorted(set(";".join(grp.participant).split(";")))
            comp_rows.append(dict(
                threshold=R, community=int(c), size=len(grp),
                n_annotated=int(grp.trex_annotated.sum()),
                n_unannotated=int((~grp.trex_annotated).sum()),
                pct_annotated=round(100 * grp.trex_annotated.mean(), 1),
                n_CMV=int(grp.trex_cmv.sum()),
                n_participants=len(parts),
                participants=";".join(parts),
                top_species=";".join(sp.value_counts().head(3).index.astype(str)) if len(sp) else ""))

        # Are annotated and unannotated clonotypes distributed across
        # communities differently from chance?
        #
        # NOT a chi-square. The communities here hold single-digit counts, so
        # the expected cell counts fall below 5 and the chi-square
        # approximation is invalid - on the uncollapsed first run it returned
        # chi2 exactly equal to n with p ~ 1e-6, which is the signature of a
        # perfectly separated table, not of evidence. A label permutation test
        # makes no distributional assumption: the annotated/unannotated labels
        # are shuffled across the clonotypes that sit in communities >= 5 and
        # the observed segregation is compared with its own null.
        #
        # Statistic: mean over communities of the squared deviation of the
        # community's annotated fraction from the overall annotated fraction,
        # weighted by community size (i.e. how purely communities segregate).
        if len(big) and big["_c"].nunique() > 1:
            lab = big.trex_annotated.to_numpy().astype(float)
            grp_id = big["_c"].to_numpy()
            uniq = np.unique(grp_id)
            sizes_arr = np.array([(grp_id == u).sum() for u in uniq], dtype=float)

            def stat(v):
                fr = np.array([v[grp_id == u].mean() for u in uniq])
                return float(np.sum(sizes_arr * (fr - v.mean()) ** 2) / sizes_arr.sum())

            obs = stat(lab)
            rng = np.random.default_rng(SEED)
            null = np.array([stat(rng.permutation(lab)) for _ in range(10000)])
            p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
            print(f"\nR={R}: {len(big)} of {len(nodes)} clonotypes sit in a community of "
                  f">= {MIN_COMMUNITY}; segregation statistic={obs:.3f}, "
                  f"permutation p={p:.4f} (10,000 shuffles)", flush=True)
            sweep.loc[sweep.threshold == R, "segregation_stat"] = round(obs, 4)
            sweep.loc[sweep.threshold == R, "segregation_perm_p"] = p
            sweep.loc[sweep.threshold == R, "n_in_communities_ge5"] = len(big)

    pd.DataFrame(comp_rows).to_csv(OUT / "tables" / "network_community_composition.csv", index=False)
    sweep.to_csv(OUT / "tables" / "network_threshold_sweep.csv", index=False)
    nodes.drop(columns=["_c"]).to_csv(OUT / "network_nodes.csv", index=False)
    print(f"\nwrote {OUT / 'network_nodes.csv'} and tables/", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
