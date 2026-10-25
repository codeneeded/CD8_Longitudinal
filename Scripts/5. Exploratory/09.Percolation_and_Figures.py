#!/usr/bin/env python
"""
EXPLORATORY - percolation guard + figures for the tcrdist / Trex comparison.

NOT FOR THE PAPER. Trex remains the published specificity call. Trex is NOT TCRex.

PERCOLATION GUARD
  Any distance cutoff has a value above which the graph collapses into one
  giant component; "communities" found above it are an artefact of
  over-connection rather than structure. The threshold is swept finely and the
  largest-component fraction recorded, so the chosen cutoff can be shown to
  sit BELOW that point. If the ROC-derived cutoff sits above it, the conflict
  is reported rather than resolved by quietly preferring one.

FIGURES
  Fig 1  tcrdist vs Trex: where a beta-chain threshold would sit, and whether
         the Trex CMV calls survive independent scoring.
  Fig 2  how the cutoff was chosen: ROC on database ground truth + percolation.
  Fig 3  the network at the calibrated cutoff and at roughly +/- 30 around it,
         with 92 included so the authors can see what it would have given.
  Fig 4  the network at the chosen cutoff under three node colourings.

Out: Projects/Analysis_Exploratory/TCRdist/*.png + tables/percolation_sweep.csv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

OUT = Path("/home/ammasakshay/Desktop/Projects/Analysis_Exploratory/TCRdist")
os.environ.setdefault("TCR_EXPLORER_DATA", str(OUT / "tcr_explorer_data"))
sys.path.insert(0, "/home/ammasakshay/Desktop/Projects/Packages/conga-master")

import igraph as ig  # noqa: E402
import leidenalg  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

mpl.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "savefig.facecolor": "white"})
SEED = 0
MIN_COMMUNITY = 5
CMV_COL = "#C44E52"
GREY = "#BBBBBB"


def graph_at(D, R):
    iu = np.triu_indices(D.shape[0], k=1)
    m = (D[iu] < R) & ~np.isnan(D[iu])
    g = ig.Graph(n=D.shape[0], edges=list(zip(iu[0][m], iu[1][m])))
    return g


# ---------------------------------------------------------------- percolation
def percolation(D, nodes):
    rows = []
    for R in range(20, 401, 4):
        g = graph_at(D, R)
        comps = g.connected_components()
        sizes = sorted((len(c) for c in comps), reverse=True)
        part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                                        seed=SEED, n_iterations=5)
        cs = pd.Series(part.membership).value_counts()
        rows.append(dict(threshold=R, n_edges=g.ecount(), n_components=len(comps),
                         largest_component=sizes[0],
                         largest_frac=sizes[0] / D.shape[0],
                         n_communities_ge5=int((cs >= MIN_COMMUNITY).sum())))
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "tables" / "percolation_sweep.csv", index=False)
    giant = s.loc[s.largest_frac >= 0.5, "threshold"]
    perc_50 = int(giant.iloc[0]) if len(giant) else None
    jump = int(s.loc[s.largest_frac.diff().idxmax(), "threshold"])
    return s, perc_50, jump


# ---------------------------------------------------------------- figures
def fig1(d, beta_marks):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    ax = axes[0]
    ann = d.loc[d.trex_annotated & d.dist_beta.notna(), "dist_beta"]
    una = d.loc[~d.trex_annotated & d.dist_beta.notna(), "dist_beta"]
    bins = np.arange(0, 105, 6)
    ax.hist([ann, una], bins=bins, stacked=False, label=["Trex-annotated", "Trex-unannotated"],
            color=["#4C72B0", GREY], edgecolor="white", linewidth=.4)
    # Marks are the PPV measured on the labelled reference set in script 08 -
    # NOT Youden's J or an FPR-derived cutoff, both of which land far outside
    # this axis (beta 127 and 74 respectively) because AUC is only 0.56.
    top = ax.get_ylim()[1]
    for k, (x, ppv) in enumerate(beta_marks):
        ax.axvline(x, color="k", ls="--", lw=.8)
        ax.text(x + 1, top * (0.97 - 0.1 * k), f"{x:g} tdus\nPPV {ppv:.2f}",
                fontsize=5.5, va="top", ha="left")
    ax.set_xlabel("tcrdist to nearest database record (beta chain, tdus)")
    ax.set_ylabel("expanded clonotypes")
    ax.set_title("Distance to the nearest annotated receptor")
    ax.legend(frameon=False, loc="upper right")
    ax.margins(0.04)

    ax = axes[1]
    radii = [6, 12, 18, 24, 30, 36]
    cmv = d[d.trex_cmv & d.beta_min_CMV.notna()]
    oth = d[~d.trex_cmv & d.beta_min_CMV.notna()]
    yc = [100 * (cmv.beta_min_CMV <= R).mean() for R in radii]
    yo = [100 * (oth.beta_min_CMV <= R).mean() for R in radii]
    ax.plot(radii, yc, "-o", color=CMV_COL, ms=4, label=f"Trex CMV calls (n={len(cmv)})")
    ax.plot(radii, yo, "-o", color=GREY, ms=4, label=f"all other clonotypes (n={len(oth)})")
    ax.set_xlabel("search radius around the clonotype (tdus)")
    ax.set_ylabel("% with a CMV record within radius")
    ax.set_title("Trex CMV calls against independent tcrdist scoring")
    ax.legend(frameon=False, loc="upper left")
    ax.set_ylim(-4, 108)
    ax.margins(0.04)
    fig.tight_layout()
    fig.savefig(OUT / "Fig1_tcrdist_vs_Trex.png")
    return fig


def fig2(roc, sweep, perc_50, jump, cuts, you, f1):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    ax = axes[0]
    for cfg, col in [("beta", "#4C72B0"), ("paired", "#55A868")]:
        f = OUT / "tables" / f"threshold_roc_curve_{cfg}.csv"
        if not f.exists():
            continue
        c = pd.read_csv(f)
        a = roc.loc[roc.config == cfg, "auc"].iloc[0]
        ax.plot(c.fpr, c.tpr, color=col, lw=1.2, label=f"{cfg} (AUC {a:.3f})")
    ax.plot([0, 1], [0, 1], color="k", lw=.6, ls=":")
    ax.set_xlabel("false-positive rate (different epitope, called same)")
    ax.set_ylabel("true-positive rate")
    ax.set_title("Threshold calibrated on database ground truth")
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1]
    ax.plot(sweep.threshold, 100 * sweep.largest_frac, color="#333333", lw=1.2)
    cp = cuts["paired"]
    # The precision-calibrated cutoff came out at exactly the value the authors
    # suggested, so they are ONE line here, not two overlaid at the same x.
    lab = (f"calibrated cutoff = suggested ({cp:.0f})" if abs(cp - 92) < 1e-9
           else f"calibrated cutoff ({cp:.0f})")
    ax.axvline(cp, color="#55A868", ls="--", lw=1.2, label=lab)
    if abs(cp - 92) >= 1e-9:
        ax.axvline(92, color="#8172B2", ls=":", lw=1, label="suggested 92")
    for v, c, l in [(you["paired"], "#8172B2", f"Youden J ({you['paired']:.0f})"),
                    (f1["paired"], "#CCB974", f"1% FPR ({f1['paired']:.0f})")]:
        ax.axvline(v, color=c, ls=":", lw=1, label=l)
    if perc_50:
        ax.axvline(perc_50, color=CMV_COL, lw=1.2,
                   label=f"giant component forms ({perc_50})")
    ax.set_xlabel("paired tcrdist edge threshold (tdus)")
    ax.set_ylabel("% of clonotypes in largest component")
    ax.set_title("Percolation: where the graph collapses")
    ax.legend(frameon=False, loc="upper left", fontsize=5.5)
    fig.tight_layout()
    fig.savefig(OUT / "Fig2_threshold_calibration.png")
    return fig


def _layout(D, R, n):
    g = graph_at(D, R)
    lay = np.array(g.layout_fruchterman_reingold(niter=800, seed=None).coords)
    return g, lay


def _draw(ax, g, lay, colors, sizes=9):
    for e in g.es:
        i, j = e.tuple
        ax.plot([lay[i, 0], lay[j, 0]], [lay[i, 1], lay[j, 1]],
                color="#999999", lw=.4, zorder=1)
    ax.scatter(lay[:, 0], lay[:, 1], c=colors, s=sizes, linewidths=.2,
               edgecolors="white", zorder=2)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def fig3(D, nodes, thresholds):
    cols = np.where(nodes.trex_cmv, CMV_COL,
                    np.where(nodes.trex_annotated, "#4C72B0", GREY))
    fig, axes = plt.subplots(1, len(thresholds), figsize=(7.2, 2.2))
    for ax, R in zip(axes, thresholds):
        g, lay = _layout(D, R, len(nodes))
        _draw(ax, g, lay, cols, sizes=7)
        comps = g.connected_components()
        big = max(len(c) for c in comps)
        ax.set_title(f"threshold {R}\n{g.ecount()} edges, largest group {big}")
    h = [plt.Line2D([], [], marker="o", ls="", color=CMV_COL, label="Trex CMV"),
         plt.Line2D([], [], marker="o", ls="", color="#4C72B0", label="Trex other"),
         plt.Line2D([], [], marker="o", ls="", color=GREY, label="unannotated")]
    fig.legend(handles=h, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(.5, -.06))
    fig.suptitle("Paired tcrdist network of expanded clonotypes across edge thresholds", y=1.04)
    fig.tight_layout()
    fig.savefig(OUT / "Fig3_network_by_threshold.png")
    return fig


def fig4(D, nodes, R):
    g, lay = _layout(D, R, len(nodes))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))

    # Legends sit BELOW their panel, not inside it: these graphs fill the whole
    # canvas with nodes, so an inset legend covers data.
    def under(ax, handles, ncol):
        ax.legend(handles=handles, frameon=False, fontsize=5.5, ncol=ncol,
                  loc="upper center", bbox_to_anchor=(0.5, -0.02),
                  handletextpad=.3, columnspacing=.8)

    c1 = np.where(nodes.trex_cmv, CMV_COL,
                  np.where(nodes.trex_annotated, "#4C72B0", GREY))
    _draw(axes[0], g, lay, c1)
    axes[0].set_title("Trex epitope call")
    h1 = [plt.Line2D([], [], marker="o", ls="", color=CMV_COL, label="CMV"),
          plt.Line2D([], [], marker="o", ls="", color="#4C72B0", label="other species"),
          plt.Line2D([], [], marker="o", ls="", color=GREY, label="unannotated")]
    under(axes[0], h1, 3)

    ann = nodes.CD8_Annotation.fillna("not in CD8 object")
    top = ann.value_counts().index[:6]
    pal = dict(zip(top, ["#4C72B0", "#55A868", GREY, "#8172B2", "#CCB974", "#64B5CD"]))
    c2 = [pal.get(a, GREY) for a in ann]
    _draw(axes[1], g, lay, c2)
    axes[1].set_title("CD8 subset")
    h2 = [plt.Line2D([], [], marker="o", ls="", color=pal[t], label=str(t)) for t in top]
    under(axes[1], h2, 2)

    pid = nodes.participant.str.split(";").str[0]
    tp = pid.value_counts().index[:8]
    pal2 = dict(zip(tp, plt.cm.tab10.colors))
    c3 = [pal2.get(p, GREY) for p in pid]
    _draw(axes[2], g, lay, c3)
    axes[2].set_title("Participant")
    h3 = [plt.Line2D([], [], marker="o", ls="", color=pal2[t], label=str(t)) for t in tp]
    under(axes[2], h3, 4)

    fig.suptitle(f"Paired tcrdist network of expanded clonotypes at threshold {R}", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "Fig4_network_colourings.png")
    return fig


def main() -> int:
    d = pd.read_csv(OUT / "tcrdist_group_min_distances.csv")
    nodes = pd.read_csv(OUT / "network_nodes.csv")
    D = np.load(OUT / "paired_tcrdist_matrix_conga.npy")
    roc = pd.read_csv(OUT / "tables" / "threshold_calibration_roc.csv")

    # CHOOSING THE CUTOFF BY PRECISION, NOT BY YOUDEN'S J.
    # The ROC AUC is poor (0.56 both configurations) and Youden's J therefore
    # lands at an absurdly permissive value (beta 127, paired 272 - around the
    # median of the whole distance distribution). That is not a failure of the
    # metric but of the question AUC asks: most pairs that share an epitope are
    # genuinely dissimilar, so no cutoff can have high recall, and a criterion
    # that rewards balanced accuracy is dragged outwards. For an epitope claim
    # the operating characteristic that matters is positive predictive value -
    # given two receptors this close, how often do they really share an
    # epitope. The cutoff used here is therefore the LARGEST one whose PPV on
    # the labelled reference set is >= PPV_TARGET.
    PPV_TARGET = 0.90
    cuts, ppv_at = {}, {}
    for cfg in ("beta", "paired"):
        pr = pd.read_csv(OUT / "tables" / f"threshold_precision_{cfg}.csv")
        ok = pr[(pr.ppv >= PPV_TARGET) & (pr.n_pairs > 0)]
        cuts[cfg] = float(ok.cutoff.max()) if len(ok) else float(pr.cutoff.min())
        ppv_at[cfg] = float(pr.loc[pr.cutoff == cuts[cfg], "ppv"].iloc[0])
    you = {r.config: float(r.youden_cutoff) for r in roc.itertuples(index=False)}
    f1 = {r.config: float(r.cutoff_fpr1) for r in roc.itertuples(index=False)}
    auc = {r.config: float(r.auc) for r in roc.itertuples(index=False)}
    print(f"AUC: {auc}", flush=True)
    print(f"Youden-J cutoffs (NOT used): {you}", flush=True)
    print(f"1% FPR cutoffs:              {f1}", flush=True)
    print(f"PPV>={PPV_TARGET} cutoffs (used): {cuts}  with PPV {ppv_at}", flush=True)

    sweep, perc_50, jump = percolation(D, nodes)
    print(f"\npercolation: giant component (>=50% of nodes) first at threshold "
          f"{perc_50}; steepest rise at {jump}", flush=True)
    chosen = cuts["paired"]
    print(f"calibrated paired cutoff {chosen:.0f} is "
          f"{'BELOW' if perc_50 and chosen < perc_50 else 'ABOVE'} the percolation point",
          flush=True)

    prb = pd.read_csv(OUT / "tables" / "threshold_precision_beta.csv")
    beta_marks = [(float(r.cutoff), float(r.ppv)) for r in prb.itertuples(index=False)
                  if r.cutoff in (6, 24, 36)]
    fig1(d, beta_marks)
    fig2(roc, sweep, perc_50, jump, cuts, you, f1)
    thresholds = sorted({int(max(20, chosen - 30)), int(chosen), 92, int(chosen + 30)})
    print("network thresholds plotted:", thresholds, flush=True)
    fig3(D, nodes, thresholds)
    fig4(D, nodes, int(chosen))
    print(f"\nwrote figures to {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
