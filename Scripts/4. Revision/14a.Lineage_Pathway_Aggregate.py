#!/usr/bin/env python3
"""
Revision - aggregate the EXISTING per-lineage pathway enrichment into one table.
                                                    (Reviewer 1, major point 2)

WHAT THIS IS AND IS NOT
    This runs no new statistics. Differential expression (MAST, with latent
    covariates) and Enrichr pathway enrichment were already run for every
    annotated population across five contrasts; the outputs are 146 enrichment workbooks
    under Analysis/Pathway_Analysis/ (alongside 278 plot files) and 148
    per-population CSVs under Analysis/Differential_Expression/. Nothing was analysed here that was not
    already analysed. This script collects those results into one table so the
    non-CD8 lineages can be read at a glance.

TWO FILTERS, AND WHY EACH IS NEEDED
    1. Ribosomal and translation terms (rRNA Processing, Peptide Chain
       Elongation, NMD, Cap-dependent Translation Initiation, Influenza
       Infection and the like) reach significance in 26 of 30 populations.
       A term significant almost everywhere does not distinguish a lineage, so
       terms significant in more than half the populations are dropped. The
       cut is recorded in the output, not applied silently.
    2. Ranking genes by fold change is not used anywhere here. In these tables
       the largest fold changes are the same handful of low-expression features
       in every population (ENSG00000257322, ENSG00000280441, ANK2, RBFOX2,
       BOLA2), which is a depth artefact rather than shared biology. Enrichment
       uses the significance-filtered gene sets, and the driver genes written
       out are the ones each term actually matched.

OUT
    Reviewer_Comments/Tables/R1M2_07_lineage_pathway_summary.csv    - per term
    Reviewer_Comments/Tables/R1M2_08_lineage_pathway_drivers.csv    - genes
    Reviewer_Comments/Tables/R1M2_09_dropped_ubiquitous_terms.csv   - the cut
"""
import csv
import pathlib
import collections

import pandas as pd

REPO = pathlib.Path("/home/ammasakshay/Desktop/Projects/Repositories/CD8_Longitudinal")
PA = REPO / "Analysis" / "Pathway_Analysis"
DE = REPO / "Analysis" / "Differential_Expression"
OUT = REPO / "Reviewer_Comments" / "Tables"

# Interpretable collections only. TF and PWM sheets are omitted: a single
# overlapping gene can carry a TF term here, which is not a lineage claim.
SHEETS = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023",
          "Reactome_2022", "KEGG_2021_Human"]

FAMILY = [("CD8", ["_CD8", "CD8_"]),
          ("CD4", ["_CD4", "CD4dim", "Treg", "Th2_Th17"]),
          ("NK", ["NK"]),
          ("B cell", ["B_cell"]),
          ("Monocyte", ["Monocyte"]),
          ("pDC", ["pDC"]),
          ("gd T", ["gd_T", "Vd1", "Vg9"]),
          ("MAIT", ["MAIT"])]


def family_of(lineage):
    for name, keys in FAMILY:
        if any(k in lineage for k in keys):
            return name
    return "other"


def load_contrast(folder):
    """Every significant enrichment row for one contrast."""
    rows = []
    csvdir = folder / "CSVs"
    if not csvdir.is_dir():
        return rows
    for f in sorted(csvdir.glob("*.xlsx")):
        lineage = f.stem.replace("_Enrichment", "")
        for suffix in ("_Suppressed_vs_Unsuppressed", "_Suppressed_vs_PreART",
                       "_Unsuppressed_vs_PreART", "_HEI_vs_HEU", "_HEU_vs_HUU"):
            lineage = lineage.replace(suffix, "")
        try:
            xl = pd.ExcelFile(f)
        except Exception as e:                       # a corrupt workbook must
            print(f"  [skip] {f.name}: {e}")         # be named, not swallowed
            continue
        for sheet in SHEETS:
            if sheet not in xl.sheet_names:
                continue
            d = xl.parse(sheet)
            if d.empty or "Adjusted.P.value" not in d.columns:
                continue
            for _, r in d[d["Adjusted.P.value"] < 0.05].iterrows():
                rows.append(dict(lineage=lineage, family=family_of(lineage),
                                 collection=sheet, term=str(r["Term"]),
                                 padj=float(r["Adjusted.P.value"]),
                                 overlap=str(r["Overlap"]),
                                 odds_ratio=float(r.get("Odds.Ratio", float("nan"))),
                                 genes=str(r["Genes"])))
    return rows


def de_direction(contrast_dir, genes):
    """Sign convention: FindMarkers(ident.1 = <first group in the folder name>),
    so a NEGATIVE avg_log2FC means higher in the SECOND group. Verified against
    Scripts/2. Analysis/2. WNN_Annotation.R line 765."""
    up_first = up_second = 0
    per_gene = collections.Counter()
    for f in sorted(contrast_dir.glob("*.csv")):
        if f.stem.startswith("Summary"):
            continue
        for r in csv.DictReader(f.open()):
            g = r.get("") or r.get("gene")
            padj = r.get("p_val_adj")
            if g not in genes or padj in (None, "", "NA"):
                continue
            if float(padj) >= 0.05:
                continue
            fc = float(r["avg_log2FC"])
            if fc > 0:
                up_first += 1
            else:
                up_second += 1
                per_gene[g] += 1
    return up_first, up_second, per_gene


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows, dropped_rows = [], []

    for folder in sorted(p for p in PA.iterdir() if p.is_dir()):
        rows = load_contrast(folder)
        if not rows:
            print(f"  {folder.name}: no enrichment rows")
            continue
        n_lineage = len({r["lineage"] for r in rows})
        spread = collections.Counter()
        for term in {(r["lineage"], r["term"]) for r in rows}:
            spread[term[1]] += 1
        ubiquitous = {t for t, n in spread.items() if n > 0.5 * n_lineage}
        for t in sorted(ubiquitous):
            dropped_rows.append(dict(contrast=folder.name, term=t,
                                     n_lineages_significant=spread[t],
                                     n_lineages_tested=n_lineage))
        kept = [r for r in rows if r["term"] not in ubiquitous]
        for r in kept:
            r["contrast"] = folder.name
        all_rows += kept
        print(f"  {folder.name}: {len(rows)} significant rows across {n_lineage} "
              f"populations; dropped {len(ubiquitous)} ubiquitous terms; "
              f"{len(kept)} kept")

    if not all_rows:
        raise SystemExit("no enrichment rows found - check Analysis/Pathway_Analysis")

    df = pd.DataFrame(all_rows)
    cols = ["contrast", "family", "lineage", "collection", "term", "padj",
            "overlap", "odds_ratio", "genes"]
    df[cols].sort_values(["contrast", "family", "padj"]).to_csv(
        OUT / "R1M2_07_lineage_pathway_summary.csv", index=False)
    pd.DataFrame(dropped_rows).to_csv(
        OUT / "R1M2_09_dropped_ubiquitous_terms.csv", index=False)

    # ---- the pan-lineage result, with its direction ------------------------
    key = "Suppressed_vs_Unsuppressed"
    tnf = df[(df.contrast == key) & df.term.str.startswith("TNF-alpha Signaling")]
    genes = sorted({g.strip() for s in tnf.genes for g in s.split(";") if g.strip()})
    up_sup, up_unsup, per_gene = de_direction(DE / key, set(genes))

    print(f"\n  TNF-alpha/NF-kB, {key}:")
    print(f"    populations significant: {tnf.lineage.nunique()}"
          f"  families: {sorted(tnf.family.unique())}")
    print(f"    driver genes: {len(genes)}")
    print(f"    significant gene-by-population instances: {up_sup + up_unsup}")
    print(f"      higher in SUPPRESSED:   {up_sup}")
    print(f"      higher in UNSUPPRESSED: {up_unsup}")
    print(f"    most recurrent drivers: "
          f"{', '.join(g for g, _ in per_gene.most_common(10))}")

    pd.DataFrame([dict(gene=g, n_populations_higher_in_unsuppressed=n)
                  for g, n in per_gene.most_common()]).to_csv(
        OUT / "R1M2_08_lineage_pathway_drivers.csv", index=False)

    # Per-family table for the figure: strongest term per family per contrast.
    import math
    fig = (df.sort_values("padj")
             .groupby(["contrast", "family", "term"], as_index=False)
             .first())
    fig["neglog10padj"] = fig["padj"].apply(lambda x: -math.log10(x) if x > 0 else 300.0)
    fig.to_csv(OUT / "R1M2_07_lineage_pathway_figure_data.csv", index=False)
    print(f"\n  wrote 4 tables to {OUT.name}/")


if __name__ == "__main__":
    main()
