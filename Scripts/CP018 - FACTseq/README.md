# CP018 FACT-Seq analysis pipeline

Second FACT-Seq participant, extending the CP003 proof-of-concept (Figure 3).
Mirrors the CP003 script structure so results are directly comparable.

## Why CP018 is the right second participant

CP003 **never suppressed** within its sampled scRNA-seq timepoints (VL >500,000
at both 12m and 24m). CP018 did — and its existing TARA data already shows the
manuscript's central pattern within one participant:

| TARA timepoint | ART status | Viral load | Exhaustion | Stemness |
|---|---|---|---|---|
| `CP018_1m` | Pre-ART | 176,970 | +0.006 | −0.071 |
| `CP018_25m` | Suppressed | 20 | −0.072 | +0.013 |
| `CP018_42m` | Suppressed | 20 | −0.188 | +0.047 |

*(from `Analysis/CD8_subset/09_viral_load/PerSample_ModuleScores_and_ViralLoad.csv`)*

So the FACT-Seq stim data (2m and 90m) brackets a suppression trajectory, and
the testable prediction is specific: **CP018's HIV-specific clonotypes should
carry the type I IFN memory signature (IFIT1/IFIT3/ISG15/MX1) and preserved
stemness that Figure 4/5 attribute to suppression — where CP003's, sampled
under persistent viremia, should not.**

Second advantage: CP018 was run in **duplicate at both timepoints** (2m = OB1/OB2,
90m = OB3/OB4). CP003 had replicates only at 2m. This supports a replicate
concordance analysis and a confidence tier for clonotypes seen in both
replicates of a timepoint.

## Run design

| sample_id | OCM | Tube | Timepoint | Replicate |
|---|---|---|---|---|
| `CP018_2m_A` | OB1 | 1 | 2m | A |
| `CP018_2m_B` | OB2 | 2 | 2m | B |
| `CP018_90m_A` | OB3 | 3 | 90m | A |
| `CP018_90m_B` | OB4 | 4 | 90m | B |

## Scripts

Run in order. Every script begins with `source(".../config.R")`.

| # | Script | Does | Writes |
|---|---|---|---|
| — | `config.R` | All paths, sample sheet, QC thresholds, helpers | — |
| 1 | `1.QC_Full.R` | H5 → Seurat, QC, cell cycle, Azimuth, scDblFinder | `CP018_postQC_..._singlets.qs2` |
| 2 | `2.Replicate_Concordance.R` | A vs B: recovery, QC, pseudobulk r, composition | `Analysis/CP018/Replicate_Concordance/` |
| 3 | `3.Integration_Anno_Plots.R` | CCA + FastMNN, clustree, marker plots | `CP018_RNA_integrated_CCA_MNN.qs2` |
| 4 | `4.VDJ_Seurat.R` | scRepertoire, 2m-vs-90m overlap, attach to Seurat | `CP018_..._withTCR.qs2` |
| 5 | `5.HIV_Specific_TCR.R` | Define HIV-specific clonotypes | `CP018_HIVSpecificTCR_annotated.qs2` |
| 6 | `6.HIV_Specific_TCR_TARA_bulk.R` | Cross-reference into unstimulated TARA | `Analysis/CP018/HIV_Specific/TARA_crossref/` |

### Config, not hardcoded paths

The CP003 scripts hardcode `/home/akshay-iyer/...` in ~40 places, which is why
they broke when the project moved. Everything here lives in `config.R` — one
edit repoints the whole pipeline.

## One required manual step

**Script 5 will not produce an HIV-specific call until you set
`BYSTANDER_CLUSTERS`.** For CP003 the bystander compartment was MNN clusters
0 and 1, but cluster numbering is not portable between datasets. Script 5 first
writes `CP018_cluster_profile_for_bystander_call.csv` ranking every cluster by
stemness, expansion rate, cytotoxicity, and cycling fraction. The bystander
cluster is high-stemness, near-zero-expansion, low-cytotoxicity. Set the vector
at the top of script 5 and re-run.

This is deliberate: silently reusing `c(0, 1)` from CP003 would produce a
plausible-looking but wrong answer.

## Differences from the CP003 scripts

1. **Subfolder names** — `CP018_2m_A` etc., so the CP003 regex
   `CP003_\d+m_(.+)$` is replaced by the `CP018_SAMPLES` sample sheet.
2. **No CD8+/CD8− split** — CP003 mapped OCM `004` to a CD8− fraction.
   All four CP018 libraries are the same sort, so `CellType_Sort` is constant
   and dropped as a grouping variable.
3. **`Replicate` / `Timepoint_Rep` added** as grouping variables throughout.
4. **Replicate concordance** (script 2) is entirely new.
5. **Two repertoire views** in script 4 — per-library for Seurat attachment,
   replicate-merged for `clonalOverlap`/diversity (which need unique sample
   names per biological group).
6. **Replicate-confirmed clonotype tier** in script 5.

## Inputs verified present

All exist in `saved_R_data/`:
`TARA_ALL_annotated_final.qs2`, `TARA_cd8_HEI_annotated_final.qs2`,
`seu_CP003_HIVSpecificTCR_annotated.qs2`, `TARA_TCR_Combined.RData`.

Requires Cell Ranger output at
`~/Desktop/Projects/cellranger_out/CP018_multi_out/outs/per_sample_outs/`
(see `Scripts/0. Data_Acquisition/`).

## Aligner version caveat

CP018 is being aligned with **Cell Ranger 10.1.0**; CP003 and the manuscript
used **9.0.1**. Cell calling and V(D)J assembly differ between versions, so any
CP018-vs-CP003 difference carries a pipeline-version component. State it in
Methods, or re-align CP003 with 10.1.0 if the comparison becomes central.
