# CP018 FACT-Seq pipeline

HIV-peptide-stimulated FACT-Seq (5' GEX + TCR) on CP018 at 2 and 90 months
post-ART, each timepoint in technical duplicate, plus cross-reference into the
older unstimulated eCITE-seq data and joint integration with CP003.

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

## Run everything

```bash
cd ~/Desktop/Projects/Repositories/CD8_Longitudinal
bash "Scripts/CP018 - FACTseq/run_pipeline.sh"
```

Runs every stage in order, stops at the first failure, and writes each stage's
full output to `Analysis/CP018/logs/<stage>.log`.

```bash
bash "Scripts/CP018 - FACTseq/run_pipeline.sh" --list    # show stages
bash "Scripts/CP018 - FACTseq/run_pipeline.sh" 4b        # resume from a stage
bash "Scripts/CP018 - FACTseq/run_pipeline.sh" --plots   # figures only
```

`--plots` regenerates every figure and table from the saved objects and skips
the two stages that rebuild objects rather than draw them (stage 3
re-integration, stage 9 when its checkpoint exists). Use it after a plotting
change.

If `Rscript` is not on your `PATH`, set it: `RSCRIPT=/path/to/Rscript bash ...`

## Stages

| Stage | Script | Does |
|---|---|---|
| 3a | `3a.Set_Resolution_Annotate.R` | Set resolution 0.4, apply annotation, annotated plots. **No re-integration** (~2 GB). |
| 3 | `3.Integration_Anno_Plots.R` | Full CCA + FastMNN re-integration. ~25 GB peak; auto-skipped when the object exists. |
| 3b | `3b.Annotation_Evidence.R` | Annotation justification table, marker panels, CP003 correspondence mapping. |
| 4 | `4.VDJ_Seurat.R` | Attach TCR, repertoire plots, clonal overlap. |
| **0** | `0.Normalise_Clonotypes.R` | **Stamp the canonical clonotype column.** Runs after 4 (needs `CTstrict`) and before anything cross-experiment. |
| 4b | `4b.Clonotype_Persistence.R` | Persistence alluvials, detection-power controls, persistence by clone size. |
| 5 | `5.HIV_Specific_TCR.R` | HIV-specific call (clonal expansion >2 within a timepoint). |
| 6 | `6.HIV_Specific_TCR_TARA_bulk.R` | Cross-reference into unstimulated eCITE-seq. |
| 6b | `6b.Lineage_Timeline.R` | Lineage timeline pre-ART → 90m, resting phenotype, Tscm test. |
| 8 | `8.CP003_vs_CP018_Comparison.R` | CP003 vs CP018 comparison, matched-rule recall. |
| 9 | `9.Integrate_CP003_CP018.R` | Joint integration, Tex-survival check. |

`config.R` holds all paths, the sample sheet, palettes, plotting helpers and the
annotation map. Every script sources it; nothing else hardcodes a path.

## The clonotype-format problem, and how it is handled

Cell Ranger 10.1.0 writes an **empty** D-gene field where 9.0.1 and the older
eCITE-seq alignment write a literal `NA`:

```
CP003 stim (CR 9.0.1)      TRBV6-1.NA.TRBJ2-5.TRBC2;TGTGCC...
unstim eCITE-seq (older)   TRBV15.NA.TRBJ2-3.TRBC2;TGTGCC...
CP018 stim (CR 10.1.0)     TRBV29-1..TRBJ1-2.TRBC1;TGCAGC...   <- ".."
```

**Consequence if ignored:** intersecting raw `CTstrict` across these experiments
returns **zero** matches by construction — indistinguishable from real clonal
absence. This is not hypothetical; it silently produced a false negative before
it was caught.

**Fix.** `0.Normalise_Clonotypes.R` stamps a stored column `CTstrict_norm`
(referenced everywhere as `CT_CANON`) onto every object, rewriting `..` → `.NA.`.
Scripts read the stored column rather than recomputing, so a stale object cannot
be compared against a normalised one. `require_ct_canon()` makes any
cross-experiment script fail loudly instead of returning zero.

**CP003 is provably unaffected — this does not change published numbers.**
Normalisation only ever rewrites `..`, so objects already using `.NA.` pass
through untouched. Measured on the actual objects:

| Object | Clonotype strings | Normalised | % |
|---|---:|---:|---:|
| CP018 stim (CR 10.1.0) | 4,645 | 3,911 | 84.2% |
| CP003 stim (CR 9.0.1) | 3,436 | **0** | **0.0%** |
| unstim eCITE-seq | 9,394 | **0** | **0.0%** |

Distinct clonotype counts are unchanged in every object (CP003 2,012 → 2,012),
so no clonotypes merge. Stage 0 **refuses to write** a legacy-flagged object if
normalisation would alter it, and aborts if the transform ever increases the
distinct-clonotype count. It is idempotent — re-running changes nothing.

Effect on recovery in the unstimulated data:

| Matching level | HIV-specific recovered |
|---|---:|
| `CTstrict` raw | 0 / 72 |
| `CTstrict_norm` *(used)* | 3 / 72 |
| `CTaa` (looser) | 5 / 72 |

Normalised `CTstrict` is preferred over `CTaa`: it requires nucleotide sequence
**and** gene calls to agree, whereas `CTaa` collapses distinct nucleotide
sequences encoding the same CDR3 peptide.

**For Methods:** state that cross-experiment clonotype matching used
`CTstrict` with the D-gene field normalised across Cell Ranger versions.

## Two decisions recorded in the code

*Historical note: earlier versions required a manual `BYSTANDER_CLUSTERS` step
before script 5 would produce a call. That is gone — the expansion-based
definition below removed the dependency on cluster numbers entirely, which is
what had broken when the saved object's resolution and the exclusion list
diverged.*

**Clustering resolution 0.4** — chosen from the clustree plus marker profiles,
with the evidence in `3b`'s output rather than in a comment.

**HIV-specific definition** — a clonotype expanded to >2 cells *within a single
timepoint*, with MAIT excluded by annotation label. Deliberately not by cluster
number: the earlier cluster-number version broke when the saved object's
resolution and the exclusion list diverged.

## Outputs

| Directory | Contents |
|---|---|
| `Analysis/CP018/Integration/` | Clustree, markers, annotation evidence, correspondence mapping |
| `Analysis/CP018/Annotation_Plots/` | Annotated UMAPs, composition |
| `Analysis/CP018/VDJ/TCR/` | Repertoire plots, clonal overlap, persistence alluvials |
| `Analysis/CP018/HIV_Specific/` | Per-cell calls with every decision input, cross-reference, lineage timeline |
| `Analysis/CP018/CP003_vs_CP018/` | Comparison tables, joint integration |
| `Analysis/CP018/logs/` | Per-stage run logs |

Every figure ships the table behind it. Start with
`Analysis/CP018/CP018_FINDINGS.md`.

## Caveats that travel with the numbers

1. **n = 1 participant per suppression group.** Any CP003/CP018 difference is
   inseparable from individual identity. Integration gives a shared coordinate
   system, not statistical power.
2. **Cell Ranger versions differ** (9.0.1 vs 10.1.0), affecting cell calling and
   V(D)J assembly. Needs a Methods sentence.
3. **Sampling depth bounds every persistence claim.** Technical replicates of one
   blood draw share only a small fraction of clonotypes, and most clonotypes are
   single cells — so non-detection is not clonal loss.
4. **Do not read cell counts across the two experiments as fold-change.**
   Different assays, cell inputs and sampling depths.

## Aligner version caveat

CP018 is being aligned with **Cell Ranger 10.1.0**; CP003 and the manuscript
used **9.0.1**. Cell calling and V(D)J assembly differ between versions, so any
CP018-vs-CP003 difference carries a pipeline-version component. State it in
Methods, or re-align CP003 with 10.1.0 if the comparison becomes central.
