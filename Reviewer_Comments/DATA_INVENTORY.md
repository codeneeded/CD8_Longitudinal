# TARA / CD8_Longitudinal — data and code inventory

Internal reference. What each saved object is, which script writes it, which
script reads it, and which manuscript figure depends on it. Built by reading
every object in `saved_R_data/` (47 files, 34.3 GB) and cross-referencing
every `qs_save`/`saveRDS`/read call in `Scripts/`.

Regenerate with `Scripts/2. Analysis/inventory_objects.R`.

## The two clustering lineages — read this first

The manuscript uses **two different CD8 clusterings**, and this is the source of
Reviewer 2's Major 1d and Reviewer 1's Minor 4:

- **Figures 1 and 2** use `TARA_ALL_annotated_final.qs2` — CD8 clusters as they
  fall out of **whole-PBMC WNN clustering**. This is what the reviewers looked at.
- **Figures 4, 5 and Supplementary 4-9** use `TARA_cd8_HEI_annotated_final.qs2` —
  a **dedicated CD8-only re-analysis** (re-normalised, re-scaled, RNA PCA + MNN
  integration, ADT PCA + MNN integration, `FindMultiModalNeighbors`, WNN
  clustering swept over resolutions 0.2-2.0 with clustree).

So the dedicated CD8 sub-clustering the reviewer asks for **already exists** — it
is just not shown, and Figure 1 does not use it. That reframes Major 1d from new
analysis to presenting analysis already done, plus deciding whether Figure 1
should be rebuilt on the CD8-only object for consistency.

## Main pipeline chain

| Step | Script | Writes | Cells |
|---|---|---|---|
| Protein normalisation | `1. QC/1.DSB_Norm.R` | — | — |
| QC | `1. QC/2.QC.R` | — | — |
| Doublets + cell cycle | `1. QC/3.QC2_Doublet_Cell_Cycle.R` | — | — |
| Isotype correction | `1. QC/4.Isotype_Correction.R` | — | — |
| Integration | `1. QC/5. Integration.R` | — | — |
| TCR/BCR attach | `2. Analysis/1. V(D)JscRepitoire.R` | — | — |
| WNN + annotation | `2. Analysis/2. WNN_Annotation.R` | `TARA_ALL_checkpoint.qs2` -> `TARA_ALL_annotated_final.qs2`, `TARA_cd8_final.qs2` | 126435 |
| Module scoring | `2. Analysis/3.Module_Scoring.R` | `TARA_ALL_ModuleScored.qs2` | 126435 |
| **Dedicated CD8 re-analysis** | `2. Analysis/4.CD8_Subseting.R` | `TARA_cd8_HEI_integrated.qs2` -> `_annotated.qs2` -> `_with_trajectory.qs2` -> `_annotated_final.qs2` | 18843 |
| CMV vs expanding | `2. Analysis/5.CMV_vs_Expanding_Clones.R` | `TARA_TRB_ED1_CMVgroup.qs2`, `TARA_CMV_vs_OtherExpanding_ED1_markers.qs2` | 142481 |
| BCR | `2. Analysis/BCR/2.Clonal_Analysis.R` | `bcr_results_full.qs2`, `bcr_results_heavy.qs2` | — |
| CP003 FACT-Seq | `CP003 - Longitudinal HIV Stim/7.*.R` | `cp003_fig3_annotated.qs2`, `tara_cp003_fig3.qs2` | 4032 / 14719 |
| CP018 FACT-Seq | `CP018 - FACTseq/1-11` | `CP018_*.qs2`, `tara_cp018_crossref.qs2` | 13880 |

## Objects the manuscript figures actually depend on

| Object | Cells | Assays | Used by |
|---|---|---|---|
| `TARA_ALL_annotated_final.qs2` | 126435 | RNA+ADT+mnn.reconstructed | Figure 1, Figure 2, and as input to 3.Module_Scoring, 4.CD8_Subseting, both FACT-Seq cross-references |
| `TARA_ALL_ModuleScored.qs2` | 126435 | RNA+ADT+mnn.reconstructed | Figure 1 |
| `TARA_cd8_HEI_annotated_final.qs2` | 18843 | RNA+ADT+mnn.reconstructed | Figure 4, Figure 5, Supplementary 4-9, `v12_Manuscript_Answers.R`, and CP018 `config.R` |
| `cp003_fig3_annotated.qs2` | 4032 | RNA+prediction.score.celltype.l1+prediction.score.celltype.l2+prediction.score.celltype.l3+mnn.reconstructed | Figure 3 |
| `tara_cp003_fig3.qs2` | 14719 | RNA+ADT+mnn.reconstructed | Figure 3 |

## Lineage subsets that already exist (Reviewer 1, Major 2)

These are on disk but **no script reads or writes them** — provenance unknown, so
confirm how they were made before using them in the revision.

| Object | Cells | Referenced in any script |
|---|---|---|
| `TARA_CD4_subset.qs2` | 51660 | no |
| `TARA_B_cell_subset.qs2` | 17671 | no |
| `TARA_NK_subset.qs2` | 16877 | no |
| `TARA_CD8_subset.qs2` | 29978 | no |

The safer route for Reviewer 1 Major 2 is to subset from
`TARA_ALL_annotated_final.qs2` with a tracked script, since its Azimuth calls are
already there (CD4 T 56,963 / CD8 T 26,908 / B 17,018 / NK 16,146 / Mono 5,407 /
other T 3,440 / DC 366 / other 187).

## Cell-count discrepancy — needs resolving before the revision

Object cell counts do not agree, and the revision requires a defensible per-sample
N for both reviewers:

| Object | Cells |
|---|---|
| `TARA_TRB_ED1_CMVgroup.qs2` | 142481 |
| `TARA_ALL_post_annotation_comparisonsplit.rds` | 142434 |
| `TARA_ALL_pass1.qs2` / `_sorted*.qs2` / `_post_annotation.rds` | 135765 |
| `TARA_ALL_annotated_final.qs2` | 126435 |

The final object is the smallest, consistent with clusters being dropped during
annotation — but which clusters, and why, has to be stated in Methods now that
both reviewers want per-participant counts.

## Untraced objects

45 of 47 are Seurat objects; 2 are marker `data.frame`s. Only 11 have a
`qs_save` call in a script; the rest are earlier iterations
(`_v2`/`_v3`/`_v4`/`_refined`/`_pass1`) or were made interactively. **None of the
untraced objects feed a manuscript figure** — every figure input is in the table
above — so they are safe to leave alone, but they should not be used in the
revision without re-deriving them from a script.

Full machine-readable detail, one row per object, in `object_inventory.csv`.

## Draft Code Availability

> All analysis code is available at https://github.com/codeneeded/CD8_Longitudinal.
> The repository is organised as sequential pipeline stages: quality control and
> normalisation (`Scripts/1. QC`), V(D)J processing, whole-PBMC WNN clustering and
> annotation, module scoring, dedicated CD8+ T cell sub-clustering, and CMV
> comparison (`Scripts/2. Analysis`), figure generation (`Scripts/3. Manuscript`),
> and the two HIV peptide-stimulation (FACT-Seq) analyses
> (`Scripts/CP003 - Longitudinal HIV Stim`, `Scripts/CP018 - FACTseq`). Each
> figure panel is accompanied by a source-data CSV containing the values plotted.
> Package versions are recorded in the repository.

## Draft Data Availability — NEEDS INPUT

Four GEO accessions are referenced in
`Scripts/0. Data_Acquisition/README_CP018_FACTseq.md`: **GSE254645, GSE306003,
GSE328561, GSE328716**. Which samples each covers is not recorded anywhere in the
repo, and CP018 FACT-Seq is not yet deposited. Confirm the mapping before writing
this statement, because Reviewer 1's HIV-alignment request depends on the raw data
being retrievable per sample.

> Raw sequencing data (scRNA-seq, scTCR-seq and ADT) generated in this study are
> available in the Gene Expression Omnibus under accessions GSE254645, GSE306003,
> GSE328561 and GSE328716 [CONFIRM SAMPLE MAPPING]. Data for the HIV
> peptide-stimulation (FACT-Seq) experiments are deposited under [PENDING].
> Processed Seurat objects are available from the corresponding author on request.
> Clinical metadata are provided in Supplementary Table 1; individual-level
> clinical data are restricted to protect participant privacy in this small
> cohort and are available on request under a data use agreement.
> Source Data are provided with this paper.
