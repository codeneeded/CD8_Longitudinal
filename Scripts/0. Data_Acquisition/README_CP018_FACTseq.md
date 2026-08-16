# CP018 FACT-Seq — data acquisition & alignment

Second FACT-Seq participant, extending the CP003 proof-of-concept (Figure 3).
Data from L. de Armas, sequenced 2026-08-07.

## Run design

One 10x chip, on-chip multiplexing (OCM). Both timepoints in technical duplicate
— note this differs from CP003, which had duplicates only at 2m.

| sample_id | OCM barcode | Tube label | Timepoint |
|---|---|---|---|
| `CP018_2m_A` | OB1 | 1 | 2 months, rep A |
| `CP018_2m_B` | OB2 | 2 | 2 months, rep B |
| `CP018_90m_A` | OB3 | 3 | 90 months, rep A |
| `CP018_90m_B` | OB4 | 4 | 90 months, rep B |

**Why CP018 matters:** CP003 never suppressed within the sampled scRNA-seq
timepoints (VL >500,000 at both 12m and 24m). CP018 suppressed at 65 days
(pre-ART VL 176,970) and has longitudinal timepoints at 1m/25m/42m. So this run
gives the HIV-specific compartment in a *suppressed* participant — the
comparison the Figure 3 analysis could not make.

## Libraries as delivered

OCM pools all four samples into one library per modality; `cellranger multi`
deconvolves by barcode. Both modalities are present:

| Library | fastq_id | Size |
|---|---|---|
| Gene expression | `DeArmas-34125-001_OCM-GEX5` | ~30 GB (R1+R2) |
| TCR V(D)J | `DeArmas-34125-001_OCM-TCR` | ~7 GB (R1+R2) |

Chemistry verified from FASTQ headers: **10x 5' v2** — R1 28 bp (16 bp barcode +
12 bp UMI), R2 90 bp, NovaSeq X flowcell `2532FJLT4` lane 8. Matches CP003.

## Layout

Neither raw data nor references are committed (`fastq_raw/` is gitignored;
references and Cell Ranger output live outside the repo entirely).

```
~/Desktop/Projects/
├── Repositories/CD8_Longitudinal/      # git repo
│   ├── fastq_raw/CP018_FACTseq/        # BaseSpace download (gitignored)
│   ├── saved_R_data/                   # Seurat objects (gitignored)
│   └── Scripts/0. Data_Acquisition/    # these scripts (committed)
├── Tools/cellranger-10.1.0/            # aligner
├── refdata/
│   ├── refdata-gex-GRCh38-2024-A/                         # 16 GB
│   └── refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/  # 12 MB
├── Manuscripts/CD8_Longitudinal/       # manuscript + supplementary
└── cellranger_out/CP018_multi_out/     # multi output
```

## References — staged

Both match the manuscript STAR Methods exactly:

- **GEX**: `refdata-gex-GRCh38-2024-A` — GRCh38, GENCODE v44, mkref 8.0.0
- **VDJ**: `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` — Ensembl 94

Tarballs are kept alongside the extracted dirs; delete
`refdata-gex-GRCh38-2024-A.tar.gz` (11 GB) once alignment succeeds.

## Steps

```bash
cd ~/Desktop/Projects/Repositories/CD8_Longitudinal

# 1. Download FASTQs (resumable; skips complete files)
PROJECT_ID=512350847 bash "Scripts/0. Data_Acquisition/download_CP018_FACTseq.sh"

# 2. Preflight — checks binary, references, FASTQ integrity, config paths, disk
bash "Scripts/0. Data_Acquisition/run_CP018_cellranger_multi.sh" --check

# 3. Align (~6-14 h on 8 cores)
bash "Scripts/0. Data_Acquisition/run_CP018_cellranger_multi.sh"
```

**Cell Ranger version note.** `Tools/cellranger-10.1.0/` is installed, but the
CP003 FACT-Seq data and everything in the manuscript were processed with
**v9.0.1** (STAR Methods). OCM config syntax (`ocm_barcode_ids`, `OB1`-`OB4`) is
unchanged between the two, so 10.1.0 runs this config as-is — but see
"Version decision" below before aligning.

## Output layout

`cellranger multi` writes the same structure the CP003 R scripts already read:

```
cellranger_out/CP018_multi_out/outs/per_sample_outs/<sample_id>/
├── count/sample_filtered_feature_bc_matrix.h5
└── vdj_t/filtered_contig_annotations.csv
```

## Downstream

Adapt `Scripts/CP003 - Longitudinal HIV Stim/1. QC_Full.R`, setting `base_dir`
to the `per_sample_outs` path above. Three differences from CP003 to handle:

1. **Subfolder names** are `CP018_2m_A` / `CP018_90m_B`, so the CP003 regex
   `CP003_\\d+m_(.+)$` needs the PID and barcode pattern updated (`A`/`B`
   rather than `001A`/`003`).
2. **No CD8+/CD8− sort split.** CP003 mapped OCM `004` to `CD8-`; here all four
   samples are the same sort, so `CellType_Sort` is constant and the
   `case_when` on barcode should be dropped.
3. **True duplicates at both timepoints** — quantify replicate concordance
   (cell recovery, per-gene correlation, clonotype overlap A vs B) before
   merging, which CP003 could not do at 101m.

## Version decision

CP018 will be compared directly against CP003 (Figure 3). Aligner version is a
batch variable: 10.1.0 vs 9.0.1 differ in cell-calling and V(D)J assembly, so
any CP018-vs-CP003 difference would be confounded with pipeline version.

Two defensible options:

1. **Install 9.0.1 alongside** and align CP018 with it — keeps CP018 directly
   comparable to CP003 and to the manuscript's stated methods. Preferred if
   CP018 is going into the existing manuscript.
2. **Use 10.1.0 for CP018 and re-align CP003 with 10.1.0** — both on the newer
   pipeline, consistent, but means re-running CP003 and regenerating Figure 3.

Do *not* align CP018 with 10.1.0 and compare against 9.0.1-processed CP003
without stating it as a limitation.

## Provenance

- BaseSpace project `DeArmas-30-06-34125` (id `512350847`), 72.79 GB, NovaSeq
- Share link: <https://basespace.illumina.com/s/fXOEy9uevGcc>
- Raw data goes to GEO on publication, alongside GSE254645 / GSE306003 /
  GSE328561 / GSE328716.
