# CD8 Longitudinal — Single-Cell Multi-Omic Analysis of CD8+ T Cells in Perinatal HIV

> Longitudinal CITE-seq (paired scRNA-seq + TCR-seq + surface protein) analysis of CD8+ T cell stemness, exhaustion, clonal dynamics, and ART-mediated transcriptional remodeling in perinatally HIV-infected infants from the TARA cohort. Code repository for a manuscript in preparation.

---

## Manuscript

**Single-cell multi-omic profiling reveals ART-mediated remodeling of CD8+ T cell stemness, exhaustion, and interferon memory in perinatally HIV-infected infants**

Akshay Iyer\*, Lesley de Armas\*, Benjamin Bone, Paula Vaz, Maria Grazia Lain, Rajendra N. Pahwa, and Savita G. Pahwa

\* These authors contributed equally

*Manuscript in preparation. Raw sequencing data will be deposited at GEO upon publication. All analysis code is available in this repository.*

---

## Overview

This repository contains the complete analysis code for a longitudinal single-cell multi-omic study of CD8+ T cells in the **TARA (Towards AIDS Remission Approaches)** cohort — perinatally HIV-infected (HEI) and HIV-exposed uninfected (HEU) infants from Maputo, Mozambique. The study applies **CITE-seq** (simultaneous scRNA-seq + TCR-seq + antibody-derived tag [ADT] surface protein profiling) at three longitudinal timepoints (pre-ART, ~12 months on ART, ~24 months on ART) to resolve the transcriptional, phenotypic, and clonal dynamics of the CD8+ T cell compartment during early HIV infection and antiretroviral therapy.

### Key Findings

- **10 transcriptionally distinct αβ CD8+ T cell populations** resolved, including 4 naïve subsets with divergent metabolic and activation states
- **TEMRA/CTL CD8+ T cells** are significantly expanded in HEI infants (p = 0.022 vs HEU; p = 0.026 vs HUU) and near-absent in uninfected infants
- **Tex frequency correlates with viral burden** (Spearman ρ = 0.74, p = 0.013), providing quantitative evidence of antigen-driven exhaustion
- **24% of HEI CD8+ T cells are clonally expanded** (4,910 of 20,418 cells; 821 unique clonotypes) vs 0.6% in HEU and 0% in HUU
- **ART suppression remodels the transcriptional program of expanding clones**: actively reduces exhaustion below pre-treatment levels, preserves stemness (*TCF7*, *SELL*, *BACH2*), and installs a durable **type I IFN memory signature** (*IFIT1*, *IFIT3*, *ISG15*, *MX1*)
- **ART-dependent naïve CD8 fate divergence**: viremia drives a metabolic stress trajectory (*Naïve CD8 2*: *IGFBP2*, *HILPDA*, *PFKFB4*), while ART enables an IFN-primed quiescent state (*Naïve CD8 3*: *BACH2*, *BCL2*, *IFIT1*)
- **A CD16+ Effector CD8 subset** (ADCC-competent; *FCGR3A*+, 41% clonally expanded) displays the strongest type I IFN memory signature under suppression
- **HIV-specific CD8+ T cell clonotypes** established within the first month of life persist as functionally armed effectors across 8+ years of childhood (validated by HIV peptide stimulation in participant CP003)
- **Intact proviral reservoir burden** correlates directionally with HIV-specific clonal persistence, suggesting antigen availability sustains the effector repertoire

---

## Cohort

| Feature | Details |
|---|---|
| **Cohort** | TARA (Towards AIDS Remission Approaches) |
| **Location** | Maputo Province, Mozambique |
| **HEI participants** | n = 11 (cross-sectional); n = 5 with complete longitudinal scRNA-seq data |
| **HEU participants** | n = 3 |
| **HUU participants** | n = 2 |
| **Enrollment** | Age 1–2 months (pre-ART initiation, median ART start 33 days) |
| **Longitudinal timepoints** | Pre-ART (~1 month), ~12 months on ART, ~24 months on ART |
| **Extended follow-up** | CP018 (42 months), CP020 (44 months); CP003 functional validation at 101 months |
| **Funding** | NIH/NIAID R01 AI127347 (S. Pahwa) |

---

## Repository Structure

```
CD8_Longitudinal/
│
├── Scripts/                          # R analysis scripts
├── Analysis/                         # Processed outputs, figures, tables
├── Integration/                      # Batch correction and multimodal integration
│                                     # (FastMNN, Seurat CCA, WNN)
├── Manuscript/                       # Publication-ready figures and tables
├── QC/                               # Quality control outputs (doublet removal,
│                                     # filtering thresholds, library metrics)
│
├── TARA_Longitudinal_Trends_v1.xlsx  # CD8+ T cell phenotype longitudinal data
├── TARA_VL_CD4.xlsx                  # Viral load and CD4 count data per timepoint
│
├── CD8_Longitudinal.Rproj            # RStudio project file
└── README.md
```

---

## Analysis Pipeline

### 1. Sequencing & Alignment
- **Platform**: 10x Genomics Chromium 5' Immune Profiling V2
- **Libraries**: Gene expression (mRNA), V(D)J (TCR), and surface protein (ADT — BioLegend TotalSeq-C Human Universal Cocktail v3)
- **Sequencing**: Illumina NovaSeq X Plus
- **Alignment**: Cell Ranger v9.0.1; GRCh38 reference (refdata-gex-GRCh38-2024-A); VDJ GRCh38 reference (refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0)
- **Multiplexing**: 10x Genomics on-chip multiplexing (OCM) for the CP003 peptide stimulation experiment

### 2. Quality Control (`QC/`)
- Cells with <200 detected genes or >20% mitochondrial RNA removed
- Doublets identified and removed with **scDblFinder**

### 3. Normalization & Integration (`Integration/`)
- RNA: log-normalization and scaling
- ADT: **DSB** (denoised and scaled by background) normalization with isotype controls
- Batch correction: **FastMNN** and **Seurat CCA**
- Multimodal integration: **Seurat WNN** (weighted nearest neighbor) combining RNA and ADT modalities
- Optimal clustering resolution (0.4) selected using **clustree**

### 4. Cell Type Annotation
- Broad PBMC annotation: **Azimuth** (automated) + manual curation
- CD8+ T cell populations annotated using curated panels of surface proteins and transcripts
- Sub-clustering of αβ CD8+ T cells at resolution 0.4 resolved 10 populations
- ADT sub-gating resolved Tscm and Naïve Intermediate populations by CD45RA/CD45RO/FAS co-expression

### 5. Core Analysis (`Scripts/`, `Analysis/`)

**Differential abundance and expression:**
- CD8 cluster frequencies compared across HEI, HEU, HUU and across ART status (pre-ART, suppressed, unsuppressed)
- Differential gene expression with **MAST**, Wilcoxon tests; Bonferroni/BH correction
- Seven curated functional gene modules scored per cell: Exhaustion, Stemness/Naïve, Cytotoxicity, Stress/IFN (acute), Type I IFN memory, Inflammatory chemokines, Terminal differentiation (Table 1 in manuscript)

**Clonal dynamics (TCR repertoire):**
- Clonotype analysis with **scRepertoire**
- Longitudinal clonal tracking with alluvial plots (top 20 clonotypes per participant)
- Cross-cluster clonotype sharing quantified per participant
- Clonal overlap across timepoints: Jaccard index
- Epitope annotation of expanded TRB CDR3 sequences using **Trex** against VDJdb, McPAS-TCR, and IEDB (Levenshtein edit distance ≤ 2 aa)

**Trajectory analysis:**
- Pseudotime inference with **Monocle3** across all αβ CD8+ T cells
- Three branches resolved: effector trajectory + two ART-status-divergent naïve branches

**Viral load correlation:**
- Spearman correlations between CD8 subset frequency/module scores and plasma HIV-1 RNA
- Per-participant mean exhaustion score vs log₁₀ VL: ρ = 0.73, p_adj = 0.002

### 6. HIV Peptide Stimulation & Functional Validation (`Scripts/`)
- PBMCs from CP003 at 2 months and 101 months stimulated with HIV Consensus C peptide pools
- Proliferating CD8+ T cells sorted by CellTrace dye dilution and profiled by scRNA-seq + TCR-seq (10x OCM)
- 4,032 CD8+ T cells retained post-QC; six functional populations resolved
- HIV-specific clonotypes cross-referenced against unstimulated TARA longitudinal data
- Clonal persistence tracked from 1 month to 101 months of age

### 7. HIV-1 Reservoir Quantification (`Scripts/`)
- Total HIV-1 DNA: LTR-targeted real-time PCR (TWC assay) across all 5 HEI scRNA-seq participants
- Intact proviral DNA: Cross-Subtype IPDA (CS-IPDA) on nanoplate digital PCR (QIAcuity) for CP006 and CS-IPDA/IPDA from prior runs for CP003, CP013, CP020 (CP018 failed QC)
- Reservoir metrics correlated with clonal persistence (Spearman, descriptive)

---

## CD8+ T Cell Populations Resolved

| Population | Key Markers | ART Association |
|---|---|---|
| Naive CD8 1 | CCR7, CD62L, TCF7, LEF1, CD73hi, CD127hi | — |
| Naive CD8 2 | CD38hi, FAS, CD45RO, *IGFBP2*, *HILPDA*, *PFKFB4* | Pre-ART (82%; 1.67×) |
| Naive CD8 3 | *BACH2*, *BCL2*, *IFIT1*, *ISG15*, CD69, ITGB7 | Suppressed (52%; 1.81×) |
| Naive CD8 4 | *GZMK*+, proliferative | — |
| Tscm CD8 | CD45RA+, CD45RO−, FAS+ | Pre-ART enriched |
| Naive Intermediate CD8 | CD45RA+CD45RO+ (double-positive), CCR7, TCF7 | Pre-ART enriched |
| TEM CD8 | Effector memory markers | Clonally expanding |
| TEMRA CD8 | CD57, KLRG1, KIR3DL1, GZMB, GNLY, PRF1 | Clonally expanding |
| CD16+ Effector CD8 | FCGR3A+, ADCC-competent | Strongest IFN memory under ART |
| KIR+ Innate-like CD8 | KIR3DL1 protein, CD94, CD57 | Suppressed enriched (49.6%) |

---

## Functional Gene Modules (Table 1)

| Module | Key Genes |
|---|---|
| Exhaustion | *TOX, PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, ENTPD1* |
| Stemness / Naïve | *TCF7, SELL, CCR7, LEF1, BACH2, BCL2, IL7R, S1PR1, KLF2* |
| Cytotoxicity | *GZMB, GNLY, PRF1, NKG7, GZMA, GZMH, FGFBP2* |
| Stress / IFN (acute) | *HSPA1A, HSPA1B, DNAJB1, IFI27, IFI44L* |
| Type I IFN memory | *IFIT1, IFIT3, ISG15, MX1* |
| Inflammatory chemokines | *CCL3, CCL3L3, CCL4, CCL4L2, CCL5* |
| Terminal differentiation | *ZEB2, PRDM1, TBX21, CX3CR1, S1PR5, ID2* |

---

## Key Data Files

| File | Description |
|---|---|
| `TARA_Longitudinal_Trends_v1.xlsx` | Longitudinal CD8+ T cell phenotype data (flow cytometry) for TARA cohort participants |
| `TARA_VL_CD4.xlsx` | Plasma HIV-1 viral load (copies/mL) and CD4+ T cell count per participant per timepoint |

> Raw scRNA-seq, TCR-seq, and ADT data will be deposited at NCBI GEO upon publication.

---

## Dependencies

All scripts are written in **R (v4.4)**. Key packages:

| Package | Purpose |
|---|---|
| `Seurat` (v5) | WNN multimodal clustering, DGE, visualization |
| `scDblFinder` | Doublet detection and removal |
| `dsb` | ADT normalization |
| `FastMNN` / `batchelor` | Batch correction |
| `scRepertoire` | TCR clonotype analysis and alluvial plots |
| `Monocle3` | Pseudotime trajectory inference |
| `MAST` | Differential gene expression |
| `clustree` | Clustering resolution selection |
| `readxl` / `dplyr` / `ggplot2` | Data ingestion and visualization |
| `pheatmap` / `corrplot` | Heatmaps and correlation matrices |

Install Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scDblFinder", "batchelor", "MAST"))

# scRepertoire and Monocle3 from GitHub/CRAN:
install.packages("scRepertoire")
# Monocle3: follow installation at https://cole-trapnell-lab.github.io/monocle3/
```

---

## Related Publication

This work extends the following published study from the same cohort:

- **de Armas LR, Dinh V, Iyer A, Pallikkuth S, Pahwa R, Cotugno N, Rinaldi S, Palma P, Vaz P, Lain MG, Pahwa SG.** (2024). Accelerated CD8+ T cell maturation in infants with perinatal HIV infection. *iScience*, 27(5), 109720. [DOI: 10.1016/j.isci.2024.109720](https://doi.org/10.1016/j.isci.2024.109720)
  > Cross-sectional flow cytometry and scRNA-seq showing accelerated CD8+ T cell maturation and reduced self-renewal markers pre-ART in the same TARA cohort. The present manuscript extends those findings longitudinally at single-cell multi-omic resolution under ART.

---

## Citation

A full citation will be added upon publication. In the interim, please reference this repository and acknowledge the TARA cohort (NIH/NIAID R01 AI127347).

---

## Contact

For questions, please open an [issue](https://github.com/codeneeded/CD8_Longitudinal/issues) or contact the corresponding author: spahwa@med.miami.edu
