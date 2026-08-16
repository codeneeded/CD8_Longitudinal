# CP018 FACT-Seq — summary of findings

HIV-peptide-stimulated FACT-Seq (5' GEX + TCR) on CP018 at two timepoints,
**2 months and 90 months** post-ART-initiation — an 88-month interval.
Both timepoints were run in **technical duplicate** (4 libraries, OCM-barcoded,
one Cell Ranger `multi` run).

CP018 is the **ART-suppressed** participant: pre-ART VL 176,970, suppressed to
VL 20 by day 65. The comparator, **CP003, never suppressed** (VL >500,000 at
both sampled timepoints). This is the contrast the analysis is built around.

**5,129 cells** passed QC; **4,645 (90.6%)** carry a productive TCR.

---

## 1. What we annotated

Nine clusters (MNN integration, resolution 0.4) collapse to **seven states**.
Labels were assigned from marker panels spanning stemness, activation,
exhaustion, metabolism, effector function and cycling — then cross-checked
against CP003 (section 4) so shared names mean the same thing in both.

| State | Cells | Defining evidence |
|---|---:|---|
| **Naive/Bystander** | 2,885 | SELL 12.4, CCR7 6.6, IL7R 14.1, TCF7 5.2; no cytotoxicity, not cycling |
| **Transitional Tem** | 496 | GZMK 8.1, XCL1 3.7; SELL down to 2.3 — losing naive surface phenotype |
| **Activated Stem-like** | 437 | TCF7 10.0 (highest), ID3 5.66, CD27 6.02 **with** CD38 2.81 (10× naive) |
| **TEMRA/Effector** | 426 | GNLY 76.5, FGFBP2 4.8, GZMB 5.2; TCF7 down to 1.3 |
| **Cycling Effector** | 408 | MKI67 3.4, GZMB 13.1; **97% in S/G2M** |
| **MAIT** | 246 | KLRB1 15.3, TRAV1-2 2.7, SLC4A10 1.1, ZBTB16 |
| **Naive (hypoxic)** | 231 | Naive surface intact, but hypoxia z = +2.6 (MIR210HG, DDIT4, P4HA1, BNIP3) |

**Two annotation calls worth flagging to reviewers:**

- **"Activated Stem-like", not Tpex.** This cluster keeps the full stem
  programme (TCF7, ID3, CD27, SELL, CCR7) while showing recent TCR activation
  (CD38 up 10-fold, IL7R down). But **PDCD1 is 0.09 and TOX 0.62** — the
  exhaustion programme Tpex requires is not engaged. CP003's equivalent
  stem-like cluster *does* engage it (ENTPD1 0.98, LAG3 4.07, HAVCR2 1.59) and
  is a genuine Tpex. Calling CP018's cluster Tpex would assert exhaustion that
  is not there. It is also not classic Tscm (IL7R is low) — though low IL7R is
  expected after peptide stimulation, so the resting-state phenotype is an open
  question testable in the unstimulated data.

- **"Naive (hypoxic)" is a real state, not debris.** Naive by surface markers,
  but dominated by a HIF1A-target programme. QC is clean (mito 2.6%, complexity
  0.865, doublet score 0.042), so this is most likely hypoxic stress during
  handling — labelled separately so the signal is not read as biology.

Per-cluster evidence, including technical confounders (depth, doublet score,
heat-shock, hypoxia) and a rule-based check that flags any cluster whose
evidence contradicts its label:
`Integration/Annotation_Evidence/CP018_annotation_justification.csv`
(all nine clusters pass).

---

## 2. Antigen-specific clonotypes

**Definition used here: a clonotype expanded to >2 cells within a single
timepoint, excluding MAIT.** FACT-Seq stimulates with HIV peptide and reads out
which clonotypes proliferated — expansion *is* the assay readout. MAIT cells are
excluded because their semi-invariant TCR makes clonotype sharing a property of
receptor biology rather than antigen response.

This differs from the published CP003 definition, which additionally filtered on
hand-picked cluster numbers. Cluster numbering is not portable between
resolutions or datasets, so specificity is now defined by expansion alone.
**Section 4 re-calls CP003 under this same rule** so the comparison is like for
like.

| | |
|---|---:|
| HIV-specific clonotypes | **72** |
| HIV-specific cells | 311 (99 at 2m, 212 at 90m) |
| **Confirmed in both technical replicates** | **57 / 72 (79%)** |
| MAIT cells excluded | 246 |

The 79% replicate confirmation is the key quality number: these are not
library artefacts.

---

## 3. Persistence — clonotypes ARE maintained at 88 months

**Three clonotypes were detected at both 2m and 90m**, one of them
HIV-specific — **21 cells in total** (18 at 2m, 3 at 90m).

| Clonotype (CDR3α_CDR3β) | 2m | 90m |
|---|---|---|
| CALSDLTGGSYIPTF_CASSYSSWGAGANVLTF (2β) | 15 cells, Cycling Effector | 1 cell, Cycling Effector |
| CAVPMEYGNKLVF_CAWRENTDTQYF | 1 TEMRA/Effector + 1 Cycling Effector | 1 cell, Cycling Effector |
| CAASGSSYKLIF_CASSQLGDEQFF | 1 cell, Cycling Effector | 1 cell, Cycling Effector |

All three have their 90m cell in **Cycling Effector** — the maintained cells are
proliferating, not quiescent. Across both timepoints **20 of the 21 cells are
Cycling Effector**; one clonotype also had a single TEMRA/Effector cell at 2m,
so these clonotypes are not strictly single-state.

Per-cell, per-timepoint breakdown:
`VDJ/TCR/Persistence/tables/CP018_persistent_clonotype_states.csv`

### Why three clonotypes is a positive result, not a negative one

The raw persistence rate is 0.11% of 2m clonotypes — which sounds like near-total
loss. **It is not, and the number must not be reported without this context.**

**Detection ceiling.** Splitting *one* blood draw into two technical replicates
and sequencing both recovers only **3.1% (2m) and 15.2% (90m)** of the same
clonotypes. Same cells, same day. **95.8%** of clonotypes are seen as a single
cell. Cross-timepoint persistence of 0.11% must be read against a within-sample
ceiling of 3–15%, not against 100%. Non-detection at 90m is the *expected*
outcome for any clone that is not large; it is not evidence of clonal loss.

**Scale of what a single detection implies.** We sampled 1,920 TCR+ cells at
90m. A clonotype observed even once therefore sits at a frequency of at least
1/1,920 = 5.2 × 10⁻⁴ of the circulating CD8 pool. At a typical 400 CD8 T
cells/µL over ~5 L of blood (~2 × 10⁹ CD8 cells):

> **A clonotype detected once at 90m represents on the order of
> 10⁶ cells in circulation.**

(Range across plausible CD8 counts: ~5 × 10⁵ at 200/µL to ~2 × 10⁶ at 800/µL.)

So the finding is not "a handful of cells survived." It is that **at least three
clonal lineages — one of them HIV-specific — were still present at
million-cell scale in the blood 7.3 years after ART initiation, and were
still proliferating on antigen re-encounter.** Sparse sampling means three is a
floor on the number of maintained clonotypes, not an estimate of it.

**Expansion drives detectability**, as expected if this is a sampling limit
rather than biology:

| Clone size at 2m | Clonotypes | Re-detected at 90m |
|---|---:|---:|
| 1 cell | 2,590 | 0.04% |
| 2 cells | 18 | 5.6% |
| 3–9 cells | 16 | 0% |
| **10+ cells** | **2** | **50%** |

---

## 4. CP003 vs CP018 — key differences

Both participants' clusters were matched by correlating pseudobulk profiles over
2,189 shared variable genes, so shared labels are earned rather than assumed.

### 4a. The exhausted population is absent in CP018

Every CP003 population is recovered in CP018 — **except Tex**.

| CP003 population | Best CP018 match | r |
|---|---|---:|
| Naive/Bystander | Naive/Bystander | 0.948 |
| Tpex | Activated Stem-like | 0.934 |
| Cycling Effector | Cycling Effector | 0.927 |
| Transitional Tem | Activated Stem-like | 0.919 |
| TEMRA/Effector | TEMRA/Effector | 0.892 |
| **Tex** | *(no comparable cluster)* | **0.695** |

Tex sits **0.20 below every other population**. Supporting this at the
signature level, no CP018 cluster reaches CP003's exhaustion profile:

| | CP003 (never suppressed) | CP018 (suppressed) |
|---|---:|---:|
| Mean exhaustion score | 0.340 | 0.137 |
| Mean stemness score | 1.053 | 1.337 |
| Highest cluster exhaustion | 0.815 | 0.517 |
| Worst exhaustion : stemness ratio | **3.05** | **0.56** |

CP003 has three clusters with an exhaustion:stemness ratio above 1.5 — high
inhibitory receptors *with* collapsed stemness, i.e. terminal exhaustion.
CP018 has none; its maximum is 0.75.

### 4b. Repertoire structure differs sharply

| | CP003 | CP018 |
|---|---:|---:|
| Cells with TCR | 3,436 (85.2%) | 4,645 (90.6%) |
| Clonotypes | 2,012 | **4,295** |
| **Cells in expanded clones** | **47.4%** | **11.4%** |
| Largest clone | **84 cells** | 16 cells |
| Shannon diversity | 6.81 | **8.31** |
| HIV-specific clonotypes (matched rule) | 119 | 72 |

CP003's repertoire is oligoclonal — nearly half its cells sit in expanded
clones, with one clone at 84 cells. CP018 has twice the clonotypes, a quarter
the clonal dominance, a five-fold smaller top clone, and higher diversity.

**Zero clonotypes are shared between the two participants** — expected for
unrelated individuals, and a useful control against index hopping or barcode
bleed-through between datasets.

---

## 5. Implications

**Two independent measurements agree.** Transcriptional (no Tex population,
r = 0.695 vs ≥0.87 for everything else) and clonal (11% vs 47% clonal
dominance) both separate the suppressed participant from the unsuppressed one,
in the direction the Figure 4/5 model predicts: chronic antigen drives
oligoclonal expansion and terminal exhaustion; suppression preserves a diverse
repertoire and a stem-like compartment that has *not* engaged the exhaustion
programme.

**The reservoir-relevant point.** HIV-specific clonotypes are maintained at
million-cell scale for 7.3 years under suppression, and they still proliferate
on antigen re-encounter — all three persistent clonotypes were captured in a
cycling state. Immunological memory against HIV is not lost during suppression.
For cure strategies depending on recall responses, that is the encouraging
direction: the responder population is still there and still functional.

**The open question this raises.** CP018's stem-like compartment retains TCF7,
ID3 and CD27 without engaging PD-1/TOX. Whether these are *bona fide* Tscm at
rest — or only look stem-like because peptide stimulation transiently lowered
IL7R — is directly testable in the unstimulated longitudinal data, and is the
natural next analysis.

---

## 6. Limitations — please read before citing any number

1. **n = 1 per group.** All CP003-vs-CP018 differences are descriptive. The
   per-cell statistics in the output tables describe the observed cell
   populations; they are **not** tests of a participant-level effect and must
   not be reported as such.
2. **Cell Ranger version differs** (CP003 9.0.1, CP018 10.1.0). Versions differ
   in cell calling and V(D)J assembly, so some of the clonotype-count gap may be
   technical. The *dominance* difference (47% vs 11%) is much harder to explain
   this way than the raw counts. This belongs in Methods.
3. **Timepoints are not identical** (CP003 2m/101m, CP018 2m/90m).
4. **Sampling depth bounds every persistence statement.** See section 3.
   Report persistence of expanded clones; do not claim clonal loss from these
   data.
5. **Cluster 2 (Activated Stem-like) carries a depth flag** — mean nFeature
   5,218, about 2× other non-cycling clusters, doublet score 0.236. Lineage
   checks rule out CD4/CD8 and T/myeloid doublets (CD4 = 0.05, normal
   CD3E/CD8A/CD8B), and high transcriptional output is expected in recently
   activated blasts — but it is flagged rather than silently accepted, and
   deserves a second look if it drives a headline result.

---

## 7. Where the numbers live

Every figure ships the table behind it.

| Directory | Contents |
|---|---|
| `Annotation_Plots/Annotated/` | Annotated UMAPs, composition, marker dot plot |
| `Integration/Annotation_Evidence/` | Per-cluster programme scores, QC confounders, justification table, CP003 mapping |
| `VDJ/TCR/Persistence/` | Alluvials, detection-power controls, persistent clonotypes |
| `HIV_Specific/` | Per-cell calls with every decision input, clonotype tables, module scores |
| `CP003_vs_CP018/` | Signature comparison, repertoire comparison, matched recall |

Pipeline: `Scripts/CP018 - FACTseq/`, run via `run_pipeline.sh`.
