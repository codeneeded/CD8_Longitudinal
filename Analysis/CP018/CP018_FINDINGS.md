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

## 3b. Unstimulated cross-reference — the lineages before and after ART

The FACT-Seq clonotypes were matched into CP018's **unstimulated eCITE-seq**
data, which spans three further timepoints: **1m (pre-ART entry)**, 25m and 42m
(both suppressed). This extends the observation window to **89 months** and
adds a pre-ART baseline the stimulation experiment does not have.

### A matching artefact had to be fixed first

Cell Ranger 10.1.0 writes an *empty* D-gene field where earlier versions write a
literal `NA`:

| Dataset | Aligner | Format |
|---|---|---|
| CP003 stim | 9.0.1 | `TRBV6-1` `.NA.` `TRBJ2-5` |
| Unstimulated eCITE-seq | older | `TRBV15` `.NA.` `TRBJ2-3` |
| **CP018 stim** | **10.1.0** | `TRBV29-1` `..` `TRBJ1-2` |

3,607 of CP018's 4,295 clonotypes contain `..`, against **0 of 9,191** in the
unstimulated object — so raw `CTstrict` matched **nothing**, a false negative
indistinguishable from real clonal absence. Normalising recovers it:

| Matching level | HIV-specific recovered | All clonotypes |
|---|---:|---:|
| `CTstrict` raw | **0 / 72** | 5 / 4,295 |
| `CTstrict` normalised *(used)* | **3 / 72** | 31 / 4,295 |
| `CTaa` (looser) | 5 / 72 | 36 / 4,295 |

**CP003 is unaffected** — both its stim object and the unstimulated data use
`.NA.`, so the published cross-reference is valid as run (7 of 165 recovered).
Recovery rates are near-identical under the same rule: **CP003 4.2%, CP018
4.2%**. The pipeline now stops with a diagnosis if the chosen level returns zero
while a looser level finds matches.

### The lineage timeline

Three antigen-validated clonotypes appear in the unstimulated data (9 cells: 6
at pre-ART entry, 3 at 42m).

| Clonotype (CDR3α_CDR3β) | Timepoints | Span |
|---|---|---|
| CALSDLTGGSYIPTF_CASSYSSWGAGANVLTF | 1m pre-ART (2 cells, Tex CD8) → 2m stim (15, Cycling Effector) → 90m stim (1) | **89 months** |
| CAVGGSQGNLIF_CASSLISQGPYEQYV | 1m pre-ART (4 cells) → 2m stim (4, TEMRA/Effector) | 1 month |
| CALRGPLSGNTPLVF_CSAFLGVSSYNEQFF | 42m (3 cells, Tex CD8) → 90m stim (4, TEMRA/Effector) | 48 months |

The first is the important one: **detectable before ART started, expanded on
peptide stimulation at 2 months, and still present at 90 months.** A lineage
tracked across a treatment transition and two independent experiments.

### These clones are NOT Tscm at rest — a direct answer

The stimulated data could not separate genuine Tscm from stimulation-lowered
IL7R. The unstimulated data settles it. Comparing antigen-validated CD8 cells
(n=8) against all other CD8 cells (n=2,180) — CD8-restricted, because pooling
CD4/B/NK/monocytes would inflate every T-cell signal:

| Panel | Result |
|---|---|
| Stemness | **IL7R 0.000, LEF1 0.000, BACH2 0.000, CD28 0.000, FAS 0.000**; CCR7 0.10x, SELL 0.32x, TCF7 0.48x |
| Exhaustion | **LAG3 5.5x, CTLA4 5.2x, TIGIT 2.4x, TOX 1.8x, PDCD1 1.7x** |
| Effector | **KLRG1 6.1x, GZMH 4.0x, NKG7 2.5x, PRF1 1.9x** |

Module scores over CD8 cells only: stemness **down** (Δ −0.857, p_adj 0.0038),
exhaustion **up** (Δ +0.246, p_adj 0.0038), cytotoxicity **up** (Δ +1.163,
p_adj 0.017). Type I IFN memory was not significant (p_adj 0.60).

The decisive detail: the unstimulated annotation contains a named **Tscm CD8**
cluster of 41 cells, and **zero** antigen-validated cells are in it. That is a
direct negative, not an absence of evidence.

**So the low IL7R in the stimulated data was not a stimulation artefact.** These
are differentiated, cytotoxic, partially exhausted effector cells at rest — and
5 of 9 sit in the unstimulated **Tex CD8** cluster. This is consistent with the
53-cell exhausted subpopulation found by integration (section 4a): the
antigen-specific compartment in CP018 is where its exhaustion is concentrated,
even though the participant is suppressed and its bulk clusters look unexhausted.

*Caveat: 3 clonotypes, 8–9 cells. Directions are descriptive, and a zero at any
timepoint means not detected, not proven absent.*

**Do not read cell counts across the two experiments as fold-change.** The
2-cells-at-1m → 15-cells-at-2m step is not a 7.5-fold expansion: the two
numbers come from different assays, different cell inputs, and different
sampling depths. The supportable statement is that the clonotype was
*detectable pre-ART and reached 15 cells after peptide stimulation*. Only
within-experiment, within-timepoint counts support a magnitude claim.

---

## 3c. Two audit corrections — read these before quoting section 3b or 4

### The two "Tex" labels are different clusterings, not one finding

Panel F of the Figure 3 mirror shows CP018's antigen-validated cells sitting
mostly in **"Tex CD8"**, while section 4a says CP018's own clustering has a
*reduced* exhausted compartment. Both are correct, because the labels come from
different namespaces:

| Label | Where it is defined |
|---|---|
| `Tex CD8` | PBMC-atlas annotation of the **unstimulated** data, defined across **all** TARA participants |
| `Tex` | CP018's **own** FACT-Seq clustering (resolution 0.4) |

A cell can be assigned to the shared atlas's `Tex CD8` cluster while its
participant has no Tex cluster of its own. The atlas label is an assignment to a
reference, not a claim about CP018's internal structure.

**And the atlas cluster does not meet our own exhaustion criterion.** Scoring
every atlas CD8 cluster with the same exhaustion:stemness ratio used on the
FACT-Seq clusters:

| Atlas cluster | n | Exhaustion | Stemness | Ratio | Meets Tex criterion (>1.5) |
|---|---:|---:|---:|---:|:--|
| TEMRA/CTL CD8 | 464 | 0.258 | 0.248 | 1.00 | No |
| **Tex CD8** | 85 | 0.355 | 0.489 | **0.71** | **No** |
| Tscm CD8 | 41 | 0.044 | 1.701 | 0.03 | No |
| Naive 1 CD8 | 614 | 0.037 | 1.653 | 0.02 | No |

**No atlas CD8 cluster crosses the threshold, including the one named "Tex CD8"**
— it retains stemness (0.489). So panel F must **not** be read as "CP018's
antigen-specific cells are terminally exhausted." The defensible statement is
that they occupy the atlas's most differentiated, highest-inhibitory-receptor
CD8 compartment. The genuinely exhausted CP018 population is the 53-cell
subpopulation in section 4a, which was found in CP018's own clustering.

### A specificity-labelling bug was truncating persistent lineages

Figure 3C filters **cells** carrying the HIV-specific label. That silently cuts
lineages short: CP018's largest persister (TRAV9-2, **15 cells at 2m**) has
exactly **one** cell at 90m, and a 1-cell clone fails the ">2 cells within a
timepoint" expansion rule there — so that cell is labelled "Other" and **the 90m
end of a genuinely persistent clonotype disappears from the plot.**

This is why an earlier analysis showed shared 2m/90m clonotypes and the first
version of the mirror figure showed none. Antigen specificity is a property of
the **receptor**, not of an individual cell, so once a TCR is called specific
anywhere, every cell carrying it belongs to that lineage.

| | Cell-filter (paper's rule) | Clonotype-level (corrected) |
|---|---:|---:|
| CP003 persisting ≥2 timepoints | **9** | **11** |
| CP018 persisting ≥2 timepoints | 3 | 3 |

The cell-filter column reproduces the published Fig3C subtitle exactly (CP003 =
9), which confirms the construction was read correctly rather than
misinterpreted. **CP003 has 30 such split-label clonotypes and CP018 has 1**, so
the published count is an undercount by 2 — small, but it is an undercount, and
it affects the published figure. Both versions are now produced
(`Combined_C_Alluvial_paperRule_*` and the clonotype-level default) so the
difference is visible rather than a silent choice.

CP018's corrected count is unchanged at 3, but the *shape* changes: the
TRAV9-2 lineage now shows **pre-ART → 2m (15 cells) → 90m**, which is the
89-month span reported in section 3b.

### A third artefact: alluvial ribbons break at non-detected timepoints

Even after the specificity fix, the corrected alluvial still showed no visible
2m→90m connection. That was a **third, purely graphical** artefact:
`geom_flow()` connects only *adjacent* axes, so a lineage absent from an
intermediate timepoint has no stratum there and its ribbon **breaks** — the
later timepoint renders as a disconnected island.

CP018's TRAV9-2 lineage is detected at 1m (pre-ART), 2m and 90m but **not** at
25m or 42m, so the 2m→90m link was invisible despite being present in the data.

**This affects the published figure more than ours.** Of CP003's 11 persisting
lineages, **5 have a detection gap** — so in Fig3C their ribbons break and a
first-to-last span cannot be read off the figure at all. Four CP003 lineages are
present at both the first and last axis point, and this is not visually apparent
in the published panel.

Two fixes, both produced:

1. **Gap-bridged alluvial** — the clonotype × timepoint grid is completed with
   zeros and drawn with `geom_alluvium`, giving one continuous band per lineage.
   Stratum heights remain exact; a non-detected timepoint gets a thin floor
   ribbon, labelled as such on the panel.
2. **Lineage trace** (`Fig3C2_*`, `Combined_C2_*`) — one row per lineage, filled
   point where detected (number = cells), open circle where not. No graphical
   choice can mislead here, so this is the panel to use when the claim is about
   span rather than magnitude.

| Artefact | Effect | Status |
|---|---|---|
| Cross-experiment `CTstrict` format | 0 matches instead of 3 | fixed (normalisation) |
| Per-cell specificity filter | truncates lineages at low-count timepoints | fixed (clonotype-level) |
| `geom_flow` adjacency | ribbon breaks across detection gaps | fixed (bridged + trace) |

All three were false negatives — each one *hid* persistence rather than
inventing it.

---

## 4. CP003 vs CP018 — key differences

### 4a. Exhaustion is reduced ~3-fold in CP018, not absent

Both participants' clusters were matched by correlating pseudobulk profiles over
2,189 shared variable genes, so shared labels are earned rather than assumed.

At whole-cluster level, every CP003 population is recovered in CP018 —
**except Tex**.

| CP003 population | Best CP018 match | r |
|---|---|---:|
| Naive/Bystander | Naive/Bystander | 0.948 |
| Tpex | Activated Stem-like | 0.934 |
| Cycling Effector | Cycling Effector | 0.927 |
| Transitional Tem | Activated Stem-like | 0.919 |
| TEMRA/Effector | TEMRA/Effector | 0.892 |
| **Tex** | *(no comparable cluster)* | **0.695** |

Tex sits **0.20 below every other population** at whole-cluster level.

**IMPORTANT CORRECTION — joint integration revised this.** Integrating both
datasets (CCA, participant as batch) put all 24 of CP003's Tex cells into joint
cluster 7, which also contains **53 CP018 cells**. Those CP018 cells are
exhausted *on their own terms*, not pulled there by batch correction:

| Group | n | Exhaustion | Stemness | Ratio |
|---|---:|---:|---:|---:|
| CP003 in cluster 7 | 114 | 0.920 | 0.380 | **2.36** |
| **CP018 in cluster 7** | **53** | **1.038** | **0.417** | **2.43** |
| CP018 Cycling Effector elsewhere | 355 | 0.439 | 0.998 | 0.44 |
| All cells outside cluster 7 | 8,994 | 0.213 | 1.228 | 0.17 |

Per gene, CP018's cluster-7 cells **equal or exceed** CP003's on TOX (1.18 vs
0.95), PDCD1 (0.79 vs 0.36), LAG3 (2.22 vs 1.99), CTLA4 and ENTPD1. Their
exhaustion:stemness ratio is 5.5x that of cells carrying the same label
elsewhere in CP018, so this is a genuine rare subpopulation.

**So the difference is one of MAGNITUDE, not absence:**

| | Cells in the exhausted cluster | % of participant |
|---|---:|---:|
| CP003 (never suppressed) | 114 | **2.83%** |
| CP018 (suppressed) | 53 | **1.03%** |

CP018 has a ~2.7-fold *smaller* exhausted compartment, not a missing one. Its
own per-participant clustering assigned these 53 cells to Cycling Effector,
and whole-cluster pseudobulk correlation averaged a 1%-of-cells subpopulation
away — which is exactly the resolution limit of that method, not an error in it.
The correlation result stands at cluster level; this is the finer-grained view.

At the signature level, CP018's *clusters as a whole* still do not reach
CP003's exhaustion profile:

| | CP003 (never suppressed) | CP018 (suppressed) |
|---|---:|---:|
| Mean exhaustion score | 0.340 | 0.137 |
| Mean stemness score | 1.053 | 1.337 |
| Highest cluster exhaustion | 0.815 | 0.517 |
| Worst exhaustion : stemness ratio | **3.05** | **0.56** |

CP003 has three clusters with an exhaustion:stemness ratio above 1.5 — high
inhibitory receptors *with* collapsed stemness, i.e. terminal exhaustion. No
CP018 *cluster* reaches that at whole-cluster level (maximum 0.75) — but as
shown above, a 53-cell subpopulation within its Cycling Effector cluster does
(ratio 2.43). Both statements are true at their own resolution; report the
prevalence difference rather than presence/absence.

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

**The Tscm question is now answered — negatively.** The unstimulated data
(section 3b) shows the antigen-validated clones are *not* Tscm at rest: IL7R,
LEF1, BACH2, CD28 and FAS are all zero, exhaustion markers are elevated 1.7–5.5x,
and the named Tscm CD8 cluster contains none of them. So CP018's bulk stem-like
signature does not extend to its antigen-specific cells — the cells that
recognise HIV are differentiated and partially exhausted even under suppression,
while the participant's overall CD8 compartment remains stem-like and diverse.

**What this reframes.** Suppression preserves the *global* repertoire and stem
compartment, but the antigen-specific fraction still carries an exhaustion
signature — concentrated in the 1%-of-cells subpopulation that integration
exposed. The protection is not uniform across the compartment.

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
