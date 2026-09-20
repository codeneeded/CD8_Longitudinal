# Revision — build the revised-figures PDF that the response letter cites
#
# WHAT THIS IS
#   One PDF, ONE FIGURE (or table) PER PAGE, so the point-by-point response can
#   say "see revised figure PDF, page N" instead of asking the reviewer to hunt
#   through the manuscript. Figures themselves stay PNG only (per project
#   convention); this PDF just embeds them.
#
#   Page 1 is a contents page listing every item with its page number, and the
#   page ORDER IS FIXED by the ITEMS list below so a citation cannot drift.
#
#   The manifest CSV it writes is the source of truth for the page numbers used
#   in Reviewer_Comments/Response_Text.csv — never type a page number by hand.
#
# Output: Deliverables/.../TARA_Revised_Figures.pdf + Revised_figure_pages.csv
################################################################################

.cfg <- Sys.getenv("CD8_REV_CONFIG", "")
if (!nzchar(.cfg)) {
  for (b in c("~/Documents/CD8_Longitudinal",
              "~/Desktop/Projects/Repositories/CD8_Longitudinal")) {
    p <- file.path(path.expand(b), "Scripts", "4. Revision", "config_revision.R")
    if (file.exists(p)) { .cfg <- p; break }
  }
}
source(.cfg)

suppressPackageStartupMessages({library(png); library(grid); library(gridExtra)})

OUTDIR <- path.expand("~/Desktop/Projects/Deliverables/CP018/collaborator_package")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUTDIR, "TARA_Revised_Figures.pdf")

# ── The page order. One entry = one page. ────────────────────────────────────
# kind: "fig" (a PNG) or "tab" (a CSV rendered as a table)
ITEMS <- list(
  list(kind = "fig", file = "R1m6_03_DSB_distributions.png",
       id = "R1-m6",
       title = "DSB-normalised protein value distributions, per marker",
       caption = paste0(
         "Answers Reviewer 1, minor point 6. Distribution of DSB-normalised values over cells ",
         "positive for each marker, with all seven isotype controls shown across all cells as ",
         "the internal negative control. Dashed lines mark the 5-20 range the reviewer expects ",
         "of positive populations. Proposed as a new Supplementary Figure S1B. The -1 to 2 range ",
         "in the existing S1A is a property of that panel, which is a Seurat DotPlot: it z-scores ",
         "each feature across groups and clips the colour scale to +/-2.5 SD, so its scale was ",
         "never in DSB units.")),
  list(kind = "tab", file = "R1m6_01_DSB_marker_distributions_positive_cells.csv",
       id = "R1-m6",
       title = "DSB values: lineage markers versus isotype controls",
       caption = paste0(
         "Source data for the preceding figure. Isotype controls are detected above zero in ",
         "roughly 1% of cells at medians of 0.20-0.32, while lineage markers are positive in ",
         "40-98% of their populations at medians up to 17.44 - confirming that isotype-based ",
         "denoising centred the background rather than compressing real signal.")),
  list(kind = "tab", file = "R1m6_02_DSB_panel_summary.csv",
       id = "R1-m6",
       title = "Panel-wide DSB summary",
       caption = paste0(
         "Of 137 antibodies, 61 have a 99th-percentile value above 5, 29 above 10 and 4 above ",
         "20. The stored minimum of 0 and the 81.7% exact zeros reflect sparse storage of the ",
         "normalised matrix, which represents DSB's below-background negatives as zero; this is ",
         "now disclosed in the Methods and affects no reported result, as every protein gate ",
         "used is a positive threshold.")),
  list(kind = "tab", file = "R1m1_01_CD8_central_memory_check.csv",
       id = "R1-m1",
       title = "Central-memory CD8 cells in the dedicated CD8 object",
       caption = paste0(
         "Answers Reviewer 1, minor point 1. CCR7 is absent from the TotalSeq-C panel, so Tcm ",
         "cells cannot be gated on surface CCR7. Scoring the transcriptional profile instead ",
         "identifies 450 of 18,843 CD8 cells (2.4%), spread across 9 clusters with no cluster ",
         "exceeding 7.2% - i.e. central-memory cells are present as a minor profile within the ",
         "naive and effector-memory compartments rather than as a separable population. ",
         "Proposed as Supplementary Table S3.")),
  list(kind = "fig", file = "R2M1d_03_ADT_subgating_crosswalk.png",
       id = "R2-M1d",
       title = "Surface-protein sub-gating of the naive CD8 clusters",
       caption = paste0(
         "Answers Reviewer 2, major point 1d, and corrects a mislabelling in our own earlier ",
         "draft of this panel. The dedicated CD8 clustering was performed on the CD8 subset ",
         "(re-normalisation, separate RNA and ADT integrations, WNN graph, resolution 0.4); this ",
         "panel shows the subsequent surface-protein sub-gating step, which redistributes the ",
         "four naive WNN clusters into Naive Intermediate (CD45RO+) and Tscm (CD45RO-/FAS+).")),
  list(kind = "fig", file = "R2m13_02_CP003_LineageTrace.png",
       id = "R2-m13",
       title = "CP003 antigen-specific lineages, one row per clonotype",
       caption = paste0(
         "Answers Reviewer 2, minor point 13. One row per clonotype, filled markers where the ",
         "clonotype was detected and labelled with the number of cells, open markers where it ",
         "was not. Four clonotypes are detected at both 2 and 101 months, and all four are ",
         "undetected at 1, 12 and 24 months. That is why the published panel shows no ",
         "2-to-101-month connection: an alluvial geometry links only adjacent axes, so a ",
         "clonotype absent at every intermediate timepoint has no ribbon to carry across. The ",
         "persistence is in the data; the published rendering cannot show it.")),
  list(kind = "fig", file = "R2m13_01_CP003_Alluvial_bridged.png",
       id = "R2-m13",
       title = "CP003 clonotype alluvial with detection gaps bridged",
       caption = paste0(
         "The same 11 clonotypes in alluvial form, with the clonotype-by-timepoint grid ",
         "completed so a ribbon spans a non-detected timepoint as a thin bridge instead of ",
         "breaking. Per-timepoint markers are omitted because several clonotypes share a ",
         "timepoint and their markers cannot be attributed; the preceding trace panel carries ",
         "that per-clonotype detail.")),
  list(kind = "tab", file = "R2m13_01_CP003_alluvial_source.csv",
       id = "R2-m13",
       title = "CP003 clonotype detections behind both panels",
       caption = paste0(
         "One row per clonotype and timepoint with its cell count, so every marker and ribbon ",
         "in the two panels above can be checked directly.")),
  # The editor's additional-data point: CP018 FACT-Seq, generated after
  # submission as an independent second participant, run in technical duplicate.
  list(kind = "fig", file = "E1_02_CP018_LineageTrace.png",
       id = "E-1",
       title = "CP018 antigen-specific lineages, one row per clonotype",
       caption = paste0(
         "Additional data generated after submission. CP018 is a second participant, processed ",
         "with the same FACT-Seq protocol in technical duplicate at both timepoints, and its ",
         "unstimulated course adds a pre-treatment baseline the published participant's design ",
         "does not have. Same construction as the CP003 panel.")),
  list(kind = "fig", file = "E1_03_CP018_Alluvial.png",
       id = "E-1",
       title = "CP018 clonotype alluvial with detection gaps bridged",
       caption = "The CP018 clonotypes in alluvial form, constructed as for CP003."),
  list(kind = "fig", file = "E1_01_Combined_C2_LineageTrace_CP003_vs_CP018.png",
       id = "E-1",
       title = "Both participants on a shared layout",
       caption = paste0(
         "CP003 above, CP018 below. The published participant's timepoint axis is not ",
         "chronological because the two stimulated draws bracket the unstimulated course; ",
         "CP018's runs in elapsed time and starts before treatment. Shown so the two ",
         "participants are compared on construction, not on axis order.")),
  # Panels for the asks that previously said "requires re-rendering". They are
  # OPTIONS: none is applied to the manuscript, and each response states what
  # the manuscript change would be if the panel is adopted.
  list(kind = "fig", file = "R2M1b_03_Fig1C_bubble_ADT.png",
       id = "R2-M1b",
       title = "Figure 1C as a bubble plot: surface protein",
       caption = paste0(
         "Answers Reviewer 2, major point 1b, for the panel the comment names. Same markers, ",
         "same six whole-PBMC clusters and same marker grouping as the published Figure 1C, ",
         "drawn so each point carries both quantities the comment asks for: area is the ",
         "percentage of cells expressing, colour is the mean scaled across clusters. Protein ",
         "is DSB-normalised and stored with a floor at zero, so expression above zero is what ",
         "the area counts; that is why a marker can be high in level and modest in area.")),
  list(kind = "fig", file = "R2M1b_04_Fig1C_bubble_RNA.png",
       id = "R2-M1b",
       title = "Figure 1C as a bubble plot: transcript",
       caption = paste0(
         "The transcript half of the same panel, constructed identically. Dot area is the ",
         "percentage of cells with a detected transcript, so dropout reduces area without ",
         "moving colour.")),
  list(kind = "tab", file = "R2M1b_03_Fig1C_bubble_source.csv",
       id = "R2-M1b",
       title = "Values behind both bubble panels",
       caption = paste0(
         "One row per feature, cluster and modality, with percentage expressing, mean and ",
         "scaled mean, so every point can be checked.")),
  list(kind = "fig", file = "R2M1e_02_Fig1E_pct_of_CD8.png",
       id = "R2-M1e",
       title = "Cluster frequency as a percentage of each sample's CD8 compartment",
       caption = paste0(
         "Answers Reviewer 2, major point 1e, for the manuscript's own clusters. Each point ",
         "is one sample and the denominator is that sample's own CD8 count, so the value no ",
         "longer depends on how large the CD8 fraction was. Both denominators are in the ",
         "source table. No pairwise comparison between exposure groups reaches ",
         "Benjamini-Hochberg adjusted p < 0.05.")),
  list(kind = "tab", file = "R2M1e_02_cluster_frequency_both_denominators.csv",
       id = "R2-M1e",
       title = "Cluster frequency under both denominators",
       caption = paste0(
         "Per sample and cluster: cell count, the sample's CD8 count, its PBMC count, and the ",
         "percentage under each denominator.")),
  list(kind = "fig", file = "R2M2a_02_CD8only_clonotype_UMAP.png",
       id = "R2-M2a",
       title = "Clonotype expansion on the CD8-only embedding",
       caption = paste0(
         "Answers Reviewer 2, major point 2a. Clonotype size drawn on the dedicated CD8 ",
         "embedding instead of the whole-PBMC one. One limit is stated on the panel: the only ",
         "CD8-only embedding in this study was built on HEI participants, so this cannot ",
         "carry the HEU and HUU facets the published Figure 2A has - producing those would ",
         "mean re-integrating CD8 cells across all three groups, a new analysis rather than a ",
         "re-render. Faceting is by treatment stage, the comparison these cells support.")),
  list(kind = "fig", file = "R2m2_02_UMAP_counts_on_panel.png",
       id = "R2-m2",
       title = "Cell numbers annotated on the UMAP panels",
       caption = paste0(
         "Answers Reviewer 2, minor point 2. Each population is labelled with its own cell ",
         "count on the panel and the total sits in the corner, so the number does not have to ",
         "be carried in the legend. Whole-PBMC embedding on the left, dedicated CD8 on the ",
         "right.")),
  list(kind = "tab", file = "R1M2_01_Lineage_annotation_findings.csv",
       id = "R1-M2",
       title = "Cluster counts and annotation calls, per lineage",
       caption = paste0(
         "Answers Reviewer 1, major point 2. Each non-CD8 lineage was subset from the published ",
         "object and re-analysed with the same method used for CD8 (re-normalisation, separate ",
         "RNA and ADT FastMNN integrations, WNN graph, resolution sweep). Clusters whose ",
         "evidence does not support a confident label are reported as ambiguous or unresolved ",
         "rather than being given a forced name.")),
  list(kind = "tab", file = "R1m4_01_CD3_detection_by_compartment.csv",
       id = "R1-m4",
       title = "CD3 transcript versus surface protein, by compartment",
       caption = paste0(
         "Answers Reviewer 1, minor point 4 (CD8 Tn cells clustering with CD4 T cells). CD3 ",
         "transcript is detected in a large minority of NK cells with essentially no ",
         "corresponding surface protein, whereas true T cells show both. Surface protein, not ",
         "transcript, is therefore the arbiter of lineage assignment - with CD4 and CD8 ",
         "compartments included as positive controls.")),
  list(kind = "tab", file = "R2M1a_01_CD8_cell_count_reconciliation.csv",
       id = "R2-M1a",
       title = "CD8 cell counts at each analysis stage",
       caption = paste0(
         "Answers Reviewer 2, major point 1a. An exhaustive partition by cell origin summing ",
         "exactly to the 18,843 cells in the dedicated CD8 object. Rows D-G partition the final ",
         "object; the residual is decomposed into marker-based gd/MAIT removal and cells ",
         "dropped with contaminant clusters and subset QC. Note that removal was marker-based ",
         "rather than label-based, so subtracting annotation labels does not reproduce the total.")),
  list(kind = "tab", file = "R2M1d_01_CD8_annotation_justification.csv",
       id = "R2-M1d",
       title = "CD8 cluster annotation: evidence per population",
       caption = paste0(
         "Marker-programme scores behind each CD8 label, including technical and cross-lineage ",
         "panels, with a rule-based verdict flagging any population whose evidence conflicts ",
         "with its assigned name. Proposed as a supplementary table supporting the annotation.")),
  list(kind = "tab", file = "R2M1d_02_CD8_flagged_population_resolution.csv",
       id = "R2-M1d",
       title = "Identity evidence for CD16+ Effector and Naive CD8 4",
       caption = paste0(
         "The two populations the evidence check flagged, each resolved with specific marker ",
         "evidence rather than reassigned by assumption. Included because a reviewer repeating ",
         "this check would raise the same two.")),
  list(kind = "fig", file = "R1m3_01_TEMRA.png",
       id = "R1-m3",
       title = "Figure 1E: TEMRA/CTL frequency by HIV exposure group",
       caption = paste0(
         "Answers Reviewer 1 minor 3 and Reviewer 2 minor 4. The y-axis title ",
         "'% of total PBMC' was cropped because a 16 pt bold rotated title does not fit ",
         "the 5-inch canvas the panel was saved at. A left plot margin and a slightly ",
         "wider canvas resolve it; the data and the statistics are unchanged.")),
  list(kind = "fig", file = "R2m8_01_ClonalSharing.png",
       id = "R2-m8",
       title = "Supplementary 2C: clonotypes shared across CD8 clusters",
       caption = paste0(
         "Answers Reviewer 2, minor point 8. Only the FIRST bar's count label was ",
         "missing. The y-axis had been set with no upper expansion and a limit equal to ",
         "the tallest bar, so that bar's label alone fell outside the panel while the ",
         "other two rendered normally. With 10% upper headroom all three render: 654 ",
         "clonotypes confined to a single cluster, 148 spanning two, 19 spanning three.")),
  list(kind = "fig", file = "R2m16_01_ADT_gating.png",
       id = "R2-m16",
       title = "Supplementary 4C: surface-protein sub-gating of the naive CD8 clusters",
       caption = paste0(
         "Answers Reviewer 2, minor point 16. The facet strips previously read ",
         "'Source cluster 0 / 2 / 4 / 6'. Cluster numbers appear nowhere else in the ",
         "manuscript text, so each facet is now labelled with the naive subset name(s) ",
         "that source cluster resolves to, with the number retained in parentheses for ",
         "traceability to the underlying clustering. As in the published panel, only ",
         "the 2,649 gated cells (Tscm and Intermediate) of 11,245 naive cells are ",
         "plotted; point size, font and colours are unchanged.")),
  list(kind = "fig", file = "R2M2a_01_CD8_dedicated_UMAP.png",
       id = "R2-M2a",
       title = "CD8-only UMAP from the dedicated CD8 clustering",
       caption = "Answers Rev 2 major 2a. The CD8-only embedding the reviewer asks for already exists as the dedicated CD8 re-clustering, shown here. Adopting it for the Figure 2 clonal panels would put the clonal and phenotypic views in one coordinate system; that change is pending the decision on figure layout."),
  list(kind = "fig", file = "R2M1b_01_bubble_RNA.png",
       id = "R2-M1b",
       title = "Bubble plot, transcript: fraction expressing and level",
       caption = "Answers Rev 2 major 1b. Point size is the fraction of cells expressing the gene and colour the mean level in those cells, separating the two quantities a heatmap of scaled means conflates. Produced for the dedicated CD8 clustering during this revision."),
  list(kind = "fig", file = "R2M1b_02_bubble_ADT.png",
       id = "R2-M1b",
       title = "Bubble plot, surface protein: fraction expressing and level",
       caption = "Answers Rev 2 major 1b, protein counterpart. The same encoding on DSB-normalised surface protein, so transcript and protein can be read side by side."),
  list(kind = "fig", file = "R2M1e_01_composition_pctCD8.png",
       id = "R2-M1e",
       title = "Cluster frequency as a percentage of total CD8",
       caption = "Answers Rev 2 major 1e. Each cluster as a percentage of that sample's CD8 compartment rather than of live cells, which removes the dependence on how much of the sample was CD8."),
  list(kind = "fig", file = "R2M1d_04_naive2_identity.png",
       id = "R2-M1d",
       title = "Naive CD8 2: identity in the dedicated clustering",
       caption = "Answers Rev 2 major 1d. Naive CD8 2 persists as a distinct population in the dedicated clustering and is unambiguously CD8 by both modalities."),
  list(kind = "fig", file = "R1M2_02_lineage_identity_RNA_ADT.png",
       id = "R2-M1d",
       title = "Lineage identity by transcript and surface protein",
       caption = "Answers Rev 1 major 2. Every lineage in the dataset was subset, re-clustered and annotated for this revision; this panel shows the identity evidence per population at both transcript and protein level."),
  list(kind = "fig", file = "R1M2_03_qc_vs_CD8.png",
       id = "R2-M1d",
       title = "Quality metrics against CD8 signal, per population",
       caption = "Answers Rev 1 major 2. Quality metrics per population against CD8 signal, so a population is not read as biology when it is a quality artefact."),
  list(kind = "fig", file = "R1M2_04_crosslineage_doublet.png",
       id = "R2-M1d",
       title = "Cross-lineage and doublet checks per population",
       caption = "Answers Rev 1 major 2. Cross-lineage marker detection and doublet score per population. A flag requires both an elevated module score and real transcript detection, because a relative score alone ranks some population top at everything."),
  list(kind = "fig", file = "R1M2_05_participant_composition.png",
       id = "R2-M1d",
       title = "Population composition by participant",
       caption = "Answers Rev 1 major 2. Composition by participant, showing which populations are shared across the cohort and which are participant-restricted."),
  list(kind = "fig", file = "R1M2_06_programme_heatmap.png",
       id = "R2-M1d",
       title = "Programme scores per population",
       caption = "Answers Rev 1 major 2. Transcriptional programme scores per population across all re-clustered lineages."),
  list(kind = "fig", file = "R2m2_01_cells_per_sample.png",
       id = "R2-m2",
       title = "Cells recovered per sample, all lineages",
       caption = "Answers Rev 2 minor 2 and the editor's point 2. Cells recovered per sample by lineage, the per-sample view behind the totals now stated in the Figure 1 and Figure 4 legends. Full counts are in Supplementary Table S1."),
  list(kind = "fig", file = "R2M1c_01_all_CD8_clusters_by_group.png",
       id = "R2-M1c",
       title = "CD8 population frequencies by HIV exposure group",
       caption = "Answers Rev 2 major 1c. All seven CD8 populations in the entry-visit cohort, on the same denominator and the same pairwise Wilcoxon test as the published Figure 1E. TEMRA/CTL reproduces the published values exactly (HEI vs HEU p = 0.022, HEI vs HUU p = 0.026). Benjamini-Hochberg adjustment across the seven populations is given in the source table; no population reaches q < 0.05, which we would rather show than omit."),
  list(kind = "fig", file = "R2m15_01_S_CD8_composition_alluvial.png",
       id = "R2-m15",
       title = "CD8 subset composition across the sampled course",
       caption = "Answers Rev 2 minor 15. Each band is one CD8 population's share of the compartment at that treatment stage, averaged across participants so a single deeply sequenced sample cannot set the ribbon shape. Timepoint groups are treatment stages rather than equally spaced times, and participants contribute to different stages, so the panel should be read as composition per stage and not as a rate of change."),
    list(kind = "fig", file = "E5_01_sex_inference.png",
         id = "E-5",
         title = "Genetic sex inferred from sex-chromosome expression",
         caption = paste0(
           "Answers the editor's point 5. Sex was absent from the clinical dataset, so it was ",
           "inferred from sex-chromosome expression: XIST against eight Y-linked genes, one ",
           "point per library. The two groups separate completely with nothing in between - ",
           "female libraries detect Y genes in exactly 0% of cells, male libraries detect XIST ",
           "in exactly 0% - and all 11 participants sampled more than once gave the same call ",
           "at every timepoint. The cohort is 12 female and 4 male. This establishes genetic ",
           "sex only; gender was not collected and is not inferable from transcriptomes.")),
    list(kind = "fig", file = "E5_02_composition_by_sex.png",
         id = "E-5",
         title = "CD8 composition by inferred sex, HEI participants",
         caption = paste0(
           "Answers the editor's point 5. Descriptive only, and deliberately so: the sex by ",
           "exposure-group table has an empty cell (HUU male, n = 0) and a cell of one ",
           "participant (HEU male), so a sex-stratified comparison cannot be estimated in 16 ",
           "participants. Every participant is shown as a point rather than summarised, ",
           "because with 8 female and 3 male HEI participants a bar of medians would imply ",
           "precision that is not present. No test is reported.")),
    list(kind = "fig", file = "R2M5c_01_CD16_signature_match.png",
         id = "R2-M5c",
         title = "CD16+ effector CD8 against the congenital-CMV signature",
         caption = paste0(
           "Answers Rev 2 major 5c. Our CD16+ effector CD8 population (860 cells, 4.6% of ",
           "CD8) tested against the phenotype reported for FcgammaRIII+ CD8+ T cells in ",
           "congenital CMV infection (Semmes et al., J Clin Invest 2024, PMID 39531313), ",
           "feature by feature. It is consistent on all 11 testable features, including the ",
           "two transcription-factor directions that paper emphasises: T-bet higher (46.4 vs ",
           "15.8% of cells) and BCL11B lower (45.1 vs 59.9%). A phenotypic match does not ",
           "establish CMV causation: CMV serostatus was never measured in this cohort.")),
    list(kind = "tab", file = "R2M5c_02_CD16_signature_comparison.csv",
         id = "R2-M5c",
         title = "Feature-by-feature comparison with PMID 39531313",
         caption = paste0(
           "Answers Rev 2 major 5c. Source data for the panel above, with the direction each ",
           "feature is reported to take in the published population and whether ours agrees.")),
    list(kind = "fig", file = "R2M3d_01_CMV_enrichment.png",
         id = "R2-M3d",
         title = "CMV-annotated clonotypes by CD8 population",
         caption = paste0(
           "Answers Rev 2 major 3d. Database-predicted CMV-specific clonotypes are not ",
           "enriched in the CD16+ effector population but depleted (4.65% vs 10.88% of ",
           "TCR-recovered cells; Fisher OR 0.39, p = 5.2e-11). The obvious artefact was ",
           "tested: epitope databases are built from alpha-beta TCRs and this population is ",
           "37% gamma-delta, but the depletion persists among alpha-beta cells alone (OR ",
           "0.47, p = 1.5e-05), and TEMRA CD8 is equally gamma-delta-rich (37.1%) without ",
           "being depleted. We read this as consistent with that population's described ",
           "mechanism being antibody-dependent and TCR-independent, not as evidence against ",
           "a CMV contribution - database annotation covers public clonotypes only, and CMV ",
           "serostatus was never measured here.")),
    list(kind = "tab", file = "R2M3d_04_CMV_enrichment_gd_controlled.csv",
         id = "R2-M3d",
         title = "CMV enrichment before and after excluding gamma-delta cells",
         caption = paste0(
           "Answers Rev 2 major 3d. Both analyses are shown together; the depletion is ",
           "not an artefact of database blindness to gamma-delta TCRs.")),
    list(kind = "fig", file = "R2M2b_02_clonotype_frequency.png",
         id = "R2-M2b",
         title = "Cells in expanded clonotypes, as a percentage of TCR-recovered cells",
         caption = paste0(
           "Answers Rev 2 major 2b. Cells in expanded clonotypes as a percentage of cells ",
           "with a recovered TCR. The reviewer's concern is well founded in this dataset: ",
           "TCR-recovered cells range 262-1438 across samples, a 5.5-fold spread, so a count ",
           "partly tracks recovery depth rather than repertoire structure, and the count and ",
           "the frequency do not rank the samples in the same order. The denominator is ",
           "TCR-recovered cells rather than all cells, because dividing by all cells would ",
           "rescale each timepoint by its TCR recovery rate and reintroduce the artefact ",
           "being removed.")),
    list(kind = "tab", file = "R2M2b_01_clonotype_frequency_by_timepoint.csv",
         id = "R2-M2b",
         title = "Clonotype frequency per sample",
         caption = paste0(
           "Answers Rev 2 major 2b. Per-sample counts and frequencies with the TCR-recovered ",
           "denominator shown, so the count-versus-frequency difference is auditable.")),
    list(kind = "fig", file = "R1M2_07_lineage_pathway_heatmap.png",
         id = "R1-M2",
         title = "Shared transcriptional programmes across all annotated lineages",
         caption = paste0(
           "Answers Rev 1 major 2. Pathway enrichment already run for all 30 annotated ",
           "populations across five contrasts, collected into one panel; no new statistics. ",
           "Infection and viraemia give different answers. HEI vs HEU before treatment: ",
           "complement is the most widespread term (15 populations, five families, ",
           "p_adj 2.5e-14), with IL-6/JAK/STAT3, unfolded protein response and PI3K/AKT/mTOR; ",
           "TNF-alpha/NF-kB is significant in no population. Suppressed vs unsuppressed on ",
           "treatment: TNF-alpha/NF-kB is dominant (14 populations, all six families, ",
           "4.9e-29), with 254 of 334 significant gene-by-population instances higher in the ",
           "unsuppressed participant. Ribosomal and translation terms are significant in 26 of ",
           "30 populations and are excluded as non-discriminating.")),
    list(kind = "tab", file = "R1M2_08_lineage_pathway_drivers.csv",
         id = "R1-M2",
         title = "NF-kB driver genes by number of populations",
         caption = paste0(
         "The 70 TNF-alpha/NF-kB genes that reach significance and are higher in the ",
         "unsuppressed participant in at least one population, ranked by how many populations ",
         "show that. The Hallmark term matched 90 genes in total across the 14 populations; the ",
         "20 not listed here reach significance only in the opposite direction. SOCS3 is the ",
         "most recurrent at 16 populations, followed by IL7R and BTG2 at 10.")),
    list(kind = "fig", file = "R1M2_11_nonCD8_lineage_brief.png",
         id = "R1-M2",
         title = "Differential expression outside the CD8 compartment",
         caption = paste0(
           "Answers Rev 1 major 2 for the non-CD8 lineages. Significant ",
           "gene-by-population calls (MAST, Benjamini-Hochberg p < 0.05) summed over the ",
           "populations in each lineage, from analysis already run. CD4 carries the largest ",
           "burden; pDC did not meet the minimum cells-per-group threshold in any contrast and ",
           "was not assessed. The dot gives our assessment of whether the lineage holds ",
           "something worth adding to the paper.")),
    list(kind = "tab", file = "R1M2_10_nonCD8_lineage_brief.csv",
         id = "R1-M2",
         title = "Non-CD8 lineages: what was found and whether it is worth pursuing",
         caption = paste0(
           "Answers Rev 1 major 2. One row per lineage outside CD8, with the ",
           "differential-expression burden in each contrast, the Hallmark terms significant in ",
           "that lineage and no other, and an assessment. Offered for the authors to decide ",
           "what belongs in the paper; nothing from it has been added to the manuscript.")),
    list(kind = "fig", file = "R1M3_01_HIV_cohort_detection.png",
         id = "R1-M3",
         title = "HIV-1 transcript detection at the entry visit",
         caption = paste0(
           "Answers Rev 1 major 3. Nine entry-visit libraries realigned to a combined ",
           "GRCh38 + HIV-1 reference carrying subtype B (HXB2) and subtype C contigs. ",
           "Bars give HIV-1 UMIs per library, split into UMIs recovered from called cells ",
           "(cell-associated) and from empty droplets (ambient). Three libraries carry ",
           "cell-associated vRNA above the ambient floor; the remainder are negative at ",
           "this detection ceiling. Per-cell burden, lineage and the three controls are in ",
           "HIV_cohort_positive_cells.csv and HIV_cohort_controls.csv.")),
  list(kind = "tab", file = "HIV_cohort_per_sample.csv",
       id = "R1-M3",
       title = "HIV-1 transcript detection per library",
       caption = paste0(
         "Nine entry-visit libraries plus the pilot, realigned with Cell Ranger 10.1.0 to a ",
         "combined GRCh38 + HIV-1 reference carrying the subtype B (HXB2, K03455) and subtype ",
         "C contigs. Per library: cells called, HIV-1 UMIs recovered from called cells, UMIs ",
         "in all barcodes, the ambient fraction, vRNA-positive cells, and the split between ",
         "the two viral references."))
)

# ── Rendering helpers ────────────────────────────────────────────────────────
FIGDIR <- REV_FIG
TABDIR <- REV_TAB

wrap_lines <- function(txt, width) {
  # strwrap() drops embedded newlines, so wrap each existing line separately.
  unlist(lapply(strsplit(txt, "\n", fixed = TRUE)[[1]],
                function(l) if (nzchar(l)) strwrap(l, width = width) else ""))
}

page_no <- 0
manifest <- data.frame()

header <- function(id, title) {
  grid.text(paste0("Revised figure PDF  |  ", id), x = 0.06, y = 0.965,
            just = c("left", "top"), gp = gpar(fontsize = 9.5, col = "grey40"))
  grid.text(paste(wrap_lines(title, 78), collapse = "\n"), x = 0.06, y = 0.935,
            just = c("left", "top"), gp = gpar(fontsize = 15, fontface = "bold"))
}

footer <- function(n) {
  grid.text(sprintf("page %d", n), x = 0.94, y = 0.028, just = c("right", "bottom"),
            gp = gpar(fontsize = 9, col = "grey45"))
  grid.lines(x = c(0.06, 0.94), y = c(0.055, 0.055), gp = gpar(col = "grey85"))
}

caption_block <- function(txt, y_top) {
  lines <- wrap_lines(txt, 116)
  grid.text(paste(lines, collapse = "\n"), x = 0.06, y = y_top,
            just = c("left", "top"), gp = gpar(fontsize = 10, lineheight = 1.35))
  length(lines)
}

render_fig <- function(it, n) {
  path <- file.path(FIGDIR, it$file)
  grid.newpage()
  # The PNG carries its own title; repeating it here reads as a duplicate.
  grid.text(paste0("Revised figure PDF  |  ", it$id, "  |  ", it$title),
            x = 0.06, y = 0.965, just = c("left", "top"),
            gp = gpar(fontsize = 9.5, col = "grey40"))
  if (!file.exists(path)) {
    grid.text(paste0("MISSING: ", it$file), gp = gpar(col = "red", fontsize = 13))
    footer(n); return(FALSE)
  }
  img <- readPNG(path)
  ar <- dim(img)[1] / dim(img)[2]           # height / width
  # Reserve room for the caption, then fit the image inside what is left.
  cap_lines <- length(wrap_lines(it$caption, 116))
  cap_h <- cap_lines * 0.0165 + 0.02
  # ASPECT RATIO: npc units are relative to the PAGE, and the page is 11 x 8.5,
  # so setting w = h / ar in npc distorts every image by 11/8.5 = 1.29. Size in
  # absolute inches instead, which is the only way the ratio survives.
  PAGE_W <- 11; PAGE_H <- 8.5
  avail_h <- 0.945 - cap_h - 0.075
  avail_h_in <- avail_h * PAGE_H
  avail_w_in <- 0.88 * PAGE_W
  h_in <- min(avail_h_in, avail_w_in * ar)
  w_in <- h_in / ar
  if (w_in > avail_w_in) { w_in <- avail_w_in; h_in <- w_in * ar }
  stopifnot(abs((h_in / w_in) - ar) < 1e-6)      # ratio preserved, in inches
  grid.raster(img, x = 0.5, y = 0.075 + cap_h + avail_h / 2 - 0.012,
              width = unit(w_in, "in"), height = unit(h_in, "in"))
  caption_block(it$caption, 0.075 + cap_h)
  footer(n)
  TRUE
}

render_tab <- function(it, n) {
  path <- file.path(TABDIR, it$file)
  grid.newpage(); header(it$id, it$title)
  if (!file.exists(path)) {
    grid.text(paste0("MISSING: ", it$file), gp = gpar(col = "red", fontsize = 13))
    footer(n); return(FALSE)
  }
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  # Long free-text columns make an unreadable table; wrap them.
  for (j in seq_along(df)) {
    if (is.character(df[[j]]))
      df[[j]] <- vapply(df[[j]], function(s)
        paste(wrap_lines(s, 46), collapse = "\n"), character(1))
  }
  cap_lines <- length(wrap_lines(it$caption, 116))
  cap_h <- cap_lines * 0.0165 + 0.02
  fs <- if (nrow(df) > 12 || ncol(df) > 5) 6.4 else 7.6
  tg <- tableGrob(df, rows = NULL,
                  theme = ttheme_minimal(
                    base_size = fs,
                    core = list(fg_params = list(hjust = 0, x = 0.02, lineheight = 1.05)),

                    colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02))))
  # Top-align: a short table centred in the box leaves a large gap under the
  # heading and reads as a rendering fault.
  th <- min(convertHeight(sum(tg$heights), "npc", valueOnly = TRUE),
            0.905 - cap_h - 0.075)
  pushViewport(viewport(x = 0.5, y = 0.905 - th / 2,
                        width = 0.9, height = th))
  grid.draw(tg); popViewport()
  caption_block(it$caption, 0.075 + cap_h)
  footer(n)
  TRUE
}

# ── Build ────────────────────────────────────────────────────────────────────
pdf(OUT, width = 11, height = 8.5, onefile = TRUE)

# Contents page is page 1, so every item's page number is its index + 1.
grid.newpage()
grid.text("TARA manuscript revision - revised figures and tables",
          x = 0.06, y = 0.94, just = c("left", "top"),
          gp = gpar(fontsize = 17, fontface = "bold"))
grid.text(paste(wrap_lines(paste0(
  "One figure or table per page, cited by page number in the point-by-point response. ",
  "Figure placement in the manuscript is NOT yet decided - these are the revised panels and ",
  "the new supplementary items we propose, not an assumed final figure set. Every number in a ",
  "caption is read from the source CSV named in the corresponding response entry."), 112),
  collapse = "\n"), x = 0.06, y = 0.885, just = c("left", "top"),
  gp = gpar(fontsize = 10, col = "grey30", lineheight = 1.35))

# Start the list BELOW the intro rather than at a fixed y: the intro wraps to a
# variable number of lines and a fixed start overlaps the first item.
.intro_lines <- length(wrap_lines(paste0(
  "One figure or table per page, cited by page number in the point-by-point response. ",
  "Figure placement in the manuscript is NOT yet decided - these are the revised panels and ",
  "the new supplementary items we propose, not an assumed final figure set. Every number in a ",
  "caption is read from the source CSV named in the corresponding response entry."), 112))
y <- 0.885 - .intro_lines * 0.019 - 0.035
for (k in seq_along(ITEMS)) {
  it <- ITEMS[[k]]
  grid.text(sprintf("%2d.", k + 1), x = 0.07, y = y, just = c("left", "top"),
            gp = gpar(fontsize = 10, col = "grey45"))
  grid.text(it$id, x = 0.115, y = y, just = c("left", "top"),
            gp = gpar(fontsize = 10, fontface = "bold"))
  grid.text(paste(wrap_lines(it$title, 82), collapse = "\n"),
            x = 0.215, y = y, just = c("left", "top"), gp = gpar(fontsize = 10))
  y <- y - max(0.030, 0.021 * length(wrap_lines(it$title, 82)))
}
footer(1)

for (k in seq_along(ITEMS)) {
  it <- ITEMS[[k]]
  n <- k + 1
  ok <- if (it$kind == "fig") render_fig(it, n) else render_tab(it, n)
  manifest <- rbind(manifest, data.frame(
    page = n, comment_id = it$id, kind = it$kind, file = it$file,
    title = it$title, present = ok, stringsAsFactors = FALSE))
}
invisible(dev.off())

save_tab("Revised_figure_pages", manifest)
cat("\nWrote:", OUT, "\n")
cat("Pages:", nrow(manifest) + 1, "(1 contents +", nrow(manifest), "items)\n")
miss <- manifest[!manifest$present, ]
if (nrow(miss)) {
  cat("\nMISSING SOURCE FILES - these pages are placeholders:\n")
  print(miss[, c("page", "comment_id", "file")], row.names = FALSE)
} else cat("every source figure and table was found\n")
print(manifest[, c("page", "comment_id", "kind", "file")], row.names = FALSE)
