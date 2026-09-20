#!/usr/bin/env bash
################################################################################
# Revision step 1 — build a combined GRCh38 + HIV-1 reference for cellranger
#
# ANSWERS: Reviewer 1, major point 3 (detect HIV-infected / vRNA+ cells).
#
# TWO DESIGN DECISIONS THAT AFFECT SENSITIVITY — read before running
#
# 1. SUBTYPE. The TARA cohort is from Maputo, Mozambique, where subtype C
#    predominates. Most published HIV scRNA-seq aligns to HXB2, a subtype B
#    reference. Divergence between B and C is roughly 10-15% at the nucleotide
#    level in env and less in gag/pol, which costs alignment sensitivity exactly
#    where we can least afford it. This script therefore builds the reference
#    with BOTH a subtype C consensus and HXB2 as separate contigs, so reads can
#    map to whichever is closer and the pilot reports the split. If subtype C
#    dominates the hits, drop HXB2 for the full run.
#
# 2. TRANSCRIPT MODEL. HIV produces overlapping spliced transcripts from a
#    single provirus, so annotating individual ORFs makes reads multi-map and
#    cellranger discards them. Instead each viral genome is annotated as ONE
#    exon spanning the whole contig, giving a single count per cell. That is
#    what we want: presence/absence per cell, not per-gene HIV expression.
#
# A KNOWN SENSITIVITY CEILING (from the literature, not from this data)
#   HIV transcripts are frequently NOT polyadenylated, so poly-dT priming in
#   standard 10x chemistry misses much of the viral RNA. Published work that
#   needed reliable detection added HIV-specific capture primers to the RT step
#   ("HIV-seq") and reported roughly a doubling of detected HIV RNA+ cells over
#   standard chemistry. Our libraries were made with standard chemistry, so
#   whatever this pipeline finds is a FLOOR, not the true infected fraction.
#   That is a Methods sentence and a caveat in the response letter.
#
# USAGE
#   bash "Scripts/4. Revision/1.Build_HIV_Reference.sh"           # build
#   bash "Scripts/4. Revision/1.Build_HIV_Reference.sh" --check   # verify only
################################################################################
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Human reference already on disk from the CP018 alignment
GEX_REF="${GEX_REF:-$HOME/Desktop/Projects/refdata/human/gex/refdata-gex-GRCh38-2024-A}"

# Output: the combined reference. NVMe - cellranger reads the index constantly.
REF_OUT="${REF_OUT:-$HOME/Desktop/Projects/refdata/human_hiv}"
WORK="$REF_OUT/_build"

CR="${CR:-$HOME/Desktop/Projects/Tools/cellranger-10.1.0/cellranger}"
THREADS="${THREADS:-8}"
MEM="${MEM:-48}"

# HIV-1 reference genomes from NCBI Nucleotide.
#   K03455  = HXB2, subtype B, the field's canonical coordinate system
#   AY772699 = subtype C consensus (South Africa, 04ZASK), close to Mozambique
HIV_ACCESSIONS=(K03455 AY772699)
declare -A HIV_LABEL=( [K03455]="HIV1_HXB2_subB" [AY772699]="HIV1_subC" )

MODE="${1:-build}"

say() { printf '%s\n' "$*"; }
hr()  { printf '%s\n' "------------------------------------------------------------"; }

hr; say "PREFLIGHT"; hr

say "tools:"
for t in "$CR" efetch samtools; do
  if command -v "$t" >/dev/null 2>&1 || [[ -x "$t" ]]; then
    say "  $(basename "$t")  ok"
  else
    say "  $(basename "$t")  MISSING"
  fi
done
say ""
say "  cellranger: $CR"
[[ -x "$CR" ]] || { say "  NOT EXECUTABLE - set CR=/path/to/cellranger"; exit 1; }

say ""
say "human reference:"
if [[ -d "$GEX_REF" ]]; then
  say "  $GEX_REF"
  say "  fasta: $(ls -la "$GEX_REF/fasta/genome.fa" 2>/dev/null | awk '{printf "%.1f GB", $5/1073741824}')"
  say "  gtf:   $(ls -la "$GEX_REF/genes/genes.gtf"* 2>/dev/null | head -1 | awk '{printf "%.0f MB", $5/1048576}')"
else
  say "  MISSING: $GEX_REF"
  say "  Set GEX_REF to the refdata-gex-GRCh38 directory."
  exit 1
fi

say ""
say "output:"
say "  $REF_OUT"
df -h "$(dirname "$REF_OUT")" 2>/dev/null | tail -1 | awk '{print "  "$4" free"}'
say "  mkref needs roughly 3x the genome size in scratch and produces ~16 GB."

if [[ "$MODE" == "--check" ]]; then
  say ""; say "Preflight only."; exit 0
fi

mkdir -p "$WORK" || exit 1

################################################################################
# 1. Fetch the HIV genomes
################################################################################
hr; say "FETCHING HIV GENOMES"; hr

for acc in "${HIV_ACCESSIONS[@]}"; do
  lab="${HIV_LABEL[$acc]}"
  fa="$WORK/${lab}.fa"
  if [[ -s "$fa" ]]; then
    say "  $lab  already present ($(grep -v '^>' "$fa" | tr -d '\n' | wc -c) bp)"
    continue
  fi
  printf '  %s (%s) ... ' "$lab" "$acc"
  if ! efetch -db nuccore -id "$acc" -format fasta < /dev/null > "$fa" 2>/dev/null; then
    say "FAILED"; say "  Could not fetch $acc. Check network access to NCBI."; exit 1
  fi
  # Rewrite the header to the contig name we want, on one line.
  bp=$(grep -v '^>' "$fa" | tr -d '\n' | wc -c)
  if [[ "$bp" -lt 8000 ]]; then
    say "SUSPECT ($bp bp - an HIV-1 genome is ~9.7 kb)"; exit 1
  fi
  { printf '>%s\n' "$lab"; grep -v '^>' "$fa" | tr -d '\n' | fold -w 60; echo; } > "$fa.tmp"
  mv "$fa.tmp" "$fa"
  say "ok ($bp bp)"
done

################################################################################
# 2. Build the combined FASTA and GTF
################################################################################
hr; say "BUILDING COMBINED REFERENCE"; hr

COMBI_FA="$WORK/genome_plus_hiv.fa"
COMBI_GTF="$WORK/genes_plus_hiv.gtf"

say "  concatenating genomes..."
HUMAN_FA="$GEX_REF/fasta/genome.fa"
cat "$HUMAN_FA" > "$COMBI_FA"
for acc in "${HIV_ACCESSIONS[@]}"; do
  cat "$WORK/${HIV_LABEL[$acc]}.fa" >> "$COMBI_FA"
done
say "  fasta: $(ls -la "$COMBI_FA" | awk '{printf "%.1f GB", $5/1073741824}')"

say "  concatenating annotations..."
HUMAN_GTF="$GEX_REF/genes/genes.gtf"
[[ -f "$HUMAN_GTF" ]] || HUMAN_GTF="$GEX_REF/genes/genes.gtf.gz"
if [[ "$HUMAN_GTF" == *.gz ]]; then
  zcat "$HUMAN_GTF" > "$COMBI_GTF"
else
  cat "$HUMAN_GTF" > "$COMBI_GTF"
fi

# ONE exon per viral genome. Overlapping ORF annotations would make reads
# multi-map and cellranger would drop them; a single feature gives one
# unambiguous count per cell, which is the presence/absence call we need.
for acc in "${HIV_ACCESSIONS[@]}"; do
  lab="${HIV_LABEL[$acc]}"
  len=$(grep -v '^>' "$WORK/${lab}.fa" | tr -d '\n' | wc -c)
  printf '%s\tHIV\tgene\t1\t%s\t.\t+\t.\tgene_id "%s"; gene_name "%s"; gene_type "protein_coding";\n' \
    "$lab" "$len" "$lab" "$lab" >> "$COMBI_GTF"
  printf '%s\tHIV\ttranscript\t1\t%s\t.\t+\t.\tgene_id "%s"; transcript_id "%s_tx"; gene_name "%s"; gene_type "protein_coding"; transcript_type "protein_coding";\n' \
    "$lab" "$len" "$lab" "$lab" "$lab" >> "$COMBI_GTF"
  printf '%s\tHIV\texon\t1\t%s\t.\t+\t.\tgene_id "%s"; transcript_id "%s_tx"; gene_name "%s"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number "1";\n' \
    "$lab" "$len" "$lab" "$lab" "$lab" >> "$COMBI_GTF"
  say "  added $lab as a single $len bp exon"
done

################################################################################
# 3. cellranger mkref
################################################################################
hr; say "RUNNING mkref (30-60 min)"; hr

cd "$REF_OUT"
rm -rf "$REF_OUT/GRCh38_HIV1"
"$CR" mkref \
  --genome=GRCh38_HIV1 \
  --fasta="$COMBI_FA" \
  --genes="$COMBI_GTF" \
  --nthreads="$THREADS" \
  --memgb="$MEM" \
  2>&1 | tail -20

################################################################################
# 4. Verify the viral contigs actually made it into the index
################################################################################
hr; say "VERIFYING"; hr

IDX="$REF_OUT/GRCh38_HIV1"
if [[ ! -d "$IDX" ]]; then
  say "  mkref did not produce $IDX - check the log above"; exit 1
fi

say "  index size: $(du -sh "$IDX" | cut -f1)"
say ""
say "  viral contigs in the index:"
for acc in "${HIV_ACCESSIONS[@]}"; do
  lab="${HIV_LABEL[$acc]}"
  n=$(grep -c "^$lab" "$IDX/star/chrNameLength.txt" 2>/dev/null || echo 0)
  bp=$(grep "^$lab" "$IDX/star/chrNameLength.txt" 2>/dev/null | awk '{print $2}')
  if [[ "$n" -ge 1 ]]; then
    say "    $lab  present  ${bp} bp"
  else
    say "    $lab  MISSING FROM INDEX - mkref dropped it"
  fi
done

say ""
say "  viral features in the gene list:"
grep -iE "HIV1_" "$IDX/star/geneInfo.tab" 2>/dev/null | sed 's/^/    /' || \
  grep -icE "HIV1_" "$IDX/star/exonInfo.tab" 2>/dev/null | sed 's/^/    features: /' || \
  say "    (could not read geneInfo.tab; check manually)"

say ""
say "Reference ready:"
say "  $IDX"
say ""
say "Next: pilot on ONE sample before committing to all 26."
say "  bash \"Scripts/4. Revision/2.Pilot_HIV_Detection.sh\""
