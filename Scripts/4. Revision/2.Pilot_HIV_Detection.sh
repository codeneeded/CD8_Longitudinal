#!/usr/bin/env bash
################################################################################
# Revision step 2 — PILOT: is HIV detectable at all in one sample?
#
# WHY A PILOT
#   Realigning 26 libraries is 3-5 weeks of compute and disk. Before spending
#   that, establish on ONE sample whether vRNA+ cells are detectable with the
#   chemistry we actually used. Reviewer 1's request is only answerable if they
#   are; if not, the honest answer is a negative result with a power argument,
#   which costs days instead of weeks.
#
# WHICH SAMPLE, AND WHY
#   CP006_1m. Highest viral load in the cohort (10,000,000 copies/mL, at the
#   assay ceiling), pre-ART, 6,233 cells of which 2,750 are CD4 T. It also sits
#   in GSE254645, which is PUBLIC - no embargo to work around. If HIV is
#   detectable anywhere in this cohort it is detectable here, so a negative
#   here is close to decisive.
#
# THE THREE CONTROLS THIS RUNS, because a bare count is not evidence
#   1. EMPTY DROPLETS. Ambient viral RNA in the cell suspension gets into
#      droplets that contain no cell. If HIV counts in empty droplets are
#      comparable to cell-containing droplets, the "infected cells" are ambient
#      contamination. This uses the raw (unfiltered) matrix to measure it.
#   2. LINEAGE PLAUSIBILITY. HIV infects CD4 T cells. If the HIV+ barcodes are
#      spread evenly across B cells, monocytes and CD8s, that is contamination
#      or mismapping, not infection.
#   3. HOST TRANSCRIPTOME UNCHANGED. Adding viral contigs must not perturb the
#      human alignment. This compares per-cell gene and UMI counts against the
#      existing GRCh38-only alignment for the same sample.
#
# USAGE
#   bash "Scripts/4. Revision/2.Pilot_HIV_Detection.sh" --check
#   nohup bash "Scripts/4. Revision/2.Pilot_HIV_Detection.sh" \
#     > ~/Desktop/Projects/pilot_hiv.log 2>&1 &
################################################################################
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

PILOT_SAMPLE="${PILOT_SAMPLE:-CP006_1m}"
PILOT_GSE="${PILOT_GSE:-GSE254645}"

RAW_ROOT="${RAW_ROOT:-/mnt/data}"
FASTQ_DIR="$RAW_ROOT/TARA_GEO_raw/$PILOT_GSE/$PILOT_SAMPLE"

REF="${REF:-$HOME/Desktop/Projects/refdata/human_hiv/GRCh38_HIV1}"
OUT_PARENT="${OUT_PARENT:-$HOME/Desktop/Projects/cellranger_out}"
OUT_ID="pilot_hiv_${PILOT_SAMPLE}"

CR="${CR:-$HOME/Desktop/Projects/Tools/cellranger-10.1.0/cellranger}"
CORES="${CORES:-24}"
MEM="${MEM:-100}"

MODE="${1:-run}"
say() { printf '%s\n' "$*"; }
hr()  { printf '%s\n' "------------------------------------------------------------"; }

hr; say "PILOT: HIV DETECTION IN $PILOT_SAMPLE"; hr
say "  viral load:  10,000,000 copies/mL (assay ceiling)"
say "  ART status:  pre-ART, entry visit"
say "  cells:       6,233 total / 2,750 CD4 T in the published alignment"
say "  series:      $PILOT_GSE (public)"
say ""

################################################################################
# Preflight
################################################################################
say "cellranger:"
[[ -x "$CR" ]] && say "  $CR" || { say "  MISSING: $CR"; exit 1; }

say ""
say "combined reference:"
if [[ -d "$REF" ]]; then
  say "  $REF"
  for lab in HIV1_HXB2_subB HIV1_subC; do
    bp=$(grep "^$lab" "$REF/star/chrNameLength.txt" 2>/dev/null | awk '{print $2}')
    [[ -n "$bp" ]] && say "    $lab  ${bp} bp" || say "    $lab  NOT IN INDEX"
  done
else
  say "  MISSING: $REF"
  say "  Run: bash \"Scripts/4. Revision/1.Build_HIV_Reference.sh\""
  exit 1
fi

say ""
say "pilot FASTQs:"
if [[ -d "$FASTQ_DIR" ]]; then
  say "  $FASTQ_DIR"
  ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | sed 's#.*/#    #'
  n=$(ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
  say "    $n files, $(du -sh "$FASTQ_DIR" 2>/dev/null | cut -f1)"
  if [[ "$n" -lt 2 ]]; then
    say "  Need at least R1 and R2. Download this sample first:"
    say "    RAW_ROOT=$RAW_ROOT bash \"Scripts/4. Revision/0.Download_GEO_for_HIV_alignment.sh\""
    exit 1
  fi
  say ""
  say "  gzip integrity:"
  bad=0
  for f in "$FASTQ_DIR"/*.fastq.gz; do gzip -t "$f" 2>/dev/null || { say "    CORRUPT: $(basename "$f")"; bad=1; }; done
  [[ "$bad" -eq 0 ]] && say "    all pass" || exit 1
else
  say "  MISSING: $FASTQ_DIR"
  say "  Download just this sample:"
  say "    grep $PILOT_SAMPLE Reviewer_Comments/Tables/gex_download_list.txt"
  say "    then run step 0 - it skips completed samples, so you can start with one."
  exit 1
fi

say ""
say "compute:"
say "  cores requested: $CORES  (machine has $(nproc 2>/dev/null || echo '?'))"
say "  mem requested:   ${MEM} GB"
df -h "$OUT_PARENT" 2>/dev/null | tail -1 | awk '{print "  scratch free:    "$4}'

if [[ "$MODE" == "--check" ]]; then
  say ""; say "Preflight only. Re-run without --check to align."; exit 0
fi

################################################################################
# Align
################################################################################
hr; say "ALIGNING (expect 2-6 h for one library)"; hr

mkdir -p "$OUT_PARENT"; cd "$OUT_PARENT"
rm -rf "$OUT_ID"

"$CR" count \
  --id="$OUT_ID" \
  --transcriptome="$REF" \
  --fastqs="$FASTQ_DIR" \
  --sample="$PILOT_SAMPLE" \
  --create-bam=true \
  --localcores="$CORES" \
  --localmem="$MEM" \
  2>&1 | tail -25

OUTS="$OUT_PARENT/$OUT_ID/outs"
if [[ ! -d "$OUTS" ]]; then
  say ""; say "Alignment did not produce outs/. See the log above."; exit 1
fi

hr; say "ALIGNMENT COMPLETE"; hr
say "  outs: $OUTS"
ls "$OUTS" | sed 's/^/    /'

say ""
say "Next: score detection and run the three controls."
say "  Rscript \"Scripts/4. Revision/3.Score_HIV_Detection.R\""
say ""
say "  DO NOT read the raw HIV count as an infected-cell count until that"
say "  script has run. Ambient viral RNA in empty droplets is the dominant"
say "  false-positive source and it is what the controls measure."
