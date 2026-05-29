#!/bin/bash
# ─────────────────────────────────────────────────────────────
# Stage 2: BCR Pipeline - Clone definition + Germline reconstruction
#
# Run AFTER threshold_finder.R has written thresholds.tsv
#
# Inputs (mounted in Docker):
#   /scripts/sample_manifest.csv
#   /output/thresholds.tsv               <- written by threshold_finder.R
#   /output/<SAMPLE>/<SAMPLE>_heavy_parse-select.tsv   <- from Stage 1
#
# Outputs per sample:
#   /output/<SAMPLE>/<SAMPLE>_clone-pass.tsv
#   /output/<SAMPLE>/<SAMPLE>_germ-pass.tsv          <- final input for R
# ─────────────────────────────────────────────────────────────

set -euo pipefail

DATA_DIR="/data"
OUTPUT_DIR="/output"
MANIFEST="/scripts/sample_manifest.csv"
THRESHOLDS="$OUTPUT_DIR/thresholds.tsv"
GERMLINE_DB="/usr/local/share/germlines/imgt/human/vdj"

DEFAULT_THRESHOLD="0.16"

# ── Sanity checks ───────────────────────────────────────────
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found at $MANIFEST"
    exit 1
fi

if [[ ! -f "$THRESHOLDS" ]]; then
    echo "ERROR: thresholds.tsv not found at $THRESHOLDS"
    echo "Did you run threshold_finder.R yet?"
    exit 1
fi

# ── Logging ─────────────────────────────────────────────────
LOG="$OUTPUT_DIR/stage2_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "════════════════════════════════════════════════════════"
echo "BCR PIPELINE — STAGE 2"
echo "Started: $(date)"
echo "Thresholds: $THRESHOLDS"
echo "════════════════════════════════════════════════════════"
echo ""
echo "Threshold table:"
cat "$THRESHOLDS"
echo ""

# ── Read sample list ────────────────────────────────────────
SAMPLES=($(tail -n +2 "$MANIFEST" | cut -d',' -f1))
echo "Found ${#SAMPLES[@]} samples in manifest"
echo ""

SUCCESS=()
FAILED=()

# ── Loop through samples ────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do

    echo "────────────────────────────────────────────────────────"
    echo "Processing: $SAMPLE"
    echo "────────────────────────────────────────────────────────"

    OUTDIR="$OUTPUT_DIR/$SAMPLE"
    HEAVY="$OUTDIR/${SAMPLE}_heavy_parse-select.tsv"

    if [[ ! -f "$HEAVY" ]]; then
        echo "  ✗ Stage 1 output missing — skipping"
        FAILED+=("$SAMPLE (no heavy file)")
        continue
    fi

    # Read threshold from R-generated table (column 2, tab-separated)
    THRESHOLD=$(awk -v s="$SAMPLE" '$1==s {print $2}' "$THRESHOLDS")

    if [[ -z "$THRESHOLD" ]]; then
        echo "  ! No threshold found for $SAMPLE, using default $DEFAULT_THRESHOLD"
        THRESHOLD="$DEFAULT_THRESHOLD"
    fi

    echo "  Threshold: $THRESHOLD"

    # ── STEP 4: Define clones ────────────────────────────────
    echo "  [1/2] DefineClones..."
    if ! DefineClones.py \
        -d "$HEAVY" \
        --act set \
        --model ham \
        --norm len \
        --dist "$THRESHOLD" \
        --outname "$SAMPLE" \
        --outdir "$OUTDIR" 2>&1 | sed 's/^/      /'; then
        echo "  ✗ DefineClones failed"
        FAILED+=("$SAMPLE (DefineClones)")
        continue
    fi

    # ── STEP 5: Reconstruct germlines (needed for SHM analysis) ──
    echo "  [2/2] CreateGermlines..."
    if ! CreateGermlines.py \
        -d "$OUTDIR/${SAMPLE}_clone-pass.tsv" \
        -g dmask \
        --cloned \
        -r "$GERMLINE_DB" \
        --outname "$SAMPLE" \
        --outdir "$OUTDIR" 2>&1 | sed 's/^/      /'; then
        echo "  ✗ CreateGermlines failed"
        FAILED+=("$SAMPLE (CreateGermlines)")
        continue
    fi

    # Report counts
    N_CLONES=$(tail -n +2 "$OUTDIR/${SAMPLE}_germ-pass.tsv" | \
               awk -F'\t' '{print $NF}' | sort -u | wc -l)
    N_SEQ=$(($(wc -l < "$OUTDIR/${SAMPLE}_germ-pass.tsv") - 1))

    echo "  ✓ Done — $N_SEQ sequences in $N_CLONES clones"
    SUCCESS+=("$SAMPLE")

done

# ── Summary ─────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════"
echo "STAGE 2 COMPLETE"
echo "════════════════════════════════════════════════════════"
echo "Successful: ${#SUCCESS[@]} / ${#SAMPLES[@]}"
for s in "${SUCCESS[@]}"; do echo "  ✓ $s"; done

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo "Failed: ${#FAILED[@]}"
    for s in "${FAILED[@]}"; do echo "  ✗ $s"; done
fi

echo ""
echo "Final output per sample: /output/<SAMPLE>/<SAMPLE>_germ-pass.tsv"
echo "These are the AIRR-format files for the R analysis pipeline."
echo "Log: $LOG"
