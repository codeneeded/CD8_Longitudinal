#!/bin/bash
# ─────────────────────────────────────────────────────────────
# Stage 1: BCR Pipeline - IgBLAST + Change-O database build
# Runs inside the Immcantation Docker container
#
# Usage (from inside container):
#   bash /scripts/run_bcr_stage1.sh
#
# Inputs:
#   /data         - mounted 10x_Genomics directory
#   /scripts      - mounted directory containing this script + manifest
#   /output       - mounted output directory
# ─────────────────────────────────────────────────────────────

set -euo pipefail

DATA_DIR="/data"
OUTPUT_DIR="/output"
MANIFEST="/scripts/sample_manifest.csv"
GERMLINE_DB="/usr/local/share/germlines/imgt/human/vdj"
IGBLAST_DB="/usr/local/share/igblast"

# ── Sanity checks ───────────────────────────────────────────
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found at $MANIFEST"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ── Log file for tracking progress ──────────────────────────
LOG="$OUTPUT_DIR/stage1_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "════════════════════════════════════════════════════════"
echo "BCR PIPELINE — STAGE 1"
echo "Started: $(date)"
echo "Manifest: $MANIFEST"
echo "════════════════════════════════════════════════════════"

# ── Read sample names from manifest (skip header) ───────────
SAMPLES=($(tail -n +2 "$MANIFEST" | cut -d',' -f1))
echo "Found ${#SAMPLES[@]} samples in manifest"
echo ""

# ── Track success/failure ───────────────────────────────────
SUCCESS=()
FAILED=()

# ── Loop through samples ────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do

    echo "────────────────────────────────────────────────────────"
    echo "Processing: $SAMPLE"
    echo "────────────────────────────────────────────────────────"

    BCR_DIR="$DATA_DIR/$SAMPLE/per_sample_outs/BCR"
    FASTA="$BCR_DIR/filtered_contig.fasta"
    CONTIGS="$BCR_DIR/filtered_contig_annotations.csv"
    OUTDIR="$OUTPUT_DIR/$SAMPLE"

    # Check inputs
    if [[ ! -f "$FASTA" ]] || [[ ! -f "$CONTIGS" ]]; then
        echo "  ✗ MISSING INPUT FILES — skipping"
        FAILED+=("$SAMPLE (missing input)")
        continue
    fi

    mkdir -p "$OUTDIR"

    # ── STEP 1: IgBLAST V(D)J gene assignment ───────────────
    echo "  [1/3] IgBLAST..."
    if ! AssignGenes.py igblast \
        -s "$FASTA" \
        -b "$IGBLAST_DB" \
        --organism human \
        --loci ig \
        --format blast \
        --outdir "$OUTDIR" \
        --outname "$SAMPLE" \
        --nproc 4 2>&1 | sed 's/^/      /'; then
        echo "  ✗ IgBLAST failed"
        FAILED+=("$SAMPLE (IgBLAST)")
        continue
    fi

    # ── STEP 2: Build Change-O AIRR database ────────────────
    echo "  [2/3] MakeDb..."
    if ! MakeDb.py igblast \
        -i "$OUTDIR/${SAMPLE}_igblast.fmt7" \
        -s "$FASTA" \
        -r "$GERMLINE_DB" \
        --10x "$CONTIGS" \
        --extended \
        --outdir "$OUTDIR" \
        --outname "$SAMPLE" 2>&1 | sed 's/^/      /'; then
        echo "  ✗ MakeDb failed"
        FAILED+=("$SAMPLE (MakeDb)")
        continue
    fi

    # ── STEP 3a: Filter productive heavy chains ─────────────
    echo "  [3/3] ParseDb filtering..."
    if ! ParseDb.py select \
        -d "$OUTDIR/${SAMPLE}_db-pass.tsv" \
        -f productive locus \
        -u T IGH \
        --outname "${SAMPLE}_heavy" \
        --outdir "$OUTDIR" 2>&1 | sed 's/^/      /'; then
        echo "  ✗ ParseDb (heavy) failed"
        FAILED+=("$SAMPLE (ParseDb heavy)")
        continue
    fi

    # ── STEP 3b: Filter productive light chains ─────────────
    ParseDb.py select \
        -d "$OUTDIR/${SAMPLE}_db-pass.tsv" \
        -f productive locus \
        -u T IGK IGL \
        --outname "${SAMPLE}_light" \
        --outdir "$OUTDIR" 2>&1 | sed 's/^/      /' || true   # don't fail if no light chains

    # Count sequences for reporting
    N_HEAVY=$(($(wc -l < "$OUTDIR/${SAMPLE}_heavy_parse-select.tsv") - 1))
    N_LIGHT=0
    [[ -f "$OUTDIR/${SAMPLE}_light_parse-select.tsv" ]] && \
        N_LIGHT=$(($(wc -l < "$OUTDIR/${SAMPLE}_light_parse-select.tsv") - 1))

    echo "  ✓ Done — heavy: $N_HEAVY, light: $N_LIGHT"
    SUCCESS+=("$SAMPLE")

done

# ── Summary ─────────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════════════════════"
echo "STAGE 1 COMPLETE"
echo "════════════════════════════════════════════════════════"
echo "Successful: ${#SUCCESS[@]} / ${#SAMPLES[@]}"
for s in "${SUCCESS[@]}"; do echo "  ✓ $s"; done

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo "Failed: ${#FAILED[@]}"
    for s in "${FAILED[@]}"; do echo "  ✗ $s"; done
fi

echo ""
echo "Next step: Run threshold_finder.R in R, then run_bcr_stage2.sh"
echo "Log saved to: $LOG"
