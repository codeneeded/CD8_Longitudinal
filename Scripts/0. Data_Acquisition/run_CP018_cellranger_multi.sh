#!/usr/bin/env bash
#===============================================================================
# CP018 FACT-Seq — Cell Ranger multi (GEX + TCR, OCM demultiplexing)
#
# Mirrors the CP003 FACT-Seq run: one pooled GEX library + one pooled TCR
# library, four on-chip-multiplexed (OCM) samples deconvolved by barcode.
#
#   sample_id      ocm_barcode   tube   timepoint
#   CP018_2m_A     OB1           1      2 months, replicate A
#   CP018_2m_B     OB2           2      2 months, replicate B
#   CP018_90m_A    OB3           3      90 months, replicate A
#   CP018_90m_B    OB4           4      90 months, replicate B
#
# Chemistry (verified from FASTQ headers): 10x 5' v2, R1 28bp (16bc+12umi),
# R2 90bp, NovaSeq X flowcell 2532FJLT4 lane 8.
#
# Output: <OUT_PARENT>/CP018_multi_out/outs/per_sample_outs/<sample_id>/
#           count/sample_filtered_feature_bc_matrix.h5
#           vdj_t/filtered_contig_annotations.csv
#         — the exact layout the CP003 R scripts expect.
#
# Usage:
#   bash "Scripts/0. Data_Acquisition/run_CP018_cellranger_multi.sh" --check
#   bash "Scripts/0. Data_Acquisition/run_CP018_cellranger_multi.sh"
#===============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_DIR="$HOME/Desktop/Projects/Repositories/CD8_Longitudinal"
REF_DIR="$HOME/Desktop/Projects/refdata/human"
OUT_PARENT="$HOME/Desktop/Projects/cellranger_out"
CONFIG="$REPO_DIR/Scripts/0. Data_Acquisition/CP018_multi_config.csv"
FASTQ_DIR="$REPO_DIR/fastq_raw/CP018_FACTseq"

GEX_REF="$REF_DIR/gex/refdata-gex-GRCh38-2024-A"
VDJ_REF="$REF_DIR/vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"

# Resources. Machine has 8 usable cores / 30 GB RAM, but only ~24 GB is
# typically FREE (browser, RStudio, desktop). Martian treats localmem as a
# scheduling budget and will over-commit if told it has more than it does —
# which risks an OOM kill during STAR alignment (GRCh38 index alone ~16 GB).
# 22 leaves headroom; raise only if the machine is otherwise idle.
LOCAL_CORES="${LOCAL_CORES:-7}"
LOCAL_MEM="${LOCAL_MEM:-22}"

CHECK_ONLY=0
[[ "${1:-}" == "--check" ]] && CHECK_ONLY=1

echo "==============================================================="
echo " CP018 FACT-Seq — Cell Ranger multi"
echo "==============================================================="

# ── 1. cellranger binary ──────────────────────────────────────────────────────
# Cell Ranger is licence-gated: it cannot be fetched non-interactively.
# Download the v9.0.1 tarball from your 10x account, then extract it once:
#   https://www.10xgenomics.com/support/software/cell-ranger/downloads
#   tar xf cellranger-<version>.tar -C ~/Desktop/Projects/Tools/
if ! command -v cellranger >/dev/null 2>&1; then
  for cand in "$HOME/Desktop/Projects/Tools"/cellranger-*/cellranger \
              "$HOME/Desktop/Projects"/cellranger-*/cellranger \
              "$HOME"/cellranger-*/cellranger \
              /opt/cellranger-*/cellranger; do
    [[ -x "$cand" ]] && { export PATH="$(dirname "$cand"):$PATH"; break; }
  done
fi

if ! command -v cellranger >/dev/null 2>&1; then
  cat <<'MSG'
[FAIL] cellranger not found.

  Cell Ranger requires accepting the 10x licence, so it cannot be downloaded
  automatically. One-time setup:

    1. https://www.10xgenomics.com/support/software/cell-ranger/downloads
    2. Download Cell Ranger (10.1.0 installed; 9.0.1 was used for CP003)
    3. tar xf cellranger-<version>.tar -C ~/Desktop/Projects/Tools/

  Then re-run this script — it finds the binary automatically.
MSG
  exit 1
fi
echo "[ok] cellranger: $(command -v cellranger)  ($(cellranger --version 2>&1 | head -1))"

# ── 2. References ─────────────────────────────────────────────────────────────
for ref in "$GEX_REF" "$VDJ_REF"; do
  if [[ -f "$ref/reference.json" ]]; then
    echo "[ok] reference: $(basename "$ref")"
  else
    echo "[FAIL] missing/incomplete reference: $ref"
    echo "       expected \$ref/reference.json — re-run the reference download"
    exit 1
  fi
done

# ── 3. FASTQs ─────────────────────────────────────────────────────────────────
n_fq=$(find "$FASTQ_DIR" -name "*.fastq.gz" 2>/dev/null | wc -l)
echo "[ok] FASTQ files found: $n_fq"
if [[ "$n_fq" -lt 4 ]]; then
  echo "[FAIL] expected at least 4 (GEX R1/R2 + TCR R1/R2); download may still be running"
  exit 1
fi

# Full gzip test on ~73 GB takes ~10 min. Set SKIP_FASTQ_CHECK=1 to skip once
# the files have already passed (e.g. relaunching after a verified --check).
bad=0
if [[ "${SKIP_FASTQ_CHECK:-0}" == "1" ]]; then
  echo "     [skipped] SKIP_FASTQ_CHECK=1 — assuming FASTQs already verified"
else
  echo "     verifying gzip integrity (catches truncated/in-progress transfers)..."
  while IFS= read -r f; do
    gzip -t "$f" 2>/dev/null || { echo "     INCOMPLETE/CORRUPT: $(basename "$f")"; bad=1; }
  done < <(find "$FASTQ_DIR" -name "*.fastq.gz")
fi
if [[ "$bad" -eq 0 ]]; then
  echo "     all FASTQs pass"
else
  cat <<'MSG'
[FAIL] One or more FASTQs failed the gzip test.

  If the BaseSpace download is STILL RUNNING this is expected — a partially
  written .gz always fails gzip -t. Wait for it to finish and re-run --check.

  If the download has finished, those files are genuinely truncated: re-run
    PROJECT_ID=512350847 bash "Scripts/0. Data_Acquisition/download_CP018_FACTseq.sh"
  which skips complete files and re-fetches the bad ones.
MSG
  exit 1
fi

# ── 4. Config sanity ──────────────────────────────────────────────────────────
[[ -f "$CONFIG" ]] || { echo "[FAIL] config missing: $CONFIG"; exit 1; }
echo "[ok] config: $CONFIG"
# every fastqs path in the config must exist
while IFS=, read -r fid fdir rest; do
  [[ "$fid" == "fastq_id" || -z "${fdir:-}" ]] && continue
  [[ -d "$fdir" ]] || { echo "[FAIL] config fastqs dir not found: $fdir"; exit 1; }
done < <(sed -n '/^\[libraries\]/,/^\[samples\]/p' "$CONFIG")
echo "     all library paths resolve"

# ── 5. Disk ───────────────────────────────────────────────────────────────────
avail=$(df -BG --output=avail "$HOME/Desktop" | tail -1 | tr -dc '0-9')
echo "[ok] disk available: ${avail} GB"
[[ "$avail" -lt 300 ]] && echo "     WARNING: cellranger multi needs ~200-300GB scratch for a run this size"

if [[ "$CHECK_ONLY" -eq 1 ]]; then
  echo
  echo "[--check] Preflight complete; nothing run."
  exit 0
fi

# ── 6. Run ────────────────────────────────────────────────────────────────────
mkdir -p "$OUT_PARENT"
cd "$OUT_PARENT"

if [[ -d CP018_multi_out ]]; then
  echo
  echo "NOTE: CP018_multi_out exists — cellranger will resume from the last"
  echo "      completed stage. To start clean: rm -rf $OUT_PARENT/CP018_multi_out"
fi

echo
echo "[run] starting cellranger multi — ${LOCAL_CORES} cores, ${LOCAL_MEM} GB"
echo "[run] expect ~6-14 h for a ~30GB GEX + ~7GB TCR library on 8 cores"
echo "[run] log: $OUT_PARENT/CP018_multi_out/_log"
echo

cellranger multi \
  --id=CP018_multi_out \
  --csv="$CONFIG" \
  --localcores="$LOCAL_CORES" \
  --localmem="$LOCAL_MEM"

# ── 7. Verify outputs ─────────────────────────────────────────────────────────
echo
echo "=== per-sample outputs ==="
PSO="$OUT_PARENT/CP018_multi_out/outs/per_sample_outs"
for s in CP018_2m_A CP018_2m_B CP018_90m_A CP018_90m_B; do
  h5="$PSO/$s/sample_filtered_feature_bc_matrix.h5"
  tcr="$PSO/$s/vdj_t/filtered_contig_annotations.csv"
  printf '%-14s  GEX:%s  TCR:%s\n' "$s" \
    "$([[ -f "$h5"  ]] && echo OK || echo MISSING)" \
    "$([[ -f "$tcr" ]] && echo OK || echo MISSING)"
done

echo
echo "Done. Next: adapt 'Scripts/CP003 - Longitudinal HIV Stim/1. QC_Full.R',"
echo "setting base_dir to:"
echo "  $PSO"
