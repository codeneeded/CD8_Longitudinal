#!/usr/bin/env bash
################################################################################
# Revision — align the entry-visit libraries against GRCh38 + HIV-1
#
# Runs cellranger count with EXACTLY the pilot's invocation, so CP006_1m's
# result stays comparable to the nine new ones:
#   --transcriptome=refdata/human_hiv/GRCh38_HIV1  --create-bam=true
#   --localcores=24  --localmem=90
# (read from the pilot's own _cmdline rather than retyped, and asserted below.)
#
# ONE SAMPLE AT A TIME, deliberately. Each run is 4-8 h and holds the 24 GB
# STAR index resident; two at once would swap. The loop is resumable — a
# sample whose outs/ carries the success marker is skipped.
#
# USAGE
#   bash "Scripts/4. Revision/2b.Align_Entry_HIV.sh" --check
#   nohup bash "Scripts/4. Revision/2b.Align_Entry_HIV.sh" > ~/Desktop/Projects/entry_align.log 2>&1 &
#   bash "Scripts/4. Revision/2b.Align_Entry_HIV.sh" --only=CP022_1m
#
# Launch it from YOUR terminal with nohup: a multi-hour job started from the
# agent's shell dies when that cell ends (confirmed the hard way on CP018).
################################################################################
set -uo pipefail

P="${P:-$HOME/Desktop/Projects}"
CR="${CR:-$P/Tools/cellranger-10.1.0/cellranger}"
IDX="${IDX:-$P/refdata/human_hiv/GRCh38_HIV1}"
RAW_ROOT="${RAW_ROOT:-/mnt/data}"
DEST="$RAW_ROOT/TARA_GEO_raw"
OUT_PARENT="${OUT_PARENT:-$P/cellranger_out}"
LOG="${LOG:-$P/entry_align.log}"

# Cores/memory: match the pilot. nproc reports the CALLING shell's allowance,
# which is 8 in a sandbox and 32 on the host — so it must not set this.
CORES="${CORES:-24}"
MEM="${MEM:-90}"

MODE="run"; ONLY=""
for a in "$@"; do
  case "$a" in
    --check)  MODE="check" ;;
    --only=*) ONLY="${a#--only=}" ;;
    *) echo "unknown flag: $a"; exit 1 ;;
  esac
done

say() { echo "$*" | tee -a "$LOG"; }
hr()  { say "------------------------------------------------------------------"; }

# CP020_1m is last because its FASTQs are not on SRA (GSE328561 is embargoed)
# and are being placed by hand. The loop already skips any sample lacking a
# .complete marker, so listing it here is harmless until it appears - and
# means no edit is needed on the day it does.
SAMPLES="CP022_1m CP002_1m CP017_2m CP016_2m CP003_1m CP018_1m CP011_2m CP042_1m CP013_1m CP020_1m"

hr; say "ENTRY-VISIT HIV ALIGNMENT   $(date '+%F %T')"; hr

fail=0
[[ -x "$CR" ]] && say "  ok    cellranger $CR" || { say "  FAIL  cellranger not executable: $CR"; fail=1; }

# Reference: both viral contigs must be present, or every zero is meaningless.
if [[ -s "$IDX/genes/genes.gtf.gz" ]]; then
  nvir=$(zcat "$IDX/genes/genes.gtf.gz" | awk '$1 ~ /^HIV1_/{print $1}' | sort -u | wc -l)
  [[ "$nvir" -eq 2 ]] && say "  ok    reference, 2 viral contigs" \
                      || { say "  FAIL  reference has $nvir viral contigs, expected 2"; fail=1; }
else
  say "  FAIL  no GTF at $IDX/genes/genes.gtf.gz"; fail=1
fi

# Assert we are using the pilot's settings, by reading the pilot's own record.
PILOT_CMD="$OUT_PARENT/pilot_hiv_CP006_1m/_cmdline"
if [[ -f "$PILOT_CMD" ]]; then
  p_idx=$(grep -oE '\-\-transcriptome=[^ ]+' "$PILOT_CMD" | head -1 | cut -d= -f2)
  p_cores=$(grep -oE '\-\-localcores=[0-9]+' "$PILOT_CMD" | head -1 | cut -d= -f2)
  p_mem=$(grep -oE '\-\-localmem=[0-9]+' "$PILOT_CMD" | head -1 | cut -d= -f2)
  [[ "$p_idx" == "$IDX" ]] && say "  ok    reference matches the pilot's" \
    || say "  WARN  pilot used a DIFFERENT reference: $p_idx"
  [[ "$p_cores" == "$CORES" && "$p_mem" == "$MEM" ]] \
    && say "  ok    cores/mem match the pilot ($CORES/$MEM)" \
    || say "  WARN  pilot used --localcores=$p_cores --localmem=$p_mem; this run uses $CORES/$MEM"
else
  say "  WARN  no pilot _cmdline to compare against"
fi

avail=$(df -BG --output=avail "$OUT_PARENT" 2>/dev/null | tail -1 | tr -dc '0-9')
say "  ok    output dir $OUT_PARENT (${avail:-?} GB free; ~22 GB per sample)"
if [[ -n "${avail:-}" && "$avail" -lt 220 ]]; then
  say "  WARN  9 samples need ~200 GB of output; only ${avail} GB free"
fi

ready=0; missing=""
for s in $SAMPLES; do
  [[ -n "$ONLY" && "$s" != "$ONLY" ]] && continue
  d=$(find "$DEST" -maxdepth 2 -type d -name "$s" 2>/dev/null | head -1)
  if [[ -n "$d" && -f "$d/.complete" ]]; then ready=$((ready+1)); else missing="$missing $s"; fi
done
# Zero alignable samples must be a FAILURE, not a pass.
#
# This is the defect that wasted a night: with the data drive unmounted (or
# RAW_ROOT pointing elsewhere) every sample is skipped, the preflight still
# printed "passed", the loop then emitted nine [skip] lines and exited in under
# a second. From the outside a skipped-everything run looks like a long quiet
# one: no output directory, no error, normal RAM.
if [[ "$ready" -eq 0 ]]; then
  say "  FAIL  0 samples have complete FASTQs under $DEST"
  say "        the data drive is probably not mounted. Check:"
  say "          ls $DEST"
  say "          bash \"Scripts/4. Revision/0b.Download_Entry_HIV.sh\" --status"
  say "        or point RAW_ROOT at the right place:"
  say "          RAW_ROOT=/path/to/drive bash \"Scripts/4. Revision/2b.Align_Entry_HIV.sh\""
  fail=1
else
  say "  ok    FASTQs ready: $ready sample(s)"
fi
[[ -n "$missing" ]] && say "  WARN  not downloaded yet:$missing"

hr
[[ "$fail" -ne 0 ]] && { say "PREFLIGHT FAILED — nothing was aligned."; exit 1; }
say "preflight passed"
# Record the run's parameters in the log. say() already tees to $LOG, so the
# file exists from the banner onward - which is itself a useful diagnostic: if
# $LOG is absent, this script never ran, whatever the shell history suggests.
# What that tee does NOT capture is which inputs were used, so pin them here.
mkdir -p "$(dirname "$LOG")"
{ echo "=============================================================="
  echo "run started $(date '+%F %T')  ready=$ready  only=${ONLY:-all}"
  echo "  DEST=$DEST"
  echo "  IDX=$IDX"
  echo "  cores=$CORES mem=$MEM"
} >>"$LOG"
[[ "$MODE" == "check" ]] && { say "(--check: stopping here)"; exit 0; }

################################################################################
mkdir -p "$OUT_PARENT"
cd "$OUT_PARENT" || exit 1

for s in $SAMPLES; do
  [[ -n "$ONLY" && "$s" != "$ONLY" ]] && continue
  id="hiv_${s}"

  if [[ -f "$OUT_PARENT/$id/outs/metrics_summary.csv" ]]; then
    say "  [skip] $s already aligned"
    continue
  fi
  d=$(find "$DEST" -maxdepth 2 -type d -name "$s" 2>/dev/null | head -1)
  if [[ -z "$d" || ! -f "$d/.complete" ]]; then
    say "  [skip] $s FASTQs not complete — run 0b.Download_Entry_HIV.sh first"
    continue
  fi

  # A stale _lock from a killed run blocks cellranger with a confusing error.
  [[ -f "$OUT_PARENT/$id/_lock" ]] && { say "  clearing stale lock on $id"; rm -f "$OUT_PARENT/$id/_lock"; }

  hr; say "$s   started $(date '+%F %T')"; hr
  "$CR" count \
    --id="$id" \
    --transcriptome="$IDX" \
    --fastqs="$d" \
    --sample="$s" \
    --create-bam=true \
    --localcores="$CORES" \
    --localmem="$MEM" >>"$LOG" 2>&1

  if [[ -f "$OUT_PARENT/$id/outs/metrics_summary.csv" ]]; then
    # Gate: the viral features must be in the matrix, or a zero means nothing.
    nvf=$(zcat "$OUT_PARENT/$id/outs/filtered_feature_bc_matrix/features.tsv.gz" 2>/dev/null \
          | grep -c "^HIV1_" || true)
    say "  DONE $s   $(date '+%H:%M:%S')   viral features in matrix: ${nvf:-0}"
    [[ "${nvf:-0}" -ne 2 ]] && say "  ** WARNING: expected 2 viral features, found ${nvf:-0} — do NOT trust counts from this sample"
  else
    say "  FAILED $s — see $LOG"
  fi
done

hr; say "SUMMARY   $(date '+%F %T')"; hr
for s in $SAMPLES; do
  [[ -n "$ONLY" && "$s" != "$ONLY" ]] && continue
  m="$OUT_PARENT/hiv_${s}/outs/metrics_summary.csv"
  if [[ -f "$m" ]]; then
    cells=$(awk -F'","' 'NR==2{gsub(/[",]/,"",$1); print $1}' "$m" 2>/dev/null)
    printf '  %-11s aligned   %s cells\n' "$s" "${cells:-?}"
  else
    printf '  %-11s not aligned\n' "$s"
  fi
done
say ""
say "Next: score detection with the three controls, across all samples:"
say "  Rscript \"Scripts/4. Revision/3.Score_HIV_Detection.R\""
