#!/usr/bin/env bash
# CP018 FACT-Seq pipeline runner
#
# Runs the scripts IN ORDER and STOPS at the first failure. Each script's full
# output goes to logs/<script>.log so a crash leaves a readable trace instead of
# scrolling off the terminal.
#
# Order matters: script 3 rebuilds the integrated object and attaches
# Fig_Annotation. Scripts 4, 5, 7 refuse to run without it, so running them
# against a stale object fails immediately rather than producing wrong numbers.
#
#   bash "Scripts/CP018 - FACTseq/run_pipeline.sh"           # from script 3
#   bash "Scripts/CP018 - FACTseq/run_pipeline.sh" 4         # resume from 4
#   bash "Scripts/CP018 - FACTseq/run_pipeline.sh" --list    # show stages

set -uo pipefail

SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SDIR/../.." && pwd)"
LOGS="$REPO/Analysis/CP018/logs"
mkdir -p "$LOGS"

RSCRIPT="${RSCRIPT:-Rscript}"
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "Rscript not found. Set RSCRIPT=/path/to/Rscript"; exit 1; }

# Stage 3 re-runs CCA + FastMNN integration and needs ~25 GB peak. If the
# integrated object already exists, stage 3a does the cheap remainder instead
# (set resolution, annotate, plot) - no integration, ~2 GB.
STAGES=(
  "3a:3a.Set_Resolution_Annotate.R:set res 0.4 + annotate + annotated plots (NO re-integration)"
  "3:3.Integration_Anno_Plots.R:FULL re-integration (only if the object is missing)"
  "3b:3b.Annotation_Evidence.R:annotation evidence + CP003 correspondence"
  "4:4.VDJ_Seurat.R:attach TCR, repertoire plots"
  "4b:4b.Clonotype_Persistence.R:persistence alluvial + detection controls"
  "5:5.HIV_Specific_TCR.R:HIV-specific call (expansion-based)"
  "8:8.CP003_vs_CP018_Comparison.R:CP003 vs CP018 comparison"
)

if [[ "${1:-}" == "--list" ]]; then
  printf '%-5s %-34s %s\n' "STAGE" "SCRIPT" "DOES"
  for s in "${STAGES[@]}"; do IFS=: read -r id f d <<< "$s"; printf '%-5s %-34s %s\n' "$id" "$f" "$d"; done
  exit 0
fi

START_AT="${1:-3a}"
started=0
FAILED=""

echo "=============================================="
echo "CP018 pipeline | $(date '+%Y-%m-%d %H:%M:%S')"
echo "logs: $LOGS"
echo "=============================================="

for s in "${STAGES[@]}"; do
  IFS=: read -r id f desc <<< "$s"
  [[ "$started" -eq 0 && "$id" != "$START_AT" ]] && continue
  started=1

  # Never auto-run full re-integration when the object is already present -
  # it costs ~25 GB and 3a achieves the same end state from the saved object.
  if [[ "$id" == "3" && "$START_AT" != "3" ]]; then
    if [[ -f "$REPO/saved_R_data/CP018_RNA_integrated_CCA_MNN.qs2" ]]; then
      echo
      echo ">>> [3] SKIPPED - integrated object exists; 3a already did the work."
      echo "    To force full re-integration: bash \"$0\" 3"
      continue
    fi
  fi

  echo
  echo ">>> [$id] $f"
  echo "    $desc"
  t0=$(date +%s)
  log="$LOGS/${id}.log"

  if "$RSCRIPT" "$SDIR/$f" > "$log" 2>&1; then
    printf '    OK  (%ds)\n' "$(( $(date +%s) - t0 ))"
    # surface the lines worth seeing without dumping the whole log
    grep -E "^(Dropping tiny|Annotation applied|Reusing marker|Using:)" "$log" | sed 's/^/    /' || true
    grep -E "\*\*\*" "$log" | sed 's/^/    /' || true
  else
    printf '    FAILED (%ds)\n' "$(( $(date +%s) - t0 ))"
    echo "    ---- last 25 lines of $log ----"
    tail -25 "$log" | sed 's/^/    /'
    FAILED="$id"
    break
  fi
done

echo
if [[ -n "$FAILED" ]]; then
  echo "STOPPED at stage $FAILED. Full log: $LOGS/${FAILED}.log"
  echo "Fix, then resume with:  bash \"$0\" $FAILED"
  exit 1
fi
echo "All stages completed. Logs in $LOGS"
