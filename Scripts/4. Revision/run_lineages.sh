#!/usr/bin/env bash
################################################################################
# Re-run the per-lineage subset -> recluster -> annotate pipeline, SAFELY.
#
#   bash "Scripts/4. Revision/run_lineages.sh"              # all four
#   bash "Scripts/4. Revision/run_lineages.sh" CD4          # one
#   INTEGRATION=cca bash "Scripts/4. Revision/run_lineages.sh"   # deliberate CCA
#
# WHY THIS EXISTS
#   A previous run cleared four checkpoints and only THEN discovered the
#   interpreter could not load SeuratWrappers/batchelor, destroying ~1.5 GB of
#   computed objects and replacing them with nothing. This wrapper makes that
#   sequence impossible:
#
#     1. find an Rscript that can actually load every required package
#     2. only then touch checkpoints - and MOVE them aside, never delete
#     3. run; on success the moved copies are removable, on failure restore them
#
#   Derived CSVs and PNGs under Analysis/<LIN>_subset/ are never touched: they
#   are what the response document quotes, and they survive an object loss.
################################################################################
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
S="$REPO/Scripts/4. Revision"
SAVED="$REPO/saved_R_data"
INTEGRATION="${INTEGRATION:-mnn}"
LINS=("${@:-}")
[[ -z "${LINS[0]:-}" ]] && LINS=(Myeloid NK B CD4)

say() { printf '  %s\n' "$*"; }
hr()  { printf '%s\n' "------------------------------------------------------------"; }

################################################################################
# 1. Find a working interpreter — BEFORE anything is destroyed
################################################################################
hr; echo "  PREFLIGHT"; hr

NEED='c("Seurat","SeuratObject","qs2","clustree","dplyr","tidyr","ggplot2")'
[[ "$INTEGRATION" == "mnn" ]] && NEED="c($NEED, \"SeuratWrappers\", \"batchelor\")"

CANDIDATES=()
[[ -n "${RSCRIPT:-}" ]] && CANDIDATES+=("$RSCRIPT")
CANDIDATES+=("$(command -v Rscript 2>/dev/null)")
CANDIDATES+=(/usr/lib/R/bin/Rscript /usr/local/bin/Rscript)
for e in "$HOME"/.claude-science/conda/envs/*/bin/Rscript; do
  [[ -x "$e" ]] && CANDIDATES+=("$e")
done

RS=""
for c in "${CANDIDATES[@]}"; do
  [[ -z "$c" || ! -x "$c" ]] && continue
  out=$("$c" -e "n <- $NEED; a <- n[!vapply(n, requireNamespace, logical(1), quietly=TRUE)];
        cat(if (length(a)) paste('MISSING:', paste(a, collapse=',')) else 'ALL_OK')" 2>&1 | tail -1)
  if [[ "$out" == *ALL_OK* ]]; then
    RS="$c"; say "usable: $c"; break
  else
    say "unusable: $c  ($(echo "$out" | cut -c1-70))"
  fi
done

if [[ -z "$RS" ]]; then
  echo
  say "NO interpreter can load the required packages for INTEGRATION=$INTEGRATION."
  say "NOTHING has been deleted or modified."
  say ""
  say "Run this from a terminal where R works (not from a sandboxed tool), or set"
  say "RSCRIPT=/path/to/Rscript. For a look without FastMNN: INTEGRATION=cca -"
  say "but CCA is NOT comparable to the published FastMNN CD8 analysis."
  exit 1
fi

################################################################################
# 2. Move checkpoints aside — never delete
################################################################################
STASH="$SAVED/.superseded_$(date +%Y%m%d_%H%M%S)"
moved=0
for lin in "${LINS[@]}"; do
  low=$(echo "$lin" | tr '[:upper:]' '[:lower:]')
  for k in reclustered annotated; do
    f="$SAVED/rev_${low}_${k}.qs2"
    if [[ -f "$f" ]]; then
      mkdir -p "$STASH"; mv "$f" "$STASH/"; moved=$((moved+1))
    fi
  done
done
echo
if (( moved )); then
  say "moved $moved existing checkpoint(s) to $(basename "$STASH")/ (NOT deleted)"
else
  say "no existing checkpoints to move"
fi

################################################################################
# 3. Run
################################################################################
fail=0
for lin in "${LINS[@]}"; do
  low=$(echo "$lin" | tr '[:upper:]' '[:lower:]')
  mkdir -p "$REPO/Analysis/${lin}_subset"
  echo; hr; echo "  $lin  $(date '+%H:%M:%S')"; hr
  INTEGRATION="$INTEGRATION" timeout 14400 "$RS" "$S/6.Lineage_Subset_Recluster.R" "$lin" \
    > "$REPO/Analysis/${lin}_subset/_run.log" 2>&1
  if [[ ! -f "$SAVED/rev_${low}_reclustered.qs2" ]]; then
    say "CLUSTERING FAILED - see Analysis/${lin}_subset/_run.log"
    grep -aE "^Error|PREFLIGHT|Execution halted" "$REPO/Analysis/${lin}_subset/_run.log" \
      | head -4 | sed 's/^/      /'
    fail=$((fail+1)); continue
  fi
  grep -aE "Integration method|^Subset:|provisional resolution" \
    "$REPO/Analysis/${lin}_subset/_run.log" | sed 's/^/    /'

  INTEGRATION="$INTEGRATION" timeout 7200 "$RS" "$S/7.Lineage_Annotation.R" "$lin" \
    > "$REPO/Analysis/${lin}_subset/_annot.log" 2>&1
  if [[ -f "$SAVED/rev_${low}_annotated.qs2" ]]; then
    grep -aE "confident:" "$REPO/Analysis/${lin}_subset/_annot.log" | tail -1 | sed 's/^/    /'
  else
    say "ANNOTATION FAILED - see Analysis/${lin}_subset/_annot.log"
    grep -aE "^Error|Execution halted" "$REPO/Analysis/${lin}_subset/_annot.log" \
      | head -3 | sed 's/^/      /'
    fail=$((fail+1))
  fi
done

################################################################################
# 4. Restore on failure, report on success
################################################################################
echo; hr
if (( fail )); then
  say "$fail lineage(s) failed."
  if [[ -d "$STASH" ]]; then
    restored=0
    for f in "$STASH"/*.qs2; do
      [[ -f "$f" ]] || continue
      b=$(basename "$f")
      [[ -f "$SAVED/$b" ]] || { mv "$f" "$SAVED/"; restored=$((restored+1)); }
    done
    say "restored $restored superseded checkpoint(s) so nothing is lost"
    rmdir "$STASH" 2>/dev/null
  fi
  say "Derived CSVs and PNGs under Analysis/*_subset/ were never touched."
else
  say "all ${#LINS[@]} lineage(s) completed with INTEGRATION=$INTEGRATION"
  [[ -d "$STASH" ]] && say "superseded copies kept in $(basename "$STASH")/ - remove when satisfied"
fi
hr
exit $(( fail > 0 ))
