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
# Stage 0 (clonotype normalisation) deliberately sits AFTER stage 4: it needs
# CTstrict to exist, which only happens once TCR data is attached. Everything
# cross-experiment (6, 6b, 8, 9) runs after it and refuses to start without it.
#
#   bash "Scripts/CP018 - FACTseq/run_pipeline.sh"           # everything, in order
#   bash "Scripts/CP018 - FACTseq/run_pipeline.sh" 4         # resume from stage 4
#   bash "Scripts/CP018 - FACTseq/run_pipeline.sh" --plots   # figures only, skip
#                                                            # the heavy rebuilds
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
  "0:0.Normalise_Clonotypes.R:stamp canonical clonotype column on every object"
  "4b:4b.Clonotype_Persistence.R:persistence alluvial + detection controls"
  "5:5.HIV_Specific_TCR.R:HIV-specific call (expansion-based)"
  "6:6.HIV_Specific_TCR_TARA_bulk.R:cross-reference into unstimulated eCITE-seq"
  "6b:6b.Lineage_Timeline.R:lineage timeline pre-ART -> 90m + resting phenotype"
  "8:8.CP003_vs_CP018_Comparison.R:CP003 vs CP018 comparison"
  "9:9.Integrate_CP003_CP018.R:joint integration + Tex-survival check"
  "10:10.Figure3_Mirror_CP018.R:manuscript Figure 3 mirror + CP003/CP018 combined panels"
  "11:11.Story_Figures.R:story figures - within-participant effect sizes, IFN inversion, lineage traces"
)

if [[ "${1:-}" == "--list" ]]; then
  printf '%-5s %-34s %s\n' "STAGE" "SCRIPT" "DOES"
  for s in "${STAGES[@]}"; do IFS=: read -r id f d <<< "$s"; printf '%-5s %-34s %s\n' "$id" "$f" "$d"; done
  exit 0
fi

# ── Preflight: check every R package the pipeline needs, BEFORE running 90s of
# stage 3b only to die at stage 4 on a missing dependency. Names the stages that
# would fail rather than just the package.
if [[ "${SKIP_PREFLIGHT:-0}" != "1" ]]; then
  echo "Checking R packages..."
  missing=$("$RSCRIPT" -e '
    need <- list(
      "Seurat"="all", "qs2"="all", "dplyr"="all", "ggplot2"="all", "tidyr"="all",
      "scRepertoire"="4", "ggalluvial"="4b, 6b, 10", "RColorBrewer"="10", "scCustomize"="3a, 3b, 4",
      "clustree"="3", "SeuratWrappers"="9 (FastMNN; falls back to CCA)", "patchwork"="several", "scales"="several",
      "Matrix"="all", "cowplot"="4")
    miss <- names(need)[!vapply(names(need), requireNamespace, logical(1), quietly=TRUE)]
    if (length(miss)) for (m in miss) cat(m, " (needed by stage: ", need[[m]], ")\n", sep="")
  ' 2>/dev/null)

  if [[ -n "$missing" ]]; then
    echo
    echo "MISSING R PACKAGES:"
    echo "$missing" | sed 's/^/  /'
    echo
    echo "  Install them, then re-run. For Bioconductor packages:"
    echo "    R -e 'BiocManager::install(c(\"scRepertoire\",\"GenomeInfoDbData\"))'"
    echo
    echo "  To run anyway and skip the stages that need them:"
    echo "    SKIP_PREFLIGHT=1 bash \"$0\" <stage>"
    exit 1
  fi
  echo "  all present"
fi


# --plots: regenerate every figure and table from the SAVED objects, skipping
# the two expensive rebuild stages (3 = re-integration, 9 = joint integration
# when its checkpoint already exists). Use this after a plotting change.
PLOTS_ONLY=0
if [[ "${1:-}" == "--plots" ]]; then
  PLOTS_ONLY=1
  shift || true
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

  # In --plots mode skip the stages that rebuild objects rather than draw them.
  if [[ "$PLOTS_ONLY" -eq 1 && "$id" == "3" ]]; then
    echo; echo ">>> [3] SKIPPED (--plots: re-integration is a rebuild, not a plot)"
    continue
  fi
  if [[ "$PLOTS_ONLY" -eq 1 && "$id" == "9" \
        && -f "$REPO/saved_R_data/CP003_CP018_integrated.qs2" ]]; then
    echo; echo ">>> [9] SKIPPED (--plots: joint object exists; re-run without --plots to rebuild)"
    continue
  fi

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
    grep -aE "^(Dropping tiny|Annotation applied|Reusing marker|Using:|Integration method|VERDICT|Legacy objects)" "$log" | sed 's/^/    /' || true
    grep -acE "^  \[csv\]" "$log" | sed 's/^/    tables written: /' || true
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
# Optional: verify the write-up numbers against the tables this run produced.
# The verifier and the write-ups live OUTSIDE this repo, under
# ~/Desktop/Projects/Deliverables/CP018/ — the repo holds the manuscript
# (code, figures, tables), not correspondence or reporting tooling. Skipped
# silently when that folder is absent, so the pipeline is self-contained.
CLAIMS="${DELIVERABLES:-$HOME/Desktop/Projects/Deliverables/CP018}/tooling/verify_claims.py"
if [[ -f "$CLAIMS" ]]; then
  echo
  echo ">>> verifying write-up claims against source tables"
  python3 "$CLAIMS" | sed 's/^/    /' || \
    echo "    *** claims above did not match - check before circulating the docs"
fi

echo "All stages completed. Logs in $LOGS"
echo
echo "Outputs:"
echo "  figures: $(find "$REPO/Analysis/CP018" -name '*.png' 2>/dev/null | wc -l)"
echo "  tables : $(find "$REPO/Analysis/CP018" -name '*.csv' 2>/dev/null | wc -l)"
