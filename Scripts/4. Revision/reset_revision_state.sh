#!/usr/bin/env bash
################################################################################
# Revision — reset state so a step can be re-run cleanly
#
# Each target is opt-in. Nothing is deleted unless you name it, and every
# deletion is printed with its size first so you can see what you are about to
# lose. Downloaded FASTQs are NEVER touched unless you pass --fastq explicitly,
# because re-downloading them is hours.
#
# USAGE
#   bash "Scripts/4. Revision/reset_revision_state.sh" --status     # show only
#   bash "Scripts/4. Revision/reset_revision_state.sh" --ref        # reference build
#   bash "Scripts/4. Revision/reset_revision_state.sh" --pilot      # pilot alignment
#   bash "Scripts/4. Revision/reset_revision_state.sh" --sra-cache  # prefetch cache
#   bash "Scripts/4. Revision/reset_revision_state.sh" --outputs    # figures/tables
#   bash "Scripts/4. Revision/reset_revision_state.sh" --fastq      # RAW DATA (asks)
#
# Combine flags freely:  --ref --pilot
################################################################################
set -uo pipefail

RAW_ROOT="${RAW_ROOT:-/mnt/data}"
DEST="$RAW_ROOT/TARA_GEO_raw"
REF_OUT="${REF_OUT:-$HOME/Desktop/Projects/refdata/human_hiv}"
OUT_PARENT="${OUT_PARENT:-$HOME/Desktop/Projects/cellranger_out}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PILOT_SAMPLE="${PILOT_SAMPLE:-CP006_1m}"

say() { printf '%s\n' "$*"; }
hr()  { printf '%s\n' "------------------------------------------------------------"; }
sz()  { du -sh "$1" 2>/dev/null | cut -f1; }

DO_REF=0; DO_PILOT=0; DO_CACHE=0; DO_OUT=0; DO_FASTQ=0; STATUS_ONLY=0
[[ $# -eq 0 ]] && STATUS_ONLY=1
for a in "$@"; do
  case "$a" in
    --status)    STATUS_ONLY=1 ;;
    --ref)       DO_REF=1 ;;
    --pilot)     DO_PILOT=1 ;;
    --sra-cache) DO_CACHE=1 ;;
    --outputs)   DO_OUT=1 ;;
    --fastq)     DO_FASTQ=1 ;;
    --all)       DO_REF=1; DO_PILOT=1; DO_CACHE=1; DO_OUT=1 ;;
    *) say "unknown flag: $a"; exit 1 ;;
  esac
done

################################################################################
# Running processes — stop these BEFORE deleting anything underneath them
################################################################################
hr; say "RUNNING PROCESSES"; hr
found=0
for pat in "cellranger" "mkref" "prefetch" "fasterq-dump" "0.Download_GEO" "2.Pilot_HIV"; do
  pids=$(pgrep -f "$pat" 2>/dev/null | tr '\n' ' ')
  if [[ -n "${pids// /}" ]]; then
    say "  $pat  pids: $pids"; found=1
  fi
done
if [[ "$found" -eq 0 ]]; then
  say "  none"
else
  say ""
  say "  Stop them first, or deleting their working directory leaves the"
  say "  process writing into nowhere:"
  say "    pkill -f 'cellranger|mkref|prefetch|fasterq-dump'"
  say ""
  say "  A martian pipeline (mkref/count) also leaves a _lock file; removing"
  say "  the whole output directory as below clears it."
fi

################################################################################
# Current state
################################################################################
hr; say "CURRENT STATE"; hr

say "raw FASTQs        $DEST"
if [[ -d "$DEST" ]]; then
  nsamp=$(find "$DEST" -name ".complete" 2>/dev/null | wc -l)
  nfq=$(find "$DEST" -name "*.fastq.gz" 2>/dev/null | wc -l)
  say "                  $(sz "$DEST"), $nsamp samples complete, $nfq fastq.gz"
else
  say "                  (nothing downloaded yet)"
fi

say ""
say "SRA cache         $DEST/_sra_cache"
[[ -d "$DEST/_sra_cache" ]] && say "                  $(sz "$DEST/_sra_cache")" \
                            || say "                  (empty)"

say ""
say "HIV reference     $REF_OUT"
if [[ -d "$REF_OUT" ]]; then
  say "                  $(sz "$REF_OUT") total"
  [[ -d "$REF_OUT/_build" ]] && say "                    _build (intermediates): $(sz "$REF_OUT/_build")"
  if [[ -d "$REF_OUT/GRCh38_HIV1/star" ]]; then
    say "                    GRCh38_HIV1: $(sz "$REF_OUT/GRCh38_HIV1")  COMPLETE"
  elif [[ -d "$REF_OUT/GRCh38_HIV1" ]] || [[ -d "$REF_OUT/mkref_GRCh38_HIV1" ]]; then
    say "                    GRCh38_HIV1: INCOMPLETE - mkref did not finish"
  fi
else
  say "                  (not built)"
fi

say ""
say "pilot alignment   $OUT_PARENT/pilot_hiv_$PILOT_SAMPLE"
if [[ -d "$OUT_PARENT/pilot_hiv_$PILOT_SAMPLE" ]]; then
  say "                  $(sz "$OUT_PARENT/pilot_hiv_$PILOT_SAMPLE")"
  [[ -d "$OUT_PARENT/pilot_hiv_$PILOT_SAMPLE/outs" ]] && say "                  outs/ present - alignment finished" \
                                                      || say "                  no outs/ - incomplete"
  [[ -f "$OUT_PARENT/pilot_hiv_$PILOT_SAMPLE/_lock" ]] && say "                  _lock present - martian thinks it is still running"
else
  say "                  (not run)"
fi

say ""
say "revision outputs  $REPO/Reviewer_Comments"
say "                  $(find "$REPO/Reviewer_Comments" -name '*.png' 2>/dev/null | wc -l) figures, $(find "$REPO/Reviewer_Comments" -name '*.csv' 2>/dev/null | wc -l) tables"

if [[ "$STATUS_ONLY" -eq 1 ]]; then
  hr; say "Status only - nothing deleted."; say ""
  say "To reset a step, pass a flag:"
  say "  --ref        reference build (re-run 1.Build_HIV_Reference.sh, 30-60 min)"
  say "  --pilot      pilot alignment (re-run 2.Pilot_HIV_Detection.sh, 2-6 h)"
  say "  --sra-cache  prefetch's .sra staging area (safe - FASTQs already converted)"
  say "  --outputs    figures and tables (cheap to regenerate)"
  say "  --fastq      THE DOWNLOAD ITSELF - asks for confirmation"
  exit 0
fi

################################################################################
# Delete what was asked for
################################################################################
hr; say "RESETTING"; hr

if [[ "$DO_REF" -eq 1 ]]; then
  if [[ -d "$REF_OUT" ]]; then
    say "  removing reference build ($(sz "$REF_OUT"))"
    rm -rf "$REF_OUT"
    say "    done - re-run: bash \"Scripts/4. Revision/1.Build_HIV_Reference.sh\""
  else
    say "  reference: nothing to remove"
  fi
fi

if [[ "$DO_PILOT" -eq 1 ]]; then
  d="$OUT_PARENT/pilot_hiv_$PILOT_SAMPLE"
  if [[ -d "$d" ]]; then
    say "  removing pilot alignment ($(sz "$d"))"
    rm -rf "$d"
    say "    done - re-run: bash \"Scripts/4. Revision/2.Pilot_HIV_Detection.sh\""
  else
    say "  pilot: nothing to remove"
  fi
fi

if [[ "$DO_CACHE" -eq 1 ]]; then
  if [[ -d "$DEST/_sra_cache" ]]; then
    say "  clearing SRA cache ($(sz "$DEST/_sra_cache"))"
    # Only the staged .sra files, not the directory itself - vdb-config points here.
    find "$DEST/_sra_cache" -mindepth 1 -delete 2>/dev/null
    say "    done (directory kept so vdb-config stays valid)"
  else
    say "  SRA cache: nothing to remove"
  fi
fi

if [[ "$DO_OUT" -eq 1 ]]; then
  n=$(find "$REPO/Reviewer_Comments/Figures" "$REPO/Reviewer_Comments/Tables" \
        -type f \( -name '*.png' -o -name '*.csv' \) 2>/dev/null | wc -l)
  say "  removing $n revision figures/tables"
  # Keep the manifests the download depends on.
  find "$REPO/Reviewer_Comments/Figures" -type f -name '*.png' -delete 2>/dev/null
  find "$REPO/Reviewer_Comments/Figures" -type f -name '*_source_data.csv' -delete 2>/dev/null
  find "$REPO/Reviewer_Comments/Tables" -type f -name 'Pilot_HIV_*.csv' -delete 2>/dev/null
  find "$REPO/Reviewer_Comments/Tables" -type f -name 'TableS1_*.csv' -delete 2>/dev/null
  say "    done (kept geo_sample_map.csv and gex_download_list.txt - the download needs them)"
fi

if [[ "$DO_FASTQ" -eq 1 ]]; then
  say ""
  say "  ** THIS DELETES DOWNLOADED RAW DATA **"
  say "     $DEST  ($(sz "$DEST"))"
  say "     Re-downloading is hours. The SRA cache is separate (--sra-cache)."
  printf '     Type DELETE to confirm: '
  read -r ans
  if [[ "$ans" == "DELETE" ]]; then
    find "$DEST" -mindepth 1 -maxdepth 1 -name 'GSE*' -exec rm -rf {} + 2>/dev/null
    rm -f "$DEST/_runs.tsv" "$DEST/_failures.log" 2>/dev/null
    say "     removed (cache and _logs kept)"
  else
    say "     cancelled - nothing deleted"
  fi
fi

hr; say "DONE"; hr
say "Correct order from here:"
say "  1. download at least the pilot sample (CP006_1m is 6th in the manifest)"
say "       nohup bash \"Scripts/4. Revision/0.Download_GEO_for_HIV_alignment.sh\" \\"
say "         > ~/Desktop/Projects/geo_download.log 2>&1 &"
say "  2. build the reference (independent - can run at the same time)"
say "       bash \"Scripts/4. Revision/1.Build_HIV_Reference.sh\""
say "  3. pilot, once BOTH the above are done"
say "       bash \"Scripts/4. Revision/2.Pilot_HIV_Detection.sh\" --check"
say "  4. score it"
say "       Rscript \"Scripts/4. Revision/3.Score_HIV_Detection.R\""
