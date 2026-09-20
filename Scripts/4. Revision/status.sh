#!/usr/bin/env bash
################################################################################
# Revision — one command that answers "is it done?" for every running job.
#
#   bash "Scripts/4. Revision/status.sh"
#
# Reports, for each job: RUNNING / DONE / NOT STARTED / FAILED, with the
# evidence it used (a completion marker, an output file, a live pid) so the
# answer is never inferred from silence.
################################################################################
set -uo pipefail

RAW_ROOT="${RAW_ROOT:-/mnt/data}"
DEST="$RAW_ROOT/TARA_GEO_raw"
PROJ="$HOME/Desktop/Projects"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SAVED="$REPO/saved_R_data"
REF="$PROJ/refdata/human_hiv"

hr()  { printf '%s\n' "------------------------------------------------------------"; }
row() { printf '  %-26s %-13s %s\n' "$1" "$2" "$3"; }

alive() { pgrep -f "$1" >/dev/null 2>&1; }

hr; echo "  REVISION STATUS  $(date '+%Y-%m-%d %H:%M')"; hr

################################################################################
echo
echo "1. GEO PILOT DOWNLOAD  (CP006_1m, 2 lanes, ~55 GB)"
L="$PROJ/geo_download_pilot.log"
d="$DEST/GSE254645/CP006_1m"
if [[ -f "$d/.complete" ]]; then
  n=$(ls "$d"/CP006_1m_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  row "state" "DONE" "$n lanes, $(du -sh "$d" 2>/dev/null | cut -f1)"
elif alive "0.Download_GEO"; then
  have=$(ls "$d"/CP006_1m_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  cache=$(du -sh "$DEST/_sra_cache" 2>/dev/null | cut -f1)
  row "state" "RUNNING" "pid $(pgrep -f '0.Download_GEO' | head -1), $have/2 lanes converted"
  row "" "" "sra cache staging: ${cache:-0}"
  [[ -f "$L" ]] && row "" "" "last: $(grep -a 'L00\|complete' "$L" | tail -1 | sed 's/^ *//' | cut -c1-60)"
elif [[ -d "$d" ]]; then
  have=$(ls "$d"/CP006_1m_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  row "state" "INCOMPLETE" "$have/2 lanes, process not running - re-run to resume"
else
  row "state" "NOT STARTED" "$(if [[ -d "${RAW_ROOT:-/mnt/data}" ]]; then echo "no output directory"; else echo "${RAW_ROOT:-/mnt/data} not mounted here"; fi)"
fi

################################################################################
echo
echo "2. HIV REFERENCE BUILD  (GRCh38 + HIV-1 subtype B and C)"
if [[ -d "$REF/GRCh38_HIV1/star" ]]; then
  row "state" "DONE" "$(du -sh "$REF/GRCh38_HIV1" 2>/dev/null | cut -f1)"
  # The build is only usable if BOTH viral contigs survived mkref.
  gi="$REF/GRCh38_HIV1/star/chrName.txt"
  if [[ -f "$gi" ]]; then
    for c in HIV1_HXB2_subB HIV1_subC; do
      grep -qx "$c" "$gi" && row "" "" "$c present in index" \
                          || row "" "** MISSING" "$c NOT in index - rebuild before the pilot"
    done
  fi
elif alive "mkref"; then
  row "state" "RUNNING" "pid $(pgrep -f mkref | head -1), $(du -sh "$REF" 2>/dev/null | cut -f1) written"
  [[ -d "$REF/mkref_GRCh38_HIV1" ]] && \
    row "" "" "stage: $(ls -t "$REF/mkref_GRCh38_HIV1" 2>/dev/null | head -1)"
elif [[ -d "$REF" ]]; then
  row "state" "INCOMPLETE" "$(du -sh "$REF" 2>/dev/null | cut -f1), no index - reset_revision_state.sh --ref"
else
  row "state" "NOT STARTED" ""
fi

################################################################################
echo
echo "3. PILOT HIV ALIGNMENT  (cellranger count, GRCh38 + HIV-1)"
PD="$HOME/Desktop/Projects/cellranger_out/pilot_hiv_CP006_1m"
PL="$HOME/Desktop/Projects/pilot_hiv.log"
if [[ -f "$PD/outs/web_summary.html" ]]; then
  # finished: pull the viral counts straight out of the matrix features
  st="DONE"; det="$(du -sh "$PD" 2>/dev/null | cut -f1)"
  printf '  %-26s %-13s %s\n' "state" "$st" "$det"
  if [[ -f "$PD/outs/filtered_feature_bc_matrix/features.tsv.gz" ]]; then
    nv=$(zcat "$PD/outs/filtered_feature_bc_matrix/features.tsv.gz" 2>/dev/null \
         | grep -ci "HIV1_" || echo 0)
    printf '  %-26s %-13s %s\n' "viral features in matrix" "$nv" \
      "$([[ "$nv" -gt 0 ]] && echo "ok" || echo "** none - reference problem **")"
  fi
  printf '  %-26s %-13s %s\n' "next" "" "Rscript \"Scripts/4. Revision/3.Score_HIV_Detection.R\""
elif [[ -f "$PD/_lock" ]] \
     || pgrep -f "2.Pilot_HIV_Detection|pilot_hiv_CP006" >/dev/null 2>&1 \
     || { [[ -f "$PL" ]] && [[ $(( $(date +%s) - $(stat -c %Y "$PL" 2>/dev/null || echo 0) )) -lt 1800 ]]; }; then
  # running: name the stage, and the gzip pass so silence is not read as a hang
  stage=""
  [[ -f "$PL" ]] && stage=$(grep -aoE "gzip integrity|ALIGNING|ALIGNMENT COMPLETE" "$PL" | tail -1)
  cr=$(ls "$PD"/SC_RNA_COUNTER_CS 2>/dev/null | head -1)
  [[ -n "$cr" ]] && stage="cellranger running"
  [[ -z "$stage" ]] && stage="preflight"
  printf '  %-26s %-13s %s\n' "state" "RUNNING" "$stage"
  [[ -d "$PD" ]] && printf '  %-26s %-13s %s\n' "output so far" \
    "$(du -sh "$PD" 2>/dev/null | cut -f1)" "expect 6-10 h total"
  printf '  %-26s %-13s %s\n' "" "" "tail -f $PL"
elif [[ -f "$PL" ]] && grep -aq "ALIGNMENT FAILED\|Pipestance failed" "$PL" 2>/dev/null; then
  printf '  %-26s %-13s %s\n' "state" "FAILED" "see $PL"
else
  printf '  %-26s %-13s %s\n' "state" "NOT STARTED" \
    "bash \"Scripts/4. Revision/2.Pilot_HIV_Detection.sh\""
fi
echo

echo "4. LINEAGE RE-CLUSTERING  (answers Rev 1 major 2)"
for lin in cd4 nk b myeloid; do
  q="$SAVED/rev_${lin}_reclustered.qs2"
  up=$(echo "$lin" | tr '[:lower:]' '[:upper:]')
  if [[ -f "$q" ]]; then
    row "$up" "DONE" "$(du -h "$q" | cut -f1) checkpoint"
  elif alive "6.Lineage_Subset_Recluster.R $up"; then
    row "$up" "RUNNING" "pid $(pgrep -f "6.Lineage_Subset_Recluster.R $up" | head -1)"
  elif [[ -f "$REPO/Analysis/${up}_subset/_run.log" ]]; then
    if grep -qa "^Error\|Execution halted" "$REPO/Analysis/${up}_subset/_run.log"; then
      row "$up" "FAILED" "see Analysis/${up}_subset/_run.log"
    else
      row "$up" "INCOMPLETE" "log exists, no checkpoint"
    fi
  else
    row "$up" "NOT STARTED" ""
  fi
done

################################################################################
echo
echo "5. ANNOTATION  (script 7)"
for lin in cd4 nk b myeloid; do
  q="$SAVED/rev_${lin}_annotated.qs2"
  up=$(echo "$lin" | tr '[:lower:]' '[:upper:]')
  [[ -f "$q" ]] && row "$up" "DONE" "$(du -h "$q" | cut -f1)" \
                || row "$up" "pending" ""
done

################################################################################
echo
echo "6. OUTPUTS"
row "reviewer figures" "" "$(find "$REPO/Reviewer_Comments/Figures" -name '*.png' 2>/dev/null | wc -l) PNG"
row "reviewer tables" "" "$(find "$REPO/Reviewer_Comments/Tables" -name '*.csv' 2>/dev/null | wc -l) CSV"
row "qs2 checkpoints" "" "$(ls "$SAVED"/rev_*.qs2 2>/dev/null | wc -l) objects"
RD="$HOME/Desktop/Projects/Deliverables/CP018/collaborator_package/TARA_Response_to_Reviewers.pdf"
[[ -f "$RD" ]] && row "response PDF" "" "$(du -h "$RD" | cut -f1), $(date -r "$RD" '+%H:%M')" \
               || row "response PDF" "" "not built yet"

################################################################################
echo
echo "7. MACHINE"
free -g | awk 'NR==2{printf "  %-26s %-13s %s\n", "memory", $7" GB free", $3" GB in use"}'
HOST_CPUS=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || echo "?")
SHELL_CPUS=$(nproc 2>/dev/null || echo "?")
printf '  %-26s %-13s %s\n' "load" "$(cut -d' ' -f1-3 /proc/loadavg)" "on $HOST_CPUS cores"
if [[ "$SHELL_CPUS" != "$HOST_CPUS" ]]; then
  printf '  %-26s %-13s %s\n' "" "" \
    "(this shell is limited to $SHELL_CPUS of $HOST_CPUS - cellranger run from"
  printf '  %-26s %-13s %s\n' "" "" \
    " your own terminal gets all $HOST_CPUS; size --localcores against that)"
fi
df -h "$PROJ" 2>/dev/null | tail -1 | awk '{printf "  %-26s %-13s %s\n", "nvme free", $4, $5" used"}'
df -h "$RAW_ROOT" 2>/dev/null | tail -1 | awk '{printf "  %-26s %-13s %s\n", "data drive free", $4, $5" used"}'

echo
hr
echo "  R jobs write Analysis/<LIN>_subset/_run.log - tail that for detail."
hr
