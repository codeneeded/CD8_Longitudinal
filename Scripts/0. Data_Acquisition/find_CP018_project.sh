#!/usr/bin/env bash
#===============================================================================
# Find which BaseSpace project/datasets hold the CP018 FACT-Seq libraries.
#
# Searches biosamples + datasets across all projects for CP018 / OB1-OB4 /
# "2m" / "90m" so we download the right project rather than guessing from
# the project name (which encodes only submitter + date + ticket number).
#
# Usage:  bash "Scripts/0. Data_Acquisition/find_CP018_project.sh"
#         bash "Scripts/0. Data_Acquisition/find_CP018_project.sh 512350847"   # single project
#===============================================================================

set -uo pipefail
export PATH="$HOME/bin:$PATH"
BS_ARGS=()
[[ -n "${BS_API:-}" ]] && BS_ARGS=(--api-server "$BS_API")

PATTERN='CP018|OB[1-4]|90m'

#-------------------------------------------------------------------------------
# Mode 1 — one project given: dump everything in it
#-------------------------------------------------------------------------------
if [[ $# -ge 1 ]]; then
  PID="$1"
  echo "=== Biosamples in project $PID ==="
  bs list biosamples --project-id "$PID" "${BS_ARGS[@]}" 2>/dev/null
  echo
  echo "=== Datasets in project $PID ==="
  bs list datasets --project-id "$PID" "${BS_ARGS[@]}" 2>/dev/null
  echo
  echo "=== AppSessions (shows what was run) ==="
  bs list appsessions --project-id "$PID" "${BS_ARGS[@]}" 2>/dev/null | head -20
  exit 0
fi

#-------------------------------------------------------------------------------
# Mode 2 — scan every project for CP018-like biosamples
#-------------------------------------------------------------------------------
echo "Scanning all projects for: $PATTERN"
echo "(newest projects are most likely — scanning newest first)"
echo

# newest-first by numeric project id
mapfile -t PIDS < <(bs list projects "${BS_ARGS[@]}" --terse 2>/dev/null | sort -rn)

if [[ ${#PIDS[@]} -eq 0 ]]; then
  echo "No projects returned. Is 'bs whoami' still working?"
  exit 1
fi

HITS=()
for pid in "${PIDS[@]}"; do
  [[ -z "$pid" ]] && continue
  name=$(bs get project --id "$pid" "${BS_ARGS[@]}" --terse 2>/dev/null | head -1)

  bios=$(bs list biosamples --project-id "$pid" "${BS_ARGS[@]}" 2>/dev/null | grep -iE "$PATTERN")
  dsets=$(bs list datasets  --project-id "$pid" "${BS_ARGS[@]}" 2>/dev/null | grep -iE "$PATTERN")

  if [[ -n "$bios" || -n "$dsets" ]]; then
    echo "############################################################"
    echo "HIT  project id=$pid  ${name:+name=$name}"
    echo "############################################################"
    [[ -n "$bios"  ]] && { echo "--- biosamples ---"; echo "$bios"; }
    [[ -n "$dsets" ]] && { echo "--- datasets ---";   echo "$dsets"; }
    echo
    HITS+=("$pid")
  else
    printf '  no match  %s\n' "$pid"
  fi
done

echo
if [[ ${#HITS[@]} -eq 0 ]]; then
  cat <<'MSG'
No CP018 match in any project.

Most likely: the BaseSpace share was viewed but not ACCEPTED into your account.
A share link is browsable without importing; the CLI only sees imported data.
Re-open https://basespace.illumina.com/s/fXOEy9uevGcc and look for an explicit
"Accept" / "Add to my projects" action, then re-run this script.
MSG
else
  echo "Matched project id(s): ${HITS[*]}"
  echo
  echo "Download with:"
  echo "  PROJECT_ID=${HITS[0]} bash \"Scripts/0. Data_Acquisition/download_CP018_FACTseq.sh\""
fi
