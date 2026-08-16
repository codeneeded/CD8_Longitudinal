#!/usr/bin/env bash
#===============================================================================
# Download CP018 FACT-Seq FASTQs from Illumina BaseSpace
#
# Run:  ~2m timepoint in duplicate (OB1, OB2), 90m timepoint in duplicate (OB3, OB4)
# Source: BaseSpace share https://basespace.illumina.com/s/fXOEy9uevGcc
#         (from L. de Armas, for CP018 FACT-Seq)
#
#   Sample Name    Tube Label    5'OCM Barcode    Chip lane
#   CP018 2m_A     1             OB1              1
#   CP018 2m_B     2             OB2              2
#   CP018 90m_A    3             OB3              3
#   CP018 90m_B    4             OB4              4
#
# Destination: <repo>/fastq_raw/CP018_FACTseq/   (gitignored)
#
# Usage:
#   bash "Scripts/0. Data_Acquisition/download_CP018_FACTseq.sh"            # download
#   bash "Scripts/0. Data_Acquisition/download_CP018_FACTseq.sh" --list     # inspect only
#===============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_DIR="$HOME/Desktop/Projects/Repositories/CD8_Longitudinal"
DEST="$REPO_DIR/fastq_raw/CP018_FACTseq"
BS_BIN="$HOME/bin/bs"

LIST_ONLY=0
[[ "${1:-}" == "--list" ]] && LIST_ONLY=1

mkdir -p "$DEST" "$(dirname "$BS_BIN")"

#-------------------------------------------------------------------------------
# STEP 0 — Accept the share in your browser FIRST (one time only)
#-------------------------------------------------------------------------------
# The link Lesley sent is a *share* link. Until you accept it, the run will not
# appear in your account and the CLI cannot see it.
#   1. Open https://basespace.illumina.com/s/fXOEy9uevGcc  (you are already logged in)
#   2. Accept / import the share so the run + project land in your account.

#-------------------------------------------------------------------------------
# STEP 1 — Install the BaseSpace CLI (one time only)
#-------------------------------------------------------------------------------
if [[ ! -x "$BS_BIN" ]]; then
  echo "[setup] Installing BaseSpace CLI to $BS_BIN"
  curl -L "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -o "$BS_BIN"
  chmod u+x "$BS_BIN"
fi
export PATH="$HOME/bin:$PATH"

#-------------------------------------------------------------------------------
# STEP 2 — Authenticate (one time only; token cached in ~/.basespace/)
#-------------------------------------------------------------------------------
# Prints a URL — open it in the browser where you are already logged in and approve.
#
# NOTE: test real auth with `bs whoami`, NOT for the presence of default.cfg —
# an interrupted `bs auth` can leave a config file behind with no valid token.
#
# Region: default is US. If your account lives in another BaseSpace region,
# set BS_API before running, e.g.
#   BS_API=https://api.euc1.sh.basespace.illumina.com   (EU / Frankfurt)
BS_API="${BS_API:-}"
BS_ARGS=()
[[ -n "$BS_API" ]] && BS_ARGS=(--api-server "$BS_API")

if ! "$BS_BIN" whoami "${BS_ARGS[@]}" >/dev/null 2>&1; then
  echo "[auth] Launching device authentication — approve in your browser"
  echo "       (if it hangs >60s after you approve: Ctrl+C, then re-run — the token"
  echo "        is often already cached and this step will be skipped)"
  "$BS_BIN" auth --force "${BS_ARGS[@]}"
fi
"$BS_BIN" whoami "${BS_ARGS[@]}"

#-------------------------------------------------------------------------------
# STEP 3 — Locate the run / project
#-------------------------------------------------------------------------------
echo
echo "=== Projects in your account (newest last) ==="
"$BS_BIN" list projects "${BS_ARGS[@]}"
echo
echo "=== Runs in your account ==="
"$BS_BIN" list runs "${BS_ARGS[@]}"
echo
echo "=== Datasets matching CP018 ==="
"$BS_BIN" list datasets "${BS_ARGS[@]}" | grep -i -E "CP018|OB[1-4]" || echo "  (no CP018 match — check the share was accepted in Step 0)"

if [[ "$LIST_ONLY" -eq 1 ]]; then
  echo
  echo "[--list] Inspection only; nothing downloaded."
  echo "Set PROJECT_NAME below to the project you see above, then re-run without --list."
  exit 0
fi

#-------------------------------------------------------------------------------
# STEP 4 — Download
#-------------------------------------------------------------------------------
# Prefer PROJECT_ID (unambiguous) over PROJECT_NAME. Get the id from
# find_CP018_project.sh, which searches biosample names for CP018/OB1-4/90m.
PROJECT_ID="${PROJECT_ID:-}"
PROJECT_NAME="${PROJECT_NAME:-}"

if [[ -z "$PROJECT_ID" && -z "$PROJECT_NAME" ]]; then
  cat <<'MSG'

  ACTION NEEDED — identify the project first:

      bash "Scripts/0. Data_Acquisition/find_CP018_project.sh"

  then re-run with the id it reports:

      PROJECT_ID=<id> bash "Scripts/0. Data_Acquisition/download_CP018_FACTseq.sh"

MSG
  exit 1
fi

if [[ -n "$PROJECT_ID" ]]; then
  SELECTOR=(--id "$PROJECT_ID"); LABEL="id $PROJECT_ID"
else
  SELECTOR=(--name "$PROJECT_NAME"); LABEL="name '$PROJECT_NAME'"
fi

echo
echo "[download] Project $LABEL  ->  $DEST"
echo "[download] Large transfer — safe to re-run; completed files are skipped."
"$BS_BIN" download project \
    "${SELECTOR[@]}" \
    --output "$DEST" \
    --extension=fastq.gz \
    --no-metadata "${BS_ARGS[@]}"

# ── Alternative: download by RUN id (gives raw BCL-converted FASTQs) ──────────
# RUN_ID=$("$BS_BIN" list runs --terse | head -1)
# "$BS_BIN" download run --id "$RUN_ID" --output "$DEST"

#-------------------------------------------------------------------------------
# STEP 5 — Verify
#-------------------------------------------------------------------------------
echo
echo "=== Downloaded FASTQs ==="
find "$DEST" -name "*.fastq.gz" -printf "%10s  %P\n" | sort -k2
echo
echo "file count : $(find "$DEST" -name '*.fastq.gz' | wc -l)"
echo "total size : $(du -sh "$DEST" | cut -f1)"

echo
echo "=== gzip integrity check ==="
fail=0
while IFS= read -r f; do
  if gzip -t "$f" 2>/dev/null; then
    echo "  OK      $(basename "$f")"
  else
    echo "  CORRUPT $(basename "$f")"; fail=1
  fi
done < <(find "$DEST" -name "*.fastq.gz")
[[ "$fail" -eq 0 ]] && echo "All files passed." || { echo "Re-download the corrupt files."; exit 1; }

echo
echo "Done. Data is in $DEST (gitignored)."
