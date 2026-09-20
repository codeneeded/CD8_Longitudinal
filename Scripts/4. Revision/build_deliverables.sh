#!/usr/bin/env bash
# Rebuild both reviewer deliverables and GATE them on citation consistency.
# Refuses to leave stale citations in place: the gate exits non-zero and this
# script stops, so a bad citation cannot be shipped silently.
set -euo pipefail
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RS="${RSCRIPT:-Rscript}"
PY="${PYTHON:-python3}"
cd "$S"
echo "== revised figures PDF =="
"$RS" 98.Build_Revision_Figure_PDF.R | grep -E "Pages:|MISSING|every source" || true
echo "== response document =="
"$RS" 99.Build_Response_Document.R | grep -E "drafted so far|Pages written|missing figure" || true
echo "== citation consistency gate =="
"$PY" 97.Check_Citations.py
echo "== deliverables are consistent =="
