#!/usr/bin/env bash
################################################################################
# Revision — download the ENTRY-VISIT libraries for HIV-1 realignment
#
# ANSWERS: Reviewer 1, major point 3, at cohort scale rather than one pilot.
#
# SCOPE — 9 samples, decided deliberately. From the 16 entry-visit libraries:
#   - CP006_1m   already aligned in the pilot (22 GB of output on disk) -> skip
#   - CE021/CE025/CE037_entry (HEU) and CS005/CS015_entry (HUU) are UNINFECTED,
#     so there is no viral RNA to find -> skip. (They would have made useful
#     specificity controls, but two of them are embargoed anyway; the pilot's
#     ambient-background control already bounds the false-positive rate.)
#   - CP020_1m   in GSE328561, which has NO SRA records at all -> cannot fetch
#
#   Left: CP002_1m CP003_1m CP011_2m CP013_1m CP016_2m CP017_2m CP018_1m
#         CP022_1m CP042_1m       = 384 GB of SRA archive, 27 runs
#
# WHY THE EMBARGOED THREE CANNOT BE FETCHED WITH A TOKEN
#   Their GSMs return ZERO SRA records — not records with restricted consent.
#   A reviewer token opens a private GEO page for viewing metadata and
#   processed files; it does not create SRA runs, and prefetch can only fetch
#   runs that exist. Releasing GSE328561 (the submitting author's call) is what
#   makes them fetchable.
#
# DISK — this is the constraint, so read it before starting.
#   Archive 384 GB; FASTQ expands to roughly 2.5-3x that. Both together do not
#   fit in 1.2 TB. Default mode is therefore STREAM: one sample at a time,
#   download -> convert -> verify -> (align) -> free, so peak stays under
#   ~200 GB. Pass --keep to retain everything if you have the room.
#
# USAGE
#   bash "Scripts/4. Revision/0b.Download_Entry_HIV.sh" --check      # gates only
#   bash "Scripts/4. Revision/0b.Download_Entry_HIV.sh" --dry-run    # resolve runs
#   bash "Scripts/4. Revision/0b.Download_Entry_HIV.sh"              # stream all 9
#   bash "Scripts/4. Revision/0b.Download_Entry_HIV.sh" --only=CP022_1m
#   RAW_ROOT=/mnt/data bash "Scripts/4. Revision/0b.Download_Entry_HIV.sh"
################################################################################
set -uo pipefail

REPO="${REPO:-$HOME/Desktop/Projects/Repositories/CD8_Longitudinal}"
RAW_ROOT="${RAW_ROOT:-/mnt/data}"
DEST="$RAW_ROOT/TARA_GEO_raw"
SRA_CACHE="$DEST/_sra_cache"
LOG="${LOG:-$HOME/Desktop/Projects/entry_hiv_download.log}"

MODE="run"; ONLY=""; KEEP=0; DEEP=0
for a in "$@"; do
  case "$a" in
    --check)   MODE="check" ;;
    --status)  MODE="status" ;;
    --deep)    DEEP=1 ;;
    --dry-run) MODE="dry" ;;
    --list)    MODE="list" ;;
    --keep)    KEEP=1 ;;
    --only=*)  ONLY="${a#--only=}" ;;
    *) echo "unknown flag: $a"; exit 1 ;;
  esac
done

say() { echo "$*" | tee -a "$LOG"; }
hr()  { say "------------------------------------------------------------------"; }

# ── The 9 samples. GSM is what resolves to SRA runs; the sample name is what
#    cellranger's --sample must match, so both are carried explicitly.
#    Ordered by viral load descending: if the run is interrupted, the samples
#    most likely to carry detectable virus are already done.
read -r -d '' MANIFEST <<'EOF'
GSM9191480	GSE306003	CP022_1m	5075764
GSM8047705	GSE254645	CP002_1m	4284389
GSM9191476	GSE306003	CP017_2m	3167384
GSM9191472	GSE306003	CP016_2m	1978332
GSM8047709	GSE254645	CP003_1m	656769
GSM8047717	GSE254645	CP018_1m	176970
GSM9191464	GSE306003	CP011_2m	36965
GSM9191484	GSE306003	CP042_1m	6113
GSM9191468	GSE306003	CP013_1m	3434
EOF

if [[ "$MODE" == "list" ]]; then
  printf '%-12s %-11s %-11s %12s\n' GSM SERIES SAMPLE VIRAL_LOAD
  printf '%s\n' "$MANIFEST" | while IFS=$'\t' read -r g s n v; do
    printf '%-12s %-11s %-11s %12s\n' "$g" "$s" "$n" "$v"
  done
  exit 0
fi

mkdir -p "$(dirname "$LOG")"
if [[ "$MODE" == "status" ]]; then
  RUNS="$DEST/_entry_runs.tsv"   # defined here too: this runs before the
                                 # resolve step, which is where it normally comes from
  hr; say "DOWNLOAD STATUS   $(date '+%F %T')"; say "  root: $DEST"; hr
  say "  integrity: $([[ "$DEEP" == "1" ]] && echo 'full gzip -t (reads every byte)' \
                          || echo 'from the download-time gzip -t; pass --deep to re-check now')"
  printf '  %-11s %-11s %7s %8s  %-9s %s\n' SAMPLE SERIES LANES SIZE MARKER GZIP
  n_ok=0; n_bad=0; n_seen=0
  while IFS=$'\t' read -r gsm gse samp vl; do
    [[ -z "$samp" ]] && continue
    [[ -n "$ONLY" && "$samp" != "$ONLY" ]] && continue
    n_seen=$((n_seen+1))
    d="$DEST/$gse/$samp"
    # Expected lane count comes from the cached run table when it exists; with
    # no table we report what is present rather than inventing an expectation.
    want=0
    [[ -s "$RUNS" ]] && want=$(awk -F'\t' -v s="$samp" '$3==s' "$RUNS" | wc -l)
    have=$(ls "$d"/${samp}_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
    size=$(du -sh "$d" 2>/dev/null | cut -f1); size=${size:--}
    mark=$([[ -f "$d/.complete" ]] && echo complete || echo MISSING)
    # INTEGRITY. The download already ran `gzip -t` over every FASTQ in a sample
    # and wrote .complete only if all of them passed (see gz_ok/verify block
    # below), so the marker is not a bare claim - it records a full CRC check
    # made at write time. Re-reading 200-400 GB here would therefore re-prove
    # what the marker already asserts, and doing it while cellranger aligns off
    # the same drive steals I/O from the alignment.
    #
    # So: default reports WHEN integrity was established. --deep re-checks from
    # scratch, which is the right call only if the drive is suspect.
    #
    # A cheap middle option was tried and rejected: reading the 4-byte gzip
    # trailer. It cannot detect truncation - a file cut to 70% reported a
    # plausible non-zero uncompressed size while `gzip -t` correctly called it
    # CORRUPT. Reporting that as "ok" would have been worse than reporting
    # nothing.
    gz="-"
    if (( have > 0 )); then
      if [[ "$DEEP" == "1" ]]; then
        if gzip -t "$d"/${samp}_S1_L*_001.fastq.gz 2>/dev/null; then gz="valid"; else gz="CORRUPT"; fi
      elif [[ "$mark" == "complete" ]]; then
        gz="at-download"
      else
        gz="unverified"
      fi
    fi
    lanes=$([[ "$want" -gt 0 ]] && echo "$have/$want" || echo "$have/?")
    printf '  %-11s %-11s %7s %8s  %-9s %s\n' "$samp" "$gse" "$lanes" "$size" "$mark" "$gz"
    if [[ "$mark" == "complete" && ( "$gz" == "valid" || "$gz" == "at-download" ) ]] \
       && { [[ "$want" -eq 0 ]] || [[ "$have" -eq "$want" ]]; }; then
      n_ok=$((n_ok+1))
    else
      n_bad=$((n_bad+1))
    fi
  done <<< "$MANIFEST"
  hr
  say "  $n_ok of $n_seen samples complete and gzip-valid"
  if [[ "$n_seen" -eq 0 ]]; then
    say "  ** the manifest produced no rows - this is the bug that made SUMMARY print empty"
    exit 2
  fi
  if (( n_bad > 0 )); then
    say "  $n_bad not ready. Re-run the download; finished samples are skipped."
    exit 1
  fi
  say "  ready to align:  bash \"Scripts/4. Revision/2b.Align_Entry_HIV.sh\""
  exit 0
fi

hr; say "ENTRY-VISIT HIV REALIGNMENT — DOWNLOAD   $(date '+%F %T')"; hr

################################################################################
# Preflight. Everything that can stop the run, checked before the first byte.
################################################################################
fail=0

# --- tools -------------------------------------------------------------------
# conda install puts these in <prefix>/bin, but they will not be found if the
# shell cached its PATH lookups before the install (bash hashes command
# locations) or if the env is not active. Locate them ourselves and prepend,
# rather than making the user debug PATH.
if ! command -v prefetch >/dev/null 2>&1; then
  hash -r 2>/dev/null || true
  for pre in "${CONDA_PREFIX:-}" "$HOME/miniforge3" "$HOME/miniconda3" \
             "$HOME/anaconda3" "$HOME/mambaforge" /opt/conda; do
    [[ -n "$pre" && -x "$pre/bin/prefetch" ]] && { export PATH="$pre/bin:$PATH"; break; }
  done
fi
if ! command -v prefetch >/dev/null 2>&1; then
  # conda envs other than base
  for e in "$HOME"/miniforge3/envs/*/bin/prefetch "$HOME"/miniconda3/envs/*/bin/prefetch; do
    [[ -x "$e" ]] && { export PATH="$(dirname "$e"):$PATH"; break; }
  done
fi
# NCBI's own tarball layout, if installed that way
if ! command -v prefetch >/dev/null 2>&1; then
  for d in "$HOME"/sratoolkit*/bin "$HOME/Desktop/Projects/Tools"/sratoolkit*/bin; do
    [[ -x "$d/prefetch" ]] && { export PATH="$d:$PATH"; break; }
  done
fi

for t in prefetch fasterq-dump vdb-config; do
  if command -v "$t" >/dev/null 2>&1; then
    say "  ok    $t  $(command -v "$t")"
  elif [[ "$t" == "vdb-config" ]]; then
    # Only used to report the cache location. prefetch -O sets the output
    # directory explicitly, so a missing vdb-config is not fatal.
    say "  WARN  vdb-config not found (cache reporting disabled; download unaffected)"
  else
    say "  FAIL  $t not on PATH"
    say "        install:  conda install -c bioconda sra-tools"
    say "        then:     hash -r   (bash caches command lookups)"
    say "        found nothing in: \$CONDA_PREFIX, ~/miniforge3, ~/miniconda3,"
    say "                          ~/anaconda3, ~/mambaforge, /opt/conda,"
    say "                          ~/sratoolkit*/bin, Tools/sratoolkit*/bin"
    fail=1
  fi
done

# --- the SRA cache must NOT sit on the OS disk ------------------------------
# prefetch writes to the vdb cache, which defaults to $HOME. On a 384 GB
# download that fills the root filesystem and takes the machine down with it.
if command -v vdb-config >/dev/null 2>&1; then
  cur="$(vdb-config -o n NCBI_SETTINGS 2>/dev/null | head -1)"
  root_cache="$(grep -hoE '/repository/user/main/public/root *= *"[^"]*"' \
                  "$HOME/.ncbi/user-settings.mkfg" 2>/dev/null \
                | sed 's/.*= *"//; s/"$//' | head -1)"
  if [[ -z "$root_cache" ]]; then
    say "  WARN  no vdb cache root configured; prefetch will use \$HOME"
    say "        set it once:  vdb-config --prefetch-to-cwd"
    say "        or:           vdb-config -s /repository/user/main/public/root=\"$SRA_CACHE\""
  else
    say "  ok    vdb cache root: $root_cache"
    case "$root_cache" in
      "$RAW_ROOT"*) : ;;
      *) say "  WARN  cache is NOT under $RAW_ROOT — it will consume that filesystem instead" ;;
    esac
  fi
fi

# --- destination and space --------------------------------------------------
ARCHIVE_GB=384
if [[ -d "$RAW_ROOT" ]]; then
  avail=$(df -BG --output=avail "$RAW_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
  say "  ok    RAW_ROOT $RAW_ROOT  (${avail:-?} GB free)"
  need=$(( KEEP == 1 ? ARCHIVE_GB * 4 : 200 ))
  if [[ -n "${avail:-}" && "$avail" -lt "$need" ]]; then
    say "  FAIL  need ~${need} GB $([[ $KEEP == 1 ]] && echo '(--keep: archive + FASTQ)' || echo '(stream mode peak)'), have ${avail} GB"
    fail=1
  fi
else
  say "  FAIL  RAW_ROOT does not exist: $RAW_ROOT"
  say "        mount the data drive, or:  RAW_ROOT=/some/path bash $0"
  fail=1
fi

# --- the reference must exist AND carry both viral contigs ------------------
# A reference silently missing a viral contig produces a clean-looking zero.
IDX="${IDX:-$HOME/Desktop/Projects/refdata/human_hiv/GRCh38_HIV1}"
if [[ -s "$IDX/genes/genes.gtf.gz" ]]; then
  nvir=$(zcat "$IDX/genes/genes.gtf.gz" | awk '$1 ~ /^HIV1_/{print $1}' | sort -u | wc -l)
  if [[ "$nvir" -eq 2 ]]; then
    say "  ok    reference $IDX  (2 viral contigs in the GTF)"
  else
    say "  FAIL  reference has $nvir viral contigs in its GTF, expected 2"
    say "        a missing contig makes every zero meaningless — rebuild with 1.Build_HIV_Reference.sh"
    fail=1
  fi
else
  say "  WARN  reference not found at $IDX (needed for alignment, not download)"
fi

# --- network ----------------------------------------------------------------
if curl -s --max-time 20 -o /dev/null \
     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=GSM8047713"; then
  say "  ok    NCBI reachable"
else
  say "  FAIL  cannot reach NCBI eutils"
  fail=1
fi

hr
if [[ "$fail" -ne 0 ]]; then
  say "PREFLIGHT FAILED — nothing was downloaded. Fix the FAIL lines above."
  exit 1
fi
say "preflight passed"
[[ "$MODE" == "check" ]] && { say "(--check: stopping here)"; exit 0; }

################################################################################
# Resolve GSM -> SRA runs. One query per sample, cached to a file so a rerun
# does not re-query and cannot silently pick up a different run set.
################################################################################
mkdir -p "$DEST"
RUNS="$DEST/_entry_runs.tsv"

################################################################################
# --status: report what is on disk and stop. Nothing is downloaded or resolved.
#
# WHY THIS EXISTS
#   The end-of-run SUMMARY printed no rows on at least one invocation, so there
#   was no way to confirm the set was complete before committing several days to
#   alignment. --check exits earlier than any completeness test, so it cannot
#   answer this either. This reads the FASTQs themselves rather than trusting a
#   marker, and exits non-zero if any sample is short, so it can gate the
#   aligner:  bash 0b.Download_Entry_HIV.sh --status && bash 2b.Align_Entry_HIV.sh
################################################################################

# The expected sample count is known, so an incomplete run table can be
# detected rather than trusted. This matters: the first version cached a table
# containing only 4 of 9 samples after NCBI rate-limited the queries, and would
# then have silently downloaded 4 samples and reported success.
N_EXPECTED=$(printf '%s\n' "$MANIFEST" | grep -c .)
[[ -n "$ONLY" ]] && N_EXPECTED=1

n_cached=0
[[ -s "$RUNS" ]] && n_cached=$(awk -F'\t' '{print $3}' "$RUNS" | sort -u | wc -l)
if [[ "$n_cached" -gt 0 && "$n_cached" -lt "$N_EXPECTED" ]]; then
  say "cached run table covers only $n_cached of $N_EXPECTED samples — discarding it"
  say "  (a partial table is what a rate-limited resolve leaves behind)"
  mv "$RUNS" "$RUNS.partial.$(date +%s)" 2>/dev/null || rm -f "$RUNS"
fi

if [[ ! -s "$RUNS" ]]; then
  hr; say "RESOLVING RUNS"; hr
  # NCBI allows 3 requests/sec anonymously, 10/sec with an API key. Each sample
  # costs TWO requests (esearch + efetch), so an unpaced loop over 9 samples is
  # throttled and esearch returns an EMPTY result set rather than an error.
  # That is indistinguishable from "no such record" unless you retry - which is
  # exactly how five public samples got reported as embargoed.
  Q=""
  [[ -n "${NCBI_API_KEY:-}" ]] && Q="&api_key=${NCBI_API_KEY}"
  [[ -n "${NCBI_EMAIL:-}" ]]   && Q="${Q}&email=${NCBI_EMAIL}&tool=cd8_longitudinal_revision"
  PACE="${PACE:-1}"            # seconds between samples
  TRIES="${TRIES:-4}"

  : > "$RUNS"
  while IFS=$'\t' read -r gsm gse samp vl; do
    [[ -z "$gsm" ]] && continue
    [[ -n "$ONLY" && "$samp" != "$ONLY" ]] && continue

    uid=""
    for try in $(seq 1 "$TRIES"); do
      # </dev/null: esearch/efetch read stdin and would consume the manifest
      uid=$(curl -s --max-time 30 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${gsm}&retmax=5${Q}" \
        </dev/null | grep -oE '<Id>[0-9]+</Id>' | sed 's/<[^>]*>//g' | head -1)
      [[ -n "$uid" ]] && break
      back=$((try * 3))
      say "  $samp: empty result (attempt $try/$TRIES) — waiting ${back}s"
      sleep "$back"
    done

    if [[ -z "$uid" ]]; then
      say "  $samp: NO SRA RECORD for $gsm after $TRIES attempts"
      say "          This is either a genuine embargo or sustained throttling."
      say "          Confirm by hand before believing it:"
      say "            curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$gsm'"
      continue
    fi

    curl -s --max-time 60 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${uid}&rettype=runinfo&retmode=text${Q}" \
      </dev/null \
      | awk -F, -v g="$gsm" -v s="$samp" -v e="$gse" '
          NR==1{for(i=1;i<=NF;i++)h[$i]=i}
          NR>1 && $1!="" {print e"\t"g"\t"s"\t"$h["Run"]"\t"$h["size_MB"]"\t"$h["Consent"]}' \
      >> "$RUNS"

    n=$(awk -F'\t' -v s="$samp" '$3==s' "$RUNS" | wc -l)
    if [[ "$n" -eq 0 ]]; then
      say "  $samp: uid $uid resolved but runinfo returned no rows — retrying once"
      sleep 5
      curl -s --max-time 60 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=${uid}&rettype=runinfo&retmode=text${Q}" \
        </dev/null \
        | awk -F, -v g="$gsm" -v s="$samp" -v e="$gse" '
            NR==1{for(i=1;i<=NF;i++)h[$i]=i}
            NR>1 && $1!="" {print e"\t"g"\t"s"\t"$h["Run"]"\t"$h["size_MB"]"\t"$h["Consent"]}' \
        >> "$RUNS"
      n=$(awk -F'\t' -v s="$samp" '$3==s' "$RUNS" | wc -l)
    fi
    say "  $samp: $n run(s)"
    sleep "$PACE"
  done <<< "$MANIFEST"

else
  say "reusing resolved run table: $RUNS"
fi

# A sample with zero runs must be reported, not silently absent.
printf '%s\n' "$MANIFEST" | while IFS=$'\t' read -r gsm gse samp vl; do
  [[ -n "$ONLY" && "$samp" != "$ONLY" ]] && continue
  n=$(awk -F'\t' -v s="$samp" '$3==s' "$RUNS" | wc -l)
  [[ "$n" -eq 0 ]] && say "  ** $samp resolved to NO runs — it will not be downloaded"
done

hr
say "runs to fetch:"
awk -F'\t' '{n[$3]++; mb[$3]+=$5} END {for (k in n) printf "  %-11s %d run(s)  %6.1f GB\n", k, n[k], mb[k]/1024}' \
  "$RUNS" | sort
awk -F'\t' '{s+=$5} END {printf "  TOTAL %.0f GB archive\n", s/1024}' "$RUNS"
say ""
n_res=$(awk -F'\t' '{print $3}' "$RUNS" | sort -u | wc -l)
if [[ "$n_res" -lt "$N_EXPECTED" ]]; then
  say "  ** RESOLVED ONLY $n_res OF $N_EXPECTED SAMPLES **"
  say "  Downloading now would quietly give you a partial dataset. Either:"
  say "    - rerun this script (the partial table has been set aside, so it"
  say "      re-resolves from scratch), ideally with an API key:"
  say "        NCBI_API_KEY=<key> bash $0 --dry-run"
  say "      (free from https://www.ncbi.nlm.nih.gov/account/ -> API Key Management)"
  say "    - or accept the subset deliberately:  ACCEPT_PARTIAL=1 bash $0"
  if [[ "${ACCEPT_PARTIAL:-0}" != "1" && "$MODE" != "dry" ]]; then
    say "  refusing to download a partial set (set ACCEPT_PARTIAL=1 to override)"
    exit 1
  fi
fi

if [[ "$MODE" == "dry" ]]; then
  say "(--dry-run: resolved only, nothing downloaded)"
  say "Sanity-check the run counts above: 2-4 per sample is a lane split and"
  say "expected. A dozen would mean the query pulled in sibling libraries."
  exit 0
fi

################################################################################
# Download. Per sample: prefetch -> fasterq-dump -> rename to 10x convention
# -> gzip -> verify -> mark complete. Idempotent: a sample with .complete is
# skipped, so an interrupted run resumes rather than restarting.
################################################################################
gz_ok() { gzip -t "$1" 2>/dev/null; }

samples=$(awk -F'\t' '{print $3}' "$RUNS" | sort -u)
for samp in $samples; do
  gse=$(awk -F'\t' -v s="$samp" '$3==s{print $1; exit}' "$RUNS")
  out="$DEST/$gse/$samp"
  want=$(awk -F'\t' -v s="$samp" '$3==s' "$RUNS" | wc -l)

  if [[ -f "$out/.complete" ]]; then
    say "  [skip] $samp already complete ($(du -sh "$out" 2>/dev/null | cut -f1))"
    continue
  fi

  hr; say "$samp  ($want run(s), $gse)   $(date '+%H:%M:%S')"; hr
  mkdir -p "$out"

  lane=0
  ok=1
  while IFS=$'\t' read -r e g s run mb consent; do
    [[ "$s" != "$samp" ]] && continue
    lane=$((lane+1))
    L=$(printf 'L%03d' "$lane")

    # already converted this lane?
    if compgen -G "$out/${samp}_S1_${L}_R1_001.fastq.gz" >/dev/null; then
      say "  [have] $run -> $L"
      continue
    fi

    say "  prefetch $run (${mb} MB, consent=$consent)"
    if ! prefetch --max-size u -O "$SRA_CACHE" "$run" </dev/null >>"$LOG" 2>&1; then
      say "  FAILED prefetch $run"; ok=0; break
    fi

    say "  fasterq-dump $run -> $L"
    # --split-files gives R1/R2 separately, which is what 10x needs.
    if ! fasterq-dump --split-files --include-technical -O "$out" \
          -t "$out/_tmp" "$SRA_CACHE/$run" </dev/null >>"$LOG" 2>&1; then
      say "  FAILED fasterq-dump $run"; ok=0; break
    fi

    # 10x naming: cellranger matches --sample against the prefix.
    for r in 1 2 3; do
      src="$out/${run}_${r}.fastq"
      [[ -f "$src" ]] || continue
      case "$r" in
        1) tag="R1" ;; 2) tag="R2" ;; 3) tag="I1" ;;
      esac
      mv "$src" "$out/${samp}_S1_${L}_${tag}_001.fastq"
    done
    rm -rf "$out/_tmp"

    say "  gzip $L"
    for f in "$out/${samp}_S1_${L}"_*_001.fastq; do
      [[ -f "$f" ]] || continue
      gzip -f "$f" || { say "  FAILED gzip $f"; ok=0; break 2; }
    done

    # free the archive as we go — this is what keeps stream mode inside 200 GB
    if [[ "$KEEP" -eq 0 ]]; then
      rm -rf "$SRA_CACHE/$run"
      say "  freed archive for $run"
    fi
  done <<< "$(cat "$RUNS")"

  # ---- verify before marking complete ----
  have=$(ls "$out"/${samp}_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  bad=0
  for f in "$out"/*.fastq.gz; do [[ -f "$f" ]] && { gz_ok "$f" || { say "  CORRUPT $(basename "$f")"; bad=1; }; }; done

  if [[ "$ok" -eq 1 && "$have" -eq "$want" && "$bad" -eq 0 ]]; then
    touch "$out/.complete"
    say "  DONE $samp: $have/$want lanes, $(du -sh "$out" | cut -f1), all gzip-valid"
  else
    say "  INCOMPLETE $samp: $have/$want lanes, corrupt=$bad, errors=$((1-ok))"
    say "             rerun this script; completed lanes are not re-fetched"
  fi
done

################################################################################
hr; say "SUMMARY   $(date '+%F %T')"; hr
printf '%s\n' "$MANIFEST" | while IFS=$'\t' read -r gsm gse samp vl; do
  [[ -n "$ONLY" && "$samp" != "$ONLY" ]] && continue
  d="$DEST/$gse/$samp"
  want=$(awk -F'\t' -v s="$samp" '$3==s' "$RUNS" 2>/dev/null | wc -l)
  have=$(ls "$d"/${samp}_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  st=$([[ -f "$d/.complete" ]] && echo complete || echo INCOMPLETE)
  printf '  %-11s VL=%9s  %d/%-2d lanes  %7s  %s\n' \
    "$samp" "$vl" "$have" "$want" "$(du -sh "$d" 2>/dev/null | cut -f1)" "$st"
done
say ""
say "Next: align. One sample at a time, ~4-8 h each:"
say "  bash \"Scripts/4. Revision/2b.Align_Entry_HIV.sh\""
say ""
say "Do NOT read a viral count as an infected-cell count until"
say "3.Score_HIV_Detection.R has run its ambient, lineage and host controls."
