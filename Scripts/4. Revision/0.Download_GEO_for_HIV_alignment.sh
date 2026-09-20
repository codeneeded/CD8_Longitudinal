#!/usr/bin/env bash
################################################################################
# Revision step 0 — download raw GEX FASTQs from GEO/SRA for HIV realignment
#
# ANSWERS: Reviewer 1, major point 3 — detect HIV-infected (vRNA+) cells by
# aligning to a combined human-HIV reference. That needs raw FASTQs; only the
# CP018 FACT-Seq run (73 GB) is on this machine.
#
# WHAT IT DOWNLOADS
#   26 Gene Expression libraries covering 24 of the 26 samples in
#   TARA_ALL_annotated_final.qs2. ADT/TCR/BCR are NOT downloaded: vRNA detection
#   only needs GEX, and skipping them saves roughly three quarters of the volume.
#   SA-AH-029 and SA-TY-021 are excluded — different cohort, not in this paper.
#
#   Series (from Reviewer_Comments/Tables/geo_sample_map.csv):
#     GSE254645   7 GEX   public
#     GSE306003   6 GEX   public
#     GSE328561  13 GEX   embargoed until 2026-09-30  <- may need a token
#     GSE328716   -       embargoed until 2027-07-26   CP003 FACT-Seq, already local
#
# WHERE IT WRITES — two drives, deliberately
#   RAW_ROOT  raw FASTQs. A SATA SSD is fine: the FASTQ read path is gated by
#             gzip decompression, measured at ~104 MB/s of compressed input per
#             stream on this project's 47 GB R2 file. SATA III gives ~550 MB/s,
#             so one stream uses about a fifth of it. The whole CP018 alignment
#             averaged 5.7 MB/s of I/O over ~10 h and was CPU-bound on STAR.
#   SCRATCH   the cellranger pipestance. Keep on NVMe. CP018's was 15,811 files
#             with small scattered reads/writes plus 47 GB of per-sample BAM —
#             the random-access pattern where NVMe actually wins.
#
# USAGE
#   1. Set RAW_ROOT below (or export it) to your data drive.
#   2. bash "Scripts/4. Revision/0.Download_GEO_for_HIV_alignment.sh" --check
#   3. bash "Scripts/4. Revision/0.Download_GEO_for_HIV_alignment.sh" --dry-run
#   4. nohup bash "Scripts/4. Revision/0.Download_GEO_for_HIV_alignment.sh" \
#        > ~/Desktop/Projects/geo_download.log 2>&1 &
#
#   Resumable: already-complete samples are skipped, so re-running after an
#   interruption picks up where it stopped.
################################################################################
set -uo pipefail

# ── EDIT: your data drive ─────────────────────────────────────────────────────
# Point this at the SATA SSD. Everything lands under $RAW_ROOT/TARA_GEO_raw/.
RAW_ROOT="${RAW_ROOT:-/mnt/datadrive}"

# NVMe scratch for the alignment that follows. Not used by this script except
# for the space preflight, so you find out now rather than at hour six.
SCRATCH="${SCRATCH:-$HOME/Desktop/Projects/cellranger_out}"

THREADS="${THREADS:-8}"
# ──────────────────────────────────────────────────────────────────────────────

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MAP="$REPO/Reviewer_Comments/Tables/gex_download_list.txt"

# Folder layout: one directory per GEO series, then per sample, so a partial
# download is obvious and a re-run is idempotent.
DEST="$RAW_ROOT/TARA_GEO_raw"
SRA_CACHE="$DEST/_sra_cache"
LOGDIR="$DEST/_logs"

MODE="run"
ONLY=""
for a in "$@"; do
  case "$a" in
    --check)   MODE="check" ;;
    --dry-run) MODE="dry" ;;
    --list)    MODE="list" ;;
    --only=*)  ONLY="${a#--only=}" ;;
    *) echo "unknown flag: $a"; echo "usage: $0 [--check|--dry-run|--list] [--only=SAMPLE]"; exit 1 ;;
  esac
done

say() { printf '%s\n' "$*"; }
hr()  { printf '%s\n' "------------------------------------------------------------"; }

################################################################################
# Preflight
################################################################################
hr; say "PREFLIGHT"; hr

say "tools:"
MISSING=()
for t in prefetch fasterq-dump esearch efetch; do
  if command -v "$t" >/dev/null 2>&1; then
    say "  $t  $(command -v "$t")"
  else
    say "  $t  MISSING"; MISSING+=("$t")
  fi
done
if ((${#MISSING[@]})); then
  cat <<'EOS'

  Install:
    conda create -n geodl -c conda-forge -c bioconda sra-tools entrez-direct
    conda activate geodl

  Then run `vdb-config --interactive` ONCE and set the cache directory to a
  path on the data drive. If you skip this, prefetch fills $HOME/ncbi and can
  exhaust the boot drive partway through a 1 TB download.
EOS
  [[ "$MODE" == "run" ]] && exit 1
fi

say ""
say "sample manifest:"
if [[ -f "$MAP" ]]; then
  say "  $MAP"
  say "  $(wc -l < "$MAP") GEX libraries listed"
  if [[ -n "$ONLY" ]]; then
    n_only=$(awk -F'\t' -v s="$ONLY" '$3==s' "$MAP" | wc -l)
    if [[ "$n_only" -eq 0 ]]; then
      say ""
      say "  --only=$ONLY matches NOTHING in the manifest. Available samples:"
      awk -F'\t' '{printf "    %s\n", $3}' "$MAP" | sort | sed 's/^/  /'
      exit 1
    fi
    say "  --only=$ONLY  -> restricting to $n_only library"
  fi
else
  say "  MISSING: $MAP"
  say "  Expected the manifest written by the revision prep step."
  [[ "$MODE" == "run" ]] && exit 1
fi

say ""
say "destinations:"
if [[ -d "$RAW_ROOT" ]]; then
  avail=$(df -BG --output=avail "$RAW_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
  say "  raw FASTQ  $DEST"
  say "             $(df -h "$RAW_ROOT" | tail -1 | awk '{print $4}') free on $(df -h "$RAW_ROOT" | tail -1 | awk '{print $1}')"
  # 26 GEX libraries at roughly 25-50 GB each
  if [[ -n "$avail" && "$avail" -lt 900 ]]; then
    say "             WARNING: 26 GEX libraries need roughly 700 GB - 1.3 TB."
    say "             Only ${avail} GB free."
  fi
else
  say "  raw FASTQ  $RAW_ROOT DOES NOT EXIST"
  say "             Set RAW_ROOT at the top of this script, or:"
  say "               RAW_ROOT=/your/data/drive bash \"\$0\""
  [[ "$MODE" == "run" ]] && exit 1
fi

if [[ -d "$SCRATCH" ]]; then
  s=$(df -BG --output=avail "$SCRATCH" 2>/dev/null | tail -1 | tr -dc '0-9')
  say "  pipestance $SCRATCH"
  say "             ${s} GB free (want 700+ for 20-plus samples on NVMe)"
fi

say ""
say "SRA cache location (must be on the data drive, not \$HOME):"
# `vdb-config --cfg` output format varies between sra-tools builds, so read the
# single key directly and fall back to the config file. Getting this wrong fills
# the boot drive partway through a terabyte download.
cache="$(vdb-config -o n /repository/user/main/public/root 2>/dev/null \
         | tr -d '"' | sed 's/.*= *//' | tail -1)"
if [[ -z "$cache" ]]; then
  cache="$(grep -hoE '/repository/user/main/public/root *= *"[^"]*"' \
            "$HOME/.ncbi/user-settings.mkfg" 2>/dev/null \
            | sed 's/.*= *"//; s/"$//' | tail -1)"
fi
if [[ -n "$cache" ]]; then
  say "  $cache"
  # Warn if the cache is NOT on the same filesystem as the data drive.
  cdev=$(df --output=source "$cache" 2>/dev/null | tail -1)
  rdev=$(df --output=source "$RAW_ROOT" 2>/dev/null | tail -1)
  if [[ -n "$cdev" && -n "$rdev" && "$cdev" != "$rdev" ]]; then
    say "  WARNING: cache is on $cdev but RAW_ROOT is on $rdev."
    say "           prefetch stages the full .sra there before conversion, so"
    say "           this will fill the wrong drive. Fix with:"
    say "             vdb-config --set /repository/user/main/public/root=$DEST/_sra_cache"
  else
    say "  on the same filesystem as RAW_ROOT - correct"
  fi
else
  say "  NOT SET. prefetch will default to \$HOME/ncbi on the NVMe and can"
  say "  exhaust the boot drive. Set it with:"
  say "    vdb-config --set /repository/user/main/public/root=$DEST/_sra_cache"
fi

if [[ "$MODE" == "check" ]]; then
  say ""; say "Preflight only. Re-run with --dry-run to see the download plan."
  exit 0
fi

mkdir -p "$DEST" "$SRA_CACHE" "$LOGDIR" 2>/dev/null || {
  say "Cannot create $DEST — check permissions on the data drive."; exit 1; }

# Refuse to run twice at once. Concurrent invocations fight over the SRA cache
# and prefetch's own locks, and produce interleaved, unattributable logs.
LOCK="$DEST/.download_lock"
if [[ -f "$LOCK" ]]; then
  other=$(cat "$LOCK" 2>/dev/null)
  if [[ -n "$other" ]] && kill -0 "$other" 2>/dev/null; then
    hr
    say "ALREADY RUNNING - pid $other"
    say ""
    say "  Another invocation of this script is live. Two at once corrupt the"
    say "  shared SRA cache. Stop it first:"
    say "    kill $other"
    say "    pkill -f 'prefetch|fasterq-dump'"
    say ""
    say "  Then re-run. Completed samples are skipped, so nothing is lost."
    hr
    exit 1
  fi
  say "  (clearing stale lock from dead pid $other)"
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT INT TERM

################################################################################
# Resolve GSM -> SRR runs
################################################################################
hr; say "RESOLVING RUNS"; hr

# Per-invocation run table. A shared _runs.tsv is actively dangerous: two
# concurrent invocations truncate the SAME inode, and a running loop reading it
# on fd 3 keeps reading at its old offset - so it silently starts processing the
# other invocation's rows.
RUNS="$DEST/_runs${ONLY:+_$ONLY}.tsv"
: > "$RUNS"
PRIVATE=()

# Read the manifest on fd 3, NOT stdin. esearch reads stdin itself, so with a
# plain `while read` it swallows the rest of the manifest and the loop stops
# after one line. `< /dev/null` on esearch is the second half of the fix.
while IFS=$'\t' read -r -u 3 gsm gse sample; do
  [[ -z "${gsm:-}" ]] && continue
  # --only: skip everything but the named sample. Used to pilot one library
  # before committing to the full 400-700 GB.
  [[ -n "$ONLY" && "$sample" != "$ONLY" ]] && continue
  # GSM -> SRX -> SRR, carrying the sample label so filenames stay meaningful
  info=$(esearch -db sra -query "$gsm" < /dev/null 2>/dev/null \
         | efetch -format runinfo 2>/dev/null)
  n=0
  while IFS=, read -r run _rest; do
    [[ "$run" =~ ^[SED]RR ]] || continue
    printf '%s\t%s\t%s\t%s\n' "$gse" "$sample" "$gsm" "$run" >> "$RUNS"
    n=$((n+1))
  done <<< "$info"
  if [[ "$n" -eq 0 ]]; then
    PRIVATE+=("$gsm ($gse, $sample)")
    printf '  %-12s %-14s %s\n' "$gsm" "$sample" "no public run - embargoed?"
  else
    printf '  %-12s %-14s %s run(s)\n' "$gsm" "$sample" "$n"
  fi
done 3< "$MAP"

TOTAL=$(wc -l < "$RUNS")
say ""
say "resolved $TOTAL runs"

if ((${#PRIVATE[@]})); then
  say ""
  hr
  say "EMBARGOED — ${#PRIVATE[@]} libraries did not resolve anonymously:"
  printf '  %s\n' "${PRIVATE[@]}"
  hr
  cat <<'EOS'
  A GEO web login does not authenticate prefetch: SRA is a separate service
  with its own access path. For a series still under embargo:

    (a) SRA runs may already be public even while the GEO series page is not.
        Check one directly:
            esearch -db sra -query GSM9685302 < /dev/null | efetch -format runinfo
        If that returns an SRR, the data is fetchable now and this script will
        pick it up on the next run.

    (b) Otherwise use the GEO reviewer/private token: open the series with the
        token URL while logged in to your browser, download the supplementary
        FASTQs, and place them in the layout this script produces so the rest
        of the pipeline finds them unchanged:
            $DEST/<GSE>/<sample>/<sample>_S1_L001_R1_001.fastq.gz
                                 <sample>_S1_L001_R2_001.fastq.gz
                                 <sample>_S1_L001_R3_001.fastq.gz   (if present)
        Then verify with:
            for f in $DEST/*/*/*.fastq.gz; do gzip -t "$f" || echo "BAD $f"; done

    (c) Or release the embargo early from your GEO submitter account, which
        makes the runs public immediately.

  Public series still download below; re-run after adding the embargoed ones
  and completed samples are skipped.
EOS
fi

if [[ "$TOTAL" -eq 0 ]]; then
  say ""; say "Nothing to download. Resolve access first."; exit 1
fi

if [[ "$MODE" == "dry" ]]; then
  hr; say "DOWNLOAD PLAN"; hr
  awk -F'\t' '{printf "  %-11s %-14s %-12s %s\n", $1, $2, $3, $4}' "$RUNS"
  say ""
  say "  into: $DEST/<GSE>/<sample>/"
  say "  re-run without --dry-run to start"
  exit 0
fi

################################################################################
# Download
################################################################################
hr; say "DOWNLOADING"; hr

# Libraries are split across flowcell lanes: one GSM can have 2-4 SRR runs, all
# with the same LibraryName. They are all part of ONE library and cellranger
# needs every one, as a separate lane. So iterate SAMPLES (not runs) and give
# each run its own L00n, or the second run overwrites the first.
SAMPLES=$(awk -F'\t' '{print $1"\t"$2}' "$RUNS" | sort -u)
NSAMP=$(printf '%s\n' "$SAMPLES" | grep -c .)
say "  $NSAMP samples, $TOTAL runs total"
say ""

i=0; failed=0
while IFS=$'\t' read -r gse sample; do
  i=$((i+1))
  out="$DEST/$gse/$sample"
  done_marker="$out/.complete"

  # Every run belonging to this sample, in accession order so lane numbering
  # is deterministic across re-runs.
  mapfile -t sruns < <(awk -F'\t' -v g="$gse" -v s="$sample" \
                         '$1==g && $2==s {print $4}' "$RUNS" | sort)
  nrun=${#sruns[@]}

  if [[ -f "$done_marker" ]]; then
    printf '  [%2d/%d] %-14s already complete (%d lanes)\n' \
      "$i" "$NSAMP" "$sample" "$nrun"
    continue
  fi
  mkdir -p "$out"
  printf '  [%2d/%d] %-14s %d run(s): %s\n' \
    "$i" "$NSAMP" "$sample" "$nrun" "${sruns[*]}"

  sample_ok=1
  lane=0
  for run in "${sruns[@]}"; do
    lane=$((lane+1))
    L=$(printf 'L%03d' "$lane")
    # If this lane is already converted, skip it - makes a partial sample
    # resumable without re-downloading what is already on disk.
    if ls "$out/${sample}_S1_${L}_R1_001.fastq.gz" >/dev/null 2>&1; then
      printf '           %s %s already present\n' "$run" "$L"
      continue
    fi
    printf '           %s -> %s ... ' "$run" "$L"

    # --max-size u removes the 20 GB default cap. Without it prefetch reports
    # success on a TRUNCATED 10x run, which then fails mid-alignment.
    if ! prefetch --max-size u --output-directory "$SRA_CACHE" "$run" \
          < /dev/null > "$LOGDIR/${run}_prefetch.log" 2>&1; then
      say "prefetch FAILED (see $LOGDIR/${run}_prefetch.log)"
      sample_ok=0; failed=$((failed+1)); continue
    fi

    # --include-technical keeps the 10x barcode/UMI read, which cellranger needs.
    if ! fasterq-dump --split-files --include-technical --threads "$THREADS" \
          --temp "$SRA_CACHE" --outdir "$out" "$SRA_CACHE/$run" \
          < /dev/null > "$LOGDIR/${run}_dump.log" 2>&1; then
      say "fasterq-dump FAILED (see $LOGDIR/${run}_dump.log)"
      sample_ok=0; failed=$((failed+1)); continue
    fi

    # Rename to the Illumina convention cellranger requires, one lane per run.
    nread=0
    for f in "$out/${run}"_*.fastq; do
      [[ -e "$f" ]] || continue
      nread=$((nread+1))
      mv "$f" "$out/${sample}_S1_${L}_R${nread}_001.fastq"
    done
    if [[ "$nread" -eq 0 ]]; then
      say "no fastq produced"; sample_ok=0; failed=$((failed+1)); continue
    fi
    gzip -f "$out/${sample}_S1_${L}"_R*_001.fastq 2>/dev/null

    # Verify this lane before moving on - a truncated gzip wastes an alignment.
    lane_ok=1
    for f in "$out/${sample}_S1_${L}"_R*_001.fastq.gz; do
      gzip -t "$f" 2>/dev/null || lane_ok=0
    done
    if [[ "$lane_ok" -eq 1 ]]; then
      say "ok ($nread reads, $(du -sh "$out" | cut -f1) cumulative)"
      rm -rf "$SRA_CACHE/$run"
    else
      say "CORRUPT after conversion"; sample_ok=0; failed=$((failed+1))
    fi
  done

  # Only mark the SAMPLE complete when every one of its lanes converted cleanly.
  nlanes=$(ls "$out"/${sample}_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  if [[ "$sample_ok" -eq 1 && "$nlanes" -eq "$nrun" ]]; then
    touch "$done_marker"
    printf '           complete: %d/%d lanes, %s\n' \
      "$nlanes" "$nrun" "$(du -sh "$out" | cut -f1)"
  else
    printf '           INCOMPLETE: %d/%d lanes present - re-run to finish\n' \
      "$nlanes" "$nrun"
  fi
done <<< "$SAMPLES"

hr; say "VERIFYING"; hr

nfq=0; bad=0
while IFS= read -r f; do
  nfq=$((nfq+1))
  gzip -t "$f" 2>/dev/null || { say "  CORRUPT: $f"; bad=$((bad+1)); }
done < <(find "$DEST" -name "*.fastq.gz")

complete=$(find "$DEST" -name ".complete" | wc -l)
say "  samples complete:  $complete"
say "  fastq.gz files:    $nfq"
say "  corrupt:           $bad"
say "  failed downloads:  $failed"
say "  total size:        $(du -sh "$DEST" 2>/dev/null | cut -f1)"

# Per-sample table, so a gap is visible rather than assumed absent
say ""
say "  per sample (lanes present / lanes expected):"
for d in "$DEST"/GSE*/*/; do
  [[ -d "$d" ]] || continue
  s=$(basename "$d"); g=$(basename "$(dirname "$d")")
  have=$(ls "$d"/${s}_S1_L*_R1_001.fastq.gz 2>/dev/null | wc -l)
  want=$(awk -F'\t' -v g="$g" -v s="$s" '$1==g && $2==s' "$RUNS" | wc -l)
  m=$([[ -f "$d/.complete" ]] && echo ok || echo INCOMPLETE)
  printf '    %-11s %-14s %2d/%-2d lanes  %-9s %s\n' \
    "$g" "$s" "$have" "$want" "$(du -sh "$d" 2>/dev/null | cut -f1)" "$m"
done

say ""
if [[ "$bad" -eq 0 && "$failed" -eq 0 ]]; then
  say "All downloads verified. Next: build the combined reference."
  say "  Scripts/4. Revision/1.Build_HIV_Reference.sh"
else
  say "Finished with $failed failures and $bad corrupt files."
  say "Re-run to retry — complete samples are skipped."
  exit 1
fi
