#!/bin/bash
# Step 2 (bias): distance distribution of fragment midpoints around bias k-mer sites (raw V-plot data)
#
# For each k-mer bed in the bias directory:
#   1) (optionally) subsample -n sites, take the single-base strand endpoint (cut site), sort;
#   2) remove sites lying inside user-specified TF motifs (+/-flank) (-x empty = no exclusion);
#   3) closestBed for nearest-fragment distance -> 3 cols: site_id, fragment length, signed distance.
#
# Usage:
#   bash run_bias_scatter.sh -i <bias_bed_dir> -a <fragment_midP_fragL.bed> -o <out_dir> [options]
#
# Required:
#   -i  bias bed directory (all *.bed in it; each file a 6-col BED for one k-mer)
#   -a  fragment midpoint+length file (6-col output of step 0)
#   -o  output directory
#
# Optional (TF exclusion):
#   -x  TF motif beds to exclude, space-separated and quoted (default: empty = no exclusion)
#         e.g. -x "motif_anchor_bed/CTCF.center.uniq.bed motif_anchor_bed/NFIC.center.uniq.bed"
#   -f  bp added to each side of every motif when excluding (default: 12)
#   note: bias and motif are both stranded; exclusion is same-strand (-s).
#
# Other optional:
#   -m  fragment-length ($5) upper bound (default: no filter)
#   -n  subsample size  (default: 200000; set 0/all/none to disable subsampling, use all sites)
#   -S  sampling seed   (default: 42, reproducible; set none/off for true random)
#   -j  maximum concurrent closestBed jobs (default: 1)
#   -h  show help

set -euo pipefail

BIAS_DIR=""
FRAG_BED=""
OUT_DIR=""
EXCLUDE_LIST=""
FLANK=12
FRAGL_MAX=""
SHUF_N=200000
SEED=42
JOBS=1

usage() { awk 'NR==1{next} /^[^#]/{exit} {sub(/^#[ ]?/,"");print}' "$0"; exit "${1:-0}"; }

while getopts "i:a:o:x:f:m:n:S:j:h" opt; do
    case "$opt" in
        i) BIAS_DIR="$OPTARG" ;;
        a) FRAG_BED="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        x) EXCLUDE_LIST="$OPTARG" ;;
        f) FLANK="$OPTARG" ;;
        m) FRAGL_MAX="$OPTARG" ;;
        n) SHUF_N="$OPTARG" ;;
        S) SEED="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -z "$BIAS_DIR" || -z "$FRAG_BED" || -z "$OUT_DIR" ]] && { echo "ERROR: -i, -a, -o are required" >&2; usage 1; }
[[ -d "$BIAS_DIR" ]] || { echo "ERROR: bias directory not found: $BIAS_DIR" >&2; exit 1; }
[[ -f "$FRAG_BED" ]] || { echo "ERROR: fragment file not found: $FRAG_BED" >&2; exit 1; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -j must be a positive integer" >&2; exit 1; }
[[ "$FLANK" =~ ^[0-9]+$ ]] || { echo "ERROR: -f must be a non-negative integer" >&2; exit 1; }
[[ -z "$FRAGL_MAX" || "$FRAGL_MAX" =~ ^[0-9]+$ ]] || { echo "ERROR: -m must be a non-negative integer" >&2; exit 1; }
case "$SHUF_N" in
    0|all|ALL|none|NONE|no|No|"") ;;
    *) [[ "$SHUF_N" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -n must be a positive integer or 0/all/none" >&2; exit 1; } ;;
esac
command -v closestBed >/dev/null || { echo "ERROR: closestBed not found, please install bedtools" >&2; exit 1; }
command -v intersectBed >/dev/null || { echo "ERROR: intersectBed not found, please install bedtools" >&2; exit 1; }
command -v shuf >/dev/null || { echo "ERROR: GNU shuf not found" >&2; exit 1; }
mkdir -p "$OUT_DIR"

validate_bed6() {
    awk '
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        NF != 6 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 < $2 ||
            !($6=="+" || $6=="-" || $6==".") { bad=1; exit }
        { seen=1 }
        END { exit !(seen && !bad) }
    ' "$1"
}

# fragment-length filter
if [[ -n "$FRAGL_MAX" ]]; then FRAGL_COND="\$5<=maxl"; FRAGL_SHOW="$FRAGL_MAX"; else FRAGL_COND="1"; FRAGL_MAX=0; FRAGL_SHOW="inf"; fi

# Pre-build the exclusion zones: combine excluded TF motifs, expand by +/-flank, sort.
EXCLUDE_BED=""
EXCLUDE_FILES=()
if [[ -n "$EXCLUDE_LIST" ]]; then
    read -r -a EXCLUDE_FILES <<< "$EXCLUDE_LIST"
    EXCLUDE_BED="$(mktemp)"
    trap 'rm -f "$EXCLUDE_BED"' EXIT
    for bed in "${EXCLUDE_FILES[@]}"; do
        [[ -f "$bed" ]] || { echo "ERROR: excluded motif bed not found: $bed" >&2; exit 1; }
        validate_bed6 "$bed" || { echo "ERROR: excluded motif is not a valid BED6: $bed" >&2; exit 1; }
        awk -v fl="$FLANK" '{s=$2-fl; if(s<0)s=0; print $1,s,$3+fl,$4,$5,$6}' OFS='\t' "$bed"
    done | sort -k1,1 -k2,2n -S 80% > "$EXCLUDE_BED"
fi

# whether to subsample: SHUF_N of 0/all/none means no subsampling (use all bias sites)
DO_SHUF=1
case "$SHUF_N" in 0|all|ALL|none|NONE|no|No|"") DO_SHUF=0 ;; esac

# sampling seed: fixed by default (reproducible); -S none/off disables it (true random); irrelevant if not subsampling
SEED_ON=1
case "$SEED" in ""|none|off|NONE|OFF|None) SEED_ON=0 ;; esac
[[ "$DO_SHUF" -eq 0 ]] && SEED_ON=0
if [[ "$SEED_ON" -eq 1 ]] && ! command -v openssl >/dev/null; then
    echo "WARNING: openssl not found, cannot fix the shuf seed; using true random this run" >&2
    SEED_ON=0
fi
# deterministic byte stream from a seed, fed to shuf --random-source (the GNU coreutils way)
get_seeded_random() { openssl enc -aes-256-ctr -pass pass:"$1" -nosalt </dev/zero 2>/dev/null; }

echo "bias dir        : $BIAS_DIR"
echo "fragment file   : $FRAG_BED"
echo "output dir      : $OUT_DIR"
echo "fragLen cap     : $FRAGL_SHOW"
echo "parallel jobs   : $JOBS"
if [[ "$DO_SHUF" -eq 1 ]]; then
    echo "subsample       : yes, $SHUF_N sites"
    echo "sampling seed   : $([[ "$SEED_ON" -eq 1 ]] && echo "$SEED (fixed/reproducible)" || echo "true random")"
else
    echo "subsample       : no (use all bias sites)"
fi
if [[ -n "$EXCLUDE_BED" ]]; then
    echo "exclude TF      : ${#EXCLUDE_FILES[@]} motif(s), flank=+/-${FLANK}bp, same-strand (-s)"
else
    echo "exclude TF      : none (no sites removed)"
fi
echo "-----------------------------------------------"

shopt -s nullglob
biasfiles=("$BIAS_DIR"/*.bed)
[[ ${#biasfiles[@]} -gt 0 ]] || { echo "ERROR: no *.bed in bias directory: $BIAS_DIR" >&2; exit 1; }

for bias in "${biasfiles[@]}"; do
    validate_bed6 "$bias" || { echo "ERROR: bias file is not a valid BED6: $bias" >&2; exit 1; }
done

pids=()
labels=()
failed=0
wait_batch() {
    local i
    for i in "${!pids[@]}"; do
        if ! wait "${pids[$i]}"; then
            echo "ERROR: closestBed job failed for ${labels[$i]}" >&2
            failed=1
        fi
    done
    pids=()
    labels=()
}

for bias in "${biasfiles[@]}"; do
    name="$(basename "$bias" .bed)"
    echo "processing: $name"

    # (optional) subsample -> take single-base strand endpoint (cut site) -> sort
    # when subsampling, derive an independent but reproducible stream per k-mer from "seed_name"
    if [[ "$DO_SHUF" -eq 0 ]]; then
        cat "$bias"
    elif [[ "$SEED_ON" -eq 1 ]]; then
        shuf --random-source=<(get_seeded_random "${SEED}_${name}") -n "$SHUF_N" "$bias"
    else
        shuf -n "$SHUF_N" "$bias"
    fi | \
    awk '{if($6=="+") {print $1,$2,$2+1,$4,$5,$6} else {print $1,$3-1,$3,$4,$5,$6}}' OFS='\t' | \
    sort -k1,1 -k2,2n -S 80% > "${OUT_DIR}/${name}_sites.bed"

    # remove sites inside a TF motif (+/-flank) on the same strand (-x empty = skip)
    if [[ -n "$EXCLUDE_BED" ]]; then
        intersectBed -sorted -s -v -a "${OUT_DIR}/${name}_sites.bed" -b "$EXCLUDE_BED" -wa \
            > "${OUT_DIR}/${name}_sites_kept.bed"
    else
        cp "${OUT_DIR}/${name}_sites.bed" "${OUT_DIR}/${name}_sites_kept.bed"
    fi

    # closestBed: cols 1-6=fragment, 7-12=bias site, 13=distance; sign by strand, keep 0-300
    (
        closestBed -a "$FRAG_BED" -b "${OUT_DIR}/${name}_sites_kept.bed" -d -t first | \
        awk -v maxl="$FRAGL_MAX" \
            "{if(\$13>=0 && \$13<=300 && $FRAGL_COND){if((\$9<=\$3 && \$12==\"+\")||(\$9>\$3 && \$12==\"-\")) print \$10,\$5,\$13; else print \$10,\$5,\$13*(-1)}}" OFS='\t' \
            > "${OUT_DIR}/${name}_fragL_dist.txt"
    ) &
    pids+=("$!")
    labels+=("$name")
    if (( ${#pids[@]} >= JOBS )); then wait_batch; fi
done

wait_batch
(( failed == 0 )) || exit 1
echo "All done. output: ${OUT_DIR}/*_fragL_dist.txt"
