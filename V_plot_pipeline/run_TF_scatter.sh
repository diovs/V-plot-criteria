#!/bin/bash
# Step 2 (TF): distance distribution of fragment midpoints around TF motif sites (raw V-plot data)
#
# For each motif bed in the directory, use closestBed to find each fragment's nearest
# motif, and output 3 cols: site_id, fragment length, signed distance -> input to
# fit_vplot_apex.py.
#
# Usage:
#   bash run_TF_scatter.sh -i <motif_dir> -a <fragment_midP_fragL.bed> -o <out_dir> [options]
#
# Required:
#   -i  directory of TF motif beds (every file in it is processed; extension ignored)
#       each file must be a 6-col motif bed: chr start end name score strand
#       (strand +/-/.) ; any file not matching this format aborts the run.
#   -a  fragment midpoint+length file (6-col output of step 0: chr midP midP+1 name fragL strand)
#   -o  output directory
#
# Optional:
#   -m  fragment-length ($5) upper bound (default: no filter; keep only $5 <= value)
#   -j  maximum concurrent closestBed jobs (default: 1)
#   -h  show help

set -euo pipefail

MOTIF_DIR=""
FRAG_BED=""
OUT_DIR=""
FRAGL_MAX=""
JOBS=1

usage() { awk 'NR==1{next} /^[^#]/{exit} {sub(/^#[ ]?/,"");print}' "$0"; exit "${1:-0}"; }

while getopts "i:a:o:m:j:h" opt; do
    case "$opt" in
        i) MOTIF_DIR="$OPTARG" ;;
        a) FRAG_BED="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        m) FRAGL_MAX="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -z "$MOTIF_DIR" || -z "$FRAG_BED" || -z "$OUT_DIR" ]] && { echo "ERROR: -i, -a, -o are required" >&2; usage 1; }
[[ -d "$MOTIF_DIR" ]] || { echo "ERROR: motif directory not found: $MOTIF_DIR" >&2; exit 1; }
[[ -f "$FRAG_BED" ]]  || { echo "ERROR: fragment file not found: $FRAG_BED" >&2; exit 1; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -j must be a positive integer" >&2; exit 1; }
[[ -z "$FRAGL_MAX" || "$FRAGL_MAX" =~ ^[0-9]+$ ]] || { echo "ERROR: -m must be a non-negative integer" >&2; exit 1; }
command -v closestBed >/dev/null || { echo "ERROR: closestBed not found, please install bedtools" >&2; exit 1; }
mkdir -p "$OUT_DIR"

# Check every data row in a motif file before launching any expensive jobs.
validate_motif() {
    awk '
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        NF != 6 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 < $2 ||
            !($6=="+" || $6=="-" || $6==".") { bad=1; exit }
        { seen=1 }
        END { exit !(seen && !bad) }
    ' "$1"
}

# fragment-length filter: always true when -m is not given
if [[ -n "$FRAGL_MAX" ]]; then FRAGL_COND="\$5<=maxl"; FRAGL_SHOW="$FRAGL_MAX"; else FRAGL_COND="1"; FRAGL_MAX=0; FRAGL_SHOW="inf"; fi

echo "motif dir       : $MOTIF_DIR"
echo "fragment file   : $FRAG_BED"
echo "output dir      : $OUT_DIR"
echo "fragLen cap     : $FRAGL_SHOW"
echo "parallel jobs   : $JOBS"
echo "-----------------------------------------------"

shopt -s nullglob
files=("$MOTIF_DIR"/*)
[[ ${#files[@]} -gt 0 ]] || { echo "ERROR: motif directory is empty: $MOTIF_DIR" >&2; exit 1; }

# Validate everything first and detect output-name collisions before starting.
declare -A output_names=()
for f in "${files[@]}"; do
    [[ -f "$f" ]] || continue
    if ! validate_motif "$f"; then
        echo "ERROR: not a valid 6-col motif bed (chr start end name score strand): $f" >&2
        exit 1
    fi
    name="$(basename "$f")"; name="${name%%.*}"
    [[ -z "${output_names[$name]:-}" ]] || {
        echo "ERROR: motif filenames map to the same output name '$name'" >&2
        exit 1
    }
    output_names[$name]=1
done
(( ${#output_names[@]} > 0 )) || { echo "ERROR: motif directory contains no regular files" >&2; exit 1; }

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

for motif in "${files[@]}"; do
    [[ -f "$motif" ]] || continue
    name="$(basename "$motif")"; name="${name%%.*}"
    echo "processing: $name"

    # closestBed -a fragments -b motif -d : cols 1-6=fragment, 7-12=motif, 13=distance
    # sign the distance by motif strand; keep only distance 0-300 (and optional $5<=maxl)
    (
        closestBed -a "$FRAG_BED" -b "$motif" -d -t first | \
        awk -v maxl="$FRAGL_MAX" \
            "{if(\$13>=0 && \$13<=300 && $FRAGL_COND){if((\$9<=\$3 && \$12!=\"-\")||(\$9>\$3 && \$12==\"-\")) print \$10,\$5,\$13; else print \$10,\$5,\$13*(-1)}}" OFS='\t' \
            > "${OUT_DIR}/${name}_fragL_dist.txt"
    ) &
    pids+=("$!")
    labels+=("$name")
    if (( ${#pids[@]} >= JOBS )); then wait_batch; fi
done

wait_batch
(( failed == 0 )) || exit 1
echo "All done. output: ${OUT_DIR}/*_fragL_dist.txt"
