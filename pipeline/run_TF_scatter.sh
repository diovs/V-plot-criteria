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
#   -h  show help

set -euo pipefail

MOTIF_DIR=""
FRAG_BED=""
OUT_DIR=""
FRAGL_MAX=""

usage() { awk 'NR==1{next} /^[^#]/{exit} {sub(/^#[ ]?/,"");print}' "$0"; exit "${1:-0}"; }

while getopts "i:a:o:m:h" opt; do
    case "$opt" in
        i) MOTIF_DIR="$OPTARG" ;;
        a) FRAG_BED="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        m) FRAGL_MAX="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -z "$MOTIF_DIR" || -z "$FRAG_BED" || -z "$OUT_DIR" ]] && { echo "ERROR: -i, -a, -o are required" >&2; usage 1; }
[[ -d "$MOTIF_DIR" ]] || { echo "ERROR: motif directory not found: $MOTIF_DIR" >&2; exit 1; }
[[ -f "$FRAG_BED" ]]  || { echo "ERROR: fragment file not found: $FRAG_BED" >&2; exit 1; }
mkdir -p "$OUT_DIR"

# check whether a file is a valid 6-col motif bed (first non-empty/non-comment line)
validate_motif() {
    local f="$1" line nf
    line=$(grep -v '^#' "$f" | grep -m1 -v '^[[:space:]]*$' || true)
    [[ -z "$line" ]] && return 1
    nf=$(awk '{print NF; exit}' <<<"$line")
    [[ "$nf" -eq 6 ]] || return 1
    awk '{exit ($6=="+"||$6=="-"||$6==".")?0:1}' <<<"$line"
}

# fragment-length filter: always true when -m is not given
if [[ -n "$FRAGL_MAX" ]]; then FRAGL_COND="\$5<=maxl"; FRAGL_SHOW="$FRAGL_MAX"; else FRAGL_COND="1"; FRAGL_MAX=0; FRAGL_SHOW="inf"; fi

echo "motif dir       : $MOTIF_DIR"
echo "fragment file   : $FRAG_BED"
echo "output dir      : $OUT_DIR"
echo "fragLen cap     : $FRAGL_SHOW"
echo "-----------------------------------------------"

shopt -s nullglob
files=("$MOTIF_DIR"/*)
[[ ${#files[@]} -gt 0 ]] || { echo "ERROR: motif directory is empty: $MOTIF_DIR" >&2; exit 1; }

# validate everything first; abort on the first bad file (so no TF is silently dropped)
for f in "${files[@]}"; do
    [[ -f "$f" ]] || continue
    if ! validate_motif "$f"; then
        echo "ERROR: not a valid 6-col motif bed (chr start end name score strand): $f" >&2
        exit 1
    fi
done

for motif in "${files[@]}"; do
    [[ -f "$motif" ]] || continue
    name="$(basename "$motif")"; name="${name%%.*}"      # CTCF.center.uniq.bed -> CTCF
    echo "processing: $name"

    # closestBed -a fragments -b motif -d : cols 1-6=fragment, 7-12=motif, 13=distance
    # sign the distance by motif strand; keep only distance 0-300 (and optional $5<=maxl)
    nohup closestBed -a "$FRAG_BED" -b "$motif" -d -t first | \
    awk -v maxl="$FRAGL_MAX" \
        "{if(\$13>=0 && \$13<=300 && $FRAGL_COND){if((\$9<=\$3 && \$12==\"+\")||(\$9>\$3 && \$12==\"-\")) print \$10,\$5,\$13; else print \$10,\$5,\$13*(-1)}}" OFS='\t' \
        > "${OUT_DIR}/${name}_fragL_dist.txt" &
done

wait
echo "All done. output: ${OUT_DIR}/*_fragL_dist.txt"
