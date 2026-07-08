#!/bin/bash
# Step 0: raw fragment bed (3 cols) -> midpoint + fragment-length bed (6 cols, sorted)
#
# closestBed needs not the raw 3-col fragments but a 6-col file:
#   chr  midP  midP+1  name  fragL  strand
# where midP = int((start+end)/2) and fragL = end-start.
# This step converts a 3-col file (e.g. fragment_bed/*.bed) to that format and sorts it.
#
# Usage:
#   bash 0_prepare_fragments.sh -i <fragment_bed(3col)> -o <out_midP_fragL.bed> [options]
#
# Required:
#   -i  raw fragment bed, 3 cols (chr, start, end)
#   -o  output 6-col midP+fragL bed
#
# Optional:
#   -m  fragment-length upper bound (default: no filter; keep only fragments with end-start <= value)
#   -M  fragment-length lower bound (default: 0)
#   -h  show help

set -euo pipefail

IN=""
OUT=""
LMAX=""
LMIN=0

usage() { awk 'NR==1{next} /^[^#]/{exit} {sub(/^#[ ]?/,"");print}' "$0"; exit "${1:-0}"; }

while getopts "i:o:m:M:h" opt; do
    case "$opt" in
        i) IN="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        m) LMAX="$OPTARG" ;;
        M) LMIN="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -z "$IN" || -z "$OUT" ]] && { echo "ERROR: -i and -o are required" >&2; usage 1; }
[[ -f "$IN" ]] || { echo "ERROR: input fragment file not found: $IN" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"      # create output dir if it does not exist

# fragment-length upper bound: always true when -m is not given
if [[ -n "$LMAX" ]]; then UP="\$3-\$2<=lmax"; LMAX_SHOW="$LMAX"; else UP="1"; LMAX=0; LMAX_SHOW="inf"; fi

echo "input fragments : $IN"
echo "output file     : $OUT"
echo "fragment length : [${LMIN}, ${LMAX_SHOW}]"
echo "-----------------------------------------------"

# 6 cols: chr  mid  mid+1  rownum  fragLen  rownum   (mid=int((start+end)/2), fragL=end-start)
awk -v lmin="$LMIN" -v lmax="$LMAX" \
    "{ fl=\$3-\$2; if(fl>=lmin && $UP){ mid=int((\$2+\$3)/2); print \$1, mid, mid+1, NR, fl, NR } }" \
    OFS='\t' "$IN" | \
sort -k1,1 -k2,2n -S 80% > "$OUT"

echo "done: $(wc -l < "$OUT") fragment midpoints written to $OUT"
