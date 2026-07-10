#!/bin/bash
# Step 1: user-defined bias sequences (k-mer list) -> genome-wide coordinate beds
#
# For each sequence in the k-mer list, use seqkit locate to find all genomic
# occurrences, writing one bed per k-mer: <out_dir>/<prefix>_<kmer>.bed (6-col BED,
# stranded). This step needs seqkit and a genome fasta (provided via -g).
#
# Usage:
#   bash 1_locate_bias.sh -k <kmer_list.txt> -g <genome.fa> -o <out_dir> [options]
#
# Required:
#   -k  bias k-mer list, one sequence per line (e.g. ATAC_Bias_3bp_top10.txt)
#   -g  genome fasta (e.g. hg38XY.fa)
#   -o  output directory
#
# Optional:
#   -p  output filename prefix  (default: Bias)
#   -j  seqkit threads          (default: 40)
#   -h  show help

set -euo pipefail

KMER=""
GENOME=""
OUT_DIR=""
PREFIX="Bias"
THREADS=40
TMP_BED=""

trap '[[ -z "$TMP_BED" ]] || rm -f "$TMP_BED"' EXIT

usage() { awk 'NR==1{next} /^[^#]/{exit} {sub(/^#[ ]?/,"");print}' "$0"; exit "${1:-0}"; }

while getopts "k:g:o:p:j:h" opt; do
    case "$opt" in
        k) KMER="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -z "$KMER" || -z "$GENOME" || -z "$OUT_DIR" ]] && { echo "ERROR: -k, -g, -o are required" >&2; usage 1; }
[[ -f "$KMER" ]]   || { echo "ERROR: k-mer list not found: $KMER" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: genome fasta not found: $GENOME" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -j must be a positive integer" >&2; exit 1; }
[[ "$PREFIX" != */* && "$PREFIX" != *\\* ]] || { echo "ERROR: -p must not contain path separators" >&2; exit 1; }
command -v seqkit >/dev/null || { echo "ERROR: seqkit not found, please install it" >&2; exit 1; }
mkdir -p "$OUT_DIR"

echo "k-mer list  : $KMER"
echo "genome      : $GENOME"
echo "output dir  : $OUT_DIR"
echo "prefix/jobs : $PREFIX / $THREADS"
echo "-----------------------------------------------"

# Locate each k-mer; seqkit locate --bed gives 6-col BED
# (chr,start,end,kmer,score,strand). Blank lines and comments are ignored.
declare -A seen=()
count=0
while IFS= read -r id || [[ -n "$id" ]]; do
    id="${id//$'\r'/}"
    id="${id//[[:space:]]/}"
    [[ -z "$id" || "$id" == \#* ]] && continue
    id="${id^^}"
    [[ "$id" =~ ^[ACGT]+$ ]] || {
        echo "ERROR: invalid k-mer '$id' (only A/C/G/T sequences are supported)" >&2
        exit 1
    }
    [[ -z "${seen[$id]:-}" ]] || {
        echo "ERROR: duplicate k-mer in list: $id" >&2
        exit 1
    }
    seen[$id]=1
    ((count+=1))
    echo "locating: $id"
    out_bed="${OUT_DIR}/${PREFIX}_${id}.bed"
    TMP_BED="${out_bed}.tmp.$$"
    if ! seqkit locate -j "$THREADS" --bed -p "$id" "$GENOME" > "$TMP_BED"; then
        rm -f "$TMP_BED"
        TMP_BED=""
        echo "ERROR: seqkit locate failed for $id" >&2
        exit 1
    fi
    mv -f "$TMP_BED" "$out_bed"
    TMP_BED=""
done < "$KMER"

(( count > 0 )) || { echo "ERROR: k-mer list contains no sequences" >&2; exit 1; }

echo "done: $count bias sequence(s) written to $OUT_DIR"
