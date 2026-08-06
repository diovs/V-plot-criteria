#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="$(cd "$HERE/.." && pwd)"
OUT_ROOT="${1:-$HERE/results}"
THREADS="${THREADS:-4}"

cd "$HERE"

for assay in loMNase DNase ATAC; do
    echo "===== Running ${assay} example ====="
    bash "$PIPELINE/run_pipeline.sh" \
        --genome "genome/hg38_chr22_test.fa" \
        --fragment "fragments/${assay}_chr22_test.bed" \
        --tf-dir "motifs" \
        --bias-kmers "bias_kmers/${assay}_top5.txt" \
        --exclude-tf "$(printf '%s ' motifs/*.bed)" \
        --mode "$assay" \
        --assay "${assay}_example" \
        --out "$OUT_ROOT/$assay" \
        --threads "$THREADS" \
        --shuf-n 300 \
        --shuf-seed 42 \
        -k 0
done

echo "===== All examples completed: $OUT_ROOT ====="
echo "Verify with: python3 verify_expected.py --actual '$OUT_ROOT'"
