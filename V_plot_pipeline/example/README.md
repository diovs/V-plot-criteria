# Bundled three-assay example

This directory contains small, real-data examples for a complete test run of
`V_plot_pipeline` with K562 loMNase-seq, DNase-seq, and ATAC-seq. The examples
show the required input formats and provide reference outputs for checking a
local installation.

## Quick start

From this directory:

```bash
bash run_examples.sh
python3 verify_expected.py
```

Set the number of concurrent jobs if needed:

```bash
THREADS=8 bash run_examples.sh
```

The runs sample 300 sites per bias sequence (`--shuf-n 300`) and use a fixed input
region, so they are deterministic. Generated files are written under `results/`.

## Input formats

All text inputs are tab-delimited/plain text and have no header unless noted
otherwise.

- Fragment BED3: `chrom`, `start`, `end`; one mapped fragment per row, with
  0-based, half-open genomic coordinates.
- Motif BED6: `chrom`, `start`, `end`, `name`, `score`, `strand`; the interval is
  the one-base motif centre used as the V-plot anchor.
- Bias k-mer list: one DNA sequence per line, without names or extra columns.
- Genome FASTA: the chromosome names must match column 1 of both BED inputs.

For example (`\t` denotes a literal tab character):

```text
# fragments/ATAC_chr22_test.bed (BED3)
chr22_test\t1\t51

# motifs/CTCF.bed (BED6)
chr22_test\t13596\t13597\t67931\t67931\t-

# bias_kmers/ATAC_top5.txt
GCC
```

## Test data

All coordinates were cropped from hg38 `chr22:20,000,000-22,000,000` (0-based,
half-open) and shifted to a 2-Mb test contig named `chr22_test`.

| Input | Rows/sites | Format |
|---|---:|---|
| `fragments/loMNase_chr22_test.bed` | 591,131 | BED3 fragments |
| `fragments/DNase_chr22_test.bed` | 586,134 | BED3 fragments |
| `fragments/ATAC_chr22_test.bed` | 604,111 | BED3 fragments |
| `motifs/CTCF.bed` | 46 | BED6 motif anchors |
| `motifs/CREB1.bed` | 22 | BED6 motif anchors |
| `motifs/CREM.bed` | 15 | BED6 motif anchors |
| `motifs/NFIC.bed` | 23 | BED6 motif anchors |
| `motifs/ZBTB33.bed` | 5 | BED6 motif anchors |
| `motifs/ZNF281.bed` | 19 | BED6 motif anchors |

CTCF is the primary illustrative TF. CREB1, CREM, NFIC, ZBTB33, and ZNF281 are
statistical companions because the final rank-sum gate compares groups of TF
and bias V-plots. They were retained only after confirming positive fitted
V-channel widths in all three test assays. The additional TFs should not be
interpreted as independent biological examples from this small interval.

The five assay-matched cleavage-bias sequences are:

| Assay | Bias k-mers |
|---|---|
| loMNase | AGA, AGT, TGG, AAG, ACA |
| DNase | TCC, TGT, TGG, TGA, TGC |
| ATAC | GCT, GTC, GCC, GGC, GTG |

## Data provenance

- loMNase-seq fragments are a K562 pooled dataset from HRA005744 and HRA017804.
- DNase-seq fragments are derived from K562 ENCODE experiment ENCSR000EOT.
- ATAC-seq fragments combine K562 ENCODE experiments ENCSR483RKN and
  ENCSR868FGK.
- The reference sequence is hg38/GRCh38.

The loMNase and DNase source tables stored fragment midpoint and fragment
length. BED3 intervals were reconstructed as
`start = midpoint - floor(length/2)` and `end = start + length`. ATAC fragments
were cropped directly from the merged BED3 file. Fragments crossing a test
region boundary were excluded.

## Expected results

`expected/<assay>/` contains the final output files from a validated run. The
verification script checks motif identities, fit status, sample counts, apex
parameters, enrichment, rank-sum results, cut-offs, and the presence of the
scatter plots. Small floating-point differences across supported SciPy versions
are allowed; changes large enough to alter the interpretation fail validation.

The validated run passes both rank-sum gates for every assay and produces these
test-specific two-dimensional cut-offs:

| Assay | Width cut-off (bp) | E cut-off (`log2(V-in/V-out)`) |
|---|---:|---:|
| loMNase | 9.126 | 0.293 |
| DNase | 6.726 | 0.649 |
| ATAC | 9.468 | 0.617 |

These files are intended as an installation and input-format test. They are a
small genomic subset and should not replace whole-genome data for estimating
study-level biological cut-offs.
