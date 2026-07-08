# FigS8 scripts

This folder contains the scripts and small input files used to reproduce the
logic of Fig. S8.

## Scripts

- `build_scrambled_ctcf_pwm.py`
  Builds a CTCF position-probability matrix from strand-aware CTCF motif
  anchors and creates a scrambled control motif by permuting the PWM columns.
  The script writes the original motif, scrambled motif, column permutation and
  base-count tables.

- `plot_ctcf_scramble_logos.py`
  Draws the CTCF logo and three scrambled CTCF logos from MEME-format motif
  files in `data/`.

- `plot_ctcf_scramble_scatter.py`
  Draws the loMNase, DNase and ATAC V-plot feature scatter panels comparing
  real CTCF with scrambled CTCF motifs.

## Small Input Files

- `data/ctcf_original.meme`
- `data/ctcf_scramble_1.meme`
- `data/ctcf_scramble_2.meme`
- `data/ctcf_scramble_3.meme`
- `data/ctcf_scramble_points.csv`
- `data/ctcf_scramble_stats.csv`

Large FIMO outputs, genome-wide BED files and V-plot distance tables are not
included. They can be regenerated from the CTCF motif coordinates and the genome
FASTA when needed.

## Run

Build a new scrambled CTCF motif:

```bash
python build_scrambled_ctcf_pwm.py \
  --bed Hs_CTCF.motif \
  --genome hg38.fa \
  --flank 9 \
  --seed 0 \
  --out-prefix CTCF_f9_s0
```

Plot the logos and scatter panels:

```bash
python plot_ctcf_scramble_logos.py
python plot_ctcf_scramble_scatter.py
```
