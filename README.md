# V-plot criteria

This repository contains the scripts, supporting tables, and reusable pipeline
for the manuscript **"Two essential criteria for transcription factor footprint
discovery using V-plot"**.

The project addresses a practical problem in V-plot interpretation: enzyme
sequence-cleavage bias can create V-shaped patterns in loMNase-seq, DNase-seq,
and ATAC-seq data, but these bias-derived V-shapes are not equivalent to
transcription-factor (TF) footprints. The analyses here use two graphical
criteria to distinguish TF-footprint-like V-plots from enzyme-bias controls:

- **V-channel width**, the protected region indicated by the empty channel at
  the V apex.
- **Inside-V enrichment**, reported as `log2(Vin / Vout)`, the fragment-density
  enrichment inside the fitted V region relative to the outside region.

## Repository Layout

| Path | Description |
| --- | --- |
| `Fig1/` | Scripts and input files for Figure 1 analyses. |
| `Fig2/` | Scripts and input tables for the revised Figure 2 panels, including the TF-vs-bias separation plots. |
| `FigS1/`-`FigS3/` | Scripts and input files for supplementary figures. |
| `FigS4-6/` | Single-panel V-plot plotting script for the revised Fig. S4-S6 montages. |
| `FigS8/` | Scripts and small input files for the scrambled CTCF motif control. |
| `V_plot_pipeline/` | Standalone pipeline for testing TF motif V-plots against enzyme-bias controls. |
| `data_preprocessing.sh` | Fragment preprocessing helper script. |
| `loMNase_mapping_pipeline.sh` | loMNase-seq read mapping and fragment-generation helper script. |

Generated figure images are not the focus of this repository. The included
scripts and tables are intended to regenerate the corresponding panels or to
serve as documented examples of the analysis workflow.

## V_plot_pipeline

`V_plot_pipeline/` is the recommended entry point for users who want to apply
the criteria to their own data. It starts from a fragment BED file and TF motif
coordinates, builds V-plot fragment-distance tables, fits the Vplot model, and
compares TF motifs with matched enzyme-bias k-mer controls.

In brief, the pipeline:

1. Converts raw fragment BED intervals to midpoint-plus-length BED records.
2. Locates user-specified enzyme-bias k-mers in a reference genome.
3. Builds V-plot distance tables for TF motifs and bias sites.
4. Estimates V-channel width and inside-V enrichment for each V-plot.
5. Tests TF-vs-bias separation with Wilcoxon rank-sum tests.
6. Applies two-dimensional cut-offs for V-channel width and inside-V enrichment;
   when TF and bias classes are separated, these correspond to the
   Youden-optimal thresholds.

A motif set is considered TF-footprint-like only when both graphical features
exceed their cut-offs. See
[`V_plot_pipeline/MANUAL.md`](V_plot_pipeline/MANUAL.md) for the complete input
formats, options, intermediate files, and output descriptions.

### Quick Start

```bash
cd V_plot_pipeline
cp config.example.sh config.sh
# Edit config.sh with paths to the genome FASTA, fragment BED, TF motif BED
# directory, bias k-mer list, and output directory.
bash run_pipeline.sh -c config.sh
```

The same options can also be provided directly:

```bash
bash run_pipeline.sh \
  --genome hg38.fa \
  --fragment fragments.bed \
  --tf-dir motif_beds/ \
  --bias-kmers bias_kmers.txt \
  --exclude-tf "motif_beds/CTCF.bed" \
  --mode ATAC \
  --out results_ATAC
```

Supported presets are `loMNase`, `DNase`, `ATAC`, and `custom`.

## Dependencies

The pipeline was tested with:

- `bedtools` >= 2.30.0
- `seqkit` >= 2.0.0
- `python3` >= 3.9.5
- Python packages: `numpy`, `pandas`, `scipy`, `matplotlib`

Some figure-generation scripts additionally use R packages such as `ggplot2`,
`ggseqlogo`, and `data.table`. See the individual scripts for figure-specific
requirements and local path variables.

## Data Sources

The manuscript uses loMNase-seq and CTCF/MAZ native ChIP-seq data available from
GSA-Human under accessions `HRA005744` and `HRA017804`. DNase-seq and ATAC-seq
examples use ENCODE datasets including `ENCSR000EOT`, `ENCSR483RKN`, and
`ENCSR868FGK`.

Large sequencing files, genome FASTA files, and large derived alignment files
are not stored in this repository. Users should download the required public
datasets and provide local paths in the scripts or in `V_plot_pipeline/config.sh`.

## Citation

If you use this code or pipeline, please cite the accompanying manuscript:

Zhang Q and Xu C. **Two essential criteria for transcription factor footprint
discovery using V-plot**.
