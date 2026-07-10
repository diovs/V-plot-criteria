# V-plot criteria

This repository contains scripts and small supporting tables for the manuscript
**"Two essential criteria for transcription factor footprint discovery using
V-plot"**. It also includes a reusable pipeline for comparing TF motif V-plots
with enzyme-bias controls.

## Repository Layout

| Path | Description |
| --- | --- |
| `Fig1/` | Figure 1 scripts and inputs. |
| `Fig2/` | Figure 2 scripts and input tables. |
| `FigS2/`-`FigS4/` | Supplementary figure scripts and inputs. |
| `FigS5-7/` | Single-panel V-plot plotting script. |
| `FigS8/` | Scrambled CTCF motif control scripts and inputs. |
| `V_plot_pipeline/` | Standalone TF-vs-bias V-plot pipeline. |
| `data_preprocessing.sh` | Fragment preprocessing helper. |
| `loMNase_mapping_pipeline.sh` | loMNase-seq mapping helper. |

Generated figure images and large sequencing files are not included.

## V_plot_pipeline

`V_plot_pipeline/` is the recommended entry point for applying the TF-vs-bias
V-plot criteria to new data. It takes fragment BED files, TF motif coordinates,
and bias k-mers, then reports V-channel width, `log2(Vin / Vout)`, statistical
tests, and cut-offs. Full options are documented in
[`V_plot_pipeline/MANUAL.md`](V_plot_pipeline/MANUAL.md).

### Quick Start

```bash
cd V_plot_pipeline
cp config.example.sh config.sh
# Edit config.sh with paths to the genome FASTA, fragment BED, TF motif BED
# directory, bias k-mer list, and output directory.
bash run_pipeline.sh -c config.sh
```

Supported presets are `loMNase`, `DNase`, `ATAC`, and `custom`. Command-line
options are described in `V_plot_pipeline/MANUAL.md`.

## Dependencies

Core tools:

- `bedtools` >= 2.30.0
- `seqkit` >= 2.0.0
- `python3` >= 3.9.5
- Python packages: `numpy`, `pandas`, `scipy`, `matplotlib`

Some figure scripts also use R packages such as `ggplot2`, `ggseqlogo`, and
`data.table`.

## Data Sources

K562 DNase-seq (`ENCSR000EOT`) and ATAC-seq (`ENCSR483RKN`, `ENCSR868FGK`) were
downloaded from ENCODE. loMNase-seq and CTCF/MAZ native ChIP-seq data are
available from GSA-Human under `HRA005744` and `HRA017804`.

Large sequencing files and genome FASTA files should be downloaded separately.

## Citation

If you use this code or pipeline, please cite the accompanying manuscript:

Zhang Q and Xu C. **Two essential criteria for transcription factor footprint
discovery using V-plot**.
