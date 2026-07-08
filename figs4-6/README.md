# FigS4-S6 V-plot panel script

This folder contains one reusable script for drawing a single V-plot panel used
in the supplementary Fig. S4-S6 montages:

- Fig. S4: loMNase
- Fig. S5: DNase
- Fig. S6: ATAC

The large `*_fragL_dist.txt` input tables are not included in this repository.
They can be regenerated with the `V_plot_pipeline` distance-table steps or
provided from an existing analysis directory.

## Input

`plot_single_vplot_panel.R` expects the three-column V-plot distance table used
by the pipeline:

```text
site_id    fragment_length    signed_distance
```

This is the same format written by `run_TF_scatter.sh` and
`run_bias_scatter.sh`.

## Examples

TF panel with fitted V-channel guides:

```bash
Rscript plot_single_vplot_panel.R \
  --input CTCF_loMNase_fragL_dist.txt \
  --assay loMNase \
  --title CTCF \
  --apex-x -1.69 \
  --apex-y 31.13 \
  --channel true \
  --output-prefix S4_CTCF
```

Bias panel without V-channel guides:

```bash
Rscript plot_single_vplot_panel.R \
  --input MNase_Bias_AAG_loMNase_fragL_dist.txt \
  --assay loMNase \
  --title AAG \
  --channel false \
  --output-prefix S4_AAG
```

The script writes PNG and PDF files by default. Use `--formats png` to write
only PNG output.
