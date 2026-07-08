# FigS3 panel scripts

This folder contains the scripts and input tables for panels C and D of the
revised supplementary V-plot model overview figure.

## Files

- `FigS3C_plot.py`: plots the loMNase apex-position likelihood profile for CTCF
  and the AAG bias control.
- `FigS3C_real_loMNase_apex_likelihood_profiles.csv`: input likelihood scan
  used by `FigS3C_plot.py`.
- `FigS3D_plot.py`: plots the two-feature cut-off derivation for loMNase, DNase
  and ATAC.
- `FigS3D_cutoff_thresholds.csv`: mid-gap cut-offs and margins for each assay
  and feature.
- `FigS3D_plot_data_best_loMNase.csv`, `FigS3D_plot_data_best_DNase.csv`,
  `FigS3D_plot_data_best_ATAC.csv`: TF and enzyme-bias values used by
  `FigS3D_plot.py`.

## Run

```bash
python FigS3C_plot.py
python FigS3D_plot.py
```

Each script writes PNG and PDF files in the same directory by default. The
AI-polished composite figure is not included here; these scripts reproduce the
same plotting logic using the archived panel data.
