# Method B

Method B fits a generative V-plot model and classifies TF-footprint-like patterns using two calibrated effect-size features:

- V-channel width
- inside-V enrichment, `log2(V-in / V-out)`

The likelihood-ratio test is used to detect V-like structure, while TF-vs-enzyme-bias classification uses the calibrated two-feature rule.

## Directory map

- `V_plot_pipeline/`  
  Main portable pipeline. Start here for rerunning Method B from fragment BED input.

- `figures/discriminant_scatter/`  
  Recreates the assay-specific TF-vs-bias scatter plots.

- `figures/cutoff_derivation/`  
  Recreates the cut-off derivation figure from the small `data/plot_data_best_*.csv` tables.

- `figures/supp_methodB_revised_overview/`  
  Recreates the final Method B supplementary overview figure.

- `figures/weblogo_vplot_grid/`  
  Recreates the selected weblogo plus V-plot montage from six curated panel PNGs.

- `source_tables/`  
  Small summary tables shared by the figure scripts.

## Notes

This is a curated GitHub-facing copy. The local backup at `new_script/method_B` still contains the full exploratory history and large intermediates.

