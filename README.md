# V-plot criteria GitHub-ready package

This folder is a curated, lightweight copy of the Method B V-plot criteria files.

The original `method_B` directory was left unchanged as a backup. Large intermediate files, old trials, local caches, raw distance distributions, genome-wide motif/bias BED files, and server-only scratch outputs were not copied.

## Contents

- `method_B/V_plot_pipeline/`  
  Portable Method B V-plot pipeline. This is the main code path for running the model from fragment BED files to TF-vs-bias cut-offs and scatter plots.

- `method_B/figures/discriminant_scatter/`  
  Script, small input tables, and final scatter outputs for the width x inside-V enrichment separation plots.

- `method_B/figures/cutoff_derivation/`  
  Self-contained cut-off derivation figure using the mid-gap rule.

- `method_B/figures/supp_methodB_revised_overview/`  
  Final revised Method B supplementary overview figure and the minimal tables needed to redraw it.

- `method_B/figures/weblogo_vplot_grid/`  
  Final weblogo-over-V-plot panel plus the selected input PNGs needed to regenerate it.

- `method_B/source_tables/`  
  Shared small summary tables copied for inspection and reuse.

## What was intentionally excluded

- `TF_vplot/*_fragL_dist.txt`
- `sequence_bias/*_fragL_dist.txt`
- genome FASTA / FAI files
- genome-wide motif and bias BED directories
- `__pycache__`, run logs, and old exploratory scripts
- full candidate-panel archives that are too large for normal GitHub use

If full raw/intermediate data are needed, regenerate them with `method_B/V_plot_pipeline/run_pipeline.sh` from the public inputs described in the manuscript or from the user's own data.
