# Revised supplementary figure notes

## Outputs

- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview.png`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview.pdf`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview.svg`


Individual panel files:

- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_A.png`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_A.pdf`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_A.svg`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_C.png`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_C.pdf`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_C.svg`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_D.png`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_D.pdf`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_overview_panel_D.svg`
- `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure\supp_methodB_revised_figure_notes.md`

## Input roots

- input root: `method_B\figures\supp_methodB_revised_overview\V_plot_pipeline`
- method root: `method_B\figures\supp_methodB_revised_overview`
- result directory: `method_B\figures\supp_methodB_revised_overview\methodB_supp_revised_figure`

## Reuse and inspection summary

- The revised script redraws the conceptual V-plot model panel, real loMNase apex-likelihood scan, workflow panel, and cut-off derivation panel.
- The previous two-dimensional scatter panel was removed from the overview.
- The figure avoids the internal project label and uses article-facing wording: V-plot model, V-channel width, inside-V enrichment, LRT, and calibrated effect-size cut-offs.

## Data sources

- cut-off derivation table: `method_B\figures\supp_methodB_revised_overview\fig_cutoff_derivation_thresholds.csv`
- separation/validation summary: `method_B\figures\supp_methodB_revised_overview\v_apex_scores\methodB_scatter_best_assay_specific\methodB_separation_stats_best_common.csv`

Cut-off derivation table:

| assay   | feature             | cutoff             | margin             | max_bias           | min_TF             | unit |
| ------- | ------------------- | ------------------ | ------------------ | ------------------ | ------------------ | ---- |
| loMNase | V-channel width     | 6.708921215961457  | 6.708471618986781  | 0.0004495969746763 | 13.417392834948238 | bp   |
| loMNase | Inside-V enrichment | 0.5067117145116604 | 0.4326199023578587 | 0.0740918121538016 | 0.939331616869519  | log2 |
| DNase   | V-channel width     | 9.154892066855576  | 9.154595157750196  | 0.0002969091053805 | 18.30948722460577  | bp   |
| DNase   | Inside-V enrichment | 0.338756749077624  | 0.2038491964576881 | 0.1349075526199359 | 0.5426059455353122 | log2 |
| ATAC    | V-channel width     | 19.998395142799787 | 19.148438893082343 | 0.849956249717442  | 39.14683403588213  | bp   |
| ATAC    | Inside-V enrichment | 0.3732815147375332 | 0.1359947540601656 | 0.2372867606773675 | 0.5092762687976988 | log2 |

## Panel provenance

- Panel A combines a V-plot model schematic with a real loMNase likelihood scan for CTCF and one sequence-bias control. The V-plot schematic uses equal data-unit aspect ratio.
- Panel C is a simplified workflow adapted from the V-plot calibration overview.
- Panel D redraws the cut-off derivation shown in `fig_cutoff_derivation.png` and adds an enlarged cut-off/validation table.

## Notes and unresolved items

- The representative TF/bias V-plot panel was removed because those examples are already shown in another supplementary figure.
- The script does not modify fitting parameters or thresholding parameters; it only redraws the figure.

## Figure legend

Supplementary Figure X. V-plot model-based classification of TF-footprint-like patterns using calibrated effect sizes.
(A) V-plot model schematic and MLE-based apex estimation. The fitted apex `(a, b)` defines the V-channel width `b`; the loMNase likelihood scan shows how candidate apex positions are supported by real CTCF and sequence-bias data.
(C) Simplified reproducible analysis workflow from fragments and anchors to fitted V-plot features and calibrated cut-offs.
(D) Derivation of the cut-offs by the mid-gap rule. For each assay and feature, the cut-off is the midpoint between the highest enzyme-bias value and the lowest TF value; the compact table summarizes the resulting cut-offs and validation performance.
