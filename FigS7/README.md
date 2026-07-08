# FigS7 README

Fig. S7 shows assay-specific TF motif logos above representative V-plot panels.
This folder intentionally contains only this README, because the panel does not
introduce a separate analysis workflow.

The figure combines two reusable operations that are already documented
elsewhere in this repository:

- Motif logo rendering: see `FigS8/plot_ctcf_scramble_logos.py` for the logo
  drawing logic used for MEME/PWM-style motifs.
- Single V-plot panel rendering: see
  `FigS4-6/plot_single_vplot_panel.R` for drawing a V-plot panel from a
  three-column fragment-distance table.

No generated figure images are stored here, and no duplicate plotting script is
included for Fig. S7. The final manuscript panel was assembled from component
logo and V-plot panels, with minor manual layout adjustment.

## Panel Content

| Assay | TF motifs shown |
| --- | --- |
| loMNase | NFIC, TFE3 |
| DNase | NFE2, REST |
| ATAC | CEBPG, CREB1 |

## Reproduction Notes

To regenerate the component panels:

1. Draw each TF motif logo using the same MEME/PWM logo-rendering approach as
   `FigS8/plot_ctcf_scramble_logos.py`.
2. Draw each corresponding V-plot using
   `FigS4-6/plot_single_vplot_panel.R` and the assay-specific
   fragment-distance table for that TF.
3. Assemble the six logo-over-V-plot pairs into the manuscript layout.

The large fragment-distance tables and generated panel images are not included
in this repository.
