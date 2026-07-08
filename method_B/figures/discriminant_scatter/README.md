# Discriminant scatter plots

Run from this directory:

```bash
python plot_methodB_best_common_TFs.py
```

The script reads the local TF apex score tables, bias score tables under `v_apex_scores/`, and motif MEME files under `weblogo_out/meme/`. Outputs are written to:

`v_apex_scores/methodB_scatter_best_assay_specific/`

Final figure files already included here:

- `methodB_discriminant_scatter_best_loMNase.png/.pdf/.svg`
- `methodB_discriminant_scatter_best_DNase.png/.pdf/.svg`
- `methodB_discriminant_scatter_best_ATAC.png/.pdf/.svg`

