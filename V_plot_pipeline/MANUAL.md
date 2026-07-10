# V_plot_pipeline — Manual

**Version 1.1 · 2026-07-10**

A pipeline that decides, for each transcription-factor (TF) motif, whether its
V-plot reflects a **TF-footprint-like pattern** or merely **enzyme sequence-cleavage
bias**, by fitting a V-plot model and applying a **two-dimensional cut-off**
(V-channel width × inside-V enrichment) learned from matched bias controls.

It runs end-to-end from a fragment BED file to a scatter plot and a table of
cut-offs, for ATAC-seq, DNase-seq or loMNase-seq.

---

## 1. Synopsis

```bash
bash run_pipeline.sh \
  --genome      GENOME.fa \
  --fragment    FRAGMENTS.bed \
  --tf-dir      TF_MOTIF_DIR/ \
  --bias-kmers  BIAS_KMERS.txt \
  --exclude-tf  "TF_MOTIF_DIR/CTCF.bed" \
  --mode        ATAC \
  --rank-alpha  0.05 \
  --out         RESULTS_DIR
```

All parameters can be passed on the command line (shown above) **or** placed in
`config.sh` as defaults; command-line values override `config.sh`. To use a config
file, copy `config.example.sh` to `config.sh` and edit the paths. Either command-line
arguments or a config file may be used.

---

## 2. Description

A **V-plot** is the 2-D distribution of sequenced fragments around a set of motif
occurrences: the x-axis is the signed distance from a fragment midpoint to the
motif centre, the y-axis is the fragment length. A bound TF protects its motif and
deflects fragment ends, producing a characteristic **"V"** whose apex sits at the
motif centre. Enzyme cleavage bias can also produce a V-like pattern, so the V shape
alone does not prove a TF-footprint-like pattern.

This pipeline distinguishes the two using two quantities extracted from a fitted
V-plot model (see [§9 Method](#9-method)):

1. **V-channel width** `b` — the apex y-coordinate, i.e. the length of the region
   protected by the bound factor;
2. **Inside-V enrichment** `E = log2(V-in / V-out)` — how much fragment density is
   concentrated inside the fitted V-shaped region versus outside.

TF-footprint-like V-plots have a resolvable channel (large `b`) **and** strong enrichment
(large `E`); enzyme-bias "V"s have `b ≈ 0` and `E ≈ 0`. The final step first uses
two-sided Wilcoxon rank-sum tests to ask whether TF motifs and bias controls differ
significantly in both features. Only when both tests pass the configured alpha and
the TF values are higher than the bias values does the pipeline calibrate cut-offs
from user-supplied bias k-mer controls and
validate them by leave-one-out cross-validation (LOO). A motif is called
TF-footprint-like only if **both** features exceed their cut-offs.

---

## 3. Dependencies

| Tool | Version |
|------|---------|
| `bash` | >= 4.0 |
| `bedtools` | >= 2.30.0 |
| `seqkit` | >= 2.0.0 |
| GNU `sort` / `shuf` | coreutils >= 8.25 |
| `python3` | >= 3.9.5 |
| `numpy` | >= 1.22.3 |
| `pandas` | >= 1.4.2 |
| `scipy` | >= 1.8.0 |
| `matplotlib` | >= 3.5.2 |

Versions above are those tested on our server; equal or newer should also work.
Check yours with:
```bash
command -v closestBed intersectBed seqkit shuf python3
```

---

## 4. Input files

### 4.1 Fragment BED (`--fragment`) — **required**
Raw sequenced fragments, **BED3** (tab-separated):
```
chr    start    end
chr1   10015    10118
chr1   10019    10090
```
- One row per fragment; `start`/`end` are the 5′ and 3′ fragment ends (0-based, BED).
- The pipeline derives the fragment **midpoint** and **length** itself (step 0).
- Chromosome names must match the genome FASTA and the motif BEDs.

### 4.2 TF motif directory (`--tf-dir`) — **required**
A directory; **every file in it** is treated as one motif set (filename/extension
is ignored — `CTCF.bed`, `CTCF.motif`, … all accepted). Each file must be a
**6-column BED**:
```
chr    start    end    name    score    strand
chr1   268007   268008  12      12       +
```
- Column 6 (`strand`) must be `+`, `-`, or `.`. Unstranded (`.`) sites use genomic
  left-to-right orientation when signed distances are calculated.
- Typically each interval is the single-base **motif centre/anchor**.
- **Any file that is not 6-column triggers an error and stops the run** — keep the
  directory clean (motif BEDs only).

### 4.3 Bias k-mer list (`--bias-kmers`) — **required for step 1**
Plain text, **one sequence per line** — the enzyme's preferred cleavage k-mers
(negative controls):
```
GCC
GGC
GCT
GTG
GTC
```
Sequences are case-insensitive and may contain only `A`, `C`, `G`, and `T`.
Blank lines and lines beginning with `#` are ignored; duplicate sequences cause an
error rather than silently overwriting an output BED.

### 4.4 Genome FASTA (`--genome`) — **required for step 1**
This is the genome assembly that the fragments were mapped to (e.g. `hg38.fa`).
`seqkit` uses it to locate every genomic occurrence of each bias k-mer.

### 4.5 TF motifs to exclude (`--exclude-tf`) — optional
One or more **6-column BED** files (space-separated, quoted), whose regions are
removed from the bias sites so that a bias k-mer accidentally lying inside a real
TF motif site cannot contaminate the bias control. Exclusion is **same-strand** and
requires a strand column.
Because this option is parsed as a space-separated list, BED paths used here must
not themselves contain spaces.

---

## 5. Basic usage

```bash
# foreground
bash run_pipeline.sh --genome hg38.fa --fragment my_ATAC.bed \
     --tf-dir motif_dir --bias-kmers my_bias.txt \
     --exclude-tf "motif_dir/CTCF.bed" --mode ATAC --out results_myATAC

# background (recommended for large fragment files)
nohup bash run_pipeline.sh ... > run.log 2>&1 &
tail -f run.log
```

Run a subset of steps (resume / re-plot):
```bash
bash run_pipeline.sh ... -f 2     # start at step 2 (reuse step 0/1 output)
bash run_pipeline.sh ... -f 5     # only re-run the rank-sum gate, scatter, and cut-offs
```

---

## 6. Options

### Input / output
| Option | Arg | Default | Meaning |
|--------|-----|---------|---------|
| `--genome` | FA | — | genome FASTA (step 1) |
| `--fragment` | BED | — | raw fragment BED3 (step 0) |
| `--tf-dir` | DIR | — | directory of 6-col TF motif BEDs (step 2) |
| `--bias-kmers` | FILE | — | bias k-mer list, one per line (step 1) |
| `--out` | DIR | `results` | output directory (created if absent) |
| `--assay` | NAME | = `--mode` | label used in output names / plot title |
| `--rank-alpha` | FLOAT | `0.05` | step-5 Wilcoxon rank-sum gate alpha; both features must pass before cut-offs are calibrated |

### Bias-site exclusion (step 3)
| Option | Arg | Default | Meaning |
|--------|-----|---------|---------|
| `--exclude-tf` | "BED ..." | "" (none) | TF motif BEDs whose regions are removed from bias sites (space-separated, quoted) |
| `--exclude-flank` | BP | `12` | bp added to each side of every excluded motif |

### Fitting preset (step 4)
| Option | Arg | Default | Meaning |
|--------|-----|---------|---------|
| `--mode` | loMNase \| DNase \| ATAC \| custom | `ATAC` | parameter preset for the model fit |

Preset values:

| mode | apex_x | apex_y (`b`) | frag-min | frag-max | x-window | max-n | perm |
|------|--------|--------------|----------|----------|----------|-------|------|
| loMNase | [-10, 10] | [0, 60] | 20 | 100 | 150 | 0 (all) | off |
| DNase   | [-10, 10] | [0, 60] | 20 | 100 | 150 | 0 (all) | off |
| ATAC    | [-10, 10] | [0, 60] | 20 | **150** | 150 | 0 (all) | off |
| custom  | from the options below; **any left unset falls back to the ATAC value** |

`custom` fit options (passed through to `fit_vplot_apex.py`):

| Option | Meaning |
|--------|---------|
| `--apex-x-lo` / `--apex-x-hi` | search bounds for the apex distance `a` (motif-centred; preset [-10, 10]) |
| `--apex-y-lo` / `--apex-y-hi` | search bounds for the apex y = V-channel width `b` (preset [0, 60]) |
| `--frag-min` / `--frag-max` | fragment-length range kept for the fit |
| `--x-window` | maximum distance considered around each motif or bias site |
| `--max-n` | cap on fragments per V-plot used in the fit (`0` = use all) |
| `--permutations` | permutation reps for a calibrated "is there a V" p-value (`0` = off, χ² screen only) |
| `--perm-n` | fixed calibration-subsample size used for both the observed LR and each permutation (`0` = all) |

When `--mode` is loMNase/DNase/ATAC these take the preset values in the table above;
under `--mode custom` any option not given falls back to the ATAC preset.

### Bias sampling (step 3)
| Option | Arg | Default | Meaning |
|--------|-----|---------|---------|
| `--shuf-n` | N \| 0 \| all | `200000` | subsample N bias sites per k-mer; `0`/`all` = use **all** sites (no subsampling) |
| `--shuf-seed` | N \| none | `42` | fixed seed → reproducible sampling; `none` = true random |

### Flow control & misc
| Option | Arg | Default | Meaning |
|--------|-----|---------|---------|
| `-c`, `--config` | FILE | `config.sh` | config file used as defaults (ignored if absent; copy from `config.example.sh`) |
| `-f`, `--from` | 0–5 | `0` | first step to run |
| `-t`, `--to` | 0–5 | `5` | last step to run |
| `-k`, `--keep` | 0 \| 1 | `1` | keep (`1`) or delete (`0`) intermediates after the run |
| `--threads` | N | `40` | threads/jobs for `seqkit`, `closestBed`, and the model fit |
| `--fragl-max` | N | none | drop fragments longer than N before distance calc (steps 2/3) |
| `-h`, `--help` | | | print built-in help |

---

## 7. Pipeline steps

| # | Script | In → Out |
|---|--------|----------|
| 0 | `0_prepare_fragments.sh` | fragment BED3 → 6-col `chr mid mid+1 row fragL .` (sorted) |
| 1 | `1_locate_bias.sh` | k-mer list + genome → one BED per k-mer (`seqkit locate`) |
| 2 | `run_TF_scatter.sh` | TF motifs + fragments → per-TF distance file (`closestBed`) |
| 3 | `run_bias_scatter.sh` | bias sites + fragments → per-k-mer distance file (sample / exclude TF / `closestBed`) |
| 4 | `fit_vplot_apex.py` | distance files → per-motif apex score TSV (TF and bias) |
| 5 | `5_scatter_cutoffs.py` | TF + bias scores → Wilcoxon rank-sum gate; if passed, scatter plot + 2-D cut-off table |

Each individual script also runs stand-alone (`bash <script> -h`).

---

## 8. Output files

All under `--out` (intermediates under `--out/intermediate/`).

### 8.1 `<ASSAY>_TF_apex.tsv` and `<ASSAY>_bias_apex.tsv`
One row per motif (TF) or per k-mer (bias). Key columns:

| Column | Meaning |
|--------|---------|
| `motif` | motif / k-mer name |
| `n` | fragments used in the fit |
| `apex_x` | apex position on the distance axis (≈ 0 for a centred V-plot) |
| `apex_y_channel_width` | apex y = **V-channel width `b`** |
| `pi_enrichment` | fitted V-signal fraction used by the model |
| `enrichment_fold` | empirical V-in / V-out density fold (→ `E = log2(fold)`) |
| `LR`, `p_chi2`, `p_perm`, `q_chi2` | likelihood-ratio test that "a V exists" |
| `se_apex_*`, `apex_*_lo/hi` | standard errors / 95% CI of the apex |
| `status` | `ok`, or a reason the motif was skipped |

### 8.2 `<ASSAY>_rank_sum_test.csv`
Wilcoxon rank-sum gate before cut-off calibration.

| Column | Meaning |
|--------|---------|
| `feature` | `V-channel width` or `log2(V-in/V-out)` |
| `median_TF` / `median_bias` | group medians |
| `WRS_statistic`, `WRS_p` | TF-group rank sum and two-sided Wilcoxon rank-sum p-value |
| `AUC` | rank-sum AUC (1.0 = perfect separation) |
| `TF_higher` | whether the TF median and AUC both indicate larger TF values |
| `alpha` | configured significance threshold |
| `significant` | whether `WRS_p < alpha` |
| `gate_passed` | whether the feature is significant and `TF_higher = True` |

Both features must pass the gate before the pipeline calibrates cut-offs.

### 8.3 `<ASSAY>_conclusion.txt`
A short text conclusion stating whether the rank-sum gate passed. If it failed, no
cut-offs are calibrated.

### 8.4 `<ASSAY>_cutoffs.csv`
The learned thresholds and their validation. If the rank-sum gate fails, this file is
still written, but `status = no_cutoff_rank_sum_gate_failed` and `cutoff` is blank.

| Column | Meaning |
|--------|---------|
| `feature` | `V-channel width`, `log2(V-in/V-out)`, or `2D rule` |
| `cutoff` | the threshold (mid-gap = Youden when classes separate) |
| `margin` | half the gap between the highest bias and the lowest TF |
| `separable` | whether bias and TF are perfectly separated on this feature |
| `auc` | rank-sum AUC (1.0 = perfect separation) |
| `WRS_p` / `WRS_significant` / `TF_higher` | feature-level rank-sum result used by the gate |
| `loo_acc` / `loo_sens` / `loo_spec` | leave-one-out accuracy / sensitivity / specificity |

### 8.5 `<ASSAY>_scatter.png` / `.pdf`
Scatter of **V-channel width (x)** vs **log2(V-in/V-out) (y)**, TF motifs as filled
blue points, bias k-mers as open red points. If the rank-sum gate passes, dashed lines
mark the two cut-offs and the shaded **upper-right region = TF-footprint-like**
(both features above cut-off). If the rank-sum gate fails, the scatter is drawn without
cut-off lines.

### 8.6 Interpreting the result
- **First check the Wilcoxon rank-sum gate**: both features must have `WRS_p < alpha`
  and `TF_higher = True`. If not, the conclusion file reports that the required
  TF-higher result was not established and no cut-off should be used.
- **TF-footprint-like call**: a motif is called TF-footprint-like when `width ≥ width-cutoff` **and**
  `E ≥ E-cutoff` (it lands in the shaded region).
- **Good separation** looks like: bias points clustered near the origin
  (`width ≈ 0`, `E ≈ 0`), TF points up and to the right, `AUC = 1.0`, `LOO = 100%`.
- The **E cut-off (~0.5, ≈1.4-fold)** is fairly transferable across enzymes; the
  **width cut-off is assay-specific** (the channel widens with the enzyme's
  fragment-size range).
- If TF points fall in the lower-left with the bias, that motif behaves like bias in
  this assay (no TF-footprint-like V-plot), not a pipeline error.

---

## 9. Method

**V-plot model fitting.** For each TF motif set or bias-control k-mer set, the
pipeline builds a V-plot from fragment midpoint distance and fragment length. It
then fits a V-plot model to estimate:

- `apex_x`: the fitted V position on the distance axis;
- `apex_y_channel_width`: the fitted **V-channel width**;
- `enrichment_fold`: the fitted V-in / V-out density ratio, reported as
  **inside-V enrichment** `E = log2(V-in / V-out)`;
- `LR`, `p_chi2`, and optional `p_perm`: evidence that a V-like pattern is
  detectable at the input motif or bias-control coordinates.

The likelihood-ratio statistic answers a limited question: **is there evidence
for a V-like pattern?** It does not by itself distinguish TF-footprint-like patterns from
enzyme-bias patterns, because both TF motifs and bias controls can produce
detectable V-like shapes.

**TF-versus-bias discrimination.** The final classification therefore uses the
two fitted graphical features, V-channel width and inside-V enrichment:

- **Wilcoxon rank-sum gate**: two-sided Wilcoxon rank-sum tests first compare TF motifs
  with bias controls for V-channel width and inside-V enrichment separately.
  Cut-offs are calibrated only if both features have `p < alpha` (default `0.05`)
  and the TF values are higher than the bias values. This direction check is needed
  because the final classification rule assumes that larger values are more
  TF-footprint-like. If either feature fails, the pipeline does not assign cut-offs.
- **Cut-off calibration**: when the rank-sum gate passes, each feature cut-off
  is chosen as the Youden threshold, the value maximizing sensitivity + specificity
  - 1 on the training data. If the two groups are fully separated, this threshold
  is the midpoint between the highest bias-control value and the lowest TF-motif
  value.
- **LOO cross-validation**: each motif or bias control is held out once, the
  cut-offs are refit on the remaining data, and the held-out item is classified.
  The pipeline reports sensitivity, specificity, and accuracy.
- **Two-feature rule**: a motif set is classified as TF-footprint-like only when
  `V-channel width ≥ width cut-off` **and**
  `inside-V enrichment ≥ enrichment cut-off`.

---

## 10. FAQ / troubleshooting

**`ERROR: not a valid 6-col motif bed ...` (step 2 aborts).**
A file in `--tf-dir` is not 6-column. Remove non-motif files; ensure every motif BED
is `chr start end name score strand` with strand `+/-/.`.

**`bedtools` / `seqkit` / `python3: command not found`.**
Install the missing dependency or load its module; see [§3](#3-dependencies).

**`WARNING: openssl not found ... using true random this run`.**
Reproducible seeding needs `openssl`. Without it, bias sampling is still correct but
not byte-reproducible. Install openssl, or accept true-random (or use `--shuf-n 0` to
remove sampling randomness entirely).

**Results differ slightly between runs (e.g. a cut-off shifts by < 1 bp).**
Only step 3 (bias) is stochastic, via random subsampling of bias sites. The TF side
is deterministic. With the default fixed seed (`--shuf-seed 42`) repeated runs are
identical; the small jitter affects only the near-zero bias widths. Use
`--shuf-n 0`/`all` to remove sampling, or a fixed seed for reproducibility.

**Chromosome-name / empty-output problems.**
`closestBed` returns nothing if chrom names mismatch between the fragment BED, motif
BEDs, and genome FASTA (e.g. `chr1` vs `1`). Make them consistent.

**A bias control looks V-shaped but has no resolvable width.**
That is the expected, informative result: enzyme-bias "V"s lack a wide V-channel
(`width ≈ 0`), which is exactly what separates them from TF-footprint-like patterns.

---

## 11. Performance & resources

- **`closestBed` (steps 2 & 3) dominates** run time and memory; cost scales with the
  **fragment BED size** (often the largest input, multi-GB for deep ATAC). Run on a
  node with ample RAM and use `--threads`.
- **Bias k-mers can have millions of genome hits.** The default `--shuf-n 200000`
  subsamples them for speed and to balance site counts against TFs. `--shuf-n 0`
  (use all sites) is much heavier — only on a capable machine.
- **`--max-n` (fit)** caps fragments per V-plot during model fitting; the assay
  presets use `0` (all). Lower it to speed up the fit on very dense motifs.
- **Disk**: intermediates (midpoint BED, per-k-mer bias BEDs, distance files) live in
  `--out/intermediate/`. Add `-k 0` to delete them after a successful run.
- **Resume**: heavy upstream steps need not be repeated — use `-f`/`-t` to re-run only
  what changed (e.g. `-f 5` to re-plot from existing apex tables).
