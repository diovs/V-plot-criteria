#!/bin/bash
# =====================================================================
# V_plot_pipeline driver: fragment BED to TF-versus-bias cut-offs and scatter plot
# Version 1.1 (2026-07-10) - full documentation in MANUAL.md
#
#   [0] prepare_fragments  raw fragments (3 cols) -> midpoint+length (6 cols)
#   [1] locate_bias        bias k-mer + genome fa -> bias coordinate beds
#   [2] run_TF_scatter     TF motif  -> V-plot distance distribution
#   [3] run_bias_scatter   bias sites -> V-plot distance distribution (excluding given TFs)
#   [4] fit_vplot_apex     model fit -> TF/bias apex score tsv
#   [5] scatter_cutoffs    apex -> Wilcoxon rank-sum test; if significant, 2-D cut-offs + scatter plot
#
# Usage:
#   Pass everything on the command line (no need to edit config.sh), e.g.:
#     bash run_pipeline.sh \
#       --genome hg38.fa --fragment my_ATAC.bed --tf-dir motif_dir \
#       --bias-kmers my_bias.txt --exclude-tf "motif_dir/CTCF.bed" \
#       --mode ATAC --out results_myATAC
#   config.sh may also be used for defaults (CLI overrides it); copy config.example.sh
#   to config.sh and edit the paths if you prefer a config file.
#
# Flow control:
#   -c, --config FILE    config file used as defaults (default: ./config.sh; ignored if absent)
#   -f, --from N         first step 0-5 (default 0); skip finished heavy steps (e.g. [1] seqkit)
#   -t, --to N           last step 0-5 (default 5)
#   -k, --keep 0|1       intermediates: 1=keep, 0=delete after the run
#   -h, --help           show help
# Inputs / parameters (override config of the same name):
#   --genome FA          genome fasta            (step 1)
#   --fragment BED       raw fragment bed (3col) (step 0)
#   --tf-dir DIR         TF motif dir (6col bed) (step 2)
#   --bias-kmers FILE    bias k-mer list         (steps 1/3)
#   --exclude-tf "..."   TF motif beds to exclude (space-separated, quoted)
#   --exclude-flank BP   bp added to each side when excluding (default 12)
#   --mode M             loMNase|DNase|ATAC|custom (fitting preset; default ATAC)
#                        preset: apex_x[-10,10] apex_y[0,60] frag-min20 xwin150 max-n0;
#                        frag-max: loMNase/DNase=100, ATAC=150 (custom: empty fields fall back to ATAC)
#   --assay NAME         output naming / plot title (default = mode)
#   --out DIR            output dir (default results; created if absent)
#   --threads N          seqkit/fit threads (default 40)
#   --shuf-n N|0|all     bias subsample size (default 200000; 0/all = no subsampling, all sites)
#   --shuf-seed N|none   bias sampling seed (default 42, reproducible; none = true random)
#   --fragl-max N        fragment-length cap (default: no filter)
#   --rank-alpha A       Wilcoxon rank-sum gate alpha for step 5 (default 0.05; both features must pass)
#   custom fit params: --apex-x-lo --apex-x-hi --apex-y-lo --apex-y-hi
#                      --frag-min --frag-max --x-window --max-n --permutations --perm-n
#
# Input formats: fragment = BED3 (chr start end);
#           TF motif = 6-col BED (chr start end name score strand, strand +/-/.; non-6-col aborts);
#           bias-kmers = one sequence per line; exclude-tf = 6-col BED (one or more).
# Output (--out dir): <ASSAY>_TF_apex.tsv / _bias_apex.tsv,
#           <ASSAY>_rank_sum_test.csv, <ASSAY>_conclusion.txt,
#           <ASSAY>_cutoffs.csv, <ASSAY>_scatter.png/.pdf;
#           intermediates in --out/intermediate/ (add -k 0 to delete after the run).
# Dependencies: bedtools(closestBed/intersectBed), seqkit, python3(numpy/pandas/matplotlib/scipy)
# =====================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$HERE/config.sh"
FROM=0
TO=5
KEEP_CLI=""
declare -A OV=()       # command-line overrides

usage() { awk 'NR==1{next} /^[^#]/{exit} {sub(/^#[ ]?/,"");print}' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)         CONFIG="$2"; shift 2 ;;
        -f|--from)           FROM="$2"; shift 2 ;;
        -t|--to)             TO="$2"; shift 2 ;;
        -k|--keep)           KEEP_CLI="$2"; shift 2 ;;
        --genome)            OV[GENOME_FA]="$2"; shift 2 ;;
        --fragment)          OV[FRAGMENT_BED]="$2"; shift 2 ;;
        --tf-dir)            OV[TF_MOTIF_DIR]="$2"; shift 2 ;;
        --bias-kmers)        OV[BIAS_KMER_LIST]="$2"; shift 2 ;;
        --exclude-tf)        OV[EXCLUDE_TF_BEDS]="$2"; shift 2 ;;
        --exclude-flank)     OV[EXCLUDE_FLANK]="$2"; shift 2 ;;
        --mode)              OV[MODE]="$2"; shift 2 ;;
        --assay)             OV[ASSAY]="$2"; shift 2 ;;
        --threads)           OV[THREADS]="$2"; shift 2 ;;
        --shuf-n)            OV[SHUF_N]="$2"; shift 2 ;;
        --shuf-seed)         OV[SHUF_SEED]="$2"; shift 2 ;;
        --fragl-max)         OV[FRAGL_MAX]="$2"; shift 2 ;;
        --rank-alpha)        OV[RANK_ALPHA]="$2"; shift 2 ;;
        --out)               OV[OUTDIR]="$2"; shift 2 ;;
        --apex-x-lo)         OV[APEX_X_LO]="$2"; shift 2 ;;
        --apex-x-hi)         OV[APEX_X_HI]="$2"; shift 2 ;;
        --apex-y-lo)         OV[APEX_Y_LO]="$2"; shift 2 ;;
        --apex-y-hi)         OV[APEX_Y_HI]="$2"; shift 2 ;;
        --frag-min)          OV[FRAG_MIN]="$2"; shift 2 ;;
        --frag-max)          OV[FRAG_MAX]="$2"; shift 2 ;;
        --x-window)          OV[X_WINDOW]="$2"; shift 2 ;;
        --max-n)             OV[MAX_N]="$2"; shift 2 ;;
        --permutations)      OV[PERMUTATIONS]="$2"; shift 2 ;;
        --perm-n)            OV[PERM_N]="$2"; shift 2 ;;
        -h|--help)           usage 0 ;;
        *) echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

# defaults (so the pipeline runs from the command line even without config.sh)
GENOME_FA=""; FRAGMENT_BED=""; TF_MOTIF_DIR=""; BIAS_KMER_LIST=""
EXCLUDE_TF_BEDS=""; EXCLUDE_FLANK=12
MODE="ATAC"; ASSAY=""; THREADS=40; SHUF_N=200000; SHUF_SEED=42; FRAGL_MAX=""
RANK_ALPHA=""
APEX_X_LO=""; APEX_X_HI=""; APEX_Y_LO=""; APEX_Y_HI=""
FRAG_MIN=""; FRAG_MAX=""; X_WINDOW=""; MAX_N=""; PERMUTATIONS=""; PERM_N=""
OUTDIR="results"; KEEP_INTERMEDIATE=1

# config as the default source (sourced only if present); CLI overrides afterwards
if [[ -f "$CONFIG" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG"
fi
if [[ ${#OV[@]} -gt 0 ]]; then
    for k in "${!OV[@]}"; do printf -v "$k" '%s' "${OV[$k]}"; done
fi
# -k overrides KEEP_INTERMEDIATE
[[ -n "$KEEP_CLI" ]] && KEEP_INTERMEDIATE="$KEEP_CLI"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-1}"
ASSAY="${ASSAY:-$MODE}"      # naming/title defaults to MODE
RANK_ALPHA="${RANK_ALPHA:-0.05}"

[[ "$FROM" =~ ^[0-5]$ && "$TO" =~ ^[0-5]$ ]] || {
    echo "ERROR: --from and --to must be integers from 0 to 5" >&2; exit 1;
}
(( FROM <= TO )) || { echo "ERROR: --from must be less than or equal to --to" >&2; exit 1; }
[[ "$KEEP_INTERMEDIATE" =~ ^[01]$ ]] || { echo "ERROR: --keep must be 0 or 1" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be a positive integer" >&2; exit 1; }
[[ "$ASSAY" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "ERROR: --assay may contain only letters, numbers, '.', '_' and '-'" >&2; exit 1;
}
[[ -n "$OUTDIR" ]] || { echo "ERROR: --out must not be empty" >&2; exit 1; }
awk -v x="$RANK_ALPHA" 'BEGIN {
    ok = (x ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ && x > 0 && x < 1)
    exit !ok
}' || { echo "ERROR: --rank-alpha must be a number between 0 and 1" >&2; exit 1; }

run_step() { [[ "$FROM" -le "$1" && "$TO" -ge "$1" ]]; }

# required-input check (only for the steps that will run)
miss=()
if run_step 0 && [[ -z "$FRAGMENT_BED" ]];   then miss+=("--fragment"); fi
if run_step 1 && [[ -z "$GENOME_FA" ]];      then miss+=("--genome"); fi
if run_step 1 && [[ -z "$BIAS_KMER_LIST" ]]; then miss+=("--bias-kmers"); fi
if run_step 2 && [[ -z "$TF_MOTIF_DIR" ]];   then miss+=("--tf-dir"); fi
if [[ ${#miss[@]} -gt 0 ]]; then
    echo "ERROR: missing required argument(s): ${miss[*]} (pass on the command line, or set in config.sh)" >&2
    exit 1
fi

# intermediates go under the output dir; the output dir is created if absent
WORKDIR="$OUTDIR/intermediate"
mkdir -p "$WORKDIR" "$OUTDIR"

# ---- step 4 fitting parameters: pick a preset by MODE (custom: empty fields fall back to ATAC) ----
apply_mode() {
    # ATAC preset (= fallback defaults for custom fields)
    local a_xlo=-10 a_xhi=10 a_ylo=0 a_yhi=60 a_fmin=20 a_fmax=150 \
          a_xwin=150 a_maxn=0 a_perm=0 a_permn=0
    case "$MODE" in
        loMNase|DNase)
            FIT_X_LO=-10; FIT_X_HI=10; FIT_Y_LO=0; FIT_Y_HI=60
            FIT_FMIN=20; FIT_FMAX=100; FIT_XWIN=150; FIT_MAXN=0; FIT_PERM=0; FIT_PERMN=0 ;;
        ATAC)
            FIT_X_LO=$a_xlo; FIT_X_HI=$a_xhi; FIT_Y_LO=$a_ylo; FIT_Y_HI=$a_yhi
            FIT_FMIN=$a_fmin; FIT_FMAX=$a_fmax; FIT_XWIN=$a_xwin; FIT_MAXN=$a_maxn
            FIT_PERM=$a_perm; FIT_PERMN=$a_permn ;;
        custom)
            FIT_X_LO="${APEX_X_LO:-$a_xlo}"; FIT_X_HI="${APEX_X_HI:-$a_xhi}"
            FIT_Y_LO="${APEX_Y_LO:-$a_ylo}"; FIT_Y_HI="${APEX_Y_HI:-$a_yhi}"
            FIT_FMIN="${FRAG_MIN:-$a_fmin}"; FIT_FMAX="${FRAG_MAX:-$a_fmax}"
            FIT_XWIN="${X_WINDOW:-$a_xwin}"; FIT_MAXN="${MAX_N:-$a_maxn}"
            FIT_PERM="${PERMUTATIONS:-$a_perm}"; FIT_PERMN="${PERM_N:-$a_permn}" ;;
        *) echo "ERROR: unknown MODE: '$MODE' (expected loMNase/DNase/ATAC/custom)" >&2; exit 1 ;;
    esac
}
apply_mode
echo "fit MODE    : $MODE  (apex_x[$FIT_X_LO,$FIT_X_HI] apex_y[$FIT_Y_LO,$FIT_Y_HI] frag[$FIT_FMIN,$FIT_FMAX] xwin=$FIT_XWIN max-n=$FIT_MAXN perm=$FIT_PERM/$FIT_PERMN)"

# derived paths
FRAG_MIDP="$WORKDIR/fragments_midP_fragL.bed"
BIAS_BED_DIR="$WORKDIR/bias_bed"
TF_DIST_DIR="$WORKDIR/tf_dist"
BIAS_DIST_DIR="$WORKDIR/bias_dist"
TF_APEX="$OUTDIR/${ASSAY}_TF_apex.tsv"
BIAS_APEX="$OUTDIR/${ASSAY}_bias_apex.tsv"

# fragment-length filter flag (-m); not passed when empty
MFLAG=(); [[ -n "${FRAGL_MAX:-}" ]] && MFLAG=(-m "$FRAGL_MAX")

echo "########## V_plot_pipeline [$ASSAY]  steps $FROM..$TO ##########"

# ---- [0] prepare fragments ----
if run_step 0; then
    echo ">>> [0] prepare_fragments"
    bash "$HERE/0_prepare_fragments.sh" -i "$FRAGMENT_BED" -o "$FRAG_MIDP" "${MFLAG[@]}"
fi

# ---- [1] locate bias ----
if run_step 1; then
    echo ">>> [1] locate_bias (seqkit)"
    bash "$HERE/1_locate_bias.sh" -k "$BIAS_KMER_LIST" -g "$GENOME_FA" -o "$BIAS_BED_DIR" \
        -p "${ASSAY}_Bias" -j "$THREADS"
fi

# ---- [2] TF V-plot ----
if run_step 2; then
    echo ">>> [2] run_TF_scatter"
    bash "$HERE/run_TF_scatter.sh" -i "$TF_MOTIF_DIR" -a "$FRAG_MIDP" -o "$TF_DIST_DIR" \
        -j "$THREADS" "${MFLAG[@]}"
fi

# ---- [3] bias V-plot ----
if run_step 3; then
    echo ">>> [3] run_bias_scatter"
    bash "$HERE/run_bias_scatter.sh" -i "$BIAS_BED_DIR" -a "$FRAG_MIDP" -o "$BIAS_DIST_DIR" \
        -x "$EXCLUDE_TF_BEDS" -f "$EXCLUDE_FLANK" -n "$SHUF_N" -S "$SHUF_SEED" \
        -j "$THREADS" "${MFLAG[@]}"
fi

# ---- [4] fit apex ----
if run_step 4; then
    echo ">>> [4] fit_vplot_apex (TF)"
    python3 "$HERE/fit_vplot_apex.py" --input "$TF_DIST_DIR" --glob "*_fragL_dist.txt" \
        --apex-x-lo "$FIT_X_LO" --apex-x-hi "$FIT_X_HI" \
        --apex-y-lo "$FIT_Y_LO" --apex-y-hi "$FIT_Y_HI" \
        --frag-min "$FIT_FMIN" --frag-max "$FIT_FMAX" --x-window "$FIT_XWIN" \
        --max-n "$FIT_MAXN" --permutations "$FIT_PERM" --perm-n "$FIT_PERMN" \
        --cores "$THREADS" --output "$TF_APEX"
    echo ">>> [4] fit_vplot_apex (bias)"
    python3 "$HERE/fit_vplot_apex.py" --input "$BIAS_DIST_DIR" --glob "*_fragL_dist.txt" \
        --apex-x-lo "$FIT_X_LO" --apex-x-hi "$FIT_X_HI" \
        --apex-y-lo "$FIT_Y_LO" --apex-y-hi "$FIT_Y_HI" \
        --frag-min "$FIT_FMIN" --frag-max "$FIT_FMAX" --x-window "$FIT_XWIN" \
        --max-n "$FIT_MAXN" --permutations "$FIT_PERM" --perm-n "$FIT_PERMN" \
        --cores "$THREADS" --output "$BIAS_APEX"
fi

# ---- [5] 2-D cut-offs + scatter ----
if run_step 5; then
    echo ">>> [5] scatter_cutoffs"
    python3 "$HERE/5_scatter_cutoffs.py" --tf "$TF_APEX" --bias "$BIAS_APEX" \
        --out-prefix "$OUTDIR/${ASSAY}" --assay "$ASSAY" \
        --rank-alpha "$RANK_ALPHA"
fi

# ---- clean up intermediates ----
if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
    echo ">>> deleting intermediates: $WORKDIR"
    rm -rf "$WORKDIR"
else
    echo ">>> intermediates kept in: $WORKDIR"
fi

echo "########## done. results in $OUTDIR/ ##########"
