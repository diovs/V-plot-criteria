# =====================================================================
# V_plot_pipeline configuration example (ATAC)
# Copy this file to config.sh and edit the paths, or pass all paths on the command line.
# Paths may be relative to this directory or absolute.
# =====================================================================

# ---- required: inputs ----
GENOME_FA="path/to/hg38.fa"                 # genome fasta (used by seqkit in step 1)
FRAGMENT_BED="path/to/fragments.bed"        # raw fragment bed (3 cols: chr start end)
TF_MOTIF_DIR="path/to/tf_motif_beds"        # TF motif dir (each file a 6-col motif bed)
BIAS_KMER_LIST="path/to/bias_kmers.txt"     # bias k-mer list (one sequence per line)

# ---- exclude TF from bias sites (step 3) ----
# TF motif beds to remove from bias sites, space-separated; empty = no exclusion
EXCLUDE_TF_BEDS="path/to/CTCF.center.bed"
EXCLUDE_FLANK=12              # bp added to each side of every motif (bias & motif are stranded; same-strand exclusion)

# ---- general parameters ----
MODE="ATAC"                  # fitting preset: loMNase / DNase / ATAC / custom
ASSAY=""                     # plot title / output naming; empty = use the MODE name
THREADS=40                   # seqkit threads / concurrent closestBed jobs / fit processes
SHUF_N=200000               # bias-site subsample size (step 3); 0/all = no subsampling, use all sites
SHUF_SEED=42                 # bias-sampling random seed (fixed = reproducible; none = true random)
FRAGL_MAX=""                 # fragment-length cap; empty = no filter (before closest in steps 2/3)
RANK_ALPHA=0.05              # step 5 Wilcoxon rank-sum gate; both features must have p < alpha before cut-offs are calibrated

# ---- fitting parameters (step 4, fit_vplot_apex.py) ----
# loMNase / DNase / ATAC use built-in presets (defined in run_pipeline.sh):
#   all three: apex_x in [-10,10], apex_y in [0,60], frag-min 20, x-window 150, max-n 0, no permutation;
#   frag-max:  loMNase=100, DNase=100, ATAC=150.
# The values below are read only when MODE=custom; any left empty falls back to the ATAC preset.
APEX_X_LO=""
APEX_X_HI=""
APEX_Y_LO=""
APEX_Y_HI=""
FRAG_MIN=""
FRAG_MAX=""
X_WINDOW=""
MAX_N=""
PERMUTATIONS=""
PERM_N=""

# ---- output ----
OUTDIR="results"             # output dir (created if absent); holds final results + intermediates
                             # intermediates -> $OUTDIR/intermediate/
KEEP_INTERMEDIATE=1         # 1 = keep intermediates; 0 = delete $OUTDIR/intermediate/ after the run
