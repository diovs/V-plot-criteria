#!/usr/bin/env python3
"""
5_scatter_cutoffs.py -- final step: from apex scores to Wilcoxon-rank-sum-gated
2-D cut-offs + scatter.

Reads the two outputs of fit_vplot_apex.py (one TF, one bias) and uses two
discriminating quantities
  - V-channel width  width = apex_y_channel_width
  - inside-V enrichment  E = log2(enrichment_fold)
to compare TF motifs with enzyme cleavage bias controls, reporting:
  1) two-sided Wilcoxon rank-sum tests for the two features
     and confirmation that the TF values are higher than the bias values;
  2) if both features pass, an empirical cut-off per feature
     (gap midpoint / Youden; identical when fully
     separable) and its margin;
  3) leave-one-out (LOO) sensitivity/specificity/accuracy, including the
     two-dimensional rule (both features above their cut-off);
  4) a 2-D scatter (width vs E). If the rank-sum gate fails, no cut-offs are drawn.

Usage:
  python 5_scatter_cutoffs.py --tf TF_apex.tsv --bias bias_apex.tsv \
      --out-prefix results/ATAC [--assay ATAC] [--rank-alpha 0.05]
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Thresholds (rule: larger value = more TF-like)
# ---------------------------------------------------------------------------
def youden_threshold(value, is_tf):
    """Candidate cuts = midpoints of adjacent sorted values; pick the one maximising
    Youden J = sens+spec-1. Predict TF when value >= cut. Equals the gap midpoint
    when the classes are fully separable."""
    v = np.asarray(value, float)
    order = np.unique(v)
    if len(order) == 1:
        return float(order[0])
    cuts = (order[:-1] + order[1:]) / 2.0
    P = is_tf.sum(); N = (~is_tf).sum()
    best_cut, best_J = cuts[0], -np.inf
    for c in cuts:
        pred = v >= c
        tp = np.sum(pred & is_tf); fp = np.sum(pred & ~is_tf)
        sens = tp / P if P else 0.0
        spec = (N - fp) / N if N else 0.0
        J = sens + spec - 1.0
        if J > best_J:
            best_J, best_cut = J, c
    return float(best_cut)


def midgap_threshold(value, is_tf):
    """Gap-midpoint cut-off and margin (half the gap) for fully separable classes."""
    v = np.asarray(value, float)
    bias = v[~is_tf]; tf = v[is_tf]
    maxb = float(np.max(bias)); mintf = float(np.min(tf))
    thr = (maxb + mintf) / 2.0
    margin = (mintf - maxb) / 2.0
    return dict(threshold=thr, margin=margin, separable=bool(mintf > maxb),
                max_bias=maxb, min_tf=mintf)


# ---------------------------------------------------------------------------
# Leave-one-out cross-validation (LOO)
# ---------------------------------------------------------------------------
def loo_1d(value, is_tf):
    """Hold out each point, fit the Youden threshold on the rest, predict the held-out
    point. Returns acc/sens/spec."""
    v = np.asarray(value, float)
    n = len(v); correct = 0
    tp = fp = tn = fn = 0
    for i in range(n):
        mask = np.ones(n, bool); mask[i] = False
        thr = youden_threshold(v[mask], is_tf[mask])
        pred = v[i] >= thr
        truth = is_tf[i]
        correct += int(pred == truth)
        if truth and pred: tp += 1
        elif truth and not pred: fn += 1
        elif (not truth) and pred: fp += 1
        else: tn += 1
    P = tp + fn; N = tn + fp
    return dict(acc=correct / n, sens=tp / P if P else np.nan,
                spec=tn / N if N else np.nan, tp=tp, fp=fp, tn=tn, fn=fn)


def loo_2d(width, enrich, is_tf):
    """LOO for the 2-D rule: each fold fits a width cut and an E cut on the rest,
    predicting the held-out point as TF iff width>=thr_w AND E>=thr_e."""
    w = np.asarray(width, float); e = np.asarray(enrich, float)
    n = len(w); correct = 0
    tp = fp = tn = fn = 0
    for i in range(n):
        mask = np.ones(n, bool); mask[i] = False
        tw = youden_threshold(w[mask], is_tf[mask])
        te = youden_threshold(e[mask], is_tf[mask])
        pred = (w[i] >= tw) and (e[i] >= te)
        truth = is_tf[i]
        correct += int(pred == truth)
        if truth and pred: tp += 1
        elif truth and not pred: fn += 1
        elif (not truth) and pred: fp += 1
        else: tn += 1
    P = tp + fn; N = tn + fp
    return dict(acc=correct / n, sens=tp / P if P else np.nan,
                spec=tn / N if N else np.nan, tp=tp, fp=fp, tn=tn, fn=fn)


def rank_sum_auc(value, is_tf):
    """Analytic rank-sum AUC (larger value = more TF-like)."""
    v = np.asarray(value, float)
    pos = v[is_tf]; neg = v[~is_tf]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(len(v), float); ranks[order] = np.arange(1, len(v) + 1)
    # handle ties with average ranks
    s = np.concatenate([pos, neg])
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt)); np.add.at(sums, inv, ranks); avg = sums / cnt
    ranks = avg[inv]
    r_pos = ranks[:len(pos)].sum()
    return (r_pos - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def format_p(value):
    if not np.isfinite(value):
        return "NA"
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3g}"


def rank_sum_summary(df, is_tf, alpha):
    """Run two-sided Wilcoxon rank-sum tests for the two discriminating features.

    SciPy's mannwhitneyu implementation supplies the equivalent two-sample
    rank-sum p-value. The reported WRS statistic is the TF-group rank sum.
    """
    rows = []
    for feat, label in [("width", "V-channel width (bp)"),
                        ("enrich", "log2(V-in/V-out)")]:
        val = df[feat].to_numpy()
        tf_v = val[is_tf]
        bias_v = val[~is_tf]
        test = stats.mannwhitneyu(tf_v, bias_v, alternative="two-sided", method="auto")
        pvalue = float(test.pvalue)
        median_tf = float(np.median(tf_v))
        median_bias = float(np.median(bias_v))
        auc = rank_sum_auc(val, is_tf)
        tf_higher = bool(median_tf > median_bias and auc > 0.5)
        rows.append(dict(
            feature=label,
            TF_n=len(tf_v),
            bias_n=len(bias_v),
            median_TF=median_tf,
            median_bias=median_bias,
            WRS_statistic=float(test.statistic + len(tf_v) * (len(tf_v) + 1) / 2.0),
            WRS_p=pvalue,
            AUC=auc,
            TF_higher=tf_higher,
            alpha=alpha,
            significant=bool(pvalue < alpha),
            gate_passed=bool(pvalue < alpha and tf_higher),
        ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Load + clean
# ---------------------------------------------------------------------------
def load_apex(path, group):
    df = pd.read_csv(path, sep="\t")
    if "status" in df:
        df = df[df["status"] == "ok"].copy()
    df["width"] = pd.to_numeric(df["apex_y_channel_width"], errors="coerce")
    fold = pd.to_numeric(df["enrichment_fold"], errors="coerce")
    df["enrich"] = np.log2(fold.where(fold > 0))
    df["group"] = group
    df = df.dropna(subset=["width", "enrich"])
    return df[["motif", "width", "enrich", "group"]]


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_scatter(df, tw, te, out_png, out_pdf, assay, rank_sum_passed=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12})
    C_TF, C_BIAS = "#104e8b", "#FF6A6A"
    tf = df[df.group == "TF"]; bi = df[df.group == "bias"]

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    xmax = df.width.max() * 1.12 + 1
    ymax = df.enrich.max() * 1.15
    if tw is not None and te is not None:
        ymax = max(ymax, te * 1.5)
    ymin = min(df.enrich.min() * 1.1, -0.1)
    if tw is not None and te is not None:
        # shaded upper-right acceptance region
        ax.add_patch(plt.Rectangle((tw, te), xmax - tw, ymax - te,
                                   facecolor="#BFD3E6", alpha=0.30, zorder=0))
        ax.axvline(tw, ls="--", lw=1.6, color="#333333", zorder=1)
        ax.axhline(te, ls="--", lw=1.6, color="#333333", zorder=1)
    ax.scatter(bi.width, bi.enrich, s=55, facecolors="none", edgecolors=C_BIAS,
               linewidths=1.7, zorder=3, label="enzyme bias")
    ax.scatter(tf.width, tf.enrich, s=48, c=C_TF, zorder=3, label="TF motif")
    if tw is not None and te is not None:
        ax.text(tw, ymax, f"  width >= {tw:.2f} bp", color="#333333", va="top",
                ha="left", fontsize=10)
        ax.text(xmax, te, f"E >= {te:.2f}  ", color="#333333", va="bottom",
                ha="right", fontsize=10)
    elif not rank_sum_passed:
        ax.text(0.02, 0.98, "rank-sum gate not passed;\nno cut-offs calibrated",
                transform=ax.transAxes, ha="left", va="top", fontsize=10,
                color="#333333",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="#bbbbbb", alpha=0.90))
    ax.set_xlim(min(df.width.min() * 1.1, -0.5), xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("V-channel width (bp)")
    ax.set_ylabel("log2(V-in / V-out)")
    title = f"{assay}: TF motifs vs enzyme bias"
    if tw is not None and te is not None:
        title = f"{assay}: genuine footprint vs enzyme bias"
    ax.set_title(title, fontsize=12.5)
    ax.legend(frameon=False, loc="lower right", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)


def write_conclusion(path, assay, rank_sum_table, rank_sum_passed, alpha, cutoffs=None):
    lines = [
        f"Assay: {assay}",
        f"Wilcoxon rank-sum gate: {'PASS' if rank_sum_passed else 'FAIL'}",
        ("Criterion: both features must have two-sided rank-sum p "
         f"< {alpha:g}, with TF values higher than bias values."),
        "",
        "Feature-level tests:",
    ]
    for _, row in rank_sum_table.iterrows():
        lines.append(
            f"- {row['feature']}: p={format_p(float(row['WRS_p']))}, "
            f"AUC={float(row['AUC']):.3f}, "
            f"median_TF={float(row['median_TF']):.4g}, "
            f"median_bias={float(row['median_bias']):.4g}, "
            f"TF_higher={bool(row['TF_higher'])}, "
            f"gate_passed={bool(row['gate_passed'])}"
        )
    lines.append("")
    if rank_sum_passed and cutoffs is not None:
        lines.append("Conclusion: TF motifs are significantly higher than bias controls in both features; cut-offs were calibrated.")
        lines.append(f"2-D rule: V-channel width >= {cutoffs['width']:.4g} bp AND log2(V-in/V-out) >= {cutoffs['enrich']:.4g}.")
    else:
        lines.append("Conclusion: both required features did not show a significant TF-higher result; no cut-offs were calibrated.")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="apex scores -> Wilcoxon-rank-sum-gated 2D cut-offs + scatter")
    ap.add_argument("--tf", required=True, help="TF apex TSV (fit_vplot_apex.py output)")
    ap.add_argument("--bias", required=True, help="bias apex TSV")
    ap.add_argument("--out-prefix", required=True, help="output prefix")
    ap.add_argument("--assay", default="assay", help="assay name (title only)")
    ap.add_argument("--rank-alpha", type=float, default=0.05,
                    help="Wilcoxon rank-sum significance threshold; both features must pass before cut-offs are calibrated")
    args = ap.parse_args()
    if not 0 < args.rank_alpha < 1:
        ap.error("--rank-alpha must be between 0 and 1")

    tf = load_apex(args.tf, "TF")
    bi = load_apex(args.bias, "bias")
    df = pd.concat([tf, bi], ignore_index=True)
    if len(tf) < 2 or len(bi) < 2:
        sys.exit(f"ERROR: at least 2 valid TF and 2 valid bias rows are required (TF={len(tf)}, bias={len(bi)})")
    is_tf = (df.group == "TF").to_numpy()

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    csv_path = args.out_prefix + "_cutoffs.csv"
    rank_sum_path = args.out_prefix + "_rank_sum_test.csv"
    conclusion_path = args.out_prefix + "_conclusion.txt"
    scatter_png = args.out_prefix + "_scatter.png"
    scatter_pdf = args.out_prefix + "_scatter.pdf"

    rank_sum_table = rank_sum_summary(df, is_tf, args.rank_alpha)
    rank_sum_table.to_csv(rank_sum_path, index=False)
    rank_sum_passed = bool(rank_sum_table["gate_passed"].all())

    if not rank_sum_passed:
        out = rank_sum_table[["feature", "WRS_p", "AUC", "TF_higher", "alpha",
                              "significant", "gate_passed"]].copy()
        out.insert(1, "cutoff", np.nan)
        out["status"] = "no_cutoff_rank_sum_gate_failed"
        out.to_csv(csv_path, index=False)
        make_scatter(df, None, None, scatter_png, scatter_pdf, args.assay, rank_sum_passed=False)
        write_conclusion(conclusion_path, args.assay, rank_sum_table, False, args.rank_alpha)

        sys.stderr.write(f"\n[{args.assay}]  TF n={len(tf)}  bias n={len(bi)}\n")
        sys.stderr.write(f"Wilcoxon rank-sum gate FAILED: no cut-offs calibrated.\n")
        with pd.option_context("display.width", 200, "display.max_columns", 30):
            sys.stderr.write(rank_sum_table.to_string(index=False) + "\n")
        sys.stderr.write(f"wrote: {rank_sum_path}\n       {conclusion_path}\n       {csv_path}\n       {scatter_png} / .pdf\n")
        return

    rows = []
    for feat, label in [("width", "V-channel width (bp)"),
                        ("enrich", "log2(V-in/V-out)")]:
        val = df[feat].to_numpy()
        mg = midgap_threshold(val, is_tf)
        thr = youden_threshold(val, is_tf)
        loo = loo_1d(val, is_tf)
        rank_sum_row = rank_sum_table[rank_sum_table["feature"] == label].iloc[0]
        rows.append(dict(feature=label, cutoff=round(thr, 4),
                         margin=round(mg["margin"], 4), separable=mg["separable"],
                         auc=round(rank_sum_auc(val, is_tf), 4),
                         WRS_p=float(rank_sum_row["WRS_p"]),
                         TF_higher=bool(rank_sum_row["TF_higher"]),
                         WRS_significant=bool(rank_sum_row["significant"]),
                         loo_acc=round(loo["acc"], 4),
                         loo_sens=round(loo["sens"], 4),
                         loo_spec=round(loo["spec"], 4)))

    tw = youden_threshold(df["width"].to_numpy(), is_tf)
    te = youden_threshold(df["enrich"].to_numpy(), is_tf)
    j2 = loo_2d(df["width"].to_numpy(), df["enrich"].to_numpy(), is_tf)
    rows.append(dict(feature="2D rule (width AND E)",
                     cutoff=f"width>={tw:.3f} & E>={te:.3f}", margin=np.nan,
                     separable=np.nan, auc=np.nan,
                     loo_acc=round(j2["acc"], 4),
                     loo_sens=round(j2["sens"], 4),
                     loo_spec=round(j2["spec"], 4)))

    out = pd.DataFrame(rows)
    out.to_csv(csv_path, index=False)
    make_scatter(df, tw, te, scatter_png, scatter_pdf, args.assay, rank_sum_passed=True)
    write_conclusion(conclusion_path, args.assay, rank_sum_table, True, args.rank_alpha,
                     cutoffs={"width": tw, "enrich": te})

    # console summary
    sys.stderr.write(f"\n[{args.assay}]  TF n={len(tf)}  bias n={len(bi)}\n")
    sys.stderr.write(f"Wilcoxon rank-sum gate PASSED: both features have p < {args.rank_alpha:g} and TF-higher direction.\n")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        sys.stderr.write(rank_sum_table.to_string(index=False) + "\n")
        sys.stderr.write(out.to_string(index=False) + "\n")
    sys.stderr.write(f"\n2-D cut-off: width >= {tw:.3f} bp and E >= {te:.3f} (log2)\n")
    sys.stderr.write(f"wrote: {rank_sum_path}\n       {conclusion_path}\n       {csv_path}\n       {scatter_png} / .pdf\n")


if __name__ == "__main__":
    main()
