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
     (implemented with the equivalent Mann-Whitney U test);
  2) if both features are significant, an empirical cut-off per feature
     (gap midpoint / Youden; identical when fully
     separable) and its margin;
  3) leave-one-out (LOO) sensitivity/specificity/accuracy, including the
     two-dimensional rule (both features above their cut-off);
  4) a 2-D scatter (width vs E). If the rank-sum gate fails, no cut-offs are drawn.

Usage:
  python 5_scatter_cutoffs.py --tf TF_apex.tsv --bias bias_apex.tsv \
      --out-prefix results/ATAC_methodB [--assay ATAC] [--mw-alpha 0.05]
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


def auc_mw(value, is_tf):
    """Analytic rank-sum / Mann-Whitney AUC (larger value = more TF-like)."""
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


def mann_whitney_summary(df, is_tf, alpha):
    """Run two-sided Wilcoxon rank-sum tests for the two discriminating features.

    SciPy's mannwhitneyu function is used because the Wilcoxon rank-sum test and
    Mann-Whitney U test are equivalent for two independent groups. The legacy
    MW_p column is retained for backward compatibility; WRS_p is the preferred
    public-facing name.
    """
    rows = []
    for feat, label in [("width", "V-channel width (bp)"),
                        ("enrich", "log2(V-in/V-out)")]:
        val = df[feat].to_numpy()
        tf_v = val[is_tf]
        bias_v = val[~is_tf]
        mw = stats.mannwhitneyu(tf_v, bias_v, alternative="two-sided", method="auto")
        pvalue = float(mw.pvalue)
        rows.append(dict(
            feature=label,
            TF_n=len(tf_v),
            bias_n=len(bias_v),
            median_TF=float(np.median(tf_v)),
            median_bias=float(np.median(bias_v)),
            U=float(mw.statistic),
            WRS_p=pvalue,
            MW_p=pvalue,
            AUC=auc_mw(val, is_tf),
            alpha=alpha,
            significant=bool(pvalue < alpha),
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
def make_scatter(df, tw, te, out_png, out_pdf, assay, mw_passed=True):
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
    elif not mw_passed:
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


def write_conclusion(path, assay, mw_table, mw_passed, alpha, cutoffs=None):
    lines = [
        f"Assay: {assay}",
        f"Wilcoxon rank-sum gate: {'PASS' if mw_passed else 'FAIL'}",
        f"Criterion: both V-channel width and inside-V enrichment must have two-sided rank-sum p < {alpha:g}.",
        "",
        "Feature-level tests:",
    ]
    for _, row in mw_table.iterrows():
        lines.append(
            f"- {row['feature']}: p={format_p(float(row['WRS_p']))}, "
            f"AUC={float(row['AUC']):.3f}, "
            f"median_TF={float(row['median_TF']):.4g}, "
            f"median_bias={float(row['median_bias']):.4g}, "
            f"significant={bool(row['significant'])}"
        )
    lines.append("")
    if mw_passed and cutoffs is not None:
        lines.append("Conclusion: TF motifs and bias controls differ significantly in both features; cut-offs were calibrated.")
        lines.append(f"2-D rule: V-channel width >= {cutoffs['width']:.4g} bp AND inside-V enrichment >= {cutoffs['enrich']:.4g}.")
    else:
        lines.append("Conclusion: TF motifs and bias controls do not differ significantly in both required features; no cut-offs were calibrated.")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="apex scores -> Wilcoxon-rank-sum-gated 2D cut-offs + scatter")
    ap.add_argument("--tf", required=True, help="TF apex TSV (fit_vplot_apex.py output)")
    ap.add_argument("--bias", required=True, help="bias apex TSV")
    ap.add_argument("--out-prefix", required=True, help="output prefix")
    ap.add_argument("--assay", default="assay", help="assay name (title only)")
    ap.add_argument("--rank-alpha", "--mw-alpha", dest="mw_alpha", type=float, default=0.05,
                    help="Wilcoxon rank-sum significance threshold; both features must pass before cut-offs are calibrated")
    args = ap.parse_args()

    tf = load_apex(args.tf, "TF")
    bi = load_apex(args.bias, "bias")
    df = pd.concat([tf, bi], ignore_index=True)
    if len(tf) == 0 or len(bi) == 0:
        sys.exit(f"ERROR: TF or bias has 0 rows (TF={len(tf)}, bias={len(bi)})")
    is_tf = (df.group == "TF").to_numpy()

    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    csv_path = args.out_prefix + "_cutoffs.csv"
    mw_path = args.out_prefix + "_mw_test.csv"
    conclusion_path = args.out_prefix + "_conclusion.txt"
    scatter_png = args.out_prefix + "_scatter.png"
    scatter_pdf = args.out_prefix + "_scatter.pdf"

    mw_table = mann_whitney_summary(df, is_tf, args.mw_alpha)
    mw_table.to_csv(mw_path, index=False)
    mw_passed = bool(mw_table["significant"].all())

    if not mw_passed:
        out = mw_table[["feature", "WRS_p", "MW_p", "AUC", "alpha", "significant"]].copy()
        out.insert(1, "cutoff", np.nan)
        out["status"] = "no_cutoff_MW_gate_failed"
        out.to_csv(csv_path, index=False)
        make_scatter(df, None, None, scatter_png, scatter_pdf, args.assay, mw_passed=False)
        write_conclusion(conclusion_path, args.assay, mw_table, False, args.mw_alpha)

        sys.stderr.write(f"\n[{args.assay}]  TF n={len(tf)}  bias n={len(bi)}\n")
        sys.stderr.write(f"Wilcoxon rank-sum gate FAILED: no cut-offs calibrated.\n")
        with pd.option_context("display.width", 200, "display.max_columns", 30):
            sys.stderr.write(mw_table.to_string(index=False) + "\n")
        sys.stderr.write(f"wrote: {mw_path}\n       {conclusion_path}\n       {csv_path}\n       {scatter_png} / .pdf\n")
        return

    rows = []
    for feat, label in [("width", "V-channel width (bp)"),
                        ("enrich", "log2(V-in/V-out)")]:
        val = df[feat].to_numpy()
        mg = midgap_threshold(val, is_tf)
        thr = youden_threshold(val, is_tf)
        loo = loo_1d(val, is_tf)
        mw_row = mw_table[mw_table["feature"] == label].iloc[0]
        rows.append(dict(feature=label, cutoff=round(thr, 4),
                         margin=round(mg["margin"], 4), separable=mg["separable"],
                         auc=round(auc_mw(val, is_tf), 4),
                         WRS_p=float(mw_row["WRS_p"]),
                         MW_p=float(mw_row["MW_p"]),
                         MW_significant=bool(mw_row["significant"]),
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
    make_scatter(df, tw, te, scatter_png, scatter_pdf, args.assay, mw_passed=True)
    write_conclusion(conclusion_path, args.assay, mw_table, True, args.mw_alpha,
                     cutoffs={"width": tw, "enrich": te})

    # console summary
    sys.stderr.write(f"\n[{args.assay}]  TF n={len(tf)}  bias n={len(bi)}\n")
    sys.stderr.write(f"Wilcoxon rank-sum gate PASSED: both features have p < {args.mw_alpha:g}.\n")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        sys.stderr.write(mw_table.to_string(index=False) + "\n")
        sys.stderr.write(out.to_string(index=False) + "\n")
    sys.stderr.write(f"\n2-D cut-off: width >= {tw:.3f} bp and E >= {te:.3f} (log2)\n")
    sys.stderr.write(f"wrote: {mw_path}\n       {conclusion_path}\n       {csv_path}\n       {scatter_png} / .pdf\n")


if __name__ == "__main__":
    main()
