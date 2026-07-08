#!/usr/bin/env python3
"""
fit_vplot_apex.py  -- Method B: generative-model V-plot apex finder.

For each input file (= one V-plot; one TF or one bias motif), fit a mixture
model of a V-cone signal + uniform background by maximum likelihood, and report
  - the apex (apex_x, apex_y) and its standard errors / 95% CI (from the Hessian)
  - a likelihood-ratio test of "a V exists" (LR statistic, chi-square p, AIC)
  - optional permutation-calibrated p (shuffling distance breaks the V)
  - the signal fraction pi (inside-V enrichment) and apex_y (= V-channel width)

Input file: 3 columns, no header by default -> (site_id, fragment, distance).
Only fragment (col 2) and distance (col 3) are used; all rows of a file are
pooled into one V-plot.

MODEL  (conditional density of distance x given fragment length y; apex (a,b)):
    w(y)   = (y - b) / 2                          # cone half-width, floored
    f_sig  = [s((x-(a-w))/tau) - s((x-(a+w))/tau)] / (2w)   # smooth top-hat
    f_bg   = 1 / (2X)                             # uniform background on [-X, X]
    p(x|y) = pi * f_sig + (1-pi) * f_bg
  where s = logistic. The smooth top-hat integrates to 1 and is differentiable,
  so a numerical Hessian at the MLE gives well-defined standard errors.

NULL model: pi = 0  ->  p(x|y) = f_bg  ->  ll_null = -n * log(2X).
  LR = 2*(ll_full - ll_null). Reference chi-square has df = 3 (a, b, pi) but is
  only APPROXIMATE: pi=0 sits on the parameter-space boundary and (a,b) are not
  identified under the null (Davies' problem). Use --permutations for an honest,
  calibrated p-value; treat p_chi2 as a fast screen only.

Example (local quick test):
  python fit_vplot_apex.py --input CTCF_midP_fragL_dist.txt \
      --frag-min 20 --frag-max 75 --x-window 150 --max-n 200000 \
      --plot-dir plots --output ctcf_apex.tsv

Example (server, all files in a dir, full data, permutation-calibrated):
  python fit_vplot_apex.py --input data_dir --glob "*_fragL_dist.txt" \
      --frag-min 20 --frag-max 75 --x-window 150 --max-n 0 \
      --permutations 500 --perm-n 50000 --cores 40 --output apex_scores.tsv
"""
import argparse
import glob
import math
import os
import sys
import numpy as np
import pandas as pd
from scipy.special import expit          # logistic sigmoid
from scipy.optimize import minimize
from scipy import stats

# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------
def neg_loglik(params, x, y, X, tau, w_floor):
    """Negative log-likelihood of the V-cone + background mixture.
    params = (a, b, logit_pi). x, y are distance / fragment arrays."""
    a, b, lpi = params
    pi = expit(lpi)
    w = np.maximum((y - b) * 0.5, w_floor)               # cone half-width
    # smooth top-hat centred at a, plateau [a-w, a+w], integrates to 1
    f_sig = (expit((x - (a - w)) / tau) - expit((x - (a + w)) / tau)) / (2.0 * w)
    f_bg = 1.0 / (2.0 * X)
    dens = pi * f_sig + (1.0 - pi) * f_bg
    np.maximum(dens, 1e-300, out=dens)                   # guard log(0)
    return -np.sum(np.log(dens))


def loglik_null(n, X):
    """Log-likelihood of the pure-background (pi=0) model."""
    return -n * math.log(2.0 * X)


# ---------------------------------------------------------------------------
# Numerical Hessian (central differences) -> standard errors at the MLE
# ---------------------------------------------------------------------------
def numerical_hessian(fun, p, eps_scale=1e-4):
    p = np.asarray(p, float)
    k = len(p)
    h = eps_scale * (np.abs(p) + eps_scale)
    H = np.zeros((k, k))
    for i in range(k):
        for j in range(i, k):
            pp = p.copy(); pp[i] += h[i]; pp[j] += h[j]
            pm = p.copy(); pm[i] += h[i]; pm[j] -= h[j]
            mp = p.copy(); mp[i] -= h[i]; mp[j] += h[j]
            mm = p.copy(); mm[i] -= h[i]; mm[j] -= h[j]
            H[i, j] = (fun(pp) - fun(pm) - fun(mp) + fun(mm)) / (4.0 * h[i] * h[j])
            H[j, i] = H[i, j]
    return H


# ---------------------------------------------------------------------------
# Fit one V-plot
# ---------------------------------------------------------------------------
def fit_apex(x, y, X, tau, w_floor, seed=1, ax_lo=-10.0, ax_hi=10.0,
             ay_lo=None, ay_hi=None):
    """Fit the mixture model with a few restarts; return MLE + diagnostics.
    apex_x is constrained to [ax_lo, ax_hi] -- a genuine footprint apex must sit
    at the motif centre; without this bound the apex of a no-V (pi->0) motif is
    unidentified and wanders far from 0.
    apex_y is constrained to [ay_lo, ay_hi]; if either is None it defaults to the
    data-driven range (min_fragment - 10  /  max_fragment)."""
    n = len(x)
    ymin, ymax = float(np.min(y)), float(np.max(y))
    b_lo = (ymin - 10.0) if ay_lo is None else float(ay_lo)
    b_hi = ymax if ay_hi is None else float(ay_hi)
    nll = lambda p: neg_loglik(p, x, y, X, tau, w_floor)

    # restarts over apex_x and apex_y starts (b near the bottom of the cone),
    # clipped into the apex_y bounds so every start is feasible
    a0_grid = sorted({0.5 * ax_lo, 0.0, 0.5 * ax_hi})
    b0_grid = [float(np.clip(v, b_lo, b_hi))
               for v in (ymin, np.quantile(y, 0.10), np.quantile(y, 0.25))]
    bounds = [(ax_lo, ax_hi), (b_lo, b_hi), (-12.0, 8.0)]

    best = None
    for a0 in a0_grid:
        for b0 in b0_grid:
            p0 = np.array([a0, float(b0), 0.0])           # logit(pi)=0 -> pi=0.5
            try:
                res = minimize(nll, p0, method="L-BFGS-B", bounds=bounds)
            except Exception:
                continue
            if res.success or np.isfinite(res.fun):
                if best is None or res.fun < best.fun:
                    best = res
    if best is None:
        return None

    a, b, lpi = best.x
    pi = float(expit(lpi))
    ll_full = -best.fun
    ll_null = loglik_null(n, X)
    # nested models -> LR >= 0 in theory; tiny negatives are boundary/numeric
    # artifacts of pi pinned at its floor (i.e. "no V"), so clamp at 0.
    LR = max(2.0 * (ll_full - ll_null), 0.0)
    df = 3
    p_chi2 = float(stats.chi2.sf(LR, df)) if LR > 0 else 1.0
    aic_full = 2 * df - 2 * ll_full
    aic_null = 2 * 0 - 2 * ll_null

    # standard errors from the observed information (Hessian of the NLL)
    se_a = se_b = np.nan
    try:
        H = numerical_hessian(nll, best.x)
        cov = np.linalg.inv(H)
        if np.all(np.isfinite(cov)):
            d = np.diag(cov)
            if d[0] > 0:
                se_a = math.sqrt(d[0])
            if d[1] > 0:
                se_b = math.sqrt(d[1])
    except np.linalg.LinAlgError:
        pass

    # empirical inside-cone vs outside-cone density ratio at the fitted apex
    # (= the "V-in / V-out" fold; ~5 for a real footprint, ~1 for bias)
    w = np.maximum((y - b) * 0.5, w_floor)
    inside = np.abs(x - a) <= w
    width_in = np.sum(2.0 * w)
    width_out = np.sum(2.0 * X - 2.0 * w)
    dens_in = inside.sum() / width_in if width_in > 0 else np.nan
    dens_out = (~inside).sum() / width_out if width_out > 0 else np.nan
    enrichment_fold = dens_in / dens_out if (dens_out and dens_out > 0) else np.nan

    return dict(
        n=n, apex_x=a, apex_y=b, pi=pi, enrichment_fold=enrichment_fold,
        se_apex_x=se_a, se_apex_y=se_b,
        apex_x_lo=a - 1.96 * se_a, apex_x_hi=a + 1.96 * se_a,
        apex_y_lo=b - 1.96 * se_b, apex_y_hi=b + 1.96 * se_b,
        loglik_full=ll_full, loglik_null=ll_null,
        LR=LR, df=df, p_chi2=p_chi2,
        aic_full=aic_full, aic_null=aic_null, delta_aic=aic_null - aic_full,
        params=best.x,
    )


# ---------------------------------------------------------------------------
# Permutation calibration (honest p-value)
# ---------------------------------------------------------------------------
def permutation_p(x, y, X, tau, w_floor, obs_LR, n_perm, perm_n, seed,
                  ax_lo=-10.0, ax_hi=10.0, ay_lo=None, ay_hi=None):
    """Shuffle distance (x) to break the V, refit, collect LR null. Empirical p."""
    rng = np.random.default_rng(seed)
    n = len(x)
    exceed = 0
    null_LR = np.empty(n_perm)
    for p in range(n_perm):
        if perm_n > 0 and perm_n < n:
            idx = rng.integers(0, n, perm_n)
            ys = y[idx]
            xs = rng.permutation(x[idx])
            Xn = perm_n
        else:
            ys = y
            xs = rng.permutation(x)
            Xn = n
        fit = fit_apex(xs, ys, X, tau, w_floor, seed=seed + p,
                       ax_lo=ax_lo, ax_hi=ax_hi, ay_lo=ay_lo, ay_hi=ay_hi)
        lr = fit["LR"] if fit is not None else 0.0
        null_LR[p] = lr
        if lr >= obs_LR:
            exceed += 1
    emp_p = (1 + exceed) / (1 + n_perm)
    return emp_p, null_LR


# ---------------------------------------------------------------------------
# Diagnostic plot
# ---------------------------------------------------------------------------
def make_plot(x, y, fit, name, out_png, X):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a, b = fit["apex_x"], fit["apex_y"]
    fig, ax = plt.subplots(figsize=(5, 4.2))
    s = min(len(x), 20000)
    sel = np.random.default_rng(0).choice(len(x), s, replace=False)
    ax.scatter(x[sel], y[sel], s=3, c="#9bb8d3", alpha=0.4, edgecolors="none")
    xs = np.linspace(a, a + X, 50)
    ax.plot(xs, b + 2 * (xs - a), color="#104e8b", lw=2)          # right arm
    ax.plot(2 * a - xs, b + 2 * (xs - a), color="#104e8b", lw=2)  # left arm
    ax.errorbar(a, b, xerr=1.96 * fit["se_apex_x"], yerr=1.96 * fit["se_apex_y"],
                fmt="*", color="red", ms=15, capsize=3, zorder=5)
    ax.set_title(f"{name}\napex=({a:.1f}, {b:.1f})  pi={fit['pi']:.2f}  "
                 f"LR={fit['LR']:.0f}  p_chi2={fit['p_chi2']:.1e}", fontsize=9)
    ax.set_xlabel("distance (bp)"); ax.set_ylabel("fragment length (bp)")
    ax.set_xlim(-X, X)
    plt.tight_layout(); plt.savefig(out_png, dpi=150); plt.close(fig)


# ---------------------------------------------------------------------------
# Per-file driver
# ---------------------------------------------------------------------------
def load_file(path, frag_min, frag_max, X, header, max_n, seed):
    hdr = 0 if header else None
    df = pd.read_csv(path, sep="\t", header=hdr, usecols=[1, 2],
                     names=["fragment", "distance"] if not header else None)
    if header:
        df = df.iloc[:, :2]
        df.columns = ["fragment", "distance"]
    y = pd.to_numeric(df["fragment"], errors="coerce").to_numpy()
    x = pd.to_numeric(df["distance"], errors="coerce").to_numpy()
    m = np.isfinite(x) & np.isfinite(y) & (y >= frag_min) & (y <= frag_max) & (np.abs(x) <= X)
    x, y = x[m], y[m]
    if max_n > 0 and len(x) > max_n:
        idx = np.random.default_rng(seed).choice(len(x), max_n, replace=False)
        x, y = x[idx], y[idx]
    return x, y


def process(path, args):
    name = os.path.splitext(os.path.basename(path))[0]
    try:
        x, y = load_file(path, args.frag_min, args.frag_max, args.x_window,
                         args.header, args.max_n, args.seed)
    except Exception as e:
        return {"motif": name, "status": f"read_error:{e}"}
    if len(x) < args.min_n:
        return {"motif": name, "n": len(x), "status": "too_few_fragments"}

    fit = fit_apex(x, y, args.x_window, args.tau, args.w_floor, args.seed,
                   ax_lo=args.apex_x_lo, ax_hi=args.apex_x_hi,
                   ay_lo=args.apex_y_lo, ay_hi=args.apex_y_hi)
    if fit is None:
        return {"motif": name, "n": len(x), "status": "fit_failed"}

    emp_p = np.nan
    if args.permutations > 0:
        emp_p, _ = permutation_p(x, y, args.x_window, args.tau, args.w_floor,
                                 fit["LR"], args.permutations, args.perm_n,
                                 args.seed + 7,
                                 ax_lo=args.apex_x_lo, ax_hi=args.apex_x_hi,
                                 ay_lo=args.apex_y_lo, ay_hi=args.apex_y_hi)

    if args.plot_dir:
        os.makedirs(args.plot_dir, exist_ok=True)
        try:
            make_plot(x, y, fit, name, os.path.join(args.plot_dir, name + ".png"),
                      args.x_window)
        except Exception as e:
            sys.stderr.write(f"[plot warn] {name}: {e}\n")

    return {
        "motif": name, "n": fit["n"],
        "apex_x": fit["apex_x"], "apex_y_channel_width": fit["apex_y"],
        "se_apex_x": fit["se_apex_x"], "se_apex_y": fit["se_apex_y"],
        "apex_x_lo": fit["apex_x_lo"], "apex_x_hi": fit["apex_x_hi"],
        "apex_y_lo": fit["apex_y_lo"], "apex_y_hi": fit["apex_y_hi"],
        "pi_enrichment": fit["pi"], "enrichment_fold": fit["enrichment_fold"],
        "LR": fit["LR"], "p_chi2": fit["p_chi2"], "p_perm": emp_p,
        "delta_aic": fit["delta_aic"],
        "loglik_full": fit["loglik_full"], "loglik_null": fit["loglik_null"],
        "status": "ok",
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def gather_inputs(inp, pattern):
    if os.path.isdir(inp):
        return sorted(glob.glob(os.path.join(inp, pattern)))
    return [inp]


def main():
    ap = argparse.ArgumentParser(description="Method B: generative V-plot apex finder")
    ap.add_argument("--input", required=True, help="file, or directory (use --glob)")
    ap.add_argument("--glob", default="*_fragL_dist.txt", help="pattern when --input is a dir")
    ap.add_argument("--output", default="vplot_apex_scores.tsv")
    ap.add_argument("--header", action="store_true", help="input has a header row")
    ap.add_argument("--frag-min", type=float, default=20.0)
    ap.add_argument("--frag-max", type=float, default=75.0)
    ap.add_argument("--x-window", type=float, default=150.0, help="half-width X of distance window / background")
    ap.add_argument("--apex-x-lo", type=float, default=-10.0, help="lower bound for apex_x search (motif-centred)")
    ap.add_argument("--apex-x-hi", type=float, default=10.0, help="upper bound for apex_x search (motif-centred)")
    ap.add_argument("--apex-y-lo", type=float, default=None,
                    help="lower bound for apex_y / channel width (default: min_fragment - 10)")
    ap.add_argument("--apex-y-hi", type=float, default=None,
                    help="upper bound for apex_y / channel width (default: max_fragment)")
    ap.add_argument("--max-n", type=int, default=200000, help="subsample cap per file for the fit; 0 = use all")
    ap.add_argument("--min-n", type=int, default=200, help="skip files with fewer fragments")
    ap.add_argument("--tau", type=float, default=1.5, help="edge smoothing (bp) of the cone top-hat")
    ap.add_argument("--w-floor", type=float, default=0.75, help="minimum cone half-width")
    ap.add_argument("--permutations", type=int, default=0, help="permutation reps for calibrated p; 0 = off")
    ap.add_argument("--perm-n", type=int, default=50000, help="subsample size inside each permutation; 0 = all")
    ap.add_argument("--plot-dir", default="", help="write a diagnostic PNG per file here")
    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    files = gather_inputs(args.input, args.glob)
    if not files:
        sys.exit("No input files matched.")
    sys.stderr.write(f"Scoring {len(files)} file(s); objective = generative LRT.\n")

    if args.cores > 1 and len(files) > 1:
        import multiprocessing as mp
        with mp.Pool(min(args.cores, len(files))) as pool:
            rows = pool.starmap(process, [(f, args) for f in files])
    else:
        rows = [process(f, args) for f in files]

    out = pd.DataFrame(rows)
    if "p_perm" in out and out["p_perm"].notna().any():
        out["q_perm"] = _bh(out["p_perm"].values)
    out["q_chi2"] = _bh(out["p_chi2"].values) if "p_chi2" in out else np.nan
    out.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    sys.stderr.write(f"Wrote {args.output}\n")
    # short console summary
    cols = [c for c in ["motif", "n", "apex_x", "apex_y_channel_width",
                        "se_apex_x", "se_apex_y", "pi_enrichment", "LR",
                        "p_chi2", "p_perm", "status"] if c in out]
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        sys.stderr.write(out[cols].to_string(index=False) + "\n")


def _bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    idx = np.where(ok)[0]
    if len(idx) == 0:
        return q
    pp = p[idx]
    order = np.argsort(pp)
    m = len(pp)
    ranked = pp[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qq = np.empty(m); qq[order] = np.clip(ranked, 0, 1)
    q[idx] = qq
    return q


if __name__ == "__main__":
    main()
