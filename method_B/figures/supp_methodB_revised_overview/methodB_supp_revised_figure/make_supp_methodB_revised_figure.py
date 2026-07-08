#!/usr/bin/env python3
"""Build a revised supplementary overview figure for the V-plot model.

The figure separates the V-like-structure test from the calibrated effect-size
rule used for TF-vs-bias classification.
"""
from __future__ import annotations

import argparse
import glob
import importlib.util
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


BLUE = "#1f5a93"
TEAL = "#2a8f7a"
ORANGE = "#d46b32"
RED = "#bd3f3f"
DARK = "#222222"
MID = "#5b6670"
PALE_BLUE = "#eaf2f8"
PALE_GREEN = "#e9f4ef"
PALE_ORANGE = "#f7ede5"
GRID = "#d7dde3"


def find_first_existing_file(patterns: list[str | Path]) -> Path | None:
    """Return the first existing path, allowing glob patterns."""
    for pattern in patterns:
        text = str(pattern)
        matches = sorted(Path(p) for p in glob.glob(text))
        for match in matches:
            if match.exists():
                return match
        path = Path(text)
        if path.exists():
            return path
    return None


def read_table_auto(path: Path) -> pd.DataFrame:
    """Read CSV/TSV-like files with a small amount of format inference."""
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def _norm_col(name: str) -> str:
    return "".join(ch.lower() for ch in str(name) if ch.isalnum())


def infer_column(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    """Find a column by case/punctuation-insensitive candidate names."""
    lookup = {_norm_col(c): c for c in df.columns}
    for candidate in candidates:
        key = _norm_col(candidate)
        if key in lookup:
            return lookup[key]
    if required:
        raise ValueError(
            "Missing required column. Tried: "
            + ", ".join(candidates)
            + ". Available: "
            + ", ".join(map(str, df.columns))
        )
    return None


def clean_name(value: str) -> str:
    name = Path(str(value)).stem
    for suffix in [
        "_ATAC_fragL_dist",
        "_DNase_fragL_dist",
        "_loMNase_fragL_dist",
        "_fragL_dist",
    ]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def load_apex_summary(input_root: Path) -> pd.DataFrame:
    tf_path = find_first_existing_file(
        [
            input_root / "ATAC_TF_apex.tsv",
            input_root / "*_TF_apex.tsv",
            input_root / "*TF_apex_scores.tsv",
        ]
    )
    bias_path = find_first_existing_file(
        [
            input_root / "ATAC_bias_apex.tsv",
            input_root / "*_bias_apex.tsv",
            input_root / "*Bias*apex_scores.tsv",
        ]
    )
    if tf_path is None or bias_path is None:
        raise FileNotFoundError("Could not find both TF and bias apex summary tables.")

    frames: list[pd.DataFrame] = []
    for path, group in [(tf_path, "TF motifs"), (bias_path, "matched enzyme-bias controls")]:
        df = read_table_auto(path)
        motif_col = infer_column(df, ["motif", "name", "label", "id"])
        width_col = infer_column(
            df,
            [
                "apex_y_channel_width",
                "width",
                "v_width",
                "V_width",
                "apex_y",
                "b",
                "channel_width",
                "v_channel_width",
                "model_fitted_width",
            ],
        )
        fold_col = infer_column(df, ["enrichment_fold", "V_in_V_out", "Vin_Vout"], required=False)
        enrich_col = infer_column(
            df,
            [
                "enrichment",
                "inside_v_enrichment",
                "insideV_enrichment",
                "log2_Vin_Vout",
                "log2_v_in_v_out",
                "log2_enrichment",
                "v_enrichment",
                "enrich",
            ],
            required=False,
        )
        p_col = infer_column(df, ["p_chi2", "lrt_p", "LRT_p", "p_lrt", "lrt_pvalue", "pvalue", "p_value"], required=False)
        lr_col = infer_column(df, ["LR", "lrt", "likelihood_ratio"], required=False)
        ax_col = infer_column(df, ["apex_x", "a"], required=False)
        pi_col = infer_column(df, ["pi_enrichment", "pi", "signal_fraction"], required=False)
        status_col = infer_column(df, ["status"], required=False)

        out = pd.DataFrame(
            {
                "motif": df[motif_col].astype(str),
                "name": df[motif_col].astype(str).map(clean_name),
                "width": pd.to_numeric(df[width_col], errors="coerce"),
                "group": group,
                "source_table": str(path),
            }
        )
        if fold_col is not None:
            fold = pd.to_numeric(df[fold_col], errors="coerce")
            out["enrichment"] = np.log2(fold)
            out["enrichment_fold"] = fold
        elif enrich_col is not None:
            out["enrichment"] = pd.to_numeric(df[enrich_col], errors="coerce")
            out["enrichment_fold"] = np.nan
        else:
            raise ValueError(f"Missing inside-V enrichment field in {path}")

        out["lrt_p"] = pd.to_numeric(df[p_col], errors="coerce") if p_col else np.nan
        out["LR"] = pd.to_numeric(df[lr_col], errors="coerce") if lr_col else np.nan
        out["apex_x"] = pd.to_numeric(df[ax_col], errors="coerce") if ax_col else 0.0
        out["pi"] = pd.to_numeric(df[pi_col], errors="coerce") if pi_col else np.nan
        out["status"] = df[status_col].astype(str) if status_col else "unknown"
        frames.append(out)

    combined = pd.concat(frames, ignore_index=True)
    combined = combined[np.isfinite(combined["width"]) & np.isfinite(combined["enrichment"])]
    return combined


def load_cutoffs(input_root: Path) -> tuple[float, float, pd.DataFrame, Path]:
    cutoff_path = find_first_existing_file(
        [
            input_root / "ATAC_methodB_cutoffs.csv",
            input_root / "*methodB_cutoffs.csv",
            input_root / "*cutoffs.csv",
        ]
    )
    if cutoff_path is None:
        raise FileNotFoundError("Could not find a cutoff table.")
    table = read_table_auto(cutoff_path)
    feature_col = infer_column(table, ["feature"])
    cutoff_col = infer_column(table, ["cutoff", "threshold"])

    width_row = table[table[feature_col].astype(str).str.contains("width", case=False, na=False)]
    enrich_row = table[
        table[feature_col].astype(str).str.contains("log2|enrich|V-in", case=False, regex=True, na=False)
    ]
    if width_row.empty or enrich_row.empty:
        raise ValueError(f"Could not identify width/enrichment cutoff rows in {cutoff_path}")

    width_cutoff = float(pd.to_numeric(width_row.iloc[0][cutoff_col], errors="coerce"))
    enrich_cutoff = float(pd.to_numeric(enrich_row.iloc[0][cutoff_col], errors="coerce"))
    return width_cutoff, enrich_cutoff, table, cutoff_path


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.120,
        1.010,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=DARK,
    )


def style_mini_vplot(ax: plt.Axes) -> None:
    ax.set_xlim(-90, 90)
    ax.set_ylim(20, 150)
    ax.set_xticks([-80, 0, 80])
    ax.set_yticks([40, 90, 140])
    ax.tick_params(labelsize=6.8, length=2.5, pad=1)
    ax.grid(True, color=GRID, lw=0.45, alpha=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
        spine.set_color("#6c7882")


def draw_schematic_vplot(ax: plt.Axes, b: float, strength: float, color: str, seed: int) -> None:
    rng = np.random.default_rng(seed)
    n_bg = 420
    x_bg = rng.uniform(-85, 85, n_bg)
    y_bg = rng.uniform(25, 145, n_bg)
    ax.scatter(x_bg, y_bg, s=3, c="#c9d1d8", alpha=0.35, edgecolors="none")

    n_sig = int(260 * strength)
    y_sig = rng.uniform(max(25, b + 6), 145, n_sig)
    half = np.maximum((y_sig - b) / 2, 2)
    x_sig = rng.normal(0, half / 2.8)
    x_sig = np.clip(x_sig, -half, half)
    ax.scatter(x_sig, y_sig, s=4, c=color, alpha=0.55, edgecolors="none")

    xs = np.linspace(0, 80, 100)
    y = b + 2 * xs
    keep = y <= 150
    ax.plot(xs[keep], y[keep], color=color, lw=1.6)
    ax.plot(-xs[keep], y[keep], color=color, lw=1.6)
    ax.scatter([0], [b], marker="*", s=38, color=RED, zorder=5)
    style_mini_vplot(ax)


_FIT_CORE_CACHE = {}


def load_fit_core(method_root: Path):
    core_path = method_root / "fit_vplot_apex.py"
    key = str(core_path.resolve())
    if key not in _FIT_CORE_CACHE:
        spec = importlib.util.spec_from_file_location("methodb_fit_vplot_apex_for_panel", core_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not import Method B fitter from {core_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _FIT_CORE_CACHE[key] = module
    return _FIT_CORE_CACHE[key]


def _logit(value: float) -> float:
    value = float(np.clip(value, 1e-6, 1.0 - 1e-6))
    return math.log(value / (1.0 - value))


def _read_profile_start(score_path: Path, motif: str) -> tuple[float, float]:
    if not score_path.exists():
        return 10.0, 0.0
    table = pd.read_csv(score_path, sep="\t")
    row = table[table["motif"].astype(str) == motif]
    if row.empty:
        row = table[table["motif"].astype(str).str.contains(motif, regex=False, na=False)]
    if row.empty:
        return 10.0, 0.0
    first = row.iloc[0]
    b = float(pd.to_numeric(first.get("apex_y_channel_width", 10.0), errors="coerce"))
    pi = float(pd.to_numeric(first.get("pi_enrichment", 0.5), errors="coerce"))
    if not np.isfinite(b):
        b = 10.0
    if not np.isfinite(pi):
        pi = 0.5
    return float(np.clip(b, 0.0, 60.0)), _logit(pi)


def compute_apex_profile(
    fit_core,
    path: Path,
    label: str,
    b0: float,
    lpi0: float,
    seed: int,
    grid: np.ndarray,
    max_n: int = 200000,
) -> pd.DataFrame:
    x, y = fit_core.load_file(
        str(path),
        frag_min=20.0,
        frag_max=75.0,
        X=150.0,
        header=False,
        max_n=max_n,
        seed=seed,
    )
    rows = []
    for candidate_a in grid:
        params = np.array([candidate_a, b0, lpi0], dtype=float)
        loglik = float(-fit_core.neg_loglik(params, x, y, 150.0, 1.5, 0.75))
        rows.append(
            {
                "label": label,
                "candidate_apex_x": float(candidate_a),
                "loglik": loglik,
                "fixed_b": float(b0),
                "fixed_pi": float(1.0 / (1.0 + math.exp(-lpi0))),
                "n": int(len(x)),
                "source_file": str(path),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["relative_loglik"] = out["loglik"] - out["loglik"].max()
    return out


def load_or_compute_real_apex_profiles(method_root: Path) -> pd.DataFrame | None:
    profile_path = method_root / "methodB_supp_revised_figure" / "real_loMNase_apex_likelihood_profiles.csv"
    if profile_path.exists():
        return pd.read_csv(profile_path)

    ctcf_path = find_first_existing_file(
        [
            method_root / "TF_vplot" / "from_server_three_assays" / "TF_fragL_dist_three_assays" / "loMNase" / "CTCF_fragL_dist.txt",
            method_root / "TF_vplot" / "CTCF_loMNase_fragL_dist.txt",
            method_root / "TF_vplot" / "CTCF_fragL_dist.txt",
        ]
    )
    bias_path = method_root / "sequence_bias" / "MNase_Bias_AAG_loMNase_fragL_dist.txt"
    if ctcf_path is None or not bias_path.exists():
        return None

    fit_core = load_fit_core(method_root)
    grid = np.linspace(-10.0, 10.0, 81)
    ctcf_b, ctcf_lpi = _read_profile_start(method_root / "loMNase_TF_apex_scores.tsv", "CTCF_loMNase_fragL_dist")
    bias_b, bias_lpi = _read_profile_start(
        method_root / "sequence_bias" / "loMNase_Bias_3bp_apex_scores.tsv",
        "MNase_Bias_AAG_loMNase_fragL_dist",
    )
    frames = [
        compute_apex_profile(fit_core, ctcf_path, "CTCF motif", ctcf_b, ctcf_lpi, seed=1, grid=grid),
        compute_apex_profile(fit_core, bias_path, "bias AAG", bias_b, bias_lpi, seed=2, grid=grid),
    ]
    profile = pd.concat([frame for frame in frames if frame is not None and not frame.empty], ignore_index=True)
    if profile.empty:
        return None
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    profile.to_csv(profile_path, index=False)
    return profile


def draw_vplot_model_mle(ax: plt.Axes, method_root: Path | None = None) -> None:
    """Combined V-plot model schematic and MLE profile panel."""
    ax.set_axis_off()
    ax.set_title("V-plot model and MLE-based apex estimation", loc="left", fontsize=11.0, fontweight="bold", pad=7)

    vax = ax.inset_axes([0.055, 0.15, 0.42, 0.66])
    rng = np.random.default_rng(42)
    b = 20.0
    n_bg = 900
    x_bg = rng.uniform(-60, 60, n_bg)
    y_bg = rng.uniform(0, 100, n_bg)
    vax.scatter(x_bg, y_bg, s=4, c="#bfc5cc", alpha=0.55, edgecolors="none")

    n_sig = 520
    y_sig = rng.uniform(b + 2, 100, n_sig)
    half = np.maximum((y_sig - b) * 0.5, 1.0)
    x_sig = rng.normal(0, half / 2.6)
    x_sig = np.clip(x_sig, -half, half)
    vax.scatter(x_sig, y_sig, s=5, c=BLUE, alpha=0.72, edgecolors="none")

    xs = np.linspace(0, 50, 100)
    yy = b + 2 * xs
    keep = yy <= 100
    vax.plot(xs[keep], yy[keep], color=BLUE, lw=2.0)
    vax.plot(-xs[keep], yy[keep], color=BLUE, lw=2.0)
    vax.scatter([0], [b], marker="*", s=115, color=RED, edgecolor="white", linewidth=0.5, zorder=5)
    vax.annotate(
        "apex (a, b)",
        xy=(0, b),
        xytext=(22, 10),
        arrowprops=dict(arrowstyle="->", lw=1.2, color=RED),
        fontsize=8.4,
        color=RED,
        fontweight="bold",
    )
    vax.text(-55, 30, "background\n(1 - pi)", fontsize=8.5, color=MID, ha="left")
    vax.text(0, 63, "inside V\nsignal fraction pi", fontsize=9.2, color=BLUE, ha="center", fontweight="bold")
    vax.text(-58, 19, "b = V-channel width", fontsize=8.2, color=MID, ha="left")
    vax.text(11, 79, "w(y) = (y - b) / 2", fontsize=8.0, color=DARK)
    vax.set_xlim(-60, 60)
    vax.set_ylim(0, 100)
    vax.set_aspect("equal", adjustable="box")
    vax.set_xlabel("distance from motif (bp)", fontsize=8.5, labelpad=2)
    vax.set_ylabel("fragment length (bp)", fontsize=8.5, labelpad=2)
    vax.set_xticks([-60, -30, 0, 30, 60])
    vax.set_yticks([0, 20, 40, 60, 80, 100])
    vax.tick_params(labelsize=7.5, length=2.6, pad=1)
    vax.grid(True, color=GRID, lw=0.4, alpha=0.55)
    for spine in vax.spines.values():
        spine.set_linewidth(0.85)
        spine.set_color("#4c5560")

    rax = ax.inset_axes([0.565, 0.15, 0.39, 0.66])
    profile = load_or_compute_real_apex_profiles(method_root) if method_root is not None else None
    if profile is not None and not profile.empty:
        rax.set_title("Real loMNase likelihood scan", fontsize=8.0, fontweight="bold", pad=5)
        for label, color, lw in [("CTCF motif", BLUE, 1.9), ("bias AAG", ORANGE, 1.7)]:
            sub = profile[profile["label"] == label].sort_values("candidate_apex_x")
            if sub.empty:
                continue
            y_values = sub["relative_loglik"] / sub["n"]
            rax.plot(sub["candidate_apex_x"], y_values, color=color, lw=lw, label=label)
            if label == "CTCF motif":
                best = sub.loc[sub["relative_loglik"].idxmax()]
                rax.axvline(float(best["candidate_apex_x"]), color=RED, lw=1.0, ls="--")
        note_text = "CTCF MLE"
        y_label = "mean ln L(a) -\nmax mean ln L"
    else:
        x = np.linspace(-20, 20, 140)
        tf_like = -0.028 * (x - 1.5) ** 2 + 1.05
        weak = 0.09 * np.sin(x / 5.2) - 0.0025 * x**2 - 0.10
        tf_like = tf_like - np.max(tf_like)
        weak = weak - np.max(weak)
        rax.set_title("Schematic apex likelihood profile", fontsize=8.0, fontweight="bold", pad=5)
        rax.plot(x, tf_like, color=BLUE, lw=1.9, label="clear apex")
        rax.plot(x, weak, color=ORANGE, lw=1.7, label="weak/flat")
        rax.axvline(1.5, color=RED, lw=1.0, ls="--")
        note_text = "selected apex =\nmaximum likelihood"
        y_label = "ln L(a) -\nmax ln L"
    rax.text(
        0.96,
        0.90,
        note_text,
        transform=rax.transAxes,
        ha="right",
        va="top",
        fontsize=7.1,
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.80),
    )
    rax.set_xlabel("candidate apex position a (bp)", fontsize=7.5, labelpad=1)
    rax.set_ylabel(y_label, fontsize=7.2, labelpad=1)
    rax.tick_params(labelsize=7.0, length=2.4, pad=1)
    rax.grid(True, color=GRID, lw=0.4, alpha=0.75)
    rax.legend(frameon=False, fontsize=7.0, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.12))


def load_best_assay_data(method_root: Path, assay: str = "ATAC") -> tuple[pd.DataFrame, pd.Series, Path, Path]:
    data_dir = method_root / "v_apex_scores" / "methodB_scatter_best_assay_specific"
    data_path = data_dir / f"plot_data_best_{assay}.csv"
    threshold_path = method_root / "fig_cutoff_derivation_thresholds.csv"
    if not data_path.exists():
        raise FileNotFoundError(f"Missing final scatter data: {data_path}")
    if not threshold_path.exists():
        raise FileNotFoundError(f"Missing cut-off derivation table: {threshold_path}")
    df = pd.read_csv(data_path)
    thr = pd.read_csv(threshold_path)
    assay_thr = thr[thr["assay"].astype(str).str.lower() == assay.lower()]
    width_row = assay_thr[assay_thr["feature"].astype(str).str.contains("width", case=False, na=False)].iloc[0]
    enrich_row = assay_thr[assay_thr["feature"].astype(str).str.contains("enrichment", case=False, na=False)].iloc[0]
    cutoffs = pd.Series({"width": float(width_row["cutoff"]), "enrichment": float(enrich_row["cutoff"])})
    return df, cutoffs, data_path, threshold_path


def draw_two_feature_scatter(ax: plt.Axes, method_root: Path, assay: str = "ATAC") -> tuple[Path, Path]:
    df, cutoffs, data_path, threshold_path = load_best_assay_data(method_root, assay)
    ax.set_title("Two-dimensional effect-size rule", loc="left", fontsize=10.3, fontweight="bold", pad=7)
    bias = df[df["type"].astype(str).str.lower() == "bias"].copy()
    tf = df[df["type"].astype(str).str.lower() != "bias"].copy()
    xmax = max(df["width"].max() * 1.08, cutoffs["width"] + 8)
    ymax = max(df["log2FC"].max() * 1.12, cutoffs["enrichment"] + 0.6)
    ax.axvspan(cutoffs["width"], xmax, color=PALE_GREEN, alpha=0.48, zorder=0)
    ax.axhspan(cutoffs["enrichment"], ymax, color=PALE_GREEN, alpha=0.22, zorder=0)
    ax.scatter(
        bias["width"],
        bias["log2FC"],
        s=42,
        facecolor="white",
        edgecolor="#ff6a6a",
        lw=1.3,
        label="enzyme bias",
        zorder=3,
    )
    ax.scatter(tf["width"], tf["log2FC"], s=46, facecolor=BLUE, edgecolor="white", lw=0.55, label="TF motif", zorder=4)
    ax.axvline(cutoffs["width"], color=DARK, ls="--", lw=1.1)
    ax.axhline(cutoffs["enrichment"], color=DARK, ls="--", lw=1.1)
    ax.set_xlim(-1.5, xmax)
    ax.set_ylim(-0.02, ymax)
    ax.set_xlabel("V-channel width (bp)", fontsize=8.4, labelpad=2)
    ax.set_ylabel("inside-V enrichment, log2(V-in/V-out)", fontsize=8.4, labelpad=2)
    ax.tick_params(labelsize=7.2, length=2.6, pad=1)
    ax.grid(True, color=GRID, lw=0.45, alpha=0.80)
    ax.legend(frameon=False, fontsize=7.4, loc="upper left", borderaxespad=0.2)
    ax.text(
        0.98,
        0.88,
        f"{assay} classified as TF-footprint-like\nonly above both cut-offs",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        color="#1e5c40",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#b9d7c8", alpha=0.92),
    )
    ax.text(
        0.055,
        0.18,
        "bias-like /\nweak V-like",
        transform=ax.transAxes,
        ha="left",
        fontsize=7.5,
        color=ORANGE,
        fontweight="bold",
    )
    ax.text(
        0.05,
        0.06,
        "below either cut-off",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.8,
        color=ORANGE,
    )
    return data_path, threshold_path


def draw_cutoff_derivation_panel(ax: plt.Axes, method_root: Path) -> tuple[Path, Path, Path]:
    ax.set_axis_off()
    ax.set_title("Cut-offs are derived by the mid-gap rule", loc="left", fontsize=10.8, fontweight="bold", pad=7)
    data_dir = method_root / "v_apex_scores" / "methodB_scatter_best_assay_specific"
    threshold_path = method_root / "fig_cutoff_derivation_thresholds.csv"
    stats_path = data_dir / "methodB_separation_stats_best_common.csv"
    thresholds = pd.read_csv(threshold_path)
    stats = pd.read_csv(stats_path) if stats_path.exists() else pd.DataFrame()
    assays = ["loMNase", "DNase", "ATAC"]
    lane_y = {"loMNase": 2.35, "DNase": 1.35, "ATAC": 0.35}
    lane_half_height = 0.22
    rng = np.random.default_rng(1)

    enrich_ax = ax.inset_axes([0.04, 0.20, 0.285, 0.68])
    width_ax = ax.inset_axes([0.405, 0.20, 0.315, 0.68])
    table_ax = ax.inset_axes([0.765, 0.20, 0.225, 0.68])
    table_ax.set_axis_off()

    assay_data = {}
    for assay in assays:
        path = data_dir / f"plot_data_best_{assay}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Missing cut-off plot data: {path}")
        assay_data[assay] = pd.read_csv(path)

    def threshold_row(assay: str, feature_key: str) -> pd.Series:
        sub = thresholds[
            (thresholds["assay"].astype(str) == assay)
            & thresholds["feature"].astype(str).str.contains(feature_key, case=False, na=False)
        ]
        if sub.empty:
            raise ValueError(f"Missing {feature_key} threshold for {assay}")
        return sub.iloc[0]

    def draw_metric(metric_ax: plt.Axes, metric_col: str, feature_key: str, title: str, xlabel: str, unit: str) -> None:
        for assay in assays:
            yc = lane_y[assay]
            df = assay_data[assay]
            bias = pd.to_numeric(df.loc[df["type"].astype(str).str.lower() == "bias", metric_col], errors="coerce")
            tf = pd.to_numeric(df.loc[df["type"].astype(str).str.lower() != "bias", metric_col], errors="coerce")
            thr_row = threshold_row(assay, feature_key)
            max_bias = float(thr_row["max_bias"])
            min_tf = float(thr_row["min_TF"])
            cutoff = float(thr_row["cutoff"])
            margin = float(thr_row["margin"])
            metric_ax.add_patch(
                Rectangle(
                    (max_bias, yc - lane_half_height),
                    min_tf - max_bias,
                    2 * lane_half_height,
                    facecolor=PALE_BLUE,
                    edgecolor="none",
                    alpha=0.72,
                    zorder=0,
                )
            )
            metric_ax.scatter(
                bias,
                yc + rng.uniform(-lane_half_height * 0.50, lane_half_height * 0.50, len(bias)),
                s=28,
                facecolors="none",
                edgecolors="#ff6a6a",
                linewidths=1.2,
                zorder=3,
            )
            metric_ax.scatter(
                tf,
                yc + rng.uniform(-lane_half_height * 0.50, lane_half_height * 0.50, len(tf)),
                s=26,
                c=BLUE,
                edgecolors="none",
                zorder=3,
            )
            metric_ax.plot([cutoff, cutoff], [yc - lane_half_height - 0.05, yc + lane_half_height + 0.05], color=DARK, ls="--", lw=1.2)
            metric_ax.text(
                cutoff,
                yc + lane_half_height + 0.08,
                f"cut-off={cutoff:.2f}{unit}",
                ha="center",
                va="bottom",
                fontsize=6.6,
                fontweight="bold",
                color=DARK,
            )
            metric_ax.add_patch(
                FancyArrowPatch(
                    (cutoff, yc - lane_half_height - 0.02),
                    (min_tf, yc - lane_half_height - 0.02),
                    arrowstyle="<->",
                    mutation_scale=8.0,
                    lw=0.9,
                    color=MID,
                    zorder=4,
                )
            )
            metric_ax.text((cutoff + min_tf) / 2, yc - lane_half_height - 0.12, f"margin {margin:.2f}", ha="center", va="top", fontsize=6.3, color=MID)

        metric_ax.set_yticks([lane_y[a] for a in assays])
        metric_ax.set_yticklabels(assays, fontsize=7.2)
        metric_ax.set_ylim(-0.08, 2.85)
        metric_ax.set_xlabel(xlabel, fontsize=7.4, labelpad=1)
        metric_ax.set_title(title, fontsize=8.6, fontweight="bold", loc="left", pad=4)
        metric_ax.tick_params(axis="x", labelsize=7.0, length=2.4, pad=1)
        metric_ax.grid(axis="x", color=GRID, lw=0.4, alpha=0.70)
        metric_ax.margins(x=0.08)
        for spine in metric_ax.spines.values():
            spine.set_linewidth(0.8)
            spine.set_color("#4c5560")

    draw_metric(enrich_ax, "log2FC", "enrichment", "Inside-V enrichment cut-off", "log2(V-in / V-out)", "")
    enrich_ax.set_xlim(-0.25, max(3.6, enrich_ax.get_xlim()[1]))
    draw_metric(width_ax, "width", "width", "V-channel width cut-off", "V-channel width (bp)", " bp")
    width_ax.set_xlim(-2.5, max(62.0, width_ax.get_xlim()[1]))

    legend_handles = [
        plt.Line2D([], [], marker="o", color=BLUE, ls="none", ms=5.5, label="TF motif"),
        plt.Line2D([], [], marker="o", mfc="white", mec="#ff6a6a", mew=1.2, color="#ff6a6a", ls="none", ms=5.5, label="enzyme bias"),
    ]
    enrich_ax.legend(handles=legend_handles, frameon=False, fontsize=6.9, loc="lower right", borderpad=0.1, handletextpad=0.35)

    rows = []
    for assay in assays:
        width = float(threshold_row(assay, "width")["cutoff"])
        enrich = float(threshold_row(assay, "enrichment")["cutoff"])
        sens = spec = 1.0
        if not stats.empty:
            rule = stats[(stats["assay"].astype(str) == assay) & stats["metric"].astype(str).str.contains("2D", case=False, na=False)]
            if not rule.empty:
                sens = float(rule.iloc[0].get("sens", 1.0))
                spec = float(rule.iloc[0].get("spec", 1.0))
        rows.append([assay, f"{width:.2f}", f"{enrich:.2f}", f"{sens * 100:.0f}%"])
    table = table_ax.table(
        cellText=rows,
        colLabels=["assay", "V-channel\nwidth", "Inside-V\nenrichment", "rule"],
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=[0.25, 0.27, 0.34, 0.16],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.9)
    table.scale(1.06, 2.05)
    for (r, _c), cell in table.get_celld().items():
        cell.set_linewidth(0.45)
        cell.set_edgecolor("#8b959e")
        if r == 0:
            cell.set_facecolor("#edf2f6")
            cell.set_text_props(weight="bold")
    ax.text(
        0.50,
        0.07,
        "cut-off = midpoint between the highest bias value and the lowest TF value; margin = half the gap",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.2,
        color=DARK,
    )
    return data_dir, threshold_path, stats_path


def draw_problem_schematic(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title(
        "TF footprints and enzyme bias both create V-like patterns",
        loc="left",
        fontsize=9.6,
        fontweight="bold",
        pad=7,
    )
    left = ax.inset_axes([0.02, 0.22, 0.44, 0.66])
    right = ax.inset_axes([0.54, 0.22, 0.44, 0.66])
    draw_schematic_vplot(left, b=48, strength=1.15, color=BLUE, seed=7)
    draw_schematic_vplot(right, b=8, strength=0.55, color=ORANGE, seed=13)
    left.set_title("TF-footprint-like\nwide V-channel\nstrong inside-V enrichment", fontsize=7.3)
    right.set_title("enzyme-bias-like\nnarrow V-like structure\nweak enrichment", fontsize=7.3)
    left.set_xlabel("midpoint relative to motif center", fontsize=6.9, labelpad=1)
    left.set_ylabel("fragment length", fontsize=6.9, labelpad=1)
    right.set_xlabel("midpoint relative to anchor", fontsize=6.9, labelpad=1)
    right.set_yticklabels([])
    ax.text(
        0.5,
        0.055,
        "V-like shape alone is not sufficient.",
        ha="center",
        va="center",
        fontsize=9.5,
        fontweight="bold",
        color=DARK,
    )


def draw_mle_apex_schematic(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title(
        "V-cone model and MLE apex estimation",
        loc="left",
        fontsize=9.6,
        fontweight="bold",
        pad=7,
    )

    geom = ax.inset_axes([0.02, 0.20, 0.31, 0.66])
    geom.set_xlim(-80, 80)
    geom.set_ylim(0, 125)
    geom.set_xticks([])
    geom.set_yticks([])
    for spine in geom.spines.values():
        spine.set_color("#6c7882")
        spine.set_linewidth(0.8)
    a, b = 0, 34
    poly = Polygon([[a, b], [-55, 125], [55, 125]], closed=True, facecolor=PALE_BLUE, edgecolor=BLUE, lw=1.5)
    geom.add_patch(poly)
    geom.scatter([a], [b], marker="*", s=70, color=RED, zorder=4)
    geom.plot([a, 32], [80, 80], color=DARK, lw=1.0)
    geom.plot([a, -32], [80, 80], color=DARK, lw=1.0)
    geom.text(a + 4, b + 3, "apex=(a,b)", fontsize=7.6)
    geom.text(9, 84, "w(y)=(y-b)/2", fontsize=7.4)
    geom.text(-72, 9, "b: model-fitted\nV-channel width", fontsize=7.4)

    ax.text(
        0.38,
        0.78,
        r"$p(x\mid y)=\pi f_{\rm signal}(x\mid y,a,b)+(1-\pi)f_{\rm background}(x)$",
        transform=ax.transAxes,
        fontsize=8.2,
        color=DARK,
    )
    ax.text(
        0.50,
        0.63,
        "pi: fitted signal fraction\ninside-V enrichment is a separate effect-size summary",
        transform=ax.transAxes,
        fontsize=7.4,
        color=MID,
    )

    prof = ax.inset_axes([0.42, 0.18, 0.54, 0.34])
    x = np.linspace(-20, 20, 120)
    tf_like = -0.030 * (x - 1.5) ** 2 + 1.2
    bias_like = 0.12 * np.sin(x / 4.5) - 0.003 * x**2 + 0.08
    prof.plot(x, tf_like, color=BLUE, lw=1.8, label="TF-like")
    prof.plot(x, bias_like, color=ORANGE, lw=1.6, label="bias-like")
    prof.axvline(1.5, color=RED, lw=1.0, ls="--")
    prof.set_xlabel("candidate apex position", fontsize=7.0, labelpad=1)
    prof.set_ylabel("profile\nlog-lik.", fontsize=6.6, labelpad=1)
    prof.tick_params(labelsize=6.6, length=2.5, pad=1)
    prof.grid(True, color=GRID, lw=0.45)
    prof.legend(frameon=False, fontsize=6.8, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.01))
    prof.text(
        0.38,
        0.81,
        "MLE apex = argmax log-likelihood",
        transform=prof.transAxes,
        fontsize=6.8,
        bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="none", alpha=0.75),
    )

    ax.text(
        0.50,
        0.055,
        "Apex is estimated by MLE, not manually selected.",
        transform=ax.transAxes,
        ha="center",
        fontsize=9.1,
        fontweight="bold",
        color=DARK,
    )


def draw_workflow(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Reproducible V-plot workflow", loc="left", fontsize=10.6, fontweight="bold", pad=7)
    boxes = [
        ("Inputs", "fragments\nTF sites\nbias controls"),
        ("V-plot table", "midpoint\nlength\nanchor distance"),
        ("Model fit", "apex\nwidth\ninside-V\nLRT"),
        ("Cut-offs", "mid-gap\ncalibration"),
        ("Classification", "width +\nenrichment\npass cut-offs"),
    ]
    x_positions = [0.030, 0.225, 0.420, 0.615, 0.810]
    box_w = 0.135
    box_h = 0.47
    y0 = 0.36
    for i, ((title, desc), x) in enumerate(zip(boxes, x_positions)):
        face = PALE_GREEN if i == 4 else (PALE_ORANGE if i == 3 else PALE_BLUE)
        box = FancyBboxPatch(
            (x, y0),
            box_w,
            box_h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=0.9,
            edgecolor="#6c7882",
            facecolor=face,
        )
        ax.add_patch(box)
        ax.text(
            x + box_w / 2,
            y0 + box_h - 0.075,
            title,
            ha="center",
            va="center",
            fontsize=7.9,
            fontweight="bold",
            color=DARK,
        )
        ax.text(
            x + box_w / 2,
            y0 + 0.175,
            desc,
            ha="center",
            va="center",
            fontsize=6.45,
            color=DARK,
            linespacing=1.18,
        )
        if i < len(boxes) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + box_w + 0.012, y0 + box_h / 2),
                    (x_positions[i + 1] - 0.012, y0 + box_h / 2),
                    arrowstyle="-|>",
                    mutation_scale=10,
                    lw=0.9,
                    color=MID,
                )
            )
    ax.text(
        0.5,
        0.15,
        "LRT reports V-like evidence; final calls use calibrated width + enrichment cut-offs.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.0,
        color=DARK,
    )


def find_distance_file(row: pd.Series, method_root: Path, input_root: Path) -> Path | None:
    name = str(row["name"])
    motif = str(row["motif"])
    group = str(row["group"])
    if group == "TF motifs":
        tf = clean_name(motif)
        return find_first_existing_file(
            [
                method_root / "TF_vplot" / f"{tf}_ATAC_fragL_dist.txt",
                method_root / "TF_vplot" / "from_server_three_assays" / "TF_fragL_dist_three_assays" / "ATAC" / f"{tf}_fragL_dist.txt",
                method_root / "TF_vplot" / f"{tf}_fragL_dist.txt",
            ]
        )

    kmer = name.replace("ATAC_Bias_", "")
    return find_first_existing_file(
        [
            input_root / "intermediate" / "**" / f"*{kmer}*fragL_dist*.txt",
            method_root / "sequence_bias" / "ATAC_scatter_n1000" / f"ATAC_Bias_{kmer}_ATAC_fragL_dist_n1000.txt",
            method_root / "sequence_bias" / "ATAC_scatter_n300" / f"ATAC_Bias_{kmer}_ATAC_fragL_dist_n300.txt",
            method_root / "sequence_bias" / f"ATAC_Bias_{kmer}_loMNase_fragL_dist.txt",
        ]
    )


def select_representative_rows(
    fit_summary: pd.DataFrame,
    width_cutoff: float,
    enrich_cutoff: float,
    method_root: Path,
    input_root: Path,
    use_real_vplots: str,
) -> list[dict[str, object]]:
    if use_real_vplots == "no":
        return []

    df = fit_summary.copy()
    df["passes"] = (df["width"] >= width_cutoff) & (df["enrichment"] >= enrich_cutoff)
    df["distance_file"] = [find_distance_file(row, method_root, input_root) for _, row in df.iterrows()]
    with_files = df[df["distance_file"].notna()].copy()
    examples: list[pd.Series] = []

    tf = with_files[(with_files["group"] == "TF motifs") & with_files["passes"]].copy()
    if not tf.empty:
        tf["score"] = tf["width"].rank(pct=True) + tf["enrichment"].rank(pct=True)
        examples.append(tf.sort_values("score", ascending=False).iloc[0])

    bias = with_files[(with_files["group"] == "matched enzyme-bias controls") & (~with_files["passes"])].copy()
    if not bias.empty:
        sig = bias[bias["lrt_p"].fillna(1) <= 0.05].copy()
        target = sig if not sig.empty else bias
        target["score"] = target["LR"].fillna(0)
        examples.append(target.sort_values("score", ascending=False).iloc[0])

    if use_real_vplots == "yes" and len(examples) < 2:
        missing = "TF and/or bias real V-plot distance distributions"
        raise FileNotFoundError(f"Requested --use-real-vplots yes, but missing {missing}.")

    out: list[dict[str, object]] = []
    for row in examples[:2]:
        item = row.to_dict()
        item["distance_file"] = Path(str(item["distance_file"]))
        out.append(item)
    return out


def read_vplot_points(path: Path, max_points: int = 55000, seed: int = 3) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, sep="\t", header=None, usecols=[1, 2], names=["fragment", "distance"])
    y = pd.to_numeric(df["fragment"], errors="coerce").to_numpy()
    x = pd.to_numeric(df["distance"], errors="coerce").to_numpy()
    mask = np.isfinite(x) & np.isfinite(y) & (y >= 20) & (y <= 150) & (np.abs(x) <= 150)
    x = x[mask]
    y = y[mask]
    if len(x) > max_points:
        idx = np.random.default_rng(seed).choice(len(x), max_points, replace=False)
        x = x[idx]
        y = y[idx]
    return x, y


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value <= 0:
        return "<1e-300"
    if value < 1e-3:
        return f"{value:.1e}"
    return f"{value:.3f}"


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    """Small dependency-free Markdown table writer."""
    cols = [str(c) for c in df.columns]
    rows = [[str(v) if pd.notna(v) else "" for v in row] for row in df.to_numpy()]
    widths = [len(c) for c in cols]
    for row in rows:
        for i, value in enumerate(row):
            widths[i] = max(widths[i], len(value))
    header = "| " + " | ".join(c.ljust(widths[i]) for i, c in enumerate(cols)) + " |"
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(cols))) + " |"
    body = ["| " + " | ".join(value.ljust(widths[i]) for i, value in enumerate(row)) + " |" for row in rows]
    return "\n".join([header, sep, *body])


def draw_fit_overlay(
    ax: plt.Axes,
    row: dict[str, object],
    width_cutoff: float,
    enrich_cutoff: float,
    show_ylabel: bool = True,
) -> None:
    path = Path(str(row["distance_file"]))
    x, y = read_vplot_points(path)
    group = str(row["group"])
    cmap = "Blues" if group == "TF motifs" else "Oranges"
    ax.hexbin(x, y, gridsize=46, extent=(-150, 150, 20, 150), cmap=cmap, mincnt=1, bins="log", linewidths=0)

    a = float(row.get("apex_x", 0.0))
    b = float(row.get("width", 0.0))
    xs = np.linspace(0, 150, 160)
    yy = b + 2 * xs
    keep = yy <= 150
    color = BLUE if group == "TF motifs" else ORANGE
    ax.plot(a + xs[keep], yy[keep], color=color, lw=1.8)
    ax.plot(a - xs[keep], yy[keep], color=color, lw=1.8)
    ax.scatter([a], [b], marker="*", s=70, color=RED, edgecolor="white", linewidth=0.6, zorder=4)

    ax.set_xlim(-150, 150)
    ax.set_ylim(20, 150)
    ax.set_xlabel("distance (bp)", fontsize=7.2, labelpad=1)
    if show_ylabel:
        ax.set_ylabel("fragment length (bp)", fontsize=7.2, labelpad=1)
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])
    ax.tick_params(labelsize=6.8, length=2.4, pad=1)
    ax.grid(True, color=GRID, lw=0.35, alpha=0.75)

    passes = (float(row["width"]) >= width_cutoff) and (float(row["enrichment"]) >= enrich_cutoff)
    call = "TF-footprint-like" if passes else "bias-like / weak V-like"
    label = str(row["name"]).replace("ATAC_Bias_", "")
    title = "TF motif: " + label if group == "TF motifs" else "bias control: " + label
    ax.set_title(title, fontsize=7.7, pad=2)
    ax.text(
        0.02,
        0.98,
        f"width={float(row['width']):.1f} bp\n"
        f"E={float(row['enrichment']):.2f}\n"
        f"LRT P={format_p(float(row.get('lrt_p', np.nan)))}\n"
        f"classification={call}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=6.8,
        color=DARK,
        bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="#c6ccd2", alpha=0.88),
    )


def draw_representative_vplots(
    ax: plt.Axes,
    examples: list[dict[str, object]],
    width_cutoff: float,
    enrich_cutoff: float,
) -> tuple[bool, list[str]]:
    ax.set_axis_off()
    ax.set_title("Representative fitted V-plots", loc="left", fontsize=10.6, fontweight="bold", pad=7)
    paths: list[str] = []
    if len(examples) >= 2:
        left = ax.inset_axes([0.02, 0.14, 0.46, 0.76])
        right = ax.inset_axes([0.52, 0.14, 0.46, 0.76])
        draw_fit_overlay(left, examples[0], width_cutoff, enrich_cutoff, show_ylabel=True)
        draw_fit_overlay(right, examples[1], width_cutoff, enrich_cutoff, show_ylabel=False)
        paths = [str(ex["distance_file"]) for ex in examples]
        return True, paths

    left = ax.inset_axes([0.04, 0.18, 0.42, 0.70])
    right = ax.inset_axes([0.55, 0.18, 0.42, 0.70])
    draw_schematic_vplot(left, b=50, strength=1.10, color=BLUE, seed=19)
    draw_schematic_vplot(right, b=5, strength=0.55, color=ORANGE, seed=23)
    left.set_title("schematic TF-footprint-like", fontsize=7.7)
    right.set_title("schematic enzyme-bias-like", fontsize=7.7)
    ax.text(
        0.50,
        0.05,
        "schematic panel: real distance distribution files were not found",
        transform=ax.transAxes,
        ha="center",
        fontsize=8.1,
        color=MID,
        fontweight="bold",
    )
    return False, paths


def draw_calibration_scatter(
    ax: plt.Axes,
    fit_summary: pd.DataFrame,
    cutoff_table: pd.DataFrame,
    width_cutoff: float,
    enrich_cutoff: float,
) -> None:
    ax.set_title("Calibration of the two-dimensional decision rule", loc="left", fontsize=9.8, fontweight="bold", pad=7)
    df = fit_summary.copy()
    tf = df[df["group"] == "TF motifs"]
    bias = df[df["group"] == "matched enzyme-bias controls"]

    ax.axvspan(width_cutoff, max(df["width"].max() * 1.08, width_cutoff + 5), ymin=0, ymax=1, color=PALE_GREEN, alpha=0.55, zorder=0)
    ax.axhspan(enrich_cutoff, max(df["enrichment"].max() * 1.12, enrich_cutoff + 0.5), xmin=0, xmax=1, color=PALE_GREEN, alpha=0.25, zorder=0)
    ax.scatter(bias["width"], bias["enrichment"], s=42, facecolor="white", edgecolor=ORANGE, lw=1.2, label="matched enzyme-bias controls", zorder=3)
    ax.scatter(tf["width"], tf["enrichment"], s=46, facecolor=BLUE, edgecolor="white", lw=0.55, label="TF motifs", zorder=4)
    ax.axvline(width_cutoff, color=DARK, ls="--", lw=1.1)
    ax.axhline(enrich_cutoff, color=DARK, ls="--", lw=1.1)
    ax.set_xlabel("model-fitted V-channel width, b (bp)", fontsize=8.0, labelpad=2)
    ax.set_ylabel("inside-V enrichment, log2(V-in/V-out)", fontsize=8.0, labelpad=2)
    ax.tick_params(labelsize=7.2, length=2.6, pad=1)
    ax.grid(True, color=GRID, lw=0.45, alpha=0.8)
    ax.legend(frameon=False, fontsize=7.0, loc="upper left", borderaxespad=0.2)

    xmax = max(df["width"].max() * 1.08, width_cutoff + 5)
    ymax = max(df["enrichment"].max() * 1.12, enrich_cutoff + 0.5)
    ax.set_xlim(min(-1, df["width"].min() - 1), xmax)
    ax.set_ylim(min(-0.08, df["enrichment"].min() - 0.08), ymax)
    ax.text(
        0.98,
        0.82,
        "TF-footprint-like zone\nwidth >= cutoff AND E >= cutoff",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.0,
        color="#1e5c40",
        bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="#b9d7c8", alpha=0.88),
    )
    ax.text(
        0.05,
        0.19,
        "bias-like /\nweak V-like",
        transform=ax.transAxes,
        ha="left",
        fontsize=7.3,
        color=ORANGE,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.02,
        "LRT detects V-like structure; final classification uses calibrated effect sizes.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.7,
        color=DARK,
        bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="#cbd2d8", alpha=0.90),
    )

    feature_col = infer_column(cutoff_table, ["feature"])
    cutoff_col = infer_column(cutoff_table, ["cutoff"])
    acc_col = infer_column(cutoff_table, ["loo_acc"], required=False)
    sens_col = infer_column(cutoff_table, ["loo_sens"], required=False)
    spec_col = infer_column(cutoff_table, ["loo_spec"], required=False)
    inset = ax.inset_axes([0.045, 0.58, 0.43, 0.25])
    inset.set_axis_off()
    rows = []
    for _, row in cutoff_table.iterrows():
        feature = str(row[feature_col])
        if "2d" in feature.lower():
            rows.append(["2D rule", "AND", f"{float(row[acc_col]):.2f}" if acc_col else "NA"])
        elif "width" in feature.lower():
            rows.append(["width", f"{float(row[cutoff_col]):.2f}", f"{float(row[acc_col]):.2f}" if acc_col else "NA"])
        elif "log2" in feature.lower() or "enrich" in feature.lower():
            rows.append(["E", f"{float(row[cutoff_col]):.2f}", f"{float(row[acc_col]):.2f}" if acc_col else "NA"])
    table = inset.table(cellText=rows, colLabels=["feature", "cutoff", "LOO acc"], loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(6.4)
    table.scale(1.0, 1.18)
    for (r, _c), cell in table.get_celld().items():
        cell.set_linewidth(0.45)
        cell.set_edgecolor("#8b959e")
        if r == 0:
            cell.set_facecolor("#edf2f6")
            cell.set_text_props(weight="bold")


def save_formats(fig: plt.Figure, out_prefix: Path, dpi: int, formats: list[str]) -> list[Path]:
    written: list[Path] = []
    for fmt in formats:
        path = out_prefix.with_suffix("." + fmt)
        try:
            fig.savefig(
                path,
                dpi=dpi if fmt.lower() == "png" else None,
                facecolor="white",
                bbox_inches="tight",
                pad_inches=0.04,
            )
        except PermissionError as exc:
            print(f"warning: could not write {path}: {exc}")
            continue
        written.append(path)
    return written


def write_notes(
    out_dir: Path,
    out_prefix: Path,
    input_root: Path,
    method_root: Path,
    threshold_path: Path,
    stats_path: Path,
    overview_outputs: list[Path],
    panel_outputs: list[Path] | None = None,
) -> Path:
    note_path = out_dir / "supp_methodB_revised_figure_notes.md"
    thresholds = pd.read_csv(threshold_path)
    keep_cols = [c for c in ["assay", "feature", "cutoff", "margin", "max_bias", "min_TF", "unit"] if c in thresholds.columns]
    threshold_md = dataframe_to_markdown(thresholds[keep_cols])
    outputs = "\n".join(f"- `{path}`" for path in overview_outputs)
    panel_text = ""
    if panel_outputs:
        panel_lines = "\n".join(f"- `{path}`" for path in panel_outputs)
        panel_text = f"\n\nIndividual panel files:\n\n{panel_lines}"

    text = f"""# Revised supplementary figure notes

## Outputs

{outputs}
{panel_text}
- `{note_path}`

## Input roots

- input root: `{input_root}`
- method root: `{method_root}`
- result directory: `{out_dir}`

## Reuse and inspection summary

- The revised script redraws the conceptual V-plot model panel, real loMNase apex-likelihood scan, workflow panel, and cut-off derivation panel.
- The previous two-dimensional scatter panel was removed from the overview.
- The figure avoids the internal project label and uses article-facing wording: V-plot model, V-channel width, inside-V enrichment, LRT, and calibrated effect-size cut-offs.

## Data sources

- cut-off derivation table: `{threshold_path}`
- separation/validation summary: `{stats_path}`

Cut-off derivation table:

{threshold_md}

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
"""
    note_path.write_text(text, encoding="utf-8")
    return note_path


def make_individual_panel_figures(method_root: Path, out_prefix: Path, dpi: int, formats: list[str]) -> list[Path]:
    panel_defs = [
        (
            "A",
            (8.2, 2.55),
            lambda ax: draw_vplot_model_mle(ax, method_root),
            dict(left=0.065, right=0.990, bottom=0.060, top=0.915),
        ),
        ("C", (8.6, 2.05), draw_workflow, dict(left=0.055, right=0.995, bottom=0.070, top=0.890)),
        (
            "D",
            (8.9, 3.25),
            lambda ax: draw_cutoff_derivation_panel(ax, method_root),
            dict(left=0.055, right=0.995, bottom=0.050, top=0.905),
        ),
    ]
    written: list[Path] = []
    for label, figsize, draw_func, margins in panel_defs:
        fig = plt.figure(figsize=figsize, constrained_layout=False)
        ax = fig.add_subplot(111)
        draw_func(ax)
        add_panel_label(ax, label)
        fig.subplots_adjust(**margins)
        panel_prefix = out_prefix.with_name(f"{out_prefix.name}_panel_{label}")
        written.extend(save_formats(fig, panel_prefix, dpi, formats))
        plt.close(fig)
    return written


def make_figure(
    input_root: Path,
    method_root: Path,
    out_prefix: Path,
    dpi: int,
    formats: list[str],
    use_real_vplots: str,
) -> tuple[list[Path], Path]:
    fig = plt.figure(figsize=(8.6, 10.2), constrained_layout=False)
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.05, 0.72, 1.22], hspace=0.36)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[2, 0])

    draw_vplot_model_mle(ax_a, method_root)
    draw_workflow(ax_c)
    _data_dir, threshold_path_d, stats_path = draw_cutoff_derivation_panel(ax_d, method_root)

    for axis, label in [(ax_a, "A"), (ax_c, "C"), (ax_d, "D")]:
        add_panel_label(axis, label)

    fig.suptitle(
        "Supplementary Fig. X | V-plot model-based classification of TF-footprint-like patterns\n"
        "using calibrated effect-size cut-offs.",
        fontsize=11.0,
        fontweight="bold",
        y=0.988,
    )
    fig.subplots_adjust(left=0.060, right=0.985, bottom=0.042, top=0.925)

    written = save_formats(fig, out_prefix, dpi, formats)
    plt.close(fig)

    panel_written = make_individual_panel_figures(method_root, out_prefix, dpi, formats)

    note_path = write_notes(
        out_prefix.parent,
        out_prefix,
        input_root,
        method_root,
        threshold_path_d,
        stats_path,
        written,
        panel_written,
    )
    return written + panel_written, note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    default_method_root = here.parent
    default_input_root = default_method_root / "V_plot_pipeline"
    parser.add_argument("--input-root", type=Path, default=default_input_root)
    parser.add_argument("--method-root", type=Path, default=default_method_root)
    parser.add_argument("--out-prefix", type=Path, default=here / "supp_methodB_revised_overview")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--use-real-vplots", choices=["auto", "yes", "no"], default="auto")
    parser.add_argument("--format", default="png,pdf,svg", help="Comma-separated output formats.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    formats = [part.strip().lower().lstrip(".") for part in str(args.format).split(",") if part.strip()]
    if not formats:
        raise ValueError("--format must contain at least one file format")
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    written, note_path = make_figure(
        input_root=args.input_root.resolve(),
        method_root=args.method_root.resolve(),
        out_prefix=args.out_prefix.resolve(),
        dpi=args.dpi,
        formats=formats,
        use_real_vplots=args.use_real_vplots,
    )
    for path in written:
        print(f"wrote {path}")
    print(f"wrote {note_path}")


if __name__ == "__main__":
    main()
