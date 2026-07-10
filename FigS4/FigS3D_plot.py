#!/usr/bin/env python3
"""Plot Fig. S3D: cut-off derivation for the two V-plot features.

The panel shows, for each assay, TF motif values and matched enzyme-bias values.
The shaded interval is the separation gap between the highest bias value and
the lowest TF value; the dashed line is the midpoint cut-off.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parent
DEFAULT_PREFIX = ROOT / "FigS3D_cutoff_derivation"
ASSAYS = ("loMNase", "DNase", "ATAC")

BLUE = "#1f5a93"
BIAS = "#ff6a6a"
PALE_BLUE = "#eaf2f8"
GRID = "#d7dde3"
DARK = "#222222"


def save_formats(fig: plt.Figure, prefix: Path, formats: list[str], dpi: int) -> list[Path]:
    written: list[Path] = []
    for fmt in formats:
        path = prefix.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)
    return written


def load_inputs(root: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    thresholds = pd.read_csv(root / "FigS3D_cutoff_thresholds.csv")
    assay_tables = {
        assay: pd.read_csv(root / f"FigS3D_plot_data_best_{assay}.csv")
        for assay in ASSAYS
    }
    return thresholds, assay_tables


def threshold_row(thresholds: pd.DataFrame, assay: str, feature_key: str) -> pd.Series:
    sub = thresholds[
        (thresholds["assay"].astype(str) == assay)
        & thresholds["feature"].astype(str).str.contains(feature_key, case=False, na=False)
    ]
    if sub.empty:
        raise ValueError(f"Missing {feature_key} threshold for {assay}")
    return sub.iloc[0]


def draw_metric(
    ax: plt.Axes,
    assay_tables: dict[str, pd.DataFrame],
    thresholds: pd.DataFrame,
    metric_col: str,
    feature_key: str,
    xlabel: str,
    xlim: tuple[float, float],
    seed: int,
) -> None:
    lane_y = {"loMNase": 2.0, "DNase": 1.0, "ATAC": 0.0}
    half_height = 0.22
    rng = np.random.default_rng(seed)

    for assay in ASSAYS:
        yc = lane_y[assay]
        df = assay_tables[assay]
        values = pd.to_numeric(df[metric_col], errors="coerce")
        is_bias = df["type"].astype(str).str.lower() == "bias"
        bias = values[is_bias].dropna()
        tf = values[~is_bias].dropna()
        row = threshold_row(thresholds, assay, feature_key)
        max_bias = float(row["max_bias"])
        min_tf = float(row["min_TF"])
        cutoff = float(row["cutoff"])

        ax.add_patch(
            Rectangle(
                (max_bias, yc - half_height),
                min_tf - max_bias,
                2 * half_height,
                facecolor=PALE_BLUE,
                edgecolor="none",
                zorder=0,
            )
        )
        ax.scatter(
            bias,
            yc + rng.uniform(-half_height * 0.45, half_height * 0.45, len(bias)),
            s=25,
            facecolors="none",
            edgecolors=BIAS,
            linewidths=1.1,
            zorder=3,
        )
        ax.scatter(
            tf,
            yc + rng.uniform(-half_height * 0.45, half_height * 0.45, len(tf)),
            s=24,
            color=BLUE,
            edgecolors="none",
            zorder=3,
        )
        ax.plot(
            [cutoff, cutoff],
            [yc - half_height - 0.08, yc + half_height + 0.08],
            color=DARK,
            ls="--",
            lw=1.2,
            zorder=4,
        )

    ax.set_yticks([lane_y[a] for a in ASSAYS])
    ax.set_yticklabels(ASSAYS)
    ax.set_ylim(-0.55, 2.55)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color=GRID, lw=0.5, alpha=0.7)
    ax.tick_params(axis="both", length=3, pad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


def plot_cutoffs(input_root: Path, out_prefix: Path, formats: list[str], dpi: int) -> list[Path]:
    thresholds, assay_tables = load_inputs(input_root)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    draw_metric(
        axes[0],
        assay_tables,
        thresholds,
        metric_col="log2FC",
        feature_key="enrichment",
        xlabel="log2(Vin / Vout)",
        xlim=(-0.30, 3.70),
        seed=1,
    )
    draw_metric(
        axes[1],
        assay_tables,
        thresholds,
        metric_col="width",
        feature_key="width",
        xlabel="V-channel width (bp)",
        xlim=(-2.5, 62.0),
        seed=2,
    )

    axes[1].tick_params(axis="y", labelleft=False)
    handles = [
        plt.Line2D([], [], marker="o", color=BLUE, ls="none", ms=5.3, label="TF motif"),
        plt.Line2D([], [], marker="o", mfc="white", mec=BIAS, mew=1.1, color=BIAS, ls="none", ms=5.3, label="enzyme bias"),
    ]
    axes[1].legend(handles=handles, frameon=False, loc="lower right", fontsize=8)
    fig.tight_layout(w_pad=1.2)
    written = save_formats(fig, out_prefix, formats, dpi)
    plt.close(fig)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=ROOT, help="Directory containing FigS3D CSV files.")
    parser.add_argument("--out-prefix", type=Path, default=DEFAULT_PREFIX, help="Output path without extension.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated output formats.")
    parser.add_argument("--dpi", type=int, default=300, help="Raster output DPI.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    formats = [item.strip().lstrip(".") for item in args.formats.split(",") if item.strip()]
    written = plot_cutoffs(args.input_root, args.out_prefix, formats, args.dpi)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
