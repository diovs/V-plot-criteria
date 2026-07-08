#!/usr/bin/env python3
"""Plot real CTCF versus scrambled CTCF motifs for loMNase, DNase and ATAC."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
    }
)
import matplotlib.pyplot as plt
import pandas as pd


BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
OUTDIR = BASE
POINTS = DATA / "ctcf_scramble_points.csv"
STATS = DATA / "ctcf_scramble_stats.csv"

ASSAYS = ("loMNase", "DNase", "ATAC")
REAL_COLOR = "#104e8b"
SHUFFLED_COLOR = "#D55E00"
SHUFFLED_TEXT_COLOR = "#B44A00"
THRESHOLDS = {
    "loMNase": {"width": 6.708921215961457, "log2FC": 0.6417757861672906},
    "DNase": {"width": 9.154892066855576, "log2FC": 0.6304094709936522},
    "ATAC": {"width": 19.998395142799787, "log2FC": 0.4580046081180493},
}
SCRAMBLE_LABEL_BY_SEED = {
    0: "scramble 1",
    2: "scramble 2",
    5: "scramble 3",
    1: "scramble 4",
    3: "scramble 5",
    4: "scramble 6",
    6: "scramble 7",
    7: "scramble 8",
    8: "scramble 9",
    9: "scramble 10",
    10: "scramble 11",
}

LABEL_POSITIONS = {
    "loMNase": {
        0: (2.0, 0.08), 1: (24.2, 0.28), 2: (11.8, 0.18), 3: (24.2, 0.70),
        4: (2.0, 0.36), 5: (2.0, 0.22), 6: (2.0, -0.02), 7: (24.2, 0.84),
        8: (24.2, 0.56), 9: (24.2, 0.42), 10: (24.2, 0.98),
    },
    "DNase": {
        0: (5.8, 0.10), 1: (5.8, 0.23), 2: (5.8, -0.03), 3: (23.4, 0.66),
        4: (23.4, 0.53), 5: (23.4, 0.40), 6: (23.4, 0.27), 7: (13.8, 0.60),
        8: (23.4, 0.79), 9: (5.8, 0.36), 10: (13.8, 0.47),
    },
    "ATAC": {
        0: (2.4, 0.00), 1: (2.4, 0.12), 2: (2.4, -0.12), 3: (2.4, 0.24),
        4: (2.4, 0.36), 5: (2.4, 0.48), 6: (9.0, 0.00), 7: (9.0, 0.12),
        8: (9.0, 0.24), 9: (9.0, 0.36), 10: (9.0, 0.48),
    },
}


def scramble_name(seed: int) -> str:
    return SCRAMBLE_LABEL_BY_SEED[seed]


def format_p(value: float) -> str:
    return "NA" if pd.isna(value) else f"{value:.3e}"


def draw_assay(points: pd.DataFrame, stats_df: pd.DataFrame, assay: str) -> None:
    sub = points[points["assay"] == assay].copy()
    thr = THRESHOLDS[assay]
    xmax = max(sub["width"].max() * 1.18, thr["width"] * 1.45, 1.0)
    xmin = min(-1.1, sub["width"].min() - 1.1)
    ymax = max(sub["log2FC"].max() + 0.35, thr["log2FC"] * 1.55)
    ymin = min(-0.08, sub["log2FC"].min() - 0.15)
    stat = stats_df[stats_df["assay"] == assay].iloc[0]

    real = sub[sub["group"] == "real CTCF"].copy()
    shuf = sub[sub["group"] == "shuffled CTCF"].copy()
    shuf["seed"] = shuf["seed"].astype(int)

    fig, ax = plt.subplots(figsize=(8.5, 5.25))
    ax.add_patch(
        plt.Rectangle(
            (thr["width"], thr["log2FC"]),
            xmax - thr["width"],
            ymax - thr["log2FC"],
            facecolor="#E6EEF6",
            alpha=1.0,
            zorder=0,
        )
    )
    ax.axvline(thr["width"], ls=(0, (5, 4)), lw=1.4, color="#666666", zorder=1)
    ax.axhline(thr["log2FC"], ls=(0, (5, 4)), lw=1.4, color="#666666", zorder=1)
    ax.scatter(
        shuf["width"],
        shuf["log2FC"],
        s=62,
        facecolors="white",
        edgecolors=SHUFFLED_COLOR,
        linewidths=1.7,
        zorder=3,
    )
    ax.scatter(
        real["width"],
        real["log2FC"],
        s=98,
        facecolors=REAL_COLOR,
        edgecolors=REAL_COLOR,
        linewidths=1.5,
        zorder=4,
    )

    for _, row in shuf.sort_values("seed").iterrows():
        seed = int(row["seed"])
        lx, ly = LABEL_POSITIONS[assay].get(seed, (row["width"] + 0.45, row["log2FC"] + 0.04))
        ax.annotate(
            scramble_name(seed),
            xy=(row["width"], row["log2FC"]),
            xytext=(lx, ly),
            textcoords="data",
            color=SHUFFLED_TEXT_COLOR,
            fontsize=7.7,
            ha="left",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.25},
            arrowprops={
                "arrowstyle": "-",
                "color": SHUFFLED_TEXT_COLOR,
                "lw": 0.55,
                "alpha": 0.75,
                "shrinkA": 2,
                "shrinkB": 3,
            },
            zorder=5,
        )

    real_row = real.iloc[0]
    ax.text(real_row["width"] + 0.75, real_row["log2FC"] + 0.06, "CTCF", color=REAL_COLOR, fontsize=10.5)
    ax.text(
        0.015,
        0.965,
        (
            "Wilcoxon rank-sum (CTCF vs scramble 1-11)\n"
            f"Width: p = {format_p(stat['wilcoxon_width_p'])}\n"
            f"V$_{{in}}$ / V$_{{out}}$: p = {format_p(stat['wilcoxon_Vin_Vout_p'])}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.3,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 1.0, "pad": 1.4},
    )
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("V-channel width (bp)", fontsize=14)
    ax.set_ylabel("log2(Vin / Vout)", fontsize=14)
    ax.set_title(f"{assay}: real CTCF versus scrambled CTCF motifs", fontsize=14, fontweight="bold")
    ax.grid(True, color="#E6E6E6", lw=1.0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("black")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUTDIR / f"ctcf_scramble_scatter_{assay}.{ext}", dpi=300, transparent=False)
    plt.close(fig)


def main() -> None:
    points = pd.read_csv(POINTS)
    stats_df = pd.read_csv(STATS)
    points["display_label"] = points.apply(
        lambda row: "CTCF" if row["group"] == "real CTCF" else scramble_name(int(row["seed"])),
        axis=1,
    )
    for assay in ASSAYS:
        draw_assay(points, stats_df, assay)


if __name__ == "__main__":
    main()
