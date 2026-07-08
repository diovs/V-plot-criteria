#!/usr/bin/env python3
"""Draft a compact supplementary figure for the portable Method B pipeline."""
from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTBASE = ROOT / "supp_methodB_pipeline_overview"

MODEL_IMG = ROOT / "method_schematic.png"
SCATTER_IMG = ROOT / "ATAC_methodB_scatter.png"
CUTOFFS = ROOT / "ATAC_methodB_cutoffs.csv"

BLUE = "#104e8b"
RED = "#ff6a6a"
PALE = "#eef3f8"
EDGE = "#49657a"
TEXT = "#222222"


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.04,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def draw_workflow(ax) -> None:
    ax.set_axis_off()
    ax.set_title("End-to-end workflow", loc="left", fontsize=14, fontweight="bold", pad=8)
    steps = [
        ("0", "Fragment BED", "midpoint + length"),
        ("1", "Bias k-mers", "genome-wide sites"),
        ("2", "TF motifs", "TF V-plot distances"),
        ("3", "Bias controls", "bias V-plot distances"),
        ("4", "Generative fit", "apex, width, E, LR"),
        ("5", "2D rule", "cut-offs + LOO"),
    ]
    xs = [0.05, 0.37, 0.69, 0.05, 0.37, 0.69]
    ys = [0.67, 0.67, 0.67, 0.27, 0.27, 0.27]
    w, h = 0.24, 0.22

    for i, ((num, title, desc), x, y) in enumerate(zip(steps, xs, ys)):
        box = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.018,rounding_size=0.02",
            linewidth=1.2,
            edgecolor=EDGE,
            facecolor=PALE if i < 4 else "#e9f2ea",
        )
        ax.add_patch(box)
        ax.text(x + 0.03, y + h - 0.055, num, fontsize=15, fontweight="bold", color=BLUE)
        ax.text(x + 0.085, y + h - 0.06, title, fontsize=10.5, fontweight="bold", color=TEXT)
        ax.text(
            x + 0.085,
            y + 0.052,
            "\n".join(textwrap.wrap(desc, width=20)),
            fontsize=9.5,
            color="#444444",
            va="bottom",
        )

    arrows = [
        ((0.05 + w, 0.78), (0.37, 0.78)),
        ((0.37 + w, 0.78), (0.69, 0.78)),
        ((0.69 + w / 2, 0.67), (0.05 + w / 2, 0.49)),
        ((0.05 + w, 0.38), (0.37, 0.38)),
        ((0.37 + w, 0.38), (0.69, 0.38)),
    ]
    for start, end in arrows:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.2,
                color="#555555",
                connectionstyle="arc3,rad=0.0",
            )
        )

    ax.text(
        0.05,
        0.08,
        "Output: apex-score TSVs, a width x enrichment scatter plot, and assay-calibrated cut-offs.",
        fontsize=9.6,
        color="#333333",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def show_image_panel(ax, path: Path, title: str) -> None:
    ax.imshow(mpimg.imread(path))
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=8)
    ax.set_axis_off()


def draw_cutoff_table(ax) -> None:
    ax.set_axis_off()
    ax.set_title("Example output: ATAC cut-offs and LOO validation", loc="left",
                 fontsize=14, fontweight="bold", pad=8)
    df = pd.read_csv(CUTOFFS)
    display = []
    for row in df.itertuples(index=False):
        feature = str(row.feature)
        if feature.startswith("V-channel"):
            feature = "V-channel width"
            cutoff = f"{float(row.cutoff):.2f} bp"
        elif feature.startswith("log2"):
            feature = "Inside-V enrichment"
            cutoff = f"{float(row.cutoff):.2f}"
        else:
            feature = "2D rule"
            cutoff = (
                str(row.cutoff)
                .replace("width>=", "width >= ")
                .replace("E>=", "E >= ")
                .replace(" & ", "\n")
            )
        display.append([
            feature,
            cutoff,
            f"{float(row.loo_acc):.2f}",
            f"{float(row.loo_sens):.2f}",
            f"{float(row.loo_spec):.2f}",
        ])

    table = ax.table(
        cellText=display,
        colLabels=["Feature", "Cut-off", "LOO acc", "Sens", "Spec"],
        cellLoc="center",
        colLoc="center",
        loc="center",
        colWidths=[0.31, 0.38, 0.12, 0.095, 0.095],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.3)
    table.scale(1, 1.65)
    for (r, c), cell in table.get_celld().items():
        cell.set_linewidth(0.8)
        cell.set_edgecolor("#222222")
        if r == 0:
            cell.set_facecolor("#dbe6f1")
            cell.set_text_props(weight="bold")
        elif c == 0:
            cell.set_text_props(weight="bold")

    ax.text(
        0.02,
        0.08,
        "Cut-offs are learned from the gap between TF motifs and matched enzyme-bias controls.",
        transform=ax.transAxes,
        fontsize=9.8,
        color="#333333",
        va="bottom",
    )
    ax.text(
        0.02,
        0.015,
        "A motif is called a footprint only when both features exceed their cut-offs.",
        transform=ax.transAxes,
        fontsize=9.8,
        color="#333333",
        va="bottom",
    )


def main() -> None:
    fig = plt.figure(figsize=(13.2, 9.3))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.0], height_ratios=[1, 1],
                          wspace=0.16, hspace=0.24)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    draw_workflow(ax_a)
    show_image_panel(ax_b, MODEL_IMG, "Model-derived features")
    show_image_panel(ax_c, SCATTER_IMG, "Example calibrated scatter")
    draw_cutoff_table(ax_d)

    for ax, label in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D")]:
        panel_label(ax, label)

    fig.suptitle("Portable Method B pipeline for V-plot footprint calling",
                 fontsize=18, fontweight="bold", y=0.985)
    fig.savefig(OUTBASE.with_suffix(".png"), dpi=220, facecolor="white", bbox_inches="tight")
    fig.savefig(OUTBASE.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    print(f"wrote {OUTBASE.with_suffix('.png')}")
    print(f"wrote {OUTBASE.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
