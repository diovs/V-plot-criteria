#!/usr/bin/env python3
"""Draw a four-panel CTCF/scrambled motif weblogo figure."""

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
import numpy as np
import pandas as pd
import logomaker


BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
OUT_PREFIX = BASE / "ctcf_scramble_logo_panel"
BASES = list("ACGT")
MOTIFS = [
    ("CTCF", DATA / "ctcf_original.meme"),
    ("scramble 1", DATA / "ctcf_scramble_1.meme"),
    ("scramble 2", DATA / "ctcf_scramble_2.meme"),
    ("scramble 3", DATA / "ctcf_scramble_3.meme"),
]


def read_meme_ppm(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    in_matrix = False
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if line.startswith("letter-probability"):
                in_matrix = True
                continue
            if not in_matrix:
                continue
            if len(fields) == 4:
                try:
                    rows.append([float(x) for x in fields])
                    continue
                except ValueError:
                    pass
            if rows:
                break
    if not rows:
        raise RuntimeError(f"No MEME probability matrix found in {path}")
    return np.array(rows)


def ppm_to_bits(ppm: np.ndarray) -> np.ndarray:
    p = np.clip(ppm, 1e-9, 1.0)
    entropy = -(p * np.log2(p)).sum(axis=1)
    information_content = 2.0 - entropy
    return ppm * information_content[:, None]


def draw_logo(ax: plt.Axes, title: str, meme_path: Path) -> None:
    ppm = read_meme_ppm(meme_path)
    bits = ppm_to_bits(ppm)
    df = pd.DataFrame(bits, columns=BASES)
    df.index = range(1, len(df) + 1)

    logomaker.Logo(df, ax=ax, color_scheme="classic", show_spines=False)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=4)
    ax.set_ylim(0, 2)
    ax.set_ylabel("bits", fontsize=9)
    ax.set_xlim(0.5, len(df) + 0.5)
    ax.set_xticks([1, 5, 10, 15, len(df)])
    ax.tick_params(axis="both", labelsize=8, length=2.5, width=0.8)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)


def main() -> None:
    fig = plt.figure(figsize=(10.2, 4.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], wspace=0.34, hspace=0.72)

    empty_left = fig.add_subplot(gs[0, 0])
    empty_right = fig.add_subplot(gs[0, 2])
    for ax in (empty_left, empty_right):
        ax.axis("off")

    axes = [
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[1, 2]),
    ]
    for ax, (title, meme_path) in zip(axes, MOTIFS):
        draw_logo(ax, title, meme_path)
        ax.set_xlabel("position", fontsize=9)

    fig.subplots_adjust(left=0.06, right=0.985, top=0.90, bottom=0.12)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{OUT_PREFIX}.{ext}", dpi=300, transparent=False)
    plt.close(fig)


if __name__ == "__main__":
    main()
