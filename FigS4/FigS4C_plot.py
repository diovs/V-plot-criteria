#!/usr/bin/env python3
"""Plot Fig. S4C: loMNase apex-position likelihood profile.

This panel uses a precomputed likelihood scan from the Vplot model. For each
candidate apex x-position, the curve shows the relative log-likelihood per
fragment, normalized to the maximum for that motif or bias control.
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


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "FigS4C_real_loMNase_apex_likelihood_profiles.csv"
DEFAULT_PREFIX = ROOT / "FigS4C_loMNase_apex_likelihood"

BLUE = "#1f5a93"
ORANGE = "#d46b32"
RED = "#d9534f"
GRID = "#d7dde3"


def save_formats(fig: plt.Figure, prefix: Path, formats: list[str], dpi: int) -> list[Path]:
    written: list[Path] = []
    for fmt in formats:
        path = prefix.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)
    return written


def plot_profile(input_csv: Path, out_prefix: Path, formats: list[str], dpi: int) -> list[Path]:
    df = pd.read_csv(input_csv)
    required = {"label", "candidate_apex_x", "relative_loglik", "n"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {input_csv}: {', '.join(sorted(missing))}")

    fig, ax = plt.subplots(figsize=(5.3, 3.1))
    series = [
        ("CTCF motif", "CTCF", BLUE, 1.9),
        ("bias AAG", "AAG", ORANGE, 1.7),
    ]

    y_min = 0.0
    best_x = None
    for source_label, plot_label, color, linewidth in series:
        sub = df[df["label"] == source_label].copy()
        if sub.empty:
            raise ValueError(f"No rows labeled {source_label!r} in {input_csv}")
        sub = sub.sort_values("candidate_apex_x")
        x = pd.to_numeric(sub["candidate_apex_x"], errors="coerce")
        y = pd.to_numeric(sub["relative_loglik"], errors="coerce") / pd.to_numeric(sub["n"], errors="coerce")
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        ax.plot(x, y, color=color, lw=linewidth, label=plot_label)
        y_min = min(y_min, float(y.min()))
        if source_label == "CTCF motif":
            best_x = float(x.iloc[int(np.argmax(y.to_numpy()))])

    if best_x is not None:
        ax.axvline(best_x, color=RED, ls="--", lw=1.0)

    ax.axhline(0, color="#333333", lw=0.5, alpha=0.45)
    ax.set_xlabel("candidate apex position (bp)")
    ax.set_ylabel("Relative ln(likelihood)")
    ax.set_xlim(-11, 11)
    ax.set_ylim(y_min - 0.03, 0.03)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.7)
    ax.legend(frameon=False, loc="upper right")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    fig.tight_layout()
    written = save_formats(fig, out_prefix, formats, dpi)
    plt.close(fig)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Likelihood-profile CSV.")
    parser.add_argument("--out-prefix", type=Path, default=DEFAULT_PREFIX, help="Output path without extension.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated output formats.")
    parser.add_argument("--dpi", type=int, default=300, help="Raster output DPI.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    formats = [item.strip().lstrip(".") for item in args.formats.split(",") if item.strip()]
    written = plot_profile(args.input, args.out_prefix, formats, args.dpi)
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
