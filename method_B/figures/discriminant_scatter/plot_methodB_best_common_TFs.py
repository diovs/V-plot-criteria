#!/usr/bin/env python3
"""
Redraw Method B discriminant scatter plots with assay-specific best TF sets.

The original figures in v_apex_scores/methodB_scatter were produced by
methodB_scatter_multiAssay_LOCAL.R from hand-picked TFs. This script makes the
TF choice explicit while keeping the same bias controls and visual grammar:

1. Use the full TF apex score tables under new_script/method_B/.
2. Exclude NFIX and keep NFIC when the two NFI motifs compete.
3. Rank TFs separately in loMNase, DNase and ATAC by their 2-D distance from the
   assay-matched bias controls, then add a small fixed set of common anchors.
4. Plot only two classes, TF and enzyme bias, matching the original scatter style.
5. Report the TFs common to all three selected panels and local MEME motif AT/GC
   composition where motif matrices are available.

Outputs are written to:
  new_script/method_B/v_apex_scores/methodB_scatter_best_assay_specific/
"""
from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "v_apex_scores"
MEME_DIR = ROOT / "weblogo_out" / "meme"
OUTDIR = BASE / "methodB_scatter_best_assay_specific"

ASSAYS = ("loMNase", "DNase", "ATAC")
TOP_N_PER_ASSAY = 12
COMMON_ANCHORS = ("CTCF", "NFIC", "TFE3", "ZNF281")
EXCLUDE_TFS = {"NFIX"}
AT_RICH_FAMILY_BONUS = 0.03
TF_FAMILY_GROUPS = {
    "SP": ("SP1", "SP2"),
    "KLF": ("KLF1", "KLF6", "KLF10", "KLF13", "KLF16"),
    "CEBP_DDIT": ("CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "DDIT3"),
    "CREB_ATF_CREM": (
        "ATF1",
        "ATF2",
        "ATF3",
        "ATF4",
        "ATF6",
        "ATF7",
        "CREB1",
        "CREB3",
        "CREB5",
        "CREM",
    ),
    "AP1_JUN_FOS": ("FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"),
    "MAF_SMALL": ("MAFF", "MAFG", "MAFK"),
}
TF_FAMILY = {
    tf: family for family, members in TF_FAMILY_GROUPS.items() for tf in members
}
BIAS_KEEP = {
    "loMNase": {"AAG", "ACA", "AGA", "TGG", "AGT"},
    "DNase": {"TCC", "TGA", "TGC", "TGG", "TGT"},
    "ATAC": {"GCC", "GCT", "GGC", "GTC", "GTG"},
}

TF_TABLE = {
    assay: ROOT / f"{assay}_TF_apex_scores.tsv" for assay in ASSAYS
}
BIAS_TABLE = {
    assay: BASE / f"{assay}_Bias_3bp_apex_scores.tsv" for assay in ASSAYS
}


def clean_name(name: str, kind: str) -> str:
    if kind == "bias":
        m = re.search(r"_Bias_([ACGTNacgtn]+)_", name)
        return m.group(1).upper() if m else name
    return re.sub(r"_(midP|loMNase|DNase|ATAC)_fragL_dist$", "", name)


def read_score_table(path: Path, kind: str, assay: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    ok = (
        (df["status"] == "ok")
        & np.isfinite(pd.to_numeric(df["enrichment_fold"], errors="coerce"))
        & (pd.to_numeric(df["enrichment_fold"], errors="coerce") > 0)
        & np.isfinite(pd.to_numeric(df["apex_y_channel_width"], errors="coerce"))
    )
    df = df.loc[ok].copy()
    df["TF"] = [clean_name(x, kind) for x in df["motif"]]
    df["assay"] = assay
    df["type"] = kind
    df["width"] = pd.to_numeric(df["apex_y_channel_width"], errors="coerce")
    df["FC"] = pd.to_numeric(df["enrichment_fold"], errors="coerce")
    df["log2FC"] = np.log2(df["FC"])
    return df[["assay", "TF", "type", "width", "FC", "log2FC", "motif"]]


def parse_meme(path: Path) -> dict[str, object] | None:
    rows: list[list[float]] = []
    in_matrix = False
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("letter-probability matrix"):
                in_matrix = True
                continue
            if in_matrix:
                if not line:
                    break
                parts = line.split()
                if len(parts) == 4:
                    rows.append([float(x) for x in parts])
    if not rows:
        return None
    arr = np.asarray(rows, dtype=float)
    avg = arr.mean(axis=0)
    base_range = float(avg.max() - avg.min())
    consensus = "".join("ACGT"[int(i)] for i in arr.argmax(axis=1))
    tf = path.stem
    return {
        "TF": tf,
        "motif_width": arr.shape[0],
        "A": avg[0],
        "C": avg[1],
        "G": avg[2],
        "T": avg[3],
        "GC": avg[1] + avg[2],
        "AT": avg[0] + avg[3],
        "base_range": base_range,
        "consensus": consensus,
        "base_class": (
            "GC-rich"
            if avg[1] + avg[2] >= 0.70
            else "AT-rich"
            if avg[0] + avg[3] >= 0.55
            else "balanced"
        ),
        "four_base_class": "ATCG-balanced" if base_range <= 0.15 else "skewed",
    }


def motif_composition() -> pd.DataFrame:
    rows = []
    for path in sorted(MEME_DIR.glob("*.meme")):
        row = parse_meme(path)
        if row:
            rows.append(row)
    return pd.DataFrame(rows)


def auc_mw(values: np.ndarray, is_tf: np.ndarray) -> float:
    pos = values[is_tf]
    neg = values[~is_tf]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # Mann-Whitney U with average ties, larger value = more TF-like.
    ranks = stats.rankdata(np.concatenate([pos, neg]), method="average")
    r_pos = ranks[: len(pos)].sum()
    u = r_pos - len(pos) * (len(pos) + 1) / 2.0
    return float(u / (len(pos) * len(neg)))


def wilcoxon_rank_sum_p(tf_v: np.ndarray, bias_v: np.ndarray) -> float:
    # Wilcoxon rank-sum and Mann-Whitney U are equivalent for two independent
    # groups; use one asymptotic, tie-corrected mode for both plotted features.
    return float(
        stats.mannwhitneyu(
            tf_v,
            bias_v,
            alternative="two-sided",
            method="asymptotic",
        ).pvalue
    )


def youden_threshold(values: np.ndarray, is_tf: np.ndarray) -> float:
    vals = np.unique(values[np.isfinite(values)])
    if len(vals) <= 1:
        return float(vals[0]) if len(vals) else float("nan")
    cuts = (vals[:-1] + vals[1:]) / 2.0
    best_cut = cuts[0]
    best_j = -math.inf
    p = is_tf.sum()
    n = (~is_tf).sum()
    for cut in cuts:
        pred = values >= cut
        sens = (pred & is_tf).sum() / p if p else 0.0
        spec = ((~pred) & (~is_tf)).sum() / n if n else 0.0
        j = sens + spec - 1.0
        if j > best_j:
            best_cut = cut
            best_j = j
    return float(best_cut)


def format_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    return f"{p:.3e}"


def score_and_select(tf_all: dict[str, pd.DataFrame], bias_all: dict[str, pd.DataFrame],
                     comp: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    ranked_rows = []
    selected_rows = []
    for assay in ASSAYS:
        tf = tf_all[assay][~tf_all[assay]["TF"].isin(EXCLUDE_TFS)].copy()
        bias = bias_all[assay][bias_all[assay]["TF"].isin(BIAS_KEEP[assay])].copy()
        max_bias_w = float(bias["width"].max())
        max_bias_e = float(bias["log2FC"].max())
        width_range = max(float(tf["width"].max()) - max_bias_w, 1e-9)
        enrich_range = max(float(tf["log2FC"].max()) - max_bias_e, 1e-9)
        tf["delta_width"] = tf["width"] - max_bias_w
        tf["delta_log2FC"] = tf["log2FC"] - max_bias_e
        tf["rank_score"] = 0.5 * (tf["delta_width"] / width_range) + 0.5 * (
            tf["delta_log2FC"] / enrich_range
        )
        tf["family"] = tf["TF"].map(TF_FAMILY).fillna(tf["TF"])
        tf = tf.merge(
            comp[
                [
                    "TF",
                    "GC",
                    "AT",
                    "base_range",
                    "base_class",
                    "four_base_class",
                    "consensus",
                ]
            ],
            on="TF",
            how="left",
        )
        tf["family_score"] = tf["rank_score"] + np.where(
            tf["base_class"].eq("AT-rich"), AT_RICH_FAMILY_BONUS, 0.0
        )
        tf = tf.sort_values(
            ["rank_score", "delta_width", "delta_log2FC"], ascending=False
        ).reset_index(drop=True)
        tf["assay_rank"] = np.arange(1, len(tf) + 1)
        representative_tfs = set(
            tf.sort_values(
                ["family_score", "rank_score", "delta_width", "delta_log2FC"],
                ascending=False,
            )
            .drop_duplicates("family")["TF"]
        )
        tf["family_representative"] = tf["TF"].isin(representative_tfs)
        ranked_rows.append(tf)

        family_reps = (
            tf.sort_values(
                ["family_score", "rank_score", "delta_width", "delta_log2FC"],
                ascending=False,
            )
            .drop_duplicates("family")
            .sort_values(
                ["rank_score", "delta_width", "delta_log2FC"], ascending=False
            )
            .reset_index(drop=True)
        )
        anchor_set = set(COMMON_ANCHORS)
        anchors = family_reps[family_reps["TF"].isin(anchor_set)].copy()
        remaining = family_reps[~family_reps["TF"].isin(anchor_set)].head(
            max(TOP_N_PER_ASSAY - len(anchors), 0)
        )
        take = pd.concat([anchors, remaining], ignore_index=True)
        take = take.sort_values(
            ["rank_score", "delta_width", "delta_log2FC"], ascending=False
        ).head(TOP_N_PER_ASSAY).copy()
        take["selected"] = True
        selected_rows.append(take)

    ranked = pd.concat(ranked_rows, ignore_index=True)
    selected = pd.concat(selected_rows, ignore_index=True)
    selected_sets = [set(selected.loc[selected["assay"] == a, "TF"]) for a in ASSAYS]
    common_selected = set.intersection(*selected_sets)
    selected["common_in_three_selected_panels"] = selected["TF"].isin(common_selected)
    return ranked, selected, common_selected


def stats_table(plot_df: pd.DataFrame, assay: str) -> tuple[pd.DataFrame, float, float]:
    is_tf = (plot_df["type"] == "TF").to_numpy()
    rows = []
    thresholds = {}
    for col, label in [
        ("width", "V-channel width (bp)"),
        ("log2FC", "log2(V-in/V-out)"),
    ]:
        values = plot_df[col].to_numpy(dtype=float)
        tf_v = values[is_tf]
        bias_v = values[~is_tf]
        wrs_p = wilcoxon_rank_sum_p(tf_v, bias_v)
        thr = youden_threshold(values, is_tf)
        pred = values >= thr
        sens = (pred & is_tf).sum() / is_tf.sum()
        spec = ((~pred) & (~is_tf)).sum() / (~is_tf).sum()
        rows.append(
            {
                "assay": assay,
                "metric": label,
                "median_TF": np.median(tf_v),
                "median_bias": np.median(bias_v),
                "WRS_p": wrs_p,
                "MW_p": wrs_p,
                "AUC": auc_mw(values, is_tf),
                "youden_thr": thr,
                "sens": sens,
                "spec": spec,
            }
        )
        thresholds[col] = thr
    rule = (plot_df["width"] >= thresholds["width"]) & (
        plot_df["log2FC"] >= thresholds["log2FC"]
    )
    truth = plot_df["type"] == "TF"
    rows.append(
        {
            "assay": assay,
            "metric": "2D rule (width AND E)",
            "median_TF": np.nan,
            "median_bias": np.nan,
            "WRS_p": np.nan,
            "MW_p": np.nan,
            "AUC": np.nan,
            "youden_thr": np.nan,
            "sens": (rule & truth).sum() / truth.sum(),
            "spec": ((~rule) & (~truth)).sum() / (~truth).sum(),
        }
    )
    return pd.DataFrame(rows), thresholds["width"], thresholds["log2FC"]


LABEL_OFFSETS = {
    "loMNase": {
        "CTCF": (42, 20), "REST": (44, -18), "NFIC": (48, 24),
        "SP1": (46, 28), "NRF1": (78, 24), "ZNF281": (-82, -40),
        "PATZ1": (-80, 28), "KLF16": (62, 36), "ZBTB33": (82, -28),
        "ZNF324": (68, -40), "ZNF148": (70, -52), "TFE3": (-64, 22),
        "AAG": (-22, 14), "ACA": (14, -18),
        "AGA": (20, -6), "AGT": (-22, -10), "TGG": (16, 10),
    },
    "DNase": {
        "CTCF": (48, 24), "NFIC": (56, 42), "REST": (50, 32),
        "NFE2": (62, 58), "SP1": (66, 26), "BACH1": (-70, -32),
        "ZNF281": (64, -16), "MAFG": (-66, 46), "KLF16": (60, -30),
        "KLF1": (-10, 10), "SP2": (9, 7),
        "TFE3": (-70, -40), "JUND": (54, -42),
        "TCC": (-26, 12), "TGA": (18, -12),
        "TGC": (18, 2), "TGG": (16, 14), "TGT": (-24, -12),
    },
    "ATAC": {
        "CTCF": (46, 22), "NFIC": (52, 46), "JUND": (-70, 54),
        "CEBPG": (66, 30), "CREM": (11, -4), "TFE3": (-82, 34),
        "CREB1": (68, 4), "CEBPB": (10, -17), "DDIT3": (-10, -13),
        "MLX": (-82, -10), "MECOM": (-82, -44), "GMEB1": (54, -52),
        "IRF9": (-82, -34), "STAT5B": (58, -44),
        "ATF2": (-12, -18), "ZNF281": (70, -26),
        "GCC": (-24, 12), "GCT": (-24, -11),
        "GGC": (16, 5), "GTC": (24, -10), "GTG": (16, -18),
    },
}


TF_LABEL_POSITIONS = {
    "loMNase": {
        "CTCF": (34.6, 3.68),
        "NFIC": (22.6, 2.12),
        "REST": (24.2, 1.98),
        "SP1": (16.4, 2.03),
        "KLF16": (17.0, 1.68),
        "NRF1": (21.0, 1.62),
        "PATZ1": (15.1, 1.38),
        "ZBTB33": (20.7, 1.35),
        "ZNF281": (15.6, 1.02),
        "ZNF148": (20.4, 0.82),
        "ZNF324": (22.5, 0.70),
        "TFE3": (11.0, 1.34),
    },
    "DNase": {
        "CTCF": (29.1, 3.24),
        "NFIC": (25.0, 2.40),
        "NFE2": (20.5, 2.05),
        "SP1": (21.5, 1.92),
        "MAFG": (17.6, 1.86),
        "BACH1": (17.3, 1.52),
        "JUND": (20.7, 1.30),
        "KLF16": (22.7, 1.55),
        "REST": (29.0, 1.55),
        "TFE3": (17.0, 1.15),
        "ZNF281": (24.6, 0.95),
        "ZNF740": (28.8, 0.75),
    },
    "ATAC": {
        "CTCF": (61.8, 3.18),
        "NFIC": (51.0, 1.65),
        "CEBPG": (49.2, 1.53),
        "JUND": (47.2, 1.42),
        "TFE3": (40.6, 1.38),
        "MLX": (39.0, 1.20),
        "CREB1": (50.2, 1.23),
        "MECOM": (39.2, 0.98),
        "ZNF281": (51.8, 0.88),
        "GMEB1": (43.0, 0.76),
        "IRF9": (39.0, 0.62),
        "STAT5B": (46.2, 0.53),
    },
}


TF_LABEL_CANDIDATES = [
    (18, 14),
    (18, -14),
    (-18, 14),
    (-18, -14),
    (28, 0),
    (-28, 0),
    (0, 26),
    (0, -26),
    (32, 18),
    (-32, 18),
    (32, -18),
    (-32, -18),
    (44, 0),
    (-44, 0),
    (0, 38),
    (0, -38),
    (52, 24),
    (-52, 24),
    (52, -24),
    (-52, -24),
    (64, 0),
    (-64, 0),
    (0, 52),
    (0, -52),
]

BIAS_LABEL_CANDIDATES = [
    (10, 10),
    (10, -10),
    (10, 24),
    (10, -24),
    (10, 38),
    (10, 52),
    (22, 6),
    (22, -6),
    (22, 20),
    (22, -20),
    (22, 34),
    (36, 0),
    (36, 16),
    (36, -16),
    (12, 20),
    (12, -20),
    (46, 12),
    (46, -12),
    (0, 22),
    (0, -22),
]


def unique_offsets(primary: tuple[int, int],
                   candidates: list[tuple[int, int]]) -> list[tuple[int, int]]:
    offsets = [primary] + candidates
    seen = set()
    out = []
    for off in offsets:
        if off not in seen:
            seen.add(off)
            out.append(off)
    return out


def min_offset(offset: tuple[int, int], min_points: float) -> tuple[int, int]:
    if min_points <= 0:
        return offset
    dx, dy = offset
    dist = math.hypot(dx, dy)
    if dist >= min_points:
        return offset
    if dist == 0:
        return (int(min_points), 0)
    scale = min_points / dist
    return (int(round(dx * scale)), int(round(dy * scale)))


def bbox_inside(inner, outer) -> bool:
    return (
        inner.x0 >= outer.x0
        and inner.x1 <= outer.x1
        and inner.y0 >= outer.y0
        and inner.y1 <= outer.y1
    )


def point_boxes(ax, rows: pd.DataFrame, pad: float = 8.0) -> list[Bbox]:
    boxes = []
    for _, row in rows.iterrows():
        x, y = ax.transData.transform((row["width"], row["log2FC"]))
        boxes.append(Bbox.from_extents(x - pad, y - pad, x + pad, y + pad))
    return boxes


def place_tf_labels_data(
    fig,
    ax,
    rows: pd.DataFrame,
    positions: dict[str, tuple[float, float]],
    color: str,
    fontsize: float,
    placed_boxes: list | None = None,
) -> list:
    boxes = list(placed_boxes or [])
    for _, row in rows.iterrows():
        tf_name = row["TF"]
        if tf_name not in positions:
            continue
        lx, ly = positions[tf_name]
        ha = "left" if lx >= row["width"] else "right"
        ann = ax.annotate(
            tf_name,
            xy=(row["width"], row["log2FC"]),
            xytext=(lx, ly),
            textcoords="data",
            color=color,
            fontsize=fontsize,
            ha=ha,
            va="center",
            clip_on=True,
            arrowprops={
                "arrowstyle": "-",
                "color": color,
                "lw": 0.6,
                "alpha": 0.75,
                "shrinkA": 2,
                "shrinkB": 4,
            },
        )
        fig.canvas.draw()
        boxes.append(ann.get_window_extent(fig.canvas.get_renderer()).padded(2))
    return boxes


def place_labels(
    fig,
    ax,
    rows: pd.DataFrame,
    offsets: dict[str, tuple[int, int]],
    color: str,
    fontsize: float,
    candidates: list[tuple[int, int]],
    placed_boxes: list | None = None,
    connector_color: str | None = None,
    min_offset_points: float = 0.0,
) -> list:
    boxes = list(placed_boxes or [])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer).padded(-3)

    for _, row in rows.iterrows():
        primary = min_offset(offsets.get(row["TF"], candidates[0]), min_offset_points)
        placed = False
        label_offsets = unique_offsets(
            primary,
            [min_offset(off, min_offset_points) for off in candidates],
        )
        for dx, dy in label_offsets:
            ann = ax.annotate(
                row["TF"],
                xy=(row["width"], row["log2FC"]),
                xytext=(dx, dy),
                textcoords="offset points",
                color=color,
                fontsize=fontsize,
                ha="left" if dx >= 0 else "right",
                va="center",
                clip_on=True,
                arrowprops=(
                    {
                        "arrowstyle": "-",
                        "color": connector_color,
                        "lw": 0.55,
                        "alpha": 0.75,
                        "shrinkA": 2,
                        "shrinkB": 3,
                    }
                    if connector_color
                    else None
                ),
            )
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            box = ann.get_window_extent(renderer).padded(2)
            if bbox_inside(box, axes_box) and not any(box.overlaps(b) for b in boxes):
                boxes.append(box)
                placed = True
                break
            ann.remove()

        if not placed:
            # Last resort: keep the label clipped within the axes instead of
            # letting it run outside the plotting area.
            dx, dy = min_offset(candidates[0], min_offset_points)
            ann = ax.annotate(
                row["TF"],
                xy=(row["width"], row["log2FC"]),
                xytext=(dx, dy),
                textcoords="offset points",
                color=color,
                fontsize=fontsize,
                ha="left",
                va="center",
                clip_on=True,
                arrowprops=(
                    {
                        "arrowstyle": "-",
                        "color": connector_color,
                        "lw": 0.55,
                        "alpha": 0.75,
                        "shrinkA": 2,
                        "shrinkB": 3,
                    }
                    if connector_color
                    else None
                ),
            )
            fig.canvas.draw()
            boxes.append(ann.get_window_extent(fig.canvas.get_renderer()).padded(2))
    return boxes


def plot_assay(plot_df: pd.DataFrame, assay: str, thr_w: float, thr_e: float) -> None:
    tf = plot_df[plot_df["type"] == "TF"].copy()
    bias = plot_df[plot_df["type"] == "bias"].copy()
    xmax = max(plot_df["width"].max() * 1.15, thr_w * 1.4, 1)
    xmin = min(-1.2, plot_df["width"].min() - 0.5)
    ymin = min(plot_df["log2FC"].min() - 0.2, -0.15)
    ymax = max(plot_df["log2FC"].max() + 0.35, thr_e * 1.5)

    stat, _, _ = stats_table(plot_df, assay)
    p_w = stat.loc[stat["metric"] == "V-channel width (bp)", "WRS_p"].iloc[0]
    p_e = stat.loc[stat["metric"] == "log2(V-in/V-out)", "WRS_p"].iloc[0]
    auc_w = stat.loc[stat["metric"] == "V-channel width (bp)", "AUC"].iloc[0]
    auc_e = stat.loc[stat["metric"] == "log2(V-in/V-out)", "AUC"].iloc[0]

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.add_patch(
        plt.Rectangle(
            (thr_w, thr_e),
            xmax - thr_w,
            ymax - thr_e,
            facecolor="#BFD3E6",
            alpha=0.35,
            zorder=0,
        )
    )
    ax.axvline(thr_w, ls=(0, (5, 4)), lw=1.4, color="#666666", zorder=1)
    ax.axhline(thr_e, ls=(0, (5, 4)), lw=1.4, color="#666666", zorder=1)
    ax.scatter(
        bias["width"],
        bias["log2FC"],
        s=56,
        facecolors="none",
        edgecolors="#FF6A6A",
        linewidths=1.5,
        zorder=3,
        label="enzyme bias",
    )
    ax.scatter(
        tf["width"],
        tf["log2FC"],
        s=52,
        c="#104e8b",
        marker="o",
        zorder=3,
        label="TF",
    )

    stat_text = ax.text(
        0.01,
        0.98,
        (
            "Wilcoxon rank-sum (TF vs bias)\n"
            f"width: p = {format_p(p_w)} (AUC {auc_w:.2f})\n"
            f"Vin / Vout: p = {format_p(p_e)} (AUC {auc_e:.2f})"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("V-channel width (bp)", fontsize=15)
    ax.set_ylabel("log2(Vin / Vout)", fontsize=15)
    ax.set_title(assay, fontsize=17, fontweight="bold", pad=8)
    ax.grid(True, color="#E6E6E6", lw=1.0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)
        spine.set_color("black")
    legend = ax.legend(frameon=False, loc="lower right", fontsize=8.8)
    fig.tight_layout()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    reserved = [
        stat_text.get_window_extent(renderer).padded(4),
        legend.get_window_extent(renderer).padded(4),
    ]
    reserved.extend(point_boxes(ax, plot_df, pad=8.0))
    positions = TF_LABEL_POSITIONS.get(assay, {})
    manual_tf = tf[tf["TF"].isin(positions)].copy()
    label_boxes = place_tf_labels_data(
        fig,
        ax,
        manual_tf.sort_values(["log2FC", "width"], ascending=[False, False]),
        positions,
        "#104e8b",
        7.8,
        reserved,
    )
    fallback_tf = tf[~tf["TF"].isin(positions)].copy()
    if not fallback_tf.empty:
        label_boxes = place_labels(
            fig,
            ax,
            fallback_tf.sort_values(["log2FC", "width"], ascending=[False, False]),
            LABEL_OFFSETS.get(assay, {}),
            "#104e8b",
            7.8,
            TF_LABEL_CANDIDATES,
            label_boxes,
            connector_color="#104e8b",
            min_offset_points=36.0,
        )
    place_labels(
        fig,
        ax,
        bias.sort_values(["log2FC", "width"], ascending=[False, True]),
        {},
        "#D24B4B",
        8.5,
        BIAS_LABEL_CANDIDATES,
        label_boxes,
    )
    fig.savefig(OUTDIR / f"methodB_discriminant_scatter_best_{assay}.png", dpi=300)
    fig.savefig(OUTDIR / f"methodB_discriminant_scatter_best_{assay}.pdf")
    fig.savefig(OUTDIR / f"methodB_discriminant_scatter_best_{assay}.svg")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    comp = motif_composition()
    comp = comp[~comp["TF"].isin(EXCLUDE_TFS)].copy()
    if comp.empty:
        raise RuntimeError(f"No MEME motif matrices found in {MEME_DIR}")

    tf_all = {assay: read_score_table(TF_TABLE[assay], "TF", assay) for assay in ASSAYS}
    bias_all = {
        assay: read_score_table(BIAS_TABLE[assay], "bias", assay) for assay in ASSAYS
    }
    ranked, selected, common_selected = score_and_select(tf_all, bias_all, comp)

    stats_rows = []
    for assay in ASSAYS:
        selected_tf = selected.loc[selected["assay"] == assay, "TF"].tolist()
        plot_df = pd.concat(
            [
                bias_all[assay][bias_all[assay]["TF"].isin(BIAS_KEEP[assay])],
                tf_all[assay][tf_all[assay]["TF"].isin(selected_tf)].merge(
                    selected[
                        selected["assay"] == assay
                    ][
                        [
                            "TF",
                            "family",
                            "rank_score",
                            "family_score",
                            "assay_rank",
                            "family_representative",
                            "GC",
                            "AT",
                            "base_range",
                            "base_class",
                            "four_base_class",
                        ]
                    ],
                    on="TF",
                    how="left",
                ),
            ],
            ignore_index=True,
        )
        plot_df.to_csv(OUTDIR / f"plot_data_best_{assay}.csv", index=False)
        st, thr_w, thr_e = stats_table(plot_df, assay)
        stats_rows.append(st)
        plot_assay(plot_df, assay, thr_w, thr_e)

    ranked.to_csv(OUTDIR / "ranked_candidate_TFs_by_assay.csv", index=False)
    selected.to_csv(OUTDIR / "selected_TFs_by_assay.csv", index=False)
    at_rich_sep = ranked[
        (ranked["base_class"] == "AT-rich")
        & (ranked["delta_width"] > 0)
        & (ranked["delta_log2FC"] > 0)
    ].sort_values(["assay", "rank_score"], ascending=[True, False])
    at_rich_sep.to_csv(OUTDIR / "AT_rich_separating_candidates.csv", index=False)
    comp.sort_values("GC").to_csv(OUTDIR / "motif_base_composition.csv", index=False)
    pd.concat(stats_rows, ignore_index=True).to_csv(
        OUTDIR / "methodB_separation_stats_best_common.csv", index=False
    )

    summary_path = OUTDIR / "selection_summary.txt"
    with summary_path.open("w", newline="") as out:
        out.write(f"Top N per assay: {TOP_N_PER_ASSAY}\n")
        out.write("Common anchors requested: " + ", ".join(COMMON_ANCHORS) + "\n")
        out.write("Excluded TFs: " + ", ".join(sorted(EXCLUDE_TFS)) + "\n")
        out.write(
            "Family de-duplication groups: "
            + "; ".join(
                f"{family}={','.join(members)}"
                for family, members in TF_FAMILY_GROUPS.items()
            )
            + "\n"
        )
        out.write(
            "Common TFs across all selected panels: "
            + ", ".join(sorted(common_selected))
            + f" (n={len(common_selected)})\n\n"
        )
        out.write("Selected TFs by assay:\n")
        for assay in ASSAYS:
            names = selected.loc[selected["assay"] == assay, "TF"].tolist()
            out.write(f"- {assay}: {', '.join(names)}\n")
        out.write("\nAT-rich candidates with positive separation from bias:\n")
        if at_rich_sep.empty:
            out.write("- none among TFs with local MEME motif matrices\n")
        else:
            for assay in ASSAYS:
                rows = at_rich_sep[at_rich_sep["assay"] == assay]
                if rows.empty:
                    continue
                compact = [
                    (
                        f"{r.TF}(rank={int(r.assay_rank)}, "
                        f"width={r.width:.2f}, log2={r.log2FC:.2f})"
                    )
                    for r in rows.itertuples()
                ]
                out.write(f"- {assay}: {', '.join(compact)}\n")
        out.write("\nLocal MEME motif composition of selected unique TFs:\n")
        unique_selected = sorted(set(selected["TF"]))
        comp_sel = comp[comp["TF"].isin(unique_selected)].sort_values("GC")
        atcg_balanced = comp_sel.loc[
            comp_sel["four_base_class"] == "ATCG-balanced", "TF"
        ].tolist()
        out.write(
            "ATCG-balanced selected TFs (base_range <= 0.15): "
            + (", ".join(atcg_balanced) if atcg_balanced else "none")
            + "\n"
        )
        for _, row in comp_sel.iterrows():
            out.write(
                f"- {row['TF']}: GC={row['GC']:.3f}, AT={row['AT']:.3f}, "
                f"base_range={row['base_range']:.3f}, "
                f"class={row['base_class']}/{row['four_base_class']}, "
                f"consensus={row['consensus']}\n"
            )
        missing = sorted(set(unique_selected) - set(comp["TF"]))
        if missing:
            out.write(
                "\nSelected TFs without local MEME matrices, not included in motif "
                f"composition audit: {', '.join(missing)}\n"
            )

    print(f"Wrote outputs to: {OUTDIR}")
    print("Common TFs:", ", ".join(sorted(common_selected)))


if __name__ == "__main__":
    main()
