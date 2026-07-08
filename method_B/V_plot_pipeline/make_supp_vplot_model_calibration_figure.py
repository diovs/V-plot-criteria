#!/usr/bin/env python3
"""Create the revised supplementary V-plot modeling overview figure."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


BLUE = "#1f5a93"
ORANGE = "#d46b32"
RED = "#bd3f3f"
DARK = "#222222"
MID = "#5b6670"
GRID = "#d7dde3"
PALE_BLUE = "#eaf2f8"
PALE_GREEN = "#e9f4ef"
PALE_ORANGE = "#f7ede5"
ASSAYS = ["loMNase", "DNase", "ATAC"]
P_FLOOR = np.nextafter(0, 1)
PLOT_X_CAP = 60.0


def clean_name(value: str) -> str:
    name = Path(str(value)).stem
    for suffix in ["_ATAC_fragL_dist", "_DNase_fragL_dist", "_loMNase_fragL_dist", "_fragL_dist"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def read_csv_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input file was not found: {path}")
    return pd.read_csv(path)


def read_tsv_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required input file was not found: {path}")
    return pd.read_csv(path, sep="\t")


def assay_paths(method_root: Path) -> dict[str, Path]:
    data_dir = method_root / "v_apex_scores" / "methodB_scatter_best_assay_specific"
    return {assay: data_dir / f"plot_data_best_{assay}.csv" for assay in ASSAYS}


def threshold_table(method_root: Path) -> pd.DataFrame:
    path = method_root / "fig_cutoff_derivation_thresholds.csv"
    df = read_csv_required(path)
    required = {"assay", "feature", "cutoff", "margin", "max_bias", "min_TF"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")
    return df


def get_cutoff(thr: pd.DataFrame, assay: str, feature: str) -> float:
    sub = thr[(thr["assay"].astype(str) == assay) & thr["feature"].astype(str).str.contains(feature, case=False, na=False)]
    if sub.empty:
        raise ValueError(f"Missing {feature} cutoff for {assay}")
    return float(sub.iloc[0]["cutoff"])


def get_full_margin(thr: pd.DataFrame, assay: str, feature: str) -> float:
    sub = thr[(thr["assay"].astype(str) == assay) & thr["feature"].astype(str).str.contains(feature, case=False, na=False)]
    if sub.empty:
        raise ValueError(f"Missing {feature} margin for {assay}")
    return float(sub.iloc[0]["margin"])


def load_plot_data(method_root: Path) -> pd.DataFrame:
    frames = []
    for assay, path in assay_paths(method_root).items():
        df = read_csv_required(path)
        required = {"assay", "TF", "type", "width", "log2FC", "motif"}
        missing = sorted(required.difference(df.columns))
        if missing:
            raise ValueError(f"Missing columns in {path}: {missing}")
        df = df.copy()
        df["assay"] = assay
        df["group_name"] = df["TF"].astype(str)
        df["group_type"] = np.where(df["type"].astype(str).str.lower() == "bias", "matched cleavage-bias controls", "TF motifs")
        df["width"] = pd.to_numeric(df["width"], errors="coerce")
        df["inside_v_enrichment"] = pd.to_numeric(df["log2FC"], errors="coerce")
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out = out[np.isfinite(out["width"]) & np.isfinite(out["inside_v_enrichment"])].copy()
    if out.empty:
        raise ValueError("No usable plot-data rows after parsing width and enrichment.")
    return out


def load_lrt_table(method_root: Path) -> pd.DataFrame:
    tf_frames = []
    for assay in ASSAYS:
        for path in [
            method_root / f"{assay}_TF_apex_scores.tsv",
            method_root / "v_apex_scores" / f"{assay}_TF_apex_scores_selected.tsv",
            method_root / "TF_vplot" / "from_server_three_assays" / "apex_scores" / f"{assay}_TF_apex_scores.tsv",
        ]:
            if not path.exists():
                continue
            df = read_tsv_required(path)
            required = {"motif", "p_chi2"}
            missing = sorted(required.difference(df.columns))
            if missing:
                raise ValueError(f"Missing columns in {path}: {missing}")
            group_name = df["TF"].astype(str) if "TF" in df.columns else df["motif"].astype(str).map(clean_name)
            tf_frames.append(
                pd.DataFrame(
                    {
                        "assay": assay,
                        "motif": df["motif"].astype(str),
                        "group_name": group_name,
                        "group_type": "TF motifs",
                        "lrt_pvalue": pd.to_numeric(df["p_chi2"], errors="coerce"),
                        "lrt_qvalue": pd.to_numeric(df["q_chi2"], errors="coerce") if "q_chi2" in df.columns else np.nan,
                        "source_priority": len(tf_frames),
                    }
                )
            )
    curated_path = method_root / "TF_vplot" / "from_server_three_assays" / "apex_scores" / "curated_TF_apex_scores_all.tsv"
    if curated_path.exists():
        curated = read_tsv_required(curated_path)
        if {"motif", "assay", "TF", "p_chi2"}.issubset(curated.columns):
            tf_frames.append(
                pd.DataFrame(
                    {
                        "assay": curated["assay"].astype(str),
                        "motif": curated["motif"].astype(str),
                        "group_name": curated["TF"].astype(str),
                        "group_type": "TF motifs",
                        "lrt_pvalue": pd.to_numeric(curated["p_chi2"], errors="coerce"),
                        "lrt_qvalue": pd.to_numeric(curated["q_chi2"], errors="coerce") if "q_chi2" in curated.columns else np.nan,
                        "source_priority": len(tf_frames),
                    }
                )
            )
    if not tf_frames:
        raise FileNotFoundError("No TF apex-score tables with LRT p-values were found.")
    tf_out = pd.concat(tf_frames, ignore_index=True)
    tf_out = tf_out.sort_values("source_priority").drop_duplicates(["assay", "motif"], keep="first")
    tf_out = tf_out.drop(columns=["source_priority"])

    bias_frames = []
    for assay in ASSAYS:
        path = method_root / "sequence_bias" / f"{assay}_Bias_3bp_apex_scores.tsv"
        df = read_tsv_required(path)
        required = {"motif", "p_chi2"}
        missing = sorted(required.difference(df.columns))
        if missing:
            raise ValueError(f"Missing columns in {path}: {missing}")
        out = pd.DataFrame(
            {
                "assay": assay,
                "motif": df["motif"].astype(str),
                "group_name": df["motif"].astype(str).map(lambda x: clean_name(x).replace(f"{assay}_Bias_", "")),
                "group_type": "matched cleavage-bias controls",
                "lrt_pvalue": pd.to_numeric(df["p_chi2"], errors="coerce"),
                "lrt_qvalue": pd.to_numeric(df["q_chi2"], errors="coerce") if "q_chi2" in df.columns else np.nan,
            }
        )
        bias_frames.append(out)
    return pd.concat([tf_out, *bias_frames], ignore_index=True)


def make_panel_b_table(method_root: Path, out_csv: Path) -> pd.DataFrame:
    plot_df = load_plot_data(method_root)
    lrt = load_lrt_table(method_root)
    thr = threshold_table(method_root)
    merged = plot_df.merge(lrt[["assay", "motif", "lrt_pvalue", "lrt_qvalue"]], on=["assay", "motif"], how="left")
    if merged["lrt_pvalue"].isna().any():
        missing = merged.loc[merged["lrt_pvalue"].isna(), ["assay", "motif"]].to_dict("records")
        raise ValueError(f"Could not map LRT p-values for rows: {missing[:10]}")
    use_qvalue = merged["lrt_qvalue"].notna().all()

    rows = []
    for row in merged.itertuples(index=False):
        width_cutoff = get_cutoff(thr, row.assay, "width")
        enrichment_cutoff = get_cutoff(thr, row.assay, "enrichment")
        lrt_value = float(row.lrt_qvalue) if use_qvalue else float(row.lrt_pvalue)
        p = max(lrt_value, P_FLOOR)
        width_ratio = float(row.width) / width_cutoff if width_cutoff else np.nan
        enrichment_ratio = float(row.inside_v_enrichment) / enrichment_cutoff if enrichment_cutoff else np.nan
        passes_rule = (width_ratio >= 1.0) and (enrichment_ratio >= 1.0)
        rows.append(
            {
                "assay": row.assay,
                "group_type": row.group_type,
                "group_name": row.group_name,
                "motif": row.motif,
                "width": float(row.width),
                "inside_v_enrichment": float(row.inside_v_enrichment),
                "lrt_pvalue": float(row.lrt_pvalue),
                "lrt_qvalue": float(row.lrt_qvalue) if pd.notna(row.lrt_qvalue) else np.nan,
                "lrt_statistic": "q-value" if use_qvalue else "P-value",
                "width_cutoff": width_cutoff,
                "enrichment_cutoff": enrichment_cutoff,
                "width_ratio": width_ratio,
                "enrichment_ratio": enrichment_ratio,
                "minus_log10_lrt": -np.log10(p),
                "minus_log10_lrt_display": min(-np.log10(p), PLOT_X_CAP),
                "passes_two_feature_rule": passes_rule,
                "classification": "TF-footprint-like" if passes_rule else "cleavage-bias-like / below cut-off",
            }
        )
    out = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    return out


def midgap(values: pd.DataFrame, feature: str) -> tuple[float, float, float, float]:
    bias = values.loc[values["group_type"] == "matched cleavage-bias controls", feature]
    tf = values.loc[values["group_type"] == "TF motifs", feature]
    if bias.empty or tf.empty:
        raise ValueError(f"Cannot compute mid-gap for {feature}: missing TF or bias rows.")
    max_bias = float(bias.max())
    min_tf = float(tf.min())
    cutoff = (max_bias + min_tf) / 2.0
    margin = (min_tf - max_bias) / 2.0
    return max_bias, min_tf, cutoff, margin


def make_panel_d_table(method_root: Path, out_csv: Path) -> tuple[pd.DataFrame, list[str]]:
    data = load_plot_data(method_root)
    thr = threshold_table(method_root)
    warnings = []
    rows = []
    for assay in ASSAYS:
        assay_df = data[data["assay"] == assay].copy()
        full_width = get_cutoff(thr, assay, "width")
        full_enrich = get_cutoff(thr, assay, "enrichment")
        for left_idx, held in assay_df.iterrows():
            remaining = assay_df.drop(index=left_idx)
            _wb, _wt, loo_width, width_margin = midgap(remaining, "width")
            _eb, _et, loo_enrich, enrich_margin = midgap(remaining, "inside_v_enrichment")
            if width_margin <= 0:
                warnings.append(f"{assay} width non-positive LOO margin after leaving out {held['group_name']}")
            if enrich_margin <= 0:
                warnings.append(f"{assay} enrichment non-positive LOO margin after leaving out {held['group_name']}")
            pred_pass = (float(held["width"]) >= loo_width) and (float(held["inside_v_enrichment"]) >= loo_enrich)
            expected = "TF-footprint-like" if held["group_type"] == "TF motifs" else "cleavage-bias-like"
            predicted = "TF-footprint-like" if pred_pass else "cleavage-bias-like"
            correct = expected == predicted
            common = {
                "assay": assay,
                "left_out_group": held["group_name"],
                "left_out_group_type": held["group_type"],
                "held_out_expected_class": expected,
                "held_out_predicted_class": predicted,
                "held_out_correct": bool(correct),
            }
            rows.append(
                {
                    **common,
                    "feature": "width",
                    "full_cutoff": full_width,
                    "loo_cutoff": loo_width,
                    "loo_over_full": loo_width / full_width if full_width else np.nan,
                    "margin": width_margin,
                }
            )
            rows.append(
                {
                    **common,
                    "feature": "inside_v_enrichment",
                    "full_cutoff": full_enrich,
                    "loo_cutoff": loo_enrich,
                    "loo_over_full": loo_enrich / full_enrich if full_enrich else np.nan,
                    "margin": enrich_margin,
                }
            )
    out = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    return out, warnings


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom", ha="left", color=DARK)


def draw_panel_a(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("V-cone geometry and maximum-likelihood apex estimation", loc="left", fontsize=11.0, fontweight="bold", pad=6)

    vax = ax.inset_axes([0.04, 0.16, 0.44, 0.73])
    rng = np.random.default_rng(42)
    b = 20.0
    x_bg = rng.uniform(-60, 60, 900)
    y_bg = rng.uniform(0, 100, 900)
    vax.scatter(x_bg, y_bg, s=4, c="#bfc5cc", alpha=0.55, edgecolors="none")
    y_sig = rng.uniform(b + 2, 100, 520)
    half = np.maximum((y_sig - b) * 0.5, 1.0)
    x_sig = np.clip(rng.normal(0, half / 2.6), -half, half)
    vax.scatter(x_sig, y_sig, s=5, c=BLUE, alpha=0.72, edgecolors="none")
    xs = np.linspace(0, 50, 100)
    yy = b + 2 * xs
    keep = yy <= 100
    vax.plot(xs[keep], yy[keep], color=BLUE, lw=2.0)
    vax.plot(-xs[keep], yy[keep], color=BLUE, lw=2.0)
    vax.scatter([0], [b], marker="*", s=120, color=RED, edgecolor="white", linewidth=0.5, zorder=5)
    vax.annotate(
        "apex = (a, b)",
        xy=(0, b),
        xytext=(24, 10),
        arrowprops=dict(arrowstyle="->", lw=1.2, color=RED),
        fontsize=8.5,
        color=RED,
        fontweight="bold",
    )
    vax.text(-55, 32, "local\nbackground", fontsize=8.4, color=MID, ha="left")
    vax.text(0, 63, "V-like signal\nfraction pi", fontsize=9.0, color=BLUE, ha="center", fontweight="bold")
    vax.text(-58, 19, "b = fitted V-channel width", fontsize=8.0, color=MID, ha="left")
    vax.text(12, 79, "w(y) = (y - b) / 2", fontsize=8.0, color=DARK)
    vax.set_xlim(-60, 60)
    vax.set_ylim(0, 100)
    vax.set_aspect("equal", adjustable="box")
    vax.set_xlabel("distance from motif center (bp)", fontsize=8.5, labelpad=2)
    vax.set_ylabel("fragment length (bp)", fontsize=8.5, labelpad=2)
    vax.set_xticks([-60, -30, 0, 30, 60])
    vax.set_yticks([0, 20, 40, 60, 80, 100])
    vax.tick_params(labelsize=7.5, length=2.6, pad=1)
    vax.grid(True, color=GRID, lw=0.4, alpha=0.55)
    for spine in vax.spines.values():
        spine.set_linewidth(0.85)
        spine.set_color("#4c5560")

    rax = ax.inset_axes([0.59, 0.54, 0.36, 0.30])
    x = np.linspace(-20, 20, 140)
    clear = -0.028 * (x - 1.5) ** 2 + 1.05
    weak = 0.09 * np.sin(x / 5.2) - 0.0025 * x**2 - 0.10
    rax.plot(x, clear, color=BLUE, lw=1.9, label="clear apex")
    rax.plot(x, weak, color=ORANGE, lw=1.7, label="weak / flat")
    rax.axvline(1.5, color=RED, lw=1.0, ls="--")
    rax.text(0.96, 0.86, "selected apex =\nmaximum likelihood", transform=rax.transAxes, ha="right", va="top", fontsize=7.0)
    rax.set_xlabel("candidate apex position", fontsize=7.5, labelpad=1)
    rax.set_ylabel("profile\nlog-likelihood", fontsize=7.0, labelpad=1)
    rax.tick_params(labelsize=7.0, length=2.4, pad=1)
    rax.grid(True, color=GRID, lw=0.4, alpha=0.75)
    rax.legend(frameon=False, fontsize=7.0, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.02))

    ax.text(
        0.58,
        0.25,
        "Apex and width are estimated by the V-plot model.\n"
        "LRT tests whether a V-like structure is detectable.\n"
        "Final calls require calibrated width and enrichment cut-offs.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.9,
        color=DARK,
        linespacing=1.32,
        bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor="#cbd2d8", alpha=0.95),
    )


def draw_panel_b(ax: plt.Axes, panel_b: pd.DataFrame) -> None:
    ax.set_axis_off()
    ax.set_title("LRT significance is insufficient for footprint-like classification", loc="left", fontsize=10.2, fontweight="bold", pad=6)
    left = ax.inset_axes([0.07, 0.19, 0.40, 0.68])
    right = ax.inset_axes([0.56, 0.19, 0.40, 0.68])
    lrt_label = str(panel_b["lrt_statistic"].iloc[0]) if "lrt_statistic" in panel_b.columns else "P-value"

    def draw_ratio(plot_ax: plt.Axes, ratio_col: str, ylabel: str, title: str) -> None:
        for assay, marker in zip(ASSAYS, ["o", "s", "^"]):
            sub = panel_b[panel_b["assay"] == assay]
            for group, color, face in [
                ("TF motifs", BLUE, BLUE),
                ("matched cleavage-bias controls", "#ff6a6a", "white"),
            ]:
                pts = sub[sub["group_type"] == group]
                plot_ax.scatter(
                    pts["minus_log10_lrt_display"],
                    pts[ratio_col],
                    s=26,
                    marker=marker,
                    facecolor=face,
                    edgecolor=color,
                    lw=1.1,
                    alpha=0.92,
                    label=f"{assay} {'TF' if group == 'TF motifs' else 'control'}",
                )
        plot_ax.axhline(1.0, color=DARK, ls="--", lw=1.0)
        plot_ax.set_xlabel(f"-log10(LRT {lrt_label})", fontsize=7.8, labelpad=1)
        plot_ax.set_ylabel(ylabel, fontsize=7.8, labelpad=1)
        plot_ax.set_title(title, fontsize=8.5, loc="left", fontweight="bold", pad=3)
        plot_ax.set_xlim(-1, PLOT_X_CAP + 3)
        plot_ax.tick_params(labelsize=7.0, length=2.4, pad=1)
        plot_ax.grid(True, color=GRID, lw=0.42, alpha=0.75)
        plot_ax.text(
            0.96,
            0.14,
            "LRT+ but below\ncut-off",
            transform=plot_ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.0,
            color=ORANGE,
            fontweight="bold",
        )

    draw_ratio(left, "width_ratio", "fitted width / width cut-off", "Width requirement")
    draw_ratio(right, "enrichment_ratio", "inside-V enrichment / enrichment cut-off", "Enrichment requirement")
    handles = [
        plt.Line2D([], [], marker="o", color=BLUE, ls="none", ms=5, label="TF motifs"),
        plt.Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor="#ff6a6a", color="#ff6a6a", ls="none", ms=5, label="matched controls"),
        plt.Line2D([], [], marker="o", color=MID, ls="none", ms=4.5, label="loMNase"),
        plt.Line2D([], [], marker="s", color=MID, ls="none", ms=4.5, label="DNase"),
        plt.Line2D([], [], marker="^", color=MID, ls="none", ms=4.5, label="ATAC"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7.0, ncol=5, loc="lower center", bbox_to_anchor=(0.52, 0.00), columnspacing=1.0, handletextpad=0.35)


def draw_panel_c(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Reproducible V-plot analysis workflow", loc="left", fontsize=10.2, fontweight="bold", pad=6)
    boxes = [
        ("Inputs", "fragment BED\nTF motif sites\nmatched cleavage-bias controls"),
        ("Build V-plot tables", "fragment midpoint + length\ndistance to anchor center"),
        ("Fit V-plot model", "estimate apex\nfitted width\ninside-V enrichment\nLRT evidence"),
        ("Calibrate cut-offs", "mid-gap cut-offs from\nTF motifs and matched controls"),
        ("Final classification", "TF-footprint-like if\nwidth >= cut-off\nand enrichment >= cut-off"),
    ]
    x_positions = [0.02, 0.215, 0.41, 0.605, 0.80]
    box_w = 0.16
    box_h = 0.42
    y = 0.38
    for i, ((title, desc), x) in enumerate(zip(boxes, x_positions)):
        face = PALE_GREEN if i == 4 else (PALE_ORANGE if i == 3 else PALE_BLUE)
        box = FancyBboxPatch(
            (x, y),
            box_w,
            box_h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=0.9,
            edgecolor="#6c7882",
            facecolor=face,
        )
        ax.add_patch(box)
        ax.text(x + box_w / 2, y + box_h - 0.08, title, ha="center", va="center", fontsize=8.0, fontweight="bold", color=DARK)
        ax.text(x + box_w / 2, y + 0.15, desc, ha="center", va="center", fontsize=6.8, color=DARK, linespacing=1.18)
        if i < len(boxes) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + box_w + 0.006, y + box_h / 2),
                    (x_positions[i + 1] - 0.006, y + box_h / 2),
                    arrowstyle="-|>",
                    mutation_scale=10,
                    lw=0.9,
                    color=MID,
                )
            )
    ax.text(
        0.5,
        0.14,
        "LRT is reported as evidence for V-like structure; the final call uses the calibrated two-feature rule.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.0,
        color=DARK,
    )


def draw_panel_d(ax: plt.Axes, panel_d: pd.DataFrame) -> None:
    ax.set_axis_off()
    ax.set_title("Leave-one-out calibration shows stable cut-offs across assays", loc="left", fontsize=10.2, fontweight="bold", pad=6)
    left = ax.inset_axes([0.06, 0.25, 0.39, 0.49])
    right = ax.inset_axes([0.54, 0.25, 0.39, 0.49])

    def draw_feature(plot_ax: plt.Axes, feature: str, title: str) -> None:
        sub = panel_d[panel_d["feature"] == feature].copy()
        x_base = {assay: i for i, assay in enumerate(ASSAYS)}
        rng = np.random.default_rng(9)
        for group, color, face in [
            ("TF motifs", BLUE, BLUE),
            ("matched cleavage-bias controls", "#ff6a6a", "white"),
        ]:
            pts = sub[sub["left_out_group_type"] == group]
            xs = np.array([x_base[a] for a in pts["assay"]], dtype=float) + rng.uniform(-0.10, 0.10, len(pts))
            plot_ax.scatter(xs, pts["loo_over_full"], s=24, facecolor=face, edgecolor=color, lw=1.0, alpha=0.90)
        for i, assay in enumerate(ASSAYS):
            plot_ax.scatter([i], [1.0], marker="D", s=42, color=DARK, zorder=4)
        plot_ax.axhline(1.0, color=DARK, lw=1.0, ls="--")
        plot_ax.axhspan(0.9, 1.1, color=PALE_GREEN, alpha=0.5, zorder=0)
        plot_ax.set_xticks(range(len(ASSAYS)))
        plot_ax.set_xticklabels(ASSAYS, fontsize=7.2)
        plot_ax.set_ylabel("LOO cut-off / full-data cut-off", fontsize=7.6, labelpad=1)
        plot_ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=3)
        plot_ax.set_ylim(max(0, sub["loo_over_full"].min() * 0.88), max(1.25, sub["loo_over_full"].max() * 1.12))
        plot_ax.tick_params(axis="y", labelsize=7.0, length=2.4, pad=1)
        plot_ax.grid(True, axis="y", color=GRID, lw=0.42, alpha=0.75)

    draw_feature(left, "width", "V-channel width cut-off")
    draw_feature(right, "inside_v_enrichment", "Inside-V enrichment cut-off")

    uniq = panel_d.drop_duplicates(["assay", "left_out_group"])
    rows = []
    for assay in ASSAYS:
        sub = uniq[uniq["assay"] == assay]
        acc = float(sub["held_out_correct"].mean()) if len(sub) else np.nan
        rows.append([assay, str(len(sub)), f"{acc:.2f}"])
    table_ax = ax.inset_axes([0.37, 0.02, 0.26, 0.17])
    table_ax.set_axis_off()
    table = table_ax.table(cellText=rows, colLabels=["assay", "LOO n", "accuracy"], loc="center", cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7.0)
    table.scale(1.0, 1.18)
    for (r, _c), cell in table.get_celld().items():
        cell.set_linewidth(0.45)
        cell.set_edgecolor("#8b959e")
        if r == 0:
            cell.set_facecolor("#edf2f6")
            cell.set_text_props(weight="bold")
    handles = [
        plt.Line2D([], [], marker="o", color=BLUE, ls="none", ms=5, label="left-out TF motif"),
        plt.Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor="#ff6a6a", color="#ff6a6a", ls="none", ms=5, label="left-out matched control"),
        plt.Line2D([], [], marker="D", color=DARK, ls="none", ms=5, label="full-data cut-off"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7.0, loc="upper center", ncol=3, bbox_to_anchor=(0.52, 0.94), columnspacing=1.0)


def write_notes(out_path: Path, panel_b_path: Path, panel_d_path: Path, warnings: list[str]) -> Path:
    notes = out_path.with_name("supp_vplot_model_calibration_overview_notes.md")
    text = [
        "# Supplementary V-plot calibration figure notes",
        "",
        "## Data sources",
        "",
        "- Panel B diagnostic table: `" + str(panel_b_path) + "`",
        "- Panel D leave-one-out table: `" + str(panel_d_path) + "`",
        "- Panel B column mapping: width=`width`; inside-V enrichment=`log2FC` -> `inside_v_enrichment`; LRT statistic=`q_chi2` when available, otherwise `p_chi2`; group type=`type`; assay=`assay`.",
        "- Panel D cut-offs were recomputed by leaving out each TF motif group or matched cleavage-bias control group from the final per-assay plot-data tables.",
        "",
        "## Warnings",
        "",
    ]
    text.extend(["- " + w for w in warnings] if warnings else ["- None."])
    notes.write_text("\n".join(text) + "\n", encoding="utf-8")
    return notes


def build_figure(method_root: Path, output_dir: Path, out_prefix: str, dpi: int) -> tuple[list[Path], Path, pd.DataFrame, pd.DataFrame, list[str]]:
    panel_b_path = output_dir / "panelB_lrt_effectsize_diagnostic.csv"
    panel_d_path = output_dir / "panelD_loo_cutoff_stability.csv"
    panel_b = make_panel_b_table(method_root, panel_b_path)
    panel_d, warnings = make_panel_d_table(method_root, panel_d_path)

    fig = plt.figure(figsize=(8.8, 13.5))
    gs = GridSpec(4, 1, figure=fig, height_ratios=[1.12, 1.05, 0.80, 1.10], hspace=0.48)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[3, 0])
    draw_panel_a(ax_a)
    draw_panel_b(ax_b, panel_b)
    draw_panel_c(ax_c)
    draw_panel_d(ax_d, panel_d)
    for ax, label in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D")]:
        panel_label(ax, label)
    fig.suptitle(
        "Supplementary Fig. X | Calibrated V-plot modeling separates TF-footprint-like signals\n"
        "from cleavage-bias-like V patterns.",
        fontsize=12.2,
        fontweight="bold",
        y=0.988,
    )
    fig.subplots_adjust(left=0.065, right=0.985, top=0.925, bottom=0.040)

    written = []
    for ext in ["png", "pdf", "svg"]:
        path = output_dir / f"{out_prefix}.{ext}"
        fig.savefig(path, dpi=dpi if ext == "png" else None, facecolor="white", bbox_inches="tight")
        written.append(path)
    plt.close(fig)
    notes = write_notes(output_dir / f"{out_prefix}.png", panel_b_path, panel_d_path, warnings)
    return written, notes, panel_b, panel_d, warnings


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    default_method_root = here.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method-root", type=Path, default=default_method_root)
    parser.add_argument("--output-dir", type=Path, default=here)
    parser.add_argument("--out-prefix", default="supp_vplot_model_calibration_overview")
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    written, notes, _panel_b, _panel_d, warnings = build_figure(
        method_root=args.method_root.resolve(),
        output_dir=args.output_dir.resolve(),
        out_prefix=args.out_prefix,
        dpi=args.dpi,
    )
    for path in written:
        print(f"wrote {path}")
    print(f"wrote {notes}")
    for warning in warnings:
        print(f"WARNING: {warning}")


if __name__ == "__main__":
    main()
