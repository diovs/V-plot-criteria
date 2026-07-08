# -*- coding: utf-8 -*-
"""Schematic showing HOW each cut-off in the table is derived (mid-gap rule),
using the real per-motif data. For each assay and feature:
   cut-off = midpoint of the gap between the highest bias and the lowest TF
   margin  = half the gap width.
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12})
C_TF = "#104e8b"; C_BIAS = "#FF6A6A"
ROOT = Path(__file__).resolve().parent
DATADIR = ROOT / "data"
OUTBASE = ROOT / "fig_cutoff_derivation"
OUTCSV = OUTBASE.with_name(OUTBASE.name + "_thresholds.csv")

def load(assay):
    w_b, e_b, w_t, e_t = [], [], [], []
    with (DATADIR / f"plot_data_best_{assay}.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            w = float(row["width"])
            e = float(row["log2FC"])
            if row["type"] == "bias": w_b.append(w); e_b.append(e)
            else: w_t.append(w); e_t.append(e)
    return dict(w_bias=np.array(w_b), e_bias=np.array(e_b),
               w_tf=np.array(w_t), e_tf=np.array(e_t))

assays = ["loMNase", "DNase", "ATAC"]
D = {a: load(a) for a in assays}
rng = np.random.default_rng(1)
threshold_rows = []

def midgap(assay, metric):
    bias = D[assay][f"{metric}_bias"]
    tf = D[assay][f"{metric}_tf"]
    maxb = bias.max()
    mintf = tf.min()
    mid = (maxb + mintf) / 2
    margin = (mintf - maxb) / 2
    return maxb, mintf, mid, margin

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.6, 6.2))
lane_y = {"loMNase": 2.4, "DNase": 1.4, "ATAC": 0.4}
HH = 0.26   # lane half-height

def draw(ax, metric, title, xlabel, unit):
    for a in assays:
        yc = lane_y[a]
        bias = D[a][f"{metric}_bias"]; tf = D[a][f"{metric}_tf"]
        maxb, mintf, mid, margin = midgap(a, metric)
        # gap shading
        ax.add_patch(plt.Rectangle((maxb, yc - HH), mintf - maxb, 2 * HH,
                                   facecolor="#BFD3E6", alpha=0.35, zorder=0))
        # points (jittered in lane)
        ax.scatter(bias, yc + rng.uniform(-HH*0.55, HH*0.55, len(bias)),
                   s=55, facecolors="none", edgecolors=C_BIAS, linewidths=1.7, zorder=3)
        ax.scatter(tf, yc + rng.uniform(-HH*0.55, HH*0.55, len(tf)),
                   s=45, c=C_TF, zorder=3)
        # cut-off line
        ax.plot([mid, mid], [yc - HH - 0.06, yc + HH + 0.06], color="#222222",
                ls="--", lw=1.8, zorder=4)
        ax.text(mid, yc + HH + 0.13, f"cut-off = {mid:.2f}{unit}", ha="center",
                va="bottom", fontsize=11, fontweight="bold", color="#222222")
        # margin double-arrow (from cut-off to min TF = one margin)
        ax.add_patch(FancyArrowPatch((mid, yc - HH - 0.02), (mintf, yc - HH - 0.02),
                     arrowstyle="<->", mutation_scale=11, lw=1.3, color="#555555", zorder=4))
        ax.text((mid + mintf) / 2, yc - HH - 0.16, f"margin {margin:.2f}", ha="center",
                va="top", fontsize=9.5, color="#555555")
    ax.set_yticks(list(lane_y.values())); ax.set_yticklabels(list(lane_y.keys()), fontsize=12)
    ax.set_ylim(-0.15, 3.0); ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=13.5, fontweight="bold", loc="left")
    ax.margins(x=0.08)

draw(axL, "e", "A   Inside-V enrichment cut-off", "log2(V-in / V-out)", "")
axL.set_xlim(-0.3, None)
draw(axR, "w", "B   V-channel width cut-off", "V-channel width (bp)", " bp")
axR.set_xlim(-3, None)

# legend (shared) + rule
h_tf = plt.Line2D([], [], marker="o", color=C_TF, ls="none", ms=8, label="TF motif")
h_bi = plt.Line2D([], [], marker="o", mfc="none", mec=C_BIAS, mew=1.7, ls="none", ms=8,
                  label="enzyme bias")
axL.legend(handles=[h_tf, h_bi], fontsize=11, loc="lower right", frameon=False)

fig.text(0.5, 0.025,
         "cut-off = midpoint of the gap between the highest bias and the lowest TF   ·   "
         "margin = half the gap width   (every assay perfectly separated → LOO 100%)",
         ha="center", fontsize=11.5, color="#333333")

fig.suptitle("How each cut-off is derived (mid-gap rule)", fontsize=15, fontweight="bold", y=0.99)
fig.subplots_adjust(left=0.075, right=0.985, top=0.9, bottom=0.13, wspace=0.18)

for assay in assays:
    maxb, mintf, mid, margin = midgap(assay, "w")
    threshold_rows.append({
        "assay": assay, "feature": "V-channel width", "cutoff": mid,
        "margin": margin, "max_bias": maxb, "min_TF": mintf, "unit": "bp"
    })
    maxb, mintf, mid, margin = midgap(assay, "e")
    threshold_rows.append({
        "assay": assay, "feature": "Inside-V enrichment", "cutoff": mid,
        "margin": margin, "max_bias": maxb, "min_TF": mintf, "unit": "log2"
    })

with OUTCSV.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(threshold_rows[0]))
    writer.writeheader()
    writer.writerows(threshold_rows)

fig.savefig(OUTBASE.with_suffix(".png"), dpi=200, facecolor="white")
fig.savefig(OUTBASE.with_suffix(".pdf"), facecolor="white")
print("saved:", f"{OUTBASE}.png/.pdf")
print("saved:", OUTCSV)
