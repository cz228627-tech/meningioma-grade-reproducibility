# -*- coding: utf-8 -*-
"""Fig 4 ç»ˆç‰ˆï¼ˆlog2 ä¸»è§„æ ¼ï¼Œä¿®æ­£æ ‡é¢˜/è„šæ³¨/é¢æ¿å­—æ¯ï¼‰"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\work\topic-audit\data"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "mathtext.default": "regular"})

df4 = pd.read_csv(os.path.join(OUT, "fig4_strata_v2.csv"))
colors = {"GSE16581": "#4C72B0", "GSE74385": "#DD8452"}
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)

for ax, cohort in zip(axes, ["GSE16581", "GSE74385"]):
    d = df4[df4.cohort == cohort].reset_index(drop=True)
    x = np.arange(len(d))
    delta = d.delta_pct.values
    shades = ["0.62"] + [colors[cohort]] * 4
    alphas = [1.0, 0.45, 0.65, 0.8, 1.0]
    for i in range(len(d)):
        ax.bar(x[i], delta[i], width=0.62, color=shades[i], alpha=alphas[i],
               edgecolor="black", linewidth=0.6, zorder=2)
        p = d.trend_p.iloc[i]
        star = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
        ax.text(x[i], delta[i] + 0.8, star, ha="center", va="bottom", fontsize=9)
        ax.text(x[i], delta[i] / 2, f"{delta[i]:.1f}", ha="center", va="center", fontsize=7,
                color="white" if i else "black")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"E{i+1}" for i in range(len(d))])
    ax.set_xlabel("Mean expression stratum (E1 = lowest)")
    ax.set_title(cohort, fontsize=9.5)
    ax.set_ylim(0, 42)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.03, 0.04, "gradient = concordance(Q5) \u2212 concordance(Q1)",
            transform=ax.transAxes, fontsize=6.6, color="0.35")

axes[0].set_ylabel("Concordance difference (percentage points)")
fig.suptitle("Effect-size gradients are attenuated in the lowest-expression stratum", fontsize=10, y=1.02)
fig.text(0.5, -0.08,
         "Markers indicate permutation-calibrated quintile-trend P (** P < 0.01, * P < 0.05, ns); they do not test "
         "the difference between Q5 and Q1. Unadjusted and exploratory.",
         ha="center", fontsize=6.4, color="0.35")
for ax, L in zip(axes, "ab"):
    ax.text(0.01, 1.02, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
fig.tight_layout()
for ext in ("png", "pdf"):
    p = os.path.join(OUT, f"Fig4.{ext}")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print("saved", p)

