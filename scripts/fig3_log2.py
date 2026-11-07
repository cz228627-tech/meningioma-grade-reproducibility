# -*- coding: utf-8 -*-
"""Fig 3 å¢žåŠ é¢æ¿ cï¼šä¸¤éªŒè¯é˜Ÿåˆ—åŒæ—¶åŒå‘"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

OUT = r"C:\work\topic-audit\data"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "mathtext.default": "regular"})

df3 = pd.read_csv(os.path.join(OUT, "fig3_final_v4_logmain.csv"))
bc = pd.read_csv(os.path.join(OUT, "both_cohorts_v2.csv"))
colors = {"GSE16581": "#4C72B0", "GSE74385": "#DD8452"}

fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.1))

# a, b: å·²æœ‰ä¸¤é˜Ÿåˆ—
for ax, cohort in zip(axes[:2], ["GSE16581", "GSE74385"]):
    d = df3[df3.cohort == cohort].reset_index(drop=True)
    x = np.arange(len(d))
    obs, q = d.obs_pct.values, d.p_BH.values
    lo, hi = d.ci_lo_pct.values, d.ci_hi_pct.values
    nm, nsd = d.null_mean_pct.values, d.null_sd_pct.values
    for i in range(len(d)):
        ax.fill_between([x[i] - 0.42, x[i] + 0.42], nm[i] - 1.96 * nsd[i], nm[i] + 1.96 * nsd[i],
                        color="0.87", zorder=0, linewidth=0)
    ax.plot([-0.6, len(d) - 0.4], [nm.mean()] * 2, "--", color="0.45", lw=1.0, zorder=1)
    ax.bar(x, obs, width=0.62, color=colors[cohort], alpha=0.85, edgecolor="black", linewidth=0.6, zorder=2)
    ax.errorbar(x, obs, yerr=[np.maximum(obs - lo, 0), np.maximum(hi - obs, 0)], fmt="none", ecolor="black", elinewidth=0.8, capsize=3, zorder=3)
    for i in range(len(d)):
        star = "**" if q[i] < 0.01 else ("*" if q[i] < 0.05 else "ns")
        ax.text(x[i], 93, star, ha="center", va="center", fontsize=9.5)
        ax.text(x[i], hi[i] + 1.5, f"{obs[i]:.1f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([f"Q{i+1}" for i in range(len(d))])
    nmin, nmax = int(d.n.min()), int(d.n.max())
    nlab = f"{nmin:,}" if nmin == nmax else f"{nmin:,}-{nmax:,}"
    ax.set_title(f"{cohort}  (n = {nlab}; trend P = {d.trend_p.iloc[0]:g})", fontsize=9)
    ax.set_xlabel("Discovery effect-size quintile (Q1 = weakest)")
    ax.set_ylim(0, 100)
    ax.spines[["top", "right"]].set_visible(False)

# c: åŒæ—¶åŒå‘
ax = axes[2]
x = np.arange(len(bc))
obs, q = bc.both_pct.values, bc.q.values
lo, hi = bc.ci_lo.values, bc.ci_hi.values
nm, nsd = bc.null_mean.values, bc.null_sd.values
for i in range(len(bc)):
    ax.fill_between([x[i] - 0.42, x[i] + 0.42], nm[i] - 1.96 * nsd[i], nm[i] + 1.96 * nsd[i],
                    color="0.87", zorder=0, linewidth=0)
ax.plot([-0.6, len(bc) - 0.4], [nm.mean()] * 2, "--", color="0.45", lw=1.0, zorder=1)
ax.bar(x, obs, width=0.62, color="#5F9E6E", alpha=0.9, edgecolor="black", linewidth=0.6, zorder=2)
ax.errorbar(x, obs, yerr=[np.maximum(obs - lo, 0), np.maximum(hi - obs, 0)], fmt="none", ecolor="black", elinewidth=0.8, capsize=3, zorder=3)
for i in range(len(bc)):
    star = "**" if q[i] < 0.01 else ("*" if q[i] < 0.05 else "ns")
    ax.text(x[i], 93, star, ha="center", va="center", fontsize=9.5)
    ax.text(x[i], hi[i] + 1.5, f"{obs[i]:.1f}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels([f"Q{i+1}" for i in range(len(bc))])
n = int(bc.n.iloc[0])
ax.set_title(f"Both cohorts  (n = {n:,}; trend P = {bc.trend_p.iloc[0]:g})", fontsize=9)
ax.set_xlabel("Discovery effect-size quintile (Q1 = weakest)")
ax.set_ylim(0, 100)
ax.spines[["top", "right"]].set_visible(False)

axes[0].set_ylabel("Directional concordance (%)", fontsize=8.5)
fig.legend(handles=[Patch(facecolor="#4C72B0", alpha=0.85, edgecolor="black", label="GSE16581"),
                    Patch(facecolor="#DD8452", alpha=0.85, edgecolor="black", label="GSE74385"),
                    Patch(facecolor="#5F9E6E", alpha=0.9, edgecolor="black", label="Both cohorts"),
                    Line2D([0], [0], color="0.45", ls="--", lw=1.0, label="Permutation null mean"),
                    Patch(facecolor="0.87", label=r"Null mean $\pm$ 1.96 SD")],
           loc="lower center", ncol=5, frameon=False, fontsize=7.6, bbox_to_anchor=(0.5, -0.10))
for ax, L in zip(axes, "abc"):
    ax.text(0.01, 1.16, L, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
fig.tight_layout()
for ext in ("png", "pdf"):
    p = os.path.join(OUT, f"Fig3.{ext}")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print("saved", p)



