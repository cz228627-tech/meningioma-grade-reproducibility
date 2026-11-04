# -*- coding: utf-8 -*-
"""Regenerate all four figures with corrections (reads frozen CSVs only)"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.lines import Line2D

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "mathtext.default": "regular"})


def panel_letters(fig, letters, y=1.02):
    for ax, L in zip(fig.axes, letters):
        ax.text(0.01, y, L, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    for ext in ("png", "pdf"):
        p = os.path.join(OUT, f"{name}.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print("saved:", p)


# ================= FIGURE 1 =================
comp = pd.read_csv(os.path.join(OUT, "fig1_cohorts.csv"))
qc = pd.read_csv(os.path.join(OUT, "fig1_qc_rna.csv"))

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.12], hspace=0.55, wspace=0.40)
axA = fig.add_subplot(gs[0, :]); axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")


def box(x, y, w, h, title, lines, fc):
    axA.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.010,rounding_size=0.02",
                                 fc=fc, ec="0.35", lw=0.8))
    axA.text(x + w / 2, y + h - 0.05, title, ha="center", va="top", fontsize=8.2, weight="bold")
    axA.text(x + w / 2, y + h - 0.15, "\n".join(lines), ha="center", va="top",
             fontsize=6.7, linespacing=1.45)


box(0.02, 0.56, 0.38, 0.40, "Discovery",
    ["GSE183653 (RNA-seq)", "n = 185 tumours", "WHO 1/2/3 = 86/74/25"], "#DCE6F1")
box(0.60, 0.55, 0.38, 0.40, "Validation 1",
    ["GSE16581", "Affymetrix U133 Plus 2.0", "n = 68; WHO 43/19/6"], "#FBE3D0")
box(0.60, 0.03, 0.38, 0.40, "Validation 2",
    ["GSE74385", "Illumina HumanHT-12 v4", "n = 53; WHO 17/8/28"], "#FBE3D0")
axA.add_patch(FancyArrowPatch((0.405, 0.79), (0.595, 0.79), arrowstyle="-|>", mutation_scale=9,
                              lw=0.9, color="0.35", connectionstyle="arc3,rad=-0.12"))
axA.add_patch(FancyArrowPatch((0.405, 0.66), (0.595, 0.23), arrowstyle="-|>", mutation_scale=9,
                              lw=0.9, color="0.35", connectionstyle="arc3,rad=0.15"))
axA.text(0.21, 0.47, "18,982 / 18,185 genes after effect filtering;\n5 effect-size quintiles",
         ha="center", va="top", fontsize=6.6, color="0.35")

axB = fig.add_subplot(gs[1, 0])
order = ["GSE183653", "GSE16581", "GSE74385"]
piv = comp.pivot(index="cohort", columns="grade", values="n").loc[order]
pct = piv.div(piv.sum(axis=1), axis=0) * 100
cols = {1: "#9EC3E6", 2: "#5B8FC0", 3: "#1F4E79"}
bottom = np.zeros(len(order))
for g in (1, 2, 3):
    axB.bar(np.arange(len(order)), pct[g].values, bottom=bottom, width=0.62,
            color=cols[g], edgecolor="white", linewidth=0.8, label=f"WHO {g}")
    for i, v in enumerate(pct[g].values):
        if v > 8:
            axB.text(i, bottom[i] + v / 2, f"{v:.0f}%", ha="center", va="center", fontsize=7,
                     color="white" if g > 1 else "black")
    bottom += pct[g].values
axB.set_xticks(np.arange(len(order)))
axB.set_xticklabels([f"{o}\n(n={int(piv.loc[o].sum())})" for o in order], fontsize=7)
axB.set_ylabel("Proportion of tumours (%)", fontsize=8)
axB.set_ylim(0, 100); axB.set_title("Grade composition", fontsize=9)
axB.legend(frameon=False, fontsize=6.8, ncol=3, loc="upper center",
           bbox_to_anchor=(0.5, -0.26), handlelength=1.0, columnspacing=1.2)
axB.spines[["top", "right"]].set_visible(False)

axC = fig.add_subplot(gs[1, 1])
low = qc.flag == "low"
axC.scatter(qc.loc[~low, "detected_genes"], qc.loc[~low, "median_corr"], s=11,
            c="#8FA8CC", edgecolor="none", alpha=0.85, label="QC not flagged")
axC.scatter(qc.loc[low, "detected_genes"], qc.loc[low, "median_corr"], s=30,
            c="#C0392B", marker="X", edgecolor="black", linewidth=0.4, label="QC flagged (retained)")
for _, r in qc[low].iterrows():
    axC.annotate(r["sample"], (r.detected_genes, r["median_corr"]), textcoords="offset points",
                 xytext=(5, -2), fontsize=6, color="#C0392B")
axC.set_xlabel("Detected genes (TPM > 1)", fontsize=8)
axC.set_ylabel("Median correlation with other samples", fontsize=8)
axC.set_title("Sample-level QC (discovery cohort)", fontsize=9)
lg = axC.legend(frameon=False, fontsize=6.6, loc="upper left",
                title="flagged: <12,000 genes\nor r < 0.86")
lg.get_title().set_fontsize(5.8)
axC.spines[["top", "right"]].set_visible(False)
panel_letters(fig, "abc")
save(fig, "Fig1")

# ================= FIGURE 2 =================
df2 = pd.read_csv(os.path.join(OUT, "fig2_effects.csv"))
e, a = df2.effect.values, df2.abs_effect.values
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))

ax = axes[0]
ax.hist(e, bins=140, range=(e.min(), e.max()), color="#8C8C8C", edgecolor="none")
ax.axvline(0, color="black", lw=0.9, ls="--")
ax.set_xlabel("Discovery effect (log$_2$(TPM+1) mean difference)")
ax.set_ylabel("Number of genes")
ax.set_title("Signed effects (full range)", fontsize=9.5)
ax.spines[["top", "right"]].set_visible(False)
ax.text(0.03, 0.95, f"n = {len(e):,}\nrange {e.min():.2f} to {e.max():.2f}",
        transform=ax.transAxes, fontsize=7, va="top")

ax = axes[1]
edges = [0.0, 0.0305, 0.0894, 0.1740, 0.3057, a.max()]
cols5 = ["#DCE6F1", "#AFC4DE", "#7E9CC4", "#4E75A6", "#2B5183"]
bins = np.logspace(np.log10(max(a.min(), 1e-4)), np.log10(a.max()), 70)
n, be = np.histogram(a, bins=bins)
cent = np.sqrt(be[:-1] * be[1:])
for i in range(len(n)):
    q = min(max(int(np.searchsorted(edges, cent[i], side="right")) - 1, 0), 4)
    ax.bar(be[i], n[i], width=be[i + 1] - be[i], align="edge", color=cols5[q], edgecolor="none")
for b in edges[1:-1]:
    ax.axvline(b, color="0.35", lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_xlabel("|Discovery effect| (log scale)")
ax.set_ylabel("Number of genes")
ax.set_title("Absolute effects by quintile", fontsize=9.5)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[Patch(facecolor=cols5[i], label=f"Q{i+1}  {edges[i]:.3f}-{edges[i+1]:.3f}")
                   for i in range(5)],
          loc="upper left", frameon=False, fontsize=6.4, handlelength=1.0,
          borderpad=0.2, labelspacing=0.25)
panel_letters(fig, "ab")
fig.tight_layout()
save(fig, "Fig2")

# ================= FIGURE 3 =================
df3 = pd.read_csv(os.path.join(OUT, "fig3_final_v3.csv"))
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
colors = {"GSE16581": "#4C72B0", "GSE74385": "#DD8452"}
for ax, cohort in zip(axes, ["GSE16581", "GSE74385"]):
    d = df3[df3.cohort == cohort].reset_index(drop=True)
    x = np.arange(len(d))
    obs, q = d.obs_pct.values, d.p_BH.values
    lo, hi = d.ci_lo_pct.values, d.ci_hi_pct.values
    nm, nsd = d.null_mean_pct.values, d.null_sd_pct.values
    for i in range(len(d)):
        ax.fill_between([x[i] - 0.42, x[i] + 0.42], nm[i] - 1.96 * nsd[i], nm[i] + 1.96 * nsd[i],
                        color="0.87", zorder=0, linewidth=0)
    ax.plot([-0.6, len(d) - 0.4], [nm.mean()] * 2, "--", color="0.45", lw=1.0, zorder=1)
    ax.bar(x, obs, width=0.62, color=colors[cohort], alpha=0.85, edgecolor="black",
           linewidth=0.6, zorder=2)
    ax.errorbar(x, obs, yerr=[obs - lo, hi - obs], fmt="none", ecolor="black",
                elinewidth=0.8, capsize=3, zorder=3)
    for i in range(len(d)):
        star = "**" if q[i] < 0.01 else ("*" if q[i] < 0.05 else "ns")
        ax.text(x[i], 93, star, ha="center", va="center", fontsize=9.5)
        ax.text(x[i], hi[i] + 1.5, f"{obs[i]:.1f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([f"Q{i+1}" for i in range(len(d))])
    nmin, nmax = int(d.n.min()), int(d.n.max())
    nlab = f"{nmin:,}" if nmin == nmax else f"{nmin:,}-{nmax:,}"
    ax.set_title(f"{cohort}  (n = {nlab}; trend P = {d.trend_p.iloc[0]:g})", fontsize=9)
    ax.set_xlabel("Discovery effect-size quintile (Q1 = weakest)")
    ax.set_ylim(0, 100); ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Directional concordance (%)")
fig.legend(handles=[Patch(facecolor="#4C72B0", alpha=0.85, edgecolor="black", label="Observed concordance"),
                    Line2D([0], [0], color="0.45", ls="--", lw=1.0, label="Permutation null mean"),
                    Patch(facecolor="0.87", label=r"Null mean $\pm$ 1.96 SD")],
           loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.09))
fig.tight_layout()
panel_letters(fig, "ab", y=1.14)
save(fig, "Fig3")

# ================= FIGURE 4 =================
df4 = pd.read_csv(os.path.join(OUT, "fig4_strata.csv"))
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, cohort in zip(axes, ["GSE16581", "GSE74385"]):
    d = df4[df4.cohort == cohort].reset_index(drop=True)
    x = np.arange(len(d)); delta = d.delta_pct.values
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
    ax.set_ylim(0, 42); ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.03, 0.90, "gradient = concordance(Q5) - concordance(Q1)",
            transform=ax.transAxes, fontsize=6.6, color="0.35")
axes[0].set_ylabel("Concordance difference (percentage points)")
fig.suptitle("Effect-size gradients are attenuated in the lowest-expression stratum",
             fontsize=10, y=1.02)
fig.text(0.5, -0.08,
         "Markers indicate permutation-calibrated quintile-trend P (** P < 0.01, * P < 0.05, ns); "
         "they do not test the difference between Q5 and Q1. Unadjusted and exploratory.",
         ha="center", fontsize=6.4, color="0.35")
panel_letters(fig, "ab")
fig.tight_layout()
save(fig, "Fig4")

