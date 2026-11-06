# -*- coding: utf-8 -*-
"""
G2-only 复现性分析（WHO 2 vs 1）与原口径（WHO 2+3 vs 1）对照。
"""
import sys, os, re, gzip, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", r"C:\work\topic-audit\data")
SEED, NPERM, NBOOT = 20261103, 1000, 1000
rng = np.random.default_rng(SEED)
W = np.array([-2, -1, 0, 1, 2], float)

GMT = os.path.join(OUT, "h.all.v2023.2.Hs.symbols.gmt")
sets = {}
for line in open(GMT, encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    sets[f[0]] = set(x for x in f[2:] if x)
PROG = sets["HALLMARK_E2F_TARGETS"] | sets["HALLMARK_G2M_CHECKPOINT"]

# ---------- 载入 ----------
ddf = pd.read_csv(os.path.join(OUT, "GSE183653_meningioma_tpm.csv.gz"), index_col=0)
ch, t2 = [], []
with gzip.open(os.path.join(OUT, "RNA_GSE183653.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            t2 = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            ch.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            break
lab = {c[0].split(":", 1)[0].strip(): [x.split(":", 1)[-1].strip() for x in c] for c in ch}
m2g = {}
for i, t in enumerate(t2):
    m = re.search(r"Meningioma\s+(\d+)\s*\(", t)
    if m:
        m2g["M" + m.group(1)] = int(lab["grade"][i])
cols = [c for c in ddf.columns if c in m2g]
gd = np.array([m2g[c] for c in cols])
Xd = np.log2(ddf[cols].values + 1)
didx = {g: i for i, g in enumerate(ddf.index)}
print("发现队列:", dict(zip(*np.unique(gd, return_counts=True))))


def platform(annot):
    p2g = {}
    with gzip.open(os.path.join(OUT, annot), "rt", encoding="utf-8", errors="replace") as f:
        cc = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                cc = f.readline().rstrip("\n").split("\t"); break
        i_id, i_sym = cc.index("ID"), cc.index("Gene symbol")
        for line in f:
            if line.startswith("!platform_table_end"):
                break
            p = line.rstrip("\n").split("\t")
            if len(p) > i_sym:
                s = p[i_sym].split("///")[0].strip()
                if s and s != "---":
                    p2g[p[i_id]] = s
    return p2g


def mat(name):
    chars, rows = [], []
    with gzip.open(os.path.join(OUT, name), "rt", encoding="utf-8", errors="replace") as f:
        on = False
        for line in f:
            if line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                on = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif on:
                rows.append(line.rstrip("\n").split("\t"))
    L = {c[0].split(":", 1)[0].strip(): [x.split(":", 1)[-1].strip() for x in c] for c in chars}
    V = np.array([[float(x) if x.strip('"') not in ("", "NA", "null") else np.nan for x in r[1:]]
                  for r in rows[1:]], dtype=float)
    return L, V, [r[0].strip('"') for r in rows[1:]]


def prep(mtx, annot, grader):
    p2g = platform(annot); L, V, pr = mat(mtx); gv = grader(L)
    keep = {}
    for k, x in enumerate(pr):
        g = p2g.get(x)
        if not g:
            continue
        if g not in keep or np.nanmean(V[k]) > np.nanmean(V[keep[g]]):
            keep[g] = k
    return V, keep, gv, L


V1, keep1, gv1, L1 = prep("GSE16581_matrix.txt.gz", "GPL570.annot.gz",
                          lambda L: np.array([int(x) for x in L["who grade"]]))
V2r, keep2, gva, L2 = prep("GSE74385_matrix.txt.gz", "GPL10558.annot.gz",
                           lambda L: np.array([int(x[0]) for x in L["subtype"]]))
km = np.array([not x[1:].startswith("M") for x in L2["subtype"]])
V2 = V2r[:, km]; gv2 = gva[km]
print("GSE16581:", dict(zip(*np.unique(gv1, return_counts=True))))
print("GSE74385:", dict(zip(*np.unique(gv2, return_counts=True))))

if "batch" in L2:
    b = np.array(L2["batch"])[km]
    print("\n=== GSE74385 批次 × 分级 ===")
    print(pd.crosstab(pd.Series(b, name="batch"), pd.Series(gv2, name="grade")).to_string())

results = []


def analyse(name, D, V, keep, gv, hi_fn, lo_fn, tag):
    genes = [g for g in keep if g in didx]
    ridx = np.array([keep[g] for g in genes])
    sel = np.where(hi_fn(gv) | lo_fn(gv))[0]
    sub = V[ridx][:, sel]
    lab_ = gv[sel]
    him, lom = hi_fn(lab_), lo_fn(lab_)
    e = sub[:, him].mean(1) - sub[:, lom].mean(1)
    ed = D[[didx[g] for g in genes]]
    fin = np.isfinite(e) & np.isfinite(ed) & (ed != 0)
    idx = np.where(fin)[0]
    edges = np.quantile(np.abs(ed[idx]), [0, .2, .4, .6, .8, 1.0])
    ql = np.digitize(np.abs(ed[idx]), edges[1:-1], right=False)
    groups = [idx[ql == q] for q in range(5)]
    obs = np.array([np.mean(np.sign(e[g]) == np.sign(ed[g])) for g in groups])
    perm = np.empty(NPERM)
    for b in range(NPERM):
        gp = rng.permutation(lab_)
        eb = sub[:, gp == 2].mean(1) - sub[:, gp == 1].mean(1)
        perm[b] = W @ np.array([np.mean(np.sign(eb[g]) == np.sign(ed[g])) for g in groups])
    stat = float(W @ obs)
    p = (np.sum(perm >= stat) + 1) / (NPERM + 1)
    boot = []
    for b in range(NBOOT):
        ii = rng.integers(0, len(lab_), len(lab_))
        if (lab_[ii] == 2).sum() < 3 or (lab_[ii] == 1).sum() < 3:
            continue
        eb = sub[:, ii][:, lab_[ii] == 2].mean(1) - sub[:, ii][:, lab_[ii] == 1].mean(1)
        boot.append([np.mean(np.sign(eb[g]) == np.sign(ed[g])) for g in groups])
    boot = np.array(boot) * 100
    print(f"\n--- {tag} | {name} (n_hi={int(him.sum())}, n_lo={int(lom.sum())}) ---")
    print("  一致率: " + "  ".join(f"{100*x:5.1f}" for x in obs))
    print("  95%CI : " + "  ".join(f"{np.percentile(boot[:,i],2.5):.0f}-{np.percentile(boot[:,i],97.5):.0f}" for i in range(5)))
    print(f"  趋势 {stat:+.4f} vs 零 {perm.mean():+.4f}±{perm.std():.4f} | P={p:.4f}")
    results.append(dict(cohort=name, tag=tag, n=int(fin.sum()),
                        n_hi=int(him.sum()), n_lo=int(lom.sum()),
                        obs=[round(100*x, 2) for x in obs],
                        ci=[[round(np.percentile(boot[:, i], 2.5), 1), round(np.percentile(boot[:, i], 97.5), 1)] for i in range(5)],
                        trend=round(stat, 4), null=round(float(perm.mean()), 4), sd=round(float(perm.std()), 4),
                        p=round(float(p), 4)))


d_g2 = Xd[:, gd == 2].mean(1) - Xd[:, gd == 1].mean(1)
d_mix = Xd[:, gd >= 2].mean(1) - Xd[:, gd == 1].mean(1)

analyse("GSE16581", d_g2, V1, keep1, gv1, lambda g: g == 2, lambda g: g == 1, "G2-only")
analyse("GSE74385", d_g2, V2, keep2, gv2, lambda g: g == 2, lambda g: g == 1, "G2-only")
analyse("GSE16581", d_mix, V1, keep1, gv1, lambda g: g >= 2, lambda g: g == 1, "原口径 G2+3")
analyse("GSE74385", d_mix, V2, keep2, gv2, lambda g: g >= 2, lambda g: g == 1, "原口径 G2+3")

json.dump(results, open(os.path.join(OUT, "g2only_results.json"), "w"), indent=1)
print("\nsaved g2only_results.json")
