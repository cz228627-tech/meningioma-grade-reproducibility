# -*- coding: utf-8 -*-
"""
Gate F1：67.69% 是否只是两个边际一致率的乘积？
在同一共同基因集、同一分母上计算：
  p1 = 与 GSE16581 同向率
  p2 = 与 GSE74385 同向率
  joint = 两队列同时同向率
  product = p1 * p2
检验：joint 是否超出 product（配对置换，保留边际）
"""
import sys, os, re, gzip, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

OUT = os.environ.get("DATA_DIR", r"C:\work\topic-audit\data")
rng = np.random.default_rng(20261103)
NPERM = 10000

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
    if m: m2g["M" + m.group(1)] = int(lab["grade"][i])
cols = [c for c in ddf.columns if c in m2g]
gd = np.array([m2g[c] for c in cols])
Xd = np.log2(ddf[cols].values + 1)
d_all = Xd[:, gd >= 2].mean(1) - Xd[:, gd == 1].mean(1)
didx = {g: i for i, g in enumerate(ddf.index)}


def platform(annot):
    p2g = {}
    with gzip.open(os.path.join(OUT, annot), "rt", encoding="utf-8", errors="replace") as f:
        cc = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                cc = f.readline().rstrip("\n").split("\t"); break
        i_id, i_sym = cc.index("ID"), cc.index("Gene symbol")
        for line in f:
            if line.startswith("!platform_table_end"): break
            p = line.rstrip("\n").split("\t")
            if len(p) > i_sym:
                s = p[i_sym].split("///")[0].strip()
                if s and s != "---": p2g[p[i_id]] = s
    return p2g


def mat(name):
    chars, rows = [], []
    with gzip.open(os.path.join(OUT, name), "rt", encoding="utf-8", errors="replace") as f:
        on = False
        for line in f:
            if line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"): on = True
            elif line.startswith("!series_matrix_table_end"): break
            elif on: rows.append(line.rstrip("\n").split("\t"))
    L = {}
    for c in chars:
        k = c[0].split(":", 1)[0].strip()
        if k not in L: L[k] = [x.split(":", 1)[-1].strip() for x in c]
    V = np.array([[float(x) if x.strip('"') not in ("", "NA", "null") else np.nan for x in r[1:]]
                  for r in rows[1:]], dtype=float)
    return L, V, [r[0].strip('"') for r in rows[1:]]


def prep(mtx, annot, grader):
    p2g = platform(annot); L, V, pr = mat(mtx); gv = grader(L)
    keep = {}
    for k, x in enumerate(pr):
        g = p2g.get(x)
        if not g: continue
        if g not in keep or np.nanmean(V[k]) > np.nanmean(V[keep[g]]): keep[g] = k
    return V, keep, gv, L


V1, keep1, gv1, _ = prep("GSE16581_matrix.txt.gz", "GPL570.annot.gz",
                         lambda L: np.array([int(x) for x in L["who grade"]]))
V2r, keep2, gva, L2 = prep("GSE74385_matrix.txt.gz", "GPL10558.annot.gz",
                           lambda L: np.array([int(x[0]) for x in L["subtype"]]))
km = np.array([not s[1:].startswith("M") for s in L2["subtype"]])
V2 = V2r[:, km]; gv2 = gva[km]

common = [g for g in keep1 if g in keep2 and g in didx
          and np.isfinite(d_all[didx[g]]) and d_all[didx[g]] != 0]
A1 = np.log2(V1[[keep1[g] for g in common]] + 1); A2 = V2[[keep2[g] for g in common]]
dd = d_all[[didx[g] for g in common]]
e1 = A1[:, gv1 >= 2].mean(1) - A1[:, gv1 == 1].mean(1)
e2 = A2[:, gv2 >= 2].mean(1) - A2[:, gv2 == 1].mean(1)
ok = np.isfinite(e1) & np.isfinite(e2)
idx = np.where(ok)[0]
edges = np.quantile(np.abs(dd[idx]), [0, .2, .4, .6, .8, 1.0])
ql = np.digitize(np.abs(dd[idx]), edges[1:-1], right=False)
s1 = np.sign(dd) == np.sign(e1)
s2 = np.sign(dd) == np.sign(e2)
print("共同基因:", len(common), "| 有效:", int(ok.sum()))

rows = []
for q in range(5):
    sel = idx[ql == q]
    a, b = s1[sel], s2[sel]
    p1, p2 = a.mean(), b.mean()
    joint = (a & b).mean()
    prod = p1 * p2
    # 配对置换：固定边际，打乱 b 与 a 的配对
    null = np.empty(NPERM)
    for i in range(NPERM):
        null[i] = (a & rng.permutation(b)).mean()
    p_perm = (np.sum(null >= joint) + 1) / (NPERM + 1)
    tab = [[int((a & b).sum()), int((a & ~b).sum())],
           [int((~a & b).sum()), int((~a & ~b).sum())]]
    orr, p_fisher = fisher_exact(tab)
    rows.append(dict(quintile=f"Q{q+1}", n=len(sel),
                     p1_16581=round(100 * p1, 2), p2_74385=round(100 * p2, 2),
                     joint=round(100 * joint, 2), product=round(100 * prod, 2),
                     excess=round(100 * (joint - prod), 2),
                     null_mean=round(100 * null.mean(), 2), null_sd=round(100 * null.std(), 2),
                     p_perm=round(float(p_perm), 4), OR=round(float(orr), 3), p_fisher=float(f"{p_fisher:.3g}")))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "joint_vs_marginals_v2.csv"), index=False, encoding="utf-8-sig")
print("\n" + df.to_string(index=False))
