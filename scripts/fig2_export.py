# -*- coding: utf-8 -*-
"""导出 Fig 2 的效应量数据（发现队列，GSE16581 匹配集）"""
import sys, os, re, gzip
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))


def eff(X, g):
    return X[:, g >= 2].mean(1) - X[:, g == 1].mean(1)


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
ed = eff(Xd, gd)
didx = {g: i for i, g in enumerate(ddf.index)}


def platform(annot):
    p2g = {}
    with gzip.open(os.path.join(OUT, annot), "rt", encoding="utf-8", errors="replace") as f:
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


p2g = platform("GPL570.annot.gz")
L, V, pr = mat("GSE16581_matrix.txt.gz")
keep = set()
for k, x in enumerate(pr):
    g = p2g.get(x)
    if g:
        keep.add(g)

genes = [g for g in keep if g in didx]
vals = np.array([ed[didx[g]] for g in genes])
ok = np.isfinite(vals) & (vals != 0)
genes = np.array(genes)[ok]; vals = vals[ok]
edges = np.quantile(np.abs(vals), [0, .2, .4, .6, .8, 1.0])
q = np.digitize(np.abs(vals), edges[1:-1], right=False) + 1

df = pd.DataFrame({"gene": genes, "effect": np.round(vals, 6),
                   "abs_effect": np.round(np.abs(vals), 6), "quintile": q})
p = os.path.join(OUT, "fig2_effects.csv")
df.to_csv(p, index=False, encoding="utf-8-sig")
print("n =", len(df), "| boundaries =", [round(float(x), 4) for x in edges])
print(df.groupby("quintile").size().to_string())
print("saved:", p)
