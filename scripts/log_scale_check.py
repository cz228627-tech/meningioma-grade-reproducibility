# -*- coding: utf-8 -*-
"""Robustness: GSE16581 linear vs log2 scale"""
import sys, os, re, gzip
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
SEED, NPERM = 20261103, 1000
rng = np.random.default_rng(SEED)
W = np.array([-2, -1, 0, 1, 2], float)


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
ed = eff(np.log2(ddf[cols].values + 1), gd)
didx = {g: i for i, g in enumerate(ddf.index)}

p2g = {}
with gzip.open(os.path.join(OUT, "GPL570.annot.gz"), "rt", encoding="utf-8", errors="replace") as f:
    cc = None
    for line in f:
        if line.startswith("!platform_table_begin"):
            cc = f.readline().rstrip("\n").split("\t")
            break
    i_id, i_sym = cc.index("ID"), cc.index("Gene symbol")
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        p = line.rstrip("\n").split("\t")
        if len(p) > i_sym:
            s = p[i_sym].split("///")[0].strip()
            if s and s != "---":
                p2g[p[i_id]] = s

chars, rows = [], []
with gzip.open(os.path.join(OUT, "GSE16581_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
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
gv = np.array([int(x) for x in L["who grade"]])
V = np.array([[float(x) if x.strip('"') not in ("", "NA", "null") else np.nan for x in r[1:]]
              for r in rows[1:]], dtype=float)
probes = [r[0].strip('"') for r in rows[1:]]
keep = {}
for k, pr in enumerate(probes):
    g = p2g.get(pr)
    if not g:
        continue
    if g not in keep or np.nanmean(V[k]) > np.nanmean(V[keep[g]]):
        keep[g] = k
genes = list(keep)
ridx = np.array([keep[g] for g in genes])
common = [(didx[g], j) for j, g in enumerate(genes) if g in didx]
di = np.array([c[0] for c in common])
vi = np.array([c[1] for c in common])
ed_s = ed[di]

for label, X in [("linear (as used)", V), ("log2(V+1)", np.log2(V + 1))]:
    Xv = X[ridx]
    ev = eff(Xv, gv)[vi]
    fin = np.isfinite(ev) & np.isfinite(ed_s) & (ed_s != 0)
    absd = np.abs(ed_s[fin])
    edges = np.quantile(absd, [0, .2, .4, .6, .8, 1.0])
    ql = np.digitize(absd, edges[1:-1], right=False)
    ia = np.where(fin)[0]
    qi = [ia[ql == q] for q in range(5)]
    obs = np.array([np.mean(np.sign(ev[ix]) == np.sign(ed_s[ix])) for ix in qi])
    perm = np.zeros(NPERM)
    for b in range(NPERM):
        e = eff(Xv, rng.permutation(gv))[vi]
        perm[b] = W @ np.array([np.mean(np.sign(e[ix]) == np.sign(ed_s[ix])) for ix in qi])
    stat = float(W @ obs)
    p = (np.sum(perm >= stat) + 1) / (NPERM + 1)
    print(f"{label:18s} quintiles {[round(100*x,2) for x in obs]} | trend {stat:+.4f} | p={p:.4f}")
