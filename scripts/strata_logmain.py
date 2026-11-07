# -*- coding: utf-8 -*-
"""表达量分层分析（含零效应剔除），并把完整 5x5 矩阵存成 CSV 供作图"""
import sys, os, re, gzip
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = r"C:\work\topic-audit\data"
SEED, NPERM = 20261103, 500
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
Xd = np.log2(ddf[cols].values + 1)
ed = eff(Xd, gd)
me = Xd.mean(1)
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


def prep(mtx, annot, grader):
    p2g = platform(annot); L, V, pr = mat(mtx); gv = grader(L)
    keep = {}
    for k, x in enumerate(pr):
        g = p2g.get(x)
        if not g:
            continue
        if g not in keep or np.nanmean(V[k]) > np.nanmean(V[keep[g]]):
            keep[g] = k
    return V, keep, gv


V1r, keep1, gv1 = prep("GSE16581_matrix.txt.gz", "GPL570.annot.gz",
                      lambda L: np.array([int(x) for x in L["who grade"]]))
V1 = np.log2(V1r + 1)   # 主规格：log2 变换
V2r, keep2, gva = prep("GSE74385_matrix.txt.gz", "GPL10558.annot.gz",
                       lambda L: np.array([int(x[0]) for x in L["subtype"]]))
L2, _, _ = mat("GSE74385_matrix.txt.gz")
km = np.array([not x[1:].startswith("M") for x in L2["subtype"]])
V2 = V2r[:, km]; gv2 = gva[km]

records = []
for name, X, keep, gv in [("GSE16581", V1, keep1, gv1), ("GSE74385", V2, keep2, gv2)]:
    genes = list(keep); ridx = np.array([keep[g] for g in genes]); Xv = X[ridx]
    ev = eff(Xv, gv)
    common = [(didx[g], j) for j, g in enumerate(genes) if g in didx]
    di = np.array([c[0] for c in common]); vi = np.array([c[1] for c in common])
    ev_s, ed_s = ev[vi], ed[di]
    me_s = me[di]
    fin = np.isfinite(ev_s) & np.isfinite(ed_s) & np.isfinite(me_s) & (ed_s != 0)
    idx = np.where(fin)[0]
    estr = np.quantile(me_s[idx], [0, .2, .4, .6, .8, 1.0])

    for e in range(5):
        lo, hi = estr[e], (estr[e + 1] if e < 4 else np.inf)
        ins = idx[(me_s[idx] >= lo) & (me_s[idx] < hi)]
        qb = np.quantile(np.abs(ed_s[ins]), [0, .2, .4, .6, .8, 1.0])
        ql = np.digitize(np.abs(ed_s[ins]), qb[1:-1], right=False)
        groups = [ins[ql == k] for k in range(5)]
        conc = [float(np.mean(np.sign(ev_s[g]) == np.sign(ed_s[g]))) for g in groups]
        obs = float(W @ np.array(conc))
        null = np.empty(NPERM)
        for b in range(NPERM):
            e2 = eff(Xv, rng.permutation(gv))[vi]
            null[b] = W @ np.array([np.mean(np.sign(e2[g]) == np.sign(ed_s[g])) for g in groups])
        p = (np.sum(null >= obs) + 1) / (NPERM + 1)
        rec = dict(cohort=name, stratum=f"E{e+1}", n=len(ins),
                   exp_lo=round(float(lo), 3), exp_hi=(round(float(hi), 3) if e < 4 else None),
                   trend_p=round(float(p), 4), delta_pct=round(100 * (conc[4] - conc[0]), 2))
        for k in range(5):
            rec[f"Q{k+1}_pct"] = round(100 * conc[k], 2)
        records.append(rec)
        print(f"{name} E{e+1}: n={len(ins):4d} " +
              " ".join(f"Q{k+1}={100*conc[k]:5.1f}%" for k in range(5)) +
              f" | d={100*(conc[4]-conc[0]):5.1f} | p={p:.3f}")

df = pd.DataFrame(records)
path = os.path.join(OUT, "fig4_strata_v2.csv")
df.to_csv(path, index=False, encoding="utf-8-sig")
print("\nsaved:", path)
