# -*- coding: utf-8 -*-
"""敏感性：排除 E2F/G2M 基因后，原口径梯度是否仍然存在"""
import sys, os, re, gzip, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", r"C:\work\topic-audit\data")
SEED, NPERM = 20261103, 1000
rng = np.random.default_rng(SEED)
W = np.array([-2, -1, 0, 1, 2], float)

sets = {}
for line in open(os.path.join(OUT, "h.all.v2023.2.Hs.symbols.gmt"), encoding="utf-8"):
    f = line.rstrip("\n").split("\t"); sets[f[0]] = set(x for x in f[2:] if x)
PROG = sets["HALLMARK_E2F_TARGETS"] | sets["HALLMARK_G2M_CHECKPOINT"]

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
d_mix = Xd[:, gd >= 2].mean(1) - Xd[:, gd == 1].mean(1)
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
    L = {c[0].split(":", 1)[0].strip(): [x.split(":", 1)[-1].strip() for x in c] for c in chars}
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
    return V, keep, gv


V1r, keep1, gv1 = prep("GSE16581_matrix.txt.gz", "GPL570.annot.gz", lambda L: np.array([int(x) for x in L["who grade"]]))
V1 = np.log2(V1r + 1)   # 主规格：log2 变换
V2r, keep2, gva = prep("GSE74385_matrix.txt.gz", "GPL10558.annot.gz", lambda L: np.array([int(x[0]) for x in L["subtype"]]))
L2, _, _ = mat("GSE74385_matrix.txt.gz")
km = np.array([not x[1:].startswith("M") for x in L2["subtype"]])
V2 = V2r[:, km]; gv2 = gva[km]

out = []
for name, V, keep, gv in [("GSE16581", V1, keep1, gv1), ("GSE74385", V2, keep2, gv2)]:
    genes = [g for g in keep if g in didx]
    ridx = np.array([keep[g] for g in genes])
    e = V[ridx][:, gv >= 2].mean(1) - V[ridx][:, gv == 1].mean(1)
    ed = d_mix[[didx[g] for g in genes]]
    fin = np.isfinite(e) & np.isfinite(ed) & (ed != 0)
    isprog = np.array([g in PROG for g in genes])
    for tag, mask in [("全部基因", fin), ("排除 E2F/G2M", fin & ~isprog)]:
        idx = np.where(mask)[0]
        edges = np.quantile(np.abs(ed[idx]), [0, .2, .4, .6, .8, 1.0])
        ql = np.digitize(np.abs(ed[idx]), edges[1:-1], right=False)
        groups = [idx[ql == q] for q in range(5)]
        obs = np.array([np.mean(np.sign(e[g]) == np.sign(ed[g])) for g in groups])
        perm = np.empty(NPERM)
        for b in range(NPERM):
            gp = rng.permutation(gv)
            eb = V[ridx][:, gp >= 2].mean(1) - V[ridx][:, gp == 1].mean(1)
            perm[b] = W @ np.array([np.mean(np.sign(eb[g]) == np.sign(ed[g])) for g in groups])
        stat = float(W @ obs); p = (np.sum(perm >= stat) + 1) / (NPERM + 1)
        print(f"\n{name} | {tag} | n={mask.sum()} | 排除的增殖基因 {int((fin & isprog).sum())}")
        print("  一致率: " + "  ".join(f"{100*x:5.1f}" for x in obs))
        print(f"  趋势 {stat:+.4f} vs 零 {perm.mean():+.4f}±{perm.std():.4f} | P={p:.4f}")
        out.append(dict(cohort=name, tag=tag, n=int(mask.sum()),
                        obs=[round(100*x, 2) for x in obs], trend=round(stat, 4), p=round(float(p), 4)))

json.dump(out, open(os.path.join(OUT, "prog_excluded_gradient_v2.json"), "w"), indent=1)
print("\nsaved prog_excluded_gradient.json")
