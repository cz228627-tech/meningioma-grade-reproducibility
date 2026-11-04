# -*- coding: utf-8 -*-
"""诊断：E1-Q1 的极低一致率是否由「效应恰好为 0」导致"""
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
me = Xd.mean(1)
didx = {g: i for i, g in enumerate(ddf.index)}
print("发现队列效应恰好为 0 的基因数:", int((ed == 0).sum()), "/", len(ed))


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


V1, keep1, gv1 = prep("GSE16581_matrix.txt.gz", "GPL570.annot.gz",
                      lambda L: np.array([int(x) for x in L["who grade"]]))
V2r, keep2, gva = prep("GSE74385_matrix.txt.gz", "GPL10558.annot.gz",
                       lambda L: np.array([int(x[0]) for x in L["subtype"]]))
L2, _, _ = mat("GSE74385_matrix.txt.gz")
km = np.array([not x[1:].startswith("M") for x in L2["subtype"]])
V2 = V2r[:, km]; gv2 = gva[km]


def diag(name, X, keep, gv):
    genes = list(keep); ridx = np.array([keep[g] for g in genes])
    Xv = X[ridx]; ev = eff(Xv, gv)
    common = [(g, didx[g], j) for j, g in enumerate(genes) if g in didx]
    di = np.array([c[1] for c in common]); vi = np.array([c[2] for c in common])
    ev_s, ed_s = ev[vi], ed[di]
    me_s = me[di]
    fin = np.isfinite(ev_s) & np.isfinite(ed_s)
    print(f"\n{'='*84}\n[{name}]")
    print("  效应恰好为 0：发现 %d | 验证 %d | 两者都 0: %d"
          % (int((ed_s[fin] == 0).sum()), int((ev_s[fin] == 0).sum()),
             int(((ed_s[fin] == 0) & (ev_s[fin] == 0)).sum())))

    absd = np.abs(ed_s[fin])
    edges = np.quantile(absd, [0, .2, .4, .6, .8, 1.0])
    ql = np.digitize(absd, edges[1:-1], right=False)
    idx = np.where(fin)[0]
    print(f"  {'分位':4s}{'n':>6s}{'ed==0':>7s}{'一致率(原)':>11s}{'剔除0后n':>9s}{'一致率(剔除0)':>13s}")
    for q in range(5):
        g_idx = idx[ql == q]
        zero = int((ed_s[g_idx] == 0).sum())
        c_all = np.mean(np.sign(ev_s[g_idx]) == np.sign(ed_s[g_idx]))
        keep_idx = g_idx[ed_s[g_idx] != 0]
        c_no0 = np.mean(np.sign(ev_s[keep_idx]) == np.sign(ed_s[keep_idx])) if len(keep_idx) else np.nan
        print(f"  Q{q+1:<3d}{len(g_idx):>6d}{zero:>7d}{100*c_all:>10.1f}%{len(keep_idx):>9d}{100*c_no0:>12.1f}%")

    # E1 层单独
    estr = np.quantile(me_s[fin], [0, .2, .4, .6, .8, 1.0])
    e1 = idx[me_s[idx] < estr[1]]
    q1 = np.quantile(np.abs(ed_s[e1]), [0, .2, .4, .6, .8, 1.0])
    sub = e1[np.abs(ed_s[e1]) <= q1[1]]
    zero = int((ed_s[sub] == 0).sum())
    print(f"  E1-Q1: n={len(sub)} 其中 ed==0: {zero} ({100*zero/len(sub):.1f}%)")
    if zero:
        k2 = sub[ed_s[sub] != 0]
        print(f"     剔除 ed==0 后 n={len(k2)}，一致率 {100*np.mean(np.sign(ev_s[k2])==np.sign(ed_s[k2])):.1f}%")


diag("GSE16581", V1, keep1, gv1)
diag("GSE74385", V2, keep2, gv2)
