# -*- coding: utf-8 -*-
"""
审计：匹配时同时平衡 效应方向(sign) + |Δ| + 平均表达 + SE
并检查基因集内部相关性结构
"""
import sys, os, re, gzip, json, urllib.request
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", r"C:\work\topic-audit\data")
SEED, NBOOT = 20261103, 1000
rng = np.random.default_rng(SEED)


def stats(X, g):
    hi, lo = X[:, g >= 2], X[:, g == 1]
    d = hi.mean(1) - lo.mean(1)
    se = np.sqrt(hi.var(1, ddof=1) / hi.shape[1] + lo.var(1, ddof=1) / lo.shape[1])
    return d, se


GMT = os.path.join(OUT, "h.all.v2023.2.Hs.symbols.gmt")
sets = {}
for line in open(GMT, encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    sets[f[0]] = set(x for x in f[2:] if x)
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
    if m:
        m2g["M" + m.group(1)] = int(lab["grade"][i])
cols = [c for c in ddf.columns if c in m2g]
gd = np.array([m2g[c] for c in cols])
Xd = np.log2(ddf[cols].values + 1)
d_d, d_se = stats(Xd, gd)
d_me = Xd.mean(1)
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

common = [g for g in keep1 if g in keep2 and g in didx and np.isfinite(d_d[didx[g]]) and d_d[didx[g]] != 0]
A1 = V1[[keep1[g] for g in common]]; A2 = V2[[keep2[g] for g in common]]
dd = d_d[[didx[g] for g in common]]; dse = d_se[[didx[g] for g in common]]; dme = d_me[[didx[g] for g in common]]
e1 = A1[:, gv1 >= 2].mean(1) - A1[:, gv1 == 1].mean(1)
e2 = A2[:, gv2 >= 2].mean(1) - A2[:, gv2 == 1].mean(1)
valid = np.isfinite(e1) & np.isfinite(e2)

idx_prog = np.array([i for i, g in enumerate(common) if g in PROG and valid[i]])
idx_ctrl = np.array([i for i, g in enumerate(common) if g not in PROG and valid[i]])
sign = np.sign(dd)

# ---------- 方向 + 效应量 + 表达 + SE 四维匹配 ----------
Z = np.column_stack([np.abs(dd), dme, dse])
Zs = (Z - Z.mean(0)) / Z.std(0)

order = idx_prog[np.argsort(-np.abs(dd[idx_prog]))]
used = set(); matched = {}
for i in order:
    same = idx_ctrl[(sign[idx_ctrl] == sign[i])]
    cand = np.array([j for j in same if j not in used])
    if len(cand) == 0:
        continue
    j = int(cand[np.argmin(np.sqrt(((Zs[cand] - Zs[i]) ** 2).sum(1)))])
    matched[i] = j; used.add(j)
P = np.array(sorted(matched)); C = np.array([matched[i] for i in P])
print(f"方向平衡匹配: {len(P)} 对 | 未匹配 {len(idx_prog)-len(P)}")
print(f"程序基因上调比例 {100*np.mean(sign[P]>0):.1f}% | 对照 {100*np.mean(sign[C]>0):.1f}%")
print(f"中位|Δ| {np.median(np.abs(dd[P])):.3f} vs {np.median(np.abs(dd[C])):.3f} | "
      f"中位表达 {np.median(dme[P]):.2f} vs {np.median(dme[C]):.2f} | "
      f"中位SE {np.median(dse[P]):.3f} vs {np.median(dse[C]):.3f}")


def conc(s1, s2, sel):
    ok = sel[np.isfinite(s1[sel]) & np.isfinite(s2[sel])]
    r = np.sign(dd[ok])
    return (np.mean(r == np.sign(s1[ok])), np.mean(r == np.sign(s2[ok])),
            np.mean((r == np.sign(s1[ok])) & (r == np.sign(s2[ok]))))


p1, p2, pb = conc(e1, e2, P); c1, c2, cb = conc(e1, e2, C)
print("\n=== 方向平衡后的观测 ===")
print(f"增殖程序 (n={len(P)}): {100*p1:.1f}% / {100*p2:.1f}% / 同时 {100*pb:.1f}%")
print(f"匹配对照 (n={len(C)}): {100*c1:.1f}% / {100*c2:.1f}% / 同时 {100*cb:.1f}%")
print(f"差值: {100*(p1-c1):+.1f} / {100*(p2-c2):+.1f} / {100*(pb-cb):+.1f} 百分点")

d = []
for b in range(NBOOT):
    i1 = rng.integers(0, len(gv1), len(gv1)); i2 = rng.integers(0, len(gv2), len(gv2))
    if (gv1[i1] >= 2).sum() < 3 or (gv1[i1] == 1).sum() < 3: continue
    if (gv2[i2] >= 2).sum() < 3 or (gv2[i2] == 1).sum() < 3: continue
    b1 = A1[:, i1][:, gv1[i1] >= 2].mean(1) - A1[:, i1][:, gv1[i1] == 1].mean(1)
    b2 = A2[:, i2][:, gv2[i2] >= 2].mean(1) - A2[:, i2][:, gv2[i2] == 1].mean(1)
    x = conc(b1, b2, P); y = conc(b1, b2, C)
    d.append(x[2] - y[2])
d = np.array(d) * 100
print(f"同时同向差值 95% CI: [{np.percentile(d,2.5):+.1f}, {np.percentile(d,97.5):+.1f}] | P(>0)={np.mean(d>0):.3f}")

# ---------- 基因集内部相关性 ----------
def mean_abs_corr(Xsub):
    if Xsub.shape[0] < 3: return np.nan
    C_ = np.corrcoef(Xsub)
    iu = np.triu_indices_from(C_, k=1)
    return float(np.nanmean(np.abs(C_[iu])))

r_prog = mean_abs_corr(Xd[[didx[common[i]] for i in P]])
r_ctrl = mean_abs_corr(Xd[[didx[common[i]] for i in C]])
print("\n=== 组内表达相关性（发现队列）===")
print(f"增殖程序 平均|r| = {r_prog:.3f} | 匹配对照 {r_ctrl:.3f}")

json.dump(dict(n=int(len(P)), dir_balanced=True,
               prog=[round(100*p1,2), round(100*p2,2), round(100*pb,2)],
               ctrl=[round(100*c1,2), round(100*c2,2), round(100*cb,2)],
               diff=[round(100*(p1-c1),2), round(100*(p2-c2),2), round(100*(pb-cb),2)],
               ci=[round(np.percentile(d,2.5),2), round(np.percentile(d,97.5),2)],
               corr_prog=round(r_prog,3), corr_ctrl=round(r_ctrl,3)),
          open(os.path.join(OUT, "proliferation_dir_balanced.json"), "w"), indent=1)
print("\nsaved proliferation_dir_balanced.json")
