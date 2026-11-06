# -*- coding: utf-8 -*-
"""
生物学增量检验：已知增殖程序（HALLMARK_E2F_TARGETS + G2M_CHECKPOINT）的跨平台
方向一致率，在匹配发现效应量/表达量/精度后，是否仍高于非增殖对照基因。

预注册口径（运行前固定）：
  - 程序基因 = E2F_TARGETS ∪ G2M_CHECKPOINT（MSigDB v2023.2.Hs）
  - 对照 = 同一共同基因空间中不属于上述两个集合的基因
  - 1:1 最近邻匹配，匹配变量 |Δ|、平均表达、SE，无放回
  - 结局：两验证队列各自与共同(同时同向)的方向一致率
  - 不确定性：患者级配对 bootstrap 1000 次
  - 无论正负，结果都报告
"""
import sys, os, re, gzip, json, urllib.request
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", r"C:\work\topic-audit\data")
SEED, NBOOT, NPERM = 20261103, 1000, 1000
rng = np.random.default_rng(SEED)


def stats(X, g):
    hi, lo = X[:, g >= 2], X[:, g == 1]
    d = hi.mean(1) - lo.mean(1)
    se = np.sqrt(hi.var(1, ddof=1) / hi.shape[1] + lo.var(1, ddof=1) / lo.shape[1])
    return d, se


# ---------- MSigDB ----------
GMT = os.path.join(OUT, "h.all.v2023.2.Hs.symbols.gmt")
if not os.path.exists(GMT):
    url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
    urllib.request.urlretrieve(url, GMT)
    print("downloaded", GMT)

sets = {}
for line in open(GMT, encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    sets[f[0]] = set(x for x in f[2:] if x)
E2F = sets["HALLMARK_E2F_TARGETS"]
G2M = sets["HALLMARK_G2M_CHECKPOINT"]
PROG = E2F | G2M
print(f"E2F {len(E2F)} | G2M {len(G2M)} | 合并 {len(PROG)}")

# ---------- 队列 ----------
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
    p2g = platform(annot)
    L, V, pr = mat(mtx)
    gv = grader(L)
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
V2 = V2r[:, km]
gv2 = gva[km]

common = [g for g in keep1 if g in keep2 and g in didx and np.isfinite(d_d[didx[g]]) and d_d[didx[g]] != 0]
A1 = V1[[keep1[g] for g in common]]
A2 = V2[[keep2[g] for g in common]]
dd = d_d[[didx[g] for g in common]]
dse = d_se[[didx[g] for g in common]]
dme = d_me[[didx[g] for g in common]]
e1 = A1[:, gv1 >= 2].mean(1) - A1[:, gv1 == 1].mean(1)
e2 = A2[:, gv2 >= 2].mean(1) - A2[:, gv2 == 1].mean(1)
valid = np.isfinite(e1) & np.isfinite(e2)
print("共同基因:", len(common), "| 两队列均有效:", int(valid.sum()))

idx_prog = np.array([i for i, g in enumerate(common) if g in PROG and valid[i]])
idx_ctrl_pool = np.array([i for i, g in enumerate(common) if g not in PROG and valid[i]])
print(f"程序基因（可用）: {len(idx_prog)} | 对照池: {len(idx_ctrl_pool)}")

# ---------- 1:1 最近邻匹配 ----------
Z = np.column_stack([np.abs(dd), dme, dse])
mu, sd = Z.mean(0), Z.std(0)
Zs = (Z - mu) / sd

# 按 |Δ| 降序逐个匹配
order = idx_prog[np.argsort(-np.abs(dd[idx_prog]))]
pool = set(idx_ctrl_pool.tolist())
matched = {}
used = set()
for i in order:
    cand = np.array(sorted(pool - used))
    if len(cand) == 0:
        break
    dist = np.sqrt(((Zs[cand] - Zs[i]) ** 2).sum(1))
    j = int(cand[np.argmin(dist)])
    matched[i] = j
    used.add(j)

P = np.array(sorted(matched.keys()))
C = np.array([matched[i] for i in P])
print(f"成功匹配 {len(P)} 对 | 未匹配 {len(idx_prog)-len(P)}")
print(f"中位标准化距离 {np.median([np.sqrt(((Zs[matched[i]]-Zs[i])**2).sum()) for i in P]):.3f}")
print(f"程序基因中位|Δ| {np.median(np.abs(dd[P])):.3f} vs 对照 {np.median(np.abs(dd[C])):.3f}")
print(f"程序基因中位表达 {np.median(dme[P]):.2f} vs 对照 {np.median(dme[C]):.2f}")
print(f"程序基因中位SE {np.median(dse[P]):.3f} vs 对照 {np.median(dse[C]):.3f}")


def conc(s1, s2, sel):
    ok = sel[np.isfinite(s1[sel]) & np.isfinite(s2[sel])]
    r = np.sign(dd[ok])
    return (np.mean(r == np.sign(s1[ok])),
            np.mean(r == np.sign(s2[ok])),
            np.mean((r == np.sign(s1[ok])) & (r == np.sign(s2[ok]))))


p1, p2, pb = conc(e1, e2, P)
c1, c2, cb = conc(e1, e2, C)
print("\n=== 观测 ===")
print(f"增殖程序 (n={len(P)}): GSE16581 {100*p1:.1f}% | GSE74385 {100*p2:.1f}% | 同时 {100*pb:.1f}%")
print(f"匹配对照 (n={len(C)}): GSE16581 {100*c1:.1f}% | GSE74385 {100*c2:.1f}% | 同时 {100*cb:.1f}%")
print("差值: {:+.1f} / {:+.1f} / {:+.1f} 百分点".format(100*(p1-c1), 100*(p2-c2), 100*(pb-cb)))

# ---------- 患者级配对 bootstrap ----------
d_b, d_1, d_2 = [], [], []
for b in range(NBOOT):
    i1 = rng.integers(0, len(gv1), len(gv1))
    i2 = rng.integers(0, len(gv2), len(gv2))
    if (gv1[i1] >= 2).sum() < 3 or (gv1[i1] == 1).sum() < 3:
        continue
    if (gv2[i2] >= 2).sum() < 3 or (gv2[i2] == 1).sum() < 3:
        continue
    b1 = A1[:, i1][:, gv1[i1] >= 2].mean(1) - A1[:, i1][:, gv1[i1] == 1].mean(1)
    b2 = A2[:, i2][:, gv2[i2] >= 2].mean(1) - A2[:, i2][:, gv2[i2] == 1].mean(1)
    x = conc(b1, b2, P); y = conc(b1, b2, C)
    d_1.append(x[0] - y[0]); d_2.append(x[1] - y[1]); d_b.append(x[2] - y[2])


def ci(a):
    a = np.array(a) * 100
    return np.percentile(a, 2.5), np.percentile(a, 97.5), np.mean(a > 0)


print("\n=== bootstrap 差值 95% CI ===")
for lab, a in [("GSE16581", d_1), ("GSE74385", d_2), ("同时同向", d_b)]:
    lo, hi, pg = ci(a)
    print(f"{lab:10s} {np.mean(a)*100:+.1f} [{lo:+.1f}, {hi:+.1f}] | P(差值>0)={pg:.3f}")

res = dict(n_prog=int(len(P)), n_ctrl=int(len(C)),
           prog=(round(100*p1, 2), round(100*p2, 2), round(100*pb, 2)),
           ctrl=(round(100*c1, 2), round(100*c2, 2), round(100*cb, 2)),
           diff=(round(100*(p1-c1), 2), round(100*(p2-c2), 2), round(100*(pb-cb), 2)),
           ci_both=[round(x, 2) for x in ci(d_b)[:2]], p_both=round(ci(d_b)[2], 4))
json.dump(res, open(os.path.join(OUT, "proliferation_test.json"), "w"), indent=1)
print("\nsaved proliferation_test.json")
