# -*- coding: utf-8 -*-
"""
三队列共同基因空间：两个验证队列是否同时与发现队列同向
- 基因集：三队列均有效、发现效应非零
- 分组：发现队列 |Δ| 五分位
- 结局：sign(D)==sign(V1) 且 sign(D)==sign(V2)
- 零分布：两个验证队列独立置换（保留各组样本量），10,000 次
- 敏感性：GSE74385 批次内置换
"""
import sys, os, re, gzip, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", r"C:\work\topic-audit\data")
SEED, NPERM, NBOOT = 20261103, 10000, 1000
rng = np.random.default_rng(SEED)
W = np.array([-2, -1, 0, 1, 2], float)


def stats(X, g):
    hi, lo = X[:, g >= 2], X[:, g == 1]
    d = hi.mean(1) - lo.mean(1)
    se = np.sqrt(hi.var(1, ddof=1) / hi.shape[1] + lo.var(1, ddof=1) / lo.shape[1])
    return d, se


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
d_d, d_se = stats(Xd, gd)
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


V1r, keep1, gv1, _ = prep("GSE16581_matrix.txt.gz", "GPL570.annot.gz",
                         lambda L: np.array([int(x) for x in L["who grade"]]))
V1 = np.log2(V1r + 1)   # 主规格：log2 变换
V2r, keep2, gva, L2 = prep("GSE74385_matrix.txt.gz", "GPL10558.annot.gz",
                           lambda L: np.array([int(x[0]) for x in L["subtype"]]))
km = np.array([not s[1:].startswith("M") for s in L2["subtype"]])
V2 = V2r[:, km]; gv2 = gva[km]; bt = np.array(L2["batch"])[km]

common = [g for g in keep1 if g in keep2 and g in didx
          and np.isfinite(d_d[didx[g]]) and d_d[didx[g]] != 0]
A1 = V1[[keep1[g] for g in common]]
A2 = V2[[keep2[g] for g in common]]
dd = d_d[[didx[g] for g in common]]
e1 = A1[:, gv1 >= 2].mean(1) - A1[:, gv1 == 1].mean(1)
e2 = A2[:, gv2 >= 2].mean(1) - A2[:, gv2 == 1].mean(1)
valid = np.isfinite(e1) & np.isfinite(e2)
print(f"共同基因 {len(common)} | 三队列均有效 {int(valid.sum())}")

idx = np.where(valid)[0]
edges = np.quantile(np.abs(dd[idx]), [0, .2, .4, .6, .8, 1.0])
ql = np.digitize(np.abs(dd[idx]), edges[1:-1], right=False)
groups = [idx[ql == q] for q in range(5)]


def both_rates(ev1, ev2):
    out = []
    for g in groups:
        r = np.sign(dd[g])
        out.append(np.mean((r == np.sign(ev1[g])) & (r == np.sign(ev2[g]))))
    return np.array(out)


obs = both_rates(e1, e2)
stat = float(W @ obs)
print("观测同时同向:", [round(100 * x, 2) for x in obs], f"| 趋势 {stat:+.4f}")

# ---- 置换：两队列独立置换 ----
perm_q = np.empty((NPERM, 5)); perm_t = np.empty(NPERM)
for b in range(NPERM):
    g1 = rng.permutation(gv1); g2 = rng.permutation(gv2)
    b1 = A1[:, g1 >= 2].mean(1) - A1[:, g1 == 1].mean(1)
    b2 = A2[:, g2 >= 2].mean(1) - A2[:, g2 == 1].mean(1)
    r = both_rates(b1, b2)
    perm_q[b] = r; perm_t[b] = W @ r
p_q = [(np.sum(perm_q[:, i] >= obs[i]) + 1) / (NPERM + 1) for i in range(5)]

# BH across 5
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1); return out


q_bh = bh(p_q)
p_trend = (np.sum(perm_t >= stat) + 1) / (NPERM + 1)
print(f"趋势 P={p_trend:.4f} | 零 {perm_t.mean():+.4f}±{perm_t.std():.4f}")
for i in range(5):
    print(f"  Q{i+1}: {100*obs[i]:5.2f}% | 零 {100*perm_q[:,i].mean():5.2f}±{100*perm_q[:,i].std():4.2f} | p={p_q[i]:.4f} q={q_bh[i]:.4f}")

# ---- bootstrap（specimen-level, 两队列同时重抽） ----
boot = []
for b in range(NBOOT):
    i1 = rng.integers(0, len(gv1), len(gv1)); i2 = rng.integers(0, len(gv2), len(gv2))
    if (gv1[i1] >= 2).sum() < 3 or (gv1[i1] == 1).sum() < 3: continue
    if (gv2[i2] >= 2).sum() < 3 or (gv2[i2] == 1).sum() < 3: continue
    b1 = A1[:, i1][:, gv1[i1] >= 2].mean(1) - A1[:, i1][:, gv1[i1] == 1].mean(1)
    b2 = A2[:, i2][:, gv2[i2] >= 2].mean(1) - A2[:, i2][:, gv2[i2] == 1].mean(1)
    boot.append(both_rates(b1, b2))
boot = np.array(boot) * 100
ci = [[round(np.percentile(boot[:, i], 2.5), 2), round(np.percentile(boot[:, i], 97.5), 2)] for i in range(5)]

rows = []
for i in range(5):
    rows.append(dict(quintile=f"Q{i+1}", n=int(len(groups[i])),
                     eff_lo=round(float(edges[i]), 4), eff_hi=round(float(edges[i + 1]), 4),
                     both_pct=round(100 * obs[i], 2),
                     null_mean=round(100 * float(perm_q[:, i].mean()), 2),
                     null_sd=round(100 * float(perm_q[:, i].std()), 2),
                     p=round(float(p_q[i]), 4), q=round(float(q_bh[i]), 4),
                     ci_lo=ci[i][0], ci_hi=ci[i][1]))
df = pd.DataFrame(rows)
df["trend_obs"] = round(stat, 4)
df["trend_p"] = round(float(p_trend), 4)
df.to_csv(os.path.join(OUT, "both_cohorts_v2.csv"), index=False, encoding="utf-8-sig")
print("\n" + df.to_string(index=False))

# ---- 敏感性：GSE74385 批次内置换 ----
perm_t2 = np.empty(2000)
for b in range(2000):
    g2 = gv2.copy()
    for bb in np.unique(bt):
        m = bt == bb
        g2[m] = rng.permutation(gv2[m])
    g1 = rng.permutation(gv1)
    b1 = A1[:, g1 >= 2].mean(1) - A1[:, g1 == 1].mean(1)
    b2 = A2[:, g2 >= 2].mean(1) - A2[:, g2 == 1].mean(1)
    perm_t2[b] = W @ both_rates(b1, b2)
p_t2 = (np.sum(perm_t2 >= stat) + 1) / (2000 + 1)
print(f"\n批次内敏感性（2000 次）: 趋势 {stat:+.4f} | 零 {perm_t2.mean():+.4f}±{perm_t2.std():.4f} | P={p_t2:.4f}")
json.dump(dict(n_common=int(valid.sum()), obs=[round(100 * x, 2) for x in obs],
               trend=round(stat, 4), p_trend=round(float(p_trend), 4),
               null_mean=round(float(perm_t.mean()), 4), null_sd=round(float(perm_t.std()), 4),
               ci=ci, p_batch=round(float(p_t2), 4)),
          open(os.path.join(OUT, "both_cohorts_summary_v2.json"), "w"), indent=1)
print("saved both_cohorts_concordance.csv / both_cohorts_summary.json")
