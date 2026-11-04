# -*- coding: utf-8 -*-
"""
ä¸¥æ ¼æ¯”è¾ƒï¼š|Î”| vs |t| vs ä¸‹ç½®ä¿¡ç•Œ
é…å¯¹ bootstrapï¼ˆæ‚£è€…çº§ï¼Œä¸¤é˜Ÿåˆ—åŒæ—¶é‡æŠ½æ ·ï¼‰ç»™å‡ºç­–ç•¥é—´å·®å€¼çš„ CI
"""
import sys, os, re, gzip
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
SEED, NPERM, NBOOT = 20261103, 1000, 1000
rng = np.random.default_rng(SEED)
KS = [100, 500, 1000, 2000]


def stats(X, g):
    hi, lo = X[:, g >= 2], X[:, g == 1]
    d = hi.mean(1) - lo.mean(1)
    se = np.sqrt(hi.var(1, ddof=1) / hi.shape[1] + lo.var(1, ddof=1) / lo.shape[1])
    return d, se


# discovery
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
d_t = d_d / (d_se + 1e-12)
mean_exp = Xd.mean(1)
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
dt = d_t[[didx[g] for g in common]]
dse = d_se[[didx[g] for g in common]]
me = mean_exp[[didx[g] for g in common]]
print("å…±åŒåŸºå› :", len(common))

e1 = A1[:, gv1 >= 2].mean(1) - A1[:, gv1 == 1].mean(1)
e2 = A2[:, gv2 >= 2].mean(1) - A2[:, gv2 == 1].mean(1)
valid = np.isfinite(e1) & np.isfinite(e2)
print("ä¸¤é˜Ÿåˆ—å‡æœ‰æœ‰æ•ˆæ•ˆåº”:", int(valid.sum()))

SCORES = {"A_absDelta": np.abs(dd), "B_absT": np.abs(dt),
          "C_shrunk": np.maximum(np.abs(dd) - 1.96 * dse, 0)}
sels = {}
for name, sc in SCORES.items():
    order = np.argsort(-np.where(valid, sc, -np.inf))
    for K in KS:
        sels[(name, K)] = order[:K]


def metrics(s1, s2, sel, ref):
    ok = sel[np.isfinite(s1[sel]) & np.isfinite(s2[sel])]
    if len(ok) < 10:
        return np.nan, np.nan, np.nan
    r = np.sign(ref[ok])
    c1 = np.mean(r == np.sign(s1[ok]))
    c2 = np.mean(r == np.sign(s2[ok]))
    b = np.mean((r == np.sign(s1[ok])) & (r == np.sign(s2[ok])))
    return c1, c2, b


rows = []
for K in KS:
    base = {n: metrics(e1, e2, sels[(n, K)], dd) for n in SCORES}
    # bootstrap
    boot = {n: [] for n in SCORES}
    for b in range(NBOOT):
        i1 = rng.integers(0, len(gv1), len(gv1))
        i2 = rng.integers(0, len(gv2), len(gv2))
        if (gv1[i1] >= 2).sum() < 3 or (gv1[i1] == 1).sum() < 3:
            continue
        if (gv2[i2] >= 2).sum() < 3 or (gv2[i2] == 1).sum() < 3:
            continue
        b1 = A1[:, i1][:, gv1[i1] >= 2].mean(1) - A1[:, i1][:, gv1[i1] == 1].mean(1)
        b2 = A2[:, i2][:, gv2[i2] >= 2].mean(1) - A2[:, i2][:, gv2[i2] == 1].mean(1)
        for n in SCORES:
            boot[n].append(metrics(b1, b2, sels[(n, K)], dd)[2])
    for n in SCORES:
        arr = np.array(boot[n], float)
        arr = arr[np.isfinite(arr)]
        c1, c2, bo = base[n]
        rows.append(dict(K=K, strategy=n, acc1=round(100 * c1, 1), acc2=round(100 * c2, 1),
                         both=round(100 * bo, 1),
                         ci_lo=round(100 * np.percentile(arr, 2.5), 1),
                         ci_hi=round(100 * np.percentile(arr, 97.5), 1)))
    for x, y in [("B_absT", "A_absDelta"), ("C_shrunk", "A_absDelta"), ("C_shrunk", "B_absT")]:
        d = np.array(boot[x]) - np.array(boot[y])
        d = d[np.isfinite(d)]
        rows.append(dict(K=K, strategy=f"{x}_minus_{y}", acc1=None, acc2=None,
                         both=round(100 * (base[x][2] - base[y][2]), 1),
                         ci_lo=round(100 * np.percentile(d, 2.5), 1),
                         ci_hi=round(100 * np.percentile(d, 97.5), 1),
                         p_greater=round(float(np.mean(d > 0)), 4)))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "strategy_bootstrap.csv"), index=False, encoding="utf-8-sig")
print(df.to_string(index=False))

