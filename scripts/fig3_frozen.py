# -*- coding: utf-8 -*-
"""
fig3_frozen.py â€” å”¯ä¸€æƒå¨è¾“å‡ºï¼ˆä¸€æ¬¡è¿è¡Œäº§å‡ºæ‰€æœ‰æ•°å­—ï¼‰
åŒæ—¶ä¿å­˜ï¼šäº”åˆ†ä½ä¸€è‡´çŽ‡ã€é›¶åˆ†å¸ƒã€ç½®æ¢ pã€BH qã€bootstrap CIã€è¶‹åŠ¿ç»Ÿè®¡é‡
"""
import sys, os, re, gzip, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
SEED, NPERM, NBOOT = 20261103, 10000, 1000
rng = np.random.default_rng(SEED)
W = np.array([-2, -1, 0, 1, 2], float)


def eff(X, g):
    return X[:, g >= 2].mean(1) - X[:, g == 1].mean(1)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1); return out


# --- å‘çŽ°é˜Ÿåˆ— ---
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

rows, pvals = [], []
for name, X, keep, gv in [("GSE16581", V1, keep1, gv1), ("GSE74385", V2, keep2, gv2)]:
    genes = list(keep); ridx = np.array([keep[g] for g in genes]); Xv = X[ridx]
    ev = eff(Xv, gv)
    common = [(didx[g], j) for j, g in enumerate(genes) if g in didx]
    di = np.array([c[0] for c in common]); vi = np.array([c[1] for c in common])
    ev_s, ed_s = ev[vi], ed[di]
    fin = np.isfinite(ev_s) & np.isfinite(ed_s)
    n_nonfinite = int((~fin).sum())
    # å…³é”®å®šä¹‰ï¼šå‘çŽ°é˜Ÿåˆ—æ•ˆåº”æ°å¥½ä¸º 0 çš„åŸºå› æ²¡æœ‰æ–¹å‘ä¿¡æ¯ï¼Œå…ˆå‰”é™¤å†åˆ’åˆ†äº”åˆ†ä½ï¼ˆæ–¹æ¡ˆ Aï¼‰
    # è¯¥æŽ’é™¤åªä¾èµ–å‘çŽ°é˜Ÿåˆ—ï¼Œç½®æ¢ä¸­æ’å®šä¸å˜
    n_zero = int((ed_s[fin] == 0).sum())
    fin = fin & (ed_s != 0)
    absd = np.abs(ed_s[fin]); edges = np.quantile(absd, [0, .2, .4, .6, .8, 1.0])
    ql = np.digitize(absd, edges[1:-1], right=False)
    ia = np.where(fin)[0]; qidx = [ia[ql == q] for q in range(5)]

    obs = np.array([np.mean(np.sign(ev_s[ix]) == np.sign(ed_s[ix])) for ix in qidx])
    perm = np.zeros((NPERM, 5)); trend = np.zeros(NPERM)
    bootc = np.zeros((NBOOT, 5)); boott = np.zeros(NBOOT)
    for b in range(NPERM):
        e = eff(Xv, rng.permutation(gv))[vi]
        c = np.array([np.mean(np.sign(e[ix]) == np.sign(ed_s[ix])) for ix in qidx])
        perm[b] = c; trend[b] = W @ c
    n = len(gv)
    for b in range(NBOOT):
        ii = rng.integers(0, n, n); gb = gv[ii]
        if (gb >= 2).sum() < 3 or (gb == 1).sum() < 3:
            bootc[b] = np.nan; boott[b] = np.nan; continue
        e = Xv[:, ii][:, gb >= 2].mean(1) - Xv[:, ii][:, gb == 1].mean(1)
        c = np.array([np.mean(np.sign(e[vi][ix]) == np.sign(ed_s[ix])) for ix in qidx])
        bootc[b] = c; boott[b] = W @ c
    p_q = [(np.sum(perm[:, q] >= obs[q]) + 1) / (NPERM + 1) for q in range(5)]
    obs_stat = float(W @ obs)
    p_trend = (np.sum(trend >= obs_stat) + 1) / (NPERM + 1)
    for q in range(5):
        lo, hi = np.nanpercentile(bootc[:, q], [2.5, 97.5])
        pvals.append(p_q[q])
        rows.append(dict(cohort=name, quintile=f"Q{q+1}", n=int(len(qidx[q])),
                         eff_lo=round(float(edges[q]), 4), eff_hi=round(float(edges[q + 1]), 4),
                         obs_pct=round(100 * obs[q], 2),
                         null_mean_pct=round(100 * perm[:, q].mean(), 2),
                         null_sd_pct=round(100 * perm[:, q].std(), 2),
                         p_perm=round(float(p_q[q]), 4),
                         ci_lo_pct=round(100 * lo, 2), ci_hi_pct=round(100 * hi, 2),
                         trend_obs=round(obs_stat, 4), trend_p=round(float(p_trend), 4)))
    print(f"[{name}] è¶‹åŠ¿ {obs_stat:+.4f} vs é›¶ {trend.mean():+.4f}Â±{trend.std():.4f} | p={p_trend:.4f}")

df = pd.DataFrame(rows)
df["p_BH"] = np.round(bh(pvals), 4)
df.to_csv(os.path.join(OUT, "fig3_final_v3.csv"), index=False, encoding="utf-8-sig")
print("\n" + df.to_string(index=False))
print("\nå”¯ä¸€æƒå¨è¡¨: fig3_final_v3.csv")

