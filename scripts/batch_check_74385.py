# -*- coding: utf-8 -*-
"""GSE74385: grade x batch structure check, within-batch permutation, and batch-adjusted effects."""
import sys, os, re, gzip, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))
SEED, NPERM = 20261103, 1000
rng = np.random.default_rng(SEED)
W = np.array([-2, -1, 0, 1, 2], float)


def eff(X, g):
    return X[:, g >= 2].mean(1) - X[:, g == 1].mean(1)


def platform(annot):
    p2g = {}
    with gzip.open(os.path.join(OUT, annot), "rt", encoding="utf-8", errors="replace") as f:
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


chars, rows = [], []
with gzip.open(os.path.join(OUT, "GSE74385_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
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
sub = L["subtype"]
batch = np.array(L["batch"])
grade = np.array([int(x[0]) for x in sub])
isM = np.array([x[1:].startswith("M") for x in sub])
keep = ~isM

gv = grade[keep]
bt = batch[keep]
print("=== 53 samples (M-prefixed excluded) by grade x batch ===")
tab = pd.crosstab(pd.Series(bt, name="batch"), pd.Series(gv, name="grade"))
print(tab.to_string())
print("\nGrade variation within each batch:")
for b in sorted(set(bt)):
    gs = sorted(set(gv[bt == b]))
    n1 = int(((bt == b) & (gv == 1)).sum())
    n23 = int(((bt == b) & (gv >= 2)).sum())
    print(f"  batch {b}: G1={n1:2d} G2+3={n23:2d} | grades present {gs}")

# design matrix
p2g = platform("GPL10558.annot.gz")
V = np.array([[float(x) if x.strip('"') not in ("", "NA", "null") else np.nan for x in r[1:]]
              for r in rows[1:]], dtype=float)[:, keep]
probes = [r[0].strip('"') for r in rows[1:]]
kmap = {}
for k, pr in enumerate(probes):
    g = p2g.get(pr)
    if not g:
        continue
    if g not in kmap or np.nanmean(V[k]) > np.nanmean(V[kmap[g]]):
        kmap[g] = k
genes = list(kmap)
ridx = np.array([kmap[g] for g in genes])
X = V[ridx]

# discovery cohort
ddf = pd.read_csv(os.path.join(OUT, "GSE183653_meningioma_tpm.csv.gz"), index_col=0)
ch2, t2 = [], []
with gzip.open(os.path.join(OUT, "RNA_GSE183653.txt.gz"), "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            t2 = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            ch2.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            break
lab2 = {c[0].split(":", 1)[0].strip(): [x.split(":", 1)[-1].strip() for x in c] for c in ch2}
m2g = {}
for i, t in enumerate(t2):
    m = re.search(r"Meningioma\s+(\d+)\s*\(", t)
    if m:
        m2g["M" + m.group(1)] = int(lab2["grade"][i])
cols = [c for c in ddf.columns if c in m2g]
gd = np.array([m2g[c] for c in cols])
ed = eff(np.log2(ddf[cols].values + 1), gd)
didx = {g: i for i, g in enumerate(ddf.index)}

common = [(g, didx[g], j) for j, g in enumerate(genes) if g in didx]
di = np.array([c[1] for c in common])
vi = np.array([c[2] for c in common])
ed_s = ed[di]


def run(Xin, label):
    ev = eff(Xin, gv)[vi]
    fin = np.isfinite(ev) & np.isfinite(ed_s) & (ed_s != 0)
    absd = np.abs(ed_s[fin])
    edges = np.quantile(absd, [0, .2, .4, .6, .8, 1.0])
    ql = np.digitize(absd, edges[1:-1], right=False)
    idx = np.where(fin)[0]
    qidx = [idx[ql == q] for q in range(5)]
    obs = np.array([np.mean(np.sign(ev[ix]) == np.sign(ed_s[ix])) for ix in qidx])
    stat = float(W @ obs)
    # global label permutation
    glob = np.empty(NPERM)
    for b in range(NPERM):
        e = eff(Xin, rng.permutation(gv))[vi]
        glob[b] = W @ np.array([np.mean(np.sign(e[ix]) == np.sign(ed_s[ix])) for ix in qidx])
    # within-batch stratified permutation
    strata = np.empty(NPERM)
    for b in range(NPERM):
        gp = gv.copy()
        for bb in sorted(set(bt)):
            m = np.where(bt == bb)[0]
            gp[m] = rng.permutation(gv[m])
        e = eff(Xin, gp)[vi]
        strata[b] = W @ np.array([np.mean(np.sign(e[ix]) == np.sign(ed_s[ix])) for ix in qidx])
    print(f"\n[{label}] quintile concordance {[round(100*x,2) for x in obs]}")
    print(f"  trend statistic {stat:+.4f}")
    print(f"  global permutation:  null {glob.mean():+.4f}+/-{glob.std():.4f} | p={(np.sum(glob>=stat)+1)/(NPERM+1):.4f}")
    print(f"  within-batch perm.:  null {strata.mean():+.4f}+/-{strata.std():.4f} | p={(np.sum(strata>=stat)+1)/(NPERM+1):.4f}")
    return dict(obs=[round(100 * x, 2) for x in obs], stat=round(stat, 4),
                p_global=round(float((np.sum(glob >= stat) + 1) / (NPERM + 1)), 4),
                p_strat=round(float((np.sum(strata >= stat) + 1) / (NPERM + 1)), 4))


res = {}
res["unadjusted"] = run(X, "batch-unadjusted")

# batch adjustment: subtract each gene's within-batch mean
Xa = X.copy()
for bb in sorted(set(bt)):
    m = bt == bb
    if m.sum() > 1:
        Xa[:, m] = X[:, m] - np.nanmean(X[:, m], axis=1, keepdims=True)
res["batch_adjusted"] = run(Xa, "batch-adjusted")

print("\n" + json.dumps(res, ensure_ascii=False, indent=2))

