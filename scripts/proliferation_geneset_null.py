# -*- coding: utf-8 -*-
"""补充：基因集层面的零分布——按效应量分层随机抽同等大小的基因集"""
import sys, os, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import importlib.util

sys.argv = ["x"]
spec = importlib.util.spec_from_file_location("pt", r"C:\work\topic-audit\proliferation_test.py")
# 不执行整个脚本，改为复用其输出所需的中间量：这里直接重算核心数组
exec(open(r"C:\work\topic-audit\proliferation_test.py", encoding="utf-8").read().split("# ---------- 1:1 最近邻匹配")[0])

SEED = 20261103
rng = np.random.default_rng(SEED)
NRAND = 1000

idx_prog = np.array([i for i, g in enumerate(common) if g in PROG and valid[i]])
idx_ctrl_pool = np.array([i for i, g in enumerate(common) if g not in PROG and valid[i]])

# 观测集合：用与主分析相同的匹配（降序匹配）
Z = np.column_stack([np.abs(dd), dme, dse])
Zs = (Z - Z.mean(0)) / Z.std(0)
order = idx_prog[np.argsort(-np.abs(dd[idx_prog]))]
pool = set(idx_ctrl_pool.tolist())
matched, used = {}, set()
for i in order:
    cand = np.array(sorted(pool - used))
    if not len(cand):
        break
    j = int(cand[np.argmin(np.sqrt(((Zs[cand] - Zs[i]) ** 2).sum(1)))])
    matched[i] = j
    used.add(j)
P = np.array(sorted(matched))

# 按 |Δ| 十分位分层，随机抽同大小集合
dec = pd.qcut(np.abs(dd[idx_ctrl_pool]), 10, labels=False, duplicates="drop")
bins = [idx_ctrl_pool[dec == k] for k in range(dec.max() + 1)]
prog_per_bin = np.bincount(pd.cut(np.abs(dd[P]), bins=[-np.inf] + list(pd.qcut(np.abs(dd[idx_ctrl_pool]), 10, retbins=True, duplicates="drop")[1][1:-1]) + [np.inf], labels=False), minlength=len(bins))

def both_conc(sel, s1, s2):
    ok = sel[np.isfinite(s1[sel]) & np.isfinite(s2[sel])]
    r = np.sign(dd[ok])
    return np.mean((r == np.sign(s1[ok])) & (r == np.sign(s2[ok])))

obs = both_conc(P, e1, e2)
null = np.empty(NRAND)
for b in range(NRAND):
    pick = []
    for k, arr in enumerate(bins):
        if len(pick) >= len(P):
            break
        need = int(round(len(P) / len(bins)))
        if len(arr):
            pick.extend(rng.choice(arr, size=min(need, len(arr)), replace=False))
    pick = np.array(pick[:len(P)])
    null[b] = both_conc(pick, e1, e2)

p = (np.sum(null >= obs) + 1) / (NRAND + 1)
print(f"观测（增殖程序，n={len(P)}）同时同向: {100*obs:.1f}%")
print(f"效应量分层随机集合（n={len(P)}, {NRAND} 次）: {100*null.mean():.1f}% ± {100*null.std():.1f}")
print(f"95 分位: {100*np.percentile(null,95):.1f}% | 经验 P = {p:.4f}")
json.dump(dict(obs=round(100*float(obs),2), null_mean=round(100*float(null.mean()),2),
               null_sd=round(100*float(null.std()),2), p=float(round(p,4))),
          open(os.path.join(OUT, "proliferation_geneset_null.json"), "w"), indent=1)
print("saved proliferation_geneset_null.json")
