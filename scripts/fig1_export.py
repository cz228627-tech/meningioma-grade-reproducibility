# -*- coding: utf-8 -*-
"""Fig 1 æ•°æ®å¯¼å‡ºï¼šé˜Ÿåˆ—æž„æˆ + å‘çŽ°é˜Ÿåˆ—ï¼ˆRNA-seqï¼‰æ ·æœ¬çº§ QC"""
import sys, os, re, gzip
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

OUT = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(__file__), "..", "data"))

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
m2g, age = {}, {}
for i, t in enumerate(t2):
    m = re.search(r"Meningioma\s+(\d+)\s*\(", t)
    if m:
        m2g["M" + m.group(1)] = int(lab["grade"][i])
        age["M" + m.group(1)] = float(lab["age"][i])
cols = [c for c in ddf.columns if c in m2g]
V = ddf[cols].values
X = np.log2(V + 1)
C = np.corrcoef(X.T)
mc = (C.sum(1) - 1) / (X.shape[1] - 1)
qc = pd.DataFrame({
    "sample": cols,
    "grade": [m2g[c] for c in cols],
    "detected_genes": (V > 1).sum(0),
    "median_corr": np.round(mc, 4),
})
qc["flag"] = np.where((qc.detected_genes < 12000) | (qc.median_corr < 0.86), "low", "")
qc.to_csv(os.path.join(OUT, "fig1_qc_rna.csv"), index=False, encoding="utf-8-sig")
print("QC è¡¨:", qc.shape, "| æ ‡è®°ä¸ºä½Žè´¨é‡:", int((qc.flag == "low").sum()))
print(qc[qc.flag == "low"].to_string(index=False))

comp = []
for coh, counts in [("GSE183653", {1: 86, 2: 74, 3: 25}),
                    ("GSE16581", {1: 43, 2: 19, 3: 6}),
                    ("GSE74385", {1: 17, 2: 8, 3: 28})]:
    for g, n in counts.items():
        comp.append(dict(cohort=coh, grade=g, n=n))
cd = pd.DataFrame(comp)
cd.to_csv(os.path.join(OUT, "fig1_cohorts.csv"), index=False, encoding="utf-8-sig")
print("\n", cd.pivot(index="cohort", columns="grade", values="n").to_string())

