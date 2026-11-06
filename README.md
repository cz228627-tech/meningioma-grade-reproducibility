# Meningioma grade-associated transcriptional reproducibility

Code and derived data for: *How Reproducible Are WHO Grade-Associated Transcriptional Differences in Meningioma? A Three-Cohort Cross-Platform Analysis*.

## What this repository contains

- `scripts/` — analysis scripts in execution order
- `data/` — frozen derived tables and result summaries (raw GEO matrices are not redistributed)
- `docs/data_dictionary.md` — column definitions

Raw expression matrices are **not** redistributed here. Download them from GEO:

| Cohort | Accession | Platform |
|---|---|---|
| Discovery | GSE183653 | RNA-seq (HiSeq 4000) |
| Validation 1 | GSE16581 | Affymetrix U133 Plus 2.0 (GPL570) |
| Validation 2 | GSE74385 | Illumina HumanHT-12 v4 (GPL10558) |

Platform annotations: GPL570 and GPL10558 annotation tables from GEO. Gene sets: MSigDB v2023.2.Hs (`h.all`, symbols).

Place the downloaded files in `data/raw/` and set `DATA_DIR`:

```bash
export DATA_DIR=/path/to/data
```

## Execution order

```bash
# primary analysis
python scripts/fig1_export.py             # cohort composition + per-sample QC
python scripts/fig2_export.py             # discovery effect sizes (GSE16581-matched set)
python scripts/fig3_frozen.py             # primary quintile analysis, 10,000 permutations
python scripts/strata_save.py             # expression-stratified sensitivity (500 permutations)
python scripts/batch_check_74385.py       # batch adjustment + within-batch permutation
python scripts/strategy_bootstrap.py      # exploratory ranking-strategy comparison

# grade-definition sensitivity (WHO 2 vs 1)
python scripts/g2only.py                  # G2-only quintiles, batch sensitivity (1,000 permutations)

# proliferation-programme sensitivity
python scripts/prog_excluded_gradient.py  # gradient after removing HALLMARK E2F/G2M genes
python scripts/proliferation_test.py      # matched programme-vs-control comparison
python scripts/proliferation_dir_balanced.py   # direction-balanced matching + within-set correlation
python scripts/proliferation_geneset_null.py   # gene-set level null distribution

# figures
python scripts/figs_final.py              # regenerates Figures 1-4 from frozen CSVs
```

`fig3_frozen.py` takes roughly 25-30 minutes (10,000 permutations x 2 cohorts); all other scripts finish within minutes.

## Verification

`scripts/` and `data/` correspond one-to-one: each script writes the result file that the manuscript cites.
`data/fig3_reproduced.csv` is a re-run of the primary analysis published alongside `fig3_final_v3.csv` for direct comparison;
the two files are identical in all 120 numeric cells.


## Reproducibility notes

- Random seed: `20261103` for all permutation and bootstrap procedures.
- Genes with an exactly zero discovery effect are excluded **before** quintile construction; they carry no direction and would otherwise be scored as discordant.
- Concordance is compared against an empirical permutation null, never against an assumed 50%.
- Bootstrap intervals are conditional on the fixed discovery ranking.

## Reporting conventions in `data/`

- `fig3_final_v3.csv` — one row per cohort x quintile; `obs_pct` observed concordance, `null_mean_pct`/`null_sd_pct` permutation null, `p_perm` empirical P, `p_BH` BH-adjusted q across ten tests.
- `strategy_bootstrap.csv` — pairwise strategy differences with 95% bootstrap intervals.

## License

MIT (see `LICENSE`).

## Data file checksums (SHA-256)

| file | bytes | sha256 |
|---|---|---|
| `fig1_cohorts.csv` | 155 | `00f19f08821ce2351533fe475faa24b8cd79c057054f8e72c6bac85f17a2966b` |
| `fig1_qc_rna.csv` | 4,040 | `f334a385a31aefb7abf50d44864cceee852ba6f16a0ae496c1d10f8d6f62127b` |
| `fig2_effects.csv` | 532,170 | `365cf083dd70e3a6942f939354c0c65857757ca5574a40646e07546343073cde` |
| `fig3_final_v3.csv` | 999 | `56192af78efbbc5e6d798532955adf591f4f291269d1991286803791c2849a4a` |
| `fig4_strata.csv` | 789 | `feea7fe66197e9e0a5b7e148c2f75a138b995b7f4836504a3cef169dbaa8b7f3` |
| `strategy_bootstrap.csv` | 1,141 | `bf7976b59070600723998293844186821d98b656a294b73a04f8f17ba3e5fe3a` |

Verify with `sha256sum data/*.csv` (Linux/macOS) or `Get-FileHash data\*.csv -Algorithm SHA256` (PowerShell).
