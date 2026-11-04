# Data dictionary

## fig1_cohorts.csv
| column | meaning |
|---|---|
| cohort | GEO accession |
| grade | WHO grade (1/2/3) |
| n | number of tumours |

## fig1_qc_rna.csv
| column | meaning |
|---|---|
| sample | discovery-cohort sample ID |
| grade | WHO grade |
| detected_genes | genes with TPM > 1 |
| median_corr | median Pearson correlation with other samples |
| flag | `low` if <12,000 detected genes or median correlation <0.86 |

## fig2_effects.csv
| column | meaning |
|---|---|
| gene | gene symbol |
| effect | grade 2-3 minus grade 1, mean log2(TPM+1) |
| abs_effect | absolute value |
| quintile | discovery-defined effect quintile (1-5) |

## fig3_final_v3.csv
| column | meaning |
|---|---|
| cohort, quintile, n | cohort, quintile label, genes in quintile |
| eff_lo, eff_hi | absolute-effect boundaries |
| obs_pct | observed directional concordance (%) |
| null_mean_pct, null_sd_pct | permutation null mean and SD (%) |
| p_perm | empirical one-sided permutation P (10,000 permutations) |
| ci_lo_pct, ci_hi_pct | 95% patient-bootstrap interval |
| trend_obs, trend_p | linear contrast value and its permutation P |
| p_BH | Benjamini-Hochberg q across ten quintile tests |

## fig4_strata.csv
| column | meaning |
|---|---|
| stratum | expression stratum E1-E5 |
| exp_lo, exp_hi | mean-expression boundaries |
| Q1_pct ... Q5_pct | concordance within each effect quintile |
| delta_pct | Q5 minus Q1, percentage points |
| trend_p | stratum-specific trend P (500 permutations, unadjusted) |

## strategy_bootstrap.csv
| column | meaning |
|---|---|
| K | candidate-set size |
| strategy | `A_absDelta`, `B_absT`, `C_shrunk`, or `X_minus_Y` |
| both | simultaneous two-cohort concordance (%) |
| ci_lo, ci_hi | 95% paired bootstrap interval |
| p_greater | proportion of bootstrap replicates with difference > 0 |
