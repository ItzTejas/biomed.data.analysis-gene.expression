# Audit card: 21-gene Recurrence Score (linear approximation)

Source: Paik et al., NEJM 2004, doi:10.1056/NEJMoa041588  
16 genes (weighted)

**Verdict: mixed**

- Points the right way (C > 0.5) in 3 of 3 independent cohort(s).
- Beats size-matched random gene sets in 2 of 3.
- Still predicts after adjusting for proliferation in 1 of 3.
- Pooled HR per SD 1.49 (1.18-1.87), heterogeneity I² = 76%.
- metabric: C-index 0.62 in age_band=70+ vs 0.69 in 50-69 (gap -0.07, 95% CI -0.13 to -0.01).

| Cohort | Endpoint | n / events | Genes | C-index (95% CI) | vs random | Prolif. C | HR per SD | HR adj. prolif. | HR adj. clinical |
|---|---|---|---|---|---|---|---|---|---|
| tcga:BRCA | PFI | 792 / 87 | 16/16 | 0.561 (0.497–0.631) | 83th pct, p=0.169 | 0.538 | 1.20 (0.99–1.46), p=0.064 | 1.13 (0.79–1.63), p=0.5 | 1.19 (0.98–1.46), p=0.08 (age, stage) |
| metabric | DSS | 1505 / 335 | 16/16 | 0.669 (0.642–0.695) | 100th pct, p=0.001 | 0.653 | 1.64 (1.50–1.80), p=3.2e-26 | 1.45 (1.25–1.68), p=8.9e-07 | 1.65 (1.50–1.81), p=2.7e-26 (age) |
| geo:GSE7390 | DMFS | 134 / 26 | 15/16 | 0.720 (0.651–0.801) | 100th pct, p=0.002 | 0.721 | 1.71 (1.20–2.43), p=0.0028 | 1.24 (0.67–2.30), p=0.5 | 1.79 (1.25–2.57), p=0.0016 (age) |

## Who it works for

| Cohort | Attribute | Group | n / events | C-index | Gap vs reference (95% CI) |
|---|---|---|---|---|---|
| tcga:BRCA | ancestry | AFR | 102 / 10 | 0.609 | +0.101 (-0.058 to +0.254) |
| tcga:BRCA | ancestry | AMR | 13 / 1 | – | too few events to judge |
| tcga:BRCA | ancestry | EAS | 26 / 2 | 0.259 | too few events to judge |
| tcga:BRCA | ancestry | EUR | 589 / 65 | 0.507 | reference |
| tcga:BRCA | age_band | 50-69 | 414 / 43 | 0.582 | reference |
| tcga:BRCA | age_band | 70+ | 178 / 18 | 0.527 | -0.055 (-0.246 to +0.106) |
| tcga:BRCA | age_band | <50 | 200 / 26 | 0.564 | -0.018 (-0.182 to +0.116) |
| metabric | age_band | 50-69 | 786 / 171 | 0.687 | reference |
| metabric | age_band | 70+ | 474 / 118 | 0.615 | -0.072 (-0.132 to -0.006) ⚠ |
| metabric | age_band | <50 | 245 / 46 | 0.735 | +0.048 (-0.038 to +0.120) |
| geo:GSE7390 | age_band | 50-69 | 50 / 11 | 0.710 | -0.012 (-0.176 to +0.169) |
| geo:GSE7390 | age_band | <50 | 84 / 15 | 0.722 | reference |

## Gene matching

Not measured in at least one cohort: CD68

## How to read this

- **C-index**: probability the score ranks two patients in the right order (0.5 = coin flip).
- **vs random**: where the signature's |C − 0.5| falls among 1000 random gene sets of the same size and weights (Venet et al. 2011).
- **HR adj. prolif.**: hazard ratio per SD of the score after adjusting for a meta-PCNA proliferation score. If this is ~1, the signature mostly measures how fast tumours divide.
- Subgroup gaps are flagged ⚠ only when the whole 95% bootstrap interval is below zero.
- Follow-up censored at 10.0 years.
