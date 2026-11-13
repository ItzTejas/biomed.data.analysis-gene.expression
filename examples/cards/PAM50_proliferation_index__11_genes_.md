# Audit card: PAM50 proliferation index (11 genes)

Source: Nielsen et al., Clin Cancer Res 2010, doi:10.1158/1078-0432.CCR-10-1282  
11 genes

**Verdict: mixed**

- Points the right way (C > 0.5) in 4 of 4 independent cohort(s).
- Beats size-matched random gene sets in 3 of 4.
- Still predicts after adjusting for proliferation in 0 of 4.
- Pooled HR per SD 1.50 (1.27-1.78), heterogeneity I² = 68%.

| Cohort | Endpoint | n / events | Genes | C-index (95% CI) | vs random | Prolif. C | HR per SD | HR adj. prolif. | HR adj. clinical |
|---|---|---|---|---|---|---|---|---|---|
| tcga:BRCA | PFI | 1076 / 139 | 11/11 | 0.564 (0.510–0.614) | 88th pct, p=0.119 | 0.554 | 1.23 (1.04–1.46), p=0.016 | 1.45 (0.69–3.06), p=0.33 | 1.26 (1.06–1.50), p=0.0084 (age, stage) |
| metabric | DSS | 1978 / 515 | 11/11 | 0.648 (0.625–0.668) | 100th pct, p=0.001 | 0.653 | 1.55 (1.43–1.69), p=4.9e-25 | 1.08 (0.76–1.53), p=0.67 | 1.57 (1.44–1.71), p=8.4e-25 (age) |
| geo:GSE7390 | DMFS | 198 / 50 | 10/11 | 0.628 (0.556–0.691) | 99th pct, p=0.007 | 0.645 | 1.47 (1.11–1.95), p=0.008 | 0.41 (0.12–1.38), p=0.15 | 1.50 (1.12–2.00), p=0.0057 (age) |
| geo:GSE11121 | DMFS | 200 / 41 | 10/11 | 0.678 (0.579–0.766) | 100th pct, p=0.006 | 0.689 | 2.08 (1.48–2.92), p=2.4e-05 | 0.96 (0.16–5.80), p=0.96 | – |

## Who it works for

| Cohort | Attribute | Group | n / events | C-index | Gap vs reference (95% CI) |
|---|---|---|---|---|---|
| tcga:BRCA | ancestry | AFR | 169 / 25 | 0.596 | +0.076 (-0.065 to +0.196) |
| tcga:BRCA | ancestry | AMR | 18 / 3 | 0.600 | too few events to judge |
| tcga:BRCA | ancestry | EAS | 48 / 6 | 0.468 | too few events to judge |
| tcga:BRCA | ancestry | EUR | 767 / 93 | 0.520 | reference |
| tcga:BRCA | age_band | 50-69 | 561 / 61 | 0.543 | reference |
| tcga:BRCA | age_band | 70+ | 221 / 28 | 0.573 | +0.030 (-0.119 to +0.174) |
| tcga:BRCA | age_band | <50 | 294 / 50 | 0.594 | +0.051 (-0.045 to +0.136) |
| metabric | age_band | 50-69 | 1022 / 252 | 0.644 | reference |
| metabric | age_band | 70+ | 532 / 138 | 0.600 | -0.045 (-0.096 to +0.010) |
| metabric | age_band | <50 | 424 / 125 | 0.689 | +0.044 (-0.008 to +0.095) |
| geo:GSE7390 | age_band | 50-69 | 67 / 19 | 0.672 | +0.058 (-0.082 to +0.200) |
| geo:GSE7390 | age_band | <50 | 131 / 31 | 0.614 | reference |

## Gene matching

Not measured in at least one cohort: NUF2

## How to read this

- **C-index**: probability the score ranks two patients in the right order (0.5 = coin flip).
- **vs random**: where the signature's |C − 0.5| falls among 1000 random gene sets of the same size and weights (Venet et al. 2011).
- **HR adj. prolif.**: hazard ratio per SD of the score after adjusting for a meta-PCNA proliferation score. If this is ~1, the signature mostly measures how fast tumours divide.
- Subgroup gaps are flagged ⚠ only when the whole 95% bootstrap interval is below zero.
- Follow-up censored at 10.0 years.
