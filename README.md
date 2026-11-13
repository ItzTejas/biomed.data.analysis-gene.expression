# sigaudit

**Does a published gene-expression signature actually hold up?**

Every year, over a thousand papers report a new "prognostic signature": a short
list of genes whose expression is supposed to predict how long cancer patients
survive. Few are ever checked by anyone other than their authors. `sigaudit`
runs the checks that matter on public patient data, and prints an audit card.

```bash
pip install -e .
sigaudit run --file examples/oncotype_rs_linear.json --subset "er=positive" \
    --cohort tcga:BRCA@PFI metabric@DSS geo:GSE7390@DMFS
```

## What it checks

For each cohort the signature was **not** built on:

| Check | Question | Why |
|---|---|---|
| **Discrimination** | Does the score rank patients by outcome? (C-index, bootstrap 95% CI, HR per SD) | The basic claim. |
| **Random genes** | Does it beat 1,000 random gene sets of the same size and weights? | Venet et al. (2011) showed most random sets of breast-cancer genes predict outcome. Beating chance isn't enough; it has to beat *random signatures*. |
| **Proliferation** | Does it still predict after adjusting for a meta-PCNA proliferation score? | Much of the prognostic signal in cancer transcriptomes is just how fast tumours divide. A signature that adds nothing beyond that isn't telling us anything new. |
| **Transfer** | Does it hold across cohorts, platforms and countries? (random-effects meta-analysis, I²) | A signature that only works in the data it was fitted on isn't a finding. |
| **Who it works for** | Is it reliably worse for some patients (genetic ancestry, self-reported race where that is all a cohort has, sex, age)? | Tools built mostly on European-ancestry cohorts can underperform for everyone else, e.g. the 21-gene Recurrence Score in Black women ([Hoskins et al., JAMA Oncol 2021](https://doi.org/10.1001/jamaoncol.2020.7320)). |

The verdict is **holds up**, **mixed**, **does not hold up**, or **untested**.
Subgroup gaps are flagged only when the whole bootstrap interval is below zero,
and groups with fewer than 10 events are reported but never judged.

## Example: two breast cancer signatures

Full cards are in [`examples/cards/`](examples/cards/). Open the `.html` files for plots.

**PAM50 proliferation index (11 genes)**: *mixed*. It predicts outcome in all
four cohorts, but after adjusting for proliferation it adds nothing (0 of 4),
which is the right answer for a signature that *is* a proliferation measure.
This is the sanity check that the proliferation baseline works.

**21-gene Recurrence Score, linear approximation, ER-positive patients**: *mixed*.
Strong in METABRIC (C = 0.67) and TRANSBIG (C = 0.72), and still significant in
METABRIC after the proliferation adjustment. In TCGA it is weaker (C = 0.56)
and doesn't beat random genes. In METABRIC it works less well for women 70 and
over (C 0.62 vs 0.69). TCGA has only 10 progression events among ER-positive
patients of African ancestry, too few to test the Hoskins et al. finding.
That's a data gap, and the card says so instead of guessing.

> This is an open approximation of the published formula, not the commercial
> assay (which uses RT-PCR, reference genes and score thresholds).

## Cohorts

All public, no login, downloaded on first use and cached.

| Spec | Data | Patients | Endpoints |
|---|---|---|---|
| `tcga:<CODE>` (33 cancers) | UCSC Xena RNA-seq · TCGA-CDR survival (Liu et al., Cell 2018) · genetic ancestry (Carrot-Zhang et al., Cancer Cell 2020) | e.g. BRCA 1,076 | OS, DSS, PFI, DFI |
| `metabric` | cBioPortal datahub, Illumina HT-12 | 1,979 breast | OS, DSS, RFS |
| `geo:GSE7390` | Breast, node-negative untreated (TRANSBIG) | 198 | OS, DMFS, RFS |
| `geo:GSE11121` | Breast, node-negative untreated (Mainz) | 200 | DMFS |
| `geo:GSE31210` | Lung adenocarcinoma, stage I–II, Japan | 204 | OS, RFS |
| `geo:GSE68465` | Lung adenocarcinoma, US multi-site (with self-reported race) | 442 | OS |
| `geo:GSE72094` | Lung adenocarcinoma, Moffitt | 398 | OS |
| `geo:GSE39582` | Colon, France | 573 | OS, RFS |
| `geo:GSE17536` | Colorectal, Moffitt (with self-reported race) | 177 | OS, DSS |
| `geo:GSE84437` | Gastric, Korea | 431 | OS |
| `geo:GSE65858` | Head and neck, Germany | 270 | OS, PFS |
| `geo:GSE41613` | Oral cancer (HPV-negative), US | 97 | OS |
| `geo:GSE13507` | Bladder, Korea | 165 | OS, DSS |
| `geo:GSE32894` | Urothelial, Sweden | 224 | DSS |
| `geo:GSE62452` | Pancreatic, US | 66 | OS |
| `geo:GSE57495` | Pancreatic, Moffitt | 63 | OS |
| `geo:GSE26712` | Ovarian (high-grade serous), US | 185 | DSS |
| `geo:GSE140082` | Ovarian, ICON7 trial standard arm | 380 | OS, PFS |
| `geo:GSE76427` | Liver (HCC), Singapore | 115 | OS, RFS |
| `geo:GSE65904` | Melanoma metastases, Sweden | 210 | DSS |

Set `SIGAUDIT_DATA` to choose where downloads go (default `~/.cache/sigaudit`).
Adding a GEO cohort is one entry in `loaders/geo.py::REGISTRY`.

## Signatures in, three ways

```bash
sigaudit run --genes "MKI67,TOP2A,BIRC5" ...                       # unweighted
sigaudit run --genes "CDK1:0.21,TP53:-0.40" ...                    # weighted
sigaudit run --genes "Risk score = 0.21*CDK1 - 0.40*TP53" ...      # pasted from a paper
```

Gene names are mapped to current HGNC symbols, so old names (`CDCA1` → `NUF2`),
Ensembl and Entrez ids, and Excel-mangled names (`1-Sep` → `SEPTIN1`) still match.
Use `--training tcga:BRCA` to mark the cohort a signature was fitted on; it's
reported but kept out of the verdict.

## How it relates to existing tools

The random-signature comparison is not new: [SigCheck](https://bioconductor.org/packages/SigCheck/)
(Bioconductor) does it for one dataset at a time, and large one-off studies
have re-tested published signatures for a single cancer type (e.g. 42 lung
cancer signatures, Tang et al., Ann Oncol 2017). `sigaudit` puts the random,
proliferation, cross-cohort and subgroup checks together in one open Python
tool with built-in public cohorts. It's also the scoring engine for
[`sigwatch`](https://github.com/ItzTejas/sigwatch), which audits published
signatures at scale.

## Method notes

- Scores are the weighted mean of per-gene z-scores within each cohort, so
  signatures transfer across platforms. Missing genes are dropped; cohorts
  with less than half the genes are skipped.
- C-index is Harrell's. Cox models use Breslow ties. Both are written from
  scratch in numpy and tested against `lifelines`.
- The random null keeps the signature's size and weights, shuffled onto random
  genes. The signature is compared on |C − 0.5| (two-sided, like Venet et al.)
  and must also point the right way.
- The proliferation score is the mean of the ~1% of genes most correlated with
  a 17-gene proliferation core in that cohort. No outcome data are used.
- Follow-up is censored at 10 years by default (`--horizon`).

**Caveat**: a signature that fails these checks is not evidence of error or
misconduct by its authors. It means the claim didn't transfer to these
cohorts under these settings.

## Tests

```bash
pip install -e ".[test]"
pytest
```

Synthetic cohorts with known answers check that a real signal holds up, noise
doesn't, a proliferation proxy is caught, an ancestry gap is flagged, and the
survival statistics match `lifelines`.

## License

MIT
