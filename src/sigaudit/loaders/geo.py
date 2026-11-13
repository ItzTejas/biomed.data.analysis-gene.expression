"""GEO microarray cohorts with survival, from series-matrix files.

GEO clinical annotations are free text and every submitter uses different
keys, so each cohort gets a small, hand-checked entry in REGISTRY saying which
fields hold time, event and covariates. Adding a cohort = adding an entry.
"""
from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..cache import fetch, processed
from ..cohort import Cohort
from ..genes import normalize


@dataclass
class GeoSpec:
    gse: str
    cancer: str
    description: str
    endpoints: dict[str, tuple[str, str, float]]  # name -> (time key, event key, units per year)
    age: str | None = None
    sex: str | None = None
    stage: str | None = None
    # an event value matches this (case-insensitive); "no death", "0", "alive" don't
    event_true: str = r"^(1|1\.0|yes|dead|death|deceased|relapsed?|true|dod)\b"
    exclude: dict[str, str] = field(default_factory=dict)  # key -> regex; matching samples dropped
    matrix: str | None = None  # series-matrix file name when a GSE spans platforms
    # column -> (key, value map or None to keep the lower-cased value)
    extra: dict[str, tuple[str, dict | None]] = field(default_factory=dict)


REGISTRY: dict[str, GeoSpec] = {s.gse: s for s in [
    GeoSpec("GSE7390", "BRCA", "TRANSBIG node-negative, untreated (Desmedt 2007), Affymetrix U133A",
            {"OS": ("t.os", "e.os", 365.25), "DMFS": ("t.dmfs", "e.dmfs", 365.25),
             "RFS": ("t.rfs", "e.rfs", 365.25)}, age="age",
            extra={"er": ("er", {"1": "positive", "0": "negative"})}),
    GeoSpec("GSE11121", "BRCA", "Mainz node-negative, untreated (Schmidt 2008), Affymetrix U133A",
            {"DMFS": ("t.dmfs", "e.dmfs", 12.0)}),
    GeoSpec("GSE31210", "LUAD", "Stage I-II lung adenocarcinoma, Japan (Okayama 2012), U133 Plus 2",
            {"OS": ("days before death/censor", "death", 365.25),
             "RFS": ("days before relapse/censor", "relapse", 365.25)},
            age="age (years)", sex="gender", stage="pathological stage",
            exclude={"exclude for prognosis analysis due to incomplete resection or adjuvant therapy": "exclude"}),
    GeoSpec("GSE39582", "COAD", "Colon cancer, CIT cohort France (Marisa 2013), U133 Plus 2",
            {"OS": ("os.delay (months)", "os.event", 12.0), "RFS": ("rfs.delay", "rfs.event", 12.0)},
            age="age.at.diagnosis (year)", sex="Sex", stage="tnm.stage"),
    GeoSpec("GSE68465", "LUAD", "Lung adenocarcinoma, NCI Director's Challenge, US multi-site (Shedden 2008), U133A",
            {"OS": ("months_to_last_contact_or_death", "vital_status", 12.0)},
            age="age", sex="Sex", extra={"race": ("race", None)}),
    GeoSpec("GSE72094", "LUAD", "Lung adenocarcinoma, Moffitt (Schabath 2016), Rosetta/Merck array",
            {"OS": ("survival_time_in_days", "vital_status", 365.25)},
            age="age_at_diagnosis", sex="gender", stage="Stage", extra={"race": ("race", None)}),
    GeoSpec("GSE17536", "COAD", "Colorectal cancer, Moffitt (Smith 2010), U133 Plus 2",
            {"OS": ("overall survival follow-up time", "overall_event (death from any cause)", 12.0),
             "DSS": ("dss_time", "dss_event (disease specific survival; death from cancer)", 12.0)},
            age="age", sex="gender", stage="ajcc_stage", extra={"race": ("ethnicity", None)}),
    GeoSpec("GSE84437", "STAD", "Gastric cancer, Korea (Yoon 2020), Illumina HT-12 v3",
            {"OS": ("duration overall survival", "death", 12.0)}, age="age", sex="Sex"),
    GeoSpec("GSE65858", "HNSC", "Head and neck squamous cell carcinoma, Germany (Wichmann 2015), Illumina HT-12 v4",
            {"OS": ("os", "os_event", 365.25), "PFS": ("pfs", "pfs_event", 365.25)},
            age="age", sex="gender", stage="uicc_stage"),
    GeoSpec("GSE41613", "HNSC", "HPV-negative oral cancer, US (Lohavanichbutr 2013), U133 Plus 2",
            {"OS": ("fu time", "vital", 12.0)}, sex="Sex", stage="tumor stage"),
    GeoSpec("GSE13507", "BLCA", "Bladder cancer, Korea (Kim 2010), Illumina Human-6 v2",
            {"OS": ("survival month", "overall survival", 12.0),
             "DSS": ("survival month", "cancer specific survival", 12.0)},
            age="AGE", sex="SEX"),
    GeoSpec("GSE32894", "BLCA", "Urothelial carcinoma, Sweden (Sjödahl 2012), Illumina HT-12 v3",
            {"DSS": ("time_to_dod_(months)", "dod_event_(yes/no)", 12.0)}, age="age", sex="gender"),
    GeoSpec("GSE62452", "PAAD", "Pancreatic ductal adenocarcinoma, US (Yang 2016), HuGene 1.0 ST",
            {"OS": ("survival months", "survival status", 12.0)}, stage="Stage"),
    GeoSpec("GSE57495", "PAAD", "Pancreatic ductal adenocarcinoma, Moffitt (Chen 2015), Rosetta/Merck array",
            {"OS": ("overall survival (month)", "vital.status", 12.0)}, stage="Stage"),
    GeoSpec("GSE26712", "OV", "High-grade serous ovarian cancer, US (Bonome 2008), U133A",
            {"DSS": ("survival years", "status", 1.0)}),
    GeoSpec("GSE140082", "OV", "Ovarian cancer, ICON7 trial standard arm (Kommoss 2017), Illumina WG-DASL",
            {"OS": ("final_ostm", "final_osid", 365.25), "PFS": ("final_pfstm", "final_pfsid", 365.25)},
            age="age", stage="figo_stage"),
    GeoSpec("GSE76427", "LIHC", "Hepatocellular carcinoma, Singapore (Grinchuk 2018), Illumina HT-12 v4",
            {"OS": ("duryears_os", "event_os", 1.0), "RFS": ("duryears_rfs", "event_rfs", 1.0)},
            age="age (years)", stage="tnm_staging_clinical",
            extra={"sex": ("gender (1=m, 2=f)", {"1": "male", "2": "female"})}),
    GeoSpec("GSE65904", "SKCM", "Melanoma metastases, Sweden (Cirenajwis 2015), Illumina HT-12 v4",
            {"DSS": ("disease specific survival in days", "disease specific survival (1=death, 0=alive)", 365.25)},
            age="age", sex="gender"),
]}


def _stub(acc: str) -> str:
    """GEO's directory grouping: GSE2034 -> GSE2nnn, GPL96 -> GPLnnn."""
    return re.sub(r"\d{1,3}$", "nnn", acc)


def _matrix_url(gse: str, name: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{_stub(gse)}/{gse}/matrix/{name}"


def _annot_url(gpl: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz"


def read_series_matrix(path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """-> (probe x sample values, per-sample characteristics table, platform id)."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        text = f.read()
    head, _, rest = text.partition("!series_matrix_table_begin")
    body = rest.partition("!series_matrix_table_end")[0]
    meta: dict[str, list[list[str]]] = {}
    for line in head.splitlines():
        if line.startswith("!Sample_"):
            k, *vals = line.split("\t")
            meta.setdefault(k, []).append([v.strip().strip('"') for v in vals])
    samples = meta["!Sample_geo_accession"][0]
    chars = pd.DataFrame(index=samples)
    for row in meta.get("!Sample_characteristics_ch1", []):
        for gsm, cell in zip(samples, row):
            if ":" in cell:
                k, v = cell.split(":", 1)
                chars.loc[gsm, k.strip()] = v.strip()
    chars["title"] = meta["!Sample_title"][0]
    platform = re.search(r'!Series_platform_id\t"?(GPL\d+)', head).group(1)
    values = pd.read_csv(io.StringIO(body.strip()), sep="\t", index_col=0)
    values.index = values.index.astype(str)
    return values.apply(pd.to_numeric, errors="coerce"), chars, platform


def probe_to_gene(gpl: str) -> pd.Series:
    out = processed(f"{gpl}_probe_gene.parquet")
    if out.exists():
        return pd.read_parquet(out)["gene"]
    try:
        path = fetch(_annot_url(gpl), f"{gpl}.annot.gz")
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            lines = [l for l in f if not l.startswith(("!", "#", "^"))]
    except Exception:  # no curated .annot for this platform: use the submitter's table
        path = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=data",
                     f"{gpl}.table.txt")
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = [l for l in f if not l.startswith(("!", "#", "^"))]
    a = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, on_bad_lines="skip")
    sym_col = next((c for c in ("Gene symbol", "GeneSymbol", "Symbol", "GENE_SYMBOL", "Gene Symbol",
                                "ILMN_Gene") if c in a), None)
    ent_col = next((c for c in ("Entrez_Gene_ID", "EntrezGeneID", "ENTREZ_GENE_ID", "Gene ID") if c in a), None)
    if sym_col is None and ent_col is None:
        raise ValueError(f"no gene symbol column in the {gpl} annotation")
    a["sym"] = a[sym_col] if sym_col else None
    if ent_col:  # Entrez ids are unambiguous; prefer them when present
        a["sym"] = a[ent_col].where(a[ent_col].notna() & (a[ent_col].astype(str).str.strip() != ""), a["sym"])
    a = a[["ID", "sym"]].dropna()
    a = a[~a["sym"].astype(str).str.contains("///")]     # probes hitting several genes are ambiguous
    a["gene"] = a["sym"].astype(str).str.split(".").str[0].map(normalize)
    s = a.dropna(subset=["gene"]).set_index("ID")["gene"]
    s.to_frame().to_parquet(out)
    return s


def expression(spec: GeoSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    out_e = processed(f"geo_{spec.gse}_expr.parquet")
    out_c = processed(f"geo_{spec.gse}_chars.parquet")
    if out_e.exists() and out_c.exists():
        return pd.read_parquet(out_e), pd.read_parquet(out_c)
    name = spec.matrix or f"{spec.gse}_series_matrix.txt.gz"
    vals, chars, gpl = read_series_matrix(fetch(_matrix_url(spec.gse, name), name))
    if np.nanmax(vals.to_numpy()) > 100:              # MAS5-scale data: log it
        vals = np.log2(vals.clip(lower=1))
    genes = probe_to_gene(gpl).reindex(vals.index)
    vals = vals[genes.notna()]
    genes = genes[genes.notna()]
    # one probe per gene: the one with the highest mean expression (a common, simple rule)
    means = vals.mean(axis=1)
    best = means.groupby(genes).idxmax()
    expr = vals.loc[best.values]
    expr.index = best.index
    expr = expr.astype("float32")
    expr.to_parquet(out_e)
    chars.astype(str).to_parquet(out_c)
    return expr, chars


def clinical(spec: GeoSpec, chars: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    if endpoint not in spec.endpoints:
        raise ValueError(f"{spec.gse} has endpoints {list(spec.endpoints)}")
    tkey, ekey, per_year = spec.endpoints[endpoint]
    c = chars.copy()
    for k, rx in spec.exclude.items():
        if k in c:
            c = c[~c[k].astype(str).str.contains(rx, case=False, na=False)]
    clin = pd.DataFrame(index=c.index)
    clin["time"] = pd.to_numeric(c[tkey], errors="coerce") / per_year
    ev = c[ekey].astype(str).str.strip().str.lower()
    clin["event"] = np.where(ev.str.match(spec.event_true), 1.0,
                             np.where(ev.isin(["nan", "na", "n/a", "", "none", "--"]), np.nan, 0.0))
    if spec.age and spec.age in c:
        clin["age"] = pd.to_numeric(c[spec.age], errors="coerce")
    if spec.sex and spec.sex in c:
        clin["sex"] = c[spec.sex].str.lower().str.strip()
    if spec.stage and spec.stage in c:
        st = c[spec.stage].astype(str).str.upper().str.extract(r"(IV|III|II|I|[1-4])")[0]
        clin["stage"] = st.map({"I": 1, "II": 2, "III": 3, "IV": 4, "1": 1, "2": 2, "3": 3, "4": 4})
    for col, (key, mapping) in spec.extra.items():
        if key in c:
            v = c[key].astype(str).str.strip()
            clin[col] = v.map(mapping) if mapping else v.str.lower().replace({"na": np.nan, "nan": np.nan})
    return clin


def load(gse: str, endpoint: str = "OS") -> Cohort:
    gse = gse.upper()
    if gse not in REGISTRY:
        raise ValueError(f"{gse} is not in the GEO registry yet; known: {sorted(REGISTRY)}")
    spec = REGISTRY[gse]
    if endpoint not in spec.endpoints:
        endpoint = next(iter(spec.endpoints))
    expr, chars = expression(spec)
    clin = clinical(spec, chars, endpoint)
    return Cohort(f"geo:{gse}", expr, clin, endpoint=endpoint, platform=spec.description.rsplit(", ", 1)[-1],
                  cancer=spec.cancer, notes=[spec.description])
