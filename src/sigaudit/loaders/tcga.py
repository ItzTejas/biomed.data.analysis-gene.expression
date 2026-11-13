"""TCGA cohorts: UCSC Xena RNA-seq + TCGA Clinical Data Resource survival
(Liu et al., Cell 2018) + genetic ancestry (Carrot-Zhang et al., Cancer Cell 2020).

All files are public and need no login.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from ..cache import fetch, processed
from ..cohort import Cohort
from ..genes import normalize_index

XENA_CLIN = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{c}.sampleMap%2F{c}_clinicalMatrix"
XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{c}.sampleMap%2FHiSeqV2.gz"
CDR = ("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
       "Survival_SupplementalTable_S1_20171025_xena_sp")
ANCESTRY = {  # GDC publication page CCG-AIM-2020
    "broad": "https://api.gdc.cancer.gov/data/1fee3458-14ee-4b4b-964c-a05164b68066",
    "ucsf": "https://api.gdc.cancer.gov/data/fdfa536a-c3c8-405d-99d9-bc9375b5084c",
}
# The CDR authors' recommended endpoint per cancer (Liu et al. 2018, Table 3)
# is not the same everywhere; OS is always available and is the default here.
ENDPOINTS = ("OS", "DSS", "PFI", "DFI")


def _stage(s) -> float:
    if not isinstance(s, str):
        return np.nan
    m = re.search(r"stage\s*(iv|iii|ii|i)[abc]?\d?\b", s.lower())
    return {"i": 1, "ii": 2, "iii": 3, "iv": 4}[m.group(1)] if m else np.nan


def expression(cancer: str) -> pd.DataFrame:
    """Primary-tumour RNA-seq, log2(RSEM normalised count + 1), one column per patient."""
    out = processed(f"tcga_{cancer}_expr.parquet")
    if out.exists():
        return pd.read_parquet(out)
    path = fetch(XENA.format(c=cancer), f"TCGA.{cancer}.HiSeqV2.gz")
    df = pd.read_csv(path, sep="\t", index_col=0)
    tumour_code = "03" if cancer == "LAML" else "01"
    keep = [c for c in df.columns if c[13:15] == tumour_code]
    df = df[keep]
    df.columns = [c[:12] for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    df = normalize_index(df).astype("float32")
    df.to_parquet(out)
    return df


def ancestry() -> pd.Series:
    """Consensus genetic ancestry per patient from two independent TCGA call sets."""
    out = processed("tcga_ancestry.parquet")
    if out.exists():
        return pd.read_parquet(out)["ancestry"]
    b = pd.read_csv(fetch(ANCESTRY["broad"], "tcga_ancestry_broad.txt"), sep="\t")
    b = b.set_index(b["SampleID"].str[:12])["Ancestry_assignment"]
    u = pd.read_csv(fetch(ANCESTRY["ucsf"], "tcga_ancestry_ucsf.csv"))
    u = u.set_index("Patient_ID")["pam.ancestry.cluster"].replace({"ASIAN": "ASN"})
    b = b[~b.index.duplicated()]
    u = u[~u.index.duplicated()]
    ids = b.index.union(u.index)
    calls = {}
    for pid in ids:
        x, y = b.get(pid), u.get(pid)
        x = x if isinstance(x, str) else None
        y = y if isinstance(y, str) else None
        if x and y:
            if x == y or (y == "ASN" and x in ("EAS", "SAS")):
                calls[pid] = x
            # disagreement: leave unassigned
        elif x:
            calls[pid] = x
        elif y and y != "ASN":
            calls[pid] = y
    s = pd.Series(calls, name="ancestry")
    s.to_frame().to_parquet(out)
    return s


def clinical(cancer: str, endpoint: str = "OS") -> pd.DataFrame:
    if endpoint not in ENDPOINTS:
        raise ValueError(f"endpoint must be one of {ENDPOINTS}")
    cdr = pd.read_csv(fetch(CDR, "TCGA_CDR_survival.tsv"), sep="\t", low_memory=False)
    cdr = cdr[cdr["cancer type abbreviation"] == cancer]
    if "Redaction" in cdr:
        cdr = cdr[cdr["Redaction"].isna()]
    cdr = cdr.drop_duplicates("_PATIENT").set_index("_PATIENT")
    clin = pd.DataFrame(index=cdr.index)
    clin["time"] = pd.to_numeric(cdr[f"{endpoint}.time"], errors="coerce") / 365.25
    clin["event"] = pd.to_numeric(cdr[endpoint], errors="coerce")
    clin["age"] = pd.to_numeric(cdr["age_at_initial_pathologic_diagnosis"], errors="coerce")
    clin["stage"] = cdr["ajcc_pathologic_tumor_stage"].map(_stage)
    if clin["stage"].isna().mean() > 0.5:
        clin["stage"] = cdr["clinical_stage"].map(_stage)
    clin["sex"] = cdr["gender"].str.lower()
    clin["year"] = pd.to_numeric(cdr["initial_pathologic_dx_year"], errors="coerce")
    clin["race_reported"] = cdr["race"].str.lower()
    clin["ancestry"] = ancestry().reindex(clin.index)
    if cancer == "BRCA":
        x = pd.read_csv(fetch(XENA_CLIN.format(c=cancer), f"TCGA.{cancer}.clinicalMatrix.tsv"),
                        sep="\t", low_memory=False)
        x = x.assign(pid=x["sampleID"].str[:12]).drop_duplicates("pid").set_index("pid")
        for col, src in (("er", "breast_carcinoma_estrogen_receptor_status"),
                         ("pr", "breast_carcinoma_progesterone_receptor_status"),
                         ("her2", "lab_proc_her2_neu_immunohistochemistry_receptor_status")):
            if src in x:
                v = x[src].reindex(clin.index).astype(str).str.lower()
                clin[col] = v.where(v.isin(["positive", "negative"]))
    return clin


def load(cancer: str, endpoint: str = "OS") -> Cohort:
    cancer = cancer.upper()
    return Cohort(f"tcga:{cancer}", expression(cancer), clinical(cancer, endpoint),
                  endpoint=endpoint, platform="RNA-seq (Illumina HiSeq, RSEM)", cancer=cancer)
