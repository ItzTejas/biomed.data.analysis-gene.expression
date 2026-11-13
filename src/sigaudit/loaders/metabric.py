"""METABRIC breast cancer cohort (Curtis et al., Nature 2012; Pereira et al.,
Nat Commun 2016), from the public cBioPortal datahub. Illumina HT-12 arrays,
~1,900 patients with up to 25 years of follow-up.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..cache import fetch, processed
from ..cohort import Cohort
from ..genes import normalize_index

BASE = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/"
ENDPOINTS = ("OS", "DSS", "RFS")


def expression() -> pd.DataFrame:
    out = processed("metabric_expr.parquet")
    if out.exists():
        return pd.read_parquet(out)
    path = fetch(BASE + "data_mrna_illumina_microarray.txt", "metabric_mrna_illumina_microarray.txt")
    df = pd.read_csv(path, sep="\t")
    sym = df["Hugo_Symbol"].where(df["Hugo_Symbol"].notna(), df["Entrez_Gene_Id"].astype(str))
    df = df.drop(columns=["Hugo_Symbol", "Entrez_Gene_Id"]).set_index(sym)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = normalize_index(df).astype("float32")
    df.to_parquet(out)
    return df


def clinical(endpoint: str = "OS") -> pd.DataFrame:
    if endpoint not in ENDPOINTS:
        raise ValueError(f"endpoint must be one of {ENDPOINTS}")
    p = pd.read_csv(fetch(BASE + "data_clinical_patient.txt", "metabric_clinical_patient.txt"),
                    sep="\t", comment="#").set_index("PATIENT_ID")
    s = pd.read_csv(fetch(BASE + "data_clinical_sample.txt", "metabric_clinical_sample.txt"),
                    sep="\t", comment="#").set_index("PATIENT_ID")
    clin = pd.DataFrame(index=p.index)
    if endpoint == "RFS":
        clin["time"] = p["RFS_MONTHS"] / 12
        clin["event"] = p["RFS_STATUS"].astype(str).str.startswith("1").astype(float)
    else:
        clin["time"] = p["OS_MONTHS"] / 12
        if endpoint == "OS":
            clin["event"] = p["OS_STATUS"].astype(str).str.startswith("1").astype(float)
        else:  # disease-specific: deaths from other causes are censored
            vs = p["VITAL_STATUS"]
            clin["event"] = np.where(vs == "Died of Disease", 1.0,
                                     np.where(vs.isin(["Living", "Died of Other Causes"]), 0.0, np.nan))
    clin["age"] = pd.to_numeric(p["AGE_AT_DIAGNOSIS"], errors="coerce")
    clin["stage"] = pd.to_numeric(s["TUMOR_STAGE"].reindex(clin.index), errors="coerce").replace(0, np.nan)
    clin["sex"] = "female"
    clin["er"] = s["ER_STATUS"].reindex(clin.index).str.lower()
    clin["her2"] = s["HER2_STATUS"].reindex(clin.index).str.lower()
    return clin


def load(endpoint: str = "OS") -> Cohort:
    c = Cohort("metabric", expression(), clinical(endpoint), endpoint=endpoint,
               platform="Illumina HT-12 v3 microarray", cancer="BRCA")
    c.notes.append("No genetic ancestry or diagnosis year in the public METABRIC release.")
    return c
