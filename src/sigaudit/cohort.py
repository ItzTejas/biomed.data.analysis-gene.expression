"""A cohort: an expression matrix plus survival and patient attributes."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# columns the audit understands in Cohort.clin
#   time      follow-up in years (required)
#   event     1 = event, 0 = censored (required)
#   age       years at diagnosis
#   stage     1-4 (ordinal)
#   sex       "female" / "male"
#   ancestry  EUR / AFR / EAS / AMR / SAS / ADMIX
#   year      calendar year of diagnosis


@dataclass
class Cohort:
    name: str                  # e.g. "tcga:BRCA"
    expr: pd.DataFrame         # genes (approved symbols) x samples, log scale
    clin: pd.DataFrame         # indexed by the same sample ids
    endpoint: str = "OS"
    platform: str = ""
    cancer: str = ""
    notes: list[str] = field(default_factory=list)

    def __post_init__(self):
        common = self.expr.columns.intersection(self.clin.index)
        clin = self.clin.loc[common]
        clin = clin[np.isfinite(clin["time"]) & clin["event"].notna() & (clin["time"] > 0)]
        self.clin = clin
        self.expr = self.expr[clin.index]

    @property
    def n(self) -> int:
        return self.clin.shape[0]

    @property
    def events(self) -> int:
        return int(self.clin["event"].sum())

    def truncate(self, years: float | None) -> "Cohort":
        """Administrative censoring at `years` so cohorts with long follow-up are comparable."""
        if not years:
            return self
        c = self.clin.copy()
        over = c["time"] > years
        c.loc[over, "event"] = 0
        c.loc[over, "time"] = years
        return Cohort(self.name, self.expr, c, self.endpoint, self.platform, self.cancer, list(self.notes))

    def subset(self, rule: str | None) -> "Cohort":
        """Keep patients matching e.g. "er=positive" or "er=positive;stage=1" (case-insensitive).

        Raises KeyError if the cohort has no such column, so callers can skip it.
        """
        if not rule:
            return self
        clin = self.clin
        for part in rule.split(";"):
            col, _, val = part.partition("=")
            col, val = col.strip(), val.strip().lower()
            if col not in clin:
                raise KeyError(f"{self.name} has no '{col}' annotation")
            clin = clin[clin[col].astype(str).str.lower() == val]
        c = Cohort(self.name, self.expr, clin, self.endpoint, self.platform, self.cancer, list(self.notes))
        c.notes.append(f"restricted to {rule}")
        return c

    def zscores(self) -> pd.DataFrame:
        """Per-gene z-scores across the cohort (the scale signatures are scored on)."""
        x = self.expr
        mu = x.mean(axis=1)
        sd = x.std(axis=1, ddof=1).replace(0, np.nan)
        return x.sub(mu, axis=0).div(sd, axis=0)
