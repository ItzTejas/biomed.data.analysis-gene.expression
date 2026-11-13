"""Cohort loaders. `load("tcga:BRCA")`, `load("metabric")`, `load("geo:GSE2034")`."""
from __future__ import annotations

from ..cohort import Cohort


def load(spec: str, endpoint: str = "OS") -> Cohort:
    kind, _, arg = spec.partition(":")
    kind = kind.lower()
    if kind == "tcga":
        from . import tcga
        return tcga.load(arg, endpoint)
    if kind == "metabric":
        from . import metabric
        return metabric.load(endpoint)
    if kind == "geo":
        from . import geo
        return geo.load(arg, endpoint)
    raise ValueError(f"unknown cohort {spec!r} (try tcga:BRCA, metabric, geo:GSE...)")
