"""Gene symbol normalisation against HGNC.

Published signatures use whatever symbols were current when the paper was
written (and sometimes Excel-mangled ones like "1-Mar"). Expression datasets
use their own vintage. Mapping both to today's approved HGNC symbols is what
lets a signature be scored on a cohort at all.
"""
from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from .cache import fetch, processed

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

# Excel turns SEPT1 -> 1-Sep, MARCH1 -> 1-Mar, DEC1 -> 1-Dec
_EXCEL = {"sep": "SEPTIN", "mar": "MARCHF", "dec": "DELEC"}
_EXCEL_RE = re.compile(r"^(\d{1,2})-(sep|mar|dec)$", re.I)


@lru_cache(maxsize=1)
def symbol_map() -> dict[str, str]:
    """Upper-case symbol/alias/previous symbol/Ensembl id -> approved symbol."""
    cache = processed("hgnc_map.parquet")
    if cache.exists():
        df = pd.read_parquet(cache)
        return dict(zip(df["key"], df["symbol"]))
    path = fetch(HGNC_URL, "hgnc_complete_set.txt")
    h = pd.read_csv(path, sep="\t", dtype=str, usecols=lambda c: c in {
        "symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id", "entrez_id", "status"})
    h = h[h["status"].fillna("Approved") == "Approved"]
    approved = {s.upper(): s for s in h["symbol"]}
    out: dict[str, str] = {}
    # previous symbols outrank aliases; ambiguous keys (pointing to >1 gene) are dropped
    for col in ("prev_symbol", "alias_symbol"):
        seen: dict[str, set] = {}
        for sym, field in zip(h["symbol"], h[col].fillna("")):
            for a in field.split("|"):
                a = a.strip().upper()
                if a:
                    seen.setdefault(a, set()).add(sym)
        for a, syms in seen.items():
            if a not in approved and a not in out and len(syms) == 1:
                out[a] = next(iter(syms))
    for col in ("ensembl_gene_id", "entrez_id"):
        for sym, v in zip(h["symbol"], h[col].fillna("")):
            if v:
                out[v.upper()] = sym
    out.update(approved)
    pd.DataFrame({"key": list(out), "symbol": list(out.values())}).to_parquet(cache)
    return out


def normalize(symbol: str) -> str | None:
    """Approved HGNC symbol for any symbol/alias/Ensembl/Entrez id, else None."""
    if symbol is None:
        return None
    s = str(symbol).strip().strip("'\"").split(".")[0] if str(symbol).upper().startswith("ENSG") \
        else str(symbol).strip().strip("'\"")
    if not s:
        return None
    m = _EXCEL_RE.match(s)
    if m:
        s = _EXCEL[m.group(2).lower()] + m.group(1)
    return symbol_map().get(s.upper())


def normalize_index(df: pd.DataFrame) -> pd.DataFrame:
    """Re-key a genes x samples matrix by approved symbols, averaging duplicates."""
    mapped = [normalize(g) for g in df.index]
    df = df.copy()
    df.index = mapped
    df = df[df.index.notna()]
    if df.index.has_duplicates:
        df = df.groupby(level=0).mean()
    return df
