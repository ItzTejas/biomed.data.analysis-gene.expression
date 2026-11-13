"""Command line: `sigaudit run ...` and `sigaudit fetch ...`."""
from __future__ import annotations

import argparse
import sys

from .audit import AuditSettings, run_audit
from .loaders import load
from .report import to_markdown, write
from .signature import Signature


def _cohort_specs(items: list[str]) -> list[tuple[str, str]]:
    """'tcga:BRCA@PFI' -> ('tcga:BRCA', 'PFI')."""
    out = []
    for it in items:
        for part in it.split(","):
            spec, _, ep = part.strip().partition("@")
            out.append((spec, ep or "OS"))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="sigaudit", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="audit a signature")
    src = r.add_mutually_exclusive_group(required=True)
    src.add_argument("--genes", help="'A,B,C' or 'A:0.3,B:-1.2' or a risk-score formula in quotes")
    src.add_argument("--file", help="signature JSON (see examples/)")
    r.add_argument("--name", default="signature")
    r.add_argument("--source", default="", help="DOI or URL of the paper")
    r.add_argument("--training", nargs="*", default=[],
                   help="cohorts the signature was built on, e.g. tcga:BRCA (excluded from the verdict)")
    r.add_argument("--cohort", nargs="+", required=True,
                   help="cohorts to test, with optional endpoint: tcga:BRCA@PFI metabric@DSS")
    r.add_argument("--subset", default=None,
                   help='restrict every cohort to matching patients, e.g. "er=positive"')
    r.add_argument("--random", type=int, default=1000, help="random signatures in the null")
    r.add_argument("--boot", type=int, default=200, help="bootstrap resamples for intervals")
    r.add_argument("--horizon", type=float, default=10.0, help="censor follow-up at N years (0 = off)")
    r.add_argument("--out", default="out", help="output folder for the audit card")
    r.add_argument("--seed", type=int, default=0)

    f = sub.add_parser("fetch", help="download and cache cohorts")
    f.add_argument("cohort", nargs="+")

    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):  # Windows consoles default to cp1252
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    if a.cmd == "fetch":
        for spec, ep in _cohort_specs(a.cohort):
            c = load(spec, ep)
            print(f"{c.name}: {c.n} patients, {c.events} events ({ep}), {c.expr.shape[0]} genes")
        return 0

    if a.file:
        sig = Signature.load(a.file)
    else:
        sig = Signature.from_text(a.genes, a.name, source=a.source, training_cohorts=a.training)
    if a.training and a.file:
        sig.training_cohorts = a.training
    cohorts = []
    for spec, ep in _cohort_specs(a.cohort):
        c = load(spec, ep)
        try:
            cohorts.append(c.subset(a.subset))
        except KeyError as err:
            print(f"[sigaudit] skipping {spec}: {err}", file=sys.stderr)
    settings = AuditSettings(n_random=a.random, n_boot=a.boot, horizon_years=a.horizon or None, seed=a.seed)
    res = run_audit(sig, cohorts, settings)
    paths = write(res, a.out)
    print(to_markdown(res))
    print(f"\nWrote {', '.join(str(p) for p in paths.values())}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
