"""Audit cards: JSON (machine), Markdown (README/GitHub), HTML (one self-contained page)."""
from __future__ import annotations

import base64
import html
import io
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .audit import AuditResult, CohortAudit  # noqa: E402

GRADE_COLOR = {"holds up": "#1a7f37", "mixed": "#9a6700", "does not hold up": "#cf222e",
               "untested": "#6e7781"}


def to_dict(res: AuditResult) -> dict:
    cohorts = []
    for c in res.cohorts:
        d = asdict(c)
        d["null_c_summary"] = (
            {"median": float(np.median(c.null_c)), "p95_abs": float(np.percentile(np.abs(np.array(c.null_c) - .5), 95))}
            if c.null_c else {})
        d.pop("null_c")
        d["beats_random"] = c.beats_random
        d["beats_proliferation"] = c.beats_proliferation
        d["coverage"] = c.coverage
        cohorts.append(d)
    return {"signature": asdict(res.signature), "unmapped_genes": res.unmapped,
            "settings": asdict(res.settings), "verdict": res.verdict(),
            "meta_analysis": res.meta, "cohorts": cohorts}


def _fmt_hr(h: dict) -> str:
    if not h:
        return "–"
    return f"{h['hr']:.2f} ({h['lo']:.2f}–{h['hi']:.2f}), p={h['p']:.2g}"


def to_markdown(res: AuditResult) -> str:
    v = res.verdict()
    s = res.signature
    out = [f"# Audit card: {s.name}", ""]
    if s.source:
        out.append(f"Source: {s.source}  ")
    out.append(f"{s.size} genes{' (weighted)' if s.weighted else ''}"
               + (f"; built on {', '.join(s.training_cohorts)}" if s.training_cohorts else ""))
    out += ["", f"**Verdict: {v['grade']}**", ""] + [f"- {line}" for line in v["lines"]] + [""]
    out += ["| Cohort | Endpoint | n / events | Genes | C-index (95% CI) | vs random | Prolif. C | "
            "HR per SD | HR adj. prolif. | HR adj. clinical |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for c in res.cohorts:
        tag = " (training)" if c.training else ""
        if c.skipped:
            out.append(f"| {c.cohort}{tag} | {c.endpoint} | {c.n} / {c.events} | "
                       f"{len(c.genes_found)}/{len(c.genes_found) + len(c.genes_missing)} | "
                       f"skipped: {c.skipped} | | | | | |")
            continue
        out.append(
            f"| {c.cohort}{tag} | {c.endpoint} | {c.n} / {c.events} | "
            f"{len(c.genes_found)}/{len(c.genes_found) + len(c.genes_missing)} | "
            f"{c.c_index:.3f} ({c.c_lo:.3f}–{c.c_hi:.3f}) | "
            f"{c.random_percentile:.0f}th pct, p={c.random_p:.3f} | {c.prolif_c:.3f} | "
            f"{_fmt_hr(c.hr)} | {_fmt_hr(c.hr_adj_prolif)} | "
            f"{_fmt_hr(c.hr_adj_clinical)}{' (' + ', '.join(c.adj_covariates) + ')' if c.adj_covariates else ''} |")
    sub = [r for c in res.cohorts for r in c.subgroups]
    if sub:
        out += ["", "## Who it works for", "",
                "| Cohort | Attribute | Group | n / events | C-index | Gap vs reference (95% CI) |",
                "|---|---|---|---|---|---|"]
        for r in sub:
            gap = (f"{r['gap']:+.3f} ({r['gap_lo']:+.3f} to {r['gap_hi']:+.3f})"
                   + (" ⚠" if r.get("flag") else "")) if "gap_lo" in r else r.get("note", "")
            c = f"{r['c']:.3f}" if np.isfinite(r["c"]) else "–"
            out.append(f"| {r['cohort']} | {r['attribute']} | {r['group']} | {r['n']} / {r['events']} | "
                       f"{c} | {gap} |")
    missing = sorted({g for c in res.cohorts for g in c.genes_missing})
    if res.unmapped or missing:
        out += ["", "## Gene matching", ""]
        if res.unmapped:
            out.append(f"Not recognised as HGNC symbols: {', '.join(res.unmapped)}  ")
        if missing:
            out.append(f"Not measured in at least one cohort: {', '.join(missing)}")
    out += ["", "## How to read this", "",
            "- **C-index**: probability the score ranks two patients in the right order (0.5 = coin flip).",
            f"- **vs random**: where the signature's |C − 0.5| falls among {res.settings.n_random} random "
            "gene sets of the same size and weights (Venet et al. 2011).",
            "- **HR adj. prolif.**: hazard ratio per SD of the score after adjusting for a meta-PCNA "
            "proliferation score. If this is ~1, the signature mostly measures how fast tumours divide.",
            "- Subgroup gaps are flagged ⚠ only when the whole 95% bootstrap interval is below zero.",
            f"- Follow-up censored at {res.settings.horizon_years} years.", ""]
    return "\n".join(out)


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _null_plot(c: CohortAudit) -> str:
    fig, ax = plt.subplots(figsize=(5.2, 2.6))
    ax.hist(c.null_c, bins=40, color="#8c959f", alpha=.8, label="random gene sets")
    ax.axvline(c.c_index, color="#0969da", lw=2.2, label=f"signature ({c.c_index:.3f})")
    ax.axvline(c.prolif_c, color="#cf222e", lw=1.6, ls="--", label=f"proliferation ({c.prolif_c:.3f})")
    ax.axvline(.5, color="#57606a", lw=.8, ls=":")
    ax.set_xlabel("C-index")
    ax.set_yticks([])
    ax.set_title(c.cohort, fontsize=10, loc="left")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    return _png(fig)


def _forest(res: AuditResult) -> str | None:
    rows = [c for c in res.cohorts if c.hr and not c.skipped]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(5.2, 0.5 + 0.45 * (len(rows) + 1)))
    labels = []
    for k, c in enumerate(rows):
        for off, h, col in ((-.12, c.hr, "#0969da"), (.12, c.hr_adj_prolif, "#cf222e")):
            if h:
                ax.plot([h["lo"], h["hi"]], [k + off] * 2, color=col, lw=1.5)
                ax.plot(h["hr"], k + off, "o", color=col, ms=5)
        labels.append(c.cohort + (" (training)" if c.training else ""))
    if res.meta.get("k", 0) >= 2:
        k = len(rows)
        m = res.meta
        ax.plot([m["lo"], m["hi"]], [k] * 2, color="#24292f", lw=2)
        ax.plot(m["hr"], k, "D", color="#24292f", ms=6)
        labels.append("pooled (independent)")
    ax.axvline(1, color="#57606a", lw=.8, ls=":")
    ax.set_xscale("log")
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Hazard ratio per SD   (blue: unadjusted · red: adjusted for proliferation)", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return _png(fig)


def _subgroup_plot(res: AuditResult) -> str | None:
    rows = [r for c in res.cohorts for r in c.subgroups if np.isfinite(r["c"])]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(5.2, 0.4 + 0.3 * len(rows)))
    for k, r in enumerate(rows):
        col = "#cf222e" if r.get("flag") else ("#0969da" if r["events"] >= 10 else "#8c959f")
        ax.plot(r["c"], k, "o", color=col, ms=5)
        if "gap_lo" in r:
            ax.plot([r["ref_c"] + r["gap_lo"], r["ref_c"] + r["gap_hi"]], [k] * 2, color=col, lw=1.2)
    ax.axvline(.5, color="#57606a", lw=.8, ls=":")
    ax.set_yticks(range(len(rows)),
                  [f"{r['cohort']} · {r['attribute']}={r['group']} (ev={r['events']})" for r in rows], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("C-index within group (grey: <10 events, red: reliably worse)", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return _png(fig)


def write(res: AuditResult, out_dir: str | Path, stem: str | None = None) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = stem or "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in res.signature.name)
    paths = {"json": out / f"{stem}.json", "md": out / f"{stem}.md", "html": out / f"{stem}.html"}
    paths["json"].write_text(json.dumps(to_dict(res), indent=2, default=float), encoding="utf-8")
    paths["md"].write_text(to_markdown(res), encoding="utf-8")
    paths["html"].write_text(render_html(res), encoding="utf-8")
    return paths


def render_html(res: AuditResult) -> str:
    v = res.verdict()
    s = res.signature
    e = html.escape
    imgs = []
    forest = _forest(res)
    if forest:
        imgs.append(("Effect size by cohort", forest))
    for c in res.cohorts:
        if c.null_c:
            imgs.append((f"Against random genes: {c.cohort}", _null_plot(c)))
    sub = _subgroup_plot(res)
    if sub:
        imgs.append(("Who it works for", sub))
    md = to_markdown(res)
    table_rows = []
    for c in res.cohorts:
        if c.skipped:
            table_rows.append(f"<tr><td>{e(c.cohort)}</td><td colspan=6>skipped: {e(c.skipped)}</td></tr>")
            continue
        def mark(ok):
            return '<span class="ok">yes</span>' if ok else '<span class="no">no</span>'
        table_rows.append(
            f"<tr><td>{e(c.cohort)}{' <small>(training)</small>' if c.training else ''}</td>"
            f"<td>{c.n} / {c.events}</td><td>{len(c.genes_found)}/{len(c.genes_found) + len(c.genes_missing)}</td>"
            f"<td>{c.c_index:.3f} <small>({c.c_lo:.3f}–{c.c_hi:.3f})</small></td>"
            f"<td>{mark(c.beats_random)} <small>p={c.random_p:.3f}</small></td>"
            f"<td>{mark(c.beats_proliferation)} <small>{e(_fmt_hr(c.hr_adj_prolif))}</small></td></tr>")
    figs = "".join(f"<figure><figcaption>{e(t)}</figcaption><img alt='{e(t)}' src='data:image/png;base64,{b}'></figure>"
                   for t, b in imgs)
    lines = "".join(f"<li>{e(x)}</li>" for x in v["lines"])
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Audit: {e(s.name)}</title>
<style>
:root{{--bg:#fff;--fg:#1f2328;--muted:#59636e;--line:#d1d9e0;--card:#f6f8fa}}
@media (prefers-color-scheme:dark){{:root{{--bg:#0d1117;--fg:#e6edf3;--muted:#9198a1;--line:#3d444d;--card:#151b23}}
img{{background:#fff;border-radius:6px}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;margin:0;padding:24px 16px}}
main{{max-width:860px;margin:auto}} h1{{font-size:1.5rem;margin:0 0 4px}} .muted{{color:var(--muted)}}
.grade{{display:inline-block;padding:4px 12px;border-radius:999px;color:#fff;font-weight:600;background:{GRADE_COLOR[v['grade']]}}}
table{{border-collapse:collapse;width:100%;font-size:.88rem;margin:12px 0;display:block;overflow-x:auto}}
td,th{{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;white-space:nowrap}}
.ok{{color:#1a7f37;font-weight:600}} .no{{color:#cf222e;font-weight:600}}
figure{{margin:18px 0;background:var(--card);padding:12px;border-radius:8px}} figure img{{max-width:100%}}
figcaption{{font-weight:600;margin-bottom:6px}} details{{margin-top:24px}} pre{{white-space:pre-wrap;font-size:.8rem}}
</style></head><body><main>
<h1>{e(s.name)}</h1>
<div class="muted">{s.size} genes{' · ' + e(s.source) if s.source else ''}{' · built on ' + e(', '.join(s.training_cohorts)) if s.training_cohorts else ''}</div>
<p><span class="grade">{e(v['grade'])}</span></p><ul>{lines}</ul>
<table><thead><tr><th>Cohort</th><th>n / events</th><th>Genes</th><th>C-index</th><th>Beats random</th><th>Survives proliferation</th></tr></thead>
<tbody>{''.join(table_rows)}</tbody></table>
{figs}
<details><summary>Full audit card (Markdown)</summary><pre>{e(md)}</pre></details>
<p class="muted"><small>Generated by sigaudit {__import__('sigaudit').__version__}. A signature that fails these checks
is not evidence of error or misconduct by its authors; it means the prognostic claim did not
transfer to these cohorts under these settings.</small></p>
</main></body></html>"""
