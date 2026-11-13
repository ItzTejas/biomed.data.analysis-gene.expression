"""A gene signature: genes, their weights, and where it came from."""
from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .genes import normalize


@dataclass
class Signature:
    name: str
    genes: list[str]
    weights: list[float]
    source: str = ""          # DOI / PMID / URL
    cancer: str = ""          # TCGA code the authors built it for, if any
    training_cohorts: list[str] = field(default_factory=list)  # e.g. ["tcga:BRCA"]
    notes: str = ""

    def __post_init__(self):
        if len(self.genes) != len(self.weights):
            raise ValueError("genes and weights must have the same length")
        if not self.genes:
            raise ValueError("a signature needs at least one gene")

    @property
    def size(self) -> int:
        return len(self.genes)

    @property
    def weighted(self) -> bool:
        return len({abs(w) for w in self.weights}) > 1

    def normalized(self) -> tuple["Signature", list[str]]:
        """Copy with HGNC-approved symbols; returns (signature, unmapped originals)."""
        genes, weights, unmapped = [], [], []
        merged: dict[str, float] = {}
        for g, w in zip(self.genes, self.weights):
            s = normalize(g)
            if s is None:
                unmapped.append(g)
                continue
            merged[s] = merged.get(s, 0.0) + w
        genes, weights = list(merged), list(merged.values())
        if not genes:
            raise ValueError(f"none of the genes in {self.name!r} map to HGNC symbols")
        return Signature(self.name, genes, weights, self.source, self.cancer,
                         list(self.training_cohorts), self.notes), unmapped

    # ------------------------------------------------------------ constructors

    @classmethod
    def from_genes(cls, genes, name="signature", direction=1.0, **kw) -> "Signature":
        genes = [g for g in (x.strip() for x in genes) if g]
        return cls(name, genes, [float(direction)] * len(genes), **kw)

    @classmethod
    def from_text(cls, text: str, name="signature", **kw) -> "Signature":
        """Parse 'A,B,C', 'A:0.3,B:-1.2', or a published risk-score formula."""
        text = text.strip()
        if re.search(r"[*×·]", text) or "=" in text:
            pairs = parse_formula(text)
            if pairs:
                return cls(name, [g for g, _ in pairs], [w for _, w in pairs], **kw)
        items = [t for t in re.split(r"[,\s;]+", text) if t]
        genes, weights = [], []
        for it in items:
            if ":" in it:
                g, w = it.split(":", 1)
                genes.append(g)
                weights.append(float(w))
            else:
                genes.append(it)
                weights.append(1.0)
        return cls(name, genes, weights, **kw)

    @classmethod
    def load(cls, path: str | Path) -> "Signature":
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        if "formula" in d and "genes" not in d:
            pairs = parse_formula(d.pop("formula"))
            d["genes"] = [g for g, _ in pairs]
            d["weights"] = [w for _, w in pairs]
        d.setdefault("weights", [1.0] * len(d["genes"]))
        return cls(**d)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")


# A coefficient: optional sign, digits, optional exponent (also "×10^-3" style)
_NUM = r"[-+−–]?\s*\(?\s*[-+−–]?\d*\.?\d+(?:\s*(?:[eE]|[×x]\s*10\^?)\s*[-+−–]?\d+)?\s*\)?"
# gene symbols: segments joined by "-" or "." (NKX2-1, AC012073.1), but a
# segment that is a decimal number ("CCL20-0.148") belongs to the next term
_GENE = r"[A-Za-z][A-Za-z0-9]*(?:[-.](?!\d*\.\d)[A-Za-z0-9]+)*"
_TERM_A = re.compile(rf"({_NUM})\s*[*×·]\s*(?:expression\s+(?:of|level\s+of)\s+|Exp\s*)?\(?({_GENE})\)?")
_TERM_B = re.compile(rf"\(?({_GENE})\)?\s*[*×·]\s*({_NUM})")


def _num(s: str) -> float:
    s = s.replace("−", "-").replace("–", "-").replace(" ", "").replace("(", "").replace(")", "")
    s = re.sub(r"[×x]10\^?", "e", s)
    s = s.replace("+-", "-").replace("--", "+").replace("-+", "-")
    return float(s)


_STOP = {"of", "the", "and", "for", "with", "exp", "expression", "risk", "score", "coef",
         "coefficient", "level", "value", "sum", "gene", "genes", "where", "each", "is", "was"}


_PLACEHOLDER = re.compile(r"(?i)^((m|lnc|mi)?rna|gene|lncrna|mirna|g|x|coef|beta|β)\d*[a-z]?$|^z-?scores?$")


def _looks_like_gene(g: str) -> bool:
    return (len(g) >= 2 and g.lower() not in _STOP and any(c.isupper() for c in g)
            and not _PLACEHOLDER.match(g))


def _clean_gene(g: str) -> str:
    g = g.strip(".-")
    # "expressionBOP1", "ExpNOS1", "ExpressionOfTP53" -> the gene itself
    g = re.sub(r"^(?:[Ee]xpression(?:[Oo]f)?|[Ee]xp)(?=[A-Z0-9])", "", g)
    # "RGS5exp", "TP53expression" -> the gene itself
    return re.sub(r"(?<=[A-Z0-9])(?:[Ee]xp|[Ee]xpr|[Ee]xpression)$", "", g)


def _delatex(t: str) -> str:
    """Turn LaTeX math ('\\mathrm{Exp}\\mathrm{GAS}5', '\\left(\\hbox{-} 0.49\\right)') into plain text."""
    if "\\" not in t:
        return t
    t = t.split("\\begin{document}")[-1]
    t = re.sub(r"\\hbox\{\s*-\s*\}", "-", t)
    t = re.sub(r"\\(times|cdot|ast)\b", "×", t)
    t = re.sub(r"\\(left|right|mathrm|mathit|text|textit|textrm|operatorname|mathbf|boldsymbol)\b", "", t)
    t = re.sub(r"\\[a-zA-Z]+", " ", t)
    t = t.replace("{", "").replace("}", "").replace("$$", " ")
    return t


def parse_formula(text: str) -> list[tuple[str, float]]:
    """Pull (gene, coefficient) pairs out of a risk-score formula.

    Handles the common published forms:
      risk score = 0.123*GENE1 + (-0.456)*GENE2
      Risk = (0.12 × ExpGENE1) − (0.34 × ExpGENE2)
      riskScore = GENE1*0.12 + GENE2*(-0.34)
    """
    t = _delatex(text).replace("`", "")
    t = t.replace("−", "-").replace("–", "-").replace("[", "(").replace("]", ")")
    t = t.split("=", 1)[1] if "=" in t else t
    # treat a minus sign before a bracketed term as a sign on the coefficient
    t = re.sub(r"-\s*\(\s*(\d)", r"+(-\1", t)
    # "(- 0.23)" -> "(-0.23)"
    t = re.sub(r"([-+])\s+(?=\d)", r"\1", t)
    pairs: list[tuple[str, float]] = []
    for rx, gi, ni in ((_TERM_A, 2, 1), (_TERM_B, 1, 2)):
        found = []
        for m in rx.finditer(t):
            g = _clean_gene(m.group(gi))
            try:
                w = _num(m.group(ni))
            except ValueError:
                continue
            if not _looks_like_gene(g):
                continue
            found.append((g, w))
        if len(found) > len(pairs):
            pairs = found
    # a template ("coef1 × gene1 + coef2 × gene2 ...") is not a signature
    ws = [w for _, w in pairs]
    if len(ws) >= 2 and ws == [float(k) for k in range(1, len(ws) + 1)]:
        return []
    return pairs
