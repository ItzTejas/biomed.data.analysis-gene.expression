import pytest

from sigaudit.signature import Signature, parse_formula


@pytest.mark.parametrize("text,expected", [
    ("Risk score = 0.123*CDK1 + (-0.456)*TP53 + 0.2*MKI67",
     [("CDK1", 0.123), ("TP53", -0.456), ("MKI67", 0.2)]),
    ("riskScore = (0.12 × ExpCDK1) − (0.34 × ExpTP53)",
     [("CDK1", 0.12), ("TP53", -0.34)]),
    ("risk = CDK1*0.12 + TP53*(-0.34)",
     [("CDK1", 0.12), ("TP53", -0.34)]),
    ("risk score = 0.0815 * AC012073.1 - 0.2136 * LINC01535",
     [("AC012073.1", 0.0815), ("LINC01535", -0.2136)]),
    ("Risk score = 0.104×CCL20−0.148×CD69 + 0.2×NKX2-1",
     [("CCL20", 0.104), ("CD69", -0.148), ("NKX2-1", 0.2)]),
    ("risk score = 0.2477 × expression of TP53 + 0.1 × MYC",
     [("TP53", 0.2477), ("MYC", 0.1)]),
    ("Risk score = (3.2 × expressionBOP1) − (1.9 × expressionSERPINE1)",
     [("BOP1", 3.2), ("SERPINE1", -1.9)]),
    ("index = [FURIN* 0.1089] + [PSMD14*(− 0.2344)] + (TUBB3*0.1028)",
     [("FURIN", 0.1089), ("PSMD14", -0.2344), ("TUBB3", 0.1028)]),
    ("risk score = βmRNA1×ExpressionmRNA1+βmRNA2×ExpressionmRNA2", []),
    ("Risk Score = (a1 * TNFRSF14) + (a2 * SUSD3) + (a3 * COX7A1)", []),
    (r"$$ \mathrm{RS}=(0.11)\times \mathrm{Exp}\mathrm{GAS}5+\left(\hbox{-} 0.49\right)\times \mathrm{Exp}\mathrm{HCP}5$$",
     [("GAS5", 0.11), ("HCP5", -0.49)]),
    ("Risk score=LINC01767×(0.2877)+`UBE2Q1-AS1`×(-0.6927)+`MAFG-DT`×(0.4587)",
     [("LINC01767", 0.2877), ("UBE2Q1-AS1", -0.6927), ("MAFG-DT", 0.4587)]),
    ("(-0.551 * PPARGC1A Exp) + (1.258*PPP1R17Exp)",
     [("PPARGC1A", -0.551), ("PPP1R17", 1.258)]),
    ("risk score = 1.2e-2*CDK1 + 3.5×10^-3*TP53",
     [("CDK1", 0.012), ("TP53", 0.0035)]),
])
def test_parse_formula(text, expected):
    got = parse_formula(text)
    assert [g for g, _ in got] == [g for g, _ in expected]
    for (_, a), (_, b) in zip(got, expected):
        assert a == pytest.approx(b)


def test_from_text_variants():
    s = Signature.from_text("A,B , C")
    assert s.genes == ["A", "B", "C"] and s.weights == [1, 1, 1]
    s = Signature.from_text("A:0.5,B:-2")
    assert s.weights == [0.5, -2.0] and s.weighted


def test_normalize_merges_aliases():
    s = Signature.from_text("CDCA1:1,NUF2:1,notagene:1")
    n, unmapped = s.normalized()
    assert n.genes == ["NUF2"] and n.weights == [2.0]
    assert unmapped == ["notagene"]
