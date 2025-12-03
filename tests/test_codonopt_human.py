import pytest
from Bio.Seq import Seq
from codonopt import codonopt_human

DEFAULT_PATTERNS = ["BsaI_site", "Esp3I_site", "7xC", "7xG"]

def contains_any(seq, patterns):
    return any(p in seq for p in patterns)


def test_human_translation_is_preserved():
    protein = "MPLATE*GG"
    dna = codonopt_human(protein)

    translated = str(Seq(dna).translate(table=1, to_stop=False))
    assert translated == protein


def test_human_avoids_forbidden_patterns():
    protein = "MA" * 20
    dna = codonopt_human(protein)

    assert not contains_any(dna, ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"])


def test_human_returns_string():
    dna = codonopt_human("MSTN")
    assert isinstance(dna, str)


def test_human_length_is_correct():
    protein = "MAST"
    dna = codonopt_human(protein)
    assert len(dna) == len(protein) * 3


def test_human_custom_avoid_patterns():
    protein = "MADG"
    dna = codonopt_human(
        protein,
        avoid_patterns=["AAAAAA"]
    )
    assert "AAAAAA" not in dna

