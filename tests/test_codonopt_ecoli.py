import pytest
from Bio.Seq import Seq
from codonopt import codonopt_ecoli

SD_PATTERNS = ["AGGAGG", "GGAGG", "AGGGA", "AGGAG", "GGAG", "AGGA", "GAGG"]
DEFAULT_PATTERNS = ["BsaI_site", "Esp3I_site", "7xC", "7xG"]


def contains_any(seq, patterns):
    return any(p in seq for p in patterns)


def test_ecoli_translation_is_preserved():
    protein = "MATE*GS"
    dna = codonopt_ecoli(protein)

    translated = str(Seq(dna).translate(table=1, to_stop=False))
    assert translated == protein


def test_ecoli_forced_stop_codon():
    protein = "MA*G"
    dna = codonopt_ecoli(protein)

    # position of stop: index 2 → bases 6–8
    assert dna[6:9] == "TAG"


def test_ecoli_avoids_forbidden_patterns():
    protein = "MA" * 20
    dna = codonopt_ecoli(protein)

    # Forbidden sequences should not occur
    assert not contains_any(dna, ["GGTCTC", "GAGACC", "CGTCTC", "GAGACG"])


def test_ecoli_avoids_sd_sequences():
    protein = "SGTQISTIAESEDSQESVDSVTDSQKRREILSRRPSYRKILNDLSSDAPGVPRIEE"*25
    dna = codonopt_ecoli(protein)

    assert not contains_any(dna, SD_PATTERNS)


def test_ecoli_returns_string():
    dna = codonopt_ecoli("MAST")
    assert isinstance(dna, str)


def test_ecoli_length_is_correct():
    protein = "MAGH"
    dna = codonopt_ecoli(protein)
    assert len(dna) == len(protein) * 3


def test_ecoli_custom_avoid_patterns():
    protein = "MADG"
    dna = codonopt_ecoli(
        protein,
        avoid_patterns=["AAAAAA"]  # extremely simple to detect
    )
    assert "AAAAAA" not in dna

