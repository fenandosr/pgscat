import pytest

from pgscat.common_utils import query_fasta


class FakeFasta:
    def __init__(self, seq_by_chrom: dict[str, str]):
        self.seq_by_chrom = seq_by_chrom
        self.references = list(seq_by_chrom.keys())

    def fetch(self, chrom: str, start: int, end: int) -> str:
        seq = self.seq_by_chrom[chrom]
        return seq[start:end]


def test_query_fasta_snv_effect_is_alt():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})
    # posiciones 1-based:
    # chr1:1 A
    # chr1:2 A
    # chr1:3 C
    # chr1:4 G
    # chr1:5 T

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="A",
        other_allele="G",
        flank_bp=2,
    )

    assert result["ok"] is True
    assert result["chrom"] == "chr1"
    assert result["pos"] == 4
    assert result["ref"] == "G"
    assert result["alt"] == "A"
    assert result["is_flip"] == 0
    assert result["reason"] is None
    assert result["ref_flank"] == "ACGTA"


def test_query_fasta_snv_effect_is_ref():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="G",
        other_allele="A",
        flank_bp=2,
    )

    assert result["ok"] is True
    assert result["ref"] == "G"
    assert result["alt"] == "A"
    assert result["is_flip"] == 1
    assert result["reason"] is None


def test_query_fasta_normalizes_chrom():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})

    result = query_fasta(
        fasta=fasta,
        chrom="1",
        pos=3,
        effect_allele="T",
        other_allele="C",
        flank_bp=1,
    )

    assert result["ok"] is True
    assert result["chrom"] == "chr1"
    assert result["ref"] == "C"
    assert result["alt"] == "T"
    assert result["is_flip"] == 0


def test_query_fasta_missing_other_allele():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="A",
        other_allele=None,
        flank_bp=2,
    )

    assert result["ok"] is False
    assert result["reason"] == "missing_allele"
    assert result["ref"] is None
    assert result["alt"] is None
    assert result["is_flip"] is None


def test_query_fasta_no_ref_match():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="T",
        other_allele="A",
        flank_bp=2,
    )

    assert result["ok"] is False
    assert result["ref"] == "G"
    assert result["alt"] is None
    assert result["is_flip"] is None
    assert result["reason"].startswith("no_ref_match")


def test_query_fasta_insertion_like_case_effect_is_alt():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})
    # pos=4 -> secuencia desde 4: G...
    # max_len = 2, ref = "GT"

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="AA",
        other_allele="GT",
        flank_bp=2,
    )

    assert result["ok"] is True
    assert result["ref"] == "GT"
    assert result["alt"] == "AA"
    assert result["is_flip"] == 0


def test_query_fasta_insertion_like_case_effect_is_ref():
    fasta = FakeFasta({"chr1": "AACGTACGTAAA"})
    # max_len = 2, ref = "GT"

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="GT",
        other_allele="AA",
        flank_bp=2,
    )

    assert result["ok"] is True
    assert result["ref"] == "GT"
    assert result["alt"] == "AA"
    assert result["is_flip"] == 1


def test_query_fasta_reference_out_of_range():
    fasta = FakeFasta({"chr1": "ACGT"})

    result = query_fasta(
        fasta=fasta,
        chrom="chr1",
        pos=4,
        effect_allele="TT",
        other_allele="AA",
        flank_bp=2,
    )

    assert result["ok"] is False
    assert result["reason"] == "reference_out_of_range"
