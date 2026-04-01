#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Iterator, TextIO

import pysam


# ─────────────────────────────────────────────
# IO helpers
# ─────────────────────────────────────────────

def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def iter_rows(path: str, delimiter: str) -> Iterator[list[str]]:
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            yield next(csv.reader([line], delimiter=delimiter))


# ─────────────────────────────────────────────
# Core logic
# ─────────────────────────────────────────────

def get_col_index(meta: dict, key: str) -> int | None:
    entry = meta["columns"].get(key)
    if entry and entry["found"]:
        return entry["index_0based"]
    return None


def fetch_ref(fasta: pysam.FastaFile, chrom: str, pos: int) -> str:
    # pysam uses 0-based, half-open
    return fasta.fetch(chrom, pos - 1, pos).upper()


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def process(
    input_file: str,
    meta_file: str,
    fasta_path: str,
    output: str | None = None,
):
    meta = json.loads(Path(meta_file).read_text())

    delimiter = "\t" if meta["delimiter"] == "\\t" else meta["delimiter"]

    fasta = pysam.FastaFile(fasta_path)

    # column indices
    idx_chr = get_col_index(meta, "CHROM")
    idx_pos = get_col_index(meta, "POS")
    idx_eff = get_col_index(meta, "effect_allele")
    idx_other = get_col_index(meta, "other_allele")

    if idx_chr is None or idx_pos is None:
        raise ValueError("CHROM/POS no encontrados")

    if idx_eff is None or idx_other is None:
        raise ValueError("effect_allele/other_allele no encontrados")

    out_path = output or str(Path(input_file).with_suffix(".posmap.tsv"))

    with open(out_path, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "IS_FLIP"])

        for row in iter_rows(input_file, delimiter):
            try:
                chrom = normalize_chrom(row[idx_chr])
                pos = int(row[idx_pos])

                effect = row[idx_eff].upper()
                other = row[idx_other].upper()

                ref = fetch_ref(fasta, chrom, pos)

                if effect == ref:
                    alt = other
                    is_flip = 0
                elif other == ref:
                    alt = effect
                    is_flip = 1
                else:
                    # ambiguous / mismatch
                    continue

                writer.writerow([chrom, pos, ref, alt, is_flip])

            except Exception:
                continue

    fasta.close()
    return out_path


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Build CHROM/POS/REF/ALT/IS_FLIP map from PGS score file"
    )
    p.add_argument("input", help="Score file (.tsv/.csv/.gz)")
    p.add_argument("meta", help="JSON de columnas (inspect_score_columns.py)")
    p.add_argument("--fasta", required=True, help="Referencia (e.g. hg38.fa)")
    p.add_argument("-o", "--output", help="Output TSV")

    args = p.parse_args()

    out = process(args.input, args.meta, args.fasta, args.output)
    print(out)


if __name__ == "__main__":
    main()
