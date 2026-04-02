#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path
from typing import Iterator, TextIO

import pysam


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def iter_data_rows(path: str, delimiter: str) -> Iterator[list[str]]:
    with open_text(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            yield next(csv.reader([line], delimiter=delimiter))


def get_col_index(meta: dict, key: str) -> int | None:
    entry = meta["columns"].get(key)
    if entry and entry.get("found"):
        return entry["index_0based"]
    return None


def fetch_ref_base(fasta: pysam.FastaFile, chrom: str, pos1: int) -> str:
    # pysam fetch usa coordinates 0-based half-open
    return fasta.fetch(chrom, pos1 - 1, pos1).upper()


def normalize_chrom_for_fasta(chrom: str, fasta: pysam.FastaFile) -> str:
    chrom = str(chrom).strip()

    refs = set(fasta.references)

    # match exacto
    if chrom in refs:
        return chrom

    # intentar agregar chr
    if not chrom.startswith("chr") and f"chr{chrom}" in refs:
        return f"chr{chrom}"

    # intentar quitar chr
    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]

    raise ValueError(f"Chromosome '{chrom}' no encontrado en FASTA")


def is_snv(allele: str) -> bool:
    allele = allele.strip().upper()
    return len(allele) == 1 and allele in {"A", "C", "G", "T"}


def process(
    input_file: str,
    meta_file: str,
    fasta_path: str,
    output: str | None = None,
) -> str:
    meta = json.loads(Path(meta_file).read_text(encoding="utf-8"))

    delimiter = "\t" if meta["delimiter"] == "\\t" else meta["delimiter"]

    fasta = pysam.FastaFile(fasta_path)

    idx_chr = get_col_index(meta, "CHROM")
    idx_pos = get_col_index(meta, "POS")
    idx_eff = get_col_index(meta, "effect_allele")
    idx_other = get_col_index(meta, "other_allele")

    if idx_chr is None or idx_pos is None:
        raise ValueError("CHROM/POS no encontrados en metadata de columnas")

    if idx_eff is None or idx_other is None:
        raise ValueError("effect_allele/other_allele no encontrados en metadata de columnas")

    # Intentar effect_weight primero, luego weight
    header = meta.get("header", [])
    norm_header = ["".join(ch for ch in h.strip().lower() if ch.isalnum()) for h in header]

    idx_beta = None
    for candidate in ("effectweight", "weight"):
        if candidate in norm_header:
            idx_beta = norm_header.index(candidate)
            break

    if idx_beta is None:
        raise ValueError("No se encontró columna de peso: effect_weight/weight")

    out_path = output
    if out_path is None:
        p = Path(input_file)
        if p.name.endswith(".txt.gz"):
            out_path = str(p.with_name(p.name[:-7] + ".posmap.beta.tsv"))
        elif p.name.endswith(".tsv.gz"):
            out_path = str(p.with_name(p.name[:-7] + ".posmap.beta.tsv"))
        elif p.name.endswith(".csv.gz"):
            out_path = str(p.with_name(p.name[:-7] + ".posmap.beta.tsv"))
        else:
            out_path = str(p.with_suffix(".posmap.beta.tsv"))

    n_total = 0
    n_written = 0
    n_skipped_non_snv = 0
    n_skipped_ref_mismatch = 0
    n_skipped_parse = 0

    with open(out_path, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "IS_FLIP", "BETA"])

        row_iter = iter_data_rows(input_file, delimiter)

        # Saltar header real
        try:
            next(row_iter)
        except StopIteration:
            raise ValueError("Archivo sin filas de datos")

        for row in row_iter:
            n_total += 1
            try:
                chrom_raw = row[idx_chr].strip()
                pos = int(row[idx_pos])
                effect = row[idx_eff].strip().upper()
                other = row[idx_other].strip().upper()
                beta = float(row[idx_beta])

                if not is_snv(effect) or not is_snv(other):
                    n_skipped_non_snv += 1
                    continue

                chrom = normalize_chrom_for_fasta(chrom_raw, fasta)
                ref = fetch_ref_base(fasta, chrom, pos)

                if effect == ref:
                    alt = other
                    is_flip = 0
                elif other == ref:
                    alt = effect
                    is_flip = 1
                else:
                    n_skipped_ref_mismatch += 1
                    continue

                writer.writerow([chrom, pos, ref, alt, is_flip, beta])
                n_written += 1

            except Exception:
                n_skipped_parse += 1
                continue

    fasta.close()

    stats = {
        "input_file": input_file,
        "meta_file": meta_file,
        "fasta": fasta_path,
        "output_file": out_path,
        "n_total_rows": n_total,
        "n_written": n_written,
        "n_skipped_non_snv": n_skipped_non_snv,
        "n_skipped_ref_mismatch": n_skipped_ref_mismatch,
        "n_skipped_parse": n_skipped_parse,
    }

    stats_path = out_path + ".stats.json"
    Path(stats_path).write_text(json.dumps(stats, indent=2), encoding="utf-8")

    print(json.dumps(stats, indent=2), file=sys.stderr)
    return out_path


def main() -> int:
    p = argparse.ArgumentParser(
        description="Build CHROM/POS/REF/ALT/IS_FLIP/BETA map from PGS score file"
    )
    p.add_argument("input", help="Score file (.tsv/.csv/.gz)")
    p.add_argument("meta", help="JSON de columnas generado por inspect_score_columns.py")
    p.add_argument("--fasta", required=True, help="Referencia FASTA, e.g. hg38.fa")
    p.add_argument("-o", "--output", help="Output TSV")
    args = p.parse_args()

    out = process(args.input, args.meta, args.fasta, args.output)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
