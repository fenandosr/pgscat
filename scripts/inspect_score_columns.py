#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path
from typing import Iterable, TextIO


COMMENT_PREFIXES = ("#", "//")


COLUMN_SYNONYMS = {
    "CHROM": [
        "chrom",
        "chr",
        "chromosome",
        "chr_name",
        "hm_chr",
        "contig",
    ],
    "POS": [
        "pos",
        "position",
        "bp",
        "chr_position",
        "hm_pos",
        "genomic_position",
    ],
    "ID": [
        "id",
        "snp",
        "snpid",
        "markername",
        "marker_name",
        "variant_id",
        "varid",
        "rsid",
        "rs_id",
        "rs number",
        "hm_rsid",
    ],
    "REF": [
        "ref",
        "reference",
        "reference_allele",
        "ref_allele",
        "other_allele",
        "non_effect_allele",
        "hm_inferotherallele",
        "hm_infer_other_allele",
    ],
    "ALT": [
        "alt",
        "alternate",
        "alternate_allele",
        "alt_allele",
        "effect_allele",
        "risk_allele",
        "hm_effect_allele",
    ],
    "other_allele": [
        "other_allele",
        "non_effect_allele",
        "reference_allele",
        "ref_allele",
        "hm_inferotherallele",
        "hm_infer_other_allele",
    ],
    "effect_allele": [
        "effect_allele",
        "risk_allele",
        "alt_allele",
        "alternate_allele",
        "hm_effect_allele",
    ],
}


def normalize_name(name: str) -> str:
    return "".join(ch for ch in name.strip().lower() if ch.isalnum())


def open_text_auto(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", newline="")
    return open(p, "rt", encoding="utf-8", newline="")


def iter_noncomment_lines(handle: TextIO) -> Iterable[str]:
    for line in handle:
        stripped = line.strip()
        if not stripped:
            continue
        if any(stripped.startswith(prefix) for prefix in COMMENT_PREFIXES):
            continue
        yield line


def detect_delimiter(sample: str, forced: str | None) -> str:
    if forced:
        if forced == "tab":
            return "\t"
        if forced == "comma":
            return ","
        raise ValueError(f"Unsupported delimiter mode: {forced}")

    tab_count = sample.count("\t")
    comma_count = sample.count(",")

    if tab_count > 0 and tab_count >= comma_count:
        return "\t"
    if comma_count > 0:
        return ","

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        return dialect.delimiter
    except csv.Error:
        return "\t"


def find_header_and_delimiter(path: str, forced_delimiter: str | None) -> tuple[list[str], str]:
    with open_text_auto(path) as handle:
        lines = iter_noncomment_lines(handle)
        try:
            first = next(lines)
        except StopIteration:
            raise ValueError("No non-comment, non-empty lines found")

        delim = detect_delimiter(first, forced_delimiter)
        header = next(csv.reader([first], delimiter=delim))
        return header, delim


def build_lookup(header: list[str]) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for idx, col in enumerate(header):
        norm = normalize_name(col)
        if norm and norm not in lookup:
            lookup[norm] = idx
    return lookup


def find_column(lookup: dict[str, int], candidates: list[str]) -> tuple[int | None, str | None]:
    for candidate in candidates:
        key = normalize_name(candidate)
        if key in lookup:
            return lookup[key], candidate
    return None, None


def inspect_file(path: str, forced_delimiter: str | None = None) -> dict:
    header, delim = find_header_and_delimiter(path, forced_delimiter)
    lookup = build_lookup(header)

    columns = {}
    missing = []

    for logical_name, synonyms in COLUMN_SYNONYMS.items():
        idx, matched_synonym = find_column(lookup, synonyms)
        if idx is None:
            columns[logical_name] = {
                "found": False,
                "name_in_file": None,
                "matched_on": None,
                "index_0based": None,
                "index_1based": None,
            }
            missing.append(logical_name)
        else:
            columns[logical_name] = {
                "found": True,
                "name_in_file": header[idx],
                "matched_on": matched_synonym,
                "index_0based": idx,
                "index_1based": idx + 1,
            }

    return {
        "file": path,
        "delimiter": "\\t" if delim == "\t" else delim,
        "n_columns": len(header),
        "header": header,
        "columns": columns,
        "missing": missing,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect a CSV/TSV(.gz) score file, skipping comments, and identify "
            "key columns in a machine-readable JSON output."
        )
    )
    parser.add_argument("input", help="Input file (.csv/.tsv optionally .gz), or - for stdin")
    parser.add_argument(
        "--delimiter",
        choices=["tab", "comma"],
        default=None,
        help="Force delimiter instead of auto-detect",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output JSON file, or - for stdout",
    )

    args = parser.parse_args()

    try:
        result = inspect_file(args.input, forced_delimiter=args.delimiter)
    except Exception as e:
        err = {
            "file": args.input,
            "error": str(e),
        }
        payload = json.dumps(err, ensure_ascii=False, indent=2 if args.pretty else None)
        if args.output == "-":
            print(payload, file=sys.stdout)
        else:
            Path(args.output).write_text(payload + "\n", encoding="utf-8")
        return 1

    payload = json.dumps(result, ensure_ascii=False, indent=2 if args.pretty else None)

    if args.output == "-":
        sys.stdout.write(payload)
        sys.stdout.write("\n")
    else:
        Path(args.output).write_text(payload + "\n", encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
