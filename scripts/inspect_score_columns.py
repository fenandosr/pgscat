#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
import logging

from pgscat.common_utils import (
    file_delim_meta_header_data,
    build_lookup,
    find_column,
)

logger = logging.getLogger(__name__)

COLUMN_SYNONYMS = {
    "CHROM": [
        "chrom",
        "chr",
        "chromosome",
        "chr_name",
        "contig",
    ],
    "POS": [
        "pos",
        "position",
        "bp",
        "chr_position",
        "genomic_position",
    ],
    "ID": [
        "hm_rsid",
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
    ],
    "REF": [
        "ref",
        "reference",
        "reference_allele",
        "ref_allele",
    ],
    "ALT": [
        "alt",
        "alternate",
        "alternate_allele",
        "alt_allele",
    ],
    "OTHER_ALLELE": [
        "other_allele",
        "non_effect_allele",
        "hm_inferotherallele",
        "hm_infer_other_allele",
    ],
    "EFFECT_ALLELE": [
        "effect_allele",
        "risk_allele",
    ],
    "EFFECT_WEIGHT": [
        "effect_weight",
        "beta",
        "weight",
        "effectsize",
        "effect_size",
    ],
}


def inspect_file(path: str, forced_delimiter: str | None = None) -> dict:
    delim, meta, header, data = file_delim_meta_header_data(
        path,
        forced_delimiter=forced_delimiter,
        meta_line_startswith="#",
        skip_meta_line_startswith="##",
        meta_kv_separator="=",
        skip_data=True,
    )
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
            }
            missing.append(logical_name)
        else:
            columns[logical_name] = {
                "found": True,
                "name_in_file": header[idx],
                "matched_on": matched_synonym,
                "index_0based": idx,
            }

    return {
        "file": path,
        "delimiter": "\\t" if delim == "\t" else delim,
        "meta": meta,
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
    parser.add_argument(
        "input", help="Input file (.csv/.tsv optionally .gz), or - for stdin"
    )
    parser.add_argument(
        "--delimiter",
        choices=["tab", "comma"],
        default=None,
        help="Force delimiter instead of auto-detect",
    )
    parser.add_argument(
        "--genome_build",
        choices=["GRCh37", "GRCh38"],
        default="GRCh38",
        help="Genome build.",
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

    logger.info(
        f"Inspecting file: {args.input} with forced delimiter: {args.delimiter} and genome build: {args.genome_build}"
    )

    if args.genome_build == "GRCh38":
        COLUMN_SYNONYMS["CHROM"] = ["hm_chr"]
        COLUMN_SYNONYMS["POS"] = ["hm_pos"]
        COLUMN_SYNONYMS["ID"] = ["hm_rsid"]

    try:
        result = inspect_file(args.input, forced_delimiter=args.delimiter)
    except Exception as e:
        logger.exception(e)
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
