#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import gzip
import sys
from typing import TextIO, Any
import logging
import datetime
import time
import csv

import pysam


COMMENT_PREFIXES = ("#", "//")


def setup_logging(log_dir: Path):
    log_dir.mkdir(parents=True, exist_ok=True)

    run_id = datetime.datetime.now().strftime("%Y%m%d_%H")
    run_log = log_dir / f"runs_{run_id}.log"
    error_log = log_dir / "error.log"
    miss_log = log_dir / "missing_variants.log"

    fmt = logging.Formatter(
        "%(asctime)s.%(msecs)03d | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    fmt.converter = time.gmtime  # UTC

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # evitar duplicados si se llama varias veces
    if root_logger.handlers:
        return

    # stdout
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)

    # run log
    fh_run = logging.FileHandler(run_log)
    fh_run.setFormatter(fmt)

    # error log
    fh_err = logging.FileHandler(error_log)
    fh_err.setLevel(logging.ERROR)
    fh_err.setFormatter(fmt)

    # missing
    fmt = logging.Formatter(
        "%(asctime)s.%(msecs)03d %(message)s", datefmt="%Y-%m-%dT%H:%M:%S"
    )

    mis = logging.getLogger("missing_variants")
    mis.setLevel(logging.ERROR)

    fh_mis = logging.FileHandler(miss_log)
    fh_mis.setFormatter(fmt)

    mis.addHandler(fh_mis)
    mis.propagate = False

    root_logger.addHandler(sh)
    root_logger.addHandler(fh_run)
    root_logger.addHandler(fh_err)
    root_logger.addHandler(mis)


def open_text_auto(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", newline="")
    return open(p, "rt", encoding="utf-8", newline="")


def detect_delimiter(sample: str, forced: str | None = None) -> str:
    if forced:
        if forced == "tab":
            return "\t"
        if forced == "comma":
            return ","
        else:
            return forced
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,; ")
        return dialect.delimiter
    except csv.Error:
        raise ValueError(f"Unsupported delimiter mode: {forced}")


def file_delim_meta_header_data(
    path: str,
    forced_delimiter: str | None = None,
    meta_line_startswith: str | None = None,
    skip_meta_line_startswith: str | None = None,
    meta_kv_separator: str = None,
    skip_data: bool = False,
) -> tuple[str, dict[str | None, str | None], list[str], list[list[str]]]:
    LINE_I_KEY = "line{num}"
    meta: dict[str, str] = {}
    with open_text_auto(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if skip_meta_line_startswith and stripped.startswith(
                skip_meta_line_startswith
            ):
                continue
            if meta_line_startswith and stripped.startswith(meta_line_startswith):
                payload = stripped.lstrip("#").strip()
                if meta_kv_separator and meta_kv_separator in payload:
                    key, value = payload.split(meta_kv_separator, 1)
                    meta[key.strip()] = value.strip()
                else:
                    meta[LINE_I_KEY.format(num=len(meta) + 1)] = payload
                continue
            delim = detect_delimiter(line, forced_delimiter)
            header = next(csv.reader([line], delimiter=delim))
            break
        else:
            raise ValueError("No header found!")
        if skip_data:
            return delim, meta, header, []
        reader = csv.reader(handle, delimiter=delim)
        data = [row for row in reader if row and any(cell.strip() for cell in row)]

    return delim, meta, header, data


def normalize_name(name: str) -> str:
    return "".join(ch for ch in name.strip().lower() if ch.isalnum() or ch in "_")


def build_lookup(header: list[str]) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for idx, col in enumerate(header):
        norm = normalize_name(col)
        if norm and norm not in lookup:
            lookup[norm] = idx
    return lookup


def find_column(
    lookup: dict[str, int], candidates: list[str]
) -> tuple[int | None, str | None]:
    res = None, None
    for candidate in candidates:
        key = normalize_name(candidate)
        if key in lookup:
            res = lookup[key], candidate
    return res


def normalize_chrom(value: str) -> str:
    s = value.strip()
    if not s:
        raise ValueError("Empty chrom")
    if s.startswith("chr"):
        return s
    return f"chr{s}"


def normalize_allele(value: str) -> str:
    return value.strip().upper()


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


def query_fasta(
    fasta: pysam.FastaFile,
    chrom: str,
    pos: int,
    effect_allele: str,
    other_allele: str | None,
    flank_bp: int = 5,
) -> dict[str, Any]:
    """
    Resuelve REF/ALT/IS_FLIP usando una sola consulta al FASTA.

    Convenciones:
    - pos es 1-based
    - para SNVs, REF es la base en pos
    - para indels, REF se toma como el segmento de referencia de longitud
      max(len(effect_allele), len(other_allele)) empezando en pos

    Regresa:
      {
        "ok": bool,
        "chrom": str,
        "pos": int,
        "ref": str,
        "alt": str | None,
        "is_flip": int | None,
        "ref_flank": str,
        "reason": str | None,
      }
    """
    c = normalize_chrom_for_fasta(chrom, fasta)

    if not effect_allele or not other_allele:
        return {
            "ok": False,
            "chrom": c,
            "pos": pos,
            "ref": None,
            "alt": None,
            "is_flip": None,
            "ref_flank": None,
            "reason": "missing_allele",
        }

    max_len = max(len(effect_allele), len(other_allele))

    # ventana 0-based half-open para una sola consulta
    start0 = max(0, pos - 1 - flank_bp)
    end0 = pos - 1 + max_len + flank_bp

    ref_flank = fasta.fetch(c, start0, end0).upper()

    # offset de la variante dentro de ref_flank
    offset = (pos - 1) - start0

    # segmento REF derivado del flank, sin otra consulta al FASTA
    ref = ref_flank[offset : offset + max_len]

    if len(ref) != max_len:
        return {
            "ok": False,
            "chrom": c,
            "pos": pos,
            "ref": ref,
            "alt": None,
            "is_flip": None,
            "ref_flank": ref_flank,
            "reason": "reference_out_of_range",
        }

    if effect_allele == ref:
        return {
            "ok": True,
            "chrom": c,
            "pos": pos,
            "ref": ref,
            "alt": other_allele,
            "is_flip": 1,
            "ref_flank": ref_flank,
            "reason": None,
        }

    if other_allele == ref:
        return {
            "ok": True,
            "chrom": c,
            "pos": pos,
            "ref": ref,
            "alt": effect_allele,
            "is_flip": 0,
            "ref_flank": ref_flank,
            "reason": None,
        }

    return {
        "ok": False,
        "chrom": c,
        "pos": pos,
        "ref": ref,
        "alt": None,
        "is_flip": None,
        "ref_flank": ref_flank,
        "reason": "no_ref_match",
    }
