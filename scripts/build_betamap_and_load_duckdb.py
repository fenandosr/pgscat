#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Iterator, TextIO

import duckdb


def open_text_auto(path: str | Path) -> TextIO:
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def normalize_header_name(name: str) -> str:
    return "".join(ch for ch in str(name).strip().lower() if ch.isalnum())


def parse_meta_and_rows(
    prs_path: Path,
    delimiter: str,
) -> tuple[dict[str, str], list[str], list[list[str]]]:
    meta: dict[str, str] = {}
    header: list[str] | None = None
    rows: list[list[str]] = []

    with open_text_auto(prs_path) as fh:
        reader = csv.reader(fh, delimiter=delimiter)

        for raw_line in fh:
            line = raw_line.rstrip("\n")
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("##"):
                continue

            if stripped.startswith("#") and "=" in stripped:
                payload = stripped.lstrip("#").strip()
                if "=" in payload:
                    key, value = payload.split("=", 1)
                    meta[key.strip()] = value.strip()
                continue

            header = next(csv.reader([line], delimiter=delimiter))
            break

        if header is None:
            raise ValueError(f"No se encontró header de datos en {prs_path}")

        for raw_line in fh:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("##"):
                continue
            if stripped.startswith("#") and "=" in stripped:
                payload = stripped.lstrip("#").strip()
                if "=" in payload:
                    key, value = payload.split("=", 1)
                    meta[key.strip()] = value.strip()
                continue

            row = next(csv.reader([raw_line], delimiter=delimiter))
            rows.append(row)

    return meta, header, rows


def resolve_columns(colsinfo: dict[str, Any], header: list[str]) -> dict[str, str]:
    """
    Devuelve mapping lógico -> nombre real de columna en el archivo.
    Usa columns[*].select si existe; fallback a name_in_file.
    """
    header_set = set(header)
    resolved: dict[str, str] = {}

    for logical_name, entry in colsinfo.get("columns", {}).items():
        select_name = entry.get("select")
        if select_name and select_name in header_set:
            resolved[logical_name] = select_name
            continue

        name_in_file = entry.get("name_in_file")
        if name_in_file and name_in_file in header_set:
            resolved[logical_name] = name_in_file
            continue

    return resolved


def infer_prs_id(prs_path: Path, meta: dict[str, str], explicit: str | None) -> str:
    if explicit:
        return explicit

    for key in ("pgs_id", "PGS_ID", "id", "ID"):
        if key in meta and meta[key].strip():
            return meta[key].strip()

    m = re.search(r"(PGS\d+)", prs_path.name)
    if m:
        return m.group(1)

    raise ValueError("No se pudo inferir prs_id; usa --prs-id")


def normalize_chrom(value: str) -> str:
    s = str(value).strip()
    if not s:
        return s
    if s.startswith("chr"):
        return s
    return f"chr{s}"


def normalize_allele(value: str) -> str:
    return str(value).strip().upper()


def coerce_beta(raw_weight: str, weight_type: str) -> tuple[float, str]:
    """
    Regresa (beta, beta_transform).
    Convención:
      beta / effect size / additive beta / log(OR) / log(HR) => identity
      OR / odds ratio / HR / hazard ratio                  => ln(x)
    """
    x = float(raw_weight)
    wt = normalize_header_name(weight_type or "beta")

    identity_types = {
        "beta",
        "effectweight",
        "effectsize",
        "betaeffect",
        "logor",
        "logoddsratio",
        "loghr",
        "loghazardratio",
    }

    log_types = {
        "or",
        "oddsratio",
        "odds",
        "hr",
        "hazardratio",
    }

    if wt in identity_types:
        return x, "identity"

    if wt in log_types:
        if x <= 0:
            raise ValueError(f"No se puede aplicar ln() a peso <= 0: {x}")
        return math.log(x), "ln"

    # default conservador: tratarlo ya como beta
    return x, "identity_default"


def fetch_ref_base(fasta: pysam.FastaFile, chrom: str, pos1: int, length: Optional[int] = 1) -> str:
    # pysam fetch usa coordinates 0-based half-open
    return fasta.fetch(chrom, pos1 - 1, pos1 + length - 1).upper()


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


def write_betamap(
    prs_path: Path,
    colsinfo_path: Path,
    out_path: Path,
    explicit_prs_id: str | None = None,
) -> dict[str, Any]:
    colsinfo = json.loads(colsinfo_path.read_text(encoding="utf-8"))

    fasta = pysam.FastaFile(fasta_path)

    delimiter = "\t" if colsinfo.get("delimiter") == "\\t" else colsinfo.get("delimiter", "\t")
    meta, header, rows = parse_meta_and_rows(prs_path, delimiter)
    resolved = resolve_columns(colsinfo, header)

    required = ["CHROM", "POS", "EFFECT_ALLELE", "OTHER_ALLELE", "EFFECT_WEIGHT"]
    missing = [x for x in required if x not in resolved]
    if missing:
        raise ValueError(f"Faltan columnas requeridas en colsinfo/header: {missing}")

    prs_id = infer_prs_id(prs_path, meta, explicit_prs_id)
    weight_type = meta.get("weight_type", meta.get("Weight Type", "beta")).strip() or "beta"

    idx = {name: header.index(colname) for name, colname in resolved.items()}

    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_total = 0
    n_written = 0
    n_parse_error = 0
    n_bad_allele = 0
    n_skipped_ref_mismatch = 0

    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "PRS_ID",
                "CHROM",
                "POS",
                "ID",
                "EFFECT_ALLELE",
                "OTHER_ALLELE",
                "BETA",
                "SOURCE_WEIGHT_TYPE",
                "BETA_TRANSFORM",
                "IS_FLIP",
            ]
        )

        id_col = idx.get("ID")

        for row in rows:
            n_total += 1
            try:
                chrom = normalize_chrom(row[idx["CHROM"]])
                pos = int(str(row[idx["POS"]]).strip())
                effect_allele = normalize_allele(row[idx["effect_allele"]])
                other_allele = normalize_allele(row[idx["other_allele"]])
                raw_weight = str(row[idx["effect_weight"]]).strip()
                beta, beta_transform = coerce_beta(raw_weight, weight_type)

                if len(effect_allele) == 0 or len(other_allele) == 0:
                    n_bad_allele += 1
                    continue

                c = normalize_chrom_for_fasta(chrom, fasta)
                if len(effect_allele) > len(other_allele):
                    ref = fetch_ref_base(fasta, c, pos, len(effect_allele))
                elif len(other_allele) > len(effect_allele):
                    ref = fetch_ref_base(fasta, c, pos, len(other_allele))
                else:
                    ref = fetch_ref_base(fasta, c, pos)

                if effect_allele == ref:
                    alt = other_allele
                    is_flip = 1
                elif other_allele == ref:
                    alt = effect_allele
                    is_flip = 0
                else:
                    n_skipped_ref_mismatch += 1
                    continue

                variant_id = str(row[id_col]).strip() if id_col is not None else ""

                writer.writerow(
                    [
                        prs_id,
                        chrom,
                        pos,
                        variant_id,
                        effect_allele,
                        other_allele,
                        beta,
                        weight_type,
                        beta_transform,
                        is_flip,
                    ]
                )
                n_written += 1

            except Exception:
                n_parse_error += 1
                continue

    meta_out = {
        "prs_id": prs_id,
        "prs_path": str(prs_path),
        "colsinfo_path": str(colsinfo_path),
        "betamap_path": str(out_path),
        "weight_type": weight_type,
        "selected_columns": resolved,
        "meta": meta,
        "n_total_rows": n_total,
        "n_written_rows": n_written,
        "n_parse_error": n_parse_error,
        "n_bad_allele": n_bad_allele,
        "n_skipped_ref_mismatch": n_skipped_ref_mismatch,
    }

    meta_json = out_path.with_suffix(out_path.suffix + ".meta.json")
    meta_json.write_text(json.dumps(meta_out, indent=2, ensure_ascii=False), encoding="utf-8")

    return meta_out


def ensure_duckdb_schema(con: duckdb.DuckDBPyConnection) -> None:
    con.execute("""
        create table if not exists prs_variants (
            prs_id varchar,
            chrom varchar,
            pos bigint,
            id varchar,
            effect_allele varchar,
            other_allele varchar,
            beta double,
            source_weight_type varchar,
            beta_transform varchar,
            is_flip boolean,
            source_file varchar
        )
    """)

    # índices
    con.execute("create index if not exists prs_variants_prs_idx on prs_variants(prs_id)")
    con.execute("create index if not exists prs_variants_lookup_idx on prs_variants(prs_id, chrom, pos)")


def load_betamap_into_duckdb(
    db_path: Path,
    betamap_path: Path,
    prs_id: str,
) -> dict[str, Any]:
    con = duckdb.connect(str(db_path))
    try:
        ensure_duckdb_schema(con)

        # reemplazar este prs_id
        con.execute("delete from prs_variants where prs_id = ?", [prs_id])

        # cargar betamap
        con.execute(f"""
            create or replace temp table tmp_betamap as
            select *
            from read_csv_auto(
                '{str(betamap_path).replace("'", "''")}',
                delim='\t',
                header=true,
                compression='gzip'
            )
        """)

        con.execute("""
            insert into prs_variants
            select
                PRS_ID as prs_id,
                CHROM as chrom,
                POS as pos,
                ID as id,
                EFFECT_ALLELE as effect_allele,
                OTHER_ALLELE as other_allele,
                BETA as beta,
                SOURCE_WEIGHT_TYPE as source_weight_type,
                BETA_TRANSFORM as beta_transform,
                IS_FLIP as is_flip,
                ? as source_file
            from tmp_betamap
        """, [str(betamap_path)])

        n_prs_variants = con.execute(
            "select count(*) from prs_variants where prs_id = ?",
            [prs_id],
        ).fetchone()[0]

        return {
            "duckdb_path": str(db_path),
            "prs_id": prs_id,
            "n_prs_variants_loaded": int(n_prs_variants),
        }
    finally:
        con.close()


def find_colsinfo_for_prs(prs_path: Path) -> Path:
    siblings = sorted(prs_path.parent.glob("*.columns.json"))
    if not siblings:
        raise FileNotFoundError(f"No se encontró *.columns.json junto a {prs_path}")
    if len(siblings) == 1:
        return siblings[0]

    prefix = prs_path.name
    # intenta match por prefijo base
    stem = prs_path.name
    for s in siblings:
        if stem.split(".")[0] in s.name:
            return s

    raise FileNotFoundError(
        f"Se encontraron varios *.columns.json y no se pudo elegir uno automáticamente: {siblings}"
    )


def default_betamap_path(prs_path: Path) -> Path:
    name = prs_path.name
    for suffix in (".txt.gz", ".tsv.gz", ".csv.gz", ".gz"):
        if name.endswith(suffix):
            return prs_path.with_name(name[: -len(suffix)] + ".betamap.tsv.gz")
    return prs_path.with_name(name + ".betamap.tsv.gz")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Construye *.betamap.tsv.gz desde PRS crudo + colsinfo y actualiza DuckDB"
    )
    ap.add_argument("prs_file", help="Archivo PRS crudo .gz")
    ap.add_argument("--duckdb", required=True, help="Ruta a DuckDB")
    ap.add_argument("--prs-id", default=None, help="PRS ID explícito; si no, se infiere")
    ap.add_argument("--colsinfo", default=None, help="Ruta explícita al *.columns.json")
    ap.add_argument("--out", default=None, help="Ruta explícita a *.betamap.tsv.gz")
    args = ap.parse_args()

    prs_path = Path(args.prs_file).resolve()
    if not prs_path.exists():
        raise SystemExit(f"No existe el archivo PRS: {prs_path}")

    colsinfo_path = Path(args.colsinfo).resolve() if args.colsinfo else find_colsinfo_for_prs(prs_path)
    out_path = Path(args.out).resolve() if args.out else default_betamap_path(prs_path)

    build_info = write_betamap(
        prs_path=prs_path,
        colsinfo_path=colsinfo_path,
        out_path=out_path,
        explicit_prs_id=args.prs_id,
    )

    db_info = load_betamap_into_duckdb(
        db_path=Path(args.duckdb).resolve(),
        betamap_path=out_path,
        prs_id=build_info["prs_id"],
    )

    summary = {
        "build": build_info,
        "duckdb": db_info,
    }

    json.dump(summary, sys.stdout, indent=2, ensure_ascii=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())