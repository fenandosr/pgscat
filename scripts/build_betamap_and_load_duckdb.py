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
from typing import Any
import pdb
import logging

import duckdb
import pysam

from pgscat.common_utils import (
    setup_logging,
    file_delim_meta_header_data,
    normalize_allele,
    query_fasta,
)

logger = logging.getLogger(__name__)
mis_log = logging.getLogger("missing_variants")


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


def coerce_beta(raw_weight: str, weight_type: str) -> tuple[float, str]:
    """
    Regresa (beta, beta_transform).
    Convención:
      beta / effect size / additive beta / log(OR) / log(HR) => identity
      OR / odds ratio / HR / hazard ratio                  => ln(x)
    """
    x = float(raw_weight)
    wt = weight_type or "beta"

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


def fetch_ref_base(
    fasta: pysam.FastaFile, chrom: str, pos1: int, length: int | None = 1
) -> str:
    # pysam fetch usa coordinates 0-based half-open
    return fasta.fetch(chrom, pos1 - 1, pos1 + length - 1).upper()


def write_betamap(
    prs_path: Path,
    fasta_path: Path,
    colsinfo_path: Path,
    out_path: Path,
    explicit_prs_id: str | None = None,
) -> dict[str, Any]:
    colsinfo = json.loads(colsinfo_path.read_text(encoding="utf-8"))

    fasta = pysam.FastaFile(fasta_path)

    meta = colsinfo.get("meta")

    required = ["CHROM", "POS", "EFFECT_ALLELE", "OTHER_ALLELE", "EFFECT_WEIGHT"]
    missing = [x for x in required if x in colsinfo.get("missing")]
    if missing:
        raise ValueError(
            f"En el prs {prs_path} faltan columnas requeridas en colsinfo/header: {missing}"
        )

    prs_id = infer_prs_id(prs_path, meta, explicit_prs_id)
    weight_type = (
        meta.get("weight_type", meta.get("Weight Type", "beta")).strip() or "beta"
    )

    idx = {k: v["index_0based"] for k, v in colsinfo["columns"].items()}

    delim = colsinfo.get("delimiter")
    if delim == "\\t":
        delim = "\t"

    _d, _m, _h, data = file_delim_meta_header_data(
        prs_path,
        forced_delimiter=delim,
        meta_line_startswith="#",
        skip_meta_line_startswith="##",
        meta_kv_separator="=",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_total = 0
    n_written = 0
    n_parse_error = 0
    n_bad_allele = 0
    n_not_translated = 0
    n_skipped_ref_mismatch = 0
    betamap_cols = [
        "PRS_ID",
        "CHROM",
        "POS",
        "ID",
        "EFFECT_ALLELE",
        "OTHER_ALLELE",
        "BETA",
        "IS_FLIP",
    ]

    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(betamap_cols)

        for row in data:
            n_total += 1
            hm_source = ""
            if "hm_source" in colsinfo["header"]:
                hm_source = row[colsinfo["header"].index("hm_source")].strip()
            if hm_source == "Unknown":
                n_not_translated += 1
                continue
            try:
                chrom = row[idx["CHROM"]]
                if not chrom:
                    n_bad_allele += 1
                    continue
                pos = int(row[idx["POS"]].strip())
                effect_allele = normalize_allele(row[idx["EFFECT_ALLELE"]])
                other_allele = normalize_allele(row[idx["OTHER_ALLELE"]])
                if not other_allele:
                    other_allele = normalize_allele(
                        row[colsinfo["header"].index("other_allele")]
                    )
                raw_weight = row[idx["EFFECT_WEIGHT"]].strip()
                beta, beta_transform = coerce_beta(raw_weight, weight_type)
                if len(effect_allele) == 0 or len(other_allele) == 0:
                    n_bad_allele += 1
                    continue

                variant_id = (
                    str(row[idx["ID"]]).strip() if idx["ID"] is not None else ""
                )

                ref = query_fasta(
                    fasta,
                    chrom,
                    pos,
                    effect_allele,
                    other_allele,
                    flank_bp=10,
                )

                if not ref["ok"]:
                    n_skipped_ref_mismatch += 1
                    mis_log.warning(
                        f"{prs_id} {ref['reason']} {chrom}:{pos} {effect_allele} {other_allele} {ref['ref']} {ref['alt']} {ref['ref_flank']}"
                    )
                    continue

                writer.writerow(
                    [
                        prs_id,
                        chrom,
                        pos,
                        variant_id,
                        effect_allele,
                        other_allele,
                        beta,
                        ref["is_flip"],
                    ]
                )
                n_written += 1

            except Exception as err:
                logger.error(f"Error processing row: {err} Row content: {row}")
                n_parse_error += 1
                continue

    meta_out = {
        "prs_id": prs_id,
        "prs_path": str(prs_path),
        "colsinfo_path": str(colsinfo_path),
        "betamap_path": str(out_path),
        "weight_type": weight_type,
        "selected_columns": betamap_cols,
        "n_total_rows": n_total,
        "n_written_rows": n_written,
        "n_parse_error": n_parse_error,
        "n_bad_allele": n_bad_allele,
        "n_not_translated": n_not_translated,
        "n_skipped_ref_mismatch": n_skipped_ref_mismatch,
    }

    meta_json = out_path.with_suffix(out_path.suffix + ".meta.json")
    meta_json.write_text(
        json.dumps(meta_out, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    return meta_out


def ensure_duckdb_schema(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        """
        create table if not exists prs_variants (
            prs_id varchar,
            chrom varchar,
            pos bigint,
            id varchar,
            effect_allele varchar,
            other_allele varchar,
            beta double,
            is_flip boolean,
            source_file varchar
        )
    """
    )

    # índices
    con.execute(
        "create index if not exists prs_variants_prs_idx on prs_variants(prs_id)"
    )
    con.execute(
        "create index if not exists prs_variants_lookup_idx on prs_variants(prs_id, chrom, pos)"
    )


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
        con.execute(
            f"""
            create or replace temp table tmp_betamap as
            select *
            from read_csv_auto(
                '{str(betamap_path).replace("'", "''")}',
                delim='\t',
                header=true,
                compression='gzip',
                types={{'CHROM': 'VARCHAR'}}
            )
        """
        )

        con.execute(
            """
            insert into prs_variants
            select
                PRS_ID as prs_id,
                CHROM as chrom,
                POS as pos,
                ID as id,
                EFFECT_ALLELE as effect_allele,
                OTHER_ALLELE as other_allele,
                BETA as beta,
                IS_FLIP as is_flip,
                ? as source_file
            from tmp_betamap
        """,
            [str(betamap_path)],
        )

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
    ap.add_argument("--fasta", required=True, help="Referencia FASTA, e.g. hg38.fa")
    ap.add_argument(
        "--prs-id", default=None, help="PRS ID explícito; si no, se infiere"
    )
    ap.add_argument("--colsinfo", default=None, help="Ruta explícita al *.columns.json")
    ap.add_argument("--out", default=None, help="Ruta explícita a *.betamap.tsv.gz")
    ap.add_argument("--log-dir", default=Path.cwd(), help="Ruta al directorio de logs")
    args = ap.parse_args()
    setup_logging(Path(args.log_dir))

    prs_path = Path(args.prs_file).resolve()
    if not prs_path.exists():
        raise SystemExit(f"No existe el archivo PRS: {prs_path}")
    fasta_path = Path(args.fasta).resolve()
    if not fasta_path.exists():
        raise SystemExit(f"No existe el archivo fasta: {fasta_path}")

    colsinfo_path = (
        Path(args.colsinfo).resolve()
        if args.colsinfo
        else find_colsinfo_for_prs(prs_path)
    )
    out_path = Path(args.out).resolve() if args.out else default_betamap_path(prs_path)

    build_info = write_betamap(
        prs_path=prs_path,
        fasta_path=fasta_path,
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
