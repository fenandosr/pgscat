#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import zarr


def decode_str_array(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)

    if arr.dtype.kind == "S":
        return np.char.decode(arr, "utf-8")

    if arr.dtype.kind == "U":
        return arr

    if arr.dtype.kind == "O":
        out = []
        for x in arr:
            if isinstance(x, (bytes, bytearray)):
                out.append(x.decode("utf-8"))
            else:
                out.append("" if x is None else str(x))
        return np.array(out, dtype=object)

    out = []
    for x in arr.tolist():
        if isinstance(x, (bytes, bytearray)):
            out.append(x.decode("utf-8"))
        else:
            out.append("" if x is None else str(x))
    return np.array(out, dtype=object)


def ensure_schema(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        """
        create table if not exists variants (
            chrom varchar,
            pos bigint,
            ref varchar,
            alt varchar,
            variant_idx bigint,
            zarr_path varchar
        )
    """
    )
    con.execute(
        "create index if not exists variants_chr_pos_idx on variants(chrom, pos)"
    )
    con.execute(
        "create index if not exists variants_chr_pos_ref_alt_idx on variants(chrom, pos, ref, alt)"
    )


def build_variants_df(
    zarr_path: Path,
    chrom: str,
    pos_array: str = "variant_position",
    allele_array: str = "variant_allele",
    use_alleles: bool = True,
) -> pd.DataFrame:
    root = zarr.open_group(zarr_path, mode="r")

    if pos_array not in root:
        raise KeyError(f"No existe array '{pos_array}' en {zarr_path}")

    pos = np.asarray(root[pos_array][:], dtype=np.int64)
    variant_idx = np.arange(len(pos), dtype=np.int64)

    if use_alleles and allele_array in root:
        alleles = np.asarray(root[allele_array][:, :2])
        ref = decode_str_array(alleles[:, 0])
        alt = decode_str_array(alleles[:, 1])
    else:
        ref = np.array([""] * len(pos), dtype="S")
        alt = np.array([""] * len(pos), dtype="S")

    df = pd.DataFrame(
        {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_idx": variant_idx,
            "zarr_path": str(zarr_path),
        }
    )
    return df


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Construye tabla variants en DuckDB a partir de chrN.zarr"
    )
    ap.add_argument("--zarr-base-dir", required=True)
    ap.add_argument("--duckdb", required=True)
    ap.add_argument("--pos-array", default="variant_position")
    ap.add_argument("--allele-array", default="variant_allele")
    ap.add_argument(
        "--chroms",
        nargs="*",
        default=[f"chr{i}" for i in range(1, 23)] + ["chrX"],
    )
    ap.add_argument(
        "--no-alleles",
        action="store_true",
        help="No leer variant_allele; deja ref/alt vacíos",
    )
    ap.add_argument(
        "--append",
        action="store_true",
        help="No borrar variants existente",
    )
    args = ap.parse_args()

    zarr_base_dir = Path(args.zarr_base_dir).resolve()
    duckdb_path = Path(args.duckdb).resolve()
    duckdb_path.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(str(duckdb_path))
    ensure_schema(con)

    if not args.append:
        con.execute("delete from variants")

    summary = []

    for chrom in args.chroms:
        zarr_path = zarr_base_dir / f"{chrom}.zarr"
        if not zarr_path.exists():
            summary.append(
                {"chrom": chrom, "zarr_path": str(zarr_path), "status": "missing"}
            )
            continue

        df = build_variants_df(
            zarr_path=zarr_path,
            chrom=chrom,
            pos_array=args.pos_array,
            allele_array=args.allele_array,
            use_alleles=not args.no_alleles,
        )

        con.register("tmp_variants_df", df)
        con.execute(
            """
            insert into variants
            select chrom, pos, ref, alt, variant_idx, zarr_path
            from tmp_variants_df
        """
        )
        con.unregister("tmp_variants_df")

        summary.append(
            {
                "chrom": chrom,
                "zarr_path": str(zarr_path),
                "status": "loaded",
                "n_variants": int(len(df)),
            }
        )

    total = con.execute("select count(*) from variants").fetchone()[0]
    con.close()

    print(json.dumps({"n_total_variants": int(total), "chroms": summary}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
