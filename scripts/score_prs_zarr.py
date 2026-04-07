#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import zarr


logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(name)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )


def open_text_auto(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


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

    # fallback robusto, evitando astype(str)
    out = []
    for x in arr.tolist():
        if isinstance(x, (bytes, bytearray)):
            out.append(x.decode("utf-8"))
        else:
            out.append("" if x is None else str(x))
    return np.array(out, dtype=object)


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def infer_betamap_path(prs_path: Path) -> Path:
    name = prs_path.name
    for suffix in (".txt.gz", ".tsv.gz", ".csv.gz", ".gz"):
        if name.endswith(suffix):
            candidate = prs_path.with_name(name[: -len(suffix)] + ".betamap.tsv.gz")
            if candidate.exists():
                return candidate

    alt = prs_path.with_name(prs_path.stem + ".betamap.tsv.gz")
    if alt.exists():
        return alt

    raise FileNotFoundError(f"No se encontró betamap junto a {prs_path}")


def read_betamap(betamap_path: Path) -> pd.DataFrame:
    required = [
        "PRS_ID",
        "CHROM",
        "POS",
        "ID",
        "EFFECT_ALLELE",
        "OTHER_ALLELE",
        "BETA",
        "IS_FLIP",
    ]

    df = pd.read_csv(
        betamap_path,
        sep="\t",
        compression="gzip" if str(betamap_path).endswith(".gz") else None,
        dtype={
            "PRS_ID": "string",
            "CHROM": "string",
            "POS": "Int64",
            "ID": "string",
            "EFFECT_ALLELE": "string",
            "OTHER_ALLELE": "string",
            "BETA": "float64",
            "IS_FLIP": "Int64",
        },
    )

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Faltan columnas requeridas en betamap: {missing}")

    df = df[required].copy()
    df["CHROM"] = df["CHROM"].map(normalize_chrom)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["IS_FLIP"] = pd.to_numeric(df["IS_FLIP"], errors="coerce")

    df = df.dropna(subset=["CHROM", "POS", "BETA", "IS_FLIP"])
    df["POS"] = df["POS"].astype(np.int64)
    df["IS_FLIP"] = df["IS_FLIP"].astype(np.int8)

    df = df.sort_values(["CHROM", "POS"]).reset_index(drop=True)
    return df


def get_sample_ids(root: zarr.Group, sample_array: str) -> np.ndarray:
    if sample_array not in root:
        raise KeyError(f"No existe sample array '{sample_array}' en Zarr")
    return decode_str_array(root[sample_array][:])


def get_variant_positions(root: zarr.Group, pos_array: str) -> np.ndarray:
    if pos_array not in root:
        raise KeyError(f"No existe variant position array '{pos_array}' en Zarr")
    return np.asarray(root[pos_array][:], dtype=np.int64)


def get_dosage_array(root: zarr.Group, dosage_array: str):
    if dosage_array not in root:
        raise KeyError(f"No existe dosage array '{dosage_array}' en Zarr")

    arr = root[dosage_array]

    if len(arr.shape) == 2:
        return arr

    if len(arr.shape) == 3 and arr.shape[2] == 1:
        return arr

    raise ValueError(
        f"Se esperaba un dosage array 2D o 3D con última dimensión = 1, recibido shape={arr.shape}"
    )


def read_dosage_chunk(dosage, idx_chunk: np.ndarray) -> np.ndarray:
    if len(dosage.shape) == 2:
        ds = np.asarray(dosage.oindex[idx_chunk, :], dtype=np.float32)
    elif len(dosage.shape) == 3 and dosage.shape[2] == 1:
        ds = np.asarray(dosage.oindex[idx_chunk, :, :], dtype=np.float32)
        ds = ds[:, :, 0]
    else:
        raise ValueError(f"Dosage array inesperado, shape={dosage.shape}")

    if ds.ndim != 2:
        raise ValueError(f"Dosage chunk inesperado, shape={ds.shape}")

    return ds


def match_positions_searchsorted(
    prs_chr: pd.DataFrame,
    zarr_positions: np.ndarray,
) -> pd.DataFrame:
    """
    Hace match exacto por POS usando searchsorted.
    Si hay duplicados en Zarr para una POS, toma la primera coincidencia.
    """
    qpos = prs_chr["POS"].to_numpy(dtype=np.int64)

    left = np.searchsorted(zarr_positions, qpos, side="left")
    right = np.searchsorted(zarr_positions, qpos, side="right")

    matched_rows = []

    for i, (l, r) in enumerate(zip(left, right, strict=False)):
        if l == r:
            continue

        zidx = int(l)
        row = prs_chr.iloc[i]

        matched_rows.append(
            {
                "PRS_ID": row["PRS_ID"],
                "CHROM": row["CHROM"],
                "POS": int(row["POS"]),
                "ID": row["ID"],
                "EFFECT_ALLELE": row["EFFECT_ALLELE"],
                "OTHER_ALLELE": row["OTHER_ALLELE"],
                "BETA": float(row["BETA"]),
                "IS_FLIP": int(row["IS_FLIP"]),
                "zidx": zidx,
            }
        )

    if not matched_rows:
        return pd.DataFrame(
            columns=[
                "PRS_ID",
                "CHROM",
                "POS",
                "ID",
                "EFFECT_ALLELE",
                "OTHER_ALLELE",
                "BETA",
                "IS_FLIP",
                "zidx",
            ]
        )

    out = pd.DataFrame(matched_rows)
    out = out.sort_values("zidx").reset_index(drop=True)
    return out


def score_chromosome(
    zarr_path: Path,
    prs_chr: pd.DataFrame,
    dosage_array_name: str,
    sample_array_name: str,
    pos_array_name: str,
    variant_chunk_size: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Devuelve:
      sample_ids, partial_scores, stats
    """
    root = zarr.open_group(zarr_path, mode="r")

    sample_ids = get_sample_ids(root, sample_array_name)
    zarr_positions = get_variant_positions(root, pos_array_name)
    dosage = get_dosage_array(root, dosage_array_name)

    if dosage.shape[1] != len(sample_ids):
        raise ValueError(
            f"Inconsistencia en {zarr_path}: dosage.shape[1]={dosage.shape[1]} != n_samples={len(sample_ids)}"
        )

    matched = match_positions_searchsorted(prs_chr, zarr_positions)

    scores = np.zeros(len(sample_ids), dtype=np.float64)

    if matched.empty:
        stats = {
            "chrom": prs_chr["CHROM"].iloc[0] if not prs_chr.empty else None,
            "zarr_path": str(zarr_path),
            "n_prs_variants_input_chr": int(len(prs_chr)),
            "n_prs_variants_matched_chr": 0,
        }
        return sample_ids, scores, stats

    zidx = matched["zidx"].to_numpy(dtype=np.int64)
    beta = matched["BETA"].to_numpy(dtype=np.float64)
    is_flip = matched["IS_FLIP"].to_numpy(dtype=np.int8).astype(bool)

    n = len(matched)
    for start in range(0, n, variant_chunk_size):
        end = min(start + variant_chunk_size, n)

        idx_chunk = zidx[start:end]
        beta_chunk = beta[start:end]
        flip_chunk = is_flip[start:end]

        ds = read_dosage_chunk(dosage, idx_chunk)
        ds = np.nan_to_num(ds, nan=0.0)

        contrib = ds.copy()
        if np.any(flip_chunk):
            contrib[flip_chunk, :] = 2.0 - contrib[flip_chunk, :]

        scores += (contrib * beta_chunk[:, None]).sum(axis=0)

    stats = {
        "chrom": str(prs_chr["CHROM"].iloc[0]),
        "zarr_path": str(zarr_path),
        "n_prs_variants_input_chr": int(len(prs_chr)),
        "n_prs_variants_matched_chr": int(len(matched)),
    }
    return sample_ids, scores, stats


def write_scores_tsv(out_path: Path, sample_ids: np.ndarray, prs: np.ndarray) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    compression = "gzip" if str(out_path).endswith(".gz") else None
    df = pd.DataFrame({"sample_id": sample_ids, "prs": prs})
    df.to_csv(out_path, sep="\t", index=False, compression=compression)


def build_output_paths(prs_path: Path, out_prefix: str | None) -> tuple[Path, Path]:
    if out_prefix:
        base = Path(out_prefix)
        scores_path = base.with_suffix("") if base.suffix == ".gz" else base
        if not str(scores_path).endswith(".tsv.gz"):
            scores_path = Path(str(scores_path) + ".tsv.gz")
        meta_path = Path(str(scores_path) + ".meta.json")
        return scores_path, meta_path

    name = prs_path.name
    for suffix in (".txt.gz", ".tsv.gz", ".csv.gz", ".gz"):
        if name.endswith(suffix):
            stem = name[: -len(suffix)]
            scores_path = prs_path.with_name(stem + ".prs.tsv.gz")
            meta_path = prs_path.with_name(stem + ".prs.tsv.gz.meta.json")
            return scores_path, meta_path

    scores_path = prs_path.with_name(prs_path.name + ".prs.tsv.gz")
    meta_path = prs_path.with_name(prs_path.name + ".prs.tsv.gz.meta.json")
    return scores_path, meta_path


def run_prs(
    prs_path: Path,
    betamap_path: Path,
    zarr_base_dir: Path,
    dosage_array_name: str,
    sample_array_name: str,
    pos_array_name: str,
    variant_chunk_size: int,
    out_prefix: str | None,
) -> dict:
    prs_df = read_betamap(betamap_path)
    if prs_df.empty:
        raise ValueError(f"Betamap vacío: {betamap_path}")

    prs_id_values = prs_df["PRS_ID"].dropna().unique().tolist()
    prs_id = prs_id_values[0] if prs_id_values else prs_path.stem

    chroms = prs_df["CHROM"].dropna().unique().tolist()
    chroms = sorted(chroms, key=lambda x: (x == "chrX", x))

    total_scores = None
    total_sample_ids = None
    chrom_stats = []

    for chrom in chroms:
        prs_chr = prs_df[prs_df["CHROM"] == chrom].copy()
        zarr_path = zarr_base_dir / f"{chrom}.zarr"

        if not zarr_path.exists():
            logger.warning("No existe Zarr para %s: %s", chrom, zarr_path)
            chrom_stats.append(
                {
                    "chrom": chrom,
                    "zarr_path": str(zarr_path),
                    "n_prs_variants_input_chr": int(len(prs_chr)),
                    "n_prs_variants_matched_chr": 0,
                    "missing_zarr": True,
                }
            )
            continue

        logger.info("Scoring %s en %s", chrom, zarr_path)

        sample_ids, partial_scores, stats = score_chromosome(
            zarr_path=zarr_path,
            prs_chr=prs_chr,
            dosage_array_name=dosage_array_name,
            sample_array_name=sample_array_name,
            pos_array_name=pos_array_name,
            variant_chunk_size=variant_chunk_size,
        )

        chrom_stats.append(stats)

        if total_scores is None:
            total_scores = partial_scores
            total_sample_ids = sample_ids
        else:
            if len(sample_ids) != len(total_sample_ids) or not np.array_equal(
                sample_ids, total_sample_ids
            ):
                raise ValueError(
                    f"Los sample_ids difieren entre cromosomas. Revisa {chrom}."
                )
            total_scores += partial_scores

    if total_scores is None or total_sample_ids is None:
        raise ValueError("No se pudo scorear ningún cromosoma")

    scores_path, meta_path = build_output_paths(prs_path, out_prefix)
    write_scores_tsv(scores_path, total_sample_ids, total_scores)

    summary = {
        "prs_id": prs_id,
        "prs_path": str(prs_path),
        "betamap_path": str(betamap_path),
        "zarr_base_dir": str(zarr_base_dir),
        "dosage_array_name": dosage_array_name,
        "sample_array_name": sample_array_name,
        "pos_array_name": pos_array_name,
        "variant_chunk_size": variant_chunk_size,
        "n_samples": int(len(total_sample_ids)),
        "n_chromosomes_seen": int(len(chroms)),
        "chrom_stats": chrom_stats,
        "output_scores": str(scores_path),
    }

    meta_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Calcula un PRS sobre stores Zarr separados por cromosoma"
    )
    ap.add_argument(
        "prs_path",
        help="Path del PRS crudo, e.g. /path/PRS.txt.gz",
    )
    ap.add_argument(
        "--betamap",
        default=None,
        help="Path explícito al betamap. Si no se da, se infiere junto al PRS.",
    )
    ap.add_argument(
        "--zarr-base-dir",
        required=True,
        help="Directorio base con chrN.zarr, e.g. /mnt/.../zar_files",
    )
    ap.add_argument(
        "--dosage-array",
        default="call_DS",
        help="Nombre del array de dosage en Zarr",
    )
    ap.add_argument(
        "--sample-array",
        default="sample_id",
        help="Nombre del array de sample IDs en Zarr",
    )
    ap.add_argument(
        "--pos-array",
        default="variant_position",
        help="Nombre del array de posiciones en Zarr",
    )
    ap.add_argument(
        "--variant-chunk-size",
        type=int,
        default=2048,
        help="Número de variantes por bloque al leer dosage",
    )
    ap.add_argument(
        "--out-prefix",
        default=None,
        help="Prefijo/path de salida. Si no se da, escribe junto al PRS",
    )
    ap.add_argument(
        "-v",
        "--verbose",
        action="store_true",
    )

    args = ap.parse_args()
    setup_logging(verbose=args.verbose)

    prs_path = Path(args.prs_path).resolve()
    if not prs_path.exists():
        raise SystemExit(f"No existe PRS: {prs_path}")

    betamap_path = (
        Path(args.betamap).resolve() if args.betamap else infer_betamap_path(prs_path)
    )
    if not betamap_path.exists():
        raise SystemExit(f"No existe betamap: {betamap_path}")

    zarr_base_dir = Path(args.zarr_base_dir).resolve()
    if not zarr_base_dir.exists():
        raise SystemExit(f"No existe zarr base dir: {zarr_base_dir}")

    summary = run_prs(
        prs_path=prs_path,
        betamap_path=betamap_path,
        zarr_base_dir=zarr_base_dir,
        dosage_array_name=args.dosage_array,
        sample_array_name=args.sample_array,
        pos_array_name=args.pos_array,
        variant_chunk_size=args.variant_chunk_size,
        out_prefix=args.out_prefix,
    )

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
