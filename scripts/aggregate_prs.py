#!/usr/bin/env python3
"""
aggregate_prs.py — Concatena los PRS por cromosoma en un único archivo final.

Espera los archivos generados por compute_prs.py:
    {SCORES_DIR}/{PGS_ID}/{PGS_ID}_chr{CHROM}_scores.tsv
    {SCORES_DIR}/{PGS_ID}/{PGS_ID}_chr{CHROM}_metadata.json

Produce:
    {SCORES_DIR}/{PGS_ID}/{PGS_ID}_PRS_total.tsv
    {SCORES_DIR}/{PGS_ID}/{PGS_ID}_PRS_total_metadata.json

Uso:
    python aggregate_prs.py --pgs-id PGS000001
    python aggregate_prs.py --pgs-id PGS000001 --chroms 1 2 3 X
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

SCORES_BASE = Path("/mnt/cephfs/hot_nvme/pgscatalog/scores")
ALL_CHROMS  = [str(c) for c in range(1, 23)] + ["X"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Agrega PRS por cromosoma en un total.")
    p.add_argument("--pgs-id",  required=True, help="PGS Catalog ID")
    p.add_argument(
        "--chroms", nargs="+", default=ALL_CHROMS,
        help="Cromosomas a incluir (default: 1-22 X)",
    )
    return p.parse_args()


def load_chr_scores(pgs_id: str, chrom: str) -> pd.DataFrame | None:
    path = SCORES_BASE / pgs_id / f"{pgs_id}_chr{chrom}_scores.tsv"
    if not path.exists():
        log.warning("Archivo no encontrado, se omite: %s", path)
        return None
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"PRS": f"PRS_chr{chrom}"})
    return df.set_index("sample_id")


def load_chr_metadata(pgs_id: str, chrom: str) -> dict:
    path = SCORES_BASE / pgs_id / f"{pgs_id}_chr{chrom}_metadata.json"
    if not path.exists():
        return {}
    with open(path) as fh:
        return json.load(fh)


def main() -> None:
    args   = parse_args()
    pgs_id = args.pgs_id
    chroms = args.chroms
    t_start = datetime.now()

    out_dir = SCORES_BASE / pgs_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── cargar scores por cromosoma ───────────────────────────────────────────
    frames:   list[pd.DataFrame] = []
    meta_per_chr: dict[str, dict] = {}
    missing_chrs: list[str]      = []

    for chrom in chroms:
        df = load_chr_scores(pgs_id, chrom)
        if df is None:
            missing_chrs.append(chrom)
            continue
        frames.append(df)
        meta_per_chr[chrom] = load_chr_metadata(pgs_id, chrom)
        log.info("chr%s cargado: %d muestras", chrom, len(df))

    if not frames:
        log.error("No se encontraron archivos de scores — abortando.")
        sys.exit(1)

    # ── verificar que todas las frames tienen las mismas muestras ─────────────
    reference_samples = frames[0].index
    for df in frames[1:]:
        if not df.index.equals(reference_samples):
            log.warning(
                "Diferencia en muestras entre cromosomas — "
                "se usará outer join y NaN para muestras faltantes."
            )
            break

    # ── combinar y calcular PRS total ─────────────────────────────────────────
    merged = pd.concat(frames, axis=1, join="outer")
    chr_cols = [c for c in merged.columns if c.startswith("PRS_chr")]
    merged["PRS_total"] = merged[chr_cols].sum(axis=1, skipna=True)
    merged = merged.reset_index()

    # ── reordenar columnas: sample_id, PRS_total, PRS_chr1...chrX ────────────
    ordered_cols = ["sample_id", "PRS_total"] + chr_cols
    merged = merged[ordered_cols]

    # ── estadísticas globales ─────────────────────────────────────────────────
    total = merged["PRS_total"]
    elapsed = (datetime.now() - t_start).total_seconds()

    # ── guardar scores ────────────────────────────────────────────────────────
    scores_out = out_dir / f"{pgs_id}_PRS_total.tsv"
    merged.to_csv(scores_out, sep="\t", index=False, float_format="%.8f")
    log.info("PRS total guardado en: %s", scores_out)

    # ── metadatos globales ────────────────────────────────────────────────────
    global_meta = {
        "pgs_id":              pgs_id,
        "aggregation_timestamp": t_start.isoformat(),
        "elapsed_seconds":     round(elapsed, 2),
        "chroms_requested":    chroms,
        "chroms_found":        [c for c in chroms if c not in missing_chrs],
        "chroms_missing":      missing_chrs,
        "n_samples":           int(len(merged)),
        "n_chr_scores":        len(frames),
        "prs_total_mean":      float(total.mean()),
        "prs_total_std":       float(total.std()),
        "prs_total_min":       float(total.min()),
        "prs_total_max":       float(total.max()),
        "per_chromosome":      meta_per_chr,
    }

    meta_out = out_dir / f"{pgs_id}_PRS_total_metadata.json"
    with open(meta_out, "w") as fh:
        json.dump(global_meta, fh, indent=2, default=str)
    log.info("Metadatos globales guardados en: %s", meta_out)

    # ── resumen ───────────────────────────────────────────────────────────────
    log.info(
        "Agregación completada en %.1f s | "
        "muestras=%d  cromosomas=%d  PRS_total μ=%.4f σ=%.4f",
        elapsed, len(merged), len(frames), total.mean(), total.std(),
    )
    if missing_chrs:
        log.warning("Cromosomas omitidos por archivos faltantes: %s", missing_chrs)


if __name__ == "__main__":
    main()
