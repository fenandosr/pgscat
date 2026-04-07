#!/usr/bin/env python3
"""
compute_prs.py — Cómputo de Polygenic Risk Score para un cromosoma.

Modelo aditivo:
    PRS = Σ_i  BETA_i × dosage_i

donde dosage_i es el número de copias del alelo de efecto (0, 1 ó 2):
  - IS_FLIP=0 → efecto=ALT  → dosage = suma de alelos ALT (lectura directa del zarr)
  - IS_FLIP=1 → efecto=REF  → dosage = 2 − suma_alelos_ALT

Uso:
    python compute_prs.py --pgs-id PGS000001 --chrom 1
    python compute_prs.py --pgs-id PGS000001 --chrom X
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import zarr

# ── logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

# ── rutas base ────────────────────────────────────────────────────────────────
SCORES_BASE = Path("/mnt/cephfs/hot_nvme/pgscatalog/scores")
ZARR_BASE   = Path("/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files")


# =============================================================================
# Argumentos
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Calcula PRS por cromosoma usando modelo aditivo."
    )
    p.add_argument("--pgs-id",  required=True, help="PGS Catalog ID, ej. PGS000001")
    p.add_argument("--chrom",   required=True, help="Cromosoma: 1-22 o X")
    p.add_argument(
        "--missing-strategy",
        choices=["mean", "zero", "skip"],
        default="mean",
        help="Estrategia para dosage faltante (NaN). "
             "mean=imputar con media, zero=contar como 0, skip=excluir variante. "
             "(default: mean)",
    )
    p.add_argument(
        "--chunk-size",
        type=int,
        default=500,
        help="Número de variantes leídas del zarr por iteración. "
             "Reducir si sigue habiendo OOM; aumentar para mayor velocidad. "
             "(default: 500)",
    )
    return p.parse_args()


# =============================================================================
# Carga de pesos del PRS
# =============================================================================

def load_score_weights(pgs_id: str, chrom: str) -> pd.DataFrame:
    """Lee el archivo de pesos y filtra por cromosoma."""
    path = SCORES_BASE / pgs_id / f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz"
    if not path.exists():
        log.error("Archivo de pesos no encontrado: %s", path)
        sys.exit(1)

    log.info("Leyendo pesos desde %s", path)
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype={"CHROM": str, "POS": np.int64, "IS_FLIP": np.int8, "BETA": np.float64},
    )

    subset = df[df["CHROM"] == str(chrom)].copy().reset_index(drop=True)
    log.info("  Variantes en chr%s: %d", chrom, len(subset))
    return subset


# =============================================================================
# Acceso al zarr
# =============================================================================

def open_zarr_store(chrom: str) -> zarr.Group:
    path = ZARR_BASE / f"chr{chrom}.zarr"
    if not path.exists():
        log.error("Zarr store no encontrado: %s", path)
        sys.exit(1)
    log.info("Abriendo zarr store: %s", path)
    return zarr.open(str(path), mode="r")


def get_samples(store: zarr.Group) -> np.ndarray:
    """Busca el array de IDs de muestras en las rutas más comunes."""
    candidates = ["samples", "sample_id", "sample/id", "calldata/samples"]
    for key in candidates:
        if key in store:
            return np.array(store[key])
    raise KeyError(
        f"No se encontró el array de muestras. Candidatos probados: {candidates}"
    )


def get_variant_arrays(
    store: zarr.Group,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Devuelve (posiciones, ref_alleles, alt_alleles).
    Soporta formato sgkit y formato variants/*.
    """
    # ── sgkit layout ──────────────────────────────────────────────────────────
    if "variant_position" in store:
        pos     = np.array(store["variant_position"])
        alleles = np.array(store["variant_allele"])  # (n_variants, n_alleles)
        return pos, alleles[:, 0], alleles[:, 1]

    # ── variants/* layout ─────────────────────────────────────────────────────
    if "variants/POS" in store:
        pos = np.array(store["variants/POS"])
        ref = np.array(store["variants/REF"])
        alt = np.array(store["variants/ALT"])
        return pos, ref, alt

    raise KeyError(
        "No se pudieron localizar arrays de variantes. "
        "Se esperaba 'variant_position' (sgkit) o 'variants/POS'."
    )


def _detect_dosage_key(store: zarr.Group) -> tuple[str, str]:
    """
    Devuelve (key, layout) donde layout es 'vs' (variants×samples) o 'sv'.
    Lanza KeyError si no hay array compatible.
    """
    for key in ("call_dosage", "calldata/DS", "dosage"):
        if key in store:
            shape = store[key].shape
            # Heurística: la dimensión más larga suele ser variantes
            return key, "vs" if shape[0] >= shape[1] else "sv"
    for key in ("call_genotype", "calldata/GT"):
        if key in store:
            return key, "gt"   # (variants, samples, ploidy)
    raise KeyError(
        "No se encontraron arrays de dosage o genotipo. "
        "Candidatos: call_dosage, calldata/DS, call_genotype, calldata/GT."
    )


def get_dosage_chunk(
    store: zarr.Group,
    dosage_key: str,
    layout: str,
    variant_indices: np.ndarray,
) -> np.ndarray:
    """
    Lee dosage para un subconjunto de índices de variantes.

    Estrategia de acceso eficiente en zarr:
      - Los índices llegan ORDENADOS (garantía del llamador).
      - Se usa un slice contíguo [min_idx : max_idx+1] y luego se seleccionan
        las filas exactas; esto evita fancy-indexing que zarr descomprime mal.

    Retorna array float64 de shape (len(variant_indices), n_samples).
    """
    arr   = store[dosage_key]
    lo, hi = int(variant_indices[0]), int(variant_indices[-1]) + 1

    if layout == "vs":
        raw = arr[lo:hi]                          # (hi-lo, n_samples)
        sel = variant_indices - lo                # índices locales
        chunk = raw[sel].astype(np.float64)       # (chunk, n_samples)

    elif layout == "sv":
        raw = arr[:, lo:hi]                       # (n_samples, hi-lo)
        sel = variant_indices - lo
        chunk = raw[:, sel].T.astype(np.float64)  # (chunk, n_samples)

    elif layout == "gt":
        raw = arr[lo:hi]                          # (hi-lo, n_samples, ploidy)
        sel = variant_indices - lo
        gt  = raw[sel].astype(np.float64)
        gt[gt < 0] = np.nan                       # missing = -1 en convención VCF
        chunk = gt.sum(axis=2)                    # (chunk, n_samples)

    else:
        raise ValueError(f"Layout desconocido: {layout}")

    return chunk


def _decode_allele(val) -> str:
    """
    Normaliza un valor de alelo proveniente de zarr a string en mayúsculas.
    Zarr puede devolver str nativos, numpy.str_, bytes o numpy.bytes_
    dependiendo del dtype del array (U1, S1, object, etc.).
    """
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode("ascii", errors="replace").upper()
    return str(val).upper()


# =============================================================================
# Matching de variantes
# =============================================================================

def match_variants(
    weights: pd.DataFrame,
    pos_arr: np.ndarray,
    ref_arr: np.ndarray,
    alt_arr: np.ndarray,
) -> tuple[list[int], list[int], list[bool], list[dict]]:
    """
    Empareja variantes del archivo de pesos con las del zarr por posición y alelos.

    Lógica de IS_FLIP:
      IS_FLIP=0 → efecto=ALT, otro=REF  → (zarr_REF==OTHER_ALLELE, zarr_ALT==EFFECT_ALLELE)
      IS_FLIP=1 → efecto=REF, otro=ALT  → (zarr_REF==EFFECT_ALLELE, zarr_ALT==OTHER_ALLELE)

    Retorna:
      weight_indices  — índices en `weights` que fueron emparejados
      zarr_indices    — índices correspondientes en el zarr
      flip_mask       — True si IS_FLIP=1 (hay que invertir el dosage)
      excluded        — lista de dicts con variantes excluidas y razón
    """
    # Índice posición → primera ocurrencia en zarr
    pos_to_zarr: dict[int, int] = {}
    for i, p in enumerate(pos_arr):
        if int(p) not in pos_to_zarr:
            pos_to_zarr[int(p)] = i

    weight_indices: list[int] = []
    zarr_indices:   list[int] = []
    flip_mask:      list[bool] = []
    excluded:       list[dict] = []

    for row_i, row in weights.iterrows():
        pos      = int(row["POS"])
        ea       = str(row["EFFECT_ALLELE"]).upper()
        oa       = str(row["OTHER_ALLELE"]).upper()
        is_flip  = int(row["IS_FLIP"])

        z_idx = pos_to_zarr.get(pos)
        if z_idx is None:
            excluded.append({
                "ID": row["ID"], "POS": pos,
                "reason": "posicion_no_encontrada_en_zarr",
            })
            continue

        z_ref = _decode_allele(ref_arr[z_idx])
        z_alt = _decode_allele(alt_arr[z_idx])

        # Alelos esperados según IS_FLIP
        expected_ref = oa if is_flip == 0 else ea
        expected_alt = ea if is_flip == 0 else oa

        if z_ref != expected_ref or z_alt != expected_alt:
            # Intento de strand flip (A↔T, C↔G)
            complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
            exp_ref_flip = "".join(complement.get(b, b) for b in expected_ref)
            exp_alt_flip = "".join(complement.get(b, b) for b in expected_alt)

            if z_ref == exp_ref_flip and z_alt == exp_alt_flip:
                log.debug("Strand flip detectado en POS=%d — aceptado", pos)
            else:
                excluded.append({
                    "ID": row["ID"], "POS": pos,
                    "reason": "discordancia_de_alelos",
                    "zarr_REF": z_ref, "zarr_ALT": z_alt,
                    "prs_EFFECT": ea, "prs_OTHER": oa,
                    "IS_FLIP": is_flip,
                })
                continue

        weight_indices.append(int(row_i))
        zarr_indices.append(z_idx)
        flip_mask.append(is_flip == 1)

    log.info(
        "  Emparejadas: %d / %d variantes  |  Excluidas: %d",
        len(zarr_indices), len(weights), len(excluded),
    )
    return weight_indices, zarr_indices, flip_mask, excluded


# =============================================================================
# Cómputo del PRS
# =============================================================================

def compute_prs_chunked(
    store: zarr.Group,
    dosage_key: str,
    layout: str,
    zarr_indices: np.ndarray,    # (n_matched,) — YA ORDENADOS por zarr_idx
    sort_order: np.ndarray,      # permutación para reordenar betas/flip
    flip_mask: np.ndarray,       # bool (n_matched,) — en orden original
    betas: np.ndarray,           # float64 (n_matched,) — en orden original
    n_samples: int,
    missing_strategy: str,
    chunk_size: int,
    excluded_from_missing: list,
) -> np.ndarray:
    """
    PRS aditivo acumulado en chunks para controlar el uso de memoria.

    En cada iteración se carga un bloque de `chunk_size` variantes desde zarr,
    se calcula la contribución parcial al PRS y se descarta la memoria.

    PRS_sample = Σ_i ( BETA_i × dosage_efectivo_i )
    dosage_efectivo_i:
        IS_FLIP=0 → dosage ALT  (lectura directa)
        IS_FLIP=1 → 2 − dosage ALT  (alelo de efecto es REF)
    """
    # Re-indexar betas y flip en el orden de zarr_indices
    flip_sorted  = flip_mask[sort_order]
    betas_sorted = betas[sort_order]

    prs = np.zeros(n_samples, dtype=np.float64)
    n_matched = len(zarr_indices)
    n_chunks  = (n_matched + chunk_size - 1) // chunk_size

    for chunk_i in range(n_chunks):
        lo = chunk_i * chunk_size
        hi = min(lo + chunk_size, n_matched)

        idx_chunk   = zarr_indices[lo:hi]
        flip_chunk  = flip_sorted[lo:hi]
        betas_chunk = betas_sorted[lo:hi]

        # (chunk, n_samples) — sólo este bloque en RAM
        dosage_chunk = get_dosage_chunk(store, dosage_key, layout, idx_chunk)

        for j in range(len(idx_chunk)):
            d = dosage_chunk[j]          # (n_samples,) view

            if flip_chunk[j]:
                d = 2.0 - d             # copy implícita, no modifica el chunk

            nan_mask  = np.isnan(d)
            n_missing = int(nan_mask.sum())

            if n_missing > 0:
                if missing_strategy == "mean":
                    d = d.copy()
                    d[nan_mask] = float(np.nanmean(d))
                elif missing_strategy == "zero":
                    d = d.copy()
                    d[nan_mask] = 0.0
                elif missing_strategy == "skip":
                    excluded_from_missing.append({
                        "zarr_index": int(idx_chunk[j]),
                        "n_missing_samples": n_missing,
                        "reason": "skip_por_missing",
                    })
                    continue

            prs += betas_chunk[j] * d

        if (chunk_i + 1) % 10 == 0 or chunk_i == n_chunks - 1:
            log.info(
                "  chunk %d/%d  (variantes %d–%d de %d)",
                chunk_i + 1, n_chunks, lo, hi - 1, n_matched,
            )

    return prs


# =============================================================================
# Guardar resultados
# =============================================================================

def save_results(
    pgs_id: str,
    chrom: str,
    samples: np.ndarray,
    prs: np.ndarray,
    metadata: dict,
) -> None:
    out_dir = SCORES_BASE / pgs_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Scores por muestra
    scores_path = out_dir / f"{pgs_id}_chr{chrom}_scores.tsv"
    pd.DataFrame({"sample_id": samples, "PRS": prs}).to_csv(
        scores_path, sep="\t", index=False, float_format="%.8f"
    )
    log.info("Scores guardados en: %s", scores_path)

    # Metadatos del análisis
    meta_path = out_dir / f"{pgs_id}_chr{chrom}_metadata.json"
    with open(meta_path, "w") as fh:
        json.dump(metadata, fh, indent=2, default=str)
    log.info("Metadatos guardados en: %s", meta_path)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args  = parse_args()
    pgs_id = args.pgs_id
    chrom  = args.chrom
    t_start = datetime.now()

    # ── 1. Pesos ──────────────────────────────────────────────────────────────
    weights = load_score_weights(pgs_id, chrom)
    if weights.empty:
        log.warning("Sin variantes para chr%s en %s — saltando.", chrom, pgs_id)
        sys.exit(0)

    # ── 2. Zarr ───────────────────────────────────────────────────────────────
    store   = open_zarr_store(chrom)
    samples = get_samples(store)
    pos_arr, ref_arr, alt_arr = get_variant_arrays(store)
    log.info("Zarr: %d variantes × %d muestras", len(pos_arr), len(samples))

    # Detectar key de dosage una sola vez (evita re-escanear en cada chunk)
    dosage_key, layout = _detect_dosage_key(store)
    log.info("Dosage array: '%s'  layout: %s", dosage_key, layout)

    # ── 3. Emparejamiento ─────────────────────────────────────────────────────
    w_idx, z_idx, flip_list, excluded = match_variants(
        weights, pos_arr, ref_arr, alt_arr
    )

    if not z_idx:
        log.error("Cero variantes emparejadas — abortando.")
        sys.exit(1)

    # Convertir a arrays numpy para indexación vectorizada
    w_idx_arr   = np.array(w_idx,    dtype=np.intp)
    z_idx_arr   = np.array(z_idx,    dtype=np.intp)
    flip_arr    = np.array(flip_list, dtype=bool)
    betas_arr   = weights.loc[w_idx_arr, "BETA"].to_numpy(dtype=np.float64)

    # Ordenar por zarr_index para acceso secuencial (mucho más eficiente en zarr)
    sort_order  = np.argsort(z_idx_arr)
    z_idx_sorted = z_idx_arr[sort_order]

    # ── 4. PRS por chunks ─────────────────────────────────────────────────────
    log.info(
        "Calculando PRS en chunks de %d variantes "
        "(%d muestras, estrategia missing: %s)...",
        args.chunk_size, len(samples), args.missing_strategy,
    )
    excluded_missing: list[dict] = []
    prs = compute_prs_chunked(
        store        = store,
        dosage_key   = dosage_key,
        layout       = layout,
        zarr_indices = z_idx_sorted,
        sort_order   = sort_order,
        flip_mask    = flip_arr,
        betas        = betas_arr,
        n_samples    = len(samples),
        missing_strategy = args.missing_strategy,
        chunk_size   = args.chunk_size,
        excluded_from_missing = excluded_missing,
    )

    elapsed = (datetime.now() - t_start).total_seconds()

    # ── 5. Metadatos ──────────────────────────────────────────────────────────
    metadata = {
        "pgs_id":                    pgs_id,
        "chrom":                     chrom,
        "run_timestamp":             t_start.isoformat(),
        "elapsed_seconds":           round(elapsed, 2),
        "missing_strategy":          args.missing_strategy,
        "chunk_size":                args.chunk_size,
        "dosage_array_key":          dosage_key,
        "dosage_layout":             layout,
        "n_samples":                 int(len(samples)),
        "n_variants_in_score_file":  int(len(weights)),
        "n_variants_matched":        int(len(z_idx)),
        "n_variants_excluded_allele":int(len(excluded)),
        "n_variants_excluded_missing": int(len(excluded_missing)),
        "prs_mean":                  float(np.mean(prs)),
        "prs_std":                   float(np.std(prs)),
        "prs_min":                   float(np.min(prs)),
        "prs_max":                   float(np.max(prs)),
        "excluded_allele_mismatch":  excluded,
        "excluded_missing":          excluded_missing,
    }

    # ── 6. Guardar ────────────────────────────────────────────────────────────
    save_results(pgs_id, chrom, samples, prs, metadata)
    log.info("chr%s completado en %.1f s ✓", chrom, elapsed)


if __name__ == "__main__":
    main()
