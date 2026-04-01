#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="${1:-}"
BASE_DIR="/mnt/cephfs/hot_nvme/pgscatalog/scores"

if [[ -z "$IDS_FILE" ]]; then
  echo "Uso: $0 pgs_ids.txt" >&2
  exit 1
fi

if [[ ! -f "$IDS_FILE" ]]; then
  echo "ERROR: no existe el archivo: $IDS_FILE" >&2
  exit 1
fi

if ! command -v pgscat >/dev/null 2>&1; then
  echo "ERROR: 'pgscat' no está en PATH" >&2
  exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: 'jq' no está instalado" >&2
  exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
  echo "ERROR: 'curl' no está instalado" >&2
  exit 1
fi

mkdir -p "$BASE_DIR"

download_json() {
  local pgs_id="$1"
  local out_json="$2"

  echo "[INFO] Descargando JSON para ${pgs_id}" >&2
  pgscat --json "$out_json" score "$pgs_id"
}

json_has_grch38_url() {
  local json_file="$1"
  jq -e -r '.ftp_harmonized_scoring_files.GRCh38.positions // empty' "$json_file" >/dev/null
}

extract_grch38_url() {
  local json_file="$1"
  jq -r '.ftp_harmonized_scoring_files.GRCh38.positions' "$json_file"
}

download_score_file() {
  local url="$1"
  local out_file="$2"

  echo "[INFO] Descargando score file: $(basename "$out_file")" >&2
  curl -fL --retry 3 --retry-delay 2 -o "$out_file" "$url"
}

while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
  pgs_id="$(echo "$raw_line" | tr -d '\r' | xargs)"

  # Saltar líneas vacías o comentarios
  [[ -z "$pgs_id" ]] && continue
  [[ "$pgs_id" =~ ^# ]] && continue

  sample_dir="${BASE_DIR}/${pgs_id}"
  json_file="${sample_dir}/${pgs_id}.json"

  mkdir -p "$sample_dir"

  echo "[INFO] Procesando ${pgs_id}" >&2

  # 1) JSON
  need_json=0
  if [[ ! -s "$json_file" ]]; then
    need_json=1
  elif ! jq empty "$json_file" >/dev/null 2>&1; then
    echo "[WARN] JSON corrupto o inválido, se volverá a descargar: $json_file" >&2
    need_json=1
  fi

  if [[ "$need_json" -eq 1 ]]; then
    tmp_json="${json_file}.tmp"
    rm -f "$tmp_json"
    download_json "$pgs_id" "$tmp_json"
    mv "$tmp_json" "$json_file"
  else
    echo "[INFO] JSON ya existe: $json_file" >&2
  fi

  # 2) URL GRCh38 positions
  if ! json_has_grch38_url "$json_file"; then
    echo "[ERROR] No se encontró ftp_harmonized_scoring_files.GRCh38.positions en $json_file" >&2
    continue
  fi

  score_url="$(extract_grch38_url "$json_file")"
  score_name="$(basename "$score_url")"
  score_file="${sample_dir}/${score_name}"

  # 3) Archivo score
  if [[ -s "$score_file" ]]; then
    echo "[INFO] Score file ya existe: $score_file" >&2
  else
    tmp_score="${score_file}.tmp"
    rm -f "$tmp_score"
    download_score_file "$score_url" "$tmp_score"
    mv "$tmp_score" "$score_file"
  fi

  echo "[OK] ${pgs_id}" >&2
  echo >&2

done < "$IDS_FILE"
