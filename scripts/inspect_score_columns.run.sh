#!/usr/bin/env bash
set -euo pipefail

cd /mnt/cephfs/orgs/home/gen/repos/pgscat

find /mnt/cephfs/hot_nvme/pgscatalog/scores \
  -type f \
  -name '*_hmPOS_GRCh38.txt.gz' \
| sort \
| while read -r f; do
    echo "[INFO] Processing $f"
    uv run python scripts/inspect_score_columns.py \
        ${f} --genome_build GRCh38 --output ${f%.txt.gz}.columns.json
    echo "[INFO] Done processing $f"
  done
