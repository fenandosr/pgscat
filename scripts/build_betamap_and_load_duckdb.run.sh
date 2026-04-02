cd /mnt/cephfs/orgs/home/gen/repos/pgscat

find /mnt/cephfs/hot_nvme/pgscatalog/scores \
  -type f \
  -name '*_hmPOS_GRCh38.txt.gz' \
| sort \
| while read -r f; do
    echo "[INFO] Processing $f"
    uv run python scripts/build_betamap_and_load_duckdb.py \
      "$f" \
      --duckdb /mnt/cephfs/hot_nvme/mcps/pgscatalog/prsvariants.duckdb
  done