for d in /mnt/cephfs/hot_nvme/pgscatalog/scores/*/; do
    pgs=$(basename "$d")
    sbatch --export=ALL,PGS_ID="$pgs" run_prs.slurm
done
