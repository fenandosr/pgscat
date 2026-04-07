#!/bin/bash
#SBATCH --job-name=run_prs
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-23
#SBATCH --output=/mnt/cephfs/orgs/home/gen/.logs/run_prs_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/gen/.logs/run_prs_%A_%a.err

# =============================================================================
# run_prs.slurm — SLURM array job: one task per chromosome (1-22, X)
#
# Enviar con:
#   sbatch --export=ALL,PGS_ID=PGS000001 run_prs.slurm
#
# O pasar múltiples PGS IDs en un loop:
#   for ID in PGS000001 PGS000002; do
#       sbatch --export=ALL,PGS_ID=$ID run_prs.slurm
#   done
# =============================================================================

# ── chromosome mapping: task array index [1..23] → label ─────────────────────
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHROM="${CHROMS[$SLURM_ARRAY_TASK_ID - 1]}"

# ── validar que PGS_ID fue pasado ─────────────────────────────────────────────
: "${PGS_ID:?[ERROR] PGS_ID no está definido. Usa: sbatch --export=ALL,PGS_ID=PGS000001 run_prs.slurm}"

# ── banner de inicio ──────────────────────────────────────────────────────────
echo "======================================================="
echo "  PRS Pipeline — chr${CHROM}"
echo "-------------------------------------------------------"
echo "  PGS_ID          : ${PGS_ID}"
echo "  CHROM           : ${CHROM}"
echo "  Array job ID    : ${SLURM_ARRAY_JOB_ID}"
echo "  Array task ID   : ${SLURM_ARRAY_TASK_ID}"
echo "  Node            : $(hostname)"
echo "  Start           : $(date)"
echo "======================================================="

# ── activar entorno Python ────────────────────────────────────────────────────
# Ajusta la ruta a tu instalación de conda/mamba
source /home/gen/micromamba/etc/profile.d/mamba.sh
micromamba activate vcfarr-hpc

# ── directorios y logs ────────────────────────────────────────────────────────
SCRIPT_DIR="/mnt/cephfs/orgs/home/gen/repos/pgscat/scripts"
LOG_DIR="/mnt/cephfs/orgs/home/gen/.logs"
mkdir -p "${LOG_DIR}"

# ── ejecutar el script de cómputo ────────────────────────────────────────────
python "${SCRIPT_DIR}/compute_prs.py" \
    --pgs-id "${PGS_ID}" \
    --chrom  "${CHROM}"

EXIT_CODE=$?

echo "-------------------------------------------------------"
echo "  chr${CHROM} finalizado con código: ${EXIT_CODE}"
echo "  Fin: $(date)"
echo "======================================================="

exit "${EXIT_CODE}"
