#!/bin/bash
#SBATCH --job-name=tauH_kidney_real
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/tauH_kidney_real_%j.out
#SBATCH --error=slurm_logs/tauH_kidney_real_%j.err

set -euo pipefail

PROJECT_DIR=${TAUH_PROJECT_DIR:-$PWD}
SCRIPT=${TAUH_REAL_SCRIPT:-$PROJECT_DIR/realdata_kidney_tauH_application.R}
export R_LIBS_USER=${TAUH_R_LIBS_USER:-$PROJECT_DIR/Rlibs}
export TAUH_REAL_CORES=${TAUH_REAL_CORES:-${SLURM_CPUS_PER_TASK:-64}}

cd "$PROJECT_DIR"
mkdir -p slurm_logs "$R_LIBS_USER"

module purge
module load R

echo "[kidney-real] host=$(hostname)"
echo "[kidney-real] date=$(date)"
echo "[kidney-real] script=$SCRIPT"
echo "[kidney-real] R_LIBS_USER=$R_LIBS_USER"
echo "[kidney-real] cores=$TAUH_REAL_CORES"
echo "[kidney-real] n_boot=${TAUH_REAL_N_BOOT:-2000}"

Rscript "$SCRIPT"
