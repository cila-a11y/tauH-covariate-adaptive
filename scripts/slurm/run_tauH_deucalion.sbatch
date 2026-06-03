#!/bin/bash
#SBATCH --job-name=tauH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=48:00:00
#SBATCH --partition=normal-x86
#SBATCH --exclusive
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --error=slurm_logs/%x_%A_%a.err

set -euo pipefail

# ------------------------------------------------------------------------------
# Single-node worker script for Deucalion.
# Submit this through submit_tauH_deucalion.sh rather than calling sbatch manually.
# ------------------------------------------------------------------------------

module purge
# Adjust this line if `module spider R` reports a versioned module name.
module load R

# Prevent nested BLAS/OpenMP oversubscription; the R script uses process-level workers.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

: "${TAUH_PROJECT_DIR:=$PWD}"
: "${TAUH_SCRIPT:=$TAUH_PROJECT_DIR/covariate_adaptive_tauH_simulation_deucalion.R}"
: "${TAUH_OUTPUT_ROOT:=$TAUH_PROJECT_DIR/results_covariate_adaptive_tauH}"
: "${TAUH_PROFILE:=paper}"
: "${TAUH_MODE:=all}"
: "${TAUH_CORES:=${SLURM_CPUS_PER_TASK:-64}}"
: "${TAUH_R_LIBS_USER:=$TAUH_PROJECT_DIR/Rlibs}"

export TAUH_PROJECT_DIR TAUH_SCRIPT TAUH_OUTPUT_ROOT TAUH_PROFILE TAUH_MODE TAUH_CORES TAUH_R_LIBS_USER
export R_LIBS_USER="$TAUH_R_LIBS_USER"
mkdir -p "$TAUH_OUTPUT_ROOT" "$R_LIBS_USER"

cd "$TAUH_PROJECT_DIR"

echo "[tauH] host=$(hostname)"
echo "[tauH] date=$(date)"
echo "[tauH] mode=$TAUH_MODE run_id=${TAUH_RUN_ID:-unset} array_task=${SLURM_ARRAY_TASK_ID:-none} cores=$TAUH_CORES"
echo "[tauH] script=$TAUH_SCRIPT"
echo "[tauH] output_root=$TAUH_OUTPUT_ROOT"
echo "[tauH] R_LIBS_USER=$R_LIBS_USER"

Rscript "$TAUH_SCRIPT"
