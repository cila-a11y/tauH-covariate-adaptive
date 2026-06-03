#!/bin/bash
#SBATCH --job-name=tauH_fboot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --partition=normal-x86
#SBATCH --exclusive
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --error=slurm_logs/%x_%A_%a.err

set -euo pipefail

module purge
module load R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

: "${TAUH_PROJECT_DIR:=$PWD}"
: "${TAUH_FOCUSED_SCRIPT:=$TAUH_PROJECT_DIR/focused_bootstrap_tauH_deucalion.R}"
: "${TAUH_MAIN_SCRIPT:=$TAUH_PROJECT_DIR/covariate_adaptive_tauH_simulation_deucalion.R}"
: "${TAUH_OUTPUT_ROOT:=$TAUH_PROJECT_DIR/results_focused_bootstrap_tauH}"
: "${TAUH_FOCUSED_MODE:=bootstrap}"
: "${TAUH_CORES:=${SLURM_CPUS_PER_TASK:-64}}"
: "${TAUH_R_LIBS_USER:=$TAUH_PROJECT_DIR/Rlibs}"

export TAUH_PROJECT_DIR TAUH_FOCUSED_SCRIPT TAUH_MAIN_SCRIPT TAUH_OUTPUT_ROOT TAUH_FOCUSED_MODE TAUH_CORES TAUH_R_LIBS_USER
export R_LIBS_USER="$TAUH_R_LIBS_USER"
mkdir -p "$TAUH_OUTPUT_ROOT" "$R_LIBS_USER"
cd "$TAUH_PROJECT_DIR"

echo "[tauH-focused] host=$(hostname)"
echo "[tauH-focused] date=$(date)"
echo "[tauH-focused] mode=$TAUH_FOCUSED_MODE run_id=${TAUH_RUN_ID:-unset} array_task=${SLURM_ARRAY_TASK_ID:-none} cores=$TAUH_CORES"
echo "[tauH-focused] focused_script=$TAUH_FOCUSED_SCRIPT"
echo "[tauH-focused] main_script=$TAUH_MAIN_SCRIPT"
echo "[tauH-focused] output_root=$TAUH_OUTPUT_ROOT"
echo "[tauH-focused] source_results=${TAUH_SOURCE_RESULTS:-unset}"
echo "[tauH-focused] R_LIBS_USER=$R_LIBS_USER"

Rscript "$TAUH_FOCUSED_SCRIPT"
