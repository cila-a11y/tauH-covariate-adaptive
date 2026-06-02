#!/bin/bash
set -euo pipefail

cd "$TAUH_PROJECT_DIR"

module purge
module load R

export R_LIBS_USER="$TAUH_PROJECT_DIR/Rlibs"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "[rhDNase-real] host=$(hostname)"
echo "[rhDNase-real] date=$(date)"
echo "[rhDNase-real] project_dir=$TAUH_PROJECT_DIR"
echo "[rhDNase-real] script=$TAUH_PROJECT_DIR/realdata_rhDNase_tauH_application.R"
echo "[rhDNase-real] output_root=$TAUH_REAL_OUTPUT_ROOT"
echo "[rhDNase-real] run_id=$TAUH_REAL_RUN_ID"
echo "[rhDNase-real] R_LIBS_USER=$R_LIBS_USER"
echo "[rhDNase-real] cores=$TAUH_REAL_CORES"
echo "[rhDNase-real] n_boot=$TAUH_REAL_N_BOOT"

Rscript "$TAUH_PROJECT_DIR/realdata_rhDNase_tauH_application.R"
