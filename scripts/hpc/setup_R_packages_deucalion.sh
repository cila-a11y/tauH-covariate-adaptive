#!/bin/bash
set -euo pipefail

# Run once, preferably in a suitable interactive/job environment with internet access.
# Adjust `module load R` if Deucalion reports a versioned module name via `module spider R`.

PROJECT_DIR="${TAUH_PROJECT_DIR:-$PWD}"
R_LIBS_USER="${TAUH_R_LIBS_USER:-$PROJECT_DIR/Rlibs}"
mkdir -p "$R_LIBS_USER"
export R_LIBS_USER

module purge
module load R

Rscript - <<'RS'
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("ranger", "gbm", "ggplot2")
installed <- rownames(installed.packages(lib.loc = lib))
missing <- setdiff(pkgs, installed)
if (length(missing)) {
  install.packages(missing, lib = lib, dependencies = TRUE)
}
cat("R library:", lib, "\n")
cat("Installed/available packages checked:", paste(pkgs, collapse = ", "), "\n")
RS
