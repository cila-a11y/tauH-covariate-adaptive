# Covariate-adaptive estimation of truncated Kendall association for censored sequential gap times

This archive contains scripts, data, selected numerical results, tables, figures, and diagnostics for the tau_H project.

## Final runs

- Simulation run: paper_blockC_20260531_194657
- Focused bootstrap run: focused_bootstrap_20260531_211047
- Real-data run: kidney_real_fixed_20260601_045254

## Main folders

- scripts/: R scripts and Deucalion submission scripts.
- data/: public/derived real-data inputs used in the kidney application.
- results/simulation/: summarized Monte Carlo results and simulation figures.
- results/focused_bootstrap/: focused bootstrap outputs, when available.
- results/realdata/: kidney application tables, bootstrap summaries, figures, LaTeX tables, and checks.
- manuscript_figures/: figures ready to be copied or referenced from the LaTeX manuscript.

## Large raw outputs

The full raw Monte Carlo outputs are not included in this small archive. They are packaged separately in the large raw-results archive.

Created on: Mon  1 Jun 05:19:22 WEST 2026

## Real-data application: rhDNase

The real-data application uses the public `rhDNase` dataset from the R package `survival`. The analysis constructs two sequential pulmonary-exacerbation gap times on accumulated at-risk time and estimates the normalized truncated Kendall association at the primary truncation level \(H=155\).

Main files:

- `scripts/realdata/realdata_rhDNase_tauH_application.R`
- `scripts/slurm/run_realdata_rhDNase_tauH.sbatch`
- `scripts/slurm/submit_realdata_rhDNase_tauH.sh`
- `scripts/checks/check_rhDNase_application_outputs.R`
- `results/realdata/rhDNase/tables/rhDNase_main_estimates_with_bootstrap.csv`
- `results/realdata/rhDNase/bootstrap/rhDNase_bootstrap_summary.csv`
- `figures/realdata/rhDNase_tauH_estimates_bootstrap_CI.png`

The final archived run is `rhDNase_final_20260602_055047`, based on 2000 subject-level bootstrap resamples. The consistency checker reported no red flags.
