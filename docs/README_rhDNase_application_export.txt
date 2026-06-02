rhDNase real-data application export
===================================

Run ID:
rhDNase_final_20260602_055047

Source results:
/projects/F202500010HPCVLABUMINHO/cecilia/tauH_simulation/results_realdata_rhDNase_tauH/rhDNase_final_20260602_055047

Dataset:
survival::rhDNase

Application:
Covariate-adaptive truncated Kendall association for two sequential pulmonary exacerbation gap times.

Time scale:
Accumulated at-risk time.

Primary truncation:
H = 155, empirical H_quantile = 0.90.

Core outputs:
- results/realdata/rhDNase/tables/rhDNase_main_estimates_with_bootstrap.csv
- results/realdata/rhDNase/bootstrap/rhDNase_bootstrap_summary.csv
- results/realdata/rhDNase/tables/rhDNase_nuisance_library_diagnostics.csv
- results/realdata/rhDNase/tables/rhDNase_sensitivity_lambda.csv
- results/realdata/rhDNase/tables/rhDNase_sensitivity_H.csv
- figures/realdata/rhDNase_*.png

Validation:
The final consistency checker reported: No red flags detected.

Bootstrap:
2000 subject-level bootstrap resamples for the methods included in inference.
