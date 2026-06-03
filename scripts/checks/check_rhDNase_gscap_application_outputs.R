args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript check_rhDNase_gscap_application_outputs.R <RESULTS_DIR>")
}

results <- normalizePath(args[1], mustWork = TRUE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

main <- read_csv(file.path(results, "tables", "rhDNase_main_point_estimates.csv"))
main_boot <- read_csv(file.path(results, "tables", "rhDNase_main_estimates_with_bootstrap.csv"))
weights <- read_csv(file.path(results, "tables", "rhDNase_weight_diagnostics_main_H.csv"))
sensH <- read_csv(file.path(results, "tables", "rhDNase_sensitivity_H.csv"))
boot <- read_csv(file.path(results, "bootstrap", "rhDNase_bootstrap_summary.csv"))

expected_main <- c(
  "trad",
  "ipcw_marg",
  "ipcw_x",
  "gscap_x_l050",
  "gscap_x_l100",
  "ps_marg_logit",
  "ps_x_logit",
  "ps_x_spline",
  "scap_x_l050",
  "scap_x_l075",
  "scap_x_l100"
)

expected_boot <- setdiff(expected_main, "scap_x_l100")

cat("\n=== Basic files loaded ===\n")
cat("main rows:", nrow(main), "\n")
cat("main_boot rows:", nrow(main_boot), "\n")
cat("weights rows:", nrow(weights), "\n")
cat("sensitivity H rows:", nrow(sensH), "\n")
cat("bootstrap rows:", nrow(boot), "\n")

cat("\n=== Main methods ===\n")
print(sort(main$method))

cat("\n=== Bootstrap methods ===\n")
print(sort(boot$method))

red <- character(0)

missing_main <- setdiff(expected_main, main$method)
extra_main <- setdiff(main$method, expected_main)

missing_boot <- setdiff(expected_boot, boot$method)
extra_boot <- setdiff(boot$method, expected_boot)

if (nrow(main) != length(expected_main)) {
  red <- c(red, paste0("Expected ", length(expected_main), " main rows, found ", nrow(main), "."))
}
if (nrow(main_boot) != length(expected_main)) {
  red <- c(red, paste0("Expected ", length(expected_main), " main_boot rows, found ", nrow(main_boot), "."))
}
if (nrow(weights) != length(expected_main)) {
  red <- c(red, paste0("Expected ", length(expected_main), " weight rows, found ", nrow(weights), "."))
}
if (nrow(sensH) != length(expected_main) * 4L) {
  red <- c(red, paste0("Expected ", length(expected_main) * 4L, " H-sensitivity rows, found ", nrow(sensH), "."))
}
if (nrow(boot) != length(expected_boot)) {
  red <- c(red, paste0("Expected ", length(expected_boot), " bootstrap rows, found ", nrow(boot), "."))
}

if (length(missing_main)) red <- c(red, paste("Missing main methods:", paste(missing_main, collapse = ", ")))
if (length(extra_main)) red <- c(red, paste("Unexpected main methods:", paste(extra_main, collapse = ", ")))
if (length(missing_boot)) red <- c(red, paste("Missing bootstrap methods:", paste(missing_boot, collapse = ", ")))
if (length(extra_boot)) red <- c(red, paste("Unexpected bootstrap methods:", paste(extra_boot, collapse = ", ")))

if (!all(c("gscap_x_l050", "gscap_x_l100") %in% main$method)) {
  red <- c(red, "G-sCAP methods missing from main estimates.")
}
if (!all(c("gscap_x_l050", "gscap_x_l100") %in% boot$method)) {
  red <- c(red, "G-sCAP methods missing from bootstrap summary.")
}
if (!all(c("gscap_x_l050", "gscap_x_l100") %in% sensH$method)) {
  red <- c(red, "G-sCAP methods missing from H sensitivity.")
}

cat("\n=== Main estimates ===\n")
print(main[, intersect(c("method", "tau_H", "M_hat", "ess", "H", "H_quantile"), names(main))])

cat("\n=== Main vs bootstrap estimate consistency ===\n")
mb <- merge(
  main[, c("method", "tau_H")],
  main_boot[, c("method", "estimate", "n_boot_valid")],
  by = "method",
  all.x = TRUE
)
mb$diff_tau_estimate <- mb$tau_H - mb$estimate
print(mb)

bad_mb <- mb[is.finite(mb$estimate) & abs(mb$diff_tau_estimate) > 1e-10, , drop = FALSE]
if (nrow(bad_mb)) {
  red <- c(red, paste0("Main and bootstrap point estimates disagree. Max abs diff = ",
                       signif(max(abs(bad_mb$diff_tau_estimate), na.rm = TRUE), 6)))
}

cat("\n=== Main vs weight diagnostics consistency ===\n")
mw <- merge(
  main[, c("method", "M_hat", "ess", "sum_weights")],
  weights[, c("method", "M_hat", "ess", "sum_weights")],
  by = "method",
  suffixes = c("_main", "_weights"),
  all.x = TRUE
)
mw$diff_M_hat <- mw$M_hat_main - mw$M_hat_weights
mw$diff_ess <- mw$ess_main - mw$ess_weights
mw$diff_sum_weights <- mw$sum_weights_main - mw$sum_weights_weights
print(mw)

if (max(abs(mw$diff_M_hat), na.rm = TRUE) > 1e-10 ||
    max(abs(mw$diff_ess), na.rm = TRUE) > 1e-10 ||
    max(abs(mw$diff_sum_weights), na.rm = TRUE) > 1e-10) {
  red <- c(red, "Main estimates and weight diagnostics disagree.")
}

cat("\n=== Main H and sensitivity-H consistency at H_quantile=0.90 ===\n")
s90 <- sensH[abs(sensH$H_quantile - 0.90) < 1e-12, , drop = FALSE]
mh <- merge(
  main[, c("method", "tau_H", "M_hat", "ess")],
  s90[, c("method", "tau_H", "M_hat", "ess")],
  by = "method",
  suffixes = c("_main", "_sensH"),
  all.x = TRUE
)
mh$diff_tau <- mh$tau_H_main - mh$tau_H_sensH
mh$diff_M_hat <- mh$M_hat_main - mh$M_hat_sensH
mh$diff_ess <- mh$ess_main - mh$ess_sensH
print(mh)

if (max(abs(mh$diff_tau), na.rm = TRUE) > 1e-10 ||
    max(abs(mh$diff_M_hat), na.rm = TRUE) > 1e-10 ||
    max(abs(mh$diff_ess), na.rm = TRUE) > 1e-10) {
  red <- c(red, paste0("Main estimates and H-sensitivity do not agree at H_quantile=0.90. Max abs tau diff = ",
                       signif(max(abs(mh$diff_tau), na.rm = TRUE), 6)))
}

cat("\n=== Bootstrap validity ===\n")
print(boot[, intersect(c("method", "estimate", "boot_se", "ci_perc_l", "ci_perc_u",
                         "ci_width_perc", "n_boot_valid"), names(boot))])

if (any(!is.finite(boot$boot_se))) {
  red <- c(red, "Some bootstrap standard errors are not finite.")
}
if (any(boot$n_boot_valid <= 0)) {
  red <- c(red, "Some bootstrap methods have zero valid resamples.")
}

cat("\n=== Red flags ===\n")
if (length(red)) {
  for (r in red) cat("- ", r, "\n", sep = "")
  quit(status = 1)
} else {
  cat("No red flags detected.\n")
}
