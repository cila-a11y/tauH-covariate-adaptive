args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript check_kidney_application_outputs.R <RESULTS_DIR>")
}

results <- args[1]
tables <- file.path(results, "tables")
boot <- file.path(results, "bootstrap")

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE)
}

main <- read_csv(file.path(tables, "kidney_main_point_estimates.csv"))
main_boot <- read_csv(file.path(tables, "kidney_main_estimates_with_bootstrap.csv"))
weights <- read_csv(file.path(tables, "kidney_weight_diagnostics_main_H.csv"))
sens_H <- read_csv(file.path(tables, "kidney_sensitivity_H.csv"))
sens_lam <- read_csv(file.path(tables, "kidney_sensitivity_lambda.csv"))
boot_sum <- read_csv(file.path(boot, "kidney_bootstrap_summary.csv"))

cat("\n=== Basic files loaded ===\n")
cat("main rows:", nrow(main), "\n")
cat("main_boot rows:", nrow(main_boot), "\n")
cat("weights rows:", nrow(weights), "\n")
cat("sensitivity H rows:", nrow(sens_H), "\n")
cat("sensitivity lambda rows:", nrow(sens_lam), "\n")
cat("bootstrap rows:", nrow(boot_sum), "\n")

cat("\n=== Main estimates ===\n")
print(main[, c("method", "tau_H", "M_hat", "ess", "H", "H_quantile")])

cat("\n=== Main vs bootstrap estimate consistency ===\n")
mb <- merge(
  main[, c("method", "tau_H")],
  main_boot[, c("method", "estimate", "n_boot_valid")],
  by = "method",
  all.x = TRUE
)
mb$diff_tau_estimate <- mb$tau_H - mb$estimate
print(mb)

cat("\n=== Main vs weight diagnostics consistency ===\n")
mw <- merge(
  main[, c("method", "M_hat", "ess", "sum_weights")],
  weights[, c("method", "M_hat", "ess", "sum_weights")],
  by = "method",
  all.x = TRUE,
  suffixes = c("_main", "_weights")
)
mw$diff_M_hat <- mw$M_hat_main - mw$M_hat_weights
mw$diff_ess <- mw$ess_main - mw$ess_weights
mw$diff_sum_weights <- mw$sum_weights_main - mw$sum_weights_weights
print(mw)

cat("\n=== Main H and sensitivity-H consistency at H_quantile=0.90 ===\n")
sh90 <- sens_H[abs(sens_H$H_quantile - 0.90) < 1e-12, ]
mh <- merge(
  main[, c("method", "tau_H", "M_hat", "ess")],
  sh90[, c("method", "tau_H", "M_hat", "ess")],
  by = "method",
  all = TRUE,
  suffixes = c("_main", "_sensH")
)
mh$diff_tau <- mh$tau_H_main - mh$tau_H_sensH
mh$diff_M_hat <- mh$M_hat_main - mh$M_hat_sensH
mh$diff_ess <- mh$ess_main - mh$ess_sensH
print(mh[order(abs(mh$diff_tau), decreasing = TRUE), ])

cat("\n=== Main sCAP vs lambda sensitivity consistency ===\n")
lambda_map <- data.frame(
  method = c("trad", "scap_x_l050", "scap_x_l075", "scap_x_l100"),
  lambda = c(0, 0.5, 0.75, 1.0),
  stringsAsFactors = FALSE
)

sl <- merge(lambda_map, sens_lam, by = "lambda", all.x = TRUE)
ml <- merge(
  main[, c("method", "tau_H", "M_hat", "ess")],
  sl[, c("method", "lambda", "tau_H", "M_hat", "ess")],
  by = "method",
  all = TRUE,
  suffixes = c("_main", "_sensLambda")
)
ml$diff_tau <- ml$tau_H_main - ml$tau_H_sensLambda
ml$diff_M_hat <- ml$M_hat_main - ml$M_hat_sensLambda
ml$diff_ess <- ml$ess_main - ml$ess_sensLambda
print(ml[order(abs(ml$diff_tau), decreasing = TRUE), ])

cat("\n=== Bootstrap validity ===\n")
print(boot_sum[, c("method", "estimate", "boot_se", "ci_perc_l", "ci_perc_u", "ci_width_perc", "n_boot_valid")])

cat("\n=== Red flags ===\n")
red <- list()

if (any(!is.finite(mb$diff_tau_estimate[!is.na(mb$estimate)]))) {
  red[[length(red) + 1L]] <- "Non-finite main vs bootstrap estimate differences."
}

if (max(abs(mh$diff_tau), na.rm = TRUE) > 1e-8) {
  red[[length(red) + 1L]] <- paste0(
    "Main estimates and H-sensitivity do not agree at H_quantile=0.90. Max abs diff = ",
    signif(max(abs(mh$diff_tau), na.rm = TRUE), 6)
  )
}

if (max(abs(ml$diff_tau), na.rm = TRUE) > 1e-8) {
  red[[length(red) + 1L]] <- paste0(
    "Main sCAP estimates and lambda-sensitivity do not agree. Max abs diff = ",
    signif(max(abs(ml$diff_tau), na.rm = TRUE), 6)
  )
}

if (any(boot_sum$n_boot_valid < 100, na.rm = TRUE)) {
  red[[length(red) + 1L]] <- "Some methods have fewer than 100 valid bootstrap replicates."
}

if (length(red) == 0L) {
  cat("No red flags detected.\n")
} else {
  for (r in red) cat("- ", r, "\n", sep = "")
}

