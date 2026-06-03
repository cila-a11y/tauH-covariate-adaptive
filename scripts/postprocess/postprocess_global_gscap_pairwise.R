args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript postprocess_global_gscap_pairwise.R <RESULTS_DIR>")
}

results <- normalizePath(args[1], mustWork = TRUE)
tables_dir <- file.path(results, "tables")

infile <- file.path(tables_dir, "monte_carlo_summary_by_method.csv")
if (!file.exists(infile)) {
  stop("Missing file: ", infile)
}

sm <- read.csv(infile, stringsAsFactors = FALSE)

needed <- c("block", "scenario_id", "method", "abs_bias", "mse")
miss <- setdiff(needed, names(sm))
if (length(miss)) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

pairs <- data.frame(
  method_new = c(
    "gscap_x_l050",
    "gscap_x_l100",
    "gscap_x_l050",
    "gscap_x_l100"
  ),
  method_ref = c(
    "ipcw_x",
    "ipcw_x",
    "scap_x_l050",
    "scap_x_l050"
  ),
  comparison_label = c(
    "G-sCAP(0.50) vs IPCW-X",
    "G-sCAP(1.00) vs IPCW-X",
    "G-sCAP(0.50) vs PL-sCAP(0.50)",
    "G-sCAP(1.00) vs PL-sCAP(0.50)"
  ),
  stringsAsFactors = FALSE
)

block_order <- c(
  "A_independent",
  "B_covariate_dependent",
  "C_nonlinear_observability"
)

one_pair <- function(method_new, method_ref, comparison_label) {
  a <- sm[sm$method == method_new,
          c("block", "scenario_id", "abs_bias", "mse"),
          drop = FALSE]
  b <- sm[sm$method == method_ref,
          c("block", "scenario_id", "abs_bias", "mse"),
          drop = FALSE]

  if (!nrow(a)) stop("No rows for method_new = ", method_new)
  if (!nrow(b)) stop("No rows for method_ref = ", method_ref)

  names(a)[names(a) == "abs_bias"] <- "abs_bias_new"
  names(a)[names(a) == "mse"] <- "mse_new"

  names(b)[names(b) == "abs_bias"] <- "abs_bias_ref"
  names(b)[names(b) == "mse"] <- "mse_ref"

  d <- merge(a, b, by = c("block", "scenario_id"))

  if (!nrow(d)) {
    stop("No merged rows for pair: ", method_new, " vs ", method_ref)
  }

  d$method_new <- method_new
  d$method_ref <- method_ref
  d$comparison <- comparison_label

  d$abs_bias_reduction <- d$abs_bias_ref - d$abs_bias_new
  d$mse_reduction <- d$mse_ref - d$mse_new
  d$rel_mse <- d$mse_new / d$mse_ref

  scenario_out <- d[, c(
    "block",
    "scenario_id",
    "comparison",
    "method_new",
    "method_ref",
    "abs_bias_new",
    "abs_bias_ref",
    "abs_bias_reduction",
    "mse_new",
    "mse_ref",
    "mse_reduction",
    "rel_mse"
  )]

  block_out <- do.call(rbind, lapply(split(d, d$block), function(x) {
    data.frame(
      block = x$block[1],
      comparison = comparison_label,
      method_new = method_new,
      method_ref = method_ref,
      n_scenarios = length(unique(x$scenario_id)),

      abs_bias_reduction_mean = mean(x$abs_bias_reduction, na.rm = TRUE),
      abs_bias_reduction_median = median(x$abs_bias_reduction, na.rm = TRUE),
      abs_bias_reduction_prop_positive = mean(x$abs_bias_reduction > 0, na.rm = TRUE),

      mse_reduction_mean = mean(x$mse_reduction, na.rm = TRUE),
      mse_reduction_median = median(x$mse_reduction, na.rm = TRUE),
      mse_reduction_prop_positive = mean(x$mse_reduction > 0, na.rm = TRUE),

      rel_mse_mean = mean(x$rel_mse, na.rm = TRUE),
      rel_mse_median = median(x$rel_mse, na.rm = TRUE),
      rel_mse_prop_less_than_1 = mean(x$rel_mse < 1, na.rm = TRUE),

      stringsAsFactors = FALSE
    )
  }))

  list(scenario = scenario_out, block = block_out)
}

ans <- lapply(seq_len(nrow(pairs)), function(i) {
  one_pair(
    method_new = pairs$method_new[i],
    method_ref = pairs$method_ref[i],
    comparison_label = pairs$comparison_label[i]
  )
})

scenario_tab <- do.call(rbind, lapply(ans, `[[`, "scenario"))
block_tab <- do.call(rbind, lapply(ans, `[[`, "block"))

block_tab$block_order <- match(block_tab$block, block_order)
block_tab$pair_order <- match(block_tab$comparison, pairs$comparison_label)
block_tab <- block_tab[order(block_tab$block_order, block_tab$pair_order), ]
block_tab$block_order <- NULL
block_tab$pair_order <- NULL

scenario_tab$block_order <- match(scenario_tab$block, block_order)
scenario_tab$pair_order <- match(scenario_tab$comparison, pairs$comparison_label)
scenario_tab <- scenario_tab[order(scenario_tab$block_order,
                                   scenario_tab$pair_order,
                                   scenario_tab$scenario_id), ]
scenario_tab$block_order <- NULL
scenario_tab$pair_order <- NULL

scenario_file <- file.path(tables_dir, "decision_criteria_pairwise_comparisons_gscap_global.csv")
block_file <- file.path(tables_dir, "decision_criteria_summary_gscap_global.csv")

write.csv(scenario_tab, scenario_file, row.names = FALSE)
write.csv(block_tab, block_file, row.names = FALSE)

cat("Created:\n")
cat(scenario_file, "\n")
cat(block_file, "\n\n")

cat("Block-level global G-sCAP pairwise summary:\n")
print(block_tab)

cat("\nScenario counts by block and comparison:\n")
print(with(block_tab, tapply(n_scenarios, list(block, comparison), unique)))

expected <- c(
  A_independent = 48,
  B_covariate_dependent = 96,
  C_nonlinear_observability = 48
)

bad <- merge(
  block_tab[, c("block", "comparison", "n_scenarios")],
  data.frame(block = names(expected), expected = as.integer(expected)),
  by = "block",
  all.x = TRUE
)

bad <- bad[bad$n_scenarios != bad$expected, , drop = FALSE]

if (nrow(bad)) {
  cat("\nUnexpected scenario counts:\n")
  print(bad)
  stop("Some block-level comparison counts are not 48/96/48.")
}

cat("\nAll block-level comparison counts are correct: 48/96/48.\n")
