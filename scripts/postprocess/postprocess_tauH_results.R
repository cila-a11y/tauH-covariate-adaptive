results <- Sys.getenv("RESULTS")
if (!nzchar(results)) stop("Set RESULTS environment variable.")

tables_dir <- file.path(results, "tables")
boot_dir <- file.path(results, "bootstrap")

summary_file <- file.path(tables_dir, "monte_carlo_summary_by_method.csv")
pair_file <- file.path(tables_dir, "decision_criteria_pairwise_comparisons.csv")
boot_file <- file.path(boot_dir, "bootstrap_summary.csv")

if (!file.exists(summary_file)) stop("Missing: ", summary_file)
if (!file.exists(pair_file)) stop("Missing: ", pair_file)

sm <- read.csv(summary_file, stringsAsFactors = FALSE)
pairs <- read.csv(pair_file, stringsAsFactors = FALSE)

cat("Loaded summary rows:", nrow(sm), "\n")
cat("Loaded pairwise rows:", nrow(pairs), "\n")
cat("Columns in pairwise file:\n")
print(names(pairs))

# -------------------------------------------------------------------
# Add block to pairwise file if absent
# -------------------------------------------------------------------

if (!("block" %in% names(pairs))) {
  cat("Column 'block' not found in pairwise file. Merging block from summary table by scenario_id.\n")
  if (!("scenario_id" %in% names(pairs))) {
    stop("pairwise file has no 'block' and no 'scenario_id'; cannot recover block.")
  }
  lookup <- unique(sm[, c("scenario_id", "block")])
  pairs <- merge(pairs, lookup, by = "scenario_id", all.x = TRUE, sort = FALSE)
}

if (any(is.na(pairs$block))) {
  warning("Some pairwise rows still have missing block after merge.")
}

needed_pair_cols <- c("block", "method_new", "method_ref",
                      "abs_bias_reduction", "mse_reduction", "rel_mse")

missing_pair_cols <- setdiff(needed_pair_cols, names(pairs))
if (length(missing_pair_cols) > 0) {
  stop("Missing required pairwise columns: ", paste(missing_pair_cols, collapse = ", "))
}

# -------------------------------------------------------------------
# 1. Corrected decision summary by block
# -------------------------------------------------------------------

summ_one <- function(d) {
  data.frame(
    n_scenarios = nrow(d),
    abs_bias_reduction_mean = mean(d$abs_bias_reduction, na.rm = TRUE),
    abs_bias_reduction_median = median(d$abs_bias_reduction, na.rm = TRUE),
    abs_bias_reduction_prop_positive = mean(d$abs_bias_reduction > 0, na.rm = TRUE),
    mse_reduction_mean = mean(d$mse_reduction, na.rm = TRUE),
    mse_reduction_median = median(d$mse_reduction, na.rm = TRUE),
    mse_reduction_prop_positive = mean(d$mse_reduction > 0, na.rm = TRUE),
    rel_mse_mean = mean(d$rel_mse, na.rm = TRUE),
    rel_mse_median = median(d$rel_mse, na.rm = TRUE),
    rel_mse_prop_less_than_1 = mean(d$rel_mse < 1, na.rm = TRUE)
  )
}

key <- paste(pairs$block, pairs$method_new, pairs$method_ref, sep = "||")

decision_by_block <- do.call(rbind, lapply(split(pairs, key), function(d) {
  cbind(
    block = d$block[1],
    method_new = d$method_new[1],
    method_ref = d$method_ref[1],
    summ_one(d)
  )
}))

decision_by_block <- decision_by_block[order(decision_by_block$block,
                                             decision_by_block$method_new,
                                             decision_by_block$method_ref), ]

out1 <- file.path(tables_dir, "decision_criteria_summary_by_block_corrected.csv")
write.csv(decision_by_block, out1, row.names = FALSE)

# -------------------------------------------------------------------
# 2. Overall performance summary by block and method
# -------------------------------------------------------------------

needed_sm_cols <- c("block", "method", "abs_bias", "mse", "rmse", "bias", "mean_M_hat", "fail_rate")
missing_sm_cols <- setdiff(needed_sm_cols, names(sm))
if (length(missing_sm_cols) > 0) {
  stop("Missing required summary columns: ", paste(missing_sm_cols, collapse = ", "))
}

perf_one <- function(d) {
  data.frame(
    n_scenarios = length(unique(d$scenario_id)),
    abs_bias_mean = mean(d$abs_bias, na.rm = TRUE),
    abs_bias_median = median(d$abs_bias, na.rm = TRUE),
    mse_mean = mean(d$mse, na.rm = TRUE),
    mse_median = median(d$mse, na.rm = TRUE),
    rmse_mean = mean(d$rmse, na.rm = TRUE),
    rmse_median = median(d$rmse, na.rm = TRUE),
    bias_mean = mean(d$bias, na.rm = TRUE),
    bias_median = median(d$bias, na.rm = TRUE),
    mean_M_hat_mean = mean(d$mean_M_hat, na.rm = TRUE),
    fail_rate_mean = mean(d$fail_rate, na.rm = TRUE)
  )
}

key2 <- paste(sm$block, sm$method, sep = "||")

performance_by_block_method <- do.call(rbind, lapply(split(sm, key2), function(d) {
  cbind(
    block = d$block[1],
    method = d$method[1],
    perf_one(d)
  )
}))

performance_by_block_method <- performance_by_block_method[
  order(performance_by_block_method$block,
        performance_by_block_method$mse_mean,
        performance_by_block_method$method),
]

out2 <- file.path(tables_dir, "performance_summary_by_block_method.csv")
write.csv(performance_by_block_method, out2, row.names = FALSE)

# -------------------------------------------------------------------
# 3. Method wins by block using MSE rank
# -------------------------------------------------------------------

if (!("best_by_mse" %in% names(sm))) {
  warning("Column best_by_mse not found; computing from rank_mse == 1.")
  if (!("rank_mse" %in% names(sm))) stop("No best_by_mse or rank_mse available.")
  sm$best_by_mse <- sm$rank_mse == 1
}

wins_one <- function(d) {
  data.frame(
    n_best_by_mse = sum(d$best_by_mse == TRUE, na.rm = TRUE),
    n_rows = nrow(d)
  )
}

key3 <- paste(sm$block, sm$method, sep = "||")

method_wins <- do.call(rbind, lapply(split(sm, key3), function(d) {
  z <- wins_one(d)
  cbind(
    block = d$block[1],
    method = d$method[1],
    z
  )
}))

n_scen_by_block <- aggregate(scenario_id ~ block, data = unique(sm[, c("block", "scenario_id")]), FUN = length)
names(n_scen_by_block)[names(n_scen_by_block) == "scenario_id"] <- "n_scenarios"

method_wins <- merge(method_wins, n_scen_by_block, by = "block", all.x = TRUE, sort = FALSE)
method_wins$win_rate <- method_wins$n_best_by_mse / method_wins$n_scenarios

method_wins <- method_wins[order(method_wins$block, -method_wins$win_rate, method_wins$method), ]

out3 <- file.path(tables_dir, "method_wins_by_block_corrected.csv")
write.csv(method_wins, out3, row.names = FALSE)

# -------------------------------------------------------------------
# 4. Bootstrap compact summary
# -------------------------------------------------------------------

if (file.exists(boot_file)) {
  boot <- read.csv(boot_file, stringsAsFactors = FALSE)

  needed_boot_cols <- c("block", "method", "coverage", "mean_width", "bias", "mse", "mean_n_boot_valid")
  missing_boot_cols <- setdiff(needed_boot_cols, names(boot))

  if (length(missing_boot_cols) == 0) {
    boot_one <- function(d) {
      data.frame(
        n_scenarios = length(unique(d$scenario_id)),
        coverage_mean = mean(d$coverage, na.rm = TRUE),
        coverage_median = median(d$coverage, na.rm = TRUE),
        mean_width_mean = mean(d$mean_width, na.rm = TRUE),
        mean_width_median = median(d$mean_width, na.rm = TRUE),
        bias_mean = mean(d$bias, na.rm = TRUE),
        bias_median = median(d$bias, na.rm = TRUE),
        mse_mean = mean(d$mse, na.rm = TRUE),
        mse_median = median(d$mse, na.rm = TRUE),
        mean_n_boot_valid_mean = mean(d$mean_n_boot_valid, na.rm = TRUE)
      )
    }

    key4 <- paste(boot$block, boot$method, sep = "||")

    bootstrap_by_block_method <- do.call(rbind, lapply(split(boot, key4), function(d) {
      cbind(
        block = d$block[1],
        method = d$method[1],
        boot_one(d)
      )
    }))

    bootstrap_by_block_method <- bootstrap_by_block_method[
      order(bootstrap_by_block_method$block,
            bootstrap_by_block_method$mse_mean,
            bootstrap_by_block_method$method),
    ]

    out4 <- file.path(boot_dir, "bootstrap_summary_by_block_method.csv")
    write.csv(bootstrap_by_block_method, out4, row.names = FALSE)
  } else {
    warning("Bootstrap summary missing columns: ", paste(missing_boot_cols, collapse = ", "))
  }
}

cat("Post-processing complete.\n")
cat("Created:\n")
cat(out1, "\n")
cat(out2, "\n")
cat(out3, "\n")
if (exists("out4")) cat(out4, "\n")
