args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript postprocess_gscap_final_outputs.R <RESULTS_DIR>")
}

results <- normalizePath(args[1], mustWork = TRUE)
tables_dir <- file.path(results, "tables")
boot_dir <- file.path(results, "bootstrap")

sm_file <- file.path(tables_dir, "monte_carlo_summary_by_method.csv")
boot_file <- file.path(boot_dir, "bootstrap_summary.csv")

if (!file.exists(sm_file)) stop("Missing: ", sm_file)

sm <- read.csv(sm_file, stringsAsFactors = FALSE)

need <- c("block", "scenario_id", "method", "abs_bias", "mse")
missing <- setdiff(need, names(sm))
if (length(missing)) stop("Missing columns in summary: ", paste(missing, collapse = ", "))

# ---------------------------------------------------------------------
# 1. Unique winner summary, with deterministic tie-breaking.
# ---------------------------------------------------------------------
method_order <- c(
  "trad",
  "ipcw_marg",
  "ipcw_x",
  "gscap_x_l050",
  "gscap_x_l100",
  "ps_marg_logit",
  "ps_x_logit",
  "ps_x_spline",
  "scap_marg_l050",
  "scap_x_l025",
  "scap_x_l050",
  "scap_x_l075",
  "scap_x_l100",
  "scap_x_nocf_l050",
  "scap_x_uncal_l050"
)

sm$method_order <- match(sm$method, method_order)
sm$method_order[is.na(sm$method_order)] <- 9999L

winner_rows <- do.call(rbind, lapply(split(sm, interaction(sm$block, sm$scenario_id, drop = TRUE)), function(d) {
  d <- d[is.finite(d$mse), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$mse, d$method_order, d$method), , drop = FALSE]
  d[1, c("block", "scenario_id", "method", "mse", "abs_bias"), drop = FALSE]
}))

win_tab <- as.data.frame(table(winner_rows$block, winner_rows$method), stringsAsFactors = FALSE)
names(win_tab) <- c("block", "method", "n_wins")
win_tab <- win_tab[win_tab$n_wins > 0, , drop = FALSE]

n_scen <- aggregate(scenario_id ~ block, winner_rows, function(x) length(unique(x)))
names(n_scen)[2] <- "n_scenarios"

win_tab <- merge(win_tab, n_scen, by = "block", all.x = TRUE)
win_tab$win_rate <- win_tab$n_wins / win_tab$n_scenarios
win_tab <- win_tab[order(win_tab$block, match(win_tab$method, method_order)), , drop = FALSE]

write.csv(
  win_tab,
  file.path(tables_dir, "method_win_summary_unique_with_gscap.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 2. Pairwise decision summaries involving G-sCAP.
# ---------------------------------------------------------------------
pairs <- data.frame(
  method_new = c(
    "gscap_x_l050",
    "gscap_x_l100",
    "gscap_x_l050",
    "gscap_x_l100",
    "gscap_x_l050",
    "gscap_x_l100",
    "gscap_x_l050",
    "gscap_x_l100"
  ),
  method_ref = c(
    "ipcw_x",
    "ipcw_x",
    "ipcw_marg",
    "ipcw_marg",
    "trad",
    "trad",
    "scap_x_l050",
    "scap_x_l100"
  ),
  stringsAsFactors = FALSE
)

one_pair <- function(method_new, method_ref) {
  a <- sm[sm$method == method_new, c("block", "scenario_id", "abs_bias", "mse")]
  b <- sm[sm$method == method_ref, c("block", "scenario_id", "abs_bias", "mse")]

  if (!nrow(a) || !nrow(b)) return(NULL)

  names(a)[names(a) == "abs_bias"] <- "abs_bias_new"
  names(a)[names(a) == "mse"] <- "mse_new"
  names(b)[names(b) == "abs_bias"] <- "abs_bias_ref"
  names(b)[names(b) == "mse"] <- "mse_ref"

  d <- merge(a, b, by = c("block", "scenario_id"))
  if (!nrow(d)) return(NULL)

  d$abs_bias_reduction <- d$abs_bias_ref - d$abs_bias_new
  d$mse_reduction <- d$mse_ref - d$mse_new
  d$rel_mse <- d$mse_new / d$mse_ref

  out <- do.call(rbind, lapply(split(d, d$block), function(x) {
    data.frame(
      block = x$block[1],
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

  out
}

pair_tab <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  one_pair(pairs$method_new[i], pairs$method_ref[i])
}))

write.csv(
  pair_tab,
  file.path(tables_dir, "decision_criteria_summary_gscap.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 3. Bootstrap summaries for G-sCAP, by block and method.
# ---------------------------------------------------------------------
if (file.exists(boot_file)) {
  boot <- read.csv(boot_file, stringsAsFactors = FALSE)

  if (all(c("block", "method", "coverage", "mean_width", "bias", "mse", "mean_n_boot_valid") %in% names(boot))) {
    boot_block <- do.call(rbind, lapply(split(boot, interaction(boot$block, boot$method, drop = TRUE)), function(d) {
      data.frame(
        block = d$block[1],
        method = d$method[1],
        n_scenarios = length(unique(d$scenario_id)),
        coverage_mean = mean(d$coverage, na.rm = TRUE),
        coverage_median = median(d$coverage, na.rm = TRUE),
        mean_width_mean = mean(d$mean_width, na.rm = TRUE),
        mean_width_median = median(d$mean_width, na.rm = TRUE),
        bias_mean = mean(d$bias, na.rm = TRUE),
        bias_median = median(d$bias, na.rm = TRUE),
        mse_mean = mean(d$mse, na.rm = TRUE),
        mse_median = median(d$mse, na.rm = TRUE),
        mean_n_boot_valid_mean = mean(d$mean_n_boot_valid, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))

    boot_block <- boot_block[order(boot_block$block, match(boot_block$method, method_order)), , drop = FALSE]
    write.csv(
      boot_block,
      file.path(boot_dir, "bootstrap_summary_by_block_method_with_gscap.csv"),
      row.names = FALSE
    )
  }

  if (all(c("block", "scenario_id", "method", "mse", "mean_width", "coverage") %in% names(boot))) {
    one_boot_pair <- function(method_new, method_ref) {
      a <- boot[boot$method == method_new, c("block", "scenario_id", "mse", "mean_width", "coverage")]
      b <- boot[boot$method == method_ref, c("block", "scenario_id", "mse", "mean_width", "coverage")]

      if (!nrow(a) || !nrow(b)) return(NULL)

      names(a)[names(a) == "mse"] <- "mse_new"
      names(a)[names(a) == "mean_width"] <- "width_new"
      names(a)[names(a) == "coverage"] <- "coverage_new"

      names(b)[names(b) == "mse"] <- "mse_ref"
      names(b)[names(b) == "mean_width"] <- "width_ref"
      names(b)[names(b) == "coverage"] <- "coverage_ref"

      d <- merge(a, b, by = c("block", "scenario_id"))
      if (!nrow(d)) return(NULL)

      d$mse_reduction <- d$mse_ref - d$mse_new
      d$rel_mse <- d$mse_new / d$mse_ref
      d$width_reduction <- d$width_ref - d$width_new

      do.call(rbind, lapply(split(d, d$block), function(x) {
        data.frame(
          block = x$block[1],
          method_new = method_new,
          method_ref = method_ref,
          n_scenarios = length(unique(x$scenario_id)),
          mse_reduction_mean = mean(x$mse_reduction, na.rm = TRUE),
          mse_reduction_median = median(x$mse_reduction, na.rm = TRUE),
          mse_reduction_prop_positive = mean(x$mse_reduction > 0, na.rm = TRUE),
          rel_mse_mean = mean(x$rel_mse, na.rm = TRUE),
          rel_mse_median = median(x$rel_mse, na.rm = TRUE),
          rel_mse_prop_less_than_1 = mean(x$rel_mse < 1, na.rm = TRUE),
          width_reduction_mean = mean(x$width_reduction, na.rm = TRUE),
          coverage_new_mean = mean(x$coverage_new, na.rm = TRUE),
          coverage_ref_mean = mean(x$coverage_ref, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }))
    }

    boot_pair <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
      one_boot_pair(pairs$method_new[i], pairs$method_ref[i])
    }))

    write.csv(
      boot_pair,
      file.path(boot_dir, "bootstrap_pairwise_summary_gscap.csv"),
      row.names = FALSE
    )
  }
}

cat("Created files:\n")
cat(file.path(tables_dir, "method_win_summary_unique_with_gscap.csv"), "\n")
cat(file.path(tables_dir, "decision_criteria_summary_gscap.csv"), "\n")
if (file.exists(boot_file)) {
  cat(file.path(boot_dir, "bootstrap_summary_by_block_method_with_gscap.csv"), "\n")
  cat(file.path(boot_dir, "bootstrap_pairwise_summary_gscap.csv"), "\n")
}
