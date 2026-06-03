args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript postprocess_focused_gscap_pairwise_augmented.R <RESULTS_DIR>")
}

results <- normalizePath(args[1], mustWork = TRUE)
boot_dir <- file.path(results, "bootstrap")
infile <- file.path(boot_dir, "focused_bootstrap_summary.csv")

if (!file.exists(infile)) {
  stop("Missing file: ", infile)
}

x <- read.csv(infile, stringsAsFactors = FALSE)

needed <- c(
  "scenario_id", "block", "method", "mse", "mean_width", "coverage"
)
miss <- setdiff(needed, names(x))
if (length(miss)) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

pairs <- data.frame(
  method_new = c(
    "gscap_x_l050",
    "gscap_x_l100",
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
    "scap_x_l050",
    "scap_x_l075",
    "scap_x_l075"
  ),
  stringsAsFactors = FALSE
)

one_pair <- function(method_new, method_ref) {
  a <- x[x$method == method_new,
         c("scenario_id", "block", "mse", "mean_width", "coverage"),
         drop = FALSE]
  b <- x[x$method == method_ref,
         c("scenario_id", "block", "mse", "mean_width", "coverage"),
         drop = FALSE]

  if (!nrow(a) || !nrow(b)) return(NULL)

  names(a)[names(a) == "mse"] <- "mse_new"
  names(a)[names(a) == "mean_width"] <- "width_new"
  names(a)[names(a) == "coverage"] <- "coverage_new"

  names(b)[names(b) == "mse"] <- "mse_ref"
  names(b)[names(b) == "mean_width"] <- "width_ref"
  names(b)[names(b) == "coverage"] <- "coverage_ref"

  d <- merge(a, b, by = c("scenario_id", "block"))
  if (!nrow(d)) return(NULL)

  d$method_new <- method_new
  d$method_ref <- method_ref
  d$mse_reduction <- d$mse_ref - d$mse_new
  d$rel_mse <- d$mse_new / d$mse_ref
  d$width_reduction <- d$width_ref - d$width_new

  block_sum <- do.call(rbind, lapply(split(d, d$block), function(z) {
    data.frame(
      block = z$block[1],
      method_new = method_new,
      method_ref = method_ref,
      n_scenarios = length(unique(z$scenario_id)),
      mse_reduction_mean = mean(z$mse_reduction, na.rm = TRUE),
      mse_reduction_median = median(z$mse_reduction, na.rm = TRUE),
      mse_reduction_prop_positive = mean(z$mse_reduction > 0, na.rm = TRUE),
      rel_mse_mean = mean(z$rel_mse, na.rm = TRUE),
      rel_mse_median = median(z$rel_mse, na.rm = TRUE),
      rel_mse_prop_less_than_1 = mean(z$rel_mse < 1, na.rm = TRUE),
      width_reduction_mean = mean(z$width_reduction, na.rm = TRUE),
      coverage_new_mean = mean(z$coverage_new, na.rm = TRUE),
      coverage_ref_mean = mean(z$coverage_ref, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(scenario = d, block = block_sum)
}

ans <- lapply(seq_len(nrow(pairs)), function(i) {
  one_pair(pairs$method_new[i], pairs$method_ref[i])
})

scenario_out <- do.call(rbind, lapply(ans, function(z) if (is.null(z)) NULL else z$scenario))
block_out <- do.call(rbind, lapply(ans, function(z) if (is.null(z)) NULL else z$block))

scenario_file <- file.path(boot_dir, "focused_bootstrap_pairwise_comparisons_gscap_augmented.csv")
block_file <- file.path(boot_dir, "focused_bootstrap_pairwise_summary_gscap_augmented.csv")

write.csv(scenario_out, scenario_file, row.names = FALSE)
write.csv(block_out, block_file, row.names = FALSE)

cat("Created:\n")
cat(scenario_file, "\n")
cat(block_file, "\n\n")

cat("Block-level augmented G-sCAP pairwise summary:\n")
print(block_out)
