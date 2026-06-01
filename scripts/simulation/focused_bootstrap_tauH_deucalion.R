# ==============================================================================
# Focused bootstrap for tau_H simulation results
# ------------------------------------------------------------------------------
# This script reuses the function definitions from
# covariate_adaptive_tauH_simulation_deucalion.R but DOES NOT rerun the Monte Carlo.
# It runs a higher-intensity bootstrap on selected central scenarios from a completed
# paper run, typically paper_blockC_YYYYMMDD_HHMMSS.
# ===============================================================================

# ---- Locate and load definitions from the main simulation script ----------------
PROJECT_DIR <- Sys.getenv("TAUH_PROJECT_DIR", unset = getwd())
MAIN_SCRIPT <- Sys.getenv(
  "TAUH_MAIN_SCRIPT",
  unset = file.path(PROJECT_DIR, "covariate_adaptive_tauH_simulation_deucalion.R")
)
if (!file.exists(MAIN_SCRIPT)) stop("Main script not found: ", MAIN_SCRIPT)

# Evaluate only definitions, stopping before the mode-driven execution block.
main_lines <- readLines(MAIN_SCRIPT, warn = FALSE)
stop_idx <- grep("^mode <- tolower\\(CONFIG\\$RUN_MODE\\)", main_lines)
if (length(stop_idx) != 1L) stop("Could not find execution block marker in main script.")
def_text <- paste(main_lines[seq_len(stop_idx - 1L)], collapse = "\n")
eval(parse(text = def_text), envir = .GlobalEnv)

# ---- Focused bootstrap configuration ------------------------------------------
FOCUSED_MODE <- tolower(Sys.getenv("TAUH_FOCUSED_MODE", unset = "bootstrap"))

CONFIG$BOOT_METHODS <- c(
  "trad",
  "ipcw_marg",
  "ipcw_x",
  "ps_marg_logit",
  "ps_x_logit",
  "scap_x_l050",
  "scap_x_l075"
)
CONFIG$PLOT_METHODS <- CONFIG$BOOT_METHODS
CONFIG$RUN_BOOTSTRAP <- TRUE
CONFIG$N_BOOT_DATASETS <- as.integer(Sys.getenv("TAUH_FOCUSED_N_BOOT_DATASETS", unset = "300"))
CONFIG$N_BOOT_RESAMPLES <- as.integer(Sys.getenv("TAUH_FOCUSED_N_BOOT_RESAMPLES", unset = "500"))
CONFIG$BOOT_N <- as.integer(Sys.getenv("TAUH_FOCUSED_BOOT_N", unset = "300"))
CONFIG$BOOT_TARGET_CENSOR2 <- as.numeric(Sys.getenv("TAUH_FOCUSED_TARGET_CENSOR2", unset = "0.5"))
CONFIG$BOOT_H_QUANTILE <- as.numeric(Sys.getenv("TAUH_FOCUSED_H_QUANTILE", unset = "0.9"))
CONFIG$N_BOOT_TASKS <- as.integer(Sys.getenv("TAUH_N_FOCUSED_BOOT_TASKS", unset = Sys.getenv("TAUH_N_BOOT_TASKS", unset = "8")))
CONFIG$TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = Sys.getenv("TAUH_TASK_ID", unset = "1")))
if (!is.finite(CONFIG$N_BOOT_TASKS) || CONFIG$N_BOOT_TASKS < 1L) CONFIG$N_BOOT_TASKS <- 1L
if (!is.finite(CONFIG$TASK_ID) || CONFIG$TASK_ID < 1L) CONFIG$TASK_ID <- 1L
CONFIG$MAX_CORES <- max(1L, min(64L, as.integer(CONFIG$MAX_CORES)))
CONFIG$RESUME <- TRUE

# Source completed paper run. This run provides scenario_grid_with_truth.csv.
SOURCE_RUN_ID <- Sys.getenv("TAUH_SOURCE_RUN_ID", unset = "paper_blockC_20260531_194657")
SOURCE_RESULTS <- Sys.getenv(
  "TAUH_SOURCE_RESULTS",
  unset = file.path(PROJECT_DIR, "results_covariate_adaptive_tauH", SOURCE_RUN_ID)
)
SOURCE_SCENARIO_CSV <- file.path(SOURCE_RESULTS, "tables", "scenario_grid_with_truth.csv")
if (!file.exists(SOURCE_SCENARIO_CSV)) stop("Scenario grid not found: ", SOURCE_SCENARIO_CSV)

# ---- Focused scenario selection ------------------------------------------------
parse_ids <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(integer(0))
  as.integer(strsplit(x, ",")[[1]])
}

select_focused_scenarios <- function(scenarios) {
  explicit_ids <- parse_ids(Sys.getenv("TAUH_FOCUSED_SCENARIO_IDS", unset = "20,32,86,110,163,164,175,176"))
  if (length(explicit_ids) > 0L) {
    out <- scenarios[scenarios$scenario_id %in% explicit_ids, , drop = FALSE]
    if (nrow(out) == length(explicit_ids)) {
      out <- out[match(explicit_ids, out$scenario_id), , drop = FALSE]
      return(out)
    }
    warning("Not all explicit scenario ids were found; falling back to rule-based selection.")
  }

  s <- scenarios[scenarios$n == CONFIG$BOOT_N &
                   abs(scenarios$target_censor2 - CONFIG$BOOT_TARGET_CENSOR2) < 1e-8 &
                   abs(scenarios$H_quantile - CONFIG$BOOT_H_QUANTILE) < 1e-8, , drop = FALSE]
  if (nrow(s) == 0L) stop("No scenarios match focused bootstrap filters.")

  pick_one <- function(d, block, tau, extra = NULL) {
    z <- d[d$block == block & abs(d$tau_design - tau) < 1e-8, , drop = FALSE]
    if (!is.null(extra)) z <- extra(z)
    if (nrow(z) == 0L) stop("No scenario for block=", block, " tau=", tau)
    z[order(z$scenario_id), , drop = FALSE][1, , drop = FALSE]
  }

  out <- list()
  out[[length(out) + 1L]] <- pick_one(s, "A_independent", 0.0)
  out[[length(out) + 1L]] <- pick_one(s, "A_independent", 0.5)
  out[[length(out) + 1L]] <- pick_one(s, "B_covariate_dependent", 0.0, function(z) z[z$x_type == "binary" & z$cov_effect == "strong", , drop = FALSE])
  out[[length(out) + 1L]] <- pick_one(s, "B_covariate_dependent", 0.5, function(z) z[z$x_type == "binary" & z$cov_effect == "strong", , drop = FALSE])

  for (ceff in c("nonlinear_moderate", "nonlinear_strong")) {
    out[[length(out) + 1L]] <- pick_one(s, "C_nonlinear_observability", 0.0, function(z) z[z$censor_x_effect == ceff, , drop = FALSE])
    out[[length(out) + 1L]] <- pick_one(s, "C_nonlinear_observability", 0.5, function(z) z[z$censor_x_effect == ceff, , drop = FALSE])
  }

  unique(rbind_safe(out))
}

# ---- Run modes -----------------------------------------------------------------
write_focused_config <- function(focused_scenarios) {
  saveRDS(CONFIG, file.path(DIRS$root, "focused_bootstrap_config.rds"))
  write_csv_atomic(focused_scenarios, file.path(DIRS$bootstrap, "focused_bootstrap_selected_scenarios.csv"))
  writeLines(c(
    paste("Focused bootstrap run:", RUN_ID),
    paste("Started:", Sys.time()),
    paste("Mode:", FOCUSED_MODE),
    paste("Source results:", SOURCE_RESULTS),
    paste("Scenario ids:", paste(focused_scenarios$scenario_id, collapse = ",")),
    paste("Methods:", paste(CONFIG$BOOT_METHODS, collapse = ",")),
    paste("N_BOOT_DATASETS:", CONFIG$N_BOOT_DATASETS),
    paste("N_BOOT_RESAMPLES:", CONFIG$N_BOOT_RESAMPLES),
    paste("N_BOOT_TASKS:", CONFIG$N_BOOT_TASKS),
    paste("MAX_CORES:", CONFIG$MAX_CORES)
  ), file.path(DIRS$logs, "focused_bootstrap_run_log.txt"))
}

run_focused_bootstrap <- function() {
  scenarios <- read.csv(SOURCE_SCENARIO_CSV, stringsAsFactors = FALSE)
  focused <- select_focused_scenarios(scenarios)
  write_focused_config(focused)

  idx <- task_indices(nrow(focused), CONFIG$N_BOOT_TASKS, CONFIG$TASK_ID)
  message("Focused bootstrap task ", CONFIG$TASK_ID, "/", CONFIG$N_BOOT_TASKS,
          " processing focused scenario indices: ", paste(idx, collapse = ","),
          "; scenario ids: ", paste(focused$scenario_id[idx], collapse = ","))
  if (length(idx) == 0L) return(invisible(data.frame()))
  out <- lapply(idx, function(i) run_bootstrap_for_scenario(focused[i, , drop = FALSE]))
  raw <- rbind_safe(out)
  task_file <- file.path(DIRS$checkpoints, sprintf("focused_bootstrap_task_%04d_of_%04d_done.txt", CONFIG$TASK_ID, CONFIG$N_BOOT_TASKS))
  writeLines(c(paste("completed", Sys.time()), paste("focused_scenarios", paste(focused$scenario_id[idx], collapse = ","))), task_file)
  invisible(raw)
}

summarize_by_block_method <- function(boot_summary) {
  if (nrow(boot_summary) == 0L) return(data.frame())
  key <- paste(boot_summary$block, boot_summary$method, sep = "||")
  rbind_safe(lapply(split(boot_summary, key), function(d) {
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
}

run_focused_collect <- function() {
  boot_summary <- collect_bootstrap_files()
  write_csv_atomic(boot_summary, file.path(DIRS$bootstrap, "focused_bootstrap_summary.csv"))
  bm <- summarize_by_block_method(boot_summary)
  if (nrow(bm) > 0L) {
    bm <- bm[order(bm$block, bm$mse_mean, bm$method), , drop = FALSE]
    write_csv_atomic(bm, file.path(DIRS$bootstrap, "focused_bootstrap_summary_by_block_method.csv"))
  }

  # Compact method comparison versus traditional and IPCW marginal within each selected scenario.
  if (nrow(boot_summary) > 0L) {
    cmp_rows <- list()
    for (sid in unique(boot_summary$scenario_id)) {
      d <- boot_summary[boot_summary$scenario_id == sid, , drop = FALSE]
      for (ref in c("trad", "ipcw_marg")) {
        dref <- d[d$method == ref, , drop = FALSE]
        if (nrow(dref) != 1L) next
        for (m in setdiff(unique(d$method), ref)) {
          dm <- d[d$method == m, , drop = FALSE]
          if (nrow(dm) != 1L) next
          cmp_rows[[length(cmp_rows) + 1L]] <- data.frame(
            scenario_id = sid,
            block = d$block[1],
            method_new = m,
            method_ref = ref,
            mse_new = dm$mse,
            mse_ref = dref$mse,
            mse_reduction = dref$mse - dm$mse,
            rel_mse = dm$mse / dref$mse,
            width_new = dm$mean_width,
            width_ref = dref$mean_width,
            width_reduction = dref$mean_width - dm$mean_width,
            coverage_new = dm$coverage,
            coverage_ref = dref$coverage,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    cmp <- rbind_safe(cmp_rows)
    if (nrow(cmp) > 0L) {
      write_csv_atomic(cmp, file.path(DIRS$bootstrap, "focused_bootstrap_pairwise_comparisons.csv"))
      key <- paste(cmp$block, cmp$method_new, cmp$method_ref, sep = "||")
      cmp_sum <- rbind_safe(lapply(split(cmp, key), function(d) {
        data.frame(
          block = d$block[1],
          method_new = d$method_new[1],
          method_ref = d$method_ref[1],
          n_scenarios = nrow(d),
          mse_reduction_mean = mean(d$mse_reduction, na.rm = TRUE),
          mse_reduction_median = median(d$mse_reduction, na.rm = TRUE),
          mse_reduction_prop_positive = mean(d$mse_reduction > 0, na.rm = TRUE),
          rel_mse_mean = mean(d$rel_mse, na.rm = TRUE),
          rel_mse_median = median(d$rel_mse, na.rm = TRUE),
          rel_mse_prop_less_than_1 = mean(d$rel_mse < 1, na.rm = TRUE),
          width_reduction_mean = mean(d$width_reduction, na.rm = TRUE),
          coverage_new_mean = mean(d$coverage_new, na.rm = TRUE),
          coverage_ref_mean = mean(d$coverage_ref, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }))
      write_csv_atomic(cmp_sum, file.path(DIRS$bootstrap, "focused_bootstrap_pairwise_summary.csv"))
    }
  }

  manifest <- data.frame(path = list.files(OUT_DIR, recursive = TRUE, full.names = TRUE), stringsAsFactors = FALSE)
  if (nrow(manifest) > 0L) {
    manifest$relative_path <- sub(paste0("^", normalizePath(OUT_DIR), "/?"), "", normalizePath(manifest$path, winslash = "/"))
    manifest$size_bytes <- file.info(manifest$path)$size
    write_csv_atomic(manifest, file.path(DIRS$root, "manifest.csv"))
  }
  writeLines(c(
    paste("Completed:", Sys.time()),
    paste("Run ID:", RUN_ID),
    paste("Source results:", SOURCE_RESULTS),
    paste("Bootstrap summary rows:", nrow(boot_summary))
  ), file.path(DIRS$logs, "completion_log_focused_collect.txt"))
  message("Focused collect completed. Results saved in: ", OUT_DIR)
  invisible(boot_summary)
}

message("Starting focused bootstrap. mode=", FOCUSED_MODE,
        ", run_id=", RUN_ID,
        ", out=", OUT_DIR,
        ", source=", SOURCE_RESULTS,
        ", max_cores=", CONFIG$MAX_CORES,
        ", task_id=", CONFIG$TASK_ID)

if (FOCUSED_MODE %in% c("bootstrap", "boot", "run")) {
  run_focused_bootstrap()
} else if (FOCUSED_MODE %in% c("collect", "summary")) {
  run_focused_collect()
} else if (FOCUSED_MODE == "all") {
  old_tasks <- CONFIG$N_BOOT_TASKS
  old_task <- CONFIG$TASK_ID
  CONFIG$N_BOOT_TASKS <- 1L
  CONFIG$TASK_ID <- 1L
  run_focused_bootstrap()
  CONFIG$N_BOOT_TASKS <- old_tasks
  CONFIG$TASK_ID <- old_task
  run_focused_collect()
} else {
  stop("Unknown TAUH_FOCUSED_MODE: ", FOCUSED_MODE)
}
