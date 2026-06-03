# ==============================================================================
# Covariate-adaptive estimation of normalized truncated Kendall tau for censored
# sequential gap-time data: simulation blocks A, B and C
# ------------------------------------------------------------------------------
# Block A: independent right censoring. This is the negative-control block.
#          Classical marginal methods should be strong and new adaptive methods
#          should not be advertised as universally superior.
#
# Block B: covariate-dependent right censoring. We generate X, (A,B)|X and C|X
#          so that C is dependent on (A,B) marginally but independent conditional
#          on X. This block is designed to assess whether covariate-adaptive IPCW
#          and stabilized presmoothing reduce bias/MSE relative to methods that
#          ignore X.
#
# Block C: nonlinear observability. C is still generated conditionally on observed
#          covariates, but the censoring distribution is nonlinear in X and in
#          time. This block is designed to assess whether flexible nuisance
#          learning (splines/ensembles) matters beyond linear adjustment.
#
# The target is always the normalized truncated Kendall tau
#   tau_H = E[ psi((A1,B1),(A2,B2)) | A1+B1 <= H, A2+B2 <= H ],
# where psi is concordance minus discordance with ties contributing zero.
#
# This script saves all numerical outputs and figures in a timestamped folder.
# It is designed for workstations with up to 32 cores.
# ===============================================================================

# ==============================================================================
# 0. USER CONFIGURATION
# ===============================================================================

CONFIG <- list(
  PROFILE = "paper",                      # "test" or "paper"
  OUTPUT_ROOT = file.path(getwd(), "results_covariate_adaptive_tauH"),
  MASTER_SEED = 20260505L,
  MAX_CORES = 32L,
  USE_PARALLEL = TRUE,
  RESUME = TRUE,
  INSTALL_MISSING_PACKAGES = FALSE,        # set TRUE if you want optional packages auto-installed

  # Main methods. Methods ending in _x use X; methods ending in _marg ignore X.
  METHODS = c(
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
    "scap_x_uncal_l050",
    "scap_x_nocf_l050"
  ),

  # Smaller method set for bootstrap coverage because bootstrap is expensive.
  BOOT_METHODS = c("trad", "ipcw_marg", "ipcw_x", "gscap_x_l050", "gscap_x_l100", "ps_x_spline", "scap_x_l050", "scap_x_l075"),

  # Methods highlighted in figures.
  PLOT_METHODS = c("trad", "ipcw_marg", "ipcw_x", "gscap_x_l050", "gscap_x_l100", "ps_marg_logit", "ps_x_spline", "scap_x_l050", "scap_x_l075"),

  EPS_PROB = 0.01,
  K_FOLDS = 5L,
  X_STRATA_Q = 4L,
  SPLINE_DF = 3L,

  # Optional ML learners. If ranger/gbm are absent, the ensemble still runs using
  # constant/logit/spline learners.
  ACAP_LEARNERS_X = c("constant", "logit", "logit_int", "spline", "rf", "gbm"),
  ACAP_LEARNERS_MARG = c("constant", "logit", "spline", "rf", "gbm"),
  RF_NUM_TREES = 300L,
  RF_MIN_NODE_SIZE = 10L,
  GBM_N_TREES = 250L,
  GBM_INTERACTION_DEPTH = 2L,
  GBM_SHRINKAGE = 0.03,
  GBM_BAG_FRACTION = 0.80,

  # Truth and Monte Carlo sizes.
  N_TRUTH = 200000L,
  N_REP_MAIN = 1000L,
  N_REP_SENS = 500L,

  # Bootstrap coverage study.
  RUN_BOOTSTRAP = TRUE,
  N_BOOT_DATASETS = 120L,
  N_BOOT_RESAMPLES = 250L,
  BOOT_N = 300L,
  BOOT_TARGET_CENSOR2 = 0.50,
  BOOT_H_QUANTILE = 0.90,

  # Scenario grid.
  N_GRID = c(100L, 300L, 500L),
  TAU_GRID = c(-0.50, 0.00, 0.50, 0.75),
  TARGET_CENSOR2_GRID = c(0.30, 0.50),
  H_QUANTILE_GRID_A = c(0.80, 0.90),
  H_QUANTILE_GRID_B = c(0.90),
  H_QUANTILE_GRID_C = c(0.90),
  X_TYPES_A = c("binary"),
  X_TYPES_B = c("binary", "normal"),
  X_TYPES_C = c("normal"),
  COV_EFFECTS_B = c("moderate", "strong"),
  CENSOR_EFFECTS_C = c("nonlinear_moderate", "nonlinear_strong"),
  MARGIN = "exp",

  # Effects used in conditional data generation.
  # A and B both depend on X; in block B, C also depends on X.
  EFFECT_AB_MODERATE = c(beta_A = 0.45, beta_B = -0.35),
  EFFECT_AB_STRONG   = c(beta_A = 0.70, beta_B = -0.55),
  EFFECT_C_MODERATE  = 0.90,
  EFFECT_C_STRONG    = 1.35,
  EFFECT_C_NONLINEAR_MODERATE = 0.95,
  EFFECT_C_NONLINEAR_STRONG   = 1.45
)


# ---- Environment overrides for Deucalion/SLURM --------------------------------
# These overrides are intentionally read before the PROFILE == "test" block so
# that TAUH_PROFILE=test can be used from the command line or SLURM job script.
env_profile <- Sys.getenv("TAUH_PROFILE", unset = "")
if (nzchar(env_profile)) CONFIG$PROFILE <- env_profile

if (CONFIG$PROFILE == "test") {
  CONFIG$OUTPUT_ROOT <- file.path(getwd(), "results_covariate_adaptive_tauH_TEST")
  CONFIG$MAX_CORES <- min(4L, CONFIG$MAX_CORES)
  CONFIG$N_TRUTH <- 30000L
  CONFIG$N_REP_MAIN <- 30L
  CONFIG$N_REP_SENS <- 20L
  CONFIG$N_GRID <- c(120L)
  CONFIG$TAU_GRID <- c(0.00, 0.50)
  CONFIG$TARGET_CENSOR2_GRID <- c(0.50)
  CONFIG$H_QUANTILE_GRID_A <- c(0.90)
  CONFIG$H_QUANTILE_GRID_B <- c(0.90)
  CONFIG$X_TYPES_B <- c("binary", "normal")
  CONFIG$X_TYPES_C <- c("normal")
  CONFIG$COV_EFFECTS_B <- c("strong")
  CONFIG$CENSOR_EFFECTS_C <- c("nonlinear_strong")
  CONFIG$RUN_BOOTSTRAP <- TRUE
  CONFIG$N_BOOT_DATASETS <- 8L
  CONFIG$N_BOOT_RESAMPLES <- 30L
  CONFIG$BOOT_N <- 120L
  CONFIG$RF_NUM_TREES <- 80L
  CONFIG$GBM_N_TREES <- 80L
}


# ---- Additional environment overrides, applied after test/profile defaults -----
env_output <- Sys.getenv("TAUH_OUTPUT_ROOT", unset = "")
if (nzchar(env_output)) CONFIG$OUTPUT_ROOT <- env_output

env_cores <- Sys.getenv("TAUH_CORES", unset = "")
if (!nzchar(env_cores)) env_cores <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
if (nzchar(env_cores)) CONFIG$MAX_CORES <- max(1L, min(64L, as.integer(env_cores)))

env_parallel <- Sys.getenv("TAUH_USE_PARALLEL", unset = "")
if (nzchar(env_parallel)) CONFIG$USE_PARALLEL <- as.logical(as.integer(env_parallel))

env_boot <- Sys.getenv("TAUH_RUN_BOOTSTRAP", unset = "")
if (nzchar(env_boot)) CONFIG$RUN_BOOTSTRAP <- as.logical(as.integer(env_boot))

CONFIG$RUN_MODE <- Sys.getenv("TAUH_MODE", unset = "all")
CONFIG$N_SCENARIO_TASKS <- as.integer(Sys.getenv("TAUH_N_SCENARIO_TASKS", unset = Sys.getenv("TAUH_N_TASKS", unset = "1")))
CONFIG$N_BOOT_TASKS <- as.integer(Sys.getenv("TAUH_N_BOOT_TASKS", unset = Sys.getenv("TAUH_N_TASKS", unset = "1")))
CONFIG$TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = Sys.getenv("TAUH_TASK_ID", unset = "1")))
if (!is.finite(CONFIG$N_SCENARIO_TASKS) || CONFIG$N_SCENARIO_TASKS < 1L) CONFIG$N_SCENARIO_TASKS <- 1L
if (!is.finite(CONFIG$N_BOOT_TASKS) || CONFIG$N_BOOT_TASKS < 1L) CONFIG$N_BOOT_TASKS <- 1L
if (!is.finite(CONFIG$TASK_ID) || CONFIG$TASK_ID < 1L) CONFIG$TASK_ID <- 1L

# Avoid accidental oversubscription from BLAS/OpenMP when many R workers run.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ==============================================================================
# 1. PACKAGE MANAGEMENT AND OUTPUT DIRECTORIES
# ===============================================================================

maybe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (isTRUE(CONFIG$INSTALL_MISSING_PACKAGES)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  requireNamespace(pkg, quietly = TRUE)
}

HAS_RANGER <- maybe_install("ranger")
HAS_GBM <- maybe_install("gbm")
HAS_GGPLOT2 <- maybe_install("ggplot2")
HAS_SPLINES <- requireNamespace("splines", quietly = TRUE)

RUN_ID <- Sys.getenv("TAUH_RUN_ID", unset = "")
if (!nzchar(RUN_ID)) RUN_ID <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", CONFIG$PROFILE)
OUT_DIR <- file.path(CONFIG$OUTPUT_ROOT, RUN_ID)
DIRS <- list(
  root = OUT_DIR,
  raw = file.path(OUT_DIR, "raw_by_scenario"),
  tables = file.path(OUT_DIR, "tables"),
  figures = file.path(OUT_DIR, "figures"),
  bootstrap = file.path(OUT_DIR, "bootstrap"),
  logs = file.path(OUT_DIR, "logs"),
  checkpoints = file.path(OUT_DIR, "checkpoints")
)
for (d in DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  "Covariate-adaptive tau_H simulation run",
  paste("Run ID:", RUN_ID),
  paste("Started:", as.character(Sys.time())),
  paste("Profile:", CONFIG$PROFILE),
  paste("R version:", R.version.string),
  paste("HAS_RANGER:", HAS_RANGER),
  paste("HAS_GBM:", HAS_GBM),
  paste("HAS_GGPLOT2:", HAS_GGPLOT2)
), file.path(DIRS$logs, "run_log.txt"))

saveRDS(CONFIG, file.path(DIRS$root, "config.rds"))
writeLines(capture.output(sessionInfo()), file.path(DIRS$root, "sessionInfo.txt"))

# ==============================================================================
# 2. GENERIC UTILITIES
# ===============================================================================

safe_seed_value <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) x <- CONFIG$MASTER_SEED
  as.integer((abs(floor(x)) %% 2147483646) + 1L)
}

make_seed <- function(...) {
  vals <- unlist(list(...))
  vals <- as.numeric(vals)
  vals[!is.finite(vals)] <- 0
  multipliers <- c(1000003, 9176, 104729, 1009, 7919, 65537, 31337, 27183)
  multipliers <- rep(multipliers, length.out = length(vals))
  s <- sum(vals * multipliers) + CONFIG$MASTER_SEED
  safe_seed_value(s)
}

clip_prob <- function(p, eps = CONFIG$EPS_PROB) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  pmin(1 - eps, pmax(eps, p))
}

safe_mean_binary <- function(y, eps = CONFIG$EPS_PROB) {
  y <- y[is.finite(y)]
  if (length(y) == 0L) return(0.5)
  clip_prob(mean(y), eps = eps)
}

softmax <- function(theta) {
  theta <- as.numeric(theta)
  theta <- theta - max(theta, na.rm = TRUE)
  ex <- exp(theta)
  ex / sum(ex)
}

logit_safe <- function(p, eps = CONFIG$EPS_PROB) {
  qlogis(clip_prob(p, eps = eps))
}

make_folds <- function(ids, k, seed) {
  ids <- as.integer(ids)
  if (length(ids) == 0L) return(vector("list", k))
  set.seed(safe_seed_value(seed))
  ids <- sample(ids)
  k <- max(2L, min(as.integer(k), length(ids)))
  split(ids, rep(seq_len(k), length.out = length(ids)))
}

rbind_safe <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (length(lst) == 0L) return(data.frame())
  do.call(rbind, lst)
}

method_uses_m <- function(method) {
  grepl("^ps_", method) || grepl("^scap_", method) || grepl("^gscap_", method)
}

method_uses_x <- function(method) {
  grepl("_x", method)
}

lambda_from_method <- function(method) {
  if (grepl("l025", method)) return(0.25)
  if (grepl("l050", method)) return(0.50)
  if (grepl("l075", method)) return(0.75)
  if (grepl("l100", method)) return(1.00)
  if (grepl("^ps_", method)) return(1.00)
  NA_real_
}

# ==============================================================================
# 3. DATA GENERATION: COPULA, MARGINS, X, CENSORING
# ===============================================================================

rho_from_tau_gaussian <- function(tau) {
  sin(pi * tau / 2)
}

r_positive_stable <- function(n, alpha) {
  if (alpha >= 1) return(rep(1, n))
  U <- stats::runif(n, 0, pi)
  E <- stats::rexp(n, rate = 1)
  (sin(alpha * U) / (sin(U))^(1 / alpha)) *
    (sin((1 - alpha) * U) / E)^((1 - alpha) / alpha)
}

r_copula_uv <- function(n, family = "gaussian", param = 0) {
  family <- tolower(family)
  if (family == "gaussian") {
    rho <- max(-0.999999, min(0.999999, as.numeric(param)))
    z1 <- stats::rnorm(n)
    z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(n)
    return(cbind(stats::pnorm(z1), stats::pnorm(z2)))
  }
  if (family == "clayton") {
    theta <- as.numeric(param)
    if (theta <= 1e-12) return(cbind(stats::runif(n), stats::runif(n)))
    W <- stats::rgamma(n, shape = 1 / theta, rate = 1)
    E1 <- stats::rexp(n)
    E2 <- stats::rexp(n)
    return(cbind((1 + E1 / W)^(-1 / theta), (1 + E2 / W)^(-1 / theta)))
  }
  if (family == "gumbel") {
    theta <- as.numeric(param)
    if (theta <= 1 + 1e-12) return(cbind(stats::runif(n), stats::runif(n)))
    alpha <- 1 / theta
    S <- r_positive_stable(n, alpha)
    E1 <- stats::rexp(n)
    E2 <- stats::rexp(n)
    return(cbind(exp(-(E1 / S)^alpha), exp(-(E2 / S)^alpha)))
  }
  stop("Unknown copula family: ", family)
}

copula_param_from_tau <- function(family, tau) {
  family <- tolower(family)
  if (family == "gaussian") return(rho_from_tau_gaussian(tau))
  if (family == "clayton") {
    if (tau <= 0) return(0)
    return(2 * tau / (1 - tau))
  }
  if (family == "gumbel") {
    if (tau <= 0) return(1)
    return(1 / (1 - tau))
  }
  stop("Unknown copula family: ", family)
}

copula_tau_full <- function(family, param) {
  family <- tolower(family)
  if (family == "gaussian") return(2 / pi * asin(as.numeric(param)))
  if (family == "clayton") {
    theta <- as.numeric(param)
    if (theta <= 0) return(0)
    return(theta / (theta + 2))
  }
  if (family == "gumbel") {
    theta <- as.numeric(param)
    if (theta <= 1) return(0)
    return(1 - 1 / theta)
  }
  NA_real_
}

q_margin_base <- function(u, margin = "exp", which = "A") {
  u <- pmin(1 - 1e-12, pmax(1e-12, u))
  margin <- tolower(margin)
  if (margin == "exp") {
    rate <- if (which == "A") 1.00 else 0.85
    return(stats::qexp(u, rate = rate))
  }
  if (margin == "weibull") {
    if (which == "A") return(stats::qweibull(u, shape = 1.25, scale = 1.00))
    return(stats::qweibull(u, shape = 0.85, scale = 1.25))
  }
  if (margin == "lognormal") {
    if (which == "A") return(stats::qlnorm(u, meanlog = -0.15, sdlog = 0.70))
    return(stats::qlnorm(u, meanlog = 0.05, sdlog = 0.85))
  }
  stop("Unknown margin: ", margin)
}

simulate_X <- function(n, x_type) {
  x_type <- tolower(x_type)
  if (x_type == "binary") {
    X <- stats::rbinom(n, 1L, 0.5)
    # Fixed centering avoids changing the data-generating law from sample to sample.
    Xc <- X - 0.5
    return(data.frame(X = X, Xc = Xc, x_type = x_type))
  }
  if (x_type == "normal") {
    X <- stats::rnorm(n)
    # X is already standard normal; no sample standardization is used in the DGP.
    Xc <- X
    return(data.frame(X = X, Xc = Xc, x_type = x_type))
  }
  stop("Unknown x_type: ", x_type)
}

simulate_true_gap_times_x <- function(n, scenario) {
  Xdf <- simulate_X(n, scenario$x_type)
  uv <- r_copula_uv(n, family = scenario$copula_family, param = scenario$copula_param)
  A0 <- q_margin_base(uv[, 1], margin = scenario$margin, which = "A")
  B0 <- q_margin_base(uv[, 2], margin = scenario$margin, which = "B")

  if (scenario$cov_effect == "none") {
    beta_A <- 0
    beta_B <- 0
  } else if (scenario$cov_effect == "moderate") {
    beta_A <- CONFIG$EFFECT_AB_MODERATE["beta_A"]
    beta_B <- CONFIG$EFFECT_AB_MODERATE["beta_B"]
  } else if (scenario$cov_effect == "strong") {
    beta_A <- CONFIG$EFFECT_AB_STRONG["beta_A"]
    beta_B <- CONFIG$EFFECT_AB_STRONG["beta_B"]
  } else if (scenario$cov_effect == "nonlinear") {
    beta_A <- NA_real_
    beta_B <- NA_real_
  } else {
    stop("Unknown cov_effect: ", scenario$cov_effect)
  }

  # Multiplicative accelerated gap-time effects. X affects the gap distribution.
  if (scenario$cov_effect == "nonlinear") {
    if (tolower(scenario$x_type) == "normal") {
      eta_A <- 0.45 * sin(1.25 * Xdf$Xc) + 0.18 * (Xdf$Xc^2 - 1)
      eta_B <- -0.40 * cos(1.10 * Xdf$Xc) + 0.20 * Xdf$Xc
    } else {
      eta_A <- 0.65 * Xdf$Xc
      eta_B <- -0.50 * Xdf$Xc
    }
    A <- A0 * exp(eta_A)
    B <- B0 * exp(eta_B)
  } else {
    A <- A0 * exp(beta_A * Xdf$Xc)
    B <- B0 * exp(beta_B * Xdf$Xc)
  }
  data.frame(X = Xdf$X, Xc = Xdf$Xc, A = A, B = B, Z = A + B)
}

censor_model_name <- function(scenario) {
  if ("censor_model" %in% names(scenario) && is.finite(match("censor_model", names(scenario)))) {
    cm <- as.character(scenario$censor_model)
    if (length(cm) > 0L && nzchar(cm[1]) && !is.na(cm[1])) return(cm[1])
  }
  "exp"
}

censor_rate_multiplier <- function(Xc, scenario) {
  if (scenario$block == "A_independent") return(rep(1, length(Xc)))
  eff <- as.character(scenario$censor_x_effect)
  if (eff == "moderate") {
    gamma <- CONFIG$EFFECT_C_MODERATE
    return(exp(gamma * Xc))
  }
  if (eff == "strong") {
    gamma <- CONFIG$EFFECT_C_STRONG
    return(exp(gamma * Xc))
  }
  if (eff == "nonlinear_moderate") {
    gamma <- CONFIG$EFFECT_C_NONLINEAR_MODERATE
    return(exp(gamma * (0.70 * sin(1.40 * Xc) + 0.30 * (Xc^2 - 1))))
  }
  if (eff == "nonlinear_strong") {
    gamma <- CONFIG$EFFECT_C_NONLINEAR_STRONG
    return(exp(gamma * (0.75 * sin(1.35 * Xc) + 0.35 * (Xc^2 - 1) + 0.15 * Xc)))
  }
  rep(1, length(Xc))
}

censor_shape_vector <- function(Xc, scenario) {
  cm <- censor_model_name(scenario)
  if (cm == "weibull_nonlinear") {
    # Positive, nonlinear time-shape. This makes m0(a,z,x) nonlinear in both z and X.
    return(0.70 + 0.65 * stats::plogis(1.25 * Xc) + 0.15 * sin(1.10 * Xc))
  }
  rep(1, length(Xc))
}

calibrate_censor_base_rate <- function(Z, Xc, scenario, target_censor2) {
  mult <- censor_rate_multiplier(Xc, scenario)
  shape <- censor_shape_vector(Xc, scenario)
  f <- function(base) mean(1 - exp(-base * mult * (Z^shape))) - target_censor2
  lo <- 1e-8
  hi <- 1
  while (f(hi) < 0 && hi < 1e6) hi <- hi * 2
  if (hi >= 1e6) return(hi)
  stats::uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
}

simulate_observed_data_x <- function(n, scenario) {
  true <- simulate_true_gap_times_x(n, scenario)
  mult <- censor_rate_multiplier(true$Xc, scenario)
  shape <- censor_shape_vector(true$Xc, scenario)
  rate <- scenario$censor_base_rate * mult
  Uc <- stats::runif(n)
  C <- ((-log(pmax(Uc, 1e-12))) / pmax(rate, 1e-12))^(1 / shape)
  Atilde <- pmin(true$A, C)
  Ztilde <- pmin(true$Z, C)
  d1 <- as.integer(true$A <= C)
  d2 <- as.integer(true$Z <= C)
  Aobs <- ifelse(d1 == 1L, Atilde, NA_real_)
  Bobs <- ifelse(d1 == 1L, Ztilde - Atilde, NA_real_)  # correct: retain partial censored second gap
  data.frame(
    X = true$X, Xc = true$Xc,
    A = true$A, B = true$B, Z = true$Z, C = C,
    Atilde = Atilde, Ztilde = Ztilde,
    d1 = d1, d2 = d2,
    Aobs = Aobs, Bobs = Bobs,
    censor_rate_i = rate
  )
}

# ==============================================================================
# 4. FAST WEIGHTED KENDALL AND TARGET tau_H
# ===============================================================================

weighted_kendall_num_fast <- function(A, B, w) {
  keep <- is.finite(A) & is.finite(B) & is.finite(w) & (w > 0)
  A <- as.numeric(A[keep])
  B <- as.numeric(B[keep])
  w <- as.numeric(w[keep])
  m <- length(w)
  if (m < 2L) return(0)

  # Fenwick tree over B ranks, processing by ordered A. Ties in A contribute zero.
  B_levels <- sort(unique(B))
  Br <- match(B, B_levels)
  nb <- length(B_levels)
  ord <- order(A, B)
  A <- A[ord]
  Br <- Br[ord]
  w <- w[ord]
  bit <- numeric(nb)

  bit_update <- function(idx, val) {
    while (idx <= nb) {
      bit[idx] <<- bit[idx] + val
      idx <- idx + bitwAnd(idx, -idx)
    }
    invisible(NULL)
  }
  bit_query <- function(idx) {
    if (idx <= 0L) return(0)
    s <- 0
    while (idx > 0L) {
      s <- s + bit[idx]
      idx <- idx - bitwAnd(idx, -idx)
    }
    s
  }

  total_prev <- 0
  num_unordered <- 0
  start <- 1L
  while (start <= m) {
    end <- start
    while (end < m && A[end + 1L] == A[start]) end <- end + 1L
    for (q in start:end) {
      r <- Br[q]
      less <- bit_query(r - 1L)
      leq <- bit_query(r)
      greater <- total_prev - leq
      num_unordered <- num_unordered + w[q] * (less - greater)
    }
    for (q in start:end) bit_update(Br[q], w[q])
    total_prev <- total_prev + sum(w[start:end])
    start <- end + 1L
  }
  2 * num_unordered
}

estimate_truth_tau_H <- function(A, B, H) {
  keep <- is.finite(A) & is.finite(B) & ((A + B) <= H)
  m <- sum(keep)
  if (m < 2L) return(NA_real_)
  w <- rep(1 / m, m)
  weighted_kendall_num_fast(A[keep], B[keep], w) / (sum(w)^2)
}

# ==============================================================================
# 5. PRODUCT-LIMIT WEIGHTS AND IPCW
# ===============================================================================

product_limit_weights <- function(e, Ztilde) {
  n <- length(e)
  e <- as.numeric(e)
  e[!is.finite(e)] <- 0
  e <- pmin(1, pmax(0, e))
  Ztilde <- as.numeric(Ztilde)
  ord <- order(Ztilde, seq_along(Ztilde), na.last = TRUE)
  e_ord <- e[ord]
  risk <- rev(seq_len(n))
  omega_ord <- numeric(n)
  surv_before <- 1
  for (k in seq_len(n)) {
    haz <- e_ord[k] / risk[k]
    haz <- min(1, max(0, haz))
    omega_ord[k] <- surv_before * haz
    surv_before <- surv_before * (1 - haz)
  }
  omega <- numeric(n)
  omega[ord] <- omega_ord
  omega
}

estimate_tau_pl <- function(df, H, e) {
  omega <- product_limit_weights(e = e, Ztilde = df$Ztilde)
  IH <- as.numeric(df$d1 == 1L & is.finite(df$Aobs) & is.finite(df$Bobs) & df$Ztilde <= H)
  wH <- omega * IH
  MH <- sum(wH, na.rm = TRUE)
  if (!is.finite(MH) || MH <= 0) {
    return(list(tau_hat = NA_real_, numerator = NA_real_, M_hat = MH,
                sum_omega = sum(omega, na.rm = TRUE)))
  }
  num <- weighted_kendall_num_fast(df$Aobs, df$Bobs, wH)
  tau <- max(-1, min(1, num / (MH^2)))
  list(tau_hat = tau, numerator = num, M_hat = MH, sum_omega = sum(omega, na.rm = TRUE))
}

km_fit_left <- function(time, event, eps = CONFIG$EPS_PROB) {
  # event=1 is censoring event observed for G. Returns step function G(t-).
  time <- as.numeric(time)
  event <- as.integer(event)
  keep <- is.finite(time) & is.finite(event)
  time <- time[keep]
  event <- event[keep]
  n <- length(time)
  if (n == 0L) return(list(times = numeric(0), surv_before = numeric(0), surv_after = numeric(0), eps = eps))
  ord <- order(time, seq_along(time))
  time <- time[ord]
  event <- event[ord]
  times_u <- sort(unique(time))
  surv_before <- numeric(length(times_u))
  surv_after <- numeric(length(times_u))
  surv <- 1
  risk <- n
  pos <- 1L
  for (j in seq_along(times_u)) {
    t <- times_u[j]
    start <- pos
    while (pos <= n && time[pos] == t) pos <- pos + 1L
    idx <- start:(pos - 1L)
    surv_before[j] <- surv
    d <- sum(event[idx] == 1L, na.rm = TRUE)
    if (risk > 0) surv <- surv * (1 - d / risk)
    surv <- max(surv, eps)
    surv_after[j] <- surv
    risk <- risk - length(idx)
  }
  list(times = times_u, surv_before = pmax(surv_before, eps), surv_after = pmax(surv_after, eps), eps = eps)
}

km_predict_left <- function(fit, t) {
  t <- as.numeric(t)
  if (length(fit$times) == 0L) return(rep(1, length(t)))
  # For continuous simulated times, exact ties are rare. Use survival after times < t.
  idx <- findInterval(t, fit$times)
  out <- rep(1, length(t))
  ii <- which(idx > 0L)
  out[ii] <- fit$surv_after[idx[ii]]
  pmax(out, fit$eps)
}

make_x_strata_info <- function(X, x_type, q = CONFIG$X_STRATA_Q) {
  x_type <- tolower(x_type)
  if (x_type == "binary") {
    return(list(type = "binary", breaks = NULL))
  }
  q <- max(2L, as.integer(q))
  br <- unique(stats::quantile(X, probs = seq(0, 1, length.out = q + 1L), na.rm = TRUE, type = 7))
  if (length(br) < 3L) return(list(type = "single", breaks = NULL))
  br[1] <- -Inf
  br[length(br)] <- Inf
  list(type = "quantile", breaks = br)
}

assign_x_strata <- function(X, info) {
  if (info$type == "binary") return(as.character(as.integer(X > 0.5)))
  if (info$type == "single") return(rep("all", length(X)))
  as.character(cut(X, breaks = info$breaks, include.lowest = TRUE, labels = FALSE))
}

predict_G_marginal <- function(df, eps = CONFIG$EPS_PROB) {
  fit <- km_fit_left(df$Ztilde, 1L - df$d2, eps = eps)
  km_predict_left(fit, df$Ztilde)
}

predict_G_x_crossfit <- function(df, eps = CONFIG$EPS_PROB, seed = 1L) {
  n <- nrow(df)
  folds <- make_folds(seq_len(n), CONFIG$K_FOLDS, seed = seed)
  G <- rep(NA_real_, n)
  for (f in seq_along(folds)) {
    pred_idx <- folds[[f]]
    train_idx <- setdiff(seq_len(n), pred_idx)
    info <- make_x_strata_info(df$X[train_idx], df$x_type[1], q = CONFIG$X_STRATA_Q)
    s_train <- assign_x_strata(df$X[train_idx], info)
    s_pred <- assign_x_strata(df$X[pred_idx], info)
    global_fit <- km_fit_left(df$Ztilde[train_idx], 1L - df$d2[train_idx], eps = eps)
    for (s in unique(s_pred)) {
      idp <- pred_idx[s_pred == s]
      idt <- train_idx[s_train == s]
      if (length(idt) < 15L || sum((1L - df$d2[idt]) == 1L) < 2L) {
        G[idp] <- km_predict_left(global_fit, df$Ztilde[idp])
      } else {
        fit_s <- km_fit_left(df$Ztilde[idt], 1L - df$d2[idt], eps = eps)
        G[idp] <- km_predict_left(fit_s, df$Ztilde[idp])
      }
    }
  }
  clip_prob(G, eps = eps)
}

estimate_tau_ipcw <- function(df, H, use_x = FALSE, seed = 1L) {
  G <- if (isTRUE(use_x)) predict_G_x_crossfit(df, seed = seed) else predict_G_marginal(df)
  IH <- as.numeric(df$d1 == 1L & df$d2 == 1L & is.finite(df$Aobs) & is.finite(df$Bobs) & df$Ztilde <= H)
  w <- IH / clip_prob(G)
  wH <- w / nrow(df)
  MH <- sum(wH, na.rm = TRUE)
  if (!is.finite(MH) || MH <= 0) {
    return(list(tau_hat = NA_real_, numerator = NA_real_, M_hat = MH,
                sum_omega = sum(wH, na.rm = TRUE), mean_G = mean(G, na.rm = TRUE)))
  }
  num <- weighted_kendall_num_fast(df$Aobs, df$Bobs, wH)
  tau <- max(-1, min(1, num / (MH^2)))
  list(tau_hat = tau, numerator = num, M_hat = MH, sum_omega = sum(wH, na.rm = TRUE),
       mean_G = mean(G, na.rm = TRUE))
}

# ==============================================================================
# 6. PRESMOOTHING NUISANCE MODELS m(A,Ztilde,X)
# ===============================================================================

standardize_train_pred <- function(train, pred) {
  mA <- mean(train$Aobs, na.rm = TRUE); sA <- stats::sd(train$Aobs, na.rm = TRUE)
  mZ <- mean(train$Ztilde, na.rm = TRUE); sZ <- stats::sd(train$Ztilde, na.rm = TRUE)
  mX <- mean(train$X, na.rm = TRUE); sX <- stats::sd(train$X, na.rm = TRUE)
  if (!is.finite(sA) || sA <= 0) sA <- 1
  if (!is.finite(sZ) || sZ <= 0) sZ <- 1
  if (!is.finite(sX) || sX <= 0) sX <- 1
  train$As <- (train$Aobs - mA) / sA
  train$Zs <- (train$Ztilde - mZ) / sZ
  train$Xs <- (train$X - mX) / sX
  pred$As <- (pred$Aobs - mA) / sA
  pred$Zs <- (pred$Ztilde - mZ) / sZ
  pred$Xs <- (pred$X - mX) / sX
  list(train = train, pred = pred)
}

fit_predict_constant <- function(df, train_idx, pred_idx, eps = CONFIG$EPS_PROB) {
  p0 <- safe_mean_binary(df$d2[train_idx], eps = eps)
  rep(p0, length(pred_idx))
}

fit_predict_logit <- function(df, train_idx, pred_idx, include_x = FALSE, interaction = FALSE, eps = CONFIG$EPS_PROB) {
  train <- df[train_idx, , drop = FALSE]
  pred <- df[pred_idx, , drop = FALSE]
  p0 <- safe_mean_binary(train$d2, eps = eps)
  if (nrow(train) < 10L || length(unique(train$d2)) < 2L) return(rep(p0, length(pred_idx)))
  sp <- standardize_train_pred(train, pred)
  train <- sp$train; pred <- sp$pred
  if (include_x && interaction) {
    form <- d2 ~ As + Zs + Xs + As:Xs + Zs:Xs + As:Zs
  } else if (include_x) {
    form <- d2 ~ As + Zs + Xs
  } else {
    form <- d2 ~ As + Zs
  }
  fit <- tryCatch(
    suppressWarnings(stats::glm(form, family = stats::binomial(), data = train)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(p0, length(pred_idx)))
  pr <- tryCatch(stats::predict(fit, newdata = pred, type = "response"), error = function(e) rep(p0, length(pred_idx)))
  clip_prob(pr, eps = eps)
}

fit_predict_spline <- function(df, train_idx, pred_idx, include_x = FALSE, eps = CONFIG$EPS_PROB) {
  train <- df[train_idx, , drop = FALSE]
  pred <- df[pred_idx, , drop = FALSE]
  p0 <- safe_mean_binary(train$d2, eps = eps)
  if (!HAS_SPLINES || nrow(train) < 25L || length(unique(train$d2)) < 2L) {
    return(fit_predict_logit(df, train_idx, pred_idx, include_x = include_x, eps = eps))
  }
  sp <- standardize_train_pred(train, pred)
  train <- sp$train; pred <- sp$pred
  df_use <- max(2L, CONFIG$SPLINE_DF)
  if (include_x) {
    form <- stats::as.formula(paste0(
      "d2 ~ splines::ns(As, df=", df_use, ") + splines::ns(Zs, df=", df_use, ") + ",
      "splines::ns(Xs, df=", df_use, ") + As:Zs + As:Xs + Zs:Xs"
    ))
  } else {
    form <- stats::as.formula(paste0(
      "d2 ~ splines::ns(As, df=", df_use, ") + splines::ns(Zs, df=", df_use, ") + As:Zs"
    ))
  }
  fit <- tryCatch(
    suppressWarnings(stats::glm(form, family = stats::binomial(), data = train)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(fit_predict_logit(df, train_idx, pred_idx, include_x = include_x, eps = eps))
  pr <- tryCatch(stats::predict(fit, newdata = pred, type = "response"), error = function(e) rep(p0, length(pred_idx)))
  clip_prob(pr, eps = eps)
}

fit_predict_rf <- function(df, train_idx, pred_idx, include_x = FALSE, eps = CONFIG$EPS_PROB, seed = 1L) {
  p0 <- safe_mean_binary(df$d2[train_idx], eps = eps)
  if (!HAS_RANGER) return(rep(p0, length(pred_idx)))
  train <- df[train_idx, , drop = FALSE]
  pred <- df[pred_idx, , drop = FALSE]
  if (nrow(train) < 20L || length(unique(train$d2)) < 2L) return(rep(p0, length(pred_idx)))
  sp <- standardize_train_pred(train, pred)
  train <- sp$train; pred <- sp$pred
  train$y <- factor(train$d2, levels = c(0, 1))
  form <- if (include_x) y ~ As + Zs + Xs else y ~ As + Zs
  fit <- tryCatch(
    ranger::ranger(
      form, data = train, probability = TRUE, num.trees = CONFIG$RF_NUM_TREES,
      min.node.size = CONFIG$RF_MIN_NODE_SIZE,
      mtry = if (include_x) 2L else 1L,
      seed = safe_seed_value(seed), num.threads = 1L
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(p0, length(pred_idx)))
  pr <- tryCatch({
    pp <- predict(fit, data = pred)$predictions
    if (is.matrix(pp) && "1" %in% colnames(pp)) pp[, "1"] else rep(p0, length(pred_idx))
  }, error = function(e) rep(p0, length(pred_idx)))
  clip_prob(pr, eps = eps)
}

fit_predict_gbm <- function(df, train_idx, pred_idx, include_x = FALSE, eps = CONFIG$EPS_PROB, seed = 1L) {
  p0 <- safe_mean_binary(df$d2[train_idx], eps = eps)
  if (!HAS_GBM) return(rep(p0, length(pred_idx)))
  train <- df[train_idx, , drop = FALSE]
  pred <- df[pred_idx, , drop = FALSE]
  if (nrow(train) < 25L || length(unique(train$d2)) < 2L) return(rep(p0, length(pred_idx)))
  sp <- standardize_train_pred(train, pred)
  train <- sp$train; pred <- sp$pred
  form <- if (include_x) d2 ~ As + Zs + Xs else d2 ~ As + Zs
  set.seed(safe_seed_value(seed))
  fit <- tryCatch(
    gbm::gbm(
      form, data = train, distribution = "bernoulli",
      n.trees = CONFIG$GBM_N_TREES,
      interaction.depth = CONFIG$GBM_INTERACTION_DEPTH,
      shrinkage = CONFIG$GBM_SHRINKAGE,
      bag.fraction = CONFIG$GBM_BAG_FRACTION,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(p0, length(pred_idx)))
  pr <- tryCatch(
    stats::predict(fit, newdata = pred, n.trees = CONFIG$GBM_N_TREES, type = "response"),
    error = function(e) rep(p0, length(pred_idx))
  )
  clip_prob(pr, eps = eps)
}

available_learners <- function(learners) {
  out <- learners
  if (!HAS_RANGER) out <- setdiff(out, "rf")
  if (!HAS_GBM) out <- setdiff(out, "gbm")
  if (!"constant" %in% out) out <- c("constant", out)
  unique(out)
}

fit_predict_learner <- function(learner, df, train_idx, pred_idx, include_x = FALSE, eps = CONFIG$EPS_PROB, seed = 1L) {
  learner <- tolower(learner)
  if (learner == "constant") return(fit_predict_constant(df, train_idx, pred_idx, eps))
  if (learner == "logit") return(fit_predict_logit(df, train_idx, pred_idx, include_x = include_x, interaction = FALSE, eps = eps))
  if (learner == "logit_int") return(fit_predict_logit(df, train_idx, pred_idx, include_x = include_x, interaction = TRUE, eps = eps))
  if (learner == "spline") return(fit_predict_spline(df, train_idx, pred_idx, include_x = include_x, eps = eps))
  if (learner == "rf") return(fit_predict_rf(df, train_idx, pred_idx, include_x = include_x, eps = eps, seed = seed))
  if (learner == "gbm") return(fit_predict_gbm(df, train_idx, pred_idx, include_x = include_x, eps = eps, seed = seed))
  stop("Unknown learner: ", learner)
}

fit_m_single_cf <- function(df, learner, include_x = FALSE, eps = CONFIG$EPS_PROB, seed = 1L) {
  ids <- which(df$d1 == 1L & is.finite(df$Aobs) & is.finite(df$Ztilde))
  n <- nrow(df)
  base_p <- safe_mean_binary(df$d2[ids], eps = eps)
  mhat <- rep(base_p, n)
  if (length(ids) < 8L) return(clip_prob(mhat, eps))
  folds <- make_folds(ids, CONFIG$K_FOLDS, seed)
  for (f in seq_along(folds)) {
    pred_idx <- folds[[f]]
    train_idx <- setdiff(ids, pred_idx)
    if (length(train_idx) < 5L) train_idx <- ids
    mhat[pred_idx] <- fit_predict_learner(learner, df, train_idx, pred_idx, include_x = include_x,
                                          eps = eps, seed = make_seed(seed, f, 101))
  }
  clip_prob(mhat, eps)
}

calibrate_mhat <- function(y, p, eps = CONFIG$EPS_PROB) {
  p <- clip_prob(p, eps)
  y <- as.numeric(y)
  if (length(y) < 20L || length(unique(y[is.finite(y)])) < 2L) return(p)
  q <- logit_safe(p, eps)
  if (stats::sd(q, na.rm = TRUE) <= 1e-10) return(p)
  fit <- tryCatch(
    suppressWarnings(stats::glm(y ~ q, family = stats::binomial(), data = data.frame(y = y, q = q))),
    error = function(e) NULL
  )
  if (is.null(fit)) return(p)
  pr <- tryCatch(stats::predict(fit, newdata = data.frame(q = q), type = "response"), error = function(e) p)
  clip_prob(pr, eps)
}

fit_m_ensemble <- function(df, include_x = TRUE, crossfit = TRUE, calibrate = TRUE,
                           eps = CONFIG$EPS_PROB, seed = 1L) {
  ids <- which(df$d1 == 1L & is.finite(df$Aobs) & is.finite(df$Ztilde))
  n <- nrow(df)
  base_p <- safe_mean_binary(df$d2[ids], eps = eps)
  learner_set <- if (include_x) CONFIG$ACAP_LEARNERS_X else CONFIG$ACAP_LEARNERS_MARG
  learners <- available_learners(learner_set)
  if (!include_x) learners <- setdiff(learners, "logit_int")
  if (!"constant" %in% learners) learners <- c("constant", learners)
  learners <- unique(learners)

  empty <- function() {
    alpha <- setNames(rep(0, length(learners)), learners)
    alpha["constant"] <- 1
    list(mhat = clip_prob(rep(base_p, n), eps), alpha = alpha, learners = learners,
         calibrated = FALSE, stack_loss = NA_real_)
  }
  if (length(ids) < 8L) return(empty())

  P <- matrix(NA_real_, nrow = n, ncol = length(learners))
  colnames(P) <- learners
  if (isTRUE(crossfit)) {
    folds <- make_folds(ids, CONFIG$K_FOLDS, seed)
    for (f in seq_along(folds)) {
      pred_idx <- folds[[f]]
      train_idx <- setdiff(ids, pred_idx)
      if (length(train_idx) < 5L) train_idx <- ids
      for (l in seq_along(learners)) {
        P[pred_idx, l] <- fit_predict_learner(learners[l], df, train_idx, pred_idx,
                                              include_x = include_x, eps = eps,
                                              seed = make_seed(seed, f, l, 77))
      }
    }
  } else {
    for (l in seq_along(learners)) {
      P[ids, l] <- fit_predict_learner(learners[l], df, ids, ids,
                                       include_x = include_x, eps = eps,
                                       seed = make_seed(seed, l, 88))
    }
  }
  for (l in seq_along(learners)) {
    P[ids, l] <- clip_prob(P[ids, l], eps)
    P[!is.finite(P[, l]), l] <- base_p
  }

  y <- as.numeric(df$d2[ids])
  P_ids <- P[ids, , drop = FALSE]
  valid <- which(colSums(!is.finite(P_ids)) == 0L)
  if (length(valid) == 0L) return(empty())
  P_ids <- P_ids[, valid, drop = FALSE]
  active_learners <- colnames(P_ids)

  neg_logloss <- function(theta) {
    a <- softmax(theta)
    p <- clip_prob(as.numeric(P_ids %*% a), eps)
    mean(-y * log(p) - (1 - y) * log(1 - p))
  }
  opt <- tryCatch(stats::optim(rep(0, ncol(P_ids)), neg_logloss, method = "BFGS", control = list(maxit = 200)),
                  error = function(e) NULL)
  if (is.null(opt)) {
    losses <- colMeans(-y * log(P_ids) - (1 - y) * log(1 - P_ids))
    best <- which.min(losses)
    a <- rep(0, ncol(P_ids)); a[best] <- 1
    stack_loss <- losses[best]
  } else {
    a <- softmax(opt$par)
    stack_loss <- opt$value
  }

  alpha <- setNames(rep(0, length(learners)), learners)
  alpha[active_learners] <- a
  p_ids <- clip_prob(as.numeric(P_ids %*% a), eps)
  if (isTRUE(calibrate)) p_ids <- calibrate_mhat(y, p_ids, eps)
  mhat <- rep(base_p, n)
  mhat[ids] <- p_ids
  list(mhat = clip_prob(mhat, eps), alpha = alpha, learners = learners,
       calibrated = isTRUE(calibrate), stack_loss = stack_loss)
}

calibration_metrics <- function(y, p, n_bins = 10L, eps = CONFIG$EPS_PROB) {
  keep <- is.finite(y) & is.finite(p)
  y <- as.numeric(y[keep])
  p <- clip_prob(p[keep], eps)
  if (length(y) == 0L) {
    return(list(logloss = NA_real_, brier = NA_real_, cal_intercept = NA_real_,
                cal_slope = NA_real_, ece = NA_real_, mean_mhat = NA_real_, mean_y = NA_real_, n_cal = 0L))
  }
  logloss <- mean(-y * log(p) - (1 - y) * log(1 - p))
  brier <- mean((y - p)^2)
  q <- logit_safe(p, eps)
  cal_intercept <- NA_real_; cal_slope <- NA_real_
  if (length(unique(y)) >= 2L && length(y) >= 20L && stats::sd(q) > 0) {
    fit <- tryCatch(suppressWarnings(stats::glm(y ~ q, family = stats::binomial(), data = data.frame(y = y, q = q))),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      cf <- stats::coef(fit)
      cal_intercept <- as.numeric(cf[1])
      cal_slope <- as.numeric(cf[2])
    }
  }
  ord <- order(p)
  bin_id <- cut(seq_along(ord), breaks = min(n_bins, length(ord)), labels = FALSE)
  ece <- 0
  for (b in unique(bin_id)) {
    idx <- ord[bin_id == b]
    ece <- ece + length(idx) / length(y) * abs(mean(y[idx]) - mean(p[idx]))
  }
  list(logloss = logloss, brier = brier, cal_intercept = cal_intercept, cal_slope = cal_slope,
       ece = ece, mean_mhat = mean(p), mean_y = mean(y), n_cal = length(y))
}

# ==============================================================================
# 7. METHOD ESTIMATION WRAPPER
# ===============================================================================

estimate_method <- function(df, H, method, seed = 1L) {
  method <- tolower(method)
  out_extra <- list(
    m_logloss = NA_real_, m_brier = NA_real_, m_cal_intercept = NA_real_,
    m_cal_slope = NA_real_, m_ece = NA_real_, mean_mhat = NA_real_,
    mean_d2_given_d1 = NA_real_, n_cal = NA_real_, mean_G = NA_real_,
    alpha_constant = NA_real_, alpha_logit = NA_real_, alpha_logit_int = NA_real_,
    alpha_spline = NA_real_, alpha_rf = NA_real_, alpha_gbm = NA_real_
  )

  if (method == "trad") {
    est <- estimate_tau_pl(df, H, e = df$d2)
    return(c(est, out_extra))
  }
  if (method == "ipcw_marg") {
    est <- estimate_tau_ipcw(df, H, use_x = FALSE, seed = seed)
    out_extra$mean_G <- est$mean_G
    return(c(est, out_extra))
  }
  if (method == "ipcw_x") {
    est <- estimate_tau_ipcw(df, H, use_x = TRUE, seed = seed)
    out_extra$mean_G <- est$mean_G
    return(c(est, out_extra))
  }

  # Presmoothing methods.
  ids <- which(df$d1 == 1L & is.finite(df$Aobs) & is.finite(df$Ztilde))
  include_x <- method_uses_x(method)
  lambda <- lambda_from_method(method)
  if (!is.finite(lambda)) lambda <- 1.0

  if (method == "ps_marg_logit") {
    mhat <- fit_m_single_cf(df, learner = "logit", include_x = FALSE, seed = seed)
    alpha <- NULL
  } else if (method == "ps_x_logit") {
    mhat <- fit_m_single_cf(df, learner = "logit", include_x = TRUE, seed = seed)
    alpha <- NULL
  } else if (method == "ps_x_spline") {
    mhat <- fit_m_single_cf(df, learner = "spline", include_x = TRUE, seed = seed)
    alpha <- NULL
  } else if (method == "scap_marg_l050") {
    fit <- fit_m_ensemble(df, include_x = FALSE, crossfit = TRUE, calibrate = TRUE, seed = seed)
    mhat <- fit$mhat; alpha <- fit$alpha
  } else if (method == "scap_x_uncal_l050") {
    fit <- fit_m_ensemble(df, include_x = TRUE, crossfit = TRUE, calibrate = FALSE, seed = seed)
    mhat <- fit$mhat; alpha <- fit$alpha
  } else if (method == "scap_x_nocf_l050") {
    fit <- fit_m_ensemble(df, include_x = TRUE, crossfit = FALSE, calibrate = TRUE, seed = seed)
    mhat <- fit$mhat; alpha <- fit$alpha
  } else if (grepl("^scap_x_l", method)) {
    fit <- fit_m_ensemble(df, include_x = TRUE, crossfit = TRUE, calibrate = TRUE, seed = seed)
    mhat <- fit$mhat; alpha <- fit$alpha
  } else if (grepl("^gscap_x_l", method)) {
    fit <- fit_m_ensemble(df, include_x = TRUE, crossfit = TRUE, calibrate = TRUE, seed = seed)
    mhat <- fit$mhat; alpha <- fit$alpha
  } else {
    stop("Unknown method: ", method)
  }

  e <- df$d1 * ((1 - lambda) * df$d2 + lambda * mhat)

  if (grepl("^gscap_x_l", method)) {
    G <- predict_G_x_crossfit(df, seed = seed)
    G <- clip_prob(G)

    IH <- as.numeric(df$d1 == 1L &
                     is.finite(df$Aobs) &
                     is.finite(df$Bobs) &
                     is.finite(df$Ztilde) &
                     df$Ztilde <= H)

    wH <- e * IH / G / nrow(df)
    wH[!is.finite(wH)] <- 0

    MH <- sum(wH, na.rm = TRUE)

    if (!is.finite(MH) || MH <= 0) {
      est <- list(
        tau_hat = NA_real_,
        numerator = NA_real_,
        M_hat = MH,
        sum_omega = sum(wH, na.rm = TRUE),
        mean_G = mean(G, na.rm = TRUE)
      )
    } else {
      num <- weighted_kendall_num_fast(df$Aobs, df$Bobs, wH)
      tau <- max(-1, min(1, num / (MH^2)))
      est <- list(
        tau_hat = tau,
        numerator = num,
        M_hat = MH,
        sum_omega = sum(wH, na.rm = TRUE),
        mean_G = mean(G, na.rm = TRUE)
      )
    }

    out_extra$mean_G <- est$mean_G
  } else {
    est <- estimate_tau_pl(df, H, e = e)
  }

  if (length(ids) > 0L) {
    cm <- calibration_metrics(df$d2[ids], mhat[ids])
    out_extra$m_logloss <- cm$logloss
    out_extra$m_brier <- cm$brier
    out_extra$m_cal_intercept <- cm$cal_intercept
    out_extra$m_cal_slope <- cm$cal_slope
    out_extra$m_ece <- cm$ece
    out_extra$mean_mhat <- cm$mean_mhat
    out_extra$mean_d2_given_d1 <- cm$mean_y
    out_extra$n_cal <- cm$n_cal
  }
  if (!is.null(alpha)) {
    for (nm in names(alpha)) {
      key <- paste0("alpha_", nm)
      if (key %in% names(out_extra)) out_extra[[key]] <- alpha[[nm]]
    }
  }
  c(est, out_extra)
}

as_result_row <- function(est_list) {
  data.frame(
    tau_hat = as.numeric(est_list$tau_hat),
    numerator = as.numeric(est_list$numerator),
    M_hat = as.numeric(est_list$M_hat),
    sum_omega = as.numeric(est_list$sum_omega),
    mean_G = as.numeric(est_list$mean_G),
    m_logloss = as.numeric(est_list$m_logloss),
    m_brier = as.numeric(est_list$m_brier),
    m_cal_intercept = as.numeric(est_list$m_cal_intercept),
    m_cal_slope = as.numeric(est_list$m_cal_slope),
    m_ece = as.numeric(est_list$m_ece),
    mean_mhat = as.numeric(est_list$mean_mhat),
    mean_d2_given_d1 = as.numeric(est_list$mean_d2_given_d1),
    n_cal = as.numeric(est_list$n_cal),
    alpha_constant = as.numeric(est_list$alpha_constant),
    alpha_logit = as.numeric(est_list$alpha_logit),
    alpha_logit_int = as.numeric(est_list$alpha_logit_int),
    alpha_spline = as.numeric(est_list$alpha_spline),
    alpha_rf = as.numeric(est_list$alpha_rf),
    alpha_gbm = as.numeric(est_list$alpha_gbm)
  )
}

# ==============================================================================
# 8. SCENARIO GRID AND TRUTH COMPUTATION
# ===============================================================================

build_scenario_grid <- function() {
  rows <- list()
  sid <- 0L
  for (tau in CONFIG$TAU_GRID) {
    for (n in CONFIG$N_GRID) {
      for (cens in CONFIG$TARGET_CENSOR2_GRID) {
        for (Hq in CONFIG$H_QUANTILE_GRID_A) {
          for (x_type in CONFIG$X_TYPES_A) {
            sid <- sid + 1L
            rows[[length(rows) + 1L]] <- data.frame(
              scenario_id = sid,
              block = "A_independent",
              scenario_group = "negative_control_independent_censoring",
              n = n,
              n_rep = CONFIG$N_REP_MAIN,
              copula_family = "gaussian",
              tau_design = tau,
              copula_param = copula_param_from_tau("gaussian", tau),
              tau_full = tau,
              margin = CONFIG$MARGIN,
              x_type = x_type,
              cov_effect = "moderate",          # X may affect A,B, but C ignores X
              censor_x_effect = "none",
              censor_model = "exp",
              target_censor2 = cens,
              H_quantile = Hq,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  for (tau in CONFIG$TAU_GRID) {
    for (n in CONFIG$N_GRID) {
      for (cens in CONFIG$TARGET_CENSOR2_GRID) {
        for (Hq in CONFIG$H_QUANTILE_GRID_B) {
          for (x_type in CONFIG$X_TYPES_B) {
            for (eff in CONFIG$COV_EFFECTS_B) {
              sid <- sid + 1L
              rows[[length(rows) + 1L]] <- data.frame(
                scenario_id = sid,
                block = "B_covariate_dependent",
                scenario_group = "covariate_dependent_censoring",
                n = n,
                n_rep = CONFIG$N_REP_MAIN,
                copula_family = "gaussian",
                tau_design = tau,
                copula_param = copula_param_from_tau("gaussian", tau),
                tau_full = tau,
                margin = CONFIG$MARGIN,
                x_type = x_type,
                cov_effect = eff,
                censor_x_effect = eff,
                censor_model = "exp",
                target_censor2 = cens,
                H_quantile = Hq,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  for (tau in CONFIG$TAU_GRID) {
    for (n in CONFIG$N_GRID) {
      for (cens in CONFIG$TARGET_CENSOR2_GRID) {
        for (Hq in CONFIG$H_QUANTILE_GRID_C) {
          for (x_type in CONFIG$X_TYPES_C) {
            for (ceff in CONFIG$CENSOR_EFFECTS_C) {
              sid <- sid + 1L
              rows[[length(rows) + 1L]] <- data.frame(
                scenario_id = sid,
                block = "C_nonlinear_observability",
                scenario_group = "nonlinear_covariate_dependent_censoring",
                n = n,
                n_rep = CONFIG$N_REP_MAIN,
                copula_family = "gaussian",
                tau_design = tau,
                copula_param = copula_param_from_tau("gaussian", tau),
                tau_full = tau,
                margin = CONFIG$MARGIN,
                x_type = x_type,
                cov_effect = "nonlinear",
                censor_x_effect = ceff,
                censor_model = "weibull_nonlinear",
                target_censor2 = cens,
                H_quantile = Hq,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  scenarios <- rbind_safe(rows)
  scenarios$total_scenarios <- nrow(scenarios)
  scenarios
}

compute_scenario_truth <- function(scenario) {
  set.seed(make_seed(scenario$scenario_id, 1001))
  truth <- simulate_true_gap_times_x(CONFIG$N_TRUTH, scenario)
  H <- as.numeric(stats::quantile(truth$Z, probs = scenario$H_quantile, type = 7, na.rm = TRUE))
  tau_H <- estimate_truth_tau_H(truth$A, truth$B, H)
  tau_untruncated <- estimate_truth_tau_H(truth$A, truth$B, Inf)
  base_rate <- calibrate_censor_base_rate(truth$Z, truth$Xc, scenario, scenario$target_censor2)
  mult <- censor_rate_multiplier(truth$Xc, scenario)
  p_d2 <- mean(exp(-base_rate * mult * truth$Z))
  p_d1 <- mean(exp(-base_rate * mult * truth$A))
  scenario$H <- H
  scenario$tau_H_true <- tau_H
  scenario$tau_untruncated_true <- tau_untruncated
  scenario$censor_base_rate <- base_rate
  scenario$p_d1_theory <- p_d1
  scenario$p_d2_theory <- p_d2
  scenario$target_censor2_check <- 1 - p_d2
  scenario$N_TRUTH <- CONFIG$N_TRUTH
  scenario
}


# ============================================================================== 
# 8b. DEUCALION/HPC EXECUTION HELPERS
# ============================================================================== 

hpc_ncores <- function() {
  dc <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
  max(1L, min(as.integer(CONFIG$MAX_CORES), as.integer(dc)))
}

task_indices <- function(n, n_tasks, task_id) {
  n <- as.integer(n)
  n_tasks <- max(1L, as.integer(n_tasks))
  task_id <- max(1L, as.integer(task_id))
  if (n <= 0L) return(integer(0))
  task_id <- ((task_id - 1L) %% n_tasks) + 1L
  which(((seq_len(n) - 1L) %% n_tasks) == (task_id - 1L))
}

write_csv_atomic <- function(df, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  utils::write.csv(df, tmp, row.names = FALSE)
  file.rename(tmp, path)
}

scenario_grid_path <- function() file.path(DIRS$root, "scenario_grid_with_truth.rds")

ensure_scenarios_loaded <- function() {
  rds <- scenario_grid_path()
  if (!file.exists(rds)) {
    stop("Scenario/truth file not found: ", rds, ". Run TAUH_MODE=truth first.")
  }
  readRDS(rds)
}

compute_truth_all <- function() {
  base_file <- file.path(DIRS$tables, "scenario_grid_base.csv")
  truth_file_csv <- file.path(DIRS$tables, "scenario_grid_with_truth.csv")
  truth_file_rds <- scenario_grid_path()
  if (isTRUE(CONFIG$RESUME) && file.exists(truth_file_rds)) {
    message("Truth grid already exists; loading: ", truth_file_rds)
    return(readRDS(truth_file_rds))
  }
  scenarios_base <- build_scenario_grid()
  write_csv_atomic(scenarios_base, base_file)
  message("Computing truth targets for ", nrow(scenarios_base), " scenarios using up to ", hpc_ncores(), " workers...")
  idx <- seq_len(nrow(scenarios_base))
  if (isTRUE(CONFIG$USE_PARALLEL) && hpc_ncores() > 1L && length(idx) > 1L) {
    cl <- parallel::makeCluster(hpc_ncores())
    try(parallel::clusterExport(cl, varlist = setdiff(ls(globalenv()), c("cl")), envir = globalenv()), silent = TRUE)
    scen_list <- parallel::parLapply(cl, idx, function(i) compute_scenario_truth(scenarios_base[i, , drop = FALSE]))
    try(parallel::stopCluster(cl), silent = TRUE)
  } else {
    scen_list <- lapply(idx, function(i) compute_scenario_truth(scenarios_base[i, , drop = FALSE]))
  }
  scenarios <- rbind_safe(scen_list)
  write_csv_atomic(scenarios, truth_file_csv)
  saveRDS(scenarios, truth_file_rds)
  scenarios
}

simulate_one_rep <- function(scenario, rep_id, methods = CONFIG$METHODS) {
  seed <- make_seed(scenario$scenario_id, rep_id, 123)
  set.seed(seed)
  df <- simulate_observed_data_x(scenario$n, scenario)
  df$x_type <- scenario$x_type
  rows <- list()
  for (method in methods) {
    est <- tryCatch(
      estimate_method(df, H = scenario$H, method = method, seed = make_seed(seed, match(method, CONFIG$METHODS), 999)),
      error = function(e) list(error = conditionMessage(e))
    )
    if (!is.null(est$error)) {
      row <- data.frame(
        scenario_id = scenario$scenario_id,
        rep_id = rep_id,
        method = method,
        status = paste0("error: ", est$error),
        tau_hat = NA_real_, numerator = NA_real_, M_hat = NA_real_, sum_omega = NA_real_,
        mean_G = NA_real_, m_logloss = NA_real_, m_brier = NA_real_,
        m_cal_intercept = NA_real_, m_cal_slope = NA_real_, m_ece = NA_real_,
        mean_mhat = NA_real_, mean_d2_given_d1 = NA_real_, n_cal = NA_real_,
        alpha_constant = NA_real_, alpha_logit = NA_real_, alpha_logit_int = NA_real_,
        alpha_spline = NA_real_, alpha_rf = NA_real_, alpha_gbm = NA_real_
      )
    } else {
      row <- cbind(
        data.frame(scenario_id = scenario$scenario_id, rep_id = rep_id, method = method, status = "ok"),
        as_result_row(est)
      )
    }
    rows[[length(rows) + 1L]] <- row
  }
  diag <- data.frame(
    scenario_id = scenario$scenario_id,
    rep_id = rep_id,
    censor1_rate = mean(df$d1 == 0L),
    censor2_rate = mean(df$d2 == 0L),
    n_d1 = sum(df$d1 == 1L),
    n_d2 = sum(df$d2 == 1L),
    n_in_H_observed = sum(df$d1 == 1L & df$Ztilde <= scenario$H),
    mean_X = mean(df$X),
    sd_X = stats::sd(df$X)
  )
  res <- rbind_safe(rows)
  merge(res, diag, by = c("scenario_id", "rep_id"), all.x = TRUE)
}

simulate_scenario <- function(scenario_row) {
  scenario <- scenario_row
  sid <- scenario$scenario_id
  out_file <- file.path(DIRS$raw, sprintf("scenario_%04d_raw.csv", sid))
  if (isTRUE(CONFIG$RESUME) && file.exists(out_file)) {
    return(read.csv(out_file))
  }
  message("Scenario ", sid, "/", scenario$total_scenarios, ": ", scenario$block,
          ", n=", scenario$n, ", tau=", scenario$tau_design,
          ", censor=", scenario$target_censor2, ", x=", scenario$x_type,
          ", effect=", scenario$cov_effect)
  rep_ids <- seq_len(scenario$n_rep)
  ncores <- hpc_ncores()
  if (!isTRUE(CONFIG$USE_PARALLEL) || ncores <= 1L || length(rep_ids) <= 1L) {
    reps <- lapply(rep_ids, function(r) simulate_one_rep(scenario, r, CONFIG$METHODS))
  } else {
    message("  Parallelizing replications with ", ncores, " workers on scenario ", sid)
    cl <- parallel::makeCluster(ncores)
    try(parallel::clusterExport(cl, varlist = setdiff(ls(globalenv()), c("cl")), envir = globalenv()), silent = TRUE)
    reps <- parallel::parLapply(cl, rep_ids, function(r) simulate_one_rep(scenario, r, CONFIG$METHODS))
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  raw <- rbind_safe(reps)
  raw <- merge(raw, scenario, by = "scenario_id", all.x = TRUE)
  write.csv(raw, out_file, row.names = FALSE)
  raw
}


# ==============================================================================
# 10. SUMMARIES, DECISION CRITERIA AND NUMERICAL OUTPUTS
# ============================================================================== 

summarize_mc <- function(raw) {
  ok <- raw[raw$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) return(data.frame())
  split_key <- interaction(ok$scenario_id, ok$method, drop = TRUE)
  pieces <- split(ok, split_key)
  out <- lapply(pieces, function(d) {
    tau_true <- d$tau_H_true[1]
    bias_vec <- d$tau_hat - tau_true
    data.frame(
      scenario_id = d$scenario_id[1],
      block = d$block[1],
      scenario_group = d$scenario_group[1],
      n = d$n[1],
      x_type = d$x_type[1],
      cov_effect = d$cov_effect[1],
      censor_x_effect = d$censor_x_effect[1],
      target_censor2 = d$target_censor2[1],
      H_quantile = d$H_quantile[1],
      H = d$H[1],
      tau_design = d$tau_design[1],
      tau_full = d$tau_full[1],
      tau_H_true = tau_true,
      method = d$method[1],
      n_rep_total = length(d$tau_hat),
      n_rep_valid = sum(is.finite(d$tau_hat)),
      fail_rate = mean(d$status != "ok", na.rm = TRUE),
      mean_tau = mean(d$tau_hat, na.rm = TRUE),
      sd_tau = stats::sd(d$tau_hat, na.rm = TRUE),
      mcse_mean_tau = stats::sd(d$tau_hat, na.rm = TRUE) / sqrt(sum(is.finite(d$tau_hat))),
      median_tau = stats::median(d$tau_hat, na.rm = TRUE),
      q025_tau = as.numeric(stats::quantile(d$tau_hat, 0.025, na.rm = TRUE)),
      q975_tau = as.numeric(stats::quantile(d$tau_hat, 0.975, na.rm = TRUE)),
      bias = mean(bias_vec, na.rm = TRUE),
      abs_bias = abs(mean(bias_vec, na.rm = TRUE)),
      mse = mean(bias_vec^2, na.rm = TRUE),
      rmse = sqrt(mean(bias_vec^2, na.rm = TRUE)),
      mean_M_hat = mean(d$M_hat, na.rm = TRUE),
      sd_M_hat = stats::sd(d$M_hat, na.rm = TRUE),
      mean_sum_omega = mean(d$sum_omega, na.rm = TRUE),
      mean_censor1_rate = mean(d$censor1_rate, na.rm = TRUE),
      mean_censor2_rate = mean(d$censor2_rate, na.rm = TRUE),
      mean_n_d1 = mean(d$n_d1, na.rm = TRUE),
      mean_n_d2 = mean(d$n_d2, na.rm = TRUE),
      mean_n_in_H_observed = mean(d$n_in_H_observed, na.rm = TRUE),
      mean_logloss = mean(d$m_logloss, na.rm = TRUE),
      mean_brier = mean(d$m_brier, na.rm = TRUE),
      mean_cal_intercept = mean(d$m_cal_intercept, na.rm = TRUE),
      mean_cal_slope = mean(d$m_cal_slope, na.rm = TRUE),
      mean_ece = mean(d$m_ece, na.rm = TRUE),
      mean_mhat = mean(d$mean_mhat, na.rm = TRUE),
      mean_d2_given_d1 = mean(d$mean_d2_given_d1, na.rm = TRUE),
      mean_alpha_constant = mean(d$alpha_constant, na.rm = TRUE),
      mean_alpha_logit = mean(d$alpha_logit, na.rm = TRUE),
      mean_alpha_logit_int = mean(d$alpha_logit_int, na.rm = TRUE),
      mean_alpha_spline = mean(d$alpha_spline, na.rm = TRUE),
      mean_alpha_rf = mean(d$alpha_rf, na.rm = TRUE),
      mean_alpha_gbm = mean(d$alpha_gbm, na.rm = TRUE)
    )
  })
  sm <- rbind_safe(out)
  if (nrow(sm) == 0L) return(sm)
  sm$rank_mse <- ave(sm$mse, sm$scenario_id, FUN = function(x) rank(x, ties.method = "min"))
  sm$best_by_mse <- sm$rank_mse == 1
  trad <- sm[sm$method == "trad", c("scenario_id", "mse", "rmse", "abs_bias")]
  names(trad) <- c("scenario_id", "mse_trad", "rmse_trad", "abs_bias_trad")
  sm <- merge(sm, trad, by = "scenario_id", all.x = TRUE)
  sm$rel_mse_vs_trad <- sm$mse / sm$mse_trad
  sm$rel_rmse_vs_trad <- sm$rmse / sm$rmse_trad
  sm$abs_bias_reduction_vs_trad <- sm$abs_bias_trad - sm$abs_bias
  sm
}

crit_pairs <- function(sm, method_new, method_ref) {
  if (nrow(sm) == 0L) return(data.frame())
  a <- sm[sm$method == method_new, c("scenario_id", "abs_bias", "mse", "rmse")]
  b <- sm[sm$method == method_ref, c("scenario_id", "abs_bias", "mse", "rmse")]
  if (nrow(a) == 0L || nrow(b) == 0L) return(data.frame())
  names(a) <- c("scenario_id", "abs_bias_new", "mse_new", "rmse_new")
  names(b) <- c("scenario_id", "abs_bias_ref", "mse_ref", "rmse_ref")
  m <- merge(a, b, by = "scenario_id")
  m$method_new <- method_new
  m$method_ref <- method_ref
  m$abs_bias_reduction <- m$abs_bias_ref - m$abs_bias_new
  m$mse_reduction <- m$mse_ref - m$mse_new
  m$rel_mse <- m$mse_new / m$mse_ref
  m
}

write_mc_summaries <- function(raw) {
  summary <- summarize_mc(raw)
  write_csv_atomic(summary, file.path(DIRS$tables, "monte_carlo_summary_by_method.csv"))
  if (nrow(raw) > 0L) {
    status <- aggregate(list(N = raw$status), by = list(block = raw$block, method = raw$method, status = raw$status), FUN = length)
    write_csv_atomic(status, file.path(DIRS$tables, "status_summary.csv"))
  }
  if (nrow(summary) > 0L) {
    win <- aggregate(list(n_wins = summary$best_by_mse), by = list(block = summary$block, method = summary$method), FUN = sum)
    count <- aggregate(list(n_scenarios = summary$scenario_id), by = list(block = summary$block, method = summary$method), FUN = length)
    win <- merge(win, count, by = c("block", "method"))
    win$win_rate <- win$n_wins / win$n_scenarios
    write_csv_atomic(win, file.path(DIRS$tables, "method_win_summary.csv"))

    blockA <- summary[summary$block == "A_independent", , drop = FALSE]
    blockB <- summary[summary$block == "B_covariate_dependent", , drop = FALSE]
    blockC <- summary[summary$block == "C_nonlinear_observability", , drop = FALSE]
    add_block <- function(d) {
      if (nrow(d) > 0L && !("block" %in% names(d))) d$block <- NA_character_
      d
    }
    decision <- rbind_safe(list(
      crit_pairs(blockA, "scap_x_l050", "trad"),
      crit_pairs(blockB, "ipcw_x", "ipcw_marg"),
      crit_pairs(blockB, "ps_x_spline", "ps_marg_logit"),
      crit_pairs(blockB, "scap_x_l050", "scap_marg_l050"),
      crit_pairs(blockB, "scap_x_l050", "trad"),
      crit_pairs(blockC, "ipcw_x", "ipcw_marg"),
      crit_pairs(blockC, "ps_x_spline", "ps_x_logit"),
      crit_pairs(blockC, "ps_x_spline", "ps_marg_logit"),
      crit_pairs(blockC, "scap_x_l050", "trad"),
      crit_pairs(blockC, "scap_x_l075", "trad"),
      crit_pairs(blockC, "scap_x_l100", "trad"),
      crit_pairs(blockC, "scap_x_l050", "scap_marg_l050")
    ))
    if (nrow(decision) > 0L) {
      # Recover block from summary, because crit_pairs works only at scenario_id level.
      lookup <- unique(summary[, c("scenario_id", "block")])
      decision <- merge(decision, lookup, by = "scenario_id", all.x = TRUE, sort = FALSE)
    }
    write_csv_atomic(decision, file.path(DIRS$tables, "decision_criteria_pairwise_comparisons.csv"))
    if (nrow(decision) > 0L) {
      groups <- split(decision, paste(decision$block, decision$method_new, decision$method_ref, sep = "||"))
      decision_summary <- rbind_safe(lapply(groups, function(d) {
        data.frame(
          block = d$block[1],
          method_new = d$method_new[1],
          method_ref = d$method_ref[1],
          n_scenarios = nrow(d),
          abs_bias_reduction_mean = mean(d$abs_bias_reduction, na.rm = TRUE),
          abs_bias_reduction_median = stats::median(d$abs_bias_reduction, na.rm = TRUE),
          abs_bias_reduction_prop_positive = mean(d$abs_bias_reduction > 0, na.rm = TRUE),
          mse_reduction_mean = mean(d$mse_reduction, na.rm = TRUE),
          mse_reduction_median = stats::median(d$mse_reduction, na.rm = TRUE),
          mse_reduction_prop_positive = mean(d$mse_reduction > 0, na.rm = TRUE),
          rel_mse_mean = mean(d$rel_mse, na.rm = TRUE),
          rel_mse_median = stats::median(d$rel_mse, na.rm = TRUE),
          rel_mse_prop_less_than_1 = mean(d$rel_mse < 1, na.rm = TRUE)
        )
      }))
      decision_summary <- decision_summary[order(decision_summary$block, decision_summary$method_new, decision_summary$method_ref), ]
      write_csv_atomic(decision_summary, file.path(DIRS$tables, "decision_criteria_summary.csv"))
    }
  }
  summary
}


# ==============================================================================
# 11. BOOTSTRAP COVERAGE STUDY
# ============================================================================== 

bootstrap_ci <- function(df, H, method, B, seed) {
  n <- nrow(df)
  vals <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    set.seed(make_seed(seed, b, 777))
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    db <- df[idx, , drop = FALSE]
    est <- tryCatch(estimate_method(db, H, method, seed = make_seed(seed, b, 888)), error = function(e) NULL)
    if (!is.null(est)) vals[b] <- as.numeric(est$tau_hat)
  }
  vals <- vals[is.finite(vals)]
  if (length(vals) < max(20L, B / 5L)) return(c(lower = NA_real_, upper = NA_real_, width = NA_real_, n_valid = length(vals)))
  qs <- stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, type = 7)
  c(lower = as.numeric(qs[1]), upper = as.numeric(qs[2]), width = as.numeric(qs[2] - qs[1]), n_valid = length(vals))
}

select_boot_scenarios <- function(scenarios) {
  s <- scenarios[scenarios$n == CONFIG$BOOT_N &
                   abs(scenarios$target_censor2 - CONFIG$BOOT_TARGET_CENSOR2) < 1e-8 &
                   abs(scenarios$H_quantile - CONFIG$BOOT_H_QUANTILE) < 1e-8, , drop = FALSE]
  if (nrow(s) == 0L) {
    scenarios$boot_dist <- abs(scenarios$n - CONFIG$BOOT_N) +
      1000 * abs(scenarios$target_censor2 - CONFIG$BOOT_TARGET_CENSOR2) +
      100 * abs(scenarios$H_quantile - CONFIG$BOOT_H_QUANTILE)
    s <- scenarios[order(scenarios$boot_dist), , drop = FALSE]
  }
  out <- list()
  for (bl in unique(s$block)) {
    sb <- s[s$block == bl, , drop = FALSE]
    if (bl == "B_covariate_dependent") sb <- sb[sb$cov_effect == "strong", , drop = FALSE]
    near0 <- sb[order(abs(sb$tau_design - 0)), , drop = FALSE][1, , drop = FALSE]
    pos <- sb[order(abs(sb$tau_design - 0.50)), , drop = FALSE][1, , drop = FALSE]
    out[[length(out) + 1L]] <- near0
    out[[length(out) + 1L]] <- pos
  }
  unique(rbind_safe(out))
}

run_bootstrap_dataset <- function(scenario, bds) {
  set.seed(make_seed(scenario$scenario_id, bds, 3333))
  df <- simulate_observed_data_x(scenario$n, scenario)
  df$x_type <- scenario$x_type
  rows <- list()
  for (method in CONFIG$BOOT_METHODS) {
    point <- tryCatch(estimate_method(df, scenario$H, method, seed = make_seed(scenario$scenario_id, bds, match(method, CONFIG$BOOT_METHODS), 44)), error = function(e) NULL)
    ci <- tryCatch(bootstrap_ci(df, scenario$H, method, CONFIG$N_BOOT_RESAMPLES,
                                seed = make_seed(scenario$scenario_id, bds, match(method, CONFIG$BOOT_METHODS), 55)),
                   error = function(e) c(lower = NA_real_, upper = NA_real_, width = NA_real_, n_valid = 0L))
    rows[[length(rows) + 1L]] <- data.frame(
      scenario_id = scenario$scenario_id,
      boot_dataset_id = bds,
      method = method,
      point = if (is.null(point)) NA_real_ else as.numeric(point$tau_hat),
      lower = ci["lower"],
      upper = ci["upper"],
      width = ci["width"],
      n_boot_valid = ci["n_valid"],
      covered = as.integer(is.finite(ci["lower"]) && ci["lower"] <= scenario$tau_H_true && scenario$tau_H_true <= ci["upper"])
    )
  }
  rbind_safe(rows)
}

run_bootstrap_for_scenario <- function(scenario) {
  out_file <- file.path(DIRS$bootstrap, sprintf("bootstrap_scenario_%04d_raw.csv", scenario$scenario_id))
  if (isTRUE(CONFIG$RESUME) && file.exists(out_file)) return(read.csv(out_file))
  message("Bootstrap scenario ", scenario$scenario_id, " with ", CONFIG$N_BOOT_DATASETS, " datasets; workers=", hpc_ncores())
  ids <- seq_len(CONFIG$N_BOOT_DATASETS)
  if (isTRUE(CONFIG$USE_PARALLEL) && hpc_ncores() > 1L && length(ids) > 1L) {
    cl <- parallel::makeCluster(hpc_ncores())
    try(parallel::clusterExport(cl, varlist = setdiff(ls(globalenv()), c("cl")), envir = globalenv()), silent = TRUE)
    pieces <- parallel::parLapply(cl, ids, function(bds) run_bootstrap_dataset(scenario, bds))
    try(parallel::stopCluster(cl), silent = TRUE)
  } else {
    pieces <- lapply(ids, function(bds) run_bootstrap_dataset(scenario, bds))
  }
  raw <- rbind_safe(pieces)
  raw <- merge(raw, scenario, by = "scenario_id", all.x = TRUE)
  write_csv_atomic(raw, out_file)
  raw
}

summarize_bootstrap <- function(boot_raw) {
  if (nrow(boot_raw) == 0L) return(data.frame())
  pieces <- split(boot_raw, interaction(boot_raw$scenario_id, boot_raw$method, drop = TRUE))
  rbind_safe(lapply(pieces, function(d) {
    data.frame(
      scenario_id = d$scenario_id[1], block = d$block[1], n = d$n[1], x_type = d$x_type[1],
      cov_effect = d$cov_effect[1], target_censor2 = d$target_censor2[1], H_quantile = d$H_quantile[1],
      tau_design = d$tau_design[1], tau_H_true = d$tau_H_true[1], method = d$method[1],
      n_boot_datasets = nrow(d), coverage = mean(d$covered, na.rm = TRUE),
      mean_width = mean(d$width, na.rm = TRUE), median_width = stats::median(d$width, na.rm = TRUE),
      mean_point = mean(d$point, na.rm = TRUE), bias = mean(d$point - d$tau_H_true, na.rm = TRUE),
      mse = mean((d$point - d$tau_H_true)^2, na.rm = TRUE), mean_n_boot_valid = mean(d$n_boot_valid, na.rm = TRUE)
    )
  }))
}

plot_with_gg <- function(summary_df, boot_df = NULL) {
  if (!HAS_GGPLOT2) return(FALSE)
  ggplot2 <- getNamespace("ggplot2")
  sm <- summary_df[summary_df$method %in% CONFIG$PLOT_METHODS, , drop = FALSE]
  sm$method <- factor(sm$method, levels = CONFIG$PLOT_METHODS)
  sm$block_label <- ifelse(sm$block == "A_independent", "A: independent censoring", "B: covariate-dependent censoring")

  p1 <- ggplot2$ggplot(sm, ggplot2$aes(x = factor(tau_design), y = bias, group = method, colour = method)) +
    ggplot2$geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2$geom_point() + ggplot2$geom_line() +
    ggplot2$facet_grid(block_label ~ target_censor2 + n, scales = "free_y") +
    ggplot2$labs(x = "Designed full Kendall tau", y = "Bias for tau_H", title = "Bias by block, sample size and censoring") +
    ggplot2$theme_bw()
  ggplot2$ggsave(file.path(DIRS$figures, "bias_by_block_tau_n_censor.png"), p1, width = 14, height = 8, dpi = 180)

  p2 <- ggplot2$ggplot(sm, ggplot2$aes(x = factor(tau_design), y = mse, group = method, colour = method)) +
    ggplot2$geom_point() + ggplot2$geom_line() +
    ggplot2$facet_grid(block_label ~ target_censor2 + n, scales = "free_y") +
    ggplot2$labs(x = "Designed full Kendall tau", y = "MSE for tau_H", title = "MSE by block, sample size and censoring") +
    ggplot2$theme_bw()
  ggplot2$ggsave(file.path(DIRS$figures, "mse_by_block_tau_n_censor.png"), p2, width = 14, height = 8, dpi = 180)

  p3 <- ggplot2$ggplot(sm, ggplot2$aes(x = factor(tau_design), y = rel_mse_vs_trad, group = method, colour = method)) +
    ggplot2$geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2$geom_point() + ggplot2$geom_line() +
    ggplot2$facet_grid(block_label ~ target_censor2 + n, scales = "free_y") +
    ggplot2$labs(x = "Designed full Kendall tau", y = "Relative MSE vs traditional", title = "Relative MSE vs traditional product-limit estimator") +
    ggplot2$theme_bw()
  ggplot2$ggsave(file.path(DIRS$figures, "relative_mse_vs_trad.png"), p3, width = 14, height = 8, dpi = 180)

  lam <- summary_df[grepl("scap_x_l", summary_df$method), , drop = FALSE]
  if (nrow(lam) > 0) {
    lam$lambda <- ifelse(grepl("l025", lam$method), 0.25, ifelse(grepl("l050", lam$method), 0.50, ifelse(grepl("l075", lam$method), 0.75, 1.00)))
    lam$block_label <- ifelse(lam$block == "A_independent", "A: independent censoring", "B: covariate-dependent censoring")
    p4 <- ggplot2$ggplot(lam, ggplot2$aes(x = lambda, y = mse, group = interaction(scenario_id), colour = factor(tau_design))) +
      ggplot2$geom_point(alpha = 0.7) + ggplot2$geom_line(alpha = 0.5) +
      ggplot2$facet_grid(block_label ~ target_censor2 + n, scales = "free_y") +
      ggplot2$labs(x = "Presmoothing intensity lambda", y = "MSE", colour = "tau design", title = "sCAP(lambda) sensitivity") +
      ggplot2$theme_bw()
    ggplot2$ggsave(file.path(DIRS$figures, "scap_lambda_sensitivity_mse.png"), p4, width = 14, height = 8, dpi = 180)
  }

  cal <- sm[is.finite(sm$mean_brier), , drop = FALSE]
  if (nrow(cal) > 0) {
    p5 <- ggplot2$ggplot(cal, ggplot2$aes(x = factor(tau_design), y = mean_brier, group = method, colour = method)) +
      ggplot2$geom_point() + ggplot2$geom_line() +
      ggplot2$facet_grid(block_label ~ target_censor2 + n, scales = "free_y") +
      ggplot2$labs(x = "Designed full Kendall tau", y = "Brier score for m", title = "Presmoothing probability calibration: Brier score") +
      ggplot2$theme_bw()
    ggplot2$ggsave(file.path(DIRS$figures, "presmoothing_brier_by_block.png"), p5, width = 14, height = 8, dpi = 180)
  }

  if (!is.null(boot_df) && nrow(boot_df) > 0) {
    boot_df$method <- factor(boot_df$method, levels = CONFIG$BOOT_METHODS)
    boot_df$block_label <- ifelse(boot_df$block == "A_independent", "A: independent censoring", "B: covariate-dependent censoring")
    p6 <- ggplot2$ggplot(boot_df, ggplot2$aes(x = method, y = coverage, fill = method)) +
      ggplot2$geom_hline(yintercept = 0.95, linetype = "dashed") +
      ggplot2$geom_col() +
      ggplot2$facet_wrap(~ block_label + tau_design + x_type, scales = "free_x") +
      ggplot2$labs(x = "Method", y = "Bootstrap coverage", title = "Percentile bootstrap coverage") +
      ggplot2$theme_bw() + ggplot2$theme(axis.text.x = ggplot2$element_text(angle = 45, hjust = 1))
    ggplot2$ggsave(file.path(DIRS$figures, "bootstrap_coverage.png"), p6, width = 12, height = 7, dpi = 180)
  }
  TRUE
}

plot_with_base <- function(summary_df, boot_df = NULL) {
  png(file.path(DIRS$figures, "baseplot_mse_by_method.png"), width = 1400, height = 800)
  sm <- summary_df[summary_df$method %in% CONFIG$PLOT_METHODS, , drop = FALSE]
  boxplot(mse ~ method + block, data = sm, las = 2, main = "MSE by method and block", ylab = "MSE")
  dev.off()

  png(file.path(DIRS$figures, "baseplot_bias_by_method.png"), width = 1400, height = 800)
  boxplot(bias ~ method + block, data = sm, las = 2, main = "Bias by method and block", ylab = "Bias")
  abline(h = 0, lty = 2)
  dev.off()
  TRUE
}

# ==============================================================================
# 13. MODE-DRIVEN EXECUTION FOR WORKSTATION OR DEUCALION/SLURM
# ============================================================================== 

run_truth_mode <- function() {
  scenarios <- compute_truth_all()
  message("Truth mode completed. Scenario file: ", scenario_grid_path())
  invisible(scenarios)
}

run_mc_mode <- function() {
  scenarios <- ensure_scenarios_loaded()
  idx <- task_indices(nrow(scenarios), CONFIG$N_SCENARIO_TASKS, CONFIG$TASK_ID)
  message("MC mode: task ", CONFIG$TASK_ID, "/", CONFIG$N_SCENARIO_TASKS,
          " will process scenario indices: ", paste(idx, collapse = ","))
  if (length(idx) == 0L) return(invisible(data.frame()))
  out <- lapply(idx, function(i) simulate_scenario(scenarios[i, , drop = FALSE]))
  raw <- rbind_safe(out)
  task_file <- file.path(DIRS$checkpoints, sprintf("mc_task_%04d_of_%04d_done.txt", CONFIG$TASK_ID, CONFIG$N_SCENARIO_TASKS))
  writeLines(c(paste("completed", Sys.time()), paste("scenarios", paste(scenarios$scenario_id[idx], collapse = ","))), task_file)
  invisible(raw)
}

collect_raw_files <- function() {
  files <- list.files(DIRS$raw, pattern = "^scenario_[0-9]+_raw\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No raw scenario files found in ", DIRS$raw)
  raw <- rbind_safe(lapply(files, read.csv))
  write_csv_atomic(raw, file.path(DIRS$tables, "monte_carlo_raw_all.csv"))
  saveRDS(raw, file.path(DIRS$root, "monte_carlo_raw_all.rds"))
  raw
}

run_bootstrap_mode <- function() {
  if (!isTRUE(CONFIG$RUN_BOOTSTRAP)) {
    message("Bootstrap disabled by CONFIG$RUN_BOOTSTRAP.")
    return(invisible(data.frame()))
  }
  scenarios <- ensure_scenarios_loaded()
  boot_scen <- select_boot_scenarios(scenarios)
  write_csv_atomic(boot_scen, file.path(DIRS$bootstrap, "bootstrap_selected_scenarios.csv"))
  idx <- task_indices(nrow(boot_scen), CONFIG$N_BOOT_TASKS, CONFIG$TASK_ID)
  message("Bootstrap mode: task ", CONFIG$TASK_ID, "/", CONFIG$N_BOOT_TASKS,
          " will process bootstrap scenario indices: ", paste(idx, collapse = ","))
  if (length(idx) == 0L) return(invisible(data.frame()))
  out <- lapply(idx, function(i) run_bootstrap_for_scenario(boot_scen[i, , drop = FALSE]))
  raw <- rbind_safe(out)
  task_file <- file.path(DIRS$checkpoints, sprintf("bootstrap_task_%04d_of_%04d_done.txt", CONFIG$TASK_ID, CONFIG$N_BOOT_TASKS))
  writeLines(c(paste("completed", Sys.time()), paste("boot_scenarios", paste(boot_scen$scenario_id[idx], collapse = ","))), task_file)
  invisible(raw)
}

collect_bootstrap_files <- function() {
  files <- list.files(DIRS$bootstrap, pattern = "^bootstrap_scenario_[0-9]+_raw\\.csv$", full.names = TRUE)
  if (length(files) == 0L) return(data.frame())
  boot_raw <- rbind_safe(lapply(files, read.csv))
  write_csv_atomic(boot_raw, file.path(DIRS$bootstrap, "bootstrap_raw_all.csv"))
  boot_summary <- summarize_bootstrap(boot_raw)
  write_csv_atomic(boot_summary, file.path(DIRS$bootstrap, "bootstrap_summary.csv"))
  boot_summary
}

write_manifest_and_log <- function(raw = NULL, summary = NULL, boot_summary = NULL) {
  manifest <- data.frame(path = list.files(OUT_DIR, recursive = TRUE, full.names = TRUE), stringsAsFactors = FALSE)
  if (nrow(manifest) > 0L) {
    manifest$relative_path <- sub(paste0("^", normalizePath(OUT_DIR), "/?"), "", normalizePath(manifest$path, winslash = "/"))
    manifest$size_bytes <- file.info(manifest$path)$size
    write_csv_atomic(manifest, file.path(DIRS$root, "manifest.csv"))
  }
  writeLines(c(
    paste("Completed:", as.character(Sys.time())),
    paste("Run mode:", CONFIG$RUN_MODE),
    paste("Run ID:", RUN_ID),
    paste("Output directory:", OUT_DIR),
    paste("Max cores per R process:", CONFIG$MAX_CORES),
    paste("Scenario tasks:", CONFIG$N_SCENARIO_TASKS),
    paste("Bootstrap tasks:", CONFIG$N_BOOT_TASKS),
    if (!is.null(raw)) paste("Raw rows:", nrow(raw)) else "Raw rows: NA",
    if (!is.null(summary)) paste("Summary rows:", nrow(summary)) else "Summary rows: NA",
    if (!is.null(boot_summary)) paste("Bootstrap summary rows:", nrow(boot_summary)) else "Bootstrap summary rows: NA"
  ), file.path(DIRS$logs, paste0("completion_log_", CONFIG$RUN_MODE, ".txt")))
}

run_collect_mode <- function() {
  raw <- collect_raw_files()
  summary <- write_mc_summaries(raw)
  boot_summary <- collect_bootstrap_files()
  if (!plot_with_gg(summary, boot_summary)) plot_with_base(summary, boot_summary)
  write_manifest_and_log(raw, summary, boot_summary)
  message("Collect mode completed. Results saved in: ", OUT_DIR)
  invisible(summary)
}

run_all_mode <- function() {
  run_truth_mode()
  # In all mode, use all scenarios in one process. This is intended for workstation/test runs.
  old_tasks <- CONFIG$N_SCENARIO_TASKS; old_task <- CONFIG$TASK_ID
  CONFIG$N_SCENARIO_TASKS <<- 1L; CONFIG$TASK_ID <<- 1L
  run_mc_mode()
  if (isTRUE(CONFIG$RUN_BOOTSTRAP)) {
    CONFIG$N_BOOT_TASKS <<- 1L; CONFIG$TASK_ID <<- 1L
    run_bootstrap_mode()
  }
  CONFIG$N_SCENARIO_TASKS <<- old_tasks; CONFIG$TASK_ID <<- old_task
  run_collect_mode()
}

mode <- tolower(CONFIG$RUN_MODE)
message("Starting tau_H simulation. mode=", mode, ", run_id=", RUN_ID, ", out=", OUT_DIR,
        ", max_cores=", CONFIG$MAX_CORES, ", task_id=", CONFIG$TASK_ID)

if (mode %in% c("truth", "init")) {
  run_truth_mode()
  write_manifest_and_log()
} else if (mode %in% c("mc", "montecarlo", "simulate")) {
  run_mc_mode()
  write_manifest_and_log()
} else if (mode %in% c("bootstrap", "boot")) {
  run_bootstrap_mode()
  write_manifest_and_log()
} else if (mode %in% c("collect", "summary", "figures")) {
  run_collect_mode()
} else if (mode %in% c("all", "workstation")) {
  run_all_mode()
} else {
  stop("Unknown TAUH_MODE/CONFIG$RUN_MODE: ", CONFIG$RUN_MODE)
}
