#!/usr/bin/env Rscript
# Real-data application for covariate-adaptive truncated Kendall association
# Dataset: survival::kidney, aggregated to subject-level sequential gap times.
# Outputs: estimates, bootstrap CIs, H/lambda sensitivity, nuisance diagnostics, plots, LaTeX tables.

options(warn = 1)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
CONFIG <- list(
  OUTPUT_ROOT = Sys.getenv("TAUH_REAL_OUTPUT_ROOT", "realdata_kidney_tauH_results"),
  RUN_ID = Sys.getenv("TAUH_REAL_RUN_ID", paste0("kidney_realdata_", format(Sys.time(), "%Y%m%d_%H%M%S"))),
  SEED = as.integer(Sys.getenv("TAUH_REAL_SEED", "20260531")),
  N_BOOT = as.integer(Sys.getenv("TAUH_REAL_N_BOOT", "2000")),
  MAX_CORES = as.integer(Sys.getenv("TAUH_REAL_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "8"))),
  H_QUANTILES = c(0.70, 0.80, 0.90, 0.95),
  H_MAIN_Q = as.numeric(Sys.getenv("TAUH_REAL_H_MAIN_Q", "0.90")),
  LAMBDA_GRID = seq(0, 1, by = 0.10),
  MAIN_LAMBDAS = c(0.50, 0.75),
  K_FOLDS = as.integer(Sys.getenv("TAUH_REAL_K_FOLDS", "5")),
  EPS = 1e-6,
  IPCW_TRUNC = as.numeric(Sys.getenv("TAUH_REAL_IPCW_TRUNC", "0.02"))
)
CONFIG$MAX_CORES <- max(1L, min(64L, CONFIG$MAX_CORES))
set.seed(CONFIG$SEED)

OUT <- file.path(CONFIG$OUTPUT_ROOT, CONFIG$RUN_ID)
DIRS <- list(
  root = OUT,
  tables = file.path(OUT, "tables"),
  figures = file.path(OUT, "figures"),
  bootstrap = file.path(OUT, "bootstrap"),
  latex = file.path(OUT, "latex"),
  logs = file.path(OUT, "logs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(DIRS$logs, "run_log.txt"), append = TRUE)
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE)
clip01 <- function(p, eps = CONFIG$EPS) pmin(pmax(as.numeric(p), eps), 1 - eps)
safe_div <- function(a, b) ifelse(abs(b) < 1e-12, NA_real_, a / b)

# -----------------------------------------------------------------------------
# Package loading
# -----------------------------------------------------------------------------
need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not available: ", pkg,
         ". Install it before running this application script.")
  }
}
need_pkg("survival")
need_pkg("splines")
HAS_GGPLOT <- requireNamespace("ggplot2", quietly = TRUE)
HAS_PARALLEL <- requireNamespace("parallel", quietly = TRUE)

# -----------------------------------------------------------------------------
# Data construction
# -----------------------------------------------------------------------------
load_kidney_gap_data <- function() {
  # Robust loading strategy:
  # 1. Try survival::kidney directly.
  # 2. Try data(cancer, package="survival"), which may load kidney in some survival versions.
  # 3. Fall back to a local CSV file kidney_survival.csv.
  env <- new.env(parent = emptyenv())

  raw <- NULL

  try(suppressWarnings(data("kidney", package = "survival", envir = env)), silent = TRUE)

  if (!exists("kidney", envir = env, inherits = FALSE)) {
    try(suppressWarnings(data("cancer", package = "survival", envir = env)), silent = TRUE)
  }

  if (exists("kidney", envir = env, inherits = FALSE)) {
    raw <- get("kidney", envir = env)
    message("Loaded kidney data from survival package.")
  }

  if (is.null(raw)) {
    local_csv <- file.path(getwd(), "kidney_survival.csv")
    if (!file.exists(local_csv)) {
      stop(
        "Could not load kidney data from survival and local file not found: ",
        local_csv,
        ". Create it with: curl -L -o kidney_survival.csv ",
        "https://vincentarelbundock.github.io/Rdatasets/csv/survival/kidney.csv"
      )
    }
    raw <- utils::read.csv(local_csv, stringsAsFactors = FALSE)
    message("Loaded kidney data from local CSV: ", local_csv)
  }

  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  # Harmonize column names across sources.
  if ("rownames" %in% names(raw)) raw$rownames <- NULL
  if ("patient" %in% names(raw) && !("id" %in% names(raw))) {
    names(raw)[names(raw) == "patient"] <- "id"
  }

  if (!("id" %in% names(raw))) stop("Kidney data do not contain an id column.")
  if (!("time" %in% names(raw))) stop("Kidney data do not contain a time column.")
  if (!("status" %in% names(raw))) stop("Kidney data do not contain a status column.")

  raw$id <- as.character(raw$id)
  raw$time <- as.numeric(raw$time)
  raw$status <- as.integer(raw$status)

  raw$.__row <- seq_len(nrow(raw))
  raw <- raw[order(raw$id, raw$.__row), , drop = FALSE]
  groups <- split(raw, raw$id)

  rows <- lapply(groups, function(g) {
    g <- g[order(g$.__row), , drop = FALSE]
    r1 <- g[1, , drop = FALSE]
    r2 <- if (nrow(g) >= 2) g[2, , drop = FALSE] else NULL

    delta1 <- as.integer(r1$status[1] == 1)
    A_tilde <- as.numeric(r1$time[1])

    if (delta1 == 1L && !is.null(r2)) {
      B_tilde <- as.numeric(r2$time[1])
      delta2 <- as.integer(r2$status[1] == 1)
      Z_tilde <- A_tilde + B_tilde
      A_for_tau <- A_tilde
      B_for_tau <- B_tilde
    } else if (delta1 == 1L && is.null(r2)) {
      B_tilde <- 0
      delta2 <- 0L
      Z_tilde <- A_tilde
      A_for_tau <- A_tilde
      B_for_tau <- 0
    } else {
      B_tilde <- 0
      delta2 <- 0L
      Z_tilde <- A_tilde
      A_for_tau <- 0
      B_for_tau <- 0
    }

    data.frame(
      id = as.character(r1$id[1]),
      A_tilde = A_tilde,
      B_tilde = B_tilde,
      Z_tilde = Z_tilde,
      Delta1 = delta1,
      Delta2 = delta2,
      A = A_for_tau,
      B = B_for_tau,
      n_rows_original = nrow(g),
      age = if ("age" %in% names(g)) as.numeric(r1$age[1]) else NA_real_,
      sex = if ("sex" %in% names(g)) as.character(r1$sex[1]) else NA_character_,
      disease = if ("disease" %in% names(g)) as.character(r1$disease[1]) else NA_character_,
      frail = if ("frail" %in% names(g)) as.numeric(r1$frail[1]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  dat <- do.call(rbind, rows)
  dat$sex[is.na(dat$sex) | dat$sex == ""] <- "Missing"
  dat$disease[is.na(dat$disease) | dat$disease == ""] <- "Missing"
  dat$sex <- factor(dat$sex)
  dat$disease <- factor(dat$disease)
  dat
}

add_covariate_design <- function(dat) {
  dat$age_z <- if (all(is.na(dat$age)) || sd(dat$age, na.rm = TRUE) == 0) 0 else as.numeric(scale(dat$age))
  mm <- model.matrix(~ sex + disease, data = dat)
  if (ncol(mm) > 1) {
    mm <- mm[, -1, drop = FALSE]
    colnames(mm) <- make.names(colnames(mm), unique = TRUE)
    dat <- cbind(dat, as.data.frame(mm))
  }
  dat
}

get_covariate_cols <- function(dat) {
  base <- c("age_z", grep("^(sex|disease)", names(dat), value = TRUE))
  base[base %in% names(dat)]
}

# -----------------------------------------------------------------------------
# Core estimators
# -----------------------------------------------------------------------------
rank_z_censored_last <- function(z, delta2) {
  n <- length(z)
  ord <- order(z, -as.integer(delta2 == 1), seq_len(n), na.last = TRUE)
  r <- integer(n)
  r[ord] <- seq_len(n)
  r
}

pl_weights <- function(z, event_score, delta2) {
  n <- length(z)
  event_score <- pmin(pmax(as.numeric(event_score), 0), 1)
  r <- rank_z_censored_last(z, delta2)
  ord <- order(r)
  w <- numeric(n)
  surv <- 1
  for (ii in ord) {
    denom <- n - r[ii] + 1
    haz <- if (denom > 0) event_score[ii] / denom else 0
    haz <- pmin(pmax(haz, 0), 1)
    w[ii] <- surv * haz
    surv <- surv * (1 - haz)
  }
  w
}

psi_matrix <- function(A, B) {
  DA <- outer(A, A, "-")
  DB <- outer(B, B, "-")
  sign(DA * DB)
}

tau_from_weights <- function(A, B, weights, H) {
  I <- as.numeric(A + B <= H)
  ww <- weights * I
  M <- sum(ww, na.rm = TRUE)
  if (!is.finite(M) || M <= 1e-12) {
    return(list(tau = NA_real_, M_hat = M, ess = NA_real_, sum_weights = sum(weights, na.rm = TRUE)))
  }
  P <- psi_matrix(A, B)
  tau <- sum(outer(ww, ww) * P, na.rm = TRUE) / (M * M)
  ess <- if (sum(ww^2) > 0) M^2 / sum(ww^2) else NA_real_
  list(tau = tau, M_hat = M, ess = ess, sum_weights = sum(weights, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# Cross-fitted observation-probability models
# -----------------------------------------------------------------------------
make_folds <- function(n, K, seed) {
  set.seed(seed)
  K <- max(2L, min(K, n))
  sample(rep(seq_len(K), length.out = n))
}

fit_predict_glm <- function(formula, train, test, fallback) {
  fit <- tryCatch(
    suppressWarnings(glm(formula, data = train, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(fallback, nrow(test)))
  p <- tryCatch(as.numeric(predict(fit, newdata = test, type = "response")), error = function(e) rep(fallback, nrow(test)))
  p[!is.finite(p)] <- fallback
  clip01(p)
}

crossfit_mhat <- function(dat, learner = c("constant", "marg_logit", "x_logit", "x_spline"),
                          K = CONFIG$K_FOLDS, seed = CONFIG$SEED + 17) {
  learner <- match.arg(learner)
  idx <- which(dat$Delta1 == 1)
  out <- rep(NA_real_, nrow(dat))
  if (length(idx) == 0) return(rep(0, nrow(dat)))
  y <- dat$Delta2[idx]
  global_mean <- mean(y, na.rm = TRUE)
  global_mean <- clip01(ifelse(is.finite(global_mean), global_mean, 0.5))
  if (learner == "constant" || length(unique(y)) < 2 || length(idx) < 6) {
    out[idx] <- global_mean
    out[is.na(out)] <- 0
    return(clip01(out))
  }
  folds <- make_folds(length(idx), K, seed)
  cov_cols <- get_covariate_cols(dat)

  for (k in sort(unique(folds))) {
    test_idx <- idx[folds == k]
    train_idx <- idx[folds != k]
    train <- dat[train_idx, , drop = FALSE]
    test <- dat[test_idx, , drop = FALSE]
    fallback <- clip01(mean(train$Delta2, na.rm = TRUE))
    if (!is.finite(fallback)) fallback <- global_mean

    if (learner == "marg_logit") {
      form <- Delta2 ~ A_tilde + Z_tilde
    } else if (learner == "x_logit") {
      rhs <- paste(c("A_tilde", "Z_tilde", cov_cols), collapse = " + ")
      form <- as.formula(paste("Delta2 ~", rhs))
    } else {
      rhs <- c("splines::ns(A_tilde, df = 2)", "splines::ns(Z_tilde, df = 2)")
      if ("age_z" %in% cov_cols) rhs <- c(rhs, "splines::ns(age_z, df = 2)")
      rhs <- c(rhs, setdiff(cov_cols, "age_z"))
      form <- as.formula(paste("Delta2 ~", paste(rhs, collapse = " + ")))
    }
    out[test_idx] <- fit_predict_glm(form, train, test, fallback)
  }
  out[is.na(out) & dat$Delta1 == 1] <- global_mean
  out[is.na(out)] <- 0
  clip01(out)
}

mhat_library <- function(dat, K = CONFIG$K_FOLDS, seed = CONFIG$SEED + 19) {
  learners <- c("constant", "marg_logit", "x_logit", "x_spline")
  preds <- lapply(seq_along(learners), function(j) crossfit_mhat(dat, learners[j], K, seed + 101 * j))
  names(preds) <- learners
  idx <- which(dat$Delta1 == 1)
  diag <- lapply(learners, function(l) {
    p <- preds[[l]][idx]
    y <- dat$Delta2[idx]
    brier <- mean((y - p)^2, na.rm = TRUE)
    logloss <- -mean(y * log(clip01(p)) + (1 - y) * log(clip01(1 - p)), na.rm = TRUE)
    cal <- calibration_stats(y, p)
    data.frame(learner = l, brier = brier, logloss = logloss,
               cal_intercept = cal$intercept, cal_slope = cal$slope, ece = cal$ece)
  })
  diag <- do.call(rbind, diag)
  inv <- 1 / pmax(diag$brier, 1e-8)
  alpha <- inv / sum(inv)
  names(alpha) <- learners
  ens <- Reduce("+", Map(function(p, a) p * a, preds, alpha))
  diag$alpha <- alpha[diag$learner]
  list(preds = preds, ensemble = clip01(ens), diagnostics = diag)
}

calibration_stats <- function(y, p) {
  p <- clip01(p)
  if (length(unique(y)) < 2) return(list(intercept = NA_real_, slope = NA_real_, ece = NA_real_))
  lp <- qlogis(p)
  fit <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial())), error = function(e) NULL)
  intercept <- if (!is.null(fit) && length(coef(fit)) >= 1) coef(fit)[1] else NA_real_
  slope <- if (!is.null(fit) && length(coef(fit)) >= 2) coef(fit)[2] else NA_real_
  br <- tryCatch({
    cuts <- unique(quantile(p, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
    if (length(cuts) <= 2) {
      NA_real_
    } else {
      bins <- cut(p, breaks = cuts, include.lowest = TRUE)
      tab <- aggregate(data.frame(y = y, p = p), list(bin = bins), mean)
      w <- as.numeric(table(bins)) / length(bins)
      sum(w * abs(tab$y - tab$p))
    }
  }, error = function(e) NA_real_)
  list(intercept = as.numeric(intercept), slope = as.numeric(slope), ece = as.numeric(br))
}

# -----------------------------------------------------------------------------
# IPCW censoring survival estimators
# -----------------------------------------------------------------------------
km_surv_at <- function(time, event, t_eval) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  ss <- summary(fit, times = t_eval, extend = TRUE)
  surv <- as.numeric(ss$surv)
  surv[!is.finite(surv)] <- 1
  pmax(surv, CONFIG$IPCW_TRUNC)
}

cox_censor_surv <- function(dat) {
  event_c <- as.integer(dat$Delta2 == 0)
  cov_cols <- get_covariate_cols(dat)
  Gm <- km_surv_at(dat$Z_tilde, event_c, dat$Z_tilde)
  if (length(cov_cols) == 0 || sum(event_c) < 5 || length(unique(event_c)) < 2) return(Gm)
  rhs <- paste(cov_cols, collapse = " + ")
  form <- as.formula(paste("survival::Surv(Z_tilde, event_c) ~", rhs))
  dd <- dat
  dd$event_c <- event_c
  fit <- tryCatch(suppressWarnings(survival::coxph(form, data = dd, ties = "breslow")), error = function(e) NULL)
  if (is.null(fit)) return(Gm)
  bh <- tryCatch(survival::basehaz(fit, centered = FALSE), error = function(e) NULL)
  lp <- tryCatch(as.numeric(predict(fit, newdata = dd, type = "lp")), error = function(e) rep(0, nrow(dd)))
  if (is.null(bh) || !all(is.finite(lp))) return(Gm)
  H0 <- approx(x = bh$time, y = bh$hazard, xout = dd$Z_tilde, method = "constant", f = 0,
               rule = 2, yleft = 0)$y
  G <- exp(-H0 * exp(lp))
  G[!is.finite(G)] <- Gm[!is.finite(G)]
  pmax(G, CONFIG$IPCW_TRUNC)
}

# -----------------------------------------------------------------------------
# Estimation wrapper
# -----------------------------------------------------------------------------
estimate_methods <- function(dat, H, seed = CONFIG$SEED, compute_diag = TRUE) {
  lib <- mhat_library(dat, CONFIG$K_FOLDS, seed)
  m_const <- lib$preds$constant
  m_marg <- lib$preds$marg_logit
  m_xlogit <- lib$preds$x_logit
  m_xspline <- lib$preds$x_spline
  m_sl <- lib$ensemble

  methods <- list()
  methods$trad <- pl_weights(dat$Z_tilde, dat$Delta2, dat$Delta2)
  methods$ps_marg_logit <- pl_weights(dat$Z_tilde, dat$Delta1 * m_marg, dat$Delta2)
  methods$ps_x_logit <- pl_weights(dat$Z_tilde, dat$Delta1 * m_xlogit, dat$Delta2)
  methods$ps_x_spline <- pl_weights(dat$Z_tilde, dat$Delta1 * m_xspline, dat$Delta2)
  for (lam in c(CONFIG$MAIN_LAMBDAS, 1.0)) {
    nm <- paste0("scap_x_l", sprintf("%03.0f", 100 * lam))
    escore <- dat$Delta1 * ((1 - lam) * dat$Delta2 + lam * m_sl)
    methods[[nm]] <- pl_weights(dat$Z_tilde, escore, dat$Delta2)
  }

  Gm <- km_surv_at(dat$Z_tilde, as.integer(dat$Delta2 == 0), dat$Z_tilde)
  Gx <- cox_censor_surv(dat)
  I <- as.numeric(dat$A + dat$B <= H)
  v_m <- dat$Delta2 * I / Gm
  v_x <- dat$Delta2 * I / Gx

  res <- lapply(names(methods), function(nm) {
    z <- tau_from_weights(dat$A, dat$B, methods[[nm]], H)
    data.frame(method = nm, tau_H = z$tau, M_hat = z$M_hat, ess = z$ess,
               sum_weights = z$sum_weights, H = H)
  })

  tau_ipcw <- function(v) {
    M <- sum(v, na.rm = TRUE)
    if (!is.finite(M) || M <= 1e-12) return(list(tau = NA_real_, M_hat = NA_real_, ess = NA_real_))
    P <- psi_matrix(dat$A, dat$B)
    tau <- sum(outer(v, v) * P, na.rm = TRUE) / (M * M)
    ess <- if (sum(v^2) > 0) M^2 / sum(v^2) else NA_real_
    list(tau = tau, M_hat = mean(v, na.rm = TRUE), ess = ess)
  }
  zi <- tau_ipcw(v_m)
  res[[length(res) + 1]] <- data.frame(method = "ipcw_marg", tau_H = zi$tau, M_hat = zi$M_hat,
                                       ess = zi$ess, sum_weights = sum(v_m), H = H)
  zi <- tau_ipcw(v_x)
  res[[length(res) + 1]] <- data.frame(method = "ipcw_x", tau_H = zi$tau, M_hat = zi$M_hat,
                                       ess = zi$ess, sum_weights = sum(v_x), H = H)
  est <- do.call(rbind, res)

  if (compute_diag) {
    idx <- dat$Delta1 == 1
    nuis <- lib$diagnostics
    nuis$n_train_observable <- sum(idx)
    nuis$event_rate_Delta2_given_Delta1 <- mean(dat$Delta2[idx], na.rm = TRUE)
    wdiag <- data.frame(
      method = est$method,
      H = H,
      M_hat = est$M_hat,
      ess = est$ess,
      sum_weights = est$sum_weights
    )
    return(list(estimates = est, nuisance = nuis, weights = wdiag,
                mhat = data.frame(id = dat$id, Delta1 = dat$Delta1, Delta2 = dat$Delta2,
                                  A_tilde = dat$A_tilde, Z_tilde = dat$Z_tilde,
                                  m_constant = m_const, m_marg_logit = m_marg,
                                  m_x_logit = m_xlogit, m_x_spline = m_xspline,
                                  m_ensemble = m_sl,
                                  G_marg = Gm, G_x = Gx)))
  }
  list(estimates = est)
}


estimate_methods_fixed_nuisance <- function(dat, H, nuisance_df) {
  # Sensitivity analyses should vary only H or lambda, not the nuisance fits.
  # This function reuses the nuisance predictions fitted in the main analysis.

  if (nrow(nuisance_df) != nrow(dat)) {
    stop("Fixed nuisance data frame has incompatible number of rows.")
  }

  m_marg <- clip01(nuisance_df$m_marg_logit)
  m_xlogit <- clip01(nuisance_df$m_x_logit)
  m_xspline <- clip01(nuisance_df$m_x_spline)
  m_sl <- clip01(nuisance_df$m_ensemble)

  methods <- list()
  methods$trad <- pl_weights(dat$Z_tilde, dat$Delta2, dat$Delta2)
  methods$ps_marg_logit <- pl_weights(dat$Z_tilde, dat$Delta1 * m_marg, dat$Delta2)
  methods$ps_x_logit <- pl_weights(dat$Z_tilde, dat$Delta1 * m_xlogit, dat$Delta2)
  methods$ps_x_spline <- pl_weights(dat$Z_tilde, dat$Delta1 * m_xspline, dat$Delta2)

  for (lam in c(CONFIG$MAIN_LAMBDAS, 1.0)) {
    nm <- paste0("scap_x_l", sprintf("%03.0f", 100 * lam))
    escore <- dat$Delta1 * ((1 - lam) * dat$Delta2 + lam * m_sl)
    methods[[nm]] <- pl_weights(dat$Z_tilde, escore, dat$Delta2)
  }

  res <- lapply(names(methods), function(nm) {
    z <- tau_from_weights(dat$A, dat$B, methods[[nm]], H)
    data.frame(
      method = nm,
      tau_H = z$tau,
      M_hat = z$M_hat,
      ess = z$ess,
      sum_weights = z$sum_weights,
      H = H
    )
  })

  Gm <- pmax(as.numeric(nuisance_df$G_marg), CONFIG$IPCW_TRUNC)
  Gx <- pmax(as.numeric(nuisance_df$G_x), CONFIG$IPCW_TRUNC)

  I <- as.numeric(dat$A + dat$B <= H)
  v_m <- dat$Delta2 * I / Gm
  v_x <- dat$Delta2 * I / Gx

  tau_ipcw <- function(v) {
    M <- sum(v, na.rm = TRUE)
    if (!is.finite(M) || M <= 1e-12) {
      return(list(tau = NA_real_, M_hat = NA_real_, ess = NA_real_))
    }
    P <- psi_matrix(dat$A, dat$B)
    tau <- sum(outer(v, v) * P, na.rm = TRUE) / (M * M)
    ess <- if (sum(v^2) > 0) M^2 / sum(v^2) else NA_real_
    list(tau = tau, M_hat = mean(v, na.rm = TRUE), ess = ess)
  }

  zi <- tau_ipcw(v_m)
  res[[length(res) + 1L]] <- data.frame(
    method = "ipcw_marg",
    tau_H = zi$tau,
    M_hat = zi$M_hat,
    ess = zi$ess,
    sum_weights = sum(v_m, na.rm = TRUE),
    H = H
  )

  zi <- tau_ipcw(v_x)
  res[[length(res) + 1L]] <- data.frame(
    method = "ipcw_x",
    tau_H = zi$tau,
    M_hat = zi$M_hat,
    ess = zi$ess,
    sum_weights = sum(v_x, na.rm = TRUE),
    H = H
  )

  list(estimates = do.call(rbind, res))
}

estimate_lambda_grid_fixed_nuisance <- function(dat, H, nuisance_df, lambdas) {
  if (nrow(nuisance_df) != nrow(dat)) {
    stop("Fixed nuisance data frame has incompatible number of rows.")
  }

  m_sl <- clip01(nuisance_df$m_ensemble)

  out <- lapply(lambdas, function(lam) {
    escore <- dat$Delta1 * ((1 - lam) * dat$Delta2 + lam * m_sl)
    w <- pl_weights(dat$Z_tilde, escore, dat$Delta2)
    z <- tau_from_weights(dat$A, dat$B, w, H)
    data.frame(
      lambda = lam,
      tau_H = z$tau,
      M_hat = z$M_hat,
      ess = z$ess,
      H = H
    )
  })

  do.call(rbind, out)
}

estimate_lambda_grid <- function(dat, H, lambdas = CONFIG$LAMBDA_GRID, seed = CONFIG$SEED + 999) {
  lib <- mhat_library(dat, CONFIG$K_FOLDS, seed)
  m_sl <- lib$ensemble
  out <- lapply(lambdas, function(lam) {
    escore <- dat$Delta1 * ((1 - lam) * dat$Delta2 + lam * m_sl)
    w <- pl_weights(dat$Z_tilde, escore, dat$Delta2)
    z <- tau_from_weights(dat$A, dat$B, w, H)
    data.frame(lambda = lam, tau_H = z$tau, M_hat = z$M_hat, ess = z$ess, H = H)
  })
  do.call(rbind, out)
}

# -----------------------------------------------------------------------------
# Bootstrap
# -----------------------------------------------------------------------------
bootstrap_once <- function(b, dat, H_main, methods_keep, seed) {
  set.seed(seed + b)
  idx <- sample(seq_len(nrow(dat)), size = nrow(dat), replace = TRUE)
  boot <- dat[idx, , drop = FALSE]
  boot$id <- paste0(boot$id, "_b", seq_len(nrow(boot)))
  est <- tryCatch(estimate_methods(boot, H_main, seed = seed + b, compute_diag = FALSE)$estimates,
                  error = function(e) NULL)
  if (is.null(est)) {
    return(data.frame(boot_id = b, method = methods_keep, tau_H = NA_real_, status = "error"))
  }
  est <- est[est$method %in% methods_keep, c("method", "tau_H"), drop = FALSE]
  est$boot_id <- b
  est$status <- "ok"
  est
}

run_bootstrap <- function(dat, H_main, methods_keep, n_boot, cores, seed) {
  log_msg("Starting bootstrap with", n_boot, "replicates and", cores, "cores")
  if (.Platform$OS.type != "windows" && HAS_PARALLEL && cores > 1) {
    ans <- parallel::mclapply(seq_len(n_boot), bootstrap_once, dat = dat, H_main = H_main,
                              methods_keep = methods_keep, seed = seed,
                              mc.cores = cores, mc.preschedule = FALSE)
  } else {
    ans <- lapply(seq_len(n_boot), bootstrap_once, dat = dat, H_main = H_main,
                  methods_keep = methods_keep, seed = seed)
  }
  do.call(rbind, ans)
}

summarize_bootstrap <- function(point_est, boot_raw) {
  keep <- boot_raw[boot_raw$status == "ok" & is.finite(boot_raw$tau_H), , drop = FALSE]
  out <- lapply(split(keep, keep$method), function(d) {
    method <- d$method[1]
    theta <- point_est$tau_H[match(method, point_est$method)]
    se <- sd(d$tau_H, na.rm = TRUE)
    q <- quantile(d$tau_H, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
    data.frame(method = method,
               estimate = theta,
               boot_mean = mean(d$tau_H, na.rm = TRUE),
               boot_bias = mean(d$tau_H, na.rm = TRUE) - theta,
               boot_se = se,
               ci_perc_l = q[1],
               ci_perc_median = q[2],
               ci_perc_u = q[3],
               ci_norm_l = theta - 1.96 * se,
               ci_norm_u = theta + 1.96 * se,
               ci_width_perc = q[3] - q[1],
               n_boot_valid = nrow(d))
  })
  out <- do.call(rbind, out)
  out[order(out$method), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Plotting and LaTeX helpers
# -----------------------------------------------------------------------------
plot_outputs <- function(dat, main_summary, H_sens, lambda_sens, mhat) {
  if (!HAS_GGPLOT) {
    log_msg("ggplot2 not available; skipping plots")
    return(invisible(NULL))
  }
  ggplot2 <- asNamespace("ggplot2")
  pdat <- dat
  pdat$gap_status <- ifelse(pdat$Delta1 == 0, "First gap censored",
                            ifelse(pdat$Delta2 == 1, "Second gap observed", "Second gap censored"))
  p <- ggplot2$ggplot(pdat, ggplot2$aes(x = A, y = B, color = gap_status)) +
    ggplot2$geom_point(size = 2, alpha = 0.85) +
    ggplot2$labs(x = "First gap time", y = "Second gap time or censoring gap",
                 color = "Observed status") +
    ggplot2$theme_bw()
  ggplot2$ggsave(file.path(DIRS$figures, "kidney_gap_scatter.png"), p, width = 7, height = 5, dpi = 300)

  p <- ggplot2$ggplot(main_summary, ggplot2$aes(x = reorder(method, estimate), y = estimate)) +
    ggplot2$geom_point(size = 2) +
    ggplot2$geom_errorbar(ggplot2$aes(ymin = ci_perc_l, ymax = ci_perc_u), width = 0.2) +
    ggplot2$coord_flip() + ggplot2$theme_bw() +
    ggplot2$labs(x = "Method", y = "Estimated truncated Kendall association",
                 title = "Kidney data: point estimates and percentile bootstrap intervals")
  ggplot2$ggsave(file.path(DIRS$figures, "kidney_tauH_estimates_bootstrap_CI.png"), p, width = 7.5, height = 5.5, dpi = 300)

  p <- ggplot2$ggplot(H_sens, ggplot2$aes(x = H_quantile, y = tau_H, group = method, color = method)) +
    ggplot2$geom_line() + ggplot2$geom_point(size = 1.8) + ggplot2$theme_bw() +
    ggplot2$labs(x = "Quantile used to define H", y = "Estimated tau_H", color = "Method",
                 title = "Sensitivity of tau_H to truncation threshold H")
  ggplot2$ggsave(file.path(DIRS$figures, "kidney_tauH_sensitivity_H.png"), p, width = 8, height = 5, dpi = 300)

  p <- ggplot2$ggplot(lambda_sens, ggplot2$aes(x = lambda, y = tau_H)) +
    ggplot2$geom_line() + ggplot2$geom_point(size = 1.8) + ggplot2$theme_bw() +
    ggplot2$labs(x = "lambda", y = "Estimated tau_H", title = "Sensitivity of sCAP(lambda)")
  ggplot2$ggsave(file.path(DIRS$figures, "kidney_scap_lambda_sensitivity.png"), p, width = 7, height = 5, dpi = 300)

  calib <- mhat[mhat$Delta1 == 1, , drop = FALSE]
  p <- ggplot2$ggplot(calib, ggplot2$aes(x = m_ensemble, y = Delta2)) +
    ggplot2$geom_jitter(height = 0.03, width = 0, alpha = 0.5) +
    ggplot2$geom_smooth(method = "loess", se = TRUE) + ggplot2$theme_bw() +
    ggplot2$labs(x = "Cross-fitted ensemble m-hat", y = "Observed Delta2",
                 title = "Observation-probability calibration")
  ggplot2$ggsave(file.path(DIRS$figures, "kidney_mhat_calibration.png"), p, width = 7, height = 5, dpi = 300)

  longG <- data.frame(id = rep(mhat$id, 2),
                      estimator = rep(c("G_marg", "G_x"), each = nrow(mhat)),
                      G = c(mhat$G_marg, mhat$G_x), Delta2 = rep(mhat$Delta2, 2))
  p <- ggplot2$ggplot(longG, ggplot2$aes(x = estimator, y = 1 / pmax(G, CONFIG$IPCW_TRUNC))) +
    ggplot2$geom_boxplot() + ggplot2$geom_jitter(width = 0.08, alpha = 0.5) + ggplot2$theme_bw() +
    ggplot2$labs(x = "Censoring model", y = "Inverse censoring weight", title = "IPCW weight diagnostics")
  ggplot2$ggsave(file.path(DIRS$figures, "kidney_ipcw_weight_diagnostics.png"), p, width = 6.5, height = 5, dpi = 300)
}

latex_table <- function(df, path, digits = 3) {
  d <- df
  for (j in seq_along(d)) if (is.numeric(d[[j]])) d[[j]] <- formatC(d[[j]], format = "f", digits = digits)
  lines <- c("\\begin{tabular}{" %+% paste(rep("l", ncol(d)), collapse = "") %+% "}",
             "\\hline",
             paste(names(d), collapse = " & ") %+% " \\\\",
             "\\hline")
  for (i in seq_len(nrow(d))) lines <- c(lines, paste(as.character(d[i, ]), collapse = " & ") %+% " \\\\")
  lines <- c(lines, "\\hline", "\\end{tabular}")
  writeLines(lines, path)
}
`%+%` <- function(a, b) paste0(a, b)

# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------
log_msg("Starting kidney real-data application")
log_msg("Output directory:", OUT)
log_msg("Bootstrap replicates:", CONFIG$N_BOOT, "cores:", CONFIG$MAX_CORES)

kid <- add_covariate_design(load_kidney_gap_data())
write_csv(kid, file.path(DIRS$tables, "kidney_gap_data_subject_level.csv"))

summary_tab <- data.frame(
  quantity = c("subjects", "original rows represented", "Delta1 observed", "Delta2 observed", "first-gap censoring rate", "second-gap censoring rate among all", "second-gap censoring rate among Delta1=1"),
  value = c(nrow(kid), sum(kid$n_rows_original), sum(kid$Delta1 == 1), sum(kid$Delta2 == 1),
            mean(kid$Delta1 == 0), mean(kid$Delta2 == 0), mean(kid$Delta2[kid$Delta1 == 1] == 0))
)
write_csv(summary_tab, file.path(DIRS$tables, "kidney_data_summary.csv"))

# H grid based on observable total time among subjects with first gap observed
z_ref <- kid$Z_tilde[kid$Delta1 == 1]
H_vals <- as.numeric(quantile(z_ref, probs = CONFIG$H_QUANTILES, na.rm = TRUE, names = FALSE))
H_grid <- data.frame(H_quantile = CONFIG$H_QUANTILES, H = H_vals)
H_main <- as.numeric(quantile(z_ref, probs = CONFIG$H_MAIN_Q, na.rm = TRUE, names = FALSE))
write_csv(H_grid, file.path(DIRS$tables, "kidney_H_grid.csv"))

main <- estimate_methods(kid, H_main, seed = CONFIG$SEED, compute_diag = TRUE)
main_est <- main$estimates
main_est$H_quantile <- CONFIG$H_MAIN_Q
main_est <- main_est[order(main_est$method), , drop = FALSE]
write_csv(main_est, file.path(DIRS$tables, "kidney_main_point_estimates.csv"))
write_csv(main$nuisance, file.path(DIRS$tables, "kidney_nuisance_library_diagnostics.csv"))
write_csv(main$weights, file.path(DIRS$tables, "kidney_weight_diagnostics_main_H.csv"))
write_csv(main$mhat, file.path(DIRS$tables, "kidney_subject_level_nuisance_predictions.csv"))

# H and lambda sensitivity for a compact method set.
# Important: reuse the nuisance predictions from the main analysis.
# This makes the sensitivity analysis vary only H or lambda.
fixed_nuisance <- main$mhat

H_sens <- lapply(seq_along(H_vals), function(j) {
  est <- estimate_methods_fixed_nuisance(kid, H_vals[j], fixed_nuisance)$estimates
  est$H_quantile <- CONFIG$H_QUANTILES[j]
  est
})
H_sens <- do.call(rbind, H_sens)

keep_methods <- c(
  "trad", "ipcw_marg", "ipcw_x",
  "ps_marg_logit", "ps_x_logit", "ps_x_spline",
  "scap_x_l050", "scap_x_l075", "scap_x_l100"
)

H_sens <- H_sens[H_sens$method %in% keep_methods, , drop = FALSE]
H_sens <- H_sens[order(H_sens$H_quantile, H_sens$method), , drop = FALSE]
write_csv(H_sens, file.path(DIRS$tables, "kidney_sensitivity_H.csv"))

lambda_grid_use <- sort(unique(c(CONFIG$LAMBDA_GRID, CONFIG$MAIN_LAMBDAS)))
lambda_sens <- estimate_lambda_grid_fixed_nuisance(
  kid,
  H_main,
  fixed_nuisance,
  lambdas = lambda_grid_use
)
lambda_sens$H_quantile <- CONFIG$H_MAIN_Q
lambda_sens <- lambda_sens[order(lambda_sens$lambda), , drop = FALSE]
write_csv(lambda_sens, file.path(DIRS$tables, "kidney_sensitivity_lambda.csv"))

# Bootstrap
boot_methods <- c("trad", "ipcw_marg", "ipcw_x", "ps_marg_logit", "ps_x_logit", "ps_x_spline", "scap_x_l050", "scap_x_l075")
boot_raw <- run_bootstrap(kid, H_main, boot_methods, CONFIG$N_BOOT, CONFIG$MAX_CORES, CONFIG$SEED + 3000)
write_csv(boot_raw, file.path(DIRS$bootstrap, "kidney_bootstrap_raw.csv"))
boot_summary <- summarize_bootstrap(main_est[main_est$method %in% boot_methods, , drop = FALSE], boot_raw)
write_csv(boot_summary, file.path(DIRS$bootstrap, "kidney_bootstrap_summary.csv"))

# Merge point and bootstrap for plotting/table
main_summary <- merge(main_est, boot_summary, by = "method", all.x = TRUE, suffixes = c("", "_boot"))
write_csv(main_summary, file.path(DIRS$tables, "kidney_main_estimates_with_bootstrap.csv"))

# LaTeX compact tables
tex_est <- main_summary[main_summary$method %in% boot_methods, c("method", "tau_H", "M_hat", "ess", "boot_se", "ci_perc_l", "ci_perc_u", "n_boot_valid"), drop = FALSE]
latex_table(tex_est, file.path(DIRS$latex, "table_kidney_estimates.tex"), digits = 3)
latex_table(summary_tab, file.path(DIRS$latex, "table_kidney_descriptives.tex"), digits = 3)

# Plots
plot_outputs(kid, main_summary, H_sens, lambda_sens, main$mhat)

# Manifest and session info
manifest <- data.frame(path = list.files(OUT, recursive = TRUE, full.names = FALSE))
write_csv(manifest, file.path(OUT, "manifest.csv"))
sink(file.path(OUT, "sessionInfo.txt")); print(sessionInfo()); sink()
saveRDS(CONFIG, file.path(OUT, "config.rds"))

log_msg("Completed kidney real-data application")
log_msg("Key output table:", file.path(DIRS$tables, "kidney_main_estimates_with_bootstrap.csv"))
log_msg("Key bootstrap table:", file.path(DIRS$bootstrap, "kidney_bootstrap_summary.csv"))
