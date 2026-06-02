need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not available: ", pkg)
  }
}

need_pkg("survival")
library(survival)

data(rhDNase, package = "survival")

cat("Raw rhDNase data:\n")
print(dim(rhDNase))
print(names(rhDNase))
cat("\nNumber of unique subjects:\n")
print(length(unique(rhDNase$id)))

first <- rhDNase[!duplicated(rhDNase$id), ]

dnase <- tmerge(
  first,
  first,
  id = id,
  tstop = as.numeric(end.dt - entry.dt)
)

temp_end <- with(
  rhDNase,
  pmin(ivstop + 6, as.numeric(end.dt - entry.dt))
)

dnase <- tmerge(
  dnase,
  rhDNase,
  id = id,
  infect = event(ivstart),
  end = event(temp_end)
)

dnase <- subset(
  dnase,
  infect == 1 | end == 0,
  select = c(id, inst, trt, fev, tstart, tstop, infect)
)

dnase <- dnase[order(dnase$id, dnase$tstart, dnase$tstop), ]
dnase$duration <- pmax(0, dnase$tstop - dnase$tstart)

build_one_subject <- function(d) {
  d <- d[order(d$tstart, d$tstop), ]
  d$at_risk_stop <- cumsum(d$duration)

  event_times <- d$at_risk_stop[d$infect == 1]
  total_at_risk <- sum(d$duration)

  Delta1 <- as.integer(length(event_times) >= 1)
  Delta2 <- as.integer(length(event_times) >= 2)

  A_obs <- if (Delta1 == 1) event_times[1] else total_at_risk
  Z_obs <- if (Delta2 == 1) event_times[2] else total_at_risk
  B_obs <- if (Delta1 == 1) Z_obs - A_obs else NA_real_

  data.frame(
    id = d$id[1],
    inst = d$inst[1],
    trt = d$trt[1],
    fev = d$fev[1],
    C_atrisk = total_at_risk,
    A_tilde = A_obs,
    Z_tilde = Z_obs,
    B_tilde = B_obs,
    Delta1 = Delta1,
    Delta2 = Delta2,
    n_events_atrisk = length(event_times)
  )
}

gap <- do.call(
  rbind,
  lapply(split(dnase, dnase$id), build_one_subject)
)

gap <- gap[is.finite(gap$C_atrisk) & gap$C_atrisk > 0, ]

cat("\nSubject-level at-risk gap data:\n")
print(dim(gap))

cat("\nEvent observability summary:\n")
print(data.frame(
  subjects = nrow(gap),
  Delta1_observed = sum(gap$Delta1 == 1),
  Delta2_observed = sum(gap$Delta2 == 1),
  first_gap_censoring_rate = mean(gap$Delta1 == 0),
  second_gap_censoring_rate_all = mean(gap$Delta2 == 0),
  second_gap_censoring_rate_given_Delta1 = mean(gap$Delta2[gap$Delta1 == 1] == 0)
))

cat("\nDistribution of number of at-risk exacerbations:\n")
print(table(gap$n_events_atrisk))

cat("\nTreatment table:\n")
print(table(gap$trt, useNA = "ifany"))

cat("\nFEV summary:\n")
print(summary(gap$fev))

dir.create("results_realdata_rhDNase_tauH/check_gapdata", recursive = TRUE, showWarnings = FALSE)
write.csv(
  gap,
  "results_realdata_rhDNase_tauH/check_gapdata/rhDNase_gap_data_subject_level.csv",
  row.names = FALSE
)
write.csv(
  dnase,
  "results_realdata_rhDNase_tauH/check_gapdata/rhDNase_atrisk_start_stop.csv",
  row.names = FALSE
)

cat("\nSaved:\n")
cat("results_realdata_rhDNase_tauH/check_gapdata/rhDNase_gap_data_subject_level.csv\n")
cat("results_realdata_rhDNase_tauH/check_gapdata/rhDNase_atrisk_start_stop.csv\n")
