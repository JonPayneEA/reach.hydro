# =============================================================================
# Tool:        reach.hydro test suite
# Description: Unit tests for all exported functions. Uses testthat 3e.
#              Tier 1 functions have full edge-case and regression coverage;
#              Tier 2 functions have basic smoke tests.
# Author:      Forecasting and Warning Team
# Created:     2026-02-01
# =============================================================================

library(testthat)
library(reach.hydro)

# ---- shared fixtures --------------------------------------------------------

set.seed(42)
.n    <- 365L * 2L
.tsq  <- seq_len(.n)
.seas <- 0.5 + 0.5 * cos(2 * pi * (.tsq - 30) / 365)
.rain <- pmax(0, stats::rnorm(.n, mean = 2.5 * .seas, sd = 3))
.pet  <- pmax(0, 3 * (1 - .seas) + stats::rnorm(.n, 0, 0.3))

.p_pareto <- pdm_params(dist = "pareto", cmax = 350, b = 0.4, St = 20,
                        kb = 150, m = 1, k1 = 5, k2 = 8)

# =============================================================================
# PdmParams construction and validation
# =============================================================================

test_that("pdm_params creates valid PdmParams for all distributions", {
  dists <- c("pareto", "rectangular", "exponential", "triangular", "lognormal")
  for (d in dists) {
    p <- pdm_params(dist = d)
    expect_s3_class(p, "PdmParams")
    expect_equal(p$dist, d)
  }
})

test_that("pdm_validate_params catches invalid parameters", {
  expect_error(pdm_params(cmax = -10),   "cmax must be > 0")
  expect_error(pdm_params(cmin = -1),    "cmin must be >= 0")
  expect_error(pdm_params(cmin = 400, cmax = 300), "cmin must be < cmax")
  expect_error(pdm_params(b = 0),        "b .shape. must be > 0")
  expect_error(pdm_params(sigma_c = 0),  "sigma_c must be > 0")
  expect_error(pdm_params(alpha = 1.5),  "alpha must be in")
})

test_that("print.PdmParams runs without error", {
  expect_output(print(.p_pareto), "PdmParams")
})

# =============================================================================
# Capacity distributions — Smax, CDF, c*, runoff
# =============================================================================

test_that("capacity_smax returns positive value for all distributions", {
  params <- list(cmin=0, cmax=300, b=0.5, mu_c=150, sigma_c=60,
                 mu_lnc=4.9, sigma_lnc=0.4)
  dists  <- c("pareto","rectangular","exponential","triangular","lognormal")
  for (d in dists) {
    s <- capacity_smax(d, params)
    expect_gt(s, 0, label = paste("Smax for", d))
  }
})

test_that("Pareto Smax matches Moore (2007) Appendix A closed form", {
  # Smax = cmin + (cmax-cmin)*b/(1+b)
  p <- list(cmin=10, cmax=400, b=0.5)
  expect_equal(capacity_smax("pareto", p), 10 + 390 * 0.5 / 1.5, tolerance = 1e-10)
})

test_that("Uniform Smax matches closed form", {
  p <- list(cmin=20, cmax=200)
  expect_equal(capacity_smax("rectangular", p), (20 + 200) / 2, tolerance = 1e-10)
})

test_that("Exponential Smax equals cmax", {
  p <- list(cmax=250)
  expect_equal(capacity_smax("exponential", p), 250)
})

test_that("capacity_cdf returns values in [0, 1]", {
  c_vals <- seq(0, 400, by = 50)
  params <- list(cmin=0, cmax=400, b=0.5, mu_c=200, sigma_c=80,
                 mu_lnc=5.0, sigma_lnc=0.5)
  for (d in c("pareto","rectangular","exponential","triangular","lognormal")) {
    f <- capacity_cdf(c_vals, d, params)
    expect_true(all(f >= 0 & f <= 1), label = paste("CDF bounds:", d))
    expect_true(all(diff(f) >= 0),    label = paste("CDF monotone:", d))
  }
})

test_that("capacity_runoff is bounded to [0, P] for all distributions", {
  params <- list(cmin=0, cmax=350, b=0.4, mu_c=175, sigma_c=70,
                 mu_lnc=5.0, sigma_lnc=0.5)
  for (d in c("pareto","rectangular","exponential","triangular","lognormal")) {
    for (S in c(0, 50, 150, 300)) {
      Qr <- capacity_runoff(P = 20, S = S, dist = d, params = params)
      expect_gte(Qr, 0,  label = paste("Runoff >= 0:", d, "S=", S))
      expect_lte(Qr, 20, label = paste("Runoff <= P:", d, "S=", S))
    }
  }
})

test_that("capacity_runoff returns 0 for P = 0", {
  params <- list(cmin=0, cmax=300, b=0.5)
  expect_equal(capacity_runoff(0, 100, "pareto", params), 0)
})

# =============================================================================
# PDM core
# =============================================================================

test_that("pdm() returns a ReachHydroResult with correct dimensions", {
  res <- pdm(.rain, .pet, params = .p_pareto)
  expect_s3_class(res, "ReachHydroResult")
  expect_equal(nrow(res), .n)
  expect_true(all(c("Q", "Qf", "Qb", "S", "AET", "SMD") %in% names(res)))
})

test_that("pdm() total flow = Qf + Qb at every timestep", {
  res <- pdm(.rain, .pet, params = .p_pareto)
  expect_equal(res$Q, res$Qf + res$Qb, tolerance = 1e-10)
})

test_that("pdm() soil moisture stays within [0, Smax]", {
  res  <- pdm(.rain, .pet, params = .p_pareto)
  Smax <- attr(res, "Smax")
  expect_true(all(res$S >= -1e-9))
  expect_true(all(res$S <= Smax + 1e-9))
})

test_that("pdm() SMD = Smax - S", {
  res  <- pdm(.rain, .pet, params = .p_pareto)
  Smax <- attr(res, "Smax")
  expect_equal(res$SMD, Smax - res$S, tolerance = 1e-10)
})

test_that("pdm() stores are non-negative at every timestep", {
  res <- pdm(.rain, .pet, params = .p_pareto)
  expect_true(all(res$Sg >= -1e-9))
  expect_true(all(res$Ss >= -1e-9))
  expect_true(all(res$Q  >= -1e-9))
})

test_that("pdm() accepts plain list params (backward-compatible)", {
  res <- pdm(.rain, .pet,
             params = list(cmax = 300, b = 0.4, St = 10, kb = 100, m = 1, k1 = 5, k2 = 5),
             dist = "pareto")
  expect_s3_class(res, "ReachHydroResult")
})

test_that("pdm() errors on mismatched rain/pet lengths", {
  expect_error(pdm(.rain, .pet[1:10], params = .p_pareto), "same length")
})

test_that("pdm() errors on empty input", {
  expect_error(pdm(numeric(0), numeric(0), params = .p_pareto), "empty")
})

test_that("pdm() runs successfully for all five distributions", {
  dists <- c("pareto","rectangular","exponential","triangular","lognormal")
  for (d in dists) {
    p <- pdm_params(dist = d, cmax = 300, mu_lnc = 4.9, sigma_lnc = 0.4, k1 = 5, k2 = 5)
    expect_no_error(pdm(.rain[1:100], .pet[1:100], params = p),
                    label = paste("pdm() with dist =", d))
  }
})

test_that("print and summary methods run without error", {
  res <- pdm(.rain, .pet, params = .p_pareto)
  expect_output(print(res),   "ReachHydroResult")
  expect_output(summary(res), "distribution|Smax")
})

# =============================================================================
# Goodness-of-fit metrics
# =============================================================================

test_that("nse returns 1 for perfect simulation", {
  obs <- c(1, 2, 3, 4, 5)
  expect_equal(nse(obs, obs), 1)
})

test_that("nse returns 0 when sim = mean(obs)", {
  obs <- c(1, 2, 3, 4, 5)
  sim <- rep(mean(obs), 5)
  expect_equal(nse(obs, sim), 0, tolerance = 1e-10)
})

test_that("kge returns 1 for perfect simulation", {
  obs <- c(1, 2, 3, 4, 5)
  expect_equal(kge(obs, obs), 1, tolerance = 1e-10)
})

test_that("pbias returns 0 for perfect simulation", {
  obs <- c(1, 2, 3)
  expect_equal(pbias(obs, obs), 0, tolerance = 1e-10)
})

test_that("pbias sign: over-prediction is positive", {
  obs <- c(1, 1, 1)
  sim <- c(2, 2, 2)
  expect_gt(pbias(obs, sim), 0)
})

test_that("far returns 0 when no false alarms", {
  obs <- c(TRUE,  TRUE,  FALSE)
  sim <- c(TRUE,  TRUE,  FALSE)
  expect_equal(far(obs, sim), 0)
})

test_that("far returns 1 when all alarms are false", {
  obs <- c(FALSE, FALSE, FALSE)
  sim <- c(TRUE,  TRUE,  TRUE)
  expect_equal(far(obs, sim), 1)
})

test_that("far returns NA when no simulated exceedances", {
  expect_true(is.na(far(c(TRUE, FALSE), c(FALSE, FALSE))))
})

test_that("metrics handle NA pairs gracefully", {
  obs <- c(1, NA, 3)
  sim <- c(1,  2, 3)
  expect_no_error(nse(obs, sim))
  expect_no_error(kge(obs, sim))
  expect_no_error(pbias(obs, sim))
})

test_that("gof_metrics returns named vector", {
  obs <- .rain[1:100]
  sim <- .rain[1:100] + stats::rnorm(100, 0, 0.1)
  m   <- gof_metrics(obs, sim)
  expect_named(m, c("nse", "kge", "pbias"))
})

# =============================================================================
# Flow statistics
# =============================================================================

test_that("flow_stats returns one-row data.table", {
  res <- flow_stats(c(1, 2, 3, 4, 5))
  expect_equal(nrow(res), 1L)
  expect_true(data.table::is.data.table(res))
})

test_that("flow_stats handles NA values", {
  res <- flow_stats(c(1, NA, 3))
  expect_equal(res$n_valid, 2L)
})

test_that("baseflow_index is in [0, 1]", {
  res <- pdm(.rain, .pet, params = .p_pareto)
  bfi <- baseflow_index(res$Qb, res$Q)
  expect_gte(bfi, 0)
  expect_lte(bfi, 1)
})

test_that("flow_duration_curve has monotone decreasing flow", {
  fdc <- flow_duration_curve(.rain, n_bins = 50L)
  expect_true(all(diff(fdc$flow) <= 0 + 1e-9))
})

test_that("annual_maxima extracts one value per water year", {
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 730)
  am    <- annual_maxima(.rain[1:730], dates)
  expect_equal(nrow(am), length(unique(am$water_year)))
})

test_that("peaks_over_threshold returns empty dt when no exceedances", {
  res <- peaks_over_threshold(rep(0, 100), seq.Date(as.Date("2020-01-01"),
                              by = "day", length.out = 100), threshold = 10)
  expect_equal(nrow(res), 0L)
})

# =============================================================================
# baseflow_separate
# =============================================================================

test_that("baseflow_separate (lyne_hollick): bfi is in (0, 1)", {
  Q  <- pmax(0.1, .rain[1:365])
  bf <- baseflow_separate(Q, method = "lyne_hollick")
  expect_gt(bf$bfi, 0)
  expect_lt(bf$bfi, 1)
})

test_that("baseflow_separate (boughton_eckhardt): bfi is in (0, 1)", {
  Q  <- pmax(0.1, .rain[1:365])
  bf <- baseflow_separate(Q, method = "boughton_eckhardt")
  expect_gt(bf$bfi, 0)
  expect_lt(bf$bfi, 1)
})

test_that("baseflow_separate: baseflow + quickflow == flow", {
  Q  <- pmax(0.1, .rain[1:200])
  bf <- baseflow_separate(Q)
  expect_equal(bf$baseflow + bf$quickflow, Q, tolerance = 1e-10)
})

test_that("baseflow_separate: baseflow <= flow at every timestep", {
  Q  <- pmax(0.1, .rain[1:200])
  bf <- baseflow_separate(Q)
  expect_true(all(bf$baseflow <= Q + 1e-10))
})

test_that("baseflow_separate: quickflow is non-negative", {
  Q  <- pmax(0.1, .rain[1:200])
  bf <- baseflow_separate(Q)
  expect_true(all(bf$quickflow >= -1e-10))
})

# =============================================================================
# q_n_day
# =============================================================================

test_that("q_n_day (min) is <= overall minimum", {
  Q <- pmax(0, .rain[1:200])
  expect_lte(q_n_day(Q, n = 7L, type = "min"), min(Q) + 1e-10)
})

test_that("q_n_day (max) is >= overall maximum", {
  Q <- pmax(0, .rain[1:200])
  expect_gte(q_n_day(Q, n = 7L, type = "max"), max(Q) - 1e-10)
})

test_that("q_n_day with dates returns one row per water year", {
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 730)
  Q     <- pmax(0, .rain[1:730])
  res   <- q_n_day(Q, n = 7L, dates = dates)
  expect_true(data.table::is.data.table(res))
  expect_equal(nrow(res), length(unique(res$water_year)))
})

# =============================================================================
# monthly_flow_stats
# =============================================================================

test_that("monthly_flow_stats returns 12 rows", {
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 730)
  Q     <- pmax(0, .rain[1:730])
  res   <- monthly_flow_stats(Q, dates)
  expect_equal(nrow(res), 12L)
  expect_true(data.table::is.data.table(res))
})

test_that("monthly_flow_stats: Q10 <= median <= Q90 for every month", {
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 730)
  Q     <- pmax(0.1, .rain[1:730])
  res   <- monthly_flow_stats(Q, dates)
  expect_true(all(res$Q10 <= res$median + 1e-9))
  expect_true(all(res$median <= res$Q90 + 1e-9))
})

# =============================================================================
# flow_deficit
# =============================================================================

test_that("flow_deficit returns zero rows when no values below threshold", {
  Q   <- rep(5, 50)
  res <- flow_deficit(Q, threshold = 1)
  expect_equal(nrow(res), 0L)
})

test_that("flow_deficit identifies at least one spell", {
  Q   <- c(5, 5, 1, 0.5, 0.8, 5, 5, 0.2, 0.3, 5)
  res <- flow_deficit(Q, threshold = 3)
  expect_gte(nrow(res), 1L)
})

test_that("flow_deficit: deficit_volume > 0 for all spells", {
  Q   <- c(5, 5, 1, 0.5, 0.8, 5, 5, 0.2, 0.3, 5)
  res <- flow_deficit(Q, threshold = 3)
  expect_true(all(res$deficit_volume > 0))
})

test_that("flow_deficit: start_idx <= end_idx for all spells", {
  Q   <- c(5, 5, 1, 0.5, 0.8, 5, 5, 0.2, 0.3, 5)
  res <- flow_deficit(Q, threshold = 3)
  expect_true(all(res$start_idx <= res$end_idx))
})

test_that("flow_deficit adds date columns when dates supplied", {
  Q     <- c(5, 5, 1, 0.5, 5, 5)
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 6)
  res   <- flow_deficit(Q, threshold = 3, dates = dates)
  expect_true("start_date" %in% names(res))
  expect_true("end_date"   %in% names(res))
})

# =============================================================================
# flow_recession
# =============================================================================

test_that("flow_recession returns positive k for a clean recession", {
  Q   <- 20 * exp(-seq(0, 11) / 8)   # perfect exponential decay, k = 8
  res <- flow_recession(Q, min_duration = 4L)
  expect_gte(res$n_recessions, 1L)
  expect_gt(res$k, 0)
})

test_that("flow_recession: k is NA when no recessions found", {
  Q   <- rep(5, 20)
  res <- flow_recession(Q, min_duration = 5L)
  expect_equal(res$n_recessions, 0L)
  expect_true(is.na(res$k))
})

test_that("flow_recession: recessions data.table has correct columns", {
  Q   <- c(20, 15, 11, 8, 6, 5, 4, 3.5, 3, 2.8)
  res <- flow_recession(Q, min_duration = 3L)
  expect_true(all(c("start_idx", "end_idx", "Q0", "k_fit") %in%
                    names(res$recessions)))
})

# =============================================================================
# Rainfall statistics
# =============================================================================

test_that("rainfall_events returns zero rows when no rain", {
  res <- rainfall_events(rep(0, 50))
  expect_equal(nrow(res), 0L)
})

test_that("rainfall_events: events are non-overlapping", {
  rain <- c(0, 0, 3, 5, 2, 0, 0, 0, 0, 0, 0, 4, 3, 0, 0)
  res  <- rainfall_events(rain, min_dry = 4L)
  if (nrow(res) >= 2L) {
    expect_true(all(res$start_idx[-1L] > res$end_idx[-nrow(res)]))
  }
  expect_gte(nrow(res), 1L)
})

test_that("rainfall_events: total_depth matches sum of rain in event", {
  rain <- c(0, 0, 3, 5, 2, 0, 0, 0, 0, 0, 0, 4, 3, 0, 0)
  res  <- rainfall_events(rain, min_dry = 5L)
  for (j in seq_len(nrow(res))) {
    expect_equal(
      sum(rain[res$start_idx[j]:res$end_idx[j]]),
      res$total_depth[j],
      tolerance = 1e-10
    )
  }
})

test_that("api: returns same-length vector, non-negative", {
  rain <- c(0, 0, 10, 5, 0, 0, 0, 3, 0, 0)
  out  <- api(rain, k = 0.9)
  expect_length(out, length(rain))
  expect_true(all(out >= 0))
})

test_that("api: decays to near-zero in long dry period", {
  rain <- c(10, rep(0, 50))
  out  <- api(rain, k = 0.9)
  expect_lt(tail(out, 1L), 0.1)
})

test_that("api: equals rain when k = 0 (no memory)", {
  rain <- c(2, 5, 0, 3, 1)
  expect_equal(api(rain, k = 0), rain, tolerance = 1e-10)
})

test_that("idf_empirical: returns one row per duration when no dates", {
  rain <- pmax(0, .rain[1:365])
  res  <- idf_empirical(rain, durations = c(1L, 3L, 7L))
  expect_equal(nrow(res), 3L)
  expect_true(data.table::is.data.table(res))
})

test_that("idf_empirical: longer durations have >= depth than shorter", {
  rain <- pmax(0, .rain[1:365])
  res  <- idf_empirical(rain, durations = c(1L, 7L, 30L))
  expect_gte(res$amax_depth[2L], res$amax_depth[1L])
  expect_gte(res$amax_depth[3L], res$amax_depth[2L])
})

test_that("idf_empirical with dates: returns water_year column", {
  rain  <- pmax(0, .rain[1:730])
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 730)
  res   <- idf_empirical(rain, durations = c(1L, 3L), dates = dates)
  expect_true("water_year" %in% names(res))
})

# =============================================================================
# HydroData dispatch -- flow and rainfall analysis functions
# =============================================================================

# Minimal mock builder (mirrors test-reach-io-compat.R approach)
.make_mock_flow_hd <- function(vals, start = "2020-10-01", freq = "day") {
  n     <- length(vals)
  dates <- seq.Date(as.Date(start), by = freq, length.out = n)
  dts   <- as.POSIXct(dates)
  structure(
    list(
      readings = data.table::data.table(
        dateTime = dts, date = dates,
        value    = vals, measure_notation = "m3s",
        quality  = rep("Good", n)
      ),
      parameter = "flow", period_name = "Mock", n_measures = n, n_rows = n,
      from_date = as.character(start),
      to_date   = as.character(max(dates)),
      downloaded_at = Sys.time()
    ),
    class = c("Flow_Daily", "HydroData")
  )
}

.make_mock_rain_hd <- function(vals, start = "2020-10-01", freq = "day") {
  n     <- length(vals)
  dates <- seq.Date(as.Date(start), by = freq, length.out = n)
  dts   <- as.POSIXct(dates)
  structure(
    list(
      readings = data.table::data.table(
        dateTime = dts, date = dates,
        value    = vals, measure_notation = "mm",
        quality  = rep("Good", n)
      ),
      parameter = "rainfall", period_name = "Mock", n_measures = n, n_rows = n,
      from_date = as.character(start),
      to_date   = as.character(max(dates)),
      downloaded_at = Sys.time()
    ),
    class = c("Rainfall_Daily", "HydroData")
  )
}

test_that("baseflow_separate accepts Flow_Daily HydroData", {
  hd <- .make_mock_flow_hd(pmax(0.1, .rain[1:100]))
  bf <- baseflow_separate(hd)
  expect_named(bf, c("baseflow", "quickflow", "bfi"))
  expect_length(bf$baseflow, 100L)
})

test_that("q_n_day accepts Flow_Daily HydroData and returns annual table", {
  hd  <- .make_mock_flow_hd(pmax(0, .rain[1:730]))
  res <- q_n_day(hd, n = 7L)
  expect_true(data.table::is.data.table(res))
  expect_true("water_year" %in% names(res))
})

test_that("monthly_flow_stats accepts Flow_Daily HydroData", {
  hd  <- .make_mock_flow_hd(pmax(0, .rain[1:730]))
  res <- monthly_flow_stats(hd)
  expect_equal(nrow(res), 12L)
})

test_that("flow_deficit accepts Flow_Daily HydroData and adds date columns", {
  hd  <- .make_mock_flow_hd(c(5, 5, 1, 0.5, 5, 5, 5, 1, 0.2, 5))
  res <- flow_deficit(hd, threshold = 3)
  expect_gte(nrow(res), 1L)
  expect_true("start_date" %in% names(res))
})

test_that("flow_recession accepts Flow_Daily HydroData", {
  Q  <- 20 * exp(-seq(0, 9) / 5)
  hd <- .make_mock_flow_hd(Q)
  res <- flow_recession(hd, min_duration = 4L)
  expect_true("start_date" %in% names(res$recessions))
})

test_that("rainfall_events accepts Rainfall_Daily HydroData and adds date columns", {
  hd  <- .make_mock_rain_hd(c(0, 0, 3, 5, 0, 0, 0, 0, 0, 0, 2, 4, 0, 0))
  res <- rainfall_events(hd, min_dry = 4L)
  expect_gte(nrow(res), 1L)
  expect_true("start_date" %in% names(res))
})

test_that("api accepts Rainfall_Daily HydroData", {
  hd  <- .make_mock_rain_hd(c(0, 0, 10, 5, 0, 0, 0, 3, 0, 0))
  out <- api(hd, k = 0.9)
  expect_length(out, 10L)
  expect_true(all(out >= 0))
})

test_that("idf_empirical accepts Rainfall_Daily HydroData and returns water_year", {
  hd  <- .make_mock_rain_hd(pmax(0, .rain[1:730]))
  res <- idf_empirical(hd, durations = c(1L, 3L))
  expect_true("water_year" %in% names(res))
})

test_that("flow_deficit errors when wrong HydroData type passed", {
  rain_hd <- .make_mock_rain_hd(rep(2, 20))
  expect_error(flow_deficit(rain_hd, threshold = 1), "Flow_Daily or Flow_15min")
})

test_that("rainfall_events errors when wrong HydroData type passed", {
  flow_hd <- .make_mock_flow_hd(pmax(0.1, .rain[1:20]))
  expect_error(rainfall_events(flow_hd), "Rainfall_Daily or Rainfall_15min")
})

# =============================================================================
# Unit conversions
# =============================================================================

test_that("mm_to_m3s and m3s_to_mm are inverse operations", {
  mm  <- c(1, 5, 10)
  m3s <- mm_to_m3s(mm, area_km2 = 100, dt_hours = 1)
  mm2 <- m3s_to_mm(m3s, area_km2 = 100, dt_hours = 1)
  expect_equal(mm2, mm, tolerance = 1e-10)
})

# =============================================================================
# compare_distributions
# =============================================================================

test_that("compare_distributions returns 5 rows", {
  tbl <- compare_distributions(.rain[1:100], .pet[1:100])
  expect_equal(nrow(tbl), 5L)
  expect_true(data.table::is.data.table(tbl))
})

test_that("compare_distributions adds NSE/KGE/PBIAS when obs_q supplied", {
  sim_q <- pdm(.rain[1:100], .pet[1:100], params = .p_pareto)$Q
  tbl   <- compare_distributions(.rain[1:100], .pet[1:100], obs_q = sim_q)
  expect_true(all(c("NSE", "KGE", "PBIAS") %in% names(tbl)))
})
