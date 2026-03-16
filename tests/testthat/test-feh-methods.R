# =============================================================================
# reach.hydro — FEH methods test suite
# Description: Unit tests for the FEH L-moments engine, Vol 3 statistical
#              methods, Vol 4 rainfall frequency, ReFH2, and FSR UH.
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# =============================================================================

library(testthat)
library(reach.hydro)

# ---- shared fixtures --------------------------------------------------------

set.seed(99)
# Simulate a 40-year AMAX series from a known GLO(xi=80, alpha=25, k=-0.15)
.amax <- rgev_sim(40, xi = 80, alpha = 25, k = -0.15)  # uses GEV sim; close enough
.amax <- pmax(.amax, 10)  # ensure positive

# Design storm fixture
.storm <- feh_design_storm(
  duration_hr   = 4,
  return_period = 100,
  rmed_1h       = 12,
  rmed_1d       = 38,
  saar          = 900,
  dt_min        = 15
)

# =============================================================================
# L-MOMENTS ENGINE
# =============================================================================

test_that("sample_lmom returns correct structure", {
  lm <- sample_lmom(c(1, 2, 3, 4, 5, 6, 7, 8))
  expect_named(lm, c("l1", "l2", "l3", "l4", "t2", "t3", "t4", "n"))
  expect_equal(lm$l1, mean(c(1:8)))
  expect_gt(lm$l2, 0)
})

test_that("sample_lmom l2 > 0 for any non-constant series", {
  lm <- sample_lmom(c(10, 20, 15, 30, 25, 18))
  expect_gt(lm$l2, 0)
})

test_that("sample_lmom errors with < 4 values", {
  expect_error(sample_lmom(c(1, 2, 3)), "at least 4")
})

test_that("qglo is inverse of GLO CDF at known points", {
  params <- list(xi = 100, alpha = 30, k = -0.1)
  p <- c(0.5, 0.9, 0.99)
  q <- qglo(p, params)
  expect_equal(qglo(0.5, params), params$xi, tolerance = 1e-6)  # median = xi when k=0; ~xi for small k
  expect_true(all(diff(q) > 0))  # monotone
})

test_that("qgev is monotone in p", {
  params <- list(xi = 50, alpha = 20, k = -0.1)
  p <- seq(0.01, 0.99, by = 0.05)
  q <- qgev(p, params)
  expect_true(all(diff(q) > 0))
})

test_that("lmrd_theoretical returns data.table with expected distributions", {
  dt <- lmrd_theoretical()
  expect_true(data.table::is.data.table(dt))
  expect_true(all(c("GLO", "GEV") %in% dt$dist))
})

# =============================================================================
# FEH VOL 3 — SINGLE-SITE
# =============================================================================

test_that("feh_single_site returns FehStatFit for GLO and GEV", {
  for (d in c("glo", "gev")) {
    fit <- feh_single_site(.amax, dist = d, n_boot = 0L)
    expect_s3_class(fit, "FehStatFit")
    expect_equal(fit$dist, d)
    expect_equal(fit$n_years, length(.amax))
  }
})

test_that("feh_single_site growth curve has correct dimensions", {
  rp  <- c(2, 10, 100, 200)
  fit <- feh_single_site(.amax, return_periods = rp, n_boot = 0L)
  expect_equal(nrow(fit$growth_curve), length(rp))
  expect_true(all(c("return_period_yr", "growth_factor", "flow") %in%
                    names(fit$growth_curve)))
})

test_that("feh_single_site growth factors are monotone in return period", {
  fit <- feh_single_site(.amax, n_boot = 0L)
  expect_true(all(diff(fit$growth_curve$growth_factor) > 0))
})

test_that("feh_single_site growth factor = 1 at return period ~2 yr", {
  fit <- feh_single_site(.amax, return_periods = c(2, 10, 100), n_boot = 0L)
  # T=2 growth factor should be close to 1 (median return period)
  expect_lt(abs(fit$growth_curve$growth_factor[1] - 1), 0.15)
})

test_that("feh_single_site bootstrap adds CI columns", {
  fit <- feh_single_site(.amax[1:20], n_boot = 50L, ci_level = 0.90)
  expect_true(all(c("ci_lo", "ci_hi") %in% names(fit$growth_curve)))
  expect_true(all(fit$growth_curve$ci_lo <= fit$growth_curve$flow + 1e-9))
  expect_true(all(fit$growth_curve$ci_hi >= fit$growth_curve$flow - 1e-9))
})

test_that("feh_single_site warns for short record", {
  expect_warning(feh_single_site(.amax[1:8], n_boot = 0L), "< 10 years")
})

test_that("print.FehStatFit runs without error", {
  fit <- feh_single_site(.amax, n_boot = 0L)
  expect_output(print(fit), "FehStatFit")
})

# =============================================================================
# FEH VOL 3 — POOLED
# =============================================================================

test_that("feh_pooled returns FehPooledFit with correct structure", {
  set.seed(7)
  donors <- list(
    D1 = rgev_sim(30, xi = 70,  alpha = 20, k = -0.1),
    D2 = rgev_sim(25, xi = 90,  alpha = 28, k = -0.12),
    D3 = rgev_sim(20, xi = 110, alpha = 35, k = -0.08)
  )
  fit <- feh_pooled(subject_amax = .amax, donor_list = donors,
                    n_boot = 0L)
  expect_s3_class(fit, "FehPooledFit")
  expect_equal(fit$n_donors, 3L)
  expect_gt(fit$pooled_years, 0)
  expect_true(data.table::is.data.table(fit$growth_curve))
})

test_that("feh_pooled errors without subject_amax and subject_qmed", {
  expect_error(feh_pooled(donor_list = list(D1 = .amax)), "subject_amax or subject_qmed")
})

test_that("feh_pooled accepts direct subject_qmed", {
  donors <- list(D1 = .amax)
  fit <- feh_pooled(donor_list = donors, subject_qmed = 95, n_boot = 0L)
  expect_equal(fit$subject_qmed, 95)
})

# =============================================================================
# FEH VOL 3 — POT
# =============================================================================

test_that("feh_pot returns FehPotFit", {
  set.seed(5)
  peaks <- sort(rexp(60, rate = 1/50)) + 30
  fit   <- feh_pot(peaks, threshold = 30, n_years = 20, n_boot = 0L)
  expect_s3_class(fit, "FehPotFit")
  expect_gt(fit$lambda, 0)
})

test_that("feh_pot return period flows are monotone", {
  set.seed(5)
  peaks <- sort(rexp(60, rate = 1/50)) + 30
  fit   <- feh_pot(peaks, threshold = 30, n_years = 20,
                   return_periods = c(2, 10, 50, 100), n_boot = 0L)
  expect_true(all(diff(fit$growth_curve$flow) > 0))
})

test_that("feh_pot errors with too few exceedances", {
  expect_error(feh_pot(c(1, 2, 3), threshold = 0, n_years = 5), "fewer than 5")
})

# =============================================================================
# FEH VOL 4 — RAINFALL FREQUENCY (DDF)
# =============================================================================

test_that("feh_ddf returns data.table with expected columns", {
  dt <- feh_ddf(c(1, 6, 24), return_period = 100,
                rmed_1h = 12, rmed_1d = 38, saar = 900)
  expect_true(data.table::is.data.table(dt))
  expect_true(all(c("duration_hr", "return_period", "rainfall_mm",
                    "growth_factor") %in% names(dt)))
})

test_that("feh_ddf rainfall increases with return period", {
  r10  <- feh_ddf(6, 10,  rmed_1h = 12, rmed_1d = 38, saar = 900)$rainfall_mm
  r100 <- feh_ddf(6, 100, rmed_1h = 12, rmed_1d = 38, saar = 900)$rainfall_mm
  expect_gt(r100, r10)
})

test_that("feh_ddf rainfall increases with duration (for same T)", {
  dt <- feh_ddf(c(1, 6, 24), return_period = 100,
                rmed_1h = 12, rmed_1d = 38, saar = 900)
  expect_true(all(diff(dt$rainfall_mm) > 0))
})

test_that("feh_design_storm total depth matches feh_ddf", {
  storm   <- feh_design_storm(4, 100, rmed_1h = 12, rmed_1d = 38, saar = 900,
                               dt_min = 15)
  ddf_tot <- feh_ddf(4, 100, rmed_1h = 12, rmed_1d = 38, saar = 900)$rainfall_mm
  expect_equal(sum(storm$rainfall_mm), ddf_tot, tolerance = 0.01)
})

test_that("feh_design_storm timesteps are consistent", {
  storm <- feh_design_storm(2, 50, rmed_1h = 10, rmed_1d = 30, saar = 750,
                             dt_min = 15)
  dt    <- diff(storm$time_min)
  expect_true(all(abs(dt - 15) < 1e-9))
})

test_that("feh_arf is in (0, 1] and decreases with area", {
  a1 <- feh_arf(10,  6, 100)
  a2 <- feh_arf(100, 6, 100)
  a3 <- feh_arf(500, 6, 100)
  expect_lte(a1, 1)
  expect_gt(a1, a2)
  expect_gt(a2, a3)
})

# =============================================================================
# ReFH2
# =============================================================================

test_that("refh2_params returns ReFH2Params with positive values", {
  p <- refh2_params(area = 250, bfihost = 0.45, saar = 850)
  expect_s3_class(p, "ReFH2Params")
  expect_gt(p$Cmax, 0)
  expect_gt(p$Tp,   0)
  expect_gt(p$BL,   0)
  expect_true(p$Cini < p$Cmax)
})

test_that("refh2_params Cini is smaller in winter", {
  p_sum <- refh2_params(250, 0.45, 850, season = "summer")
  p_win <- refh2_params(250, 0.45, 850, season = "winter")
  expect_lt(p_win$Cini, p_sum$Cini)
})

test_that("refh2_run returns data.table with correct columns", {
  p   <- refh2_params(250, 0.45, 850)
  res <- refh2_run(.storm, p)
  expect_true(data.table::is.data.table(res))
  expect_true(all(c("Q_mm", "Qd_mm", "Qb_mm", "net_rain_mm") %in% names(res)))
})

test_that("refh2_run Q >= 0 at all timesteps", {
  p   <- refh2_params(250, 0.45, 850)
  res <- refh2_run(.storm, p)
  expect_true(all(res$Q_mm >= -1e-9))
})

test_that("refh2_run peak_flow_mm attribute is set", {
  p   <- refh2_params(250, 0.45, 850)
  res <- refh2_run(.storm, p)
  expect_gt(attr(res, "peak_flow_mm"), 0)
})

test_that("refh2_summary returns list with expected names", {
  p   <- refh2_params(250, 0.45, 850)
  res <- refh2_run(.storm, p)
  s   <- refh2_summary(res, area_km2 = 250)
  expect_named(s, c("peak_mm", "volume_mm_hr", "time_to_peak_hr",
                    "peak_m3s", "volume_Mm3"))
  expect_gt(s$peak_m3s, 0)
})

test_that("print.ReFH2Params runs without error", {
  p <- refh2_params(250, 0.45, 850)
  expect_output(print(p), "ReFH2Params")
})

# =============================================================================
# FSR/FEH UNIT HYDROGRAPH
# =============================================================================

test_that("fsr_params returns FsrParams with plausible SPR", {
  p <- fsr_params(area = 150, bfihost = 0.4, saar = 800)
  expect_s3_class(p, "FsrParams")
  expect_true(p$spr > 0 && p$spr < 100)
  expect_gt(p$Tp, 0)
})

test_that("fsr_unit_hydrograph sums to 1", {
  uh <- fsr_unit_hydrograph(Tp = 3, dt_hr = 0.25)
  expect_equal(sum(uh), 1, tolerance = 1e-9)
})

test_that("fsr_unit_hydrograph peak is at Tp", {
  dt   <- 0.25
  Tp   <- 4
  uh   <- fsr_unit_hydrograph(Tp = Tp, dt_hr = dt)
  peak_t <- (which.max(uh) - 1) * dt
  expect_equal(peak_t, Tp, tolerance = dt)
})

test_that("fsr_percentage_runoff is in [0, 100]", {
  pr <- fsr_percentage_runoff(spr = 40, cwi = 130, storm_depth = 30, m5_60min = 12)
  expect_gte(pr, 0)
  expect_lte(pr, 100)
})

test_that("fsr_percentage_runoff increases with storm depth", {
  pr_low  <- fsr_percentage_runoff(40, 125, 10, 12)
  pr_high <- fsr_percentage_runoff(40, 125, 50, 12)
  expect_gt(pr_high, pr_low)
})

test_that("fsr_run returns data.table with Q_mm >= 0", {
  p   <- fsr_params(150, 0.4, 800)
  res <- fsr_run(.storm, p, m5_60min = 12)
  expect_true(data.table::is.data.table(res))
  expect_true(all(res$Q_mm >= -1e-9))
})

test_that("fsr_run PR attribute is between 0 and 100", {
  p   <- fsr_params(150, 0.4, 800)
  res <- fsr_run(.storm, p, m5_60min = 12)
  expect_gte(attr(res, "PR"), 0)
  expect_lte(attr(res, "PR"), 100)
})

test_that("fsr_run pr_override is respected", {
  p   <- fsr_params(150, 0.4, 800)
  res <- fsr_run(.storm, p, m5_60min = 12, pr_override = 55)
  expect_equal(attr(res, "PR"), 55)
})

# =============================================================================
# BACKWARD-COMPATIBLE WRAPPERS
# =============================================================================

test_that("fit_gev and fit_glo return FehStatFit", {
  expect_s3_class(fit_gev(.amax), "FehStatFit")
  expect_s3_class(fit_glo(.amax), "FehStatFit")
})

test_that("return_period_flow returns correct number of rows", {
  fit <- fit_glo(.amax)
  rp  <- c(10, 50, 100)
  out <- return_period_flow(fit, return_periods = rp)
  expect_equal(nrow(out), length(rp))
  expect_true(all(diff(out$flow) > 0))
})
