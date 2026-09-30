# =============================================================================
# reach.hydro — ReFH parameter class test suite
# Description: FlodeReFHParams (S7) construction, ReFH_params(), and
#              ReFH_run() parameter handling.
# Author:      Forecasting and Warning Team
# Created:     2026-09-30
# =============================================================================

library(testthat)
library(reach.hydro)

test_that("ReFH_params warns and returns a FlodeReFHParams S7 object", {
  p <- suppressWarnings(
    ReFH_params(area = 250, bfihost = 0.45, saar = 850, farl = 0.98)
  )
  expect_true(S7::S7_inherits(p, FlodeReFHParams))
  expect_true(p@Cmax > 0)
  expect_true(p@Cini <= p@Cmax)
})

test_that("ReFH_params raises the unverified-coefficients warning", {
  expect_warning(
    ReFH_params(area = 250, bfihost = 0.45, saar = 850),
    regexp = "UNVERIFIED"
  )
})

test_that("FlodeReFHParams can be constructed directly with named properties", {
  p <- FlodeReFHParams(
    Cmax = 100, Cini = 40, alpha = 0.4, Tp = 5, BL = 20, BR = 0.02, Kb = 0.98,
    area = 250, bfihost = 0.45, saar = 850, farl = 1, urbext = 0, season = "summer"
  )
  expect_true(S7::S7_inherits(p, FlodeReFHParams))
  expect_equal(p@area, 250)
  expect_equal(p@season, "summer")
})

test_that("printing a FlodeReFHParams object does not error", {
  p <- suppressWarnings(
    ReFH_params(area = 250, bfihost = 0.45, saar = 850, farl = 0.98)
  )
  expect_output(print(p), "FlodeReFHParams")
})

test_that("ReFH_run accepts a FlodeReFHParams object", {
  storm <- feh_design_storm(
    duration_hr   = 4,
    return_period = 100,
    rmed_1h       = 12,
    rmed_1d       = 38,
    saar          = 900,
    dt_min        = 15
  )
  p <- suppressWarnings(
    ReFH_params(area = 250, bfihost = 0.45, saar = 850, farl = 0.98)
  )
  result <- ReFH_run(storm, p)
  expect_true(data.table::is.data.table(result))
  expect_true(all(c("Qd_mm", "Qb_mm", "Q_mm") %in% names(result)))
  expect_true(S7::S7_inherits(attr(result, "params"), FlodeReFHParams))
})

test_that("ReFH_run still accepts a plain named list for params", {
  storm <- feh_design_storm(
    duration_hr   = 4,
    return_period = 100,
    rmed_1h       = 12,
    rmed_1d       = 38,
    saar          = 900,
    dt_min        = 15
  )
  plain_params <- list(
    Cmax = 100, Cini = 40, alpha = 0.4, Tp = 5, BL = 20, BR = 0.02, Kb = 0.98,
    area = 250, bfihost = 0.45, saar = 850, farl = 1, urbext = 0, season = "summer"
  )
  result <- ReFH_run(storm, plain_params)
  expect_true(data.table::is.data.table(result))
  expect_true(S7::S7_inherits(attr(result, "params"), FlodeReFHParams))
})
