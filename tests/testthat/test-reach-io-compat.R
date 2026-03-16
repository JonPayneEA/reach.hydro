# =============================================================================
# reach.hydro — reach.io compatibility layer tests
# Description: Tests for the reach.io compat layer. Two categories:
#              1. Tests that run without reach.io (mock objects, error paths)
#              2. Tests skipped unless reach.io is installed
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# =============================================================================

library(testthat)
library(reach.hydro)

# ---- helpers: minimal mock HydroData objects --------------------------------
# These mimic just enough of the reach.io S7 structure to test the compat
# layer without requiring reach.io to be installed.

.make_mock_hydrodata <- function(class_name, n = 100L, start = "2020-01-01",
                                 freq = "day") {
  dates <- seq.Date(as.Date(start), by = freq, length.out = n)
  dts   <- as.POSIXct(dates)
  vals  <- pmax(0, rnorm(n, mean = 5, sd = 2))
  readings <- data.table::data.table(
    dateTime          = dts,
    date              = dates,
    value             = vals,
    measure_notation  = rep("m", n),
    quality           = sample(c("Good", "Suspect", "Missing"),
                               n, replace = TRUE, prob = c(0.85, 0.1, 0.05))
  )
  structure(
    list(
      readings     = readings,
      parameter    = "flow",
      period_name  = "Test Site",
      from_date    = as.character(start),
      to_date      = as.character(max(dates)),
      n_measures   = n,
      n_rows       = n,
      downloaded_at = Sys.time()
    ),
    class = c(class_name, "HydroData")
  )
}

# Mock reach.io::as_data_table generic for tests
# (only registered if reach.io is NOT available, to avoid conflicts)
if (!requireNamespace("reach.io", quietly = TRUE)) {
  # Register a minimal as_data_table that works with our mocks
  reach.io_as_data_table_mock <- function(x) x$readings
}

# =============================================================================
# CLASS DETECTION
# =============================================================================

test_that(".is_hydrodata detects HydroData subclasses", {
  mock <- .make_mock_hydrodata("Flow_Daily")
  expect_true(reach.hydro:::.is_hydrodata(mock))
})

test_that(".is_hydrodata returns FALSE for plain objects", {
  expect_false(reach.hydro:::.is_hydrodata(c(1, 2, 3)))
  expect_false(reach.hydro:::.is_hydrodata(list(a = 1)))
  expect_false(reach.hydro:::.is_hydrodata(data.frame()))
})

test_that(".is_rainfall correctly identifies Rainfall subclasses", {
  r_daily <- .make_mock_hydrodata("Rainfall_Daily")
  r_15min <- .make_mock_hydrodata("Rainfall_15min")
  f_daily <- .make_mock_hydrodata("Flow_Daily")
  expect_true(reach.hydro:::.is_rainfall(r_daily))
  expect_true(reach.hydro:::.is_rainfall(r_15min))
  expect_false(reach.hydro:::.is_rainfall(f_daily))
})

test_that(".is_flow correctly identifies Flow subclasses", {
  f_daily <- .make_mock_hydrodata("Flow_Daily")
  f_15min <- .make_mock_hydrodata("Flow_15min")
  r_daily <- .make_mock_hydrodata("Rainfall_Daily")
  expect_true(reach.hydro:::.is_flow(f_daily))
  expect_true(reach.hydro:::.is_flow(f_15min))
  expect_false(reach.hydro:::.is_flow(r_daily))
})

test_that(".assert_same_timestep catches daily/15min mismatch", {
  f_daily <- .make_mock_hydrodata("Flow_Daily")
  r_15min <- .make_mock_hydrodata("Rainfall_15min")
  expect_error(reach.hydro:::.assert_same_timestep(f_daily, r_15min),
               "timestep mismatch")
})

test_that(".assert_same_timestep passes for matching timesteps", {
  f_daily <- .make_mock_hydrodata("Flow_Daily")
  r_daily <- .make_mock_hydrodata("Rainfall_Daily")
  expect_no_error(reach.hydro:::.assert_same_timestep(f_daily, r_daily))
})

# =============================================================================
# pdm() — HydroData input errors without reach.io
# =============================================================================

test_that("pdm() errors when only one of rain/pet is HydroData", {
  skip_if(requireNamespace("reach.io", quietly = TRUE),
          "reach.io installed; error path differs")
  mock_rain <- .make_mock_hydrodata("Rainfall_Daily")
  plain_pet <- rep(2, 100)
  # .is_hydrodata will return FALSE (reach.io absent) so this passes through
  # to the length check — just ensure no crash
  expect_error(pdm(mock_rain, plain_pet, params = pdm_params()))
})

# =============================================================================
# feh_single_site() / feh_pooled() / feh_pot() — error on wrong class
# =============================================================================

test_that(".assert_flow stops with informative message for wrong class", {
  rainfall <- .make_mock_hydrodata("Rainfall_Daily")
  expect_error(reach.hydro:::.assert_flow(rainfall, "x"),
               "Flow_Daily or Flow_15min")
})

test_that(".assert_rainfall stops with informative message for wrong class", {
  flow <- .make_mock_hydrodata("Flow_Daily")
  expect_error(reach.hydro:::.assert_rainfall(flow, "x"),
               "Rainfall_Daily or Rainfall_15min")
})

# =============================================================================
# Tests requiring reach.io
# =============================================================================

skip_if_not_installed <- function() {
  skip_if_not(requireNamespace("reach.io", quietly = TRUE),
              "reach.io not installed")
}

test_that("hydrodata_values extracts values and applies quality filter", {
  skip_if_not_installed()
  mock <- .make_mock_hydrodata("Flow_Daily", n = 50L)
  vals <- hydrodata_values(mock)
  n_bad <- sum(mock$readings$quality %in% c("Missing", "Suspect"))
  expect_equal(sum(is.na(vals)), n_bad)
  expect_equal(length(vals), 50L)
})

test_that("hydrodata_values with na_quality = character(0) keeps all values", {
  skip_if_not_installed()
  mock <- .make_mock_hydrodata("Flow_Daily", n = 50L)
  vals <- hydrodata_values(mock, na_quality = character(0))
  expect_equal(sum(is.na(vals)), 0L)
})

test_that("as_rain_input errors on non-Rainfall object", {
  skip_if_not_installed()
  flow <- .make_mock_hydrodata("Flow_Daily")
  expect_error(as_rain_input(flow), "Rainfall_Daily or Rainfall_15min")
})

test_that("as_pdm_input aligns series and returns correct structure", {
  skip_if_not_installed()
  rain <- .make_mock_hydrodata("Rainfall_Daily", n = 100L)
  pet  <- .make_mock_hydrodata("Flow_Daily",     n = 100L)  # PET as flow proxy
  inp  <- as_pdm_input(rain, pet, na_quality = character(0))
  expect_named(inp, c("rain", "pet", "dates", "datetimes", "n"))
  expect_equal(inp$n, 100L)
  expect_equal(length(inp$rain), length(inp$pet))
})

test_that("as_pdm_input errors on timestep mismatch", {
  skip_if_not_installed()
  rain <- .make_mock_hydrodata("Rainfall_Daily")
  pet  <- .make_mock_hydrodata("Flow_15min")
  expect_error(as_pdm_input(rain, pet), "timestep mismatch")
})

test_that("as_amax errors on non-Flow object", {
  skip_if_not_installed()
  rain <- .make_mock_hydrodata("Rainfall_Daily")
  expect_error(as_amax(rain), "Flow_Daily or Flow_15min")
})

test_that("as_amax returns numeric vector with water_years attribute", {
  skip_if_not_installed()
  flow <- .make_mock_hydrodata("Flow_Daily", n = 365 * 5L,
                               start = "2015-10-01")
  am   <- as_amax(flow, min_years = 3L)
  expect_true(is.numeric(am))
  expect_true(!is.null(attr(am, "water_years")))
  expect_equal(length(am), length(attr(am, "water_years")))
})

test_that("as_pot returns data.table with n_years attribute", {
  skip_if_not_installed()
  flow <- .make_mock_hydrodata("Flow_Daily", n = 365 * 3L,
                               start = "2018-10-01")
  # Set values high enough to get some exceedances
  flow$readings[, value := abs(value) * 5 + 3]
  pot <- as_pot(flow, threshold = 5)
  expect_true(data.table::is.data.table(pot))
  expect_gt(attr(pot, "n_years"), 0)
})

test_that("as_obs_flow aligns flow to supplied datetimes", {
  skip_if_not_installed()
  flow  <- .make_mock_hydrodata("Flow_Daily", n = 100L)
  dts   <- flow$readings$dateTime[1:80]  # subset of datetimes
  obs   <- as_obs_flow(flow, datetimes = dts, na_quality = character(0))
  expect_equal(length(obs), 80L)
})

test_that("hydrodata_provenance returns named list with expected fields", {
  skip_if_not_installed()
  flow <- .make_mock_hydrodata("Flow_Daily")
  prov <- hydrodata_provenance(flow)
  expect_named(prov, c("class", "parameter", "period_name", "from_date",
                        "to_date", "n_rows", "n_measures", "downloaded_at"))
  expect_equal(prov$class, "Flow_Daily")
})

test_that("pdm() accepts HydroData objects end-to-end", {
  skip_if_not_installed()
  set.seed(1)
  rain <- .make_mock_hydrodata("Rainfall_Daily", n = 200L)
  pet  <- .make_mock_hydrodata("Flow_Daily",     n = 200L)
  rain$readings[, value := pmax(0, rnorm(200, 3, 2))]
  pet$readings[,  value := pmax(0, rnorm(200, 2, 0.3))]
  p   <- pdm_params(dist = "pareto", cmax = 300)
  res <- pdm(rain, pet, params = p)
  expect_s3_class(res, "ReachHydroResult")
  expect_equal(nrow(res), 200L)
})

test_that("feh_single_site() accepts a Flow HydroData object", {
  skip_if_not_installed()
  flow <- .make_mock_hydrodata("Flow_Daily", n = 365 * 15L,
                               start = "2005-10-01")
  flow$readings[, value := pmax(5, abs(rnorm(365 * 15, 80, 30)))]
  fit <- feh_single_site(flow, n_boot = 0L)
  expect_s3_class(fit, "FehStatFit")
})
