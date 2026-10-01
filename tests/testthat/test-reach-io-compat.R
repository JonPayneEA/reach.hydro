# =============================================================================
# reach.hydro — reach.io compatibility layer tests
# Description: Tests for the reach.io compat layer. Fixtures are real reach.io
#              objects from helper-reach-io.R, so any test that builds one
#              skips when reach.io is not installed.
# Author:      Forecasting and Warning Team
# Created:     2026-02-23
# Modified:    2026-10-01 - real reach.io fixtures; Flode class names
# =============================================================================

library(testthat)
library(reach.hydro)

# =============================================================================
# CLASS DETECTION
# =============================================================================

test_that("reach.io objects carry the namespaced class names the bridge matches", {
  skip_if_not_installed("reach.io")
  flow <- .make_io_hd("Flow_Daily")
  expect_true(inherits(flow, "reach.io::FlodeFlow_Daily"))
  expect_true(inherits(flow, "reach.io::FlodeHydroData"))
})

test_that(".is_hydrodata detects reach.io objects", {
  skip_if_not_installed("reach.io")
  expect_true(reach.hydro:::.is_hydrodata(.make_io_hd("Flow_Daily")))
})

test_that(".is_hydrodata returns FALSE for plain objects", {
  expect_false(reach.hydro:::.is_hydrodata(c(1, 2, 3)))
  expect_false(reach.hydro:::.is_hydrodata(list(a = 1)))
  expect_false(reach.hydro:::.is_hydrodata(data.frame()))
})

test_that(".is_hydrodata rejects lookalikes carrying bare class names", {
  skip_if_not_installed("reach.io")
  lookalike <- structure(list(), class = c("Flow_Daily", "HydroData"))
  expect_false(reach.hydro:::.is_hydrodata(lookalike))
})

test_that(".is_rainfall correctly identifies Rainfall subclasses", {
  skip_if_not_installed("reach.io")
  expect_true(reach.hydro:::.is_rainfall(.make_io_hd("Rainfall_Daily")))
  expect_true(reach.hydro:::.is_rainfall(.make_io_hd("Rainfall_15min")))
  expect_false(reach.hydro:::.is_rainfall(.make_io_hd("Flow_Daily")))
})

test_that(".is_flow correctly identifies Flow subclasses", {
  skip_if_not_installed("reach.io")
  expect_true(reach.hydro:::.is_flow(.make_io_hd("Flow_Daily")))
  expect_true(reach.hydro:::.is_flow(.make_io_hd("Flow_15min")))
  expect_false(reach.hydro:::.is_flow(.make_io_hd("Rainfall_Daily")))
})

test_that(".assert_same_timestep catches daily/15min mismatch", {
  skip_if_not_installed("reach.io")
  expect_error(
    reach.hydro:::.assert_same_timestep(.make_io_hd("Flow_Daily"),
                                        .make_io_hd("Rainfall_15min")),
    "timestep mismatch"
  )
})

test_that(".assert_same_timestep passes for matching timesteps", {
  skip_if_not_installed("reach.io")
  expect_no_error(
    reach.hydro:::.assert_same_timestep(.make_io_hd("Flow_Daily"),
                                        .make_io_hd("Rainfall_Daily"))
  )
})

# =============================================================================
# VALIDATORS
# =============================================================================

test_that("pdm() errors when only one of rain/pet is a reach.io object", {
  skip_if_not_installed("reach.io")
  rain <- .make_io_hd("Rainfall_Daily")
  expect_error(pdm(rain, rep(2, 100), params = pdm_params()), "both must be")
})

test_that(".assert_flow names the expected classes and the class received", {
  skip_if_not_installed("reach.io")
  rainfall <- .make_io_hd("Rainfall_Daily")
  expect_error(reach.hydro:::.assert_flow(rainfall, "x"),
               "^x: expected a FlodeFlow_Daily or FlodeFlow_15min object; got reach\\.io::FlodeRainfall_Daily")
})

test_that(".assert_rainfall names the expected classes and the class received", {
  skip_if_not_installed("reach.io")
  flow <- .make_io_hd("Flow_Daily")
  expect_error(reach.hydro:::.assert_rainfall(flow, "x"),
               "^x: expected a FlodeRainfall_Daily or FlodeRainfall_15min object; got reach\\.io::FlodeFlow_Daily")
})

# =============================================================================
# EXTRACTORS AND COERCERS
# =============================================================================

test_that("hydrodata_values extracts values and applies quality filter", {
  skip_if_not_installed("reach.io")
  set.seed(42)
  q    <- sample(c("Good", "Suspect", "Missing"), 50L, replace = TRUE,
                 prob = c(0.85, 0.1, 0.05))
  flow <- .make_io_hd("Flow_Daily", n = 50L, quality = q)
  vals <- hydrodata_values(flow)
  expect_equal(sum(is.na(vals)), sum(q %in% c("Missing", "Suspect")))
  expect_equal(length(vals), 50L)
})

test_that("hydrodata_values with na_quality = character(0) keeps all values", {
  skip_if_not_installed("reach.io")
  q    <- rep(c("Good", "Missing"), 25L)
  flow <- .make_io_hd("Flow_Daily", n = 50L, quality = q)
  vals <- hydrodata_values(flow, na_quality = character(0))
  expect_equal(sum(is.na(vals)), 0L)
})

test_that("as_rain_input errors on non-Rainfall object", {
  skip_if_not_installed("reach.io")
  expect_error(as_rain_input(.make_io_hd("Flow_Daily")),
               "FlodeRainfall_Daily or FlodeRainfall_15min")
})

test_that("as_pdm_input aligns series and returns correct structure", {
  skip_if_not_installed("reach.io")
  rain <- .make_io_hd("Rainfall_Daily", n = 100L)
  pet  <- .make_io_hd("Flow_Daily",     n = 100L)  # PET as flow proxy
  inp  <- as_pdm_input(rain, pet, na_quality = character(0))
  expect_named(inp, c("rain", "pet", "dates", "datetimes", "n"))
  expect_equal(inp$n, 100L)
  expect_equal(length(inp$rain), length(inp$pet))
})

test_that("as_pdm_input errors on timestep mismatch", {
  skip_if_not_installed("reach.io")
  expect_error(as_pdm_input(.make_io_hd("Rainfall_Daily"),
                            .make_io_hd("Flow_15min")),
               "timestep mismatch")
})

test_that("as_amax errors on non-Flow object", {
  skip_if_not_installed("reach.io")
  expect_error(as_amax(.make_io_hd("Rainfall_Daily")),
               "FlodeFlow_Daily or FlodeFlow_15min")
})

test_that("as_amax returns numeric vector with water_years attribute", {
  skip_if_not_installed("reach.io")
  flow <- .make_io_hd("Flow_Daily", n = 365L * 5L, start = "2015-10-01")
  am   <- as_amax(flow, min_years = 3L)
  expect_true(is.numeric(am))
  expect_true(!is.null(attr(am, "water_years")))
  expect_equal(length(am), length(attr(am, "water_years")))
})

test_that("as_pot returns data.table with n_years attribute", {
  skip_if_not_installed("reach.io")
  set.seed(7)
  vals <- abs(stats::rnorm(365L * 3L, mean = 5, sd = 2)) * 5 + 3
  flow <- .make_io_hd("Flow_Daily", values = vals, start = "2018-10-01")
  pot  <- as_pot(flow, threshold = 5)
  expect_true(data.table::is.data.table(pot))
  expect_gt(attr(pot, "n_years"), 0)
})

test_that("as_obs_flow aligns flow to supplied datetimes", {
  skip_if_not_installed("reach.io")
  flow <- .make_io_hd("Flow_Daily", n = 100L)
  dts  <- flow@readings$dateTime[1:80]
  obs  <- as_obs_flow(flow, datetimes = dts, na_quality = character(0))
  expect_equal(length(obs), 80L)
})

test_that("hydrodata_provenance reports the class without its namespace", {
  skip_if_not_installed("reach.io")
  prov <- hydrodata_provenance(.make_io_hd("Flow_Daily"))
  expect_named(prov, c("class", "parameter", "period_name", "from_date",
                       "to_date", "n_rows", "n_measures", "downloaded_at"))
  expect_equal(prov$class, "FlodeFlow_Daily")
})

# =============================================================================
# END-TO-END
# =============================================================================

test_that("pdm() accepts reach.io objects end-to-end", {
  skip_if_not_installed("reach.io")
  set.seed(1)
  rain <- .make_io_hd("Rainfall_Daily", values = pmax(0, stats::rnorm(200, 3, 2)))
  pet  <- .make_io_hd("Flow_Daily",     values = pmax(0, stats::rnorm(200, 2, 0.3)))
  res  <- pdm(rain, pet, params = pdm_params(dist = "pareto", cmax = 300))
  expect_s3_class(res, "ReachHydroResult")
  expect_equal(nrow(res), 200L)
})

test_that("feh_single_site() accepts a reach.io Flow object", {
  skip_if_not_installed("reach.io")
  set.seed(3)
  vals <- pmax(5, abs(stats::rnorm(365L * 15L, 80, 30)))
  flow <- .make_io_hd("Flow_Daily", values = vals, start = "2005-10-01")
  fit  <- feh_single_site(flow, n_boot = 0L)
  expect_s3_class(fit, "FehStatFit")
})
