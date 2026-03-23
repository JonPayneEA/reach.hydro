# =============================================================================
# Tool:        reach.hydro test suite — PDM state updating
# Description: Unit tests for update_pdm_states(). Covers both methods,
#              edge cases, and integration with pdm().
# Author:      Forecasting and Warning Team
# Created:     2026-03-21
# =============================================================================

library(testthat)
library(reach.hydro)

# ---- shared fixtures --------------------------------------------------------

set.seed(99)
.n_win  <- 60L                            # 60-timestep assimilation window
.tsq_w  <- seq_len(.n_win)
.seas_w <- 0.5 + 0.5 * cos(2 * pi * (.tsq_w - 30) / 365)
.rain_w <- pmax(0, stats::rnorm(.n_win, mean = 2.5 * .seas_w, sd = 3))
.pet_w  <- pmax(0, 2 * (1 - .seas_w) + stats::rnorm(.n_win, 0, 0.2))

.p <- pdm_params(dist = "pareto", cmax = 300, b = 0.4, St = 15,
                 kg = 150, kb = 180, m = 1, k1 = 5, k2 = 8)

# Simulate obs_q from a "true" run with slightly different params
.p_true <- pdm_params(dist = "pareto", cmax = 320, b = 0.35, St = 12,
                      kg = 140, kb = 160, m = 1, k1 = 5, k2 = 8)
.obs_q  <- pdm(.rain_w, .pet_w, params = .p_true)$Q

# =============================================================================
# Construction and class
# =============================================================================

test_that("update_pdm_states returns a PdmStateUpdate for method = window", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = "window")
  expect_s3_class(upd, "PdmStateUpdate")
  expect_equal(upd$method, "window")
})

test_that("update_pdm_states returns a PdmStateUpdate for method = inversion", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = "inversion")
  expect_s3_class(upd, "PdmStateUpdate")
  expect_equal(upd$method, "inversion")
})

test_that("PdmStateUpdate has correct named elements", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p)
  expect_named(upd, c("S0", "Sg0", "Ss10", "Ss20",
                       "method", "hindcast", "diagnostics"))
})

# =============================================================================
# State validity
# =============================================================================

test_that("window method: S0 is in [0, Smax]", {
  upd  <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = "window")
  Smax <- capacity_smax("pareto", .p)
  expect_gte(upd$S0,  0)
  expect_lte(upd$S0,  Smax + 1e-9)
})

test_that("all returned stores are non-negative", {
  for (m in c("window", "inversion")) {
    upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = m)
    expect_gte(upd$S0,   0, label = paste(m, "S0"))
    expect_gte(upd$Sg0,  0, label = paste(m, "Sg0"))
    expect_gte(upd$Ss10, 0, label = paste(m, "Ss10"))
    expect_gte(upd$Ss20, 0, label = paste(m, "Ss20"))
  }
})

test_that("inversion with Q_obs == Q_sim leaves stores almost unchanged", {
  # Use the same params for obs, so Q_obs ≈ Q_sim
  obs_same <- pdm(.rain_w, .pet_w, params = .p)$Q
  upd <- update_pdm_states(.rain_w, .pet_w, obs_same, .p, method = "inversion")
  d   <- upd$diagnostics
  expect_equal(d$scale_factor, 1, tolerance = 1e-9)
  expect_equal(upd$Sg0,  d$Sg_before,  tolerance = 1e-9)
  expect_equal(upd$Ss10, d$Ss1_before, tolerance = 1e-9)
  expect_equal(upd$Ss20, d$Ss2_before, tolerance = 1e-9)
})

test_that("inversion increases Sg when Q_obs > Q_sim", {
  obs_high <- .obs_q * 2
  upd <- update_pdm_states(.rain_w, .pet_w, obs_high, .p, method = "inversion")
  expect_gt(upd$Sg0, upd$diagnostics$Sg_before)
})

test_that("inversion decreases Sg when Q_obs < Q_sim", {
  obs_low <- .obs_q * 0.4
  upd <- update_pdm_states(.rain_w, .pet_w, obs_low, .p, method = "inversion")
  expect_lt(upd$Sg0, upd$diagnostics$Sg_before)
})

test_that("inversion: Sg capped at Sg_max when Sg_max > 0", {
  p_cap <- pdm_params(dist = "pareto", cmax = 300, b = 0.4, St = 15,
                      recharge_type = "demand",
                      kg = 150, kb = 180, m = 1, k1 = 5, k2 = 8,
                      Sg_max = 50, alpha = 0.4, beta = 1, q_sat = 2)
  obs_high <- .obs_q * 100   # force large scale
  upd <- update_pdm_states(.rain_w, .pet_w, obs_high, p_cap,
                            method = "inversion")
  expect_lte(upd$Sg0, 50 + 1e-9)
})

# =============================================================================
# Diagnostics
# =============================================================================

test_that("diagnostics contains expected fields", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p)
  expected <- c("Q_obs", "Q_sim_window", "Q_implied_update",
                "Sg_before", "Sg_after", "Ss1_before", "Ss2_before",
                "Ss1_after", "Ss2_after", "scale_factor")
  expect_true(all(expected %in% names(upd$diagnostics)))
})

test_that("hindcast is a ReachHydroResult with correct window length", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p)
  expect_s3_class(upd$hindcast, "ReachHydroResult")
  expect_equal(nrow(upd$hindcast), .n_win)
})

test_that("window method: Q_implied_update equals Q_sim_window", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = "window")
  d   <- upd$diagnostics
  expect_equal(d$Q_implied_update, d$Q_sim_window)
})

# =============================================================================
# Integration: updated states can be passed directly to pdm()
# =============================================================================

test_that("updated states work as pdm() initial conditions without error", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = "inversion")
  n_fcast <- 48L
  rain_f  <- pmax(0, stats::rnorm(n_fcast, 2, 3))
  pet_f   <- pmax(0, rep(2, n_fcast))
  expect_no_error(
    pdm(rain_f, pet_f, params = .p,
        S0   = upd$S0,
        Sg0  = upd$Sg0,
        Ss10 = upd$Ss10,
        Ss20 = upd$Ss20)
  )
})

test_that("forecast starting from updated states produces non-negative flow", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p, method = "inversion")
  rain_f <- pmax(0, stats::rnorm(48L, 2, 3))
  pet_f  <- pmax(0, rep(2, 48L))
  fcast  <- pdm(rain_f, pet_f, params = .p,
                S0 = upd$S0, Sg0 = upd$Sg0, Ss10 = upd$Ss10, Ss20 = upd$Ss20)
  expect_true(all(fcast$Q >= -1e-9))
})

# =============================================================================
# Input validation
# =============================================================================

test_that("errors on mismatched rain/pet lengths", {
  expect_error(
    update_pdm_states(.rain_w, .pet_w[1:10], .obs_q, .p),
    "same length"
  )
})

test_that("errors on mismatched rain/obs_q lengths", {
  expect_error(
    update_pdm_states(.rain_w, .pet_w, .obs_q[1:10], .p),
    "same length"
  )
})

test_that("errors on single-element input", {
  expect_error(
    update_pdm_states(1, 1, 1, .p),
    "at least 2 timesteps"
  )
})

test_that("warns and returns window states when Q_obs is NA at forecast origin", {
  obs_na       <- .obs_q
  obs_na[.n_win] <- NA
  expect_warning(
    upd <- update_pdm_states(.rain_w, .pet_w, obs_na, .p, method = "inversion"),
    regexp = "NA or negative"
  )
  # Should still return a valid PdmStateUpdate
  expect_s3_class(upd, "PdmStateUpdate")
})

# =============================================================================
# print method
# =============================================================================

test_that("print.PdmStateUpdate runs without error", {
  upd <- update_pdm_states(.rain_w, .pet_w, .obs_q, .p)
  expect_output(print(upd), "PdmStateUpdate")
})
