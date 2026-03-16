# =============================================================================
# data-raw/prepare_example_data.R
# Description: Scripts to prepare any built-in example datasets for reach.hydro
#              Run manually; outputs saved to data/ via usethis::use_data().
# Author:      Forecasting and Warning Team
# =============================================================================

# Example: synthetic 3-year daily forcing dataset for vignettes and tests
# Uncomment and run to regenerate:

# set.seed(42)
# n      <- 365 * 3
# t_seq  <- seq_len(n)
# season <- 0.5 + 0.5 * cos(2 * pi * (t_seq - 30) / 365)
# reach_hydro_example <- data.table::data.table(
#   date = seq.Date(as.Date("2020-10-01"), by = "day", length.out = n),
#   rain = pmax(0, stats::rnorm(n, mean = 2.5 * season, sd = 3)),
#   pet  = pmax(0, 3 * (1 - season) + stats::rnorm(n, 0, 0.3))
# )
# usethis::use_data(reach_hydro_example, overwrite = TRUE)
