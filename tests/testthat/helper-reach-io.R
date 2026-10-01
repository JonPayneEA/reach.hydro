# =============================================================================
# reach.hydro — reach.io test fixtures
# Description: Builds real reach.io S7 objects. Plain list mocks reproduce
#              neither reach.io's namespaced class vector nor its S7 dispatch,
#              so tests built on them pass while real objects fail. Every
#              caller must skip_if_not_installed("reach.io").
# Author:      Forecasting and Warning Team
# Created:     2026-10-01
# =============================================================================

.make_io_hd <- function(class_name,
                        values  = NULL,
                        n       = 100L,
                        start   = "2020-01-01",
                        quality = NULL) {
  if (is.null(values)) values <- pmax(0, stats::rnorm(n, mean = 5, sd = 2))
  n     <- length(values)
  dates <- seq.Date(as.Date(start), by = "day", length.out = n)
  if (is.null(quality)) quality <- rep("Good", n)

  readings <- data.table::data.table(
    dateTime         = as.POSIXct(dates, tz = "UTC"),
    date             = dates,
    value            = values,
    measure_notation = rep("m", n),
    quality          = quality
  )

  ctor <- getExportedValue("reach.io", paste0("Flode", class_name))
  ctor(readings  = readings,
       from_date = as.character(dates[1L]),
       to_date   = as.character(dates[n]))
}
