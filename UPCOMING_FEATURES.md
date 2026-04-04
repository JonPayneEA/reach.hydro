
# Upcoming Features - reach.hydro

Planned enhancements and new capabilities for `reach.hydro`.

**Primary use: flood forecasting.** Features in §1–§5 directly support
operational and pre-operational forecasting workflows and are prioritised
accordingly.

**Secondary use: design hydrology.** Features in §6–§7 support design flood
estimation (ReFH2, FSR, DDF) and are planned but lower priority.

Items within each section run roughly top-to-bottom by priority.

---

## PRIMARY - Flood Forecasting

---

## 1. PDM Performance and Robustness

### 1.1 Rcpp acceleration of the PDM loop

The per-timestep update loop in `pdm_core.R` is the primary performance
bottleneck for operational runs and long calibration jobs. Rewrite the inner
loop in C++ via Rcpp, keeping the R interface identical. Expected 20–50×
speed-up based on profiling. Tier 1 tests must pass unchanged after migration.

### 1.2 Multi-objective calibration

Extend `calibrate_pdm()` to support Pareto-optimal calibration against multiple
objectives simultaneously (e.g. KGE on high flows + NSE on low flows). Returns
a Pareto front rather than a single parameter set, enabling uncertainty-aware
operational use.

### 1.3 Dual-zone and multi-zone PDM

The standard PDM uses a single soil moisture store. Many UK catchments exhibit
mixed responses (e.g. fast-responding urban/impervious areas alongside slower
rural areas) that a single-zone structure cannot represent well. This extension
partitions the catchment into two or more zones, each with independent PDM
parameters, and combines their outputs as a weighted sum.

```r
# Two-zone model: 30% urban, 70% rural
params_urban <- pdm_params(cmax = 20,  b = 0.3, kg = 10,  ...)
params_rural <- pdm_params(cmax = 150, b = 1.2, kg = 120, ...)

result <- pdm_multizone(
  rain, pet,
  zones = list(urban = params_urban, rural = params_rural),
  weights = c(urban = 0.30, rural = 0.70)
)

# Weights can also be calibrated jointly
cal <- calibrate_pdm_multizone(
  rain, pet, obs_flow,
  n_zones = 2,
  fixed_weights = FALSE
)
```

Scope:
- `pdm_multizone()` runs each zone's PDM loop independently then combines
  direct runoff and baseflow contributions by zone weight.
- Zone weights optionally treated as free parameters during calibration.
- Compatible with §1.1 (Rcpp loop) and §1.2 (multi-objective calibration).
- Initial implementation targets two zones (dual-zone); generalised to
  n-zones in a follow-up.

### 1.4 Monte Carlo uncertainty propagation

Add `pdm_uncertainty()` to sample from a parameter distribution (or bootstrap
calibration residuals) and return an ensemble of flow series with quantile
bands, supporting probabilistic forecasting workflows.

```r
ens <- pdm_uncertainty(rain, pet, cal$params, n = 500, method = "bootstrap")
ens$quantiles   # data.table: q05, q25, q50, q75, q95 at each timestep
```

### 1.5 State updating — removing the burn-in requirement *(implemented)*

`update_pdm_states()` runs the PDM over a recent assimilation window to produce
calibrated initial conditions (`S0`, `Sg0`, `Ss10`, `Ss20`) for a forecast run,
removing the need for a warm-up period. Two methods are provided:

- **`"window"`** — run the model over the last N timesteps of observed rain/PET;
  use the final states as initial conditions. Effective for windows ≥ 30 days.
- **`"inversion"`** — as `"window"`, then algebraically adjust the groundwater
  store `Sg` and surface routing stores `Ss1`/`Ss2` so that the implied
  instantaneous outflow matches the observed flow at the forecast origin.

```r
upd <- update_pdm_states(
  rain  = tail(rain_recent, 30L),
  pet   = tail(pet_recent,  30L),
  obs_q = tail(obs_q_recent, 30L),
  params  = cal$params,
  method  = "inversion"
)

# Pass updated states directly into the forecast run
fcast <- pdm(rain_fcast, pet_fcast, params = cal$params,
             S0 = upd$S0, Sg0 = upd$Sg0, Ss10 = upd$Ss10, Ss20 = upd$Ss20)
```

See [design/pdm-state-updating.md](design/pdm-state-updating.md) for the full
design, including the inversion algebra, window length guidelines, and planned
extensions (ensemble method, reach.validate integration, write_pdm_run() hook).

---

## 2. Operational Data Ingestion

### 2.1 NRFA API client *(high priority)*

Fetch gauged flow time series and AMAX/POT records on demand from the NRFA web
service by station number. Removes the need to download bulk archive files
during model setup and calibration. Returns reach.io-compatible output.

```r
flow_ts <- nrfa_gdf(station = 25001, from = "1990-10-01", to = "2020-09-30")
amax    <- nrfa_amax(station = 25001)
```

### 2.2 HiFlows-UK zip ingestion

The NRFA distributes the HiFlows-UK dataset as a zip archive containing one
file per gauging station. Each station file bundles AMAX and POT records with
station metadata in a fixed-format header/data layout.

Planned API:

```r
# Read all stations from the downloaded zip
hf <- read_hiflows_zip("HiFlowsUK_v4.zip")

# Returns a named list of HiFlowsStation objects
# hf[["25001"]]$amax  - data.table of annual maxima
# hf[["25001"]]$pot   - data.table of peaks over threshold
# hf[["25001"]]$meta  - list of station metadata

# Direct pass-through to FEH functions for index flood / pooled analysis
fit_pool <- feh_pooled(
  subject_amax = hf[["25001"]]$amax$flow,
  donor_list   = lapply(hf[c("25002", "25003", "25006")], \(s) s$amax$flow)
)
```

Scope:
- Parse `.am` (AMAX) and `.pt` (POT) formats within the zip without
  extracting to disk (using `unz()` connections).
- Return `HiFlowsStation` S3 objects carrying `amax`, `pot`, and `meta` slots.
- Validate water year alignment and flag incomplete records.
- Helper `hiflows_pool_candidates()` to filter donor stations by distance or
  QMED similarity.

### 2.3 WINFAP-FEH file ingestion

Read individual `.am` and `.pt` files from locally archived WINFAP-FEH exports.
Feeds directly into `feh_single_site()` and `feh_pot()` for catchment
verification during model setup.

---

## 3. Flow Statistics and Baseflow

### 3.1 Baseflow separation *(implemented)*

`baseflow_separate()` separates an observed flow series into baseflow and
quickflow components using a recursive digital filter. Works at any timestep
(daily, hourly, 15-minute).

- **Lyne-Hollick** (one-parameter recursive filter, default `alpha = 0.925`)
- **Boughton-Eckhardt** (two-parameter: `k = 0.975`, `C = 0.1`)
- Multiple passes (forward/backward) to remove phase shift (default 3 passes)

```r
bf <- baseflow_separate(flow, method = "lyne_hollick", alpha = 0.925)
bf$baseflow   # numeric vector
bf$quickflow  # numeric vector
bf$bfi        # Baseflow Index for the period
```

Provides an observed-data alternative to `baseflow_index()` (which requires
pre-separated Qb from a model run).

### 3.2 Flow percentile and deficit analysis *(implemented)*

`flow_statistics.R` now includes:

- `flow_deficit()` - volume and duration of low-flow spells below a threshold;
  works with daily or sub-daily data; returns duration in real time units when
  dates are supplied
- `flow_recession()` - automatic recession curve fitting (`Q = Q0 * exp(-t/k)`)
  across identified recession limbs; returns median k and per-event statistics
- `q_n_day()` - rolling n-day minimum/maximum flow (Q7 for low-flow indices);
  returns scalar or annual values by water year when dates supplied
- `monthly_flow_stats()` - mean, median, Q10, Q90, max by calendar month;
  works at any timestep

### 3.3 Rainfall time series analysis *(implemented)*

`R/rainfall_statistics.R` provides standalone rainfall analysis tools:

- `rainfall_events()` - storm event extraction; identifies discrete events
  separated by dry gaps; characterises each event by duration, total depth,
  peak and mean intensity; works at any timestep
- `api()` - antecedent precipitation index (exponential decay memory):
  `API[t] = k * API[t-1] + P[t]`; used for antecedent wetness estimation
- `idf_empirical()` - empirical intensity-duration-frequency table; extracts
  rolling-window maximum depths for specified durations; returns annual maxima
  by water year when dates supplied

---

## 4. Model Validation - reach.validate Integration

Provide a compatibility layer for the `reach.validate` package (analogous to
`reach_io_compat.R`) so that PDM simulation outputs can be passed directly into
`reach.validate` validation workflows without manual extraction.

```r
# Coerce a PDM simulation result to a reach.validate ModelRun object
val_run <- as_model_run(pdm_result, observed = obs_flow, warmup = 365)

# Run the full reach.validate metric suite and return a tidy report
report  <- validate(val_run)
report$metrics   # data.frame: NSE, KGE, PBIAS, FAR, RVE, …
report$plots     # list of ggplot2 objects (hydrograph, FDC, scatter)

# Compare across distributions or parameter sets
compare_runs <- lapply(dist_list, \(d) {
  res <- pdm(rain, pet, calibrate_pdm(rain, pet, obs_flow, dist = d)$params)
  as_model_run(res, observed = obs_flow)
})
validation_table(compare_runs)
```

Scope:
- Detect `reach.validate` at runtime via `requireNamespace()` - package
  remains fully functional without it (same pattern as `reach.io`).
- `as_model_run()` aligns simulated and observed series by date, applies the
  warmup mask, and attaches metadata (catchment ID, model type, run date).
- `split_validation()` - calibration/validation period split with metrics
  reported for each period separately.
- Feed FEH flood frequency fit objects into `reach.validate` benchmark
  comparisons against observed AMAX records.

---

## 5. Infrastructure

### 5.1 Rcpp-ready package structure

Scaffold the `src/` directory and `Makevars` so that Rcpp compilation is
supported cleanly, ahead of §1.1. Includes CI pipeline updates for compiled
code.

### 5.2 S7 class migration

Replace the current `PdmParams` and `ReachHydroResult` S3 stubs with full S7
class definitions once the S7 package is available on CRAN. S7 provides formal
property validation, inheritance, and method dispatch, removing boilerplate
from the existing `pdm_validate_params()` calls.

### 5.3 reach.io HydroData round-trip for PDM output

Extend `reach_io_compat.R` so that `ReachHydroResult` objects can be coerced
to and from `Flow_Daily` / `Flow_15min` HydroData objects. Enables PDM
simulations to be stored and retrieved via the reach.io data layer without
manual extraction.

### 5.4 Thiessen polygon weights

Implement `thiessen_weights()` to compute areal weighting for rain gauges from
gauge and catchment boundary coordinates, replacing the current placeholder in
`areal_rainfall()`. Supports operational areal rainfall estimation from raingauge
networks.

### 5.5 PDM run storage — Parquet + manifest convention

Implement `write_pdm_run()` and `read_pdm_run()` helpers (Tier 2) to persist
`ReachHydroResult` objects to disk in a governed, reproducible format. Uses
Apache Parquet for the time series and a companion `manifest.json` for provenance
and parameter metadata.

Addresses the scale problem: 30-year 15-minute hindcasts (~1.05M rows) and
calibration jobs that must remain reproducible and auditable under the Data &
Digital Asset Governance Framework.

See [design/pdm-run-storage.md](design/pdm-run-storage.md) for the full design,
including directory structure, manifest schema, calibration conventions, and
alignment with §1.2, §1.4, §5.2, and §5.3.

---

## SECONDARY - Design Hydrology

*Lower priority. These features support design flood estimation workflows
(ReFH2, FSR, DDF) rather than operational forecasting.*

---

## 6. FEH Design Methods

### 6.1 ReFH2 regression coefficient verification

The ReFH2 regression coefficients in `feh_refh2.R` are currently flagged as
AI-generated and unverified. Replace with values from the published ReFH2
report (Kjeldsen et al., 2008) and add validation tests against the worked
examples in that report. **Must be resolved before any design use.**

### 6.2 FEH DDF regional parameter sets

`feh_ddf()` currently uses UK national-average DDF parameters. Add the
regional coefficient sets from the FEH CD-ROM so that estimates use the
appropriate region for a given catchment.

```r
feh_ddf(..., region = "north_west")
# or derive region automatically from BNG coordinates
feh_ddf(..., easting = 358000, northing = 387000)
```

### 6.3 WINFAP benchmark suite

Compare `feh_single_site()` GLO/GEV outputs against WINFAP-FEH on a set of
reference AMAX series to confirm numerical equivalence before design use
(see Known TODOs in README).

### 6.4 L-moment ratio diagram

Add `plot_lmrd()` to visualise sample L-skewness vs L-kurtosis alongside
theoretical curves for GLO, GEV, GNO, and PE3. Uses the existing
`lmrd_theoretical()` helper and supports both base R and ggplot2 output.

### 6.5 FEH catchment descriptor lookup

Integrate with the FEH Web Service to retrieve standard catchment descriptors
(AREA, BFIHOST, SAAR, FARL, URBEXT, …) by outlet coordinates or station
number. Primarily useful for design hydrology parameter derivation.

```r
descs <- feh_descriptors(easting = 358000, northing = 387000)
p <- refh2_params(descs)
```

---

## 7. Documentation

### 7.1 Worked example vignettes

Priority order reflects the primary/secondary split above:

1. **PDM calibration and operational forecasting** - from raw inputs through
   calibration, validation via `reach.validate`, and ensemble uncertainty.
2. **FEH Vol 3 flood frequency** - single-site and pooled analysis with
   HiFlows-UK data ingestion, for catchment verification.
3. **Design flood estimation** *(secondary)* - ReFH2 and FSR with FEH Vol 4
   storms.

### 7.2 Expanded test coverage

Raise Tier 2 coverage to 70% (currently lower than the Tier 1 target). Prioritise
tests for forecasting-critical paths (PDM, calibration, flow statistics) before
design hydrology methods.

---

*Last updated: 2026-03-19 - restructured to prioritise flood forecasting over design hydrology*
