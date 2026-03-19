# Upcoming Features — reach.hydro

Planned enhancements and new capabilities for `reach.hydro`. Items are grouped
by theme. Priority within each group runs roughly top-to-bottom.

---

## 1. Data Ingestion

### 1.1 HiFlows-UK zip ingestion *(requested)*

The NRFA distributes the HiFlows-UK dataset as a zip archive containing one
file per gauging station. Each station file bundles AMAX and POT records with
station metadata in a fixed-format header/data layout.

Planned API:

```r
# Read all stations from the downloaded zip
hf <- read_hiflows_zip("HiFlowsUK_v4.zip")

# Returns a named list of HiFlowsStation objects
# hf[["25001"]]$amax  — data.table of annual maxima
# hf[["25001"]]$pot   — data.table of peaks over threshold
# hf[["25001"]]$meta  — list of station metadata

# Direct pass-through to existing FEH functions
fit <- feh_single_site(hf[["25001"]]$amax$flow, dist = "glo")
fit_pool <- feh_pooled(
  subject_amax = hf[["25001"]]$amax$flow,
  donor_list   = lapply(hf[c("25002", "25003", "25006")], \(s) s$amax$flow)
)
```

Scope:
- Parse `.am` (AMAX) and `.pt` (POT) file formats within the zip without
  extracting to disk (using `unz()` connections).
- Return `HiFlowsStation` S3 objects carrying `amax`, `pot`, and `meta` slots.
- Validate water year alignment and flag incomplete records.
- Helper `hiflows_pool_candidates()` to filter stations by distance or QMED
  similarity, ready for `feh_pooled()`.

### 1.2 WINFAP-FEH file ingestion

Read `.am` and `.pt` files produced by WINFAP-FEH directly (individual station
files, as opposed to the HiFlows zip bundle). Allows users working from locally
archived WINFAP exports to feed data straight into `feh_single_site()` and
`feh_pot()`.

### 1.3 NRFA API client

Fetch AMAX and POT series on demand from the NRFA web service by station
number, removing the need to download bulk zip archives. Wraps the NRFA JSON
API with caching support and reach.io-compatible output.

```r
amax <- nrfa_amax(station = 25001)
flow_ts <- nrfa_gdf(station = 25001, from = "1990-10-01", to = "2020-09-30")
```

---

## 2. PDM Performance and Robustness

### 2.1 Rcpp acceleration of the PDM loop

The per-timestep update loop in `pdm_core.R` is the primary performance
bottleneck for long calibration runs. Rewrite the inner loop in C++ via Rcpp,
keeping the R interface identical. Expected 20–50× speed-up based on profiling.
Tier 1 tests must pass unchanged after the migration.

### 2.2 Multi-objective calibration

Extend `calibrate_pdm()` to support Pareto-optimal calibration against multiple
objectives simultaneously (e.g. KGE on high flows + NSE on low flows). Returns
a Pareto front rather than a single parameter set, enabling uncertainty-aware
operational use.

### 2.3 Monte Carlo uncertainty propagation

Add `pdm_uncertainty()` to sample from a parameter distribution (or bootstrap
calibration residuals) and return an ensemble of flow series with quantile
bands, supporting probabilistic forecasting workflows.

---

## 3. Flow Statistics and Baseflow

### 3.1 Baseflow separation

Implement two standard digital filter algorithms for extracting baseflow from
an observed flow series:

- **Lyne-Hollick** (one-parameter recursive filter)
- **Boughton-Eckhardt** (two-parameter recursive filter)

```r
bf <- baseflow_separate(flow, method = "lyne_hollick", alpha = 0.925)
bf$baseflow   # numeric vector
bf$bfi        # Baseflow Index for the period
```

Integrates with `flow_stats()` output and provides an observed-data alternative
to `baseflow_index()` (which currently uses `bfihost`).

### 3.2 Flow percentile and deficit analysis

Extend `flow_statistics.R` with:
- `flow_deficit()` — volume and duration of low-flow spells below a threshold
- `flow_recession()` — automatic recession curve fitting (master recession curve)
- `q_n_day()` — n-day minimum/maximum flow (e.g. Q7 for low-flow indices)

---

## 4. FEH Methods

### 4.1 FEH DDF regional parameter sets

`feh_ddf()` currently uses UK national-average DDF parameters. Add the
regional coefficient sets from the FEH CD-ROM so that estimates use the
appropriate region for a given catchment.

```r
feh_ddf(..., region = "north_west")
# or derive region automatically from BNG coordinates
feh_ddf(..., easting = 358000, northing = 387000)
```

### 4.2 ReFH2 regression coefficient verification

The ReFH2 regression coefficients in `feh_refh2.R` are currently flagged as
AI-generated and unverified. Replace these with values from the published
ReFH2 report (Kjeldsen et al., 2008) and add regression-based validation tests
against the worked examples in that report.

### 4.3 WINFAP benchmark suite

Compare `feh_single_site()` GLO/GEV outputs against WINFAP-FEH on a set of
reference AMAX series to confirm numerical equivalence before operational use
(see Known TODOs in README).

### 4.4 L-moment ratio diagram

Add `plot_lmrd()` to visualise sample L-skewness vs L-kurtosis alongside
theoretical curves for GLO, GEV, GNO, and PE3. Uses the existing
`lmrd_theoretical()` helper and supports both base R and ggplot2 output.

---

## 5. Spatial / Catchment

### 5.1 Thiessen polygon weights

Implement `thiessen_weights()` to compute areal weighting for rain gauges from
gauge and catchment boundary coordinates, replacing the current placeholder in
`areal_rainfall()`.

### 5.2 FEH catchment descriptor lookup

Integrate with the FEH Web Service to retrieve standard catchment descriptors
(AREA, BFIHOST, SAAR, FARL, URBEXT, …) by providing outlet coordinates or a
station number, removing the need for manual entry.

```r
descs <- feh_descriptors(easting = 358000, northing = 387000)
p <- refh2_params(descs)
```

---

## 6. Infrastructure

### 6.1 S7 class migration

Replace the current `PdmParams` and `ReachHydroResult` S3 stubs with full S7
class definitions once the S7 package is available on CRAN. S7 provides formal
property validation, inheritance, and method dispatch that removes boilerplate
from the existing `pdm_validate_params()` calls.

### 6.2 reach.io HydroData round-trip for PDM output

Extend `reach_io_compat.R` so that `ReachHydroResult` objects can be coerced
to and from `Flow_Daily` / `Flow_15min` HydroData objects. This enables PDM
simulations to be stored and retrieved via the reach.io data layer without
manual extraction.

---

## 7. Documentation and Validation

### 7.1 Worked example vignettes

Add long-form vignettes covering the most common end-to-end workflows:

1. **PDM calibration and simulation** — from raw CSV inputs to calibrated
   parameters and goodness-of-fit diagnostics.
2. **FEH Vol 3 flood frequency** — single-site and pooled analysis with
   HiFlows-UK data ingestion.
3. **Design flood estimation** — ReFH2 and FSR methods with FEH Vol 4 storms.

### 7.2 Expanded test coverage

Raise Tier 2 coverage to 70% (currently lower than the Tier 1 target) and add
regression snapshots for FEH Vol 3 / Vol 4 outputs against published worked
examples.

---

*Last updated: 2026-03-19*
