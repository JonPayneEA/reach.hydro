# reach.hydro

<img src="man/figures/logo.svg" align="right" width=160/>

**Core hydrological calculations for the REACH framework**  
Forecasting and Warning Team | Tier 1/2 | Version 0.2.0

---

## Overview

`reach.hydro` provides a complete set of hydrological calculation tools for
the Forecasting and Warning team, covering rainfall-runoff modelling, flood
frequency estimation, and design flood methods. It is the hydrology module of
the REACH framework within the flode meta-package.

## Methods covered

| Method | Standard | Key functions |
|--------|----------|---------------|
| PDM rainfall-runoff | Moore (2007) | `pdm()`, `calibrate_pdm()` |
| Statistical flood frequency | FEH Vol 3 | `feh_single_site()`, `feh_pooled()`, `feh_pot()` |
| Rainfall frequency (DDF) | FEH Vol 4 | `feh_ddf()`, `feh_design_storm()` |
| ReFH2 rainfall-runoff | Kjeldsen (2007) | `refh2_params()`, `refh2_run()` |
| FSR/FEH unit hydrograph | NERC (1975) | `fsr_params()`, `fsr_run()` |
| Flow statistics | — | `flow_stats()`, `flow_duration_curve()` |
| Flood frequency (convenience) | FEH | `fit_glo()`, `fit_gev()` |

---

## Installation

```r
# From the team Git repository
devtools::install_git("https://github.com/forecasting-warning/reach.hydro")

# Restore the renv lockfile first in a new project
renv::restore()
```

---

## Quick start

### PDM rainfall-runoff

```r
library(reach.hydro)

p   <- pdm_params(dist = "pareto", cmax = 350, b = 0.4, St = 20, kg = 150, ks = 8)
res <- pdm(rain, pet, params = p)
summary(res)

# Performance against observed flow
gof_metrics(obs_q, res$Q, warmup = 365)

# Compare all six distributions
tbl <- compare_distributions(rain, pet, obs_q = obs_q)

# Calibrate
cal <- calibrate_pdm(rain, pet, obs_q, dist = "pareto", metric = "kge")
cal$kge
```

### FEH Vol 3 — statistical flood frequency

```r
# Single-site GLO fit with bootstrap confidence intervals
fit <- feh_single_site(amax, dist = "glo", n_boot = 1000)
fit$growth_curve

# Pooled analysis with donor sites
fit_pool <- feh_pooled(
  subject_amax = amax,
  donor_list   = list(site_a = amax_a, site_b = amax_b),
  dist         = "glo"
)
fit_pool$growth_curve

# POT analysis
peaks <- peaks_over_threshold(flow, dates, threshold = 150)
fit_pot <- feh_pot(peaks$peak_flow, threshold = 150, n_years = 30)
```

### FEH Vol 4 — rainfall frequency and design storm

```r
# T-year rainfall depth at multiple durations
feh_ddf(duration_hr = c(1, 2, 4, 6, 12, 24),
        return_period = 100,
        rmed_1h = 12, rmed_1d = 38, saar = 900)

# Design storm hyetograph (15-minute timestep, summer profile)
storm <- feh_design_storm(duration_hr = 4, return_period = 100,
                           rmed_1h = 12, rmed_1d = 38, saar = 900,
                           dt_min = 15, season = "summer")
```

### ReFH2

```r
# Derive parameters from catchment descriptors
p <- refh2_params(area = 250, bfihost = 0.45, saar = 850, farl = 0.98)

# Run model on design storm
res <- refh2_run(storm, p)
refh2_summary(res, area_km2 = 250)
```

### FSR/FEH unit hydrograph

```r
p   <- fsr_params(area = 150, bfihost = 0.4, saar = 800)
res <- fsr_run(storm, p, m5_60min = feh_m5_ratios(12, 38)["m5_60min"])
max(res$Q_mm)
```

---

## Package structure

| File | Contents |
|------|----------|
| `R/pdm_classes.R` | `PdmParams` and `ReachHydroResult` S3 classes |
| `R/pdm_numerical.R` | Numerical integration and root-finding helpers |
| `R/pdm_distributions.R` | Six capacity distributions + dispatch functions |
| `R/pdm_stores.R` | AET, soil drainage, groundwater, surface store updates |
| `R/pdm_core.R` | Main `pdm()` function |
| `R/pdm_metrics.R` | NSE, KGE, PBIAS, FAR, `gof_metrics()` |
| `R/pdm_calibration.R` | `calibrate_pdm()`, `compare_distributions()` |
| `R/feh_lmoments.R` | L-moments engine: sample L-moments, GLO/GEV/GNO/PE3 fitting |
| `R/feh_vol3_statistical.R` | FEH Vol 3: single-site, pooled, POT flood frequency |
| `R/feh_vol4_rainfall.R` | FEH Vol 4: DDF model, design storm, M5 ratios, ARF |
| `R/feh_refh2.R` | ReFH2: parameter estimation, model run, summary |
| `R/feh_fsr_uh.R` | FSR/FEH: unit hydrograph, percentage runoff, full run |
| `R/flood_frequency.R` | Convenience wrappers: `fit_glo()`, `fit_gev()`, `return_period_flow()` |
| `R/flow_statistics.R` | `flow_stats()`, BFI, FDC, annual maxima, POT extraction |
| `R/unit_conversions.R` | mm ↔ m³/s conversions |
| `R/catchment_aggregation.R` | Areal rainfall, weighted spatial means |

---

## PDM capacity distributions

| Distribution | Parameters | Closed-form Smax | Closed-form c* |
|---|---|:---:|:---:|
| Pareto (default) | cmin, cmax, b | ✓ | ✓ |
| Uniform | cmin, cmax | ✓ | ✓ |
| Exponential | cmax | ✓ | ✓ |
| Generalised Logistic | cmin, cmax, b | — | — |
| Normal | mu_c, sigma_c | — | — |
| Log-Normal | mu_lnc, sigma_lnc | ✓ | — |

---

## Governance

| Attribute | Value |
|-----------|-------|
| Version | 0.2.0 |
| Tier | 1 (PDM core, metrics) / 2 (FEH methods, calibration, flow stats) |
| Owner | Deputy Director (Technology) |
| Steward | Lead Developer (G7) |
| Parent document | Data & Digital Asset Governance Framework v1.3 |
| Dependency management | `renv` — run `renv::restore()` on first use |
| Testing | `testthat` — 70% line coverage required for Tier 1 functions |

---

## Upcoming features

See [UPCOMING_FEATURES.md](UPCOMING_FEATURES.md) for the full roadmap, including
HiFlows-UK zip ingestion, Rcpp PDM acceleration, baseflow separation, FEH
regional DDF parameters, and `reach.validate` integration.

---

## Known TODOs

- Replace `PdmParams` / `ReachHydroResult` stubs with full `S7` class
  definitions once S7 is on CRAN.
- FEH DDF regional parameter sets (currently UK national average only);
  add regional coefficients from FEH CD-ROM when available.
- Baseflow separation from observed series (Lyne-Hollick, Boughton-Eckhardt).
- Thiessen polygon weight computation for `areal_rainfall()`.
- L-moment ratio diagram plotting function using `lmrd_theoretical()`.
- Benchmark FEH Vol 3 GLO/GEV fits against WINFAP outputs on real AMAX series
  before operational use.

---

## References

Faulkner, D. (1999). *Rainfall Frequency Estimation*. FEH Vol. 2.
CEH Wallingford.

Hosking, J.R.M. (1990). L-moments: Analysis and estimation of distributions
using linear combinations of order statistics. *JRSS-B*, 52(1), 105–124.

Hosking, J.R.M. & Wallis, J.R. (1997). *Regional Frequency Analysis*.
Cambridge University Press.

Kjeldsen, T.R. (2007). *The revitalised FSR/FEH rainfall-runoff method*.
Flood Studies Report Technical Note No. 1. CEH Wallingford.

Moore, R.J. (2007). The PDM rainfall-runoff model. *Hydrology and Earth System
Sciences*, 11, 483–499. <https://doi.org/10.5194/hess-11-483-2007>

NERC (1975). *Flood Studies Report*. Vol. 1. Natural Environment Research
Council, London.

Robson, A. & Reed, D. (1999). *Statistical Procedures for Flood Frequency
Estimation*. FEH Vol. 3. CEH Wallingford.
