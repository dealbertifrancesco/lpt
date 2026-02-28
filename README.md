# lpt: Local Parallel Trends

R package for partial identification in continuous difference-in-differences
designs under Local Parallel Trends (LPT).

> Dealberti (2026). *Local Parallel Trends for Continuous
> Difference-in-Differences.*

## Installation

```r
# Install from GitHub
remotes::install_github("dealbertifrancesco/lpt")
```

## Overview

In continuous DiD, the observable dose-slope decomposes as:

$$\lambda(d) = \frac{\partial ATT(d|d)}{\partial d} + \mu'(d)$$

where $\lambda(d)$ is identified from data but the **selection slope** $\mu'(d)$
is not. Under the **Local Parallel Trends** assumption $|\mu'(d)| \leq B$,
this yields identified sets for three estimands:

| Estimand | Identified set |
|----------|---------------|
| $\partial ATT(d\|d)/\partial d$ | $[\lambda(d) - B,\ \lambda(d) + B]$ |
| $ATT(d\|d)$ | $[\Lambda(d) - Bd,\ \Lambda(d) + Bd]$ |
| $ATT^o$ (overall) | $[ATT^o - B\bar{D}_+,\ ATT^o + B\bar{D}_+]$ |

where $\Lambda(d) = E[\Delta Y \mid D=d] - E[\Delta Y \mid D=0]$ and
$ATT^o_{bin} = E[\Delta Y \mid D>0] - E[\Delta Y \mid D=0]$.

The parameter $B$ governs how far from standard parallel trends the
selection process is allowed to deviate. Setting $B = 0$ recovers point
identification under standard parallel trends. $B$ can be supplied by the
user or **calibrated from pre-treatment periods**.

## Quick Start

```r
library(lpt)
data(sru)

# Estimate with calibrated B from pre-treatment periods
fit <- lpt(sru, "commune", "year", "outcome", "dose",
           post_period = 2019, pre_periods = 1993:1999,
           B = "calibrate")
summary(fit)
```

### Plotting identified sets

```r
# Dose-response slope IS
plot(fit, type = "datt")

# ATT level IS (requires untreated units)
plot(fit, type = "att")

# Pre-period selection slopes (calibration diagnostic)
plot(fit, type = "pretrends")
```

### Sensitivity analysis

The sensitivity plot shows how an identified set widens as $B$ grows.
The `estimand` argument selects which quantity to trace:

```r
# Sensitivity of the overall ATT^o  [default — no dose needed]
plot(fit, type = "sensitivity")

# Sensitivity of ATT(d|d) at a specific dose — d0 required
plot(fit, type = "sensitivity", estimand = "att", d0 = 0.3)

# Sensitivity of the dose-response slope at a specific dose — d0 required
plot(fit, type = "sensitivity", estimand = "datt", d0 = 0.3)
```

The vertical dotted line marks $\hat{B}$ (calibrated value); the x-axis is
expressed as $B/\hat{B}$ when calibration was used, making it easy to read
off "how much larger than the pre-trend bound does B need to be for the
conclusion to overturn?"

## Key Functions

| Function | Description |
|----------|-------------|
| `lpt()` | Main estimation function |
| `calibrate_B()` | Calibrate B from pre-treatment periods |
| `estimate_dose_slope()` | Estimate $\lambda(d)$ via penalized splines |
| `plot.lpt()` | Visualize identified sets and sensitivity |
| `summary.lpt()` | Tabular summary at dose quartiles |

## Data

The package ships with the `sru` dataset: panel data on 991 French communes
affected by the 2000 Solidarity and Urban Renewal (SRU) law, with annual
observations from 1993 to 2019.

```r
data(sru)
?sru
```

## Multiple Post-Periods

Pass a vector to `post_period` to estimate across several waves:

```r
fit_multi <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2005, 2010, 2015, 2019),
                 pre_periods = 1993:1999,
                 B = "calibrate")
plot(fit_multi, type = "datt")   # faceted by period
```

## B Calibration Details

When `B = "calibrate"`, the package estimates $\mu'(d)$ in each consecutive
pre-period pair and sets $\hat{B} = \max_s \sup_d |\hat{\mu}'_s(d)|$. This
uses the fact that in pre-treatment periods there is no treatment effect, so
the entire dose-slope reflects selection. The `"pretrends"` plot visualises
these slopes with $\pm\hat{B}$ bands.

## Dependencies

- **Required:** `mgcv`, `stats`
- **Suggested:** `ggplot2` (plotting), `testthat`, `knitr`, `rmarkdown`

## References

Dealberti, F. (2026). Local Parallel Trends for Continuous
Difference-in-Differences.

## AI Use Acknowledgment

See [AI_USE.md](AI_USE.md) for details on AI-assisted development.
