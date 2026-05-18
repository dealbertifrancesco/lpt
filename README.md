# lpt: Local Parallel Trends

R package for partial identification in continuous difference-in-differences
designs under Local Parallel Trends (LPT).

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
this yields identified sets whose width grows with the post-treatment
horizon $t$ (number of periods since treatment onset):

$$IS_{\partial ATT}(d, t;\,B) = [\lambda_t(d) - (t+1)B,\ \lambda_t(d) + (t+1)B]$$

$$IS_{ATT}(d, t;\,B) = [\Lambda_t(d) - (t+1)Bd,\ \Lambda_t(d) + (t+1)Bd]$$

The identified set widens linearly with the horizon because each
additional post-period accumulates one more unit of potential
counterfactual drift.

The package also computes a **time-aggregated** overall ATT across all
requested post-periods ($ATT^0$), with identified set:

$$IS_{ATT^0}(B) = \left[\bar\Lambda^{agg} - \frac{T+2}{2}B\bar{D}_+,\ \bar\Lambda^{agg} + \frac{T+2}{2}B\bar{D}_+\right]$$

The parameter $B$ governs how far from standard parallel trends the
selection process is allowed to deviate. Setting $B = 0$ recovers point
identification under standard parallel trends. $B$ can be supplied by the
user or **calibrated from pre-treatment periods**.

## Quick Start

```r
library(lpt)
data(sru)

# Estimate with calibrated B across all post-periods
fit <- lpt(sru, "commune", "year", "outcome", "dose",
           post_period = 0:5, pre_periods = -7:-1,
           B = "calibrate")
summary(fit)
```

### Plotting

```r
# Event study: ATT^o across all periods (pre + post)
plot(fit, type = "eventstudy")

# Dose-response slope IS (faceted by period)
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
# Sensitivity of the overall ATT^o  [default -- no dose needed]
plot(fit, type = "sensitivity")

# Sensitivity of ATT(d|d) at a specific dose -- d0 required
plot(fit, type = "sensitivity", estimand = "att", d0 = 0.3)

# Sensitivity of the dose-response slope at a specific dose
plot(fit, type = "sensitivity", estimand = "datt", d0 = 0.3)
```

The vertical dotted line marks $\hat{B}$ (calibrated value); the x-axis is
expressed as $B/\hat{B}$ when calibration was used.

## Key Functions

| Function | Description |
|----------|-------------|
| `lpt()` | Main estimation function |
| `calibrate_B()` | Calibrate B from pre-treatment periods |
| `plot.lpt()` | Visualize identified sets and sensitivity |
| `summary.lpt()` | Tabular summary at dose quartiles |

## Data

The package ships with the `sru` dataset: panel data on 991 French communes
affected by the 2000 SRU French Law. Periods are re-indexed so that
treatment onset corresponds to period 0 (pre-treatment: -7 to -1;
post-treatment: 0 to 5).

```r
data(sru)
?sru
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

## AI Use Acknowledgment

See [AI_USE.md](AI_USE.md) for details on AI-assisted development.
