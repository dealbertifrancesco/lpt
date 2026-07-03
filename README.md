# lpt: Local Parallel Trends

R package for partial identification in continuous difference-in-differences
designs under local restrictions on counterfactual trends.

## Installation

```r
# Install from GitHub
remotes::install_github("dealbertifrancesco/lpt")
```

## Overview

In continuous DiD, parallel trends requires untreated counterfactual trends
to be identical across all dose groups — an assumption that is often
implausible and, even when it holds, does not identify causal responses.
`lpt` replaces the equality with **bounds** on the selection function
$\mu_t(d) = E[\Delta Y_t(0) \mid D = d]$, matching the strength of the
bound to the target parameter:

- **Level restriction (Assumption L).** A bound $M$ on trend deviations,
  $|\mu_t(d) - \mu_t(0)| \leq M$, delivers identified sets for treatment
  effect **levels**:

$$IS_{ATT}(d, t;\,M) = [\Lambda_t(d) - (t+1)M,\ \Lambda_t(d) + (t+1)M]$$

- **Slope restriction (Assumption S).** A separate bound $B$ on selection
  slopes, $|\mu'_t(d)| \leq B$ on the interior of the treated dose support,
  delivers identified sets for **dose-response derivatives**:

$$IS_{\partial ATT}(d, t;\,B) = [\lambda_t(d) - (t+1)B,\ \lambda_t(d) + (t+1)B]$$

Here $\Lambda_t(d) = E[Y_t - Y_{-1} \mid D=d] - E[Y_t - Y_{-1} \mid D=0]$ is
the horizon-$t$ DiD estimand and $\lambda_t(d)$ its dose derivative. The
identified sets widen linearly with the horizon because each additional
post-period accumulates one more unit of potential counterfactual drift.

Neither assumption implies the other: the dose support has a hole between 0
and the smallest positive dose, so a within-support Lipschitz condition
cannot control the deviation from the untreated group without extrapolating
across the gap. This is why levels get their own parameter $M$ instead of a
dose-proportional bound derived from $B$.

The package also computes the **time-aggregated** overall ATT across all
requested post-periods ($ATT^o$), with identified set

$$IS_{ATT^o}(M) = \left[\bar\Lambda^{agg} - \frac{T+2}{2}M,\ \bar\Lambda^{agg} + \frac{T+2}{2}M\right]$$

Setting $M = 0$ (resp. $B = 0$) recovers point identification under
standard parallel trends for the corresponding target. Both parameters can
be supplied by the user or **calibrated from pre-treatment periods**.

## Quick Start

```r
library(lpt)
data(sru)

# Estimate with calibrated M and B across all post-periods
fit <- lpt(sru, "commune", "year", "outcome", "dose",
           post_period = 0:5, pre_periods = -7:-1,
           M = "calibrate", B = "calibrate")
summary(fit)
```

### Plotting

```r
# Event study: ATT^o across all periods with M-based identified sets
plot(fit, type = "eventstudy")

# Dose-response slope IS (slope bound B; faceted by period)
plot(fit, type = "datt")

# ATT level IS (level bound M; requires untreated units)
plot(fit, type = "att")

# Pre-period diagnostics: trend deviations (M) and selection slopes (B)
plot(fit, type = "pretrends")
```

### Sensitivity analysis

The sensitivity plot shows how an identified set widens as the relevant
bound grows: level estimands (`"att"`, `"att_o"`) trace $M$, the slope
estimand (`"datt"`) traces $B$:

```r
# Sensitivity of the overall ATT^o in M  [default -- no dose needed]
plot(fit, type = "sensitivity")

# Sensitivity of ATT(d|d) in M at a specific dose -- d0 required
plot(fit, type = "sensitivity", estimand = "att", d0 = 0.3)

# Sensitivity of the dose-response slope in B at a specific dose
plot(fit, type = "sensitivity", estimand = "datt", d0 = 0.3)
```

The vertical dotted line marks the calibrated value ($\hat{M}$ or
$\hat{B}$); the x-axis is expressed as a ratio to it when calibration was
used.

### Alternative backend: contdid

The `method = "contdid"` option uses B-spline estimation from
[Callaway, Goodman-Bacon & Sant'Anna (2024)](https://doi.org/10.3386/w32117)
instead of GAM-based penalized splines:

```r
fit_cd <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 0:5, pre_periods = -7:-1,
              M = "calibrate", B = "calibrate", method = "contdid",
              contdid_args = list(num_knots = 1, degree = 3))
summary(fit_cd)
plot(fit_cd, type = "eventstudy")
```

All plot and summary methods work identically. The `contdid` package must be
installed separately (`remotes::install_github("bcallaway11/contdid")`).

### Alternative backend: npiv

The `method = "npiv"` option uses nonparametric B-spline sieve regression from
[Chen, Christensen & Kankanala (2024)](https://doi.org/10.1093/restud/rdae025):

```r
fit_np <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 0:5, pre_periods = -7:-1,
              M = "calibrate", B = "calibrate", method = "npiv",
              npiv_args = list(J.x.segments = 2))
summary(fit_np)
plot(fit_np, type = "datt")
```

Uses the same long-difference approach as GAM but with B-spline
sieves and data-driven tuning parameter selection. Default `J.x.segments = 2`
(5-dimensional basis) matches the GAM default of `k = 5`.

### Alternative backend: kernel

The `method = "kernel"` option uses local polynomial kernel regression from
[Hayfield & Racine (2008)](https://doi.org/10.18637/jss.v027.i05):

```r
fit_kr <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 0:5, pre_periods = -7:-1,
              M = "calibrate", B = "calibrate", method = "kernel",
              kernel_args = list(bw = "cv.ls", regtype = "ll"))
summary(fit_kr)
plot(fit_kr, type = "datt")
```

Local linear regression (`regtype = "ll"`, the default) automatically corrects
for boundary bias (Fan & Gijbels, 1996). Bandwidth can be selected via
least-squares cross-validation (`bw = "cv.ls"`, default) or set manually
(`bw = 0.1`). Kernel choices: `"gaussian"` (default), `"epanechnikov"`,
`"uniform"`.

## Key Functions

| Function | Description |
|----------|-------------|
| `lpt()` | Main estimation (`method = "gam"`, `"contdid"`, `"npiv"`, or `"kernel"`) |
| `calibrate_bounds()` | Calibrate M and B from pre-treatment periods |
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

## Calibration Details

Under no anticipation, pre-period first differences identify the selection
function directly: $E[\Delta Y_s \mid D=d] = \mu_s(d)$ for $s < 0$. Under
**stable selection** (post-treatment selection is no worse than the worst
pre-treatment selection), the package sets

$$\hat{M} = \max_{s<0}\ \sup_d |\hat{\mu}_s(d) - \hat{\mu}_s(0)|, \qquad
\hat{B} = \max_{s<0}\ \sup_d |\hat{\mu}'_s(d)|$$

where $\hat{\mu}_s$ is a penalized-spline fit of $\Delta Y_s$ on dose using
**treated units only** (no smooth fit bridges the support gap) and
$\hat{\mu}_s(0)$ is the untreated sample mean. This replaces a pass/fail
pre-trend test with an explicit degree of agnosticism, in the spirit of
Rambachan & Roth (2023), transposed from the time axis to the dose axis.
Note the practical asymmetry: $\hat{M}$ needs only conditional means, while
$\hat{B}$ needs derivatives — a harder, noisier estimation problem — making
the level route the more robust default when only ATTs are of interest.
The `"pretrends"` plot visualises both diagnostics with $\pm\hat{M}$ and
$\pm\hat{B}$ bands.

## Dependencies

- **Required:** `mgcv`, `stats`
- **Suggested:** `contdid` (alternative backend), `np` (kernel backend), `npiv` (sieve backend), `ggplot2` (plotting), `testthat`, `knitr`, `rmarkdown`

## References

- Dealberti F (2026). *Local Parallel Trends: Relaxing Parallel Trends in Continuous DiD*.
- Callaway B, Goodman-Bacon A, Sant'Anna PHC (2024). *Difference-in-differences with a continuous treatment*. NBER Working Paper.
- Rambachan A, Roth J (2023). *A More Credible Approach to Parallel Trends*. Review of Economic Studies, 90(5), 2555–2591.
- Chen X, Christensen T, Kankanala S (2024). *Adaptive Estimation and Uniform Confidence Bands for Nonparametric Structural Functions and Elasticities*. Review of Economic Studies, 91(6), 3337–3369.
- Hayfield T, Racine JS (2008). *Nonparametric Econometrics: The np Package*. Journal of Statistical Software, 27(5).
- Fan J, Gijbels I (1996). *Local Polynomial Modelling and Its Applications*. Chapman and Hall/CRC.
- Wood SN (2017). *Generalized Additive Models: An Introduction with R* (2nd ed.). Chapman and Hall/CRC.

## AI Use Acknowledgment

See [AI_USE.md](AI_USE.md) for details on AI-assisted development.
