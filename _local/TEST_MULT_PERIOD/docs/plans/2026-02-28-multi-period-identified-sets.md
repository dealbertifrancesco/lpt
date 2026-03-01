# Multi-Period Identified Sets Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the `lpt` R package to support proper multi-period identified sets with `t`-stacked LPT bounds and optional time-direction restrictions (RM-Time, SD-Time), analogous to Rambachan & Roth (2023).

**Architecture:** All work is done in `_local/TEST_MULT_PERIOD/` as standalone R files that can be sourced. Modified copies of `R/*.R` live in `_local/TEST_MULT_PERIOD/R/`. Tests live in `_local/TEST_MULT_PERIOD/tests/`. A `source_all.R` helper bootstraps the environment for interactive use.

**Tech Stack:** R, mgcv (already a dependency), testthat (for tests)

---

## Background & Key Concepts

### The `t` multiplier (most important change)

The current code computes identified sets as:
- `ATT(d|d)`: `[Lambda(d) - B*d, Lambda(d) + B*d]`
- `dATT(d)`: `[lambda(d) - B, lambda(d) + B]`
- `ATT^o`: `[att_o_bin - B*D_bar, att_o_bin + B*D_bar]`

The **correct** multi-period formula applies a `t` multiplier (# of time increments from ref to post):
- `ATT_t(d|d)`: `[Lambda(d,t) - t*B*d, Lambda(d,t) + t*B*d]`
- `dATT_t(d)`: `[lambda(d,t) - t*B, lambda(d,t) + t*B]`
- `ATT^o_t`: `[Lambda_bar_t - t*B*D_bar+, Lambda_bar_t + t*B*D_bar+]`

### How `t` is computed

For integer time indices (simulate_lpt): `t_val = post_period_index - ref_period_index`

For calendar years (SRU data): `t_val = count of distinct times in data with ref_period < time <= post_period`

Example: ref_period=1999, post_period=2005, annual data → t=6.

### Three restriction classes

**LPT alone (dose-direction, parameter B_d):**
- Bounds above, with `t` multiplier.

**RM-Time (time-direction, relative magnitudes, parameter M_bar):**
- `delta_star_pre(d) = max over pre-period pairs of |delta_tilde_s(d)|`
- `delta_tilde_s(d) = E[DeltaY_s | D=d] - E[DeltaY_s | D=0]` (estimated via GAM: `predict(gam, dose=d) - predict(gam, dose=0)`)
- Half-width for ATT: `t * M_bar * delta_star_pre(d)`
- Half-width for ATT^o: `t * M_bar * delta_bar_star_pre` where `delta_bar_star_pre = mean(delta_star_pre(D) | D > 0)` (weighted by empirical distribution of D)
- Slope NOT bounded by RM-Time alone (needs LPT for that)

**SD-Time (time-direction, smoothness, parameter B_t):**
- Uses `delta_tilde_0(d)` = last pre-period increment bias (identified from data)
- Center adjustment: `Lambda(d,t) - t * delta_tilde_0(d)`
- Half-width: `B_t * t*(t+1)/2`
- ATT^o center: `Lambda_bar_t - t * delta_tilde_0_bar` where `delta_tilde_0_bar = mean(delta_tilde_0(D) | D > 0)`
- Slope NOT bounded by SD-Time alone

**Joint (LPT + time):** Intersection = `[max(lower_lpt, lower_time), min(upper_lpt, upper_time)]`

### Output additions

Add `binding_lower` and `binding_upper` columns to `att` data frame indicating which restriction is tighter at each `(d, t)` point: `"lpt"`, `"time"`, or `"equal"`.

---

## Working Directory Convention

All paths in this plan are relative to the **package root**: `C:/Users/dealb/Research/lpt/`

Modified files live in: `_local/TEST_MULT_PERIOD/R/`
Test files live in: `_local/TEST_MULT_PERIOD/tests/`

---

## Task 1: Setup — Directory Bootstrap and Source Helper

**Files:**
- Create: `_local/TEST_MULT_PERIOD/source_all.R`
- Create: `_local/TEST_MULT_PERIOD/tests/helpers.R`

**Step 1: Create `source_all.R`**

```r
# _local/TEST_MULT_PERIOD/source_all.R
# Run this from the package root to load all modified functions.
# Usage: source("_local/TEST_MULT_PERIOD/source_all.R")

if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv required")
if (!requireNamespace("ggplot2", quietly = TRUE)) message("ggplot2 not found - plot tests will be skipped")

# Load unmodified package functions from package R/ directory
source("R/estimate_dose_slope.R")

# Load modified functions from test directory (order matters: dependencies first)
source("_local/TEST_MULT_PERIOD/R/calibrate_B.R")
source("_local/TEST_MULT_PERIOD/R/att_bounds.R")
source("_local/TEST_MULT_PERIOD/R/simulate.R")
source("_local/TEST_MULT_PERIOD/R/lpt.R")
source("_local/TEST_MULT_PERIOD/R/summary.R")
source("_local/TEST_MULT_PERIOD/R/plot.R")

# Load SRU data
load("data/sru.rda")

message("TEST_MULT_PERIOD environment loaded.")
```

**Step 2: Create `tests/helpers.R`**

```r
# _local/TEST_MULT_PERIOD/tests/helpers.R
# Common setup for all test files. Source this at top of each test file.

# Ensure we're run from package root
if (!file.exists("data/sru.rda")) {
  stop("Run tests from the package root directory (lpt/).")
}

source("_local/TEST_MULT_PERIOD/source_all.R")
library(testthat)
```

**Step 3: Verify the helpers source without error**

From the package root, run:
```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")
```
Expected: message "TEST_MULT_PERIOD environment loaded." — but also errors because the modified R files don't exist yet. That's fine for now.

**Step 4: Copy original R files as starting point**

Copy each of the 6 files that will be modified:
```bash
cp R/calibrate_B.R _local/TEST_MULT_PERIOD/R/calibrate_B.R
cp R/att_bounds.R  _local/TEST_MULT_PERIOD/R/att_bounds.R
cp R/simulate.R    _local/TEST_MULT_PERIOD/R/simulate.R
cp R/lpt.R         _local/TEST_MULT_PERIOD/R/lpt.R
cp R/summary.R     _local/TEST_MULT_PERIOD/R/summary.R
cp R/plot.R        _local/TEST_MULT_PERIOD/R/plot.R
```

**Step 5: Verify source_all.R works with copies**

```r
source("_local/TEST_MULT_PERIOD/source_all.R")
# Should print "TEST_MULT_PERIOD environment loaded." with no errors
```

**Step 6: Commit**
```bash
git add _local/TEST_MULT_PERIOD/
git commit -m "feat(multiperiod): scaffold TEST_MULT_PERIOD directory and source helper"
```

---

## Task 2: Compute `t` per Post-Period in `lpt.R`

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/lpt.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-t-computation.R`

**Background:** `t` is the number of distinct time periods in the data strictly after `ref_period` and up to (inclusive) `pp`. For the SRU data with ref=1999, post=2005: `t = sum(all_times > 1999 & all_times <= 2005)`. For simulate_lpt with integer indices and ref=n_pre, post=n_pre+2: `t = 2`.

**Step 1: Write failing test**

Create `_local/TEST_MULT_PERIOD/tests/test-t-computation.R`:

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("t computation is correct for annual calendar-year data", {
  # Use simulate_lpt as proxy: n_pre=3, n_post=5 → time indices 1..8
  # ref_period = 3 (last pre-period). Post-periods: 4,5,6,7,8.
  sim <- simulate_lpt(n = 200, n_pre = 3, n_post = 5, dgp = "linear_selection", seed = 1)
  all_times <- sort(unique(sim$time))
  ref_period <- 3L

  # t for post_period=4 should be 1; for post_period=6 should be 3
  t4 <- sum(all_times > ref_period & all_times <= 4L)
  t6 <- sum(all_times > ref_period & all_times <= 6L)

  expect_equal(t4, 1L)
  expect_equal(t6, 3L)
})

test_that("t computation is stored in lpt output", {
  sim <- simulate_lpt(n = 300, n_pre = 2, n_post = 3, dgp = "linear_selection", seed = 42)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L),
             pre_periods = 1:2,
             B = attr(sim, "true_B"))

  # t values should be stored in fit$t_values
  expect_true(!is.null(fit$t_values))
  # ref_period = 2, so t for post 3 = 1, post 4 = 2, post 5 = 3
  expect_equal(fit$t_values[["3"]], 1L)
  expect_equal(fit$t_values[["4"]], 2L)
  expect_equal(fit$t_values[["5"]], 3L)
})
```

**Step 2: Run test to confirm it fails**

```r
source("_local/TEST_MULT_PERIOD/tests/test-t-computation.R")
```
Expected: second test fails — `fit$t_values` is NULL.

**Step 3: Add `t_values` computation to `_local/TEST_MULT_PERIOD/R/lpt.R`**

In `lpt()`, after setting `ref_period`, add a block that computes `t_values` for each post-period. Find the line:
```r
  if (length(pre_period_set) < 1) {
```
Insert before it:

```r
  # Compute t (# time increments from ref to each post-period) from data
  all_times_sorted <- sort(unique(data[[time_col]]))
  t_values <- vapply(post_periods, function(pp) {
    sum(all_times_sorted > ref_period & all_times_sorted <= pp)
  }, integer(1L))
  names(t_values) <- as.character(post_periods)
```

Also add `t_values` to the returned `structure()` at the bottom:
```r
  structure(
    list(
      datt     = datt,
      att      = att,
      att_o    = att_o,
      B_hat    = B_hat,
      B_values = B_values,
      calibration = calibration_result,
      slopes   = slopes,
      call     = the_call,
      n        = nrow(ref_data),
      has_untreated = has_untreated,
      post_periods  = post_periods,
      ref_period    = ref_period,
      t_values      = t_values,          # NEW
      specifications = list(
        k = k, spline_bs = spline_bs, alpha = alpha,
        post_periods = post_periods, ref_period = ref_period
      )
    ),
    class = "lpt"
  )
```

**Step 4: Run test to confirm it passes**

```r
source("_local/TEST_MULT_PERIOD/tests/test-t-computation.R")
```
Expected: both tests PASS.

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/lpt.R _local/TEST_MULT_PERIOD/tests/test-t-computation.R
git commit -m "feat(multiperiod): compute and store t-values per post-period"
```

---

## Task 3: Apply `t` Multiplier to ATT Bounds

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/att_bounds.R`
- Modify: `_local/TEST_MULT_PERIOD/R/lpt.R` (pass `t` to `compute_att_bounds`)
- Create: `_local/TEST_MULT_PERIOD/tests/test-t-multiplier.R`

**Step 1: Write failing test**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("ATT bounds half-width scales with t (LPT alone)", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 4, dgp = "linear_selection",
                      beta = c(1, 0), rho = 0.3, seed = 1)
  B_true <- attr(sim, "true_B")  # = 0.3

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L, 6L),
             pre_periods = 1:2,
             B = B_true)

  # For a dose d, ATT half-width at post-period t should be t * B * d
  # Check at median dose
  ep <- fit$att$d[fit$att$period == 3L]
  d_med <- ep[which.min(abs(ep - stats::median(ep)))]

  for (t_idx in 1:4) {
    pp_t <- 2L + t_idx  # post-periods are 3,4,5,6
    att_pp <- fit$att[fit$att$period == pp_t & fit$att$B == B_true, ]
    row <- att_pp[which.min(abs(att_pp$d - d_med)), ]
    expected_hw <- t_idx * B_true * row$d
    actual_hw   <- (row$att_upper - row$att_lower) / 2
    expect_equal(actual_hw, expected_hw, tolerance = 1e-10,
                 label = paste("half-width at t =", t_idx))
  }
})

test_that("ATT bounds are wider for later periods", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 3, dgp = "linear_selection", seed = 2)
  B_true <- attr(sim, "true_B")
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L), pre_periods = 1:2, B = B_true)

  ep <- fit$att$d[fit$att$period == 3L][1]

  hw_t1 <- with(fit$att[fit$att$period == 3L & fit$att$B == B_true, ],
                (att_upper - att_lower)[which.min(abs(d - ep))])
  hw_t2 <- with(fit$att[fit$att$period == 4L & fit$att$B == B_true, ],
                (att_upper - att_lower)[which.min(abs(d - ep))])
  hw_t3 <- with(fit$att[fit$att$period == 5L & fit$att$B == B_true, ],
                (att_upper - att_lower)[which.min(abs(d - ep))])

  expect_lt(hw_t1, hw_t2)
  expect_lt(hw_t2, hw_t3)
})
```

**Step 2: Run to confirm tests fail**

```r
source("_local/TEST_MULT_PERIOD/tests/test-t-multiplier.R")
```
Expected: first test fails — half-width is currently `B * d` not `t * B * d`.

**Step 3: Update `compute_att_bounds()` to accept and apply `t`**

In `_local/TEST_MULT_PERIOD/R/att_bounds.R`, change the signature and formula:

```r
compute_att_bounds <- function(slope_result, B_values, dose, period = NA, t = 1L,
                               time_restriction = "none",
                               M_bar = NULL, delta_star_pre_d = NULL,
                               B_t = NULL, delta_tilde_0_d = NULL) {
  ep       <- slope_result$eval_points
  gam_fit  <- slope_result$gam_fit

  # Lambda(d) = E[DeltaY | D=d] - E[DeltaY | D=0]
  cm_at_ep   <- as.numeric(stats::predict(gam_fit, newdata = data.frame(dose = ep)))
  cm_at_zero <- as.numeric(stats::predict(gam_fit, newdata = data.frame(dose = 0)))
  Lambda_d   <- cm_at_ep - cm_at_zero

  att_list <- lapply(B_values, function(b) {
    # --- LPT bounds (with t multiplier) ---
    lpt_hw     <- t * b * ep  # half-width under LPT
    lpt_lower  <- Lambda_d - lpt_hw
    lpt_upper  <- Lambda_d + lpt_hw

    # --- Time-restriction bounds ---
    if (time_restriction == "rm") {
      if (is.null(delta_star_pre_d))
        stop("delta_star_pre_d required for time_restriction = 'rm'.")
      rm_hw    <- t * M_bar * delta_star_pre_d
      rm_lower <- Lambda_d - rm_hw
      rm_upper <- Lambda_d + rm_hw
    } else if (time_restriction == "sd") {
      if (is.null(delta_tilde_0_d) || is.null(B_t))
        stop("delta_tilde_0_d and B_t required for time_restriction = 'sd'.")
      sd_hw     <- B_t * t * (t + 1) / 2
      sd_center <- Lambda_d - t * delta_tilde_0_d
      sd_lower  <- sd_center - sd_hw
      sd_upper  <- sd_center + sd_hw
    }

    # --- Determine final bounds and which restriction binds ---
    if (time_restriction == "none") {
      final_lower <- lpt_lower
      final_upper <- lpt_upper
      bind_lower  <- rep("lpt", length(ep))
      bind_upper  <- rep("lpt", length(ep))
    } else if (time_restriction == "rm") {
      if (b > 0) {
        final_lower <- pmax(lpt_lower, rm_lower)
        final_upper <- pmin(lpt_upper, rm_upper)
      } else {
        final_lower <- rm_lower
        final_upper <- rm_upper
      }
      bind_lower <- ifelse(b > 0 & lpt_lower >= rm_lower, "lpt",
                    ifelse(b > 0 & lpt_lower < rm_lower, "rm", "rm"))
      bind_upper <- ifelse(b > 0 & lpt_upper <= rm_upper, "lpt",
                    ifelse(b > 0 & lpt_upper > rm_upper, "rm", "rm"))
    } else {  # "sd"
      if (b > 0) {
        final_lower <- pmax(lpt_lower, sd_lower)
        final_upper <- pmin(lpt_upper, sd_upper)
      } else {
        final_lower <- sd_lower
        final_upper <- sd_upper
      }
      bind_lower <- ifelse(b > 0 & lpt_lower >= sd_lower, "lpt", "sd")
      bind_upper <- ifelse(b > 0 & lpt_upper <= sd_upper, "lpt", "sd")
    }

    # Check for empty intersection
    empty_idx <- final_lower > final_upper
    if (any(empty_idx)) {
      warning(sprintf("%d dose point(s) have empty intersection of LPT and time bounds.",
                      sum(empty_idx)))
      final_lower[empty_idx] <- NA_real_
      final_upper[empty_idx] <- NA_real_
    }

    data.frame(
      period        = period,
      d             = ep,
      Lambda_d      = Lambda_d,
      B             = b,
      t             = t,
      att_lower     = final_lower,
      att_upper     = final_upper,
      binding_lower = bind_lower,
      binding_upper = bind_upper
    )
  })
  do.call(rbind, att_list)
}
```

**Step 4: Update the call site in `lpt.R`**

In the loop `for (pp in post_periods)`, the line:
```r
att_pp <- compute_att_bounds(sr, B_values, dose_vec, period = pp)
```
Becomes:
```r
t_val  <- t_values[[as.character(pp)]]
att_pp <- compute_att_bounds(sr, B_values, dose_vec, period = pp, t = t_val)
```

**Step 5: Run test to confirm it passes**

```r
source("_local/TEST_MULT_PERIOD/tests/test-t-multiplier.R")
```
Expected: both tests PASS.

**Step 6: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/att_bounds.R _local/TEST_MULT_PERIOD/R/lpt.R \
        _local/TEST_MULT_PERIOD/tests/test-t-multiplier.R
git commit -m "feat(multiperiod): apply t multiplier to ATT bounds in compute_att_bounds"
```

---

## Task 4: Apply `t` Multiplier to dATT and ATT^o Bounds in `lpt.R`

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/lpt.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-t-datt-atto.R`

**Step 1: Write failing tests**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("dATT half-width scales with t", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 3, dgp = "linear_selection",
                      beta = c(1, 0), rho = 0.3, seed = 3)
  B_true <- attr(sim, "true_B")  # 0.3

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L), pre_periods = 1:2, B = B_true)

  # dATT half-width at each t should be t * B
  for (t_idx in 1:3) {
    pp_t  <- 2L + t_idx
    datt_pp <- fit$datt[fit$datt$period == pp_t & fit$datt$B == B_true, ]
    hw <- (datt_pp$datt_upper - datt_pp$datt_lower) / 2
    expect_equal(unique(hw), t_idx * B_true, tolerance = 1e-10,
                 label = paste("dATT half-width at t =", t_idx))
  }
})

test_that("ATT^o half-width scales with t", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 3, dgp = "linear_selection",
                      beta = c(1, 0), rho = 0.3, seed = 4)
  B_true <- attr(sim, "true_B")  # 0.3

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L), pre_periods = 1:2, B = B_true)

  D_bar <- fit$att_o$D_bar[fit$att_o$period == 3L & fit$att_o$B == B_true]

  for (t_idx in 1:3) {
    pp_t  <- 2L + t_idx
    atto_row <- fit$att_o[fit$att_o$period == pp_t & fit$att_o$B == B_true, ]
    hw <- (atto_row$att_o_upper - atto_row$att_o_lower) / 2
    expect_equal(hw, t_idx * B_true * atto_row$D_bar, tolerance = 1e-10,
                 label = paste("ATT^o half-width at t =", t_idx))
  }
})
```

**Step 2: Run to confirm tests fail**

```r
source("_local/TEST_MULT_PERIOD/tests/test-t-datt-atto.R")
```

**Step 3: Update dATT bounds in `lpt.R`**

In the loop over `B_values` inside `for (pp in post_periods)`, find:
```r
    for (b in B_values) {
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period = pp,
        d = ep,
        lambda_d = lambda_d,
        se_lambda = se_lambda,
        B = b,
        datt_lower = lambda_d - b,
        datt_upper = lambda_d + b,
        ci_lower = lambda_d - b - z_alpha * se_lambda,
        ci_upper = lambda_d + b + z_alpha * se_lambda
      )
    }
```

Replace with:
```r
    t_val <- t_values[[pp_char]]  # already computed above

    for (b in B_values) {
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period    = pp,
        d         = ep,
        lambda_d  = lambda_d,
        se_lambda = se_lambda,
        B         = b,
        t         = t_val,
        datt_lower = lambda_d - t_val * b,
        datt_upper = lambda_d + t_val * b,
        ci_lower   = lambda_d - t_val * b - z_alpha * se_lambda,
        ci_upper   = lambda_d + t_val * b + z_alpha * se_lambda
      )
    }
```

**Step 4: Update ATT^o bounds in `lpt.R`**

Find:
```r
      for (b in B_values) {
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period = pp,
          att_o_bin = att_o_bin,
          D_bar = D_bar,
          B = b,
          att_o_lower = att_o_bin - b * D_bar,
          att_o_upper = att_o_bin + b * D_bar
        )
      }
```

Replace with:
```r
      for (b in B_values) {
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period    = pp,
          att_o_bin = att_o_bin,
          D_bar     = D_bar,
          B         = b,
          t         = t_val,
          att_o_lower = att_o_bin - t_val * b * D_bar,
          att_o_upper = att_o_bin + t_val * b * D_bar
        )
      }
```

**Step 5: Run tests to confirm they pass**

```r
source("_local/TEST_MULT_PERIOD/tests/test-t-datt-atto.R")
```
Expected: both PASS.

**Step 6: Verify backward compatibility (t=1 single period)**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")
data(sru)
fit_old <- lpt(sru, "commune", "year", "outcome", "dose",
               post_period = 2019, pre_periods = 1993:1999, B = 0.1)
# t for 2019 with ref=1999 in annual data = 20
# dATT half-width should be t * B = 20 * 0.1 = 2.0
unique((fit_old$datt$datt_upper - fit_old$datt$datt_lower) / 2)
# Expect: 2.0
```

**Step 7: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/lpt.R _local/TEST_MULT_PERIOD/tests/test-t-datt-atto.R
git commit -m "feat(multiperiod): apply t multiplier to dATT and ATT^o bounds"
```

---

## Task 5: Extend `calibrate_B` to Return `delta_tilde_s(d)` Levels

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/calibrate_B.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-delta-tilde.R`

**Background:** `delta_tilde_s(d) = E[DeltaY_s | D=d] - E[DeltaY_s | D=0]`. For the GAM fit of pre-period pair `(t0, t1)`, this is `predict(gam_fit, newdata=data.frame(dose=d)) - predict(gam_fit, newdata=data.frame(dose=0))`. This gives the level of differential trend at each evaluation point. We need this for RM-Time and SD-Time.

**Step 1: Write failing test**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("calibrate_B returns delta_tilde_s levels", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 10)

  cal <- calibrate_B(sim, "id", "time", "outcome", "dose",
                     pre_periods = 1:3)

  # Should have pre_slopes with delta_tilde_d column
  expect_true("delta_tilde_d" %in% names(cal$pre_slopes))

  # For linear_selection DGP, delta_tilde_s(d) = rho * d at any d
  # (because selection_fn(d,t) = rho*d*t, increment delta_tilde_s(d) = rho*d)
  # So delta_tilde_d should be approximately 0.3 * d at each eval point
  ep <- unique(cal$pre_slopes$d)
  pair1_dt <- cal$pre_slopes[cal$pre_slopes$period_pair == "1-2", ]
  expected <- 0.3 * pair1_dt$d
  expect_equal(pair1_dt$delta_tilde_d, expected, tolerance = 0.15)
})

test_that("calibrate_B returns delta_star_pre (max |delta_tilde_s| per d)", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 11)

  cal <- calibrate_B(sim, "id", "time", "outcome", "dose",
                     pre_periods = 1:3)

  expect_true(!is.null(cal$delta_star_pre))
  expect_true(all(cal$delta_star_pre >= 0))
  # delta_star_pre should equal max |delta_tilde_s(d)| over pairs
  # In this DGP all pairs have same delta_tilde ≈ 0.3*d, so delta_star_pre ≈ 0.3*d
  expect_equal(cal$delta_star_pre$delta_star_pre,
               abs(cal$delta_star_pre$d * 0.3),
               tolerance = 0.15)
})

test_that("calibrate_B returns delta_tilde_0 (last pre-period increment)", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 12)

  cal <- calibrate_B(sim, "id", "time", "outcome", "dose",
                     pre_periods = 1:3)

  # delta_tilde_0 is the delta_tilde for the last pre-period pair (2-3)
  expect_true(!is.null(cal$delta_tilde_0))
  expect_equal(cal$delta_tilde_0$delta_tilde_0, 0.3 * cal$delta_tilde_0$d,
               tolerance = 0.15)
})
```

**Step 2: Run to confirm they fail**

```r
source("_local/TEST_MULT_PERIOD/tests/test-delta-tilde.R")
```

**Step 3: Modify `calibrate_B.R`**

Inside the per-pair loop, after `slope_result <- estimate_dose_slope(...)`:

```r
    # Compute delta_tilde_s(d) = E[DeltaY_s | D=d] - E[DeltaY_s | D=0]
    gam_fit    <- slope_result$gam_fit
    ep         <- slope_result$eval_points
    pred_at_ep <- as.numeric(stats::predict(gam_fit, newdata = data.frame(dose = ep)))
    pred_at_0  <- as.numeric(stats::predict(gam_fit, newdata = data.frame(dose = 0)))
    delta_tilde_s <- pred_at_ep - pred_at_0
```

Add `delta_tilde_d = delta_tilde_s` to the `slopes_df` data frame:
```r
    slopes_df <- data.frame(
      period_pair   = pair_label,
      d             = slope_result$eval_points,
      mu_prime_d    = slope_result$lambda_d,
      se            = slope_result$se_lambda,
      delta_tilde_d = delta_tilde_s          # NEW
    )
```

After the loop, compute derived quantities. Find `pre_slopes <- do.call(rbind, all_slopes)` and add below it:

```r
  # delta_star_pre(d) = max over all pre-period pairs of |delta_tilde_s(d)|
  # Computed per eval point
  ep_vals        <- unique(pre_slopes$d)
  delta_star_vals <- vapply(ep_vals, function(dv) {
    dt_at_d <- pre_slopes$delta_tilde_d[pre_slopes$d == dv]
    max(abs(dt_at_d))
  }, numeric(1L))
  delta_star_pre <- data.frame(d = ep_vals, delta_star_pre = delta_star_vals)

  # delta_tilde_0(d) = delta_tilde for the LAST pre-period pair
  last_pair      <- pair_labels[length(pair_labels)]
  dt0_df         <- pre_slopes[pre_slopes$period_pair == last_pair,
                                c("d", "delta_tilde_d")]
  delta_tilde_0  <- data.frame(d = dt0_df$d, delta_tilde_0 = dt0_df$delta_tilde_d)
```

Add to the return list:
```r
  list(
    B_hat          = max(sup_vals),
    pre_slopes     = pre_slopes,
    sup_by_period  = sup_vals,
    delta_star_pre = delta_star_pre,  # NEW
    delta_tilde_0  = delta_tilde_0    # NEW
  )
```

**Step 4: Run tests to confirm they pass**

```r
source("_local/TEST_MULT_PERIOD/tests/test-delta-tilde.R")
```

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/calibrate_B.R _local/TEST_MULT_PERIOD/tests/test-delta-tilde.R
git commit -m "feat(multiperiod): extend calibrate_B to return delta_tilde_s levels and delta_star_pre"
```

---

## Task 6: Add `time_restriction` Argument to `lpt()` — RM-Time

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/lpt.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-rm-time.R`

**Background (RM-Time):**
- For ATT: half-width = `t * M_bar * delta_star_pre(d)`
- For ATT^o: half-width = `t * M_bar * delta_bar_star_pre` where `delta_bar_star_pre = mean(delta_star_pre(D) | D > 0)` weighted by the empirical distribution of treated doses
- Slope (dATT) is NOT bounded by RM-Time alone
- When both `B_d > 0` and `time_restriction = "rm"`, take intersection (handled in `compute_att_bounds`)

**Step 1: Write failing tests**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("lpt accepts time_restriction = 'rm' argument without error", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 20)
  expect_no_error(
    lpt(sim, "id", "time", "outcome", "dose",
        post_period = c(4L, 5L), pre_periods = 1:3,
        B = attr(sim, "true_B"),
        time_restriction = "rm", M_bar = 1)
  )
})

test_that("RM-Time ATT half-width = t * M_bar * delta_star_pre(d)", {
  sim <- simulate_lpt(n = 1000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 21)
  B_true <- attr(sim, "true_B")

  fit_rm <- lpt(sim, "id", "time", "outcome", "dose",
                post_period = c(4L, 5L), pre_periods = 1:3,
                B = 0,  # no LPT, pure RM-Time
                time_restriction = "rm", M_bar = 1)

  # For linear_selection DGP: delta_star_pre(d) ≈ rho * d = 0.3 * d
  # RM-Time half-width at t=1 (post=4): ≈ 1 * 1 * 0.3 * d = 0.3 * d
  att_t1 <- fit_rm$att[fit_rm$att$period == 4L & fit_rm$att$B == 0, ]
  hw_t1  <- (att_t1$att_upper - att_t1$att_lower) / 2
  expect_equal(hw_t1, fit_rm$calibration$delta_star_pre$delta_star_pre,
               tolerance = 0.2, label = "RM half-width at t=1 ≈ delta_star_pre")

  # At t=2, half-width should be double t=1
  att_t2 <- fit_rm$att[fit_rm$att$period == 5L & fit_rm$att$B == 0, ]
  hw_t2  <- (att_t2$att_upper - att_t2$att_lower) / 2
  expect_equal(hw_t2 / hw_t1, rep(2, length(hw_t1)), tolerance = 0.01,
               label = "RM half-width doubles from t=1 to t=2")
})

test_that("RM-Time joint with LPT produces intersection (tighter bounds)", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 22)
  B_true <- attr(sim, "true_B")

  fit_lpt <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = 4L, pre_periods = 1:3, B = B_true)
  fit_rm  <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = 4L, pre_periods = 1:3, B = 0,
                 time_restriction = "rm", M_bar = 1)
  fit_joint <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3, B = B_true,
                   time_restriction = "rm", M_bar = 1)

  # Joint set should be at least as tight as LPT or RM alone at each point
  hw_lpt   <- (fit_lpt$att$att_upper - fit_lpt$att$att_lower) / 2
  hw_rm    <- (fit_rm$att$att_upper - fit_rm$att$att_lower) / 2
  hw_joint <- (fit_joint$att$att_upper - fit_joint$att$att_lower) / 2

  expect_true(all(hw_joint <= hw_lpt + 1e-10),
              label = "Joint no wider than LPT alone")
  expect_true(all(hw_joint <= hw_rm + 1e-10),
              label = "Joint no wider than RM alone")
})

test_that("RM-Time: true ATT is inside the identified set", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 23)
  B_true     <- attr(sim, "true_B")
  true_att   <- attr(sim, "true_att")

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = B_true, time_restriction = "rm", M_bar = 1)

  for (pp in c(4L, 5L)) {
    att_pp <- fit$att[fit$att$period == pp & fit$att$B == B_true, ]
    true_vals <- true_att(att_pp$d)
    expect_true(all(att_pp$att_lower <= true_vals + 0.5 &
                      true_vals - 0.5 <= att_pp$att_upper),
                label = paste("True ATT inside IS at period", pp))
  }
})
```

**Step 2: Run to confirm tests fail**

```r
source("_local/TEST_MULT_PERIOD/tests/test-rm-time.R")
```

**Step 3: Add `time_restriction`, `M_bar`, `B_t` arguments to `lpt()`**

Change function signature in `_local/TEST_MULT_PERIOD/R/lpt.R`:
```r
lpt <- function(data, id_col, time_col, outcome_col, dose_col,
                post_period, pre_periods = NULL,
                B = "calibrate", eval_points = NULL,
                k = 20, spline_bs = "cr", alpha = 0.05,
                time_restriction = c("none", "rm", "sd"),
                M_bar = 1, B_t = 0) {
  time_restriction <- match.arg(time_restriction)
```

Add validation for RM-Time params (after existing validation):
```r
  if (time_restriction == "rm") {
    if (!is.numeric(M_bar) || length(M_bar) != 1 || M_bar <= 0)
      stop("M_bar must be a positive numeric scalar.")
  }
  if (time_restriction == "sd") {
    if (!is.numeric(B_t) || length(B_t) != 1 || B_t < 0)
      stop("B_t must be a non-negative numeric scalar.")
  }
  if (time_restriction != "none" && length(pre_period_set) < 2) {
    stop("time_restriction requires at least 2 pre-treatment periods for calibration.")
  }
```

Inside the `B == "calibrate"` block and also for numeric B, after setting `B_hat`:
When `time_restriction != "none"`, we always need `calibrate_B` to run (for `delta_tilde_s`). Add:
```r
  # Always run calibrate_B when we need pre-period delta info
  if (time_restriction != "none" && is.null(calibration_result)) {
    if (length(pre_period_set) < 2)
      stop("time_restriction requires at least 2 pre-treatment periods.")
    first_slope <- slopes[[1]]
    calibration_result <- calibrate_B(
      data        = data,
      id_col      = id_col,
      time_col    = time_col,
      outcome_col = outcome_col,
      dose_col    = dose_col,
      pre_periods = pre_period_set,
      eval_points = first_slope$eval_points,
      k           = k,
      spline_bs   = spline_bs
    )
  }
```

Update the `compute_att_bounds` call to pass time-restriction parameters. Replace:
```r
t_val  <- t_values[[as.character(pp)]]
att_pp <- compute_att_bounds(sr, B_values, dose_vec, period = pp, t = t_val)
```
With:
```r
t_val <- t_values[[pp_char]]

# Prepare time-restriction inputs
delta_star_pre_d <- NULL
delta_tilde_0_d  <- NULL
if (!is.null(calibration_result) && time_restriction != "none") {
  # Interpolate calibration results to match sr$eval_points
  cal_ep <- calibration_result$delta_star_pre$d
  sr_ep  <- sr$eval_points
  if (time_restriction == "rm") {
    delta_star_pre_d <- stats::approx(
      cal_ep,
      calibration_result$delta_star_pre$delta_star_pre,
      xout = sr_ep, rule = 2
    )$y
  } else if (time_restriction == "sd") {
    delta_tilde_0_d <- stats::approx(
      calibration_result$delta_tilde_0$d,
      calibration_result$delta_tilde_0$delta_tilde_0,
      xout = sr_ep, rule = 2
    )$y
  }
}

att_pp <- compute_att_bounds(
  slope_result     = sr,
  B_values         = B_values,
  dose             = dose_vec,
  period           = pp,
  t                = t_val,
  time_restriction = time_restriction,
  M_bar            = if (time_restriction == "rm") M_bar else NULL,
  delta_star_pre_d = delta_star_pre_d,
  B_t              = if (time_restriction == "sd") B_t else NULL,
  delta_tilde_0_d  = delta_tilde_0_d
)
```

**Update ATT^o bounds for RM-Time** — find the ATT^o block and update:
```r
      # delta_bar_star_pre = E[delta_star_pre(D) | D > 0]
      delta_bar_star_pre <- NULL
      if (time_restriction == "rm" && !is.null(calibration_result)) {
        # Weighted mean of delta_star_pre over treated dose values
        treat_d      <- d_merged[treated_idx]
        dsp          <- calibration_result$delta_star_pre
        dsp_at_treat <- stats::approx(dsp$d, dsp$delta_star_pre,
                                      xout = treat_d, rule = 2)$y
        delta_bar_star_pre <- mean(dsp_at_treat)
      }

      # delta_tilde_0_bar = E[delta_tilde_0(D) | D > 0]
      delta_tilde_0_bar <- NULL
      if (time_restriction == "sd" && !is.null(calibration_result)) {
        treat_d       <- d_merged[treated_idx]
        dt0           <- calibration_result$delta_tilde_0
        dt0_at_treat  <- stats::approx(dt0$d, dt0$delta_tilde_0,
                                       xout = treat_d, rule = 2)$y
        delta_tilde_0_bar <- mean(dt0_at_treat)
      }

      for (b in B_values) {
        # LPT bounds
        lpt_atto_hw    <- t_val * b * D_bar
        lpt_atto_lower <- att_o_bin - lpt_atto_hw
        lpt_atto_upper <- att_o_bin + lpt_atto_hw

        # Time-restriction bounds
        if (time_restriction == "rm" && !is.null(delta_bar_star_pre)) {
          rm_atto_hw    <- t_val * M_bar * delta_bar_star_pre
          rm_atto_lower <- att_o_bin - rm_atto_hw
          rm_atto_upper <- att_o_bin + rm_atto_hw
          if (b > 0) {
            final_lower <- max(lpt_atto_lower, rm_atto_lower)
            final_upper <- min(lpt_atto_upper, rm_atto_upper)
          } else {
            final_lower <- rm_atto_lower
            final_upper <- rm_atto_upper
          }
        } else if (time_restriction == "sd" && !is.null(delta_tilde_0_bar)) {
          sd_atto_hw    <- B_t * t_val * (t_val + 1) / 2
          sd_atto_ctr   <- att_o_bin - t_val * delta_tilde_0_bar
          sd_atto_lower <- sd_atto_ctr - sd_atto_hw
          sd_atto_upper <- sd_atto_ctr + sd_atto_hw
          if (b > 0) {
            final_lower <- max(lpt_atto_lower, sd_atto_lower)
            final_upper <- min(lpt_atto_upper, sd_atto_upper)
          } else {
            final_lower <- sd_atto_lower
            final_upper <- sd_atto_upper
          }
        } else {
          final_lower <- lpt_atto_lower
          final_upper <- lpt_atto_upper
        }

        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period    = pp,
          att_o_bin = att_o_bin,
          D_bar     = D_bar,
          B         = b,
          t         = t_val,
          att_o_lower = final_lower,
          att_o_upper = final_upper
        )
      }
```

Add `time_restriction`, `M_bar`, `B_t` to `specifications`:
```r
      specifications = list(
        k = k, spline_bs = spline_bs, alpha = alpha,
        post_periods = post_periods, ref_period = ref_period,
        time_restriction = time_restriction,
        M_bar = M_bar, B_t = B_t
      )
```

**Step 4: Run tests**

```r
source("_local/TEST_MULT_PERIOD/tests/test-rm-time.R")
```
Expected: all PASS.

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/lpt.R _local/TEST_MULT_PERIOD/tests/test-rm-time.R
git commit -m "feat(multiperiod): add time_restriction='rm' (RM-Time) to lpt()"
```

---

## Task 7: SD-Time Restriction

**Files:**
- `_local/TEST_MULT_PERIOD/R/lpt.R` and `R/att_bounds.R` (already updated in Task 6/3)
- Create: `_local/TEST_MULT_PERIOD/tests/test-sd-time.R`

**Background (SD-Time):**
- Center: `Lambda(d,t) - t * delta_tilde_0(d)`
- Half-width: `B_t * t * (t+1) / 2`
- When `B_t = 0`: point identification (no uncertainty from time dimension)

**Step 1: Write failing tests**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("SD-Time with B_t=0 gives point identification (linear_selection DGP)", {
  # For linear_selection: delta_tilde_s(d) = rho * d for all s (constant in time)
  # delta_tilde_0(d) = rho * d, so center = Lambda(d,t) - t * rho * d
  # With true ATT(d) = tau(d), Lambda(d,t) = tau(d) + delta_t(d) = tau(d) + t*rho*d
  # Center = tau(d) + t*rho*d - t*rho*d = tau(d) ← point identification
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3,
                      beta = c(1, 0), seed = 30)
  true_att <- attr(sim, "true_att")

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = 0, time_restriction = "sd", B_t = 0)

  for (pp in c(4L, 5L)) {
    att_pp <- fit$att[fit$att$period == pp & fit$att$B == 0, ]
    # Half-width should be zero
    hw <- (att_pp$att_upper - att_pp$att_lower) / 2
    expect_equal(hw, rep(0, length(hw)), tolerance = 1e-10)
    # Center should equal tau(d) up to estimation error
    expect_equal(att_pp$att_lower, true_att(att_pp$d), tolerance = 0.3)
  }
})

test_that("SD-Time half-width = B_t * t*(t+1)/2", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 31)
  B_t_val <- 0.1

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = 0, time_restriction = "sd", B_t = B_t_val)

  # t=1: hw = 0.1 * 1 * 2 / 2 = 0.1
  att_t1 <- fit$att[fit$att$period == 4L & fit$att$B == 0, ]
  hw_t1  <- (att_t1$att_upper - att_t1$att_lower) / 2
  expect_equal(unique(hw_t1), B_t_val * 1 * 2 / 2, tolerance = 1e-10)

  # t=2: hw = 0.1 * 2 * 3 / 2 = 0.3
  att_t2 <- fit$att[fit$att$period == 5L & fit$att$B == 0, ]
  hw_t2  <- (att_t2$att_upper - att_t2$att_lower) / 2
  expect_equal(unique(hw_t2), B_t_val * 2 * 3 / 2, tolerance = 1e-10)
})

test_that("SD-Time: true ATT inside identified set", {
  sim <- simulate_lpt(n = 1500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 32)
  true_att <- attr(sim, "true_att")
  B_t_val  <- 0.5

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = 0, time_restriction = "sd", B_t = B_t_val)

  for (pp in c(4L, 5L)) {
    att_pp <- fit$att[fit$att$period == pp & fit$att$B == 0, ]
    true_vals <- true_att(att_pp$d)
    expect_true(all(att_pp$att_lower - 0.5 <= true_vals &
                      true_vals <= att_pp$att_upper + 0.5),
                label = paste("True ATT inside SD-Time IS at period", pp))
  }
})
```

**Step 2: Run to confirm they fail**

```r
source("_local/TEST_MULT_PERIOD/tests/test-sd-time.R")
```

**Step 3: Verify SD-Time code paths were added in Task 6**

SD-Time is already handled in `compute_att_bounds()` (added in Task 3) and in `lpt.R` ATT^o block (Task 6). Run the tests — they should now pass if the implementation was correct. If not, debug the SD-Time path in `compute_att_bounds()`:

The SD-Time center in `compute_att_bounds` for `att_lower`/`att_upper`:
```r
sd_center <- Lambda_d - t * delta_tilde_0_d
sd_lower  <- sd_center - sd_hw
sd_upper  <- sd_center + sd_hw
```
Verify this is correct.

**Step 4: Run tests**

```r
source("_local/TEST_MULT_PERIOD/tests/test-sd-time.R")
```
Expected: all PASS.

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/tests/test-sd-time.R
git commit -m "test(multiperiod): add SD-Time restriction tests"
```

---

## Task 8: Add `"time_varying_selection"` DGP to `simulate_lpt()`

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/simulate.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-simulate-tvs.R`

**Background:** This DGP makes `delta_tilde_s(d)` change over time, useful for testing SD-Time. The selection function is `rho * d * t + rho2 * d * t^2` (where `t` is the integer time index). The period-specific increment is the derivative w.r.t. t, approximately `rho * d + rho2 * d * (2t - 1)`. New parameter `rho2` (default 0.05).

**Step 1: Write failing test**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("simulate_lpt supports time_varying_selection DGP", {
  expect_no_error(
    simulate_lpt(n = 200, n_pre = 3, n_post = 2,
                 dgp = "time_varying_selection", rho = 0.2, rho2 = 0.05, seed = 40)
  )
})

test_that("time_varying_selection has correct true_B and true_delta_tilde", {
  sim <- simulate_lpt(n = 200, n_pre = 3, n_post = 2,
                      dgp = "time_varying_selection", rho = 0.2, rho2 = 0.05, seed = 41)

  # true_B should be abs(rho) (LPT bound)
  expect_equal(attr(sim, "true_B"), abs(0.2))

  # true_delta_tilde_0 should be stored (last pre-period increment)
  expect_true(!is.null(attr(sim, "true_delta_tilde_0")))
})
```

**Step 2: Run to confirm failure**

```r
source("_local/TEST_MULT_PERIOD/tests/test-simulate-tvs.R")
```

**Step 3: Add DGP to `_local/TEST_MULT_PERIOD/R/simulate.R`**

Add `rho2 = 0.05` parameter to the function signature:
```r
simulate_lpt <- function(n = 1000, n_pre = 3, n_post = 1,
                          dgp = c("no_selection", "linear_selection",
                                  "time_varying_selection"),
                          beta = c(1, -0.5), rho = 0.3, rho2 = 0.05,
                          sigma_eps = 1, seed = NULL) {
```

Add to `selection_fn` switch:
```r
  selection_fn <- switch(dgp,
    "no_selection"          = function(d, t) rep(0, length(d)),
    "linear_selection"      = function(d, t) rho * d * t,
    "time_varying_selection"= function(d, t) rho * d * t + rho2 * d * t^2
  )
```

Add to `true_B` switch:
```r
  true_B <- switch(dgp,
    "no_selection"           = 0,
    "linear_selection"       = abs(rho),
    "time_varying_selection" = abs(rho)  # Lipschitz bound is still |rho| (slope w.r.t. d)
  )
```

After `attr(result, "true_B") <- true_B`, add:
```r
  # For time_varying_selection: true delta_tilde_0(d)
  # The last pre-period is t_idx = n_pre. Its increment relative to t_idx = n_pre - 1.
  # delta_tilde_{n_pre}(d) = [selection(d, n_pre) - selection(d, n_pre-1)]
  #   = rho*d + rho2*d*(2*n_pre - 1)   (for time_varying_selection)
  #   = rho*d                            (for linear_selection)
  if (dgp == "time_varying_selection") {
    true_delta_tilde_0 <- function(d) rho * d + rho2 * d * (2 * n_pre - 1)
  } else {
    true_delta_tilde_0 <- function(d) rho * d
  }
  attr(result, "true_delta_tilde_0") <- true_delta_tilde_0
  attr(result, "dgp") <- dgp
```

**Step 4: Run tests**

```r
source("_local/TEST_MULT_PERIOD/tests/test-simulate-tvs.R")
```
Expected: both PASS.

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/simulate.R _local/TEST_MULT_PERIOD/tests/test-simulate-tvs.R
git commit -m "feat(multiperiod): add time_varying_selection DGP to simulate_lpt"
```

---

## Task 9: Update `summary.R` to Report Multi-Period Info

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/summary.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-summary-multiperiod.R`

**Step 1: Write failing test**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("summary reports time_restriction and t multipliers", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 50)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = attr(sim, "true_B"),
             time_restriction = "rm", M_bar = 1)

  # Summary output should mention time_restriction
  out <- capture.output(summary(fit))
  expect_true(any(grepl("rm", out, ignore.case = TRUE)))
  expect_true(any(grepl("t\\s*=\\s*[12]", out)))  # t = 1, t = 2 mentioned
})

test_that("summary mentions binding restriction when applicable", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 51)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = 4L, pre_periods = 1:3,
             B = attr(sim, "true_B"),
             time_restriction = "rm", M_bar = 1)

  out <- capture.output(summary(fit))
  # Should mention binding_lower or binding_upper somewhere
  expect_true(any(grepl("binding|bind", out, ignore.case = TRUE)))
})
```

**Step 2: Run to confirm failure**

```r
source("_local/TEST_MULT_PERIOD/tests/test-summary-multiperiod.R")
```

**Step 3: Update `summary.lpt()` in `_local/TEST_MULT_PERIOD/R/summary.R`**

At the top of the body, after the `rule()` / header, add a "Time restriction" block:
```r
  # --- Time restriction info ---
  tr <- object$specifications$time_restriction
  if (!is.null(tr) && tr != "none") {
    if (tr == "rm") {
      cat(sprintf("\n  Time restriction: RM-Time (M_bar = %.4f)\n",
                  object$specifications$M_bar))
    } else if (tr == "sd") {
      cat(sprintf("\n  Time restriction: SD-Time (B_t = %.4f)\n",
                  object$specifications$B_t))
    }
  } else {
    cat("\n  Time restriction: none (LPT only)\n")
  }

  # Report t values
  if (!is.null(object$t_values)) {
    cat("  t multipliers per post-period:\n")
    for (nm in names(object$t_values)) {
      cat(sprintf("    period %s: t = %d\n", nm, object$t_values[[nm]]))
    }
  }
```

In the ATT identified sets section, update to show binding restriction when available. After printing `[ATT lower, upper]`, add a column if `binding_lower` exists:
```r
      has_binding <- "binding_lower" %in% names(att_pp)
      hdr_fmt <- if (has_binding) {
        "  %-10s %-12s %-22s %-12s %-10s\n"
      } else {
        "  %-10s %-12s %-22s %-10s\n"
      }

      if (has_binding) {
        cat(sprintf(hdr_fmt, "Dose", "Lambda(d)", "[ATT lower, upper]",
                    "Binding", "Excl. 0?"))
      } else {
        cat(sprintf(hdr_fmt, "Dose", "Lambda(d)", "[ATT lower, upper]", "Excl. 0?"))
      }

      for (qt in qtiles) {
        idx <- which.min(abs(dose_vals - qt))
        row <- att_pp[idx, ]
        excludes_zero <- (row$att_lower > 0) | (row$att_upper < 0)
        if (has_binding) {
          bind_str <- paste0(row$binding_lower, "/", row$binding_upper)
          cat(sprintf("  %-10.3f %-12.4f [%-8.4f, %-8.4f]  %-12s %-10s\n",
                      row$d, row$Lambda_d, row$att_lower, row$att_upper,
                      bind_str, ifelse(excludes_zero, "Yes", "No")))
        } else {
          cat(sprintf("  %-10.3f %-12.4f [%-8.4f, %-8.4f]  %-10s\n",
                      row$d, row$Lambda_d, row$att_lower, row$att_upper,
                      ifelse(excludes_zero, "Yes", "No")))
        }
      }
```

**Step 4: Run tests**

```r
source("_local/TEST_MULT_PERIOD/tests/test-summary-multiperiod.R")
```
Expected: both PASS.

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/summary.R _local/TEST_MULT_PERIOD/tests/test-summary-multiperiod.R
git commit -m "feat(multiperiod): update summary.lpt to report time_restriction and t multipliers"
```

---

## Task 10: Update `plot.R` for New Parameters

**Files:**
- Modify: `_local/TEST_MULT_PERIOD/R/plot.R`
- Create: `_local/TEST_MULT_PERIOD/tests/test-plot-multiperiod.R`

**Changes needed:**
1. `plot_att()`: if `binding_lower`/`binding_upper` columns exist, color the ribbon by which restriction binds (or add a separate visual indicator)
2. `plot_sensitivity()`: when `time_restriction = "rm"`, x-axis should be `M_bar / M_bar_hat`; when `"sd"`, x-axis = `B_t`
3. `plot_datt()`: dATT is only bounded by LPT; add note when `time_restriction != "none"` and slope is unbounded

**Step 1: Write failing test**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("plot.lpt works with time_restriction='rm'", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_lpt(n = 300, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 60)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = attr(sim, "true_B"),
             time_restriction = "rm", M_bar = 1)

  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "sensitivity", estimand = "att_o"))
})

test_that("sensitivity plot with rm uses M_bar axis", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_lpt(n = 300, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 61)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = 4L, pre_periods = 1:3,
             B = 0,
             time_restriction = "rm", M_bar = 1)

  p <- plot(fit, type = "sensitivity", estimand = "att_o")
  # The x-axis label should reference M_bar
  x_label <- as.character(p$labels$x)
  expect_true(grepl("M_bar|M", x_label, ignore.case = TRUE))
})
```

**Step 2: Run to confirm failure**

```r
source("_local/TEST_MULT_PERIOD/tests/test-plot-multiperiod.R")
```

**Step 3: Update `plot_sensitivity()` in `_local/TEST_MULT_PERIOD/R/plot.R`**

At the top of `plot_sensitivity()`, read the time_restriction from specifications:
```r
plot_sensitivity <- function(x, d0, B_grid, col_band, col_line, col_marker,
                              estimand = "datt", period = NULL) {

  tr   <- x$specifications$time_restriction %||% "none"
  # Helper for NULL coalescing
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # For RM-Time with B=0, use M_bar sensitivity axis
  if (tr == "rm" && x$B_hat == 0) {
    return(plot_sensitivity_rm(x, d0, B_grid, col_band, col_line, col_marker,
                               estimand, period))
  }
  if (tr == "sd") {
    return(plot_sensitivity_sd(x, d0, B_grid, col_band, col_line, col_marker,
                               estimand, period))
  }

  # ... existing logic for LPT-only sensitivity ...
```

Add helper `plot_sensitivity_rm()`:
```r
plot_sensitivity_rm <- function(x, d0, B_grid, col_band, col_line, col_marker,
                                 estimand, period) {
  pp     <- if (is.null(period)) x$post_periods[1] else period
  M_hat  <- x$specifications$M_bar
  M_grid <- seq(0, 3 * M_hat, length.out = 50)
  t_val  <- x$t_values[[as.character(pp)]]

  if (estimand == "att_o") {
    if (is.null(x$att_o)) stop("ATT^o not available.")
    att_o_pp   <- x$att_o[x$att_o$period == pp, ]
    center_val <- att_o_pp$att_o_bin[1]
    # delta_bar_star_pre is stored via half-width: hw = t * M_bar * delta_bar_star
    # Recover delta_bar_star from the stored row (M_bar row)
    hw_stored      <- (att_o_pp$att_o_upper[1] - att_o_pp$att_o_lower[1]) / 2
    delta_bar_star <- hw_stored / (t_val * M_hat)

    sens_df <- data.frame(
      M        = M_grid,
      M_ratio  = M_grid / M_hat,
      is_lower = center_val - t_val * M_grid * delta_bar_star,
      is_upper = center_val + t_val * M_grid * delta_bar_star
    )
    y_label    <- expression(ATT^o)
    plot_title <- sprintf("Sensitivity: ATT^o (RM-Time, t = %d)", t_val)

  } else {
    stop("For RM-Time sensitivity, currently only estimand = 'att_o' is supported.")
  }

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = .data$M_ratio)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$is_lower, ymax = .data$is_upper),
                         fill = col_band, alpha = 0.35) +
    ggplot2::geom_hline(yintercept = sens_df$is_lower[1] + (sens_df$is_upper[1] - sens_df$is_lower[1])/2,
                         color = col_line, linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                         color = col_marker, linewidth = 0.7) +
    ggplot2::annotate("text", x = 1, y = max(sens_df$is_upper),
                       label = expression(hat(M)[bar]),
                       vjust = -0.5, hjust = -0.1, color = col_marker, size = 3.5) +
    ggplot2::labs(x = expression(M[bar] / hat(M)[bar]), y = y_label, title = plot_title) +
    ggplot2::theme_minimal(base_size = 12)
  p
}
```

Add a similar `plot_sensitivity_sd()` for SD-Time (x-axis = B_t, or B_t over some reference).

**Step 4: Run tests**

```r
source("_local/TEST_MULT_PERIOD/tests/test-plot-multiperiod.R")
```
Expected: both PASS.

**Step 5: Commit**
```bash
git add _local/TEST_MULT_PERIOD/R/plot.R _local/TEST_MULT_PERIOD/tests/test-plot-multiperiod.R
git commit -m "feat(multiperiod): update plot.lpt for RM-Time and SD-Time sensitivity plots"
```

---

## Task 11: Integration Tests (Tests 1–6 from Spec)

**Files:**
- Create: `_local/TEST_MULT_PERIOD/tests/test-integration.R`

**Step 1: Write all integration tests**

```r
source("_local/TEST_MULT_PERIOD/tests/helpers.R")

# --- Test 1: LPT multiplier ---
test_that("Test 1 - LPT t multiplier: IS contains true ATT for all (d,t)", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 5,
                      dgp = "linear_selection", seed = 42)
  B_true  <- attr(sim, "true_B")
  true_att <- attr(sim, "true_att")

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = 4:8, pre_periods = 1:3, B = B_true)

  for (pp in 4:8) {
    att_pp  <- fit$att[fit$att$period == pp & fit$att$B == B_true, ]
    t_val   <- fit$t_values[[as.character(pp)]]
    true_vals <- true_att(att_pp$d)

    # True ATT inside IS (with small tolerance for finite sample)
    expect_true(all(att_pp$att_lower <= true_vals + 0.5),
                label = paste("lower bound, period", pp))
    expect_true(all(true_vals <= att_pp$att_upper + 0.5),
                label = paste("upper bound, period", pp))

    # Half-width = t * B * d
    hw_expected <- t_val * B_true * att_pp$d
    hw_actual   <- (att_pp$att_upper - att_pp$att_lower) / 2
    expect_equal(hw_actual, hw_expected, tolerance = 1e-10,
                 label = paste("half-width formula at period", pp))

    # For t > 1, bounds wider than at t=1
    if (t_val > 1) {
      att_t1  <- fit$att[fit$att$period == 4L & fit$att$B == B_true, ]
      hw_t1   <- (att_t1$att_upper - att_t1$att_lower) / 2
      d_common <- intersect(att_pp$d, att_t1$d)
      if (length(d_common) > 0) {
        hw_now  <- (att_pp$att_upper - att_pp$att_lower)[att_pp$d %in% d_common] / 2
        hw_base <- hw_t1[att_t1$d %in% d_common] / 2
        expect_true(all(hw_now >= hw_base - 1e-10),
                    label = paste("t>1 wider than t=1 at period", pp))
      }
    }
  }
})

# --- Test 2: RM-Time ---
test_that("Test 2 - RM-Time: delta_star_pre ≈ |rho| * d, IS contains true ATT", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 42)
  B_true   <- attr(sim, "true_B")  # 0.3
  true_att <- attr(sim, "true_att")

  fit_rm <- lpt(sim, "id", "time", "outcome", "dose",
                post_period = c(4L, 5L), pre_periods = 1:3,
                B = B_true, time_restriction = "rm", M_bar = 1)

  # delta_star_pre(d) ≈ |rho| * d = 0.3 * d
  dsp <- fit_rm$calibration$delta_star_pre
  expect_equal(dsp$delta_star_pre, abs(0.3 * dsp$d), tolerance = 0.15)

  # With M_bar=1, RM half-width ≈ LPT half-width for this DGP
  # True ATT inside IS
  for (pp in c(4L, 5L)) {
    att_pp    <- fit_rm$att[fit_rm$att$period == pp & fit_rm$att$B == B_true, ]
    true_vals <- true_att(att_pp$d)
    expect_true(all(att_pp$att_lower <= true_vals + 0.5),
                label = paste("lower, period", pp))
    expect_true(all(true_vals <= att_pp$att_upper + 0.5),
                label = paste("upper, period", pp))
  }
})

# --- Test 3: SD-Time ---
test_that("Test 3 - SD-Time B_t=0: point identification in linear_selection DGP", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3,
                      beta = c(1, 0), seed = 42)
  true_att <- attr(sim, "true_att")

  fit_sd <- lpt(sim, "id", "time", "outcome", "dose",
                post_period = c(4L, 5L), pre_periods = 1:3,
                B = 0, time_restriction = "sd", B_t = 0)

  for (pp in c(4L, 5L)) {
    att_pp <- fit_sd$att[fit_sd$att$period == pp & fit_sd$att$B == 0, ]
    hw <- (att_pp$att_upper - att_pp$att_lower) / 2
    expect_equal(hw, rep(0, length(hw)), tolerance = 1e-10,
                 label = paste("zero half-width at period", pp))

    # Center ≈ true ATT (up to finite-sample error)
    expect_equal(att_pp$att_lower, true_att(att_pp$d), tolerance = 0.3,
                 label = paste("center ≈ true ATT at period", pp))
  }

  # With B_t > 0, set widens
  fit_sd2 <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = c(4L, 5L), pre_periods = 1:3,
                 B = 0, time_restriction = "sd", B_t = 0.2)
  att_t2_sd2 <- fit_sd2$att[fit_sd2$att$period == 5L & fit_sd2$att$B == 0, ]
  hw_t2 <- (att_t2_sd2$att_upper - att_t2_sd2$att_lower) / 2
  expect_equal(unique(hw_t2), 0.2 * 2 * 3 / 2, tolerance = 1e-10)
})

# --- Test 4: Joint restriction tightening ---
test_that("Test 4 - Joint LPT + RM is tighter than either alone", {
  sim <- simulate_lpt(n = 1000, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 99)
  B_true <- attr(sim, "true_B")

  fit_lpt   <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3, B = B_true)
  fit_rm    <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3,
                   B = 0, time_restriction = "rm", M_bar = 1)
  fit_joint <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3,
                   B = B_true, time_restriction = "rm", M_bar = 1)

  hw_lpt   <- (fit_lpt$att$att_upper - fit_lpt$att$att_lower) / 2
  hw_rm    <- (fit_rm$att$att_upper - fit_rm$att$att_lower) / 2
  hw_joint <- (fit_joint$att$att_upper - fit_joint$att$att_lower) / 2

  expect_true(all(hw_joint <= hw_lpt + 1e-10))
  expect_true(all(hw_joint <= hw_rm + 1e-10))
})

# --- Test 5: Real SRU data ---
test_that("Test 5 - SRU data: later periods have wider IS", {
  skip_if(!file.exists("data/sru.rda"), "sru data not found")
  data(sru)

  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2005, 2010, 2015, 2019),
             pre_periods = 1993:1999,
             B = "calibrate",
             time_restriction = "rm", M_bar = 1)

  # Later post-periods should have wider ATT^o IS
  get_hw_atto <- function(pp) {
    row <- fit$att_o[fit$att_o$period == pp & fit$att_o$B == fit$B_hat, ]
    (row$att_o_upper - row$att_o_lower) / 2
  }

  hw_2005 <- get_hw_atto(2005)
  hw_2019 <- get_hw_atto(2019)

  expect_gt(hw_2019, hw_2005)

  # t values should be 6, 11, 16, 20
  expect_equal(fit$t_values[["2005"]], 6L)
  expect_equal(fit$t_values[["2019"]], 20L)
})

# --- Test 6: Backward compatibility ---
test_that("Test 6 - Backward compat: single period, time_restriction='none' works", {
  skip_if(!file.exists("data/sru.rda"), "sru data not found")
  data(sru)

  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0.1, time_restriction = "none")

  expect_s3_class(fit, "lpt")
  expect_equal(fit$t_values[["2019"]], 20L)

  # dATT half-width should be t * B = 20 * 0.1 = 2.0
  hw_datt <- unique((fit$datt$datt_upper - fit$datt$datt_lower) / 2)
  expect_equal(hw_datt, 20 * 0.1, tolerance = 1e-10)

  # ATT^o half-width should be t * B * D_bar = 20 * 0.1 * D_bar
  atto_row <- fit$att_o[fit$att_o$B == 0.1, ]
  hw_atto  <- (atto_row$att_o_upper - atto_row$att_o_lower) / 2
  expect_equal(hw_atto, 20 * 0.1 * atto_row$D_bar, tolerance = 1e-10)
})
```

**Step 2: Run all integration tests**

```r
source("_local/TEST_MULT_PERIOD/tests/test-integration.R")
```
Expected: all 6 tests PASS.

**Step 3: Fix any failures** — follow `superpowers:systematic-debugging` skill if tests fail.

**Step 4: Commit**
```bash
git add _local/TEST_MULT_PERIOD/tests/test-integration.R
git commit -m "test(multiperiod): add integration tests 1-6 from spec"
```

---

## Task 12: Full Test Suite Run and Documentation

**Files:**
- Create: `_local/TEST_MULT_PERIOD/run_all_tests.R`

**Step 1: Create `run_all_tests.R`**

```r
# _local/TEST_MULT_PERIOD/run_all_tests.R
# Run from package root: source("_local/TEST_MULT_PERIOD/run_all_tests.R")

library(testthat)
source("_local/TEST_MULT_PERIOD/source_all.R")

test_files <- list.files("_local/TEST_MULT_PERIOD/tests",
                          pattern = "^test-.*\\.R$", full.names = TRUE)

results <- lapply(test_files, function(f) {
  cat("\n=== Running", basename(f), "===\n")
  tryCatch(
    source(f, local = new.env(parent = globalenv())),
    error = function(e) cat("ERROR:", conditionMessage(e), "\n")
  )
})

cat("\n=== All test files run ===\n")
```

**Step 2: Run the full test suite**

```r
source("_local/TEST_MULT_PERIOD/run_all_tests.R")
```
Expected: all tests PASS with no errors.

**Step 3: Fix any remaining failures**, using `superpowers:systematic-debugging` for each.

**Step 4: Final commit**
```bash
git add _local/TEST_MULT_PERIOD/run_all_tests.R
git commit -m "feat(multiperiod): complete multi-period IS implementation with full test suite"
```

---

## Summary of All Modified Files

| File | Key Changes |
|------|-------------|
| `_local/TEST_MULT_PERIOD/R/lpt.R` | Add `time_restriction`, `M_bar`, `B_t` args; compute `t_values`; pass `t` to bounds functions; store pre-period `delta_tilde` via calibration |
| `_local/TEST_MULT_PERIOD/R/att_bounds.R` | Accept `t`, `time_restriction`, `M_bar`, `delta_star_pre_d`, `B_t`, `delta_tilde_0_d`; compute bounds; intersect; add `binding_lower`/`binding_upper` columns |
| `_local/TEST_MULT_PERIOD/R/calibrate_B.R` | Add `delta_tilde_d` column to `pre_slopes`; compute `delta_star_pre` and `delta_tilde_0` in return list |
| `_local/TEST_MULT_PERIOD/R/simulate.R` | Add `"time_varying_selection"` DGP; add `rho2` param; store `true_delta_tilde_0` attribute |
| `_local/TEST_MULT_PERIOD/R/plot.R` | Add `plot_sensitivity_rm()` and `plot_sensitivity_sd()` helpers; route by `time_restriction` |
| `_local/TEST_MULT_PERIOD/R/summary.R` | Report `time_restriction`, `M_bar`/`B_t`, `t_values`, `binding_lower`/`binding_upper` |

## Key Invariants to Verify

1. When `time_restriction = "none"` and single post-period: `t=1`, so all half-widths are `1*B*d = B*d` (backward compatible except for the explicit `t` column)
2. When `B = 0` and `time_restriction = "rm"`: pure RM-Time bounds, no LPT component
3. When `B = 0` and `time_restriction = "sd"` and `B_t = 0`: point identification
4. Joint intersection never gives wider bounds than either restriction alone
5. `t_values` is always a named integer vector with same names as `post_periods`
