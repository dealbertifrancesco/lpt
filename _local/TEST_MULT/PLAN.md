# Three-Type LPT (C/P/S) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend TEST_MULT to implement all three multi-period LPT assumptions (LPT-C, LPT-P, LPT-S), renaming existing "a"/"b" types to "C"/"P" and adding the new LPT-S smooth-curvature type.

**Architecture:** Approach A (minimal diff). `calibrate_B.R` gets a `"smooth"` mode that computes b_0 and Ĉ from first-difference pre-period fits. `lpt.R` renames "a"→"C", "b"→"P", adds "S" with Φ(t) = t·b₀ + Ĉ·t(t+1)/2 as the per-period multiplier. `summary.R` and `plot.R` get updated labels and LPT-S branches. All work in `_local/TEST_MULT/`; the main `R/` package is never touched.

**Tech Stack:** R, mgcv (splines), ggplot2 (plots), base R testing via `stopifnot()`.

**Reference:** Dealberti (2026), §6.3 and §7.2. Design: `docs/plans/2026-03-01-lpt-three-types-design.md`.

---

## Task 1: calibrate_B.R — Add "smooth" mode

**Files:**
- Modify: `_local/TEST_MULT/R/calibrate_B.R`
- Test: `_local/TEST_MULT/tests/test-calibrate-cumulative.R` (add smooth tests)

### Step 1: Add smooth tests to test-calibrate-cumulative.R

Append to the BOTTOM of the existing test file (after the final `cat("calibrate_B cumulative tests PASSED.\n")`):

```r
# --- Test 6: smooth mode returns b0, C_hat, b_s_sequence ---
cal_smooth <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                           pre_periods = 1993:1999, type = "smooth")
stopifnot(is.numeric(cal_smooth$B_hat), cal_smooth$B_hat > 0)   # B_hat = b0
stopifnot(cal_smooth$type == "smooth")
stopifnot(is.numeric(cal_smooth$C_hat), cal_smooth$C_hat >= 0)
stopifnot(is.numeric(cal_smooth$b_s_sequence))
stopifnot(length(cal_smooth$b_s_sequence) == 6)  # 6 consecutive pairs for 7 pre-periods
stopifnot(is.data.frame(cal_smooth$pre_slopes))
stopifnot("period_pair" %in% names(cal_smooth$pre_slopes))
cat(sprintf("  smooth: b0=%.4f, C_hat=%.4f, %d b_s values. OK\n",
            cal_smooth$B_hat, cal_smooth$C_hat, length(cal_smooth$b_s_sequence)))

# --- Test 7: b0 is the LAST entry in b_s_sequence (closest to treatment) ---
stopifnot(cal_smooth$B_hat == cal_smooth$b_s_sequence[length(cal_smooth$b_s_sequence)])
cat("  smooth: b0 == last b_s. OK\n")

# --- Test 8: C_hat is the max abs diff of consecutive b_s ---
diffs <- abs(diff(cal_smooth$b_s_sequence))
stopifnot(abs(cal_smooth$C_hat - max(diffs)) < 1e-12)
cat("  smooth: C_hat == max|diff(b_s)|. OK\n")

# --- Test 9: with 2 pre-periods, C_hat == 0 ---
cal_2_smooth <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                              pre_periods = c(1998, 1999), type = "smooth")
stopifnot(cal_2_smooth$C_hat == 0)
stopifnot(length(cal_2_smooth$b_s_sequence) == 1)
cat("  smooth, 2 pre-periods: C_hat = 0. OK\n")

# --- Test 10: non-smooth types return NULL for C_hat and b_s_sequence ---
stopifnot(is.null(cal_fd$C_hat))
stopifnot(is.null(cal_fd$b_s_sequence))
stopifnot(is.null(cal_cum$C_hat))
stopifnot(is.null(cal_cum$b_s_sequence))
cat("  non-smooth types: C_hat and b_s_sequence are NULL. OK\n")

cat("calibrate_B smooth tests PASSED.\n")
```

### Step 2: Run test to verify it fails

```bash
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-calibrate-cumulative.R', local=TRUE)"
```

Expected: FAIL with "object 'cal_smooth' not found" or similar.

### Step 3: Implement "smooth" mode in calibrate_B.R

Replace the entire `_local/TEST_MULT/R/calibrate_B.R` with:

```r
#' Calibrate the sensitivity parameter B from pre-treatment periods
#'
#' @param data Data frame in long format.
#' @param id_col Character. Unit identifier column.
#' @param time_col Character. Time period column.
#' @param outcome_col Character. Outcome column.
#' @param dose_col Character. Dose column.
#' @param pre_periods Vector. Pre-treatment period identifiers. At least 2.
#' @param eval_points Numeric vector or NULL.
#' @param k Integer. Spline basis dimension. Default: 5.
#' @param spline_bs Character. Spline basis type. Default: "cr".
#' @param type Character. One of "cumulative" (LPT-C), "first_diff" (LPT-P),
#'   "smooth" (LPT-S).
#'
#' @return A list with components B_hat, C_hat (NULL for non-smooth),
#'   b_s_sequence (NULL for non-smooth), pre_slopes, sup_by_period, type.
#' @export
calibrate_B <- function(data, id_col, time_col, outcome_col, dose_col,
                         pre_periods, eval_points = NULL,
                         k = 5, spline_bs = "cr",
                         type = c("cumulative", "first_diff", "smooth")) {

  type <- match.arg(type)
  pre_periods <- sort(pre_periods)
  if (length(pre_periods) < 2) {
    stop("Need at least 2 pre-treatment periods to calibrate B.")
  }

  all_slopes <- list()
  sup_vals   <- numeric(0)
  pair_labels <- character(0)

  if (type == "cumulative") {
    # LPT-C: ALL C(n,2) ordered pairs — estimates slope of cumulative mu(d, k)
    n <- length(pre_periods)
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1L, n)) {
        t0 <- pre_periods[i]; t1 <- pre_periods[j]
        window_len <- j - i
        pair_label <- paste0(t0, "-", t1)
        pair_labels <- c(pair_labels, pair_label)

        dat_t0 <- data[data[[time_col]] == t0, ]
        dat_t1 <- data[data[[time_col]] == t1, ]
        merged <- merge(
          dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
          dat_t1[, c(id_col, outcome_col), drop = FALSE],
          by = id_col, suffixes = c("_early", "_late")
        )
        delta_y  <- merged[[paste0(outcome_col, "_late")]] -
                    merged[[paste0(outcome_col, "_early")]]
        dose_vec <- merged[[dose_col]]

        sr <- estimate_dose_slope(delta_y, dose_vec,
                                  eval_points = eval_points,
                                  k = k, spline_bs = spline_bs)
        if (is.null(eval_points)) eval_points <- sr$eval_points

        all_slopes[[length(all_slopes) + 1]] <- data.frame(
          period_pair = pair_label, window_length = window_len,
          d = sr$eval_points, mu_prime_d = sr$lambda_d, se = sr$se_lambda
        )
        sup_vals <- c(sup_vals, max(abs(sr$lambda_d)))
      }
    }
    names(sup_vals) <- pair_labels
    pre_slopes <- do.call(rbind, all_slopes); rownames(pre_slopes) <- NULL
    return(list(B_hat = max(sup_vals), C_hat = NULL, b_s_sequence = NULL,
                pre_slopes = pre_slopes, sup_by_period = sup_vals, type = type))

  } else {
    # "first_diff" and "smooth": consecutive first-difference pairs
    b_s_vec <- numeric(0)

    for (j in seq_len(length(pre_periods) - 1)) {
      t0 <- pre_periods[j]; t1 <- pre_periods[j + 1]
      pair_label <- paste0(t0, "-", t1)
      pair_labels <- c(pair_labels, pair_label)

      dat_t0 <- data[data[[time_col]] == t0, ]
      dat_t1 <- data[data[[time_col]] == t1, ]
      merged <- merge(
        dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
        dat_t1[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_0", "_1")
      )
      delta_y  <- merged[[paste0(outcome_col, "_1")]] -
                  merged[[paste0(outcome_col, "_0")]]
      dose_vec <- merged[[dose_col]]

      sr <- estimate_dose_slope(delta_y, dose_vec,
                                eval_points = eval_points,
                                k = k, spline_bs = spline_bs)
      if (is.null(eval_points)) eval_points <- sr$eval_points

      b_s <- max(abs(sr$lambda_d))
      b_s_vec <- c(b_s_vec, b_s)

      all_slopes[[length(all_slopes) + 1]] <- data.frame(
        period_pair = pair_label, window_length = 1L,
        d = sr$eval_points, mu_prime_d = sr$lambda_d, se = sr$se_lambda
      )
      sup_vals <- c(sup_vals, b_s)
    }
    names(sup_vals)  <- pair_labels
    names(b_s_vec)   <- pair_labels
    pre_slopes <- do.call(rbind, all_slopes); rownames(pre_slopes) <- NULL

    if (type == "first_diff") {
      return(list(B_hat = max(sup_vals), C_hat = NULL, b_s_sequence = NULL,
                  pre_slopes = pre_slopes, sup_by_period = sup_vals, type = type))
    }

    # type == "smooth"
    # b0 = b_s at last pre-period (s=0, closest to treatment)
    b0    <- b_s_vec[length(b_s_vec)]
    # C_hat = max change in max-curvature across consecutive pre-periods
    C_hat <- if (length(b_s_vec) >= 2) max(abs(diff(b_s_vec))) else 0

    return(list(B_hat = b0, C_hat = C_hat, b_s_sequence = b_s_vec,
                pre_slopes = pre_slopes, sup_by_period = sup_vals, type = type))
  }
}
```

### Step 4: Run test to verify it passes

```bash
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-calibrate-cumulative.R', local=TRUE)"
```

Expected: All tests PASSED.

### Step 5: Commit

```bash
git add _local/TEST_MULT/R/calibrate_B.R _local/TEST_MULT/tests/test-calibrate-cumulative.R
git commit -m "feat: calibrate_B smooth mode computes b0 and C_hat for LPT-S"
```

---

## Task 2: lpt.R — Rename C/P and Add LPT-S

**Files:**
- Modify: `_local/TEST_MULT/R/lpt.R`
- Create: `_local/TEST_MULT/tests/test-lpt-s.R`

### Step 1: Write test-lpt-s.R

```r
# test-lpt-s.R
cat("Testing lpt() with lpt_type = 'S'...\n")

# --- Test 1: LPT-S with B=0 collapses to point identification ---
fit_s0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = c(b0 = 0, C = 0), lpt_type = "S")
stopifnot(inherits(fit_s0, "lpt"))
stopifnot(fit_s0$lpt_type == "S")
datt0 <- fit_s0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
stopifnot(all(abs(datt0$datt_upper - datt0$lambda_d) < 1e-10))
cat("  LPT-S, b0=0, C=0: point identification. OK\n")

# --- Test 2: LPT-S IS width grows with Phi(t) across periods ---
# With b0=0.1, C=0.2:
# Phi(t1) = t1*0.1 + 0.2*t1*(t1+1)/2
# Phi(t2) = t2*0.1 + 0.2*t2*(t2+1)/2
fit_s_multi <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = c(2018, 2019), pre_periods = 1993:1999,
                    B = c(b0 = 0.1, C = 0.2), lpt_type = "S")

datt_2018 <- fit_s_multi$datt[fit_s_multi$datt$period == 2018, ]
datt_2019 <- fit_s_multi$datt[fit_s_multi$datt$period == 2019, ]
t1 <- unique(datt_2018$t_index)
t2 <- unique(datt_2019$t_index)

phi_t1 <- t1 * 0.1 + 0.2 * t1 * (t1 + 1) / 2
phi_t2 <- t2 * 0.1 + 0.2 * t2 * (t2 + 1) / 2

width_2018 <- datt_2018$datt_upper[1] - datt_2018$datt_lower[1]
width_2019 <- datt_2019$datt_upper[1] - datt_2019$datt_lower[1]

stopifnot(abs(width_2018 - 2 * phi_t1) < 1e-10)
stopifnot(abs(width_2019 - 2 * phi_t2) < 1e-10)
stopifnot(width_2019 > width_2018)
cat(sprintf("  LPT-S: Phi(t1=%d)=%.3f->w=%.3f, Phi(t2=%d)=%.3f->w=%.3f. OK\n",
            t1, phi_t1, width_2018, t2, phi_t2, width_2019))

# --- Test 3: LPT-S with calibration uses smooth mode ---
fit_s_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "S")
stopifnot(fit_s_cal$calibration$type == "smooth")
stopifnot(!is.null(fit_s_cal$b0))
stopifnot(!is.null(fit_s_cal$C_hat))
stopifnot(fit_s_cal$B_hat == fit_s_cal$b0)
cat(sprintf("  LPT-S calibrated: b0=%.4f, C_hat=%.4f. OK\n",
            fit_s_cal$b0, fit_s_cal$C_hat))

# --- Test 4: ATT width under LPT-S uses Phi(t)*d ---
stopifnot(!is.null(fit_s_multi$att))
att_2018 <- fit_s_multi$att[fit_s_multi$att$period == 2018, ]
att_2019 <- fit_s_multi$att[fit_s_multi$att$period == 2019, ]
# Width at dose d: 2 * Phi(t) * d
# Check ratio: width_2019 / width_2018 = Phi(t2)/Phi(t1) at same dose
idx <- 25
w18 <- att_2018$att_upper[idx] - att_2018$att_lower[idx]
w19 <- att_2019$att_upper[idx] - att_2019$att_lower[idx]
expected_ratio <- phi_t2 / phi_t1
stopifnot(abs(w19 / w18 - expected_ratio) < 1e-6)
cat(sprintf("  LPT-S: ATT width ratio = %.3f (expected Phi(t2)/Phi(t1) = %.3f). OK\n",
            w19 / w18, expected_ratio))

# --- Test 5: ATT^o width under LPT-S ---
stopifnot(!is.null(fit_s_multi$att_o))
atto_2018 <- fit_s_multi$att_o[fit_s_multi$att_o$period == 2018, ]
atto_2019 <- fit_s_multi$att_o[fit_s_multi$att_o$period == 2019, ]
w_o18 <- atto_2018$att_o_upper - atto_2018$att_o_lower
w_o19 <- atto_2019$att_o_upper - atto_2019$att_o_lower
stopifnot(w_o19 > w_o18)
cat("  LPT-S: ATT^o width grows with t. OK\n")

# --- Test 6: LPT-S with C=0 reduces to LPT-C with B=b0 ---
fit_s_c0 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = c(b0 = 0.3, C = 0), lpt_type = "S")
fit_c_b03 <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = c(2018, 2019), pre_periods = 1993:1999,
                  B = 0.3, lpt_type = "C")
# With C=0, Phi(t) = t*b0. With LPT-C, t_mult=1 and bias=B.
# They are NOT identical in general (S gives t*b0 width, C gives B width).
# But at t=1, LPT-S with C=0 gives Phi(1)=b0 = same as LPT-C with B=b0.
datt_s_2018 <- fit_s_c0$datt[fit_s_c0$datt$period == 2018, ]
datt_c_2018 <- fit_c_b03$datt[fit_c_b03$datt$period == 2018, ]
# Width: LPT-S: 2*t1*b0; LPT-C: 2*B = 2*0.3
phi_t1_c0 <- t1 * 0.3 + 0 * t1 * (t1 + 1) / 2  # = t1 * 0.3
w_s <- datt_s_2018$datt_upper[1] - datt_s_2018$datt_lower[1]
stopifnot(abs(w_s - 2 * phi_t1_c0) < 1e-10)
cat(sprintf("  LPT-S with C=0: width = 2*t*b0 = %.3f. OK\n", 2 * phi_t1_c0))

# --- Test 7: lpt object has b0 and C_hat fields for LPT-S ---
stopifnot(!is.null(fit_s_multi$b0))
stopifnot(!is.null(fit_s_multi$C_hat))
stopifnot(fit_s_multi$b0 == 0.1)
stopifnot(fit_s_multi$C_hat == 0.2)
# LPT-C and LPT-P should have NULL for these
fit_c <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999, B = 0.1, lpt_type = "C")
stopifnot(is.null(fit_c$b0))
stopifnot(is.null(fit_c$C_hat))
cat("  b0/C_hat present for S, NULL for C/P. OK\n")

cat("lpt_type='S' tests PASSED.\n")
```

### Step 2: Run test to verify it fails

```bash
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-lpt-s.R', local=TRUE)"
```

Expected: FAIL — `lpt_type = "S"` not accepted by `match.arg`.

### Step 3: Implement updated lpt.R

Replace entire `_local/TEST_MULT/R/lpt.R` with:

```r
#' Local Parallel Trends estimation
#'
#' @export
lpt <- function(data, id_col, time_col, outcome_col, dose_col,
                post_period, pre_periods = NULL,
                B = "calibrate", eval_points = NULL,
                k = 5, spline_bs = "cr", alpha = 0.05,
                lpt_type = c("C", "P", "S")) {

  lpt_type <- match.arg(lpt_type)

  # --- Input validation ---
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  required_cols <- c(id_col, time_col, outcome_col, dose_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found in data: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  all_times <- sort(unique(data[[time_col]]))
  bad_periods <- setdiff(post_period, all_times)
  if (length(bad_periods) > 0) {
    stop(sprintf("post_period value(s) not found in '%s': %s",
                 time_col, paste(bad_periods, collapse = ", ")))
  }

  the_call <- match.call()

  # --- Period structure ---
  post_periods <- sort(post_period)

  if (!is.null(pre_periods)) {
    pre_period_set <- sort(pre_periods)
    bad_pre <- setdiff(pre_period_set, all_times)
    if (length(bad_pre) > 0) {
      stop(sprintf("pre_periods value(s) not found in '%s': %s",
                   time_col, paste(bad_pre, collapse = ", ")))
    }
    ref_period <- max(pre_period_set)
  } else {
    pre_period_set <- all_times[all_times < min(post_periods)]
    ref_period <- max(pre_period_set)
  }

  if (length(pre_period_set) < 1) {
    stop("Need at least one pre-period before the earliest post-period.")
  }

  # --- Compute t_index for each post-period ---
  ref_rank <- which(all_times == ref_period)
  t_index_map <- vapply(post_periods, function(pp) {
    as.integer(which(all_times == pp) - ref_rank)
  }, integer(1L))
  names(t_index_map) <- as.character(post_periods)

  # --- Reference-period data ---
  ref_data  <- data[data[[time_col]] == ref_period, ]
  dose_vec  <- ref_data[[dose_col]]
  has_untreated <- any(dose_vec == 0)

  # --- Estimate dose slope for each post-period ---
  slopes <- list()
  for (pp in post_periods) {
    post_data <- data[data[[time_col]] == pp, ]
    merged <- merge(
      ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
      post_data[, c(id_col, outcome_col), drop = FALSE],
      by = id_col, suffixes = c("_ref", "_post")
    )
    delta_y <- merged[[paste0(outcome_col, "_post")]] -
               merged[[paste0(outcome_col, "_ref")]]
    wd <- merged[[dose_col]]
    sr <- estimate_dose_slope(delta_y = delta_y, dose = wd,
                              eval_points = eval_points,
                              k = k, spline_bs = spline_bs)
    if (is.null(eval_points) && length(slopes) == 0) eval_points <- sr$eval_points
    slopes[[as.character(pp)]] <- sr
  }

  # --- Handle B ---
  calibration_result <- NULL
  b0_val  <- NULL   # only used for LPT-S
  C_hat   <- NULL   # only used for LPT-S

  if (is.character(B) && B == "calibrate") {
    if (length(pre_period_set) < 2) {
      stop("B = 'calibrate' requires at least 2 pre-treatment periods.")
    }
    cal_type <- switch(lpt_type,
      "C" = "cumulative",
      "P" = "first_diff",
      "S" = "smooth"
    )
    calibration_result <- calibrate_B(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      pre_periods = pre_period_set,
      eval_points = slopes[[1]]$eval_points,
      k = k, spline_bs = spline_bs, type = cal_type
    )
    if (lpt_type == "S") {
      b0_val  <- calibration_result$B_hat   # B_hat stores b0 for smooth
      C_hat   <- calibration_result$C_hat
      B_hat   <- b0_val
      B_values <- b0_val
      message(sprintf("Calibrated b0 = %.4f, C = %.4f (%s) from %d pre-period(s).",
                      b0_val, C_hat, cal_type, length(pre_period_set) - 1))
    } else {
      B_hat   <- calibration_result$B_hat
      B_values <- B_hat
      message(sprintf("Calibrated B = %.4f (%s) from %d pre-period(s).",
                      B_hat, cal_type, length(pre_period_set) - 1))
    }

  } else if (is.numeric(B)) {
    if (lpt_type == "S") {
      # B should be c(b0 = ..., C = ...) or length-2 unnamed vector
      if (!is.null(names(B)) && all(c("b0", "C") %in% names(B))) {
        b0_val <- B[["b0"]]; C_hat <- B[["C"]]
      } else if (length(B) == 2) {
        b0_val <- B[1]; C_hat <- B[2]
      } else if (length(B) == 1) {
        b0_val <- B[1]; C_hat <- 0
      } else {
        stop("For lpt_type='S', B must be c(b0=..., C=...) or length 1-2 numeric.")
      }
      if (b0_val < 0 || C_hat < 0) stop("b0 and C must be non-negative.")
      B_hat   <- b0_val
      B_values <- b0_val
    } else {
      B_values <- as.numeric(B)
      if (any(!is.finite(B_values)) || any(B_values < 0)) {
        stop("B must be non-negative finite numeric values.")
      }
      B_hat <- max(B_values)
    }
  } else {
    stop("B must be numeric or 'calibrate'.")
  }

  # --- Construct identified sets per period ---
  z_alpha  <- stats::qnorm(1 - alpha / 2)
  datt_all <- list(); att_all <- list(); att_o_all <- list()

  for (pp in post_periods) {
    pp_char <- as.character(pp)
    sr      <- slopes[[pp_char]]
    ep      <- sr$eval_points
    lambda_d  <- sr$lambda_d
    se_lambda <- sr$se_lambda
    t_idx     <- t_index_map[[pp_char]]

    # Compute per-period multiplier
    if (lpt_type == "C") {
      t_mult <- 1L
    } else if (lpt_type == "P") {
      t_mult <- t_idx
    } else {
      # LPT-S: Phi(t) = t*b0 + C*t*(t+1)/2
      t_mult <- t_idx * b0_val + C_hat * t_idx * (t_idx + 1) / 2
    }

    # For LPT-S, the IS width = t_mult (Phi_t) regardless of b.
    # Pass b_eff=1 so that bias_margin = t_mult * 1 = Phi_t.
    b_loop <- if (lpt_type == "S") 1 else B_values

    for (b in b_loop) {
      bias_margin <- t_mult * b
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period     = pp, t_index  = t_idx,
        d          = ep, lambda_d = lambda_d, se_lambda = se_lambda,
        B          = if (lpt_type == "S") B_hat else b,
        datt_lower = lambda_d - bias_margin,
        datt_upper = lambda_d + bias_margin,
        ci_lower   = lambda_d - bias_margin - z_alpha * se_lambda,
        ci_upper   = lambda_d + bias_margin + z_alpha * se_lambda
      )
    }

    if (has_untreated) {
      # ATT level bounds: bias_hw = t_mult * b_eff * d
      b_for_att <- if (lpt_type == "S") 1 else B_values
      att_pp <- compute_att_bounds(sr, b_for_att, dose_vec,
                                    period = pp, t_multiplier = t_mult)
      if (lpt_type == "S") att_pp$B <- B_hat
      att_pp$t_index <- t_idx
      att_all[[length(att_all) + 1]] <- att_pp
    }

    if (has_untreated) {
      post_data <- data[data[[time_col]] == pp, ]
      merged_atto <- merge(
        ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
        post_data[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_ref", "_post")
      )
      dy       <- merged_atto[[paste0(outcome_col, "_post")]] -
                  merged_atto[[paste0(outcome_col, "_ref")]]
      d_merged <- merged_atto[[dose_col]]
      treated_idx   <- d_merged > 0
      untreated_idx <- d_merged == 0
      att_o_bin <- mean(dy[treated_idx]) - mean(dy[untreated_idx])
      D_bar     <- mean(d_merged[treated_idx])

      b_for_atto <- if (lpt_type == "S") 1 else B_values
      for (b in b_for_atto) {
        bias_margin_o <- t_mult * b * D_bar
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period      = pp, t_index    = t_idx,
          att_o_bin   = att_o_bin, D_bar = D_bar,
          B           = if (lpt_type == "S") B_hat else b,
          att_o_lower = att_o_bin - bias_margin_o,
          att_o_upper = att_o_bin + bias_margin_o
        )
      }
    }
  }

  datt  <- do.call(rbind, datt_all);  rownames(datt)  <- NULL
  att   <- if (length(att_all)   > 0) { r <- do.call(rbind, att_all);   rownames(r) <- NULL; r } else NULL
  att_o <- if (length(att_o_all) > 0) { r <- do.call(rbind, att_o_all); rownames(r) <- NULL; r } else NULL

  if (!has_untreated) message("No untreated units (dose = 0). ATT level bounds not computed.")

  structure(
    list(
      datt = datt, att = att, att_o = att_o,
      B_hat = B_hat, B_values = B_values,
      b0 = b0_val, C_hat = C_hat,
      calibration = calibration_result,
      slopes = slopes, call = the_call,
      n = nrow(ref_data), has_untreated = has_untreated,
      post_periods = post_periods, ref_period = ref_period,
      lpt_type = lpt_type, t_index_map = t_index_map,
      specifications = list(k = k, spline_bs = spline_bs, alpha = alpha,
                            post_periods = post_periods, ref_period = ref_period,
                            lpt_type = lpt_type)
    ),
    class = "lpt"
  )
}
```

### Step 4: Run test to verify it passes

```bash
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-lpt-s.R', local=TRUE)"
```

Expected: All tests PASSED.

### Step 5: Commit

```bash
git add _local/TEST_MULT/R/lpt.R _local/TEST_MULT/tests/test-lpt-s.R
git commit -m "feat: lpt() supports lpt_type C/P/S with Phi(t) multiplier for LPT-S"
```

---

## Task 3: summary.R — Update Labels for C/P/S

**Files:**
- Modify: `_local/TEST_MULT/R/summary.R`

No new test needed — the summary output is checked in test-summary-plot.R (Task 6).

Replace entire `_local/TEST_MULT/R/summary.R` with:

```r
#' Print an lpt object
#' @method print lpt
#' @export
print.lpt <- function(x, ...) {
  b_source  <- if (!is.null(x$calibration)) "calibrated" else "user-supplied"
  lpt_label <- if (!is.null(x$lpt_type)) sprintf("LPT-%s", x$lpt_type) else "LPT"
  n_post    <- length(x$post_periods)
  post_str  <- if (n_post == 1) as.character(x$post_periods) else
    sprintf("%d periods (%s)", n_post, paste(x$post_periods, collapse = ", "))

  if (!is.null(x$lpt_type) && x$lpt_type == "S") {
    cat(sprintf("%s | n = %d | b0 = %.4f, C = %.4f (%s) | Post: %s\n",
                lpt_label, x$n, x$b0, x$C_hat, b_source, post_str))
  } else {
    cat(sprintf("%s | n = %d | B = %.4f (%s) | Post: %s\n",
                lpt_label, x$n, x$B_hat, b_source, post_str))
  }
  invisible(x)
}

#' Summarize an lpt object
#' @method summary lpt
#' @export
summary.lpt <- function(object, ...) {
  rule <- function(char = "=", width = 60)
    cat(paste(rep(char, width), collapse = ""), "\n")

  rule()
  lpt_label <- if (!is.null(object$lpt_type)) {
    sprintf("Local Parallel Trends Summary (LPT-%s)", object$lpt_type)
  } else "Local Parallel Trends Summary"
  cat(sprintf("  %s\n", lpt_label))
  rule()

  # --- Type description ---
  if (!is.null(object$lpt_type)) {
    switch(object$lpt_type,
      "C" = {
        cat("\n  LPT-C: smoothness on cumulative trend (Assumption 7).")
        cat("\n  Identified set width is CONSTANT across post-periods.\n")
      },
      "P" = {
        cat("\n  LPT-P: smoothness on period-specific increments (Assumption 8).")
        cat("\n  Identified set width GROWS linearly with horizon t.\n")
      },
      "S" = {
        cat("\n  LPT-S: smooth curvature path (Assumption 9).")
        cat("\n  IS width grows via Phi(t) = t*b0 + C*t*(t+1)/2 (~quadratic in t).\n")
      }
    )
  }

  # --- Sensitivity parameters ---
  if (!is.null(object$lpt_type) && object$lpt_type == "S") {
    b_source <- if (!is.null(object$calibration)) "calibrated" else "user-supplied"
    cat(sprintf("\n  b0 = %.4f, C = %.4f (%s)\n", object$b0, object$C_hat, b_source))
    cat("  b0 = 0, C = 0 is standard parallel trends.\n")
    if (!is.null(object$calibration)) {
      cat("  Pre-period b_s sequence (max|mu_tilde'(d)|):\n")
      for (nm in names(object$calibration$b_s_sequence)) {
        cat(sprintf("    %s: %.4f\n", nm, object$calibration$b_s_sequence[nm]))
      }
    }
  } else {
    if (!is.null(object$calibration)) {
      cal_label <- if (!is.null(object$calibration$type)) {
        sprintf("calibrated, %s", object$calibration$type)
      } else "calibrated"
      cat(sprintf("\n  Sensitivity parameter B: %.4f (%s)\n", object$B_hat, cal_label))
      cat("  Pre-period sup|mu'(d)| by pair/window:\n")
      for (nm in names(object$calibration$sup_by_period)) {
        cat(sprintf("    %s: %.4f\n", nm, object$calibration$sup_by_period[nm]))
      }
    } else {
      cat(sprintf("\n  Sensitivity parameter B: %.4f (user-supplied)\n", object$B_hat))
    }
    cat("  B = 0 is standard parallel trends (point identification).\n")
  }

  # --- Post-period horizons ---
  if (!is.null(object$t_index_map)) {
    cat("\n  Post-period horizons:\n")
    for (nm in names(object$t_index_map)) {
      t_val <- object$t_index_map[[nm]]
      if (!is.null(object$lpt_type) && object$lpt_type == "P") {
        cat(sprintf("    period %s: t = %d (bias multiplier = %d)\n", nm, t_val, t_val))
      } else if (!is.null(object$lpt_type) && object$lpt_type == "S") {
        phi_t <- t_val * object$b0 + object$C_hat * t_val * (t_val + 1) / 2
        cat(sprintf("    period %s: t = %d (Phi(t) = %.4f)\n", nm, t_val, phi_t))
      } else {
        cat(sprintf("    period %s: t = %d\n", nm, t_val))
      }
    }
  }

  # --- ATT^o summary ---
  if (!is.null(object$att_o)) {
    for (pp in object$post_periods) {
      atto_pp <- object$att_o[object$att_o$period == pp &
                                object$att_o$B == object$B_hat, ]
      if (nrow(atto_pp) == 0) next
      cat(sprintf("\n  --- ATT^o overall summary (period %s) ---\n", pp))
      cat(sprintf("  ATT^o_bin (binarized DiD): %.4f\n", atto_pp$att_o_bin))
      cat(sprintf("  Mean dose among treated:   %.4f\n", atto_pp$D_bar))
      cat(sprintf("  IS_{ATT^o}:                [%.4f, %.4f]\n",
                  atto_pp$att_o_lower, atto_pp$att_o_upper))
    }
  }

  # --- dATT results ---
  for (pp in object$post_periods) {
    datt_pp <- object$datt[object$datt$period == pp &
                             object$datt$B == object$B_hat, ]
    if (nrow(datt_pp) == 0) next
    t_str <- if ("t_index" %in% names(datt_pp))
      sprintf(", t=%d", datt_pp$t_index[1]) else ""
    cat(sprintf("\n  --- dATT identified sets (period %s%s) ---\n", pp, t_str))
    dose_vals <- datt_pp$d
    qtiles <- stats::quantile(dose_vals, probs = c(0.25, 0.5, 0.75))
    cat(sprintf("  %-10s %-12s %-22s %-10s\n",
                "Dose", "lambda(d)", "[dATT lower, upper]", "Excl. 0?"))
    for (qt in qtiles) {
      idx <- which.min(abs(dose_vals - qt))
      row <- datt_pp[idx, ]
      cat(sprintf("  %-10.3f %-12.4f [%-8.4f, %-8.4f]  %-10s\n",
                  row$d, row$lambda_d, row$datt_lower, row$datt_upper,
                  ifelse((row$datt_lower > 0) | (row$datt_upper < 0), "Yes", "No")))
    }
  }

  # --- ATT identified sets ---
  if (!is.null(object$att)) {
    for (pp in object$post_periods) {
      att_pp <- object$att[object$att$period == pp &
                             object$att$B == object$B_hat, ]
      if (nrow(att_pp) == 0) next
      t_str <- if ("t_index" %in% names(att_pp))
        sprintf(", t=%d", att_pp$t_index[1]) else ""
      cat(sprintf("\n  --- ATT identified sets (period %s%s) ---\n", pp, t_str))
      dose_vals <- att_pp$d
      qtiles <- stats::quantile(dose_vals, probs = c(0.25, 0.5, 0.75))
      cat(sprintf("  %-10s %-12s %-22s %-10s\n",
                  "Dose", "Lambda(d)", "[ATT lower, upper]", "Excl. 0?"))
      for (qt in qtiles) {
        idx <- which.min(abs(dose_vals - qt))
        row <- att_pp[idx, ]
        cat(sprintf("  %-10.3f %-12.4f [%-8.4f, %-8.4f]  %-10s\n",
                    row$d, row$Lambda_d, row$att_lower, row$att_upper,
                    ifelse((row$att_lower > 0) | (row$att_upper < 0), "Yes", "No")))
      }
    }
  }

  cat(sprintf("\n  n = %d | Spline: %s (k = %d) | alpha = %.2f\n",
              object$n, object$specifications$spline_bs,
              object$specifications$k, object$specifications$alpha))
  rule()
  invisible(object)
}
```

Commit:
```bash
git add _local/TEST_MULT/R/summary.R
git commit -m "feat: summary.R labels LPT-C/P/S and shows Phi(t) for LPT-S"
```

---

## Task 4: plot.R — Sensitivity Plot for LPT-S

**Files:**
- Modify: `_local/TEST_MULT/R/plot.R`

Replace only the `plot_sensitivity` function (lines 218–358) with the version below. Leave all other functions (`plot_datt`, `plot_att`, `plot_pretrends`, `plot.lpt`) unchanged.

```r
plot_sensitivity <- function(x, d0, B_grid, col_band, col_line, col_marker,
                              estimand = "datt", period = NULL) {

  pp      <- if (is.null(period)) x$post_periods[1] else period
  pp_char <- as.character(pp)
  sr      <- x$slopes[[pp_char]]
  ep      <- sr$eval_points

  if (estimand %in% c("datt", "att") && is.null(d0)) {
    stop("'d0' must be specified when estimand = '", estimand, "'.",
         " Supply e.g. d0 = ", round(stats::median(ep), 3), ".")
  }

  t_idx <- x$t_index_map[[pp_char]]
  lpt_t  <- x$lpt_type

  # --- Build sensitivity grid and multiplier ---
  if (lpt_t == "S") {
    # For LPT-S: vary C while holding b0 fixed.
    # B_grid is interpreted as a C-grid.
    C_ref <- if (!is.null(x$C_hat) && x$C_hat > 0) x$C_hat else 1
    if (is.null(B_grid)) B_grid <- seq(0, 3 * C_ref, length.out = 50)
    b0_x  <- if (!is.null(x$b0)) x$b0 else 0
    phi_vec <- t_idx * b0_x + B_grid * t_idx * (t_idx + 1) / 2
    use_ratio <- x$C_hat > 0
    x_var <- if (use_ratio) "B_ratio" else "B"
    x_ratio <- if (use_ratio) B_grid / x$C_hat else B_grid
    x_label <- if (use_ratio) expression(C / hat(C)) else "Sensitivity parameter C"
    vline_x <- 1
    vline_label <- expression(hat(C))
  } else {
    # LPT-C and LPT-P
    t_mult <- if (lpt_t == "P") t_idx else 1L
    if (is.null(B_grid)) {
      B_grid <- if (x$B_hat > 0) seq(0, 3, length.out = 50) * x$B_hat else
        seq(0, 1, length.out = 50)
    }
    phi_vec   <- t_mult * B_grid
    use_ratio <- x$B_hat > 0
    x_var     <- if (use_ratio) "B_ratio" else "B"
    x_ratio   <- if (use_ratio) B_grid / x$B_hat else B_grid
    x_label   <- if (use_ratio) expression(B / hat(B)) else "Sensitivity parameter B"
    vline_x   <- 1
    vline_label <- expression(hat(B))
  }

  z_alpha <- stats::qnorm(1 - x$specifications$alpha / 2)

  if (estimand == "datt") {
    idx        <- which.min(abs(ep - d0))
    center_val <- sr$lambda_d[idx]
    se_val     <- sr$se_lambda[idx]
    d0_actual  <- ep[idx]
    sens_df <- data.frame(
      B = B_grid, B_ratio = x_ratio, center = center_val,
      is_lower = center_val - phi_vec,
      is_upper = center_val + phi_vec,
      ci_lower = center_val - phi_vec - z_alpha * se_val,
      ci_upper = center_val + phi_vec + z_alpha * se_val
    )
    y_label    <- expression(partialdiff * ATT(d) / partialdiff * d ~ "at" ~ d[0])
    plot_title <- sprintf("Sensitivity: dATT at d = %.2f (t=%d)", d0_actual, t_idx)
    show_ci    <- TRUE

  } else if (estimand == "att") {
    if (is.null(x$att)) stop("ATT level bounds not available.")
    att_pp  <- x$att[x$att$period == pp, ]
    att_sub <- att_pp[att_pp$B == att_pp$B[1], ]
    idx        <- which.min(abs(att_sub$d - d0))
    center_val <- att_sub$Lambda_d[idx]
    d0_actual  <- att_sub$d[idx]
    sens_df <- data.frame(
      B = B_grid, B_ratio = x_ratio, center = center_val,
      is_lower = center_val - phi_vec * d0_actual,
      is_upper = center_val + phi_vec * d0_actual,
      ci_lower = NA_real_, ci_upper = NA_real_
    )
    y_label    <- expression(ATT(d[0]))
    plot_title <- sprintf("Sensitivity: ATT(d|d) at d = %.2f (t=%d)", d0_actual, t_idx)
    show_ci    <- FALSE

  } else {  # att_o
    if (is.null(x$att_o)) stop("Overall ATT bounds not available.")
    att_o_pp   <- x$att_o[x$att_o$period == pp, ]
    center_val <- att_o_pp$att_o_bin[1]
    D_bar      <- att_o_pp$D_bar[1]
    sens_df <- data.frame(
      B = B_grid, B_ratio = x_ratio, center = center_val,
      is_lower = center_val - phi_vec * D_bar,
      is_upper = center_val + phi_vec * D_bar,
      ci_lower = NA_real_, ci_upper = NA_real_
    )
    y_label    <- expression(ATT^o)
    plot_title <- sprintf("Sensitivity: Overall ATT (t=%d)", t_idx)
    show_ci    <- FALSE
  }

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = .data[[x_var]]))
  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = col_band, alpha = 0.2)
  }
  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$is_lower, ymax = .data$is_upper),
      fill = col_band, alpha = 0.35) +
    ggplot2::geom_hline(yintercept = sens_df$center[1],
                         color = col_line, linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

  if (use_ratio) {
    p <- p +
      ggplot2::geom_vline(xintercept = vline_x, linetype = "dotted",
                           color = col_marker, linewidth = 0.7) +
      ggplot2::annotate("text", x = vline_x,
                         y = max(sens_df$is_upper, na.rm = TRUE),
                         label = vline_label,
                         vjust = -0.5, hjust = -0.1,
                         color = col_marker, size = 3.5) +
      ggplot2::labs(x = x_label, y = y_label, title = plot_title)
  } else {
    p <- p + ggplot2::labs(x = x_label, y = y_label, title = plot_title)
  }
  p + ggplot2::theme_minimal(base_size = 12)
}
```

Commit:
```bash
git add _local/TEST_MULT/R/plot.R
git commit -m "feat: sensitivity plot uses C/C_hat x-axis for LPT-S"
```

---

## Task 5: Replace test-lpt-a.R with test-lpt-c.R

**Files:**
- Delete: `_local/TEST_MULT/tests/test-lpt-a.R`
- Create: `_local/TEST_MULT/tests/test-lpt-c.R`

Content of `test-lpt-c.R` — identical to current `test-lpt-a.R` with all `"a"` → `"C"`:

```r
# test-lpt-c.R
cat("Testing lpt() with lpt_type = 'C'...\n")

fit_c0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = 0, lpt_type = "C")
stopifnot(inherits(fit_c0, "lpt"))
stopifnot(fit_c0$lpt_type == "C")
stopifnot(fit_c0$B_hat == 0)
stopifnot("t_index" %in% names(fit_c0$datt))
datt0 <- fit_c0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
stopifnot(all(abs(datt0$datt_upper - datt0$lambda_d) < 1e-10))
cat("  LPT-C, B=0, single period: OK\n")

fit_c_multi <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = c(2018, 2019), pre_periods = 1993:1999,
                    B = 0, lpt_type = "C")
stopifnot(length(fit_c_multi$post_periods) == 2)
t_for_2018 <- unique(fit_c_multi$datt$t_index[fit_c_multi$datt$period == 2018])
t_for_2019 <- unique(fit_c_multi$datt$t_index[fit_c_multi$datt$period == 2019])
stopifnot(t_for_2019 > t_for_2018)
cat(sprintf("  LPT-C, multi-period: t(2018)=%d, t(2019)=%d. OK\n", t_for_2018, t_for_2019))

fit_c_b1 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = 0.5, lpt_type = "C")
datt_2018 <- fit_c_b1$datt[fit_c_b1$datt$period == 2018, ]
datt_2019 <- fit_c_b1$datt[fit_c_b1$datt$period == 2019, ]
width_2018 <- datt_2018$datt_upper - datt_2018$datt_lower
width_2019 <- datt_2019$datt_upper - datt_2019$datt_lower
stopifnot(all(abs(width_2018 - 1.0) < 1e-10))
stopifnot(all(abs(width_2019 - 1.0) < 1e-10))
cat("  LPT-C, B=0.5: constant dATT width across periods. OK\n")

fit_c_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "C")
stopifnot(!is.null(fit_c_cal$calibration))
stopifnot(fit_c_cal$calibration$type == "cumulative")
stopifnot(fit_c_cal$B_hat > 0)
cat(sprintf("  LPT-C, calibrated: B_hat=%.4f (cumulative). OK\n", fit_c_cal$B_hat))

att_2018 <- fit_c_b1$att[fit_c_b1$att$period == 2018, ]
att_2019 <- fit_c_b1$att[fit_c_b1$att$period == 2019, ]
att_width_2018 <- att_2018$att_upper - att_2018$att_lower
att_width_2019 <- att_2019$att_upper - att_2019$att_lower
stopifnot(all(abs(att_width_2018 - att_width_2019) < 1e-10))
cat("  LPT-C: constant ATT width across periods. OK\n")

cat("lpt_type='C' tests PASSED.\n")
```

Run:
```bash
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-lpt-c.R', local=TRUE)"
```

Then delete the old file and commit:
```bash
git rm _local/TEST_MULT/tests/test-lpt-a.R
git add _local/TEST_MULT/tests/test-lpt-c.R
git commit -m "test: rename test-lpt-a -> test-lpt-c, update type labels"
```

---

## Task 6: Replace test-lpt-b.R with test-lpt-p.R

Same pattern — create `test-lpt-p.R` from `test-lpt-b.R` with `"b"` → `"P"`:

```r
# test-lpt-p.R
cat("Testing lpt() with lpt_type = 'P'...\n")

fit_p0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = 0, lpt_type = "P")
stopifnot(inherits(fit_p0, "lpt"))
stopifnot(fit_p0$lpt_type == "P")
datt0 <- fit_p0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
cat("  LPT-P, B=0: point identification. OK\n")

fit_p_b1 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = 0.5, lpt_type = "P")
datt_2018 <- fit_p_b1$datt[fit_p_b1$datt$period == 2018, ]
datt_2019 <- fit_p_b1$datt[fit_p_b1$datt$period == 2019, ]
t_2018    <- unique(datt_2018$t_index)
t_2019    <- unique(datt_2019$t_index)
width_2018 <- datt_2018$datt_upper[1] - datt_2018$datt_lower[1]
width_2019 <- datt_2019$datt_upper[1] - datt_2019$datt_lower[1]
stopifnot(abs(width_2018 - 2 * t_2018 * 0.5) < 1e-10)
stopifnot(abs(width_2019 - 2 * t_2019 * 0.5) < 1e-10)
stopifnot(width_2019 > width_2018)
cat(sprintf("  LPT-P, B=0.5: width grows. t(2018)=%d->w=%.1f, t(2019)=%d->w=%.1f. OK\n",
            t_2018, width_2018, t_2019, width_2019))

fit_p_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "P")
stopifnot(fit_p_cal$calibration$type == "first_diff")
cat(sprintf("  LPT-P, calibrated: B_hat=%.4f (first_diff). OK\n", fit_p_cal$B_hat))

att_2018 <- fit_p_b1$att[fit_p_b1$att$period == 2018, ]
att_2019 <- fit_p_b1$att[fit_p_b1$att$period == 2019, ]
expected_ratio <- t_2019 / t_2018
actual_ratio   <- (att_2019$att_upper[25] - att_2019$att_lower[25]) /
                  (att_2018$att_upper[25] - att_2018$att_lower[25])
stopifnot(abs(actual_ratio - expected_ratio) < 1e-6)
cat(sprintf("  LPT-P: ATT width ratio = %.2f (expected = %.2f). OK\n",
            actual_ratio, expected_ratio))

atto_2018 <- fit_p_b1$att_o[fit_p_b1$att_o$period == 2018, ]
atto_2019 <- fit_p_b1$att_o[fit_p_b1$att_o$period == 2019, ]
stopifnot((atto_2019$att_o_upper - atto_2019$att_o_lower) >
          (atto_2018$att_o_upper - atto_2018$att_o_lower))
cat("  LPT-P: ATT^o width grows with t. OK\n")

cat("lpt_type='P' tests PASSED.\n")
```

```bash
git rm _local/TEST_MULT/tests/test-lpt-b.R
git add _local/TEST_MULT/tests/test-lpt-p.R
git commit -m "test: rename test-lpt-b -> test-lpt-p, update type labels"
```

---

## Task 7: Update test-backwards-compat.R and test-summary-plot.R

**Files:**
- Modify: `_local/TEST_MULT/tests/test-backwards-compat.R`
- Modify: `_local/TEST_MULT/tests/test-summary-plot.R`

### test-backwards-compat.R

Replace the full file:

```r
# test-backwards-compat.R
cat("Testing backwards compatibility...\n")

# --- Test 1: B=0 produces same results for C and P ---
fit_c <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0, lpt_type = "C")
fit_p <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0, lpt_type = "P")
stopifnot(all(abs(fit_c$datt$datt_lower - fit_p$datt$datt_lower) < 1e-10))
stopifnot(all(abs(fit_c$datt$datt_upper - fit_p$datt$datt_upper) < 1e-10))
cat("  B=0: LPT-C == LPT-P. OK\n")

# --- Test 2: Default lpt_type is "C" ---
fit_default <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = 2019, pre_periods = 1993:1999, B = 0)
stopifnot(fit_default$lpt_type == "C")
cat("  Default lpt_type is 'C'. OK\n")

# --- Test 3: t_index >= 1 ---
t_val <- fit_c$t_index_map[["2019"]]
stopifnot(t_val >= 1)
cat(sprintf("  t_index for 2019 (ref=1999): %d. OK\n", t_val))

# --- Test 4: Numeric B works with C and P ---
fit_c_num <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = 0.1, lpt_type = "C")
fit_p_num <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = 0.1, lpt_type = "P")
stopifnot(fit_c_num$B_hat == 0.1)
stopifnot(fit_p_num$B_hat == 0.1)
cat("  Numeric B works with C and P. OK\n")

# --- Test 5: lpt object has all expected fields ---
expected_fields <- c("datt", "att", "att_o", "B_hat", "B_values",
                      "b0", "C_hat", "calibration", "slopes", "call",
                      "n", "has_untreated", "post_periods", "ref_period",
                      "lpt_type", "t_index_map", "specifications")
missing <- setdiff(expected_fields, names(fit_c))
stopifnot(length(missing) == 0)
cat("  All expected fields present. OK\n")

# --- Test 6: b0 and C_hat are NULL for C and P ---
stopifnot(is.null(fit_c$b0))
stopifnot(is.null(fit_c$C_hat))
stopifnot(is.null(fit_p$b0))
stopifnot(is.null(fit_p$C_hat))
cat("  b0/C_hat NULL for C and P. OK\n")

# --- Test 7: datt has t_index column ---
stopifnot("t_index" %in% names(fit_c$datt))
cat("  datt has t_index column. OK\n")

cat("Backwards compatibility tests PASSED.\n")
```

### test-summary-plot.R

Replace the full file:

```r
# test-summary-plot.R
cat("Testing summary and plot methods...\n")

fit_c <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "C")
fit_p <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "P")
fit_s <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "S")

cat("\n--- print LPT-C ---\n"); print(fit_c)
cat("\n--- print LPT-P ---\n"); print(fit_p)
cat("\n--- print LPT-S ---\n"); print(fit_s)
cat("  print works for all types. OK\n")

cat("\n--- summary LPT-C ---\n"); summary(fit_c)
cat("\n--- summary LPT-P ---\n"); summary(fit_p)
cat("\n--- summary LPT-S ---\n"); summary(fit_s)
cat("  summary works for all types. OK\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p1 <- plot(fit_c, type = "datt"); p2 <- plot(fit_p, type = "datt")
  p3 <- plot(fit_s, type = "datt")
  cat("  plot(type='datt') works for all. OK\n")

  p4 <- plot(fit_c, type = "att"); p5 <- plot(fit_p, type = "att")
  p6 <- plot(fit_s, type = "att")
  cat("  plot(type='att') works for all. OK\n")

  p7 <- plot(fit_c, type = "pretrends"); p8 <- plot(fit_p, type = "pretrends")
  p9 <- plot(fit_s, type = "pretrends")
  cat("  plot(type='pretrends') works for all. OK\n")

  p10 <- plot(fit_c, type = "sensitivity")
  p11 <- plot(fit_p, type = "sensitivity")
  p12 <- plot(fit_s, type = "sensitivity")
  cat("  plot(type='sensitivity') works for all. OK\n")
}

cat("Summary and plot tests PASSED.\n")
```

Run both tests:
```bash
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-backwards-compat.R', local=TRUE)"
Rscript -e "source('_local/TEST_MULT/source_all.R'); source('_local/TEST_MULT/tests/test-summary-plot.R', local=TRUE)"
```

Commit:
```bash
git add _local/TEST_MULT/tests/test-backwards-compat.R _local/TEST_MULT/tests/test-summary-plot.R
git commit -m "test: update backwards-compat and summary-plot tests for C/P/S types"
```

---

## Task 8: Run Full Test Suite

```bash
Rscript -e "source('_local/TEST_MULT/run_tests.R')"
```

Expected output — all six test files run with PASSED:
```
=== Running test-backwards-compat.R ===
Backwards compatibility tests PASSED.
=== Running test-calibrate-cumulative.R ===
calibrate_B cumulative tests PASSED.
calibrate_B smooth tests PASSED.
=== Running test-lpt-c.R ===
lpt_type='C' tests PASSED.
=== Running test-lpt-p.R ===
lpt_type='P' tests PASSED.
=== Running test-lpt-s.R ===
lpt_type='S' tests PASSED.
=== Running test-summary-plot.R ===
Summary and plot tests PASSED.
=== All test files executed ===
```

If any test fails, debug the implementation — do NOT weaken the test.

Final commit:
```bash
git add -A
git commit -m "test: all six test files pass for LPT-C/P/S"
```

---

## Verification Checklist

- [ ] `Rscript -e "source('_local/TEST_MULT/run_tests.R')"` — all tests pass
- [ ] LPT-C with B=0.5, two periods: dATT width is 1.0 for BOTH periods
- [ ] LPT-P with B=0.5, two periods: dATT width grows (2*t*0.5)
- [ ] LPT-S with b0=0.1, C=0.2, two periods: width = 2*Phi(t) = 2*(t*0.1 + 0.2*t*(t+1)/2)
- [ ] `summary(lpt(..., lpt_type="S"))` prints "LPT-S" and shows Phi(t) per period
- [ ] `plot(fit_s, type="sensitivity")` x-axis label is C/Ĉ
