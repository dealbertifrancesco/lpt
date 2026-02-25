# Helper to get first differences from simulated data
get_post_diffs <- function(dat) {
  n_pre <- max(dat$time) - 1
  pre <- dat[dat$time == n_pre, ]
  post <- dat[dat$time == n_pre + 1, ]
  merged <- merge(pre[, c("id", "outcome", "dose")],
                  post[, c("id", "outcome")],
                  by = "id", suffixes = c("_pre", "_post"))
  merged$delta_y <- merged$outcome_post - merged$outcome_pre
  merged
}

test_that("estimate_dose_slope recovers truth under no selection", {
  dat <- simulate_lpt(n = 5000, n_pre = 1, dgp = "no_selection", seed = 42)
  m <- get_post_diffs(dat)
  result <- estimate_dose_slope(m$delta_y, m$dose)

  true_datt <- attr(dat, "true_datt")
  true_lambda <- true_datt(result$eval_points)

  mae <- mean(abs(result$lambda_d - true_lambda))
  expect_lt(mae, 0.2)
})

test_that("estimate_dose_slope returns correct structure", {
  dat <- simulate_lpt(n = 500, dgp = "no_selection", seed = 1)
  m <- get_post_diffs(dat)
  result <- estimate_dose_slope(m$delta_y, m$dose)

  expect_true(is.list(result))
  expect_s3_class(result, "dose_slope")
  expect_true(all(c("eval_points", "conditional_mean", "lambda_d",
                     "se_lambda", "gam_fit") %in% names(result)))
  expect_equal(length(result$eval_points), length(result$lambda_d))
  expect_true(all(result$se_lambda >= 0))
})

test_that("k adjustment warning fires for small n", {
  expect_warning(
    estimate_dose_slope(rnorm(30), runif(30), k = 20),
    "too large"
  )
})

test_that("estimate_dose_slope input validation", {
  expect_error(estimate_dose_slope("a", 1:10), "numeric")
  expect_error(estimate_dose_slope(1:5, 1:10), "same length")
  expect_error(estimate_dose_slope(1:5, 1:5), "at least 10")
})
