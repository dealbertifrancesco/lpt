test_that("estimate_dose_slope returns correct structure on sru data", {
  data(sru, package = "lpt")
  ref <- sru[sru$year == 1999, ]
  post <- sru[sru$year == 2019, ]
  m <- merge(ref[, c("commune", "outcome", "dose")],
             post[, c("commune", "outcome")],
             by = "commune", suffixes = c("_ref", "_post"))
  m$delta_y <- m$outcome_post - m$outcome_ref

  result <- estimate_dose_slope(m$delta_y, m$dose)

  expect_true(is.list(result))
  expect_s3_class(result, "dose_slope")
  expect_true(all(c("eval_points", "conditional_mean", "lambda_d",
                     "gam_fit") %in% names(result)))
  expect_equal(length(result$eval_points), length(result$lambda_d))
  expect_equal(length(result$eval_points), 50)
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
