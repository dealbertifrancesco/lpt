test_that("lpt with method='npiv' produces valid lpt object", {
  skip_if_not_installed("npiv")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1,
             B = "calibrate", method = "npiv")

  expect_s3_class(fit, "lpt")

  # All standard fields present
  expect_true(is.data.frame(fit$datt))
  expect_true(is.data.frame(fit$att))
  expect_true(is.data.frame(fit$att_o))
  expect_true(!is.null(fit$npiv_fits))
  expect_null(fit$contdid_fit)
  expect_gt(fit$B_hat, 0)

  # datt has expected columns
  expect_true(all(c("period", "horizon", "d", "lambda_d", "B",
                     "datt_lower", "datt_upper") %in% names(fit$datt)))

  # att has expected columns
  expect_true(all(c("period", "horizon", "d", "Lambda_d", "B",
                     "att_lower", "att_upper") %in% names(fit$att)))

  # att_o has expected columns
  expect_true(all(c("period", "horizon", "att_o_bin", "D_bar", "B",
                     "att_o_lower", "att_o_upper") %in% names(fit$att_o)))

  # Identified sets widen with B
  datt_B <- fit$datt[fit$datt$period == 0, ]
  widths <- datt_B$datt_upper - datt_B$datt_lower
  expect_true(all(widths > 0))

  # Multiple post-periods present
  expect_equal(length(unique(fit$datt$period)), 6)

  # att_o_agg computed for multiple periods
  expect_true(!is.null(fit$att_o_agg))

  # specifications record method
  expect_equal(fit$specifications$method, "npiv")
})

test_that("npiv calibration produces pre_slopes for pretrends plot", {
  skip_if_not_installed("npiv")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             B = "calibrate", method = "npiv")

  # Calibration output matches expected structure
  expect_true(!is.null(fit$calibration))
  expect_true(!is.null(fit$calibration$pre_slopes))
  expect_true(all(c("period_pair", "d", "mu_prime_d")
                    %in% names(fit$calibration$pre_slopes)))
  expect_equal(length(fit$calibration$sup_by_period), 6)
})

test_that("npiv plots do not error", {
  skip_if_not_installed("npiv")
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1,
             B = "calibrate", method = "npiv")

  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "pretrends"))
  expect_no_error(plot(fit, type = "eventstudy"))
  expect_no_error(plot(fit, type = "sensitivity"))
})

test_that("npiv with numeric B works", {
  skip_if_not_installed("npiv")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             B = 0, method = "npiv")

  expect_s3_class(fit, "lpt")
  expect_equal(fit$B_hat, 0)

  # With B=0, IS width should be 0
  expect_true(all(abs(fit$datt$datt_upper - fit$datt$datt_lower) < 1e-10))

  # Calibration should be NULL
  expect_null(fit$calibration)
})

test_that("npiv_args are forwarded correctly", {
  skip_if_not_installed("npiv")

  data(sru, package = "lpt")

  # Custom J.x.segments
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0, pre_periods = -7:-1,
             B = "calibrate", method = "npiv",
             npiv_args = list(J.x.segments = 3))

  expect_s3_class(fit, "lpt")

  # The npiv fit should have the custom segments
  npiv_obj <- fit$npiv_fits[["0"]]
  expect_equal(npiv_obj$J.x.segments, 3)
})

test_that("npiv slopes are compatible with dose_slope class", {
  skip_if_not_installed("npiv")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0, pre_periods = -7:-1,
             B = "calibrate", method = "npiv")

  # Slopes list should have one entry per post-period
  expect_equal(length(fit$slopes), 1)
  sl <- fit$slopes[["0"]]
  expect_s3_class(sl, "dose_slope")
  expect_true(is.numeric(sl$eval_points))
  expect_true(is.numeric(sl$lambda_d))
  expect_true(is.numeric(sl$conditional_mean))
  expect_null(sl$gam_fit)  # npiv doesn't use gam
  expect_equal(length(sl$eval_points), length(sl$lambda_d))
})
