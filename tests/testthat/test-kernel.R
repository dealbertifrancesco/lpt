test_that("lpt with method='kernel' produces valid lpt object", {
  skip_if_not_installed("np")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1,
             M = "calibrate", B = "calibrate", method = "kernel")

  expect_s3_class(fit, "lpt")

  # All standard fields present
  expect_true(is.data.frame(fit$datt))
  expect_true(is.data.frame(fit$att))
  expect_true(is.data.frame(fit$att_o))
  expect_true(!is.null(fit$kernel_fits))
  expect_null(fit$contdid_fit)
  expect_null(fit$npiv_fits)
  expect_gt(fit$B_hat, 0)
  expect_gt(fit$M_hat, 0)

  # datt has expected columns
  expect_true(all(c("period", "horizon", "d", "lambda_d", "B",
                     "datt_lower", "datt_upper") %in% names(fit$datt)))

  # att has expected columns (bounded by M)
  expect_true(all(c("period", "horizon", "d", "Lambda_d", "M",
                     "att_lower", "att_upper") %in% names(fit$att)))

  # att_o has expected columns (bounded by M)
  expect_true(all(c("period", "horizon", "att_o_bin", "M",
                     "att_o_lower", "att_o_upper") %in% names(fit$att_o)))

  # Identified sets widen with B
  datt_B <- fit$datt[fit$datt$period == 0, ]
  widths <- datt_B$datt_upper - datt_B$datt_lower
  expect_true(all(widths > 0))

  # ATT IS width is constant across doses: 2*(t+1)*M
  att_h0 <- fit$att[fit$att$period == 0, ]
  att_widths <- att_h0$att_upper - att_h0$att_lower
  expect_true(all(abs(att_widths - 2 * fit$M_hat) < 1e-10))

  # Multiple post-periods present
  expect_equal(length(unique(fit$datt$period)), 6)

  # att_o_agg computed for multiple periods
  expect_true(!is.null(fit$att_o_agg))

  # specifications record method
  expect_equal(fit$specifications$method, "kernel")
})

test_that("kernel calibration produces deviations and slopes", {
  skip_if_not_installed("np")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0, pre_periods = -7:-1,
             M = "calibrate", B = "calibrate", method = "kernel")

  # Calibration should exist
  expect_true(!is.null(fit$calibration))
  expect_true(!is.null(fit$calibration$pre_slopes))
  expect_true(!is.null(fit$calibration$pre_deviations))

  # pre_slopes / pre_deviations have expected structure
  ps <- fit$calibration$pre_slopes
  expect_true(all(c("period_pair", "d", "mu_prime_d") %in% names(ps)))
  pd <- fit$calibration$pre_deviations
  expect_true(all(c("period_pair", "d", "deviation") %in% names(pd)))

  # 6 consecutive pairs from periods -7 to -1
  expect_equal(length(unique(ps$period_pair)), 6)

  # sup vectors have named entries
  sup_b <- fit$calibration$sup_slope_by_period
  sup_m <- fit$calibration$sup_dev_by_period
  expect_equal(length(sup_b), 6)
  expect_equal(length(sup_m), 6)
  expect_true(all(nchar(names(sup_b)) > 0))

  # B_hat / M_hat = max of sup values
  expect_equal(fit$B_hat, max(sup_b))
  expect_equal(fit$M_hat, max(sup_m))
})

test_that("all plot types work with kernel method", {
  skip_if_not_installed("np")
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:2, pre_periods = -7:-1,
             M = "calibrate", B = "calibrate", method = "kernel")

  for (ptype in c("datt", "att", "pretrends", "eventstudy", "sensitivity")) {
    p <- plot(fit, type = ptype)
    expect_s3_class(p, "gg")
  }
})

test_that("kernel with M = 0 and B = 0 gives point identification", {
  skip_if_not_installed("np")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0, pre_periods = -7:-1,
             M = 0, B = 0, method = "kernel")

  # Point identification: IS width = 0
  expect_true(all(fit$datt$datt_lower == fit$datt$datt_upper))
  expect_true(all(fit$att$att_lower == fit$att$att_upper))
  expect_equal(fit$B_hat, 0)
  expect_equal(fit$M_hat, 0)

  # Calibration should be NULL when both bounds are numeric
  expect_null(fit$calibration)
})

test_that("kernel_args forwarding works", {
  skip_if_not_installed("np")

  data(sru, package = "lpt")

  # Manual bandwidth
  fit_manual <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = 0, pre_periods = -7:-1,
                    M = 0, B = 0, method = "kernel",
                    kernel_args = list(bw = 0.1))
  expect_s3_class(fit_manual, "lpt")

  # Check that the specified bandwidth was used
  bws <- sapply(fit_manual$kernel_fits, function(x) x$bw)
  expect_true(all(abs(bws - 0.1) < 1e-10))

  # Epanechnikov kernel
  fit_ep <- lpt(sru, "commune", "year", "outcome", "dose",
                post_period = 0, pre_periods = -7:-1,
                M = 0, B = 0, method = "kernel",
                kernel_args = list(bw = 0.1, ckertype = "epanechnikov"))
  expect_s3_class(fit_ep, "lpt")
})

test_that("kernel slopes are dose_slope compatible", {
  skip_if_not_installed("np")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0, pre_periods = -7:-1,
             M = 0, B = 0, method = "kernel")

  # First slope should be a dose_slope object
  slope1 <- fit$slopes[[1]]
  expect_s3_class(slope1, "dose_slope")
  expect_true(!is.null(slope1$eval_points))
  expect_true(!is.null(slope1$conditional_mean))
  expect_true(!is.null(slope1$lambda_d))
  expect_equal(length(slope1$eval_points), length(slope1$lambda_d))
})
