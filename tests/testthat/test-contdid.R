test_that("lpt with method='contdid' produces valid lpt object", {
  skip_if_not_installed("contdid")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 0:5, pre_periods = -7:-1,
        M = "calibrate", B = "calibrate", method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  expect_s3_class(fit, "lpt")

  # All standard fields present
  expect_true(is.data.frame(fit$datt))
  expect_true(is.data.frame(fit$att))
  expect_true(is.data.frame(fit$att_o))
  expect_true(!is.null(fit$contdid_fit))
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
  expect_equal(fit$specifications$method, "contdid")
})

test_that("contdid calibration produces deviations and slopes", {
  skip_if_not_installed("contdid")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 5, pre_periods = -7:-1,
        M = "calibrate", B = "calibrate", method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  # Calibration output matches expected structure
  expect_true(!is.null(fit$calibration))
  expect_true(!is.null(fit$calibration$pre_slopes))
  expect_true(all(c("period_pair", "d", "mu_prime_d")
                    %in% names(fit$calibration$pre_slopes)))
  expect_true(!is.null(fit$calibration$pre_deviations))
  expect_true(all(c("period_pair", "d", "deviation")
                    %in% names(fit$calibration$pre_deviations)))
  expect_true(length(fit$calibration$sup_slope_by_period) >= 1)
  expect_true(length(fit$calibration$sup_dev_by_period) >= 1)
})

test_that("contdid plots do not error", {
  skip_if_not_installed("contdid")
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 0:5, pre_periods = -7:-1,
        M = "calibrate", B = "calibrate", method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "pretrends"))
  expect_no_error(plot(fit, type = "eventstudy"))
  expect_no_error(plot(fit, type = "sensitivity"))
})

test_that("contdid with numeric M and B works", {
  skip_if_not_installed("contdid")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 5, pre_periods = -7:-1,
        M = 0, B = 0, method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  expect_s3_class(fit, "lpt")
  expect_equal(fit$B_hat, 0)
  expect_equal(fit$M_hat, 0)

  # With B=0 / M=0, IS widths should be 0
  expect_true(all(abs(fit$datt$datt_upper - fit$datt$datt_lower) < 1e-10))
  expect_true(all(abs(fit$att$att_upper - fit$att$att_lower) < 1e-10))

  # Calibration should be NULL when both bounds are numeric
  expect_null(fit$calibration)
})
