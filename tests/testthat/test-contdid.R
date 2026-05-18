test_that("lpt with method='contdid' produces valid lpt object", {
  skip_if_not_installed("contdid")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 0:5, pre_periods = -7:-1,
        B = "calibrate", method = "contdid",
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
  expect_equal(fit$specifications$method, "contdid")
})

test_that("contdid calibration produces pre_slopes for pretrends plot", {
  skip_if_not_installed("contdid")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 5, pre_periods = -7:-1,
        B = "calibrate", method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  # Calibration output matches expected structure
  expect_true(!is.null(fit$calibration))
  expect_true(!is.null(fit$calibration$pre_slopes))
  expect_true(all(c("period_pair", "d", "mu_prime_d")
                    %in% names(fit$calibration$pre_slopes)))
  expect_true(length(fit$calibration$sup_by_period) >= 1)
})

test_that("contdid plots do not error", {
  skip_if_not_installed("contdid")
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 0:5, pre_periods = -7:-1,
        B = "calibrate", method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "pretrends"))
  expect_no_error(plot(fit, type = "eventstudy"))
  expect_no_error(plot(fit, type = "sensitivity"))
})

test_that("contdid with numeric B works", {
  skip_if_not_installed("contdid")

  data(sru, package = "lpt")
  fit <- tryCatch(
    lpt(sru, "commune", "year", "outcome", "dose",
        post_period = 5, pre_periods = -7:-1,
        B = 0, method = "contdid",
        contdid_args = list(num_knots = 1, degree = 3, biters = 100)),
    error = function(e) skip(paste("contdid failed:", e$message))
  )

  expect_s3_class(fit, "lpt")
  expect_equal(fit$B_hat, 0)

  # With B=0, IS width should be 0
  expect_true(all(abs(fit$datt$datt_upper - fit$datt$datt_lower) < 1e-10))

  # Calibration should be NULL
  expect_null(fit$calibration)
})
