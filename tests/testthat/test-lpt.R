test_that("lpt works on sru data with B=0", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, B = 0)

  expect_s3_class(fit, "lpt")
  expect_true(is.data.frame(fit$datt))
  expect_true(is.data.frame(fit$att))
  expect_equal(fit$B_hat, 0)
  expect_true(all(c("period", "horizon", "d", "lambda_d",
                     "datt_lower", "datt_upper") %in% names(fit$datt)))
  expect_true(all(c("period", "horizon", "d", "Lambda_d",
                     "att_lower", "att_upper") %in% names(fit$att)))
  # With B=0, IS width should be 0
  expect_true(all(abs(fit$datt$datt_upper - fit$datt$datt_lower) < 1e-10))
})

test_that("lpt works with B=calibrate on sru data", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             B = "calibrate")

  expect_s3_class(fit, "lpt")
  expect_gt(fit$B_hat, 0)
  expect_true(!is.null(fit$calibration))
  datt <- fit$datt[fit$datt$B == fit$B_hat, ]
  expect_true(all(datt$datt_upper - datt$datt_lower > 0))
})

test_that("lpt computes ATT^o", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, B = 0)

  expect_true(!is.null(fit$att_o))
  expect_true(is.data.frame(fit$att_o))
  expect_true(all(c("att_o_bin", "D_bar", "att_o_lower", "att_o_upper")
                    %in% names(fit$att_o)))
  # With B=0, bounds should equal att_o_bin
  expect_equal(fit$att_o$att_o_lower, fit$att_o$att_o_bin)
  expect_equal(fit$att_o$att_o_upper, fit$att_o$att_o_bin)
})

test_that("lpt works with multiple post-periods", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(4, 5), pre_periods = -7:-1,
             B = 0)

  expect_equal(length(fit$post_periods), 2)
  expect_true(all(c(4, 5) %in% fit$datt$period))
})

test_that("print and summary methods work", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, B = 0)

  expect_output(print(fit), "lpt")
  expect_output(summary(fit), "Local Parallel Trends")
})

test_that("plot methods do not error", {
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             B = "calibrate")

  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "pretrends"))
  expect_no_error(plot(fit, type = "sensitivity"))
})

test_that("multi-period IS width scales with (t+1)", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1, B = 0.1)

  # At horizon 0 (period 0): IS width = 2 * (0+1) * 0.1 = 0.2
  datt_h0 <- fit$datt[fit$datt$period == 0 & fit$datt$B == 0.1, ]
  widths_h0 <- datt_h0$datt_upper - datt_h0$datt_lower
  expect_equal(widths_h0[1], 2 * 1 * 0.1, tolerance = 1e-10)

  # At horizon 5 (period 5): IS width = 2 * (5+1) * 0.1 = 1.2
  datt_h5 <- fit$datt[fit$datt$period == 5 & fit$datt$B == 0.1, ]
  widths_h5 <- datt_h5$datt_upper - datt_h5$datt_lower
  expect_equal(widths_h5[1], 2 * 6 * 0.1, tolerance = 1e-10)
})

test_that("att_o_agg is computed for multiple post-periods", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1, B = 0.1)

  expect_true(!is.null(fit$att_o_agg))
  expect_true(is.data.frame(fit$att_o_agg))
  expect_true(all(c("Lambda_agg", "D_bar", "n_periods",
                     "att_o_agg_lower", "att_o_agg_upper")
                    %in% names(fit$att_o_agg)))
  expect_equal(fit$att_o_agg$n_periods[1], 6)
})

test_that("att_o_agg is NULL for single post-period", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, B = 0.1)

  expect_null(fit$att_o_agg)
})

test_that("horizon column exists in output data frames", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(0, 5), pre_periods = -7:-1, B = 0)

  expect_true("horizon" %in% names(fit$datt))
  expect_true("horizon" %in% names(fit$att))
  expect_true("horizon" %in% names(fit$att_o))
  expect_equal(fit$datt$horizon[fit$datt$period == 0][1], 0L)
  expect_equal(fit$datt$horizon[fit$datt$period == 5][1], 5L)
})
