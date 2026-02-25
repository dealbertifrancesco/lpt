test_that("lpt works on bundled sru data with B=0", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999, B = 0)

  expect_s3_class(fit, "lpt")
  expect_true(is.data.frame(fit$datt))
  expect_true(is.data.frame(fit$att))
  expect_equal(fit$B_hat, 0)
  expect_true(all(c("period", "d", "lambda_d", "datt_lower", "datt_upper")
                    %in% names(fit$datt)))
  expect_true(all(c("period", "d", "Lambda_d", "att_lower", "att_upper")
                    %in% names(fit$att)))
})

test_that("lpt works with B=calibrate on sru data", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = "calibrate")

  expect_s3_class(fit, "lpt")
  expect_gt(fit$B_hat, 0)
  expect_true(!is.null(fit$calibration))
  # Identified set should be wider than point estimate
  datt <- fit$datt[fit$datt$B == fit$B_hat, ]
  expect_true(all(datt$datt_upper - datt$datt_lower > 0))
})

test_that("lpt computes ATT^o (Corollary 1)", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999, B = 0)

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
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = 0)

  expect_equal(length(fit$post_periods), 2)
  expect_true(all(c(2018, 2019) %in% fit$datt$period))
})

test_that("print and summary methods work", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999, B = 0)

  expect_output(print(fit), "lpt")
  expect_output(summary(fit), "Lipschitz Parallel Trends")
})

test_that("plot methods do not error", {
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = "calibrate")

  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "pretrends"))
  expect_no_error(plot(fit, type = "sensitivity"))
})
