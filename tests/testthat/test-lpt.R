test_that("lpt works on sru data with M=0, B=0 (point identification)", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, M = 0, B = 0)

  expect_s3_class(fit, "lpt")
  expect_true(is.data.frame(fit$datt))
  expect_true(is.data.frame(fit$att))
  expect_equal(fit$B_hat, 0)
  expect_equal(fit$M_hat, 0)
  expect_true(all(c("period", "horizon", "d", "lambda_d", "B",
                     "datt_lower", "datt_upper") %in% names(fit$datt)))
  expect_true(all(c("period", "horizon", "d", "Lambda_d", "M",
                     "att_lower", "att_upper") %in% names(fit$att)))
  # With B=0, dATT IS width should be 0
  expect_true(all(abs(fit$datt$datt_upper - fit$datt$datt_lower) < 1e-10))
  # With M=0, ATT IS width should be 0
  expect_true(all(abs(fit$att$att_upper - fit$att$att_lower) < 1e-10))
})

test_that("lpt calibrates M and B on sru data", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             M = "calibrate", B = "calibrate")

  expect_s3_class(fit, "lpt")
  expect_gt(fit$B_hat, 0)
  expect_gt(fit$M_hat, 0)
  expect_equal(fit$M_source, "calibrated")
  expect_equal(fit$B_source, "calibrated")
  expect_true(!is.null(fit$calibration))
  datt <- fit$datt[fit$datt$B == fit$B_hat, ]
  expect_true(all(datt$datt_upper - datt$datt_lower > 0))
  att <- fit$att[fit$att$M == fit$M_hat, ]
  expect_true(all(att$att_upper - att$att_lower > 0))
})

test_that("ATT identified-set width is constant across doses", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, M = 0.2, B = 0)

  # IS_ATT(d, t; M) = Lambda_t(d) -/+ (t+1)M: width 2*(t+1)*M at every dose
  widths <- fit$att$att_upper - fit$att$att_lower
  h <- fit$att$horizon[1]
  expect_true(all(abs(widths - 2 * (h + 1) * 0.2) < 1e-10))
})

test_that("lpt computes ATT^o with M-based bounds", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, M = 0, B = 0)

  expect_true(!is.null(fit$att_o))
  expect_true(is.data.frame(fit$att_o))
  expect_true(all(c("att_o_bin", "M", "att_o_lower", "att_o_upper")
                    %in% names(fit$att_o)))
  # With M=0, bounds should equal att_o_bin
  expect_equal(fit$att_o$att_o_lower, fit$att_o$att_o_bin)
  expect_equal(fit$att_o$att_o_upper, fit$att_o$att_o_bin)
})

test_that("lpt works with multiple post-periods", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(4, 5), pre_periods = -7:-1,
             M = 0, B = 0)

  expect_equal(length(fit$post_periods), 2)
  expect_true(all(c(4, 5) %in% fit$datt$period))
})

test_that("print and summary methods work", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, M = 0, B = 0)

  expect_output(print(fit), "lpt")
  expect_output(print(fit), "M = ")
  expect_output(summary(fit), "Local Parallel Trends")
})

test_that("plot methods do not error", {
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             M = "calibrate", B = "calibrate")

  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "pretrends"))
  expect_no_error(plot(fit, type = "sensitivity"))
  expect_no_error(plot(fit, type = "sensitivity", estimand = "att", d0 = 0.3))
  expect_no_error(plot(fit, type = "sensitivity", estimand = "datt", d0 = 0.3))
  expect_no_error(plot(fit, type = "eventstudy"))
})

test_that("pre_att_o is computed for event study", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1, M = 0, B = 0)

  expect_true(!is.null(fit$pre_att_o))
  expect_true(is.data.frame(fit$pre_att_o))
  expect_true(all(c("period", "att_o_bin") %in% names(fit$pre_att_o)))
  # Should have entries for pre-periods -7 to -2 (not ref period -1)
  expect_equal(nrow(fit$pre_att_o), 6)
  expect_true(all(fit$pre_att_o$period %in% -7:-2))
  expect_false(-1 %in% fit$pre_att_o$period)
})

test_that("eventstudy plot returns ggplot object", {
  skip_if_not_installed("ggplot2")

  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1,
             M = "calibrate", B = "calibrate")

  p <- plot(fit, type = "eventstudy")
  expect_s3_class(p, "ggplot")
})

test_that("multi-period IS widths scale with (t+1)", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1, M = 0.2, B = 0.1)

  # dATT at horizon 0: width = 2 * (0+1) * B
  datt_h0 <- fit$datt[fit$datt$period == 0 & fit$datt$B == 0.1, ]
  expect_equal((datt_h0$datt_upper - datt_h0$datt_lower)[1],
               2 * 1 * 0.1, tolerance = 1e-10)

  # dATT at horizon 5: width = 2 * (5+1) * B
  datt_h5 <- fit$datt[fit$datt$period == 5 & fit$datt$B == 0.1, ]
  expect_equal((datt_h5$datt_upper - datt_h5$datt_lower)[1],
               2 * 6 * 0.1, tolerance = 1e-10)

  # ATT at horizon 0: width = 2 * (0+1) * M
  att_h0 <- fit$att[fit$att$period == 0 & fit$att$M == 0.2, ]
  expect_equal((att_h0$att_upper - att_h0$att_lower)[1],
               2 * 1 * 0.2, tolerance = 1e-10)

  # ATT at horizon 5: width = 2 * (5+1) * M
  att_h5 <- fit$att[fit$att$period == 5 & fit$att$M == 0.2, ]
  expect_equal((att_h5$att_upper - att_h5$att_lower)[1],
               2 * 6 * 0.2, tolerance = 1e-10)

  # ATT^o at horizon 5: width = 2 * (5+1) * M
  atto_h5 <- fit$att_o[fit$att_o$period == 5 & fit$att_o$M == 0.2, ]
  expect_equal(atto_h5$att_o_upper - atto_h5$att_o_lower,
               2 * 6 * 0.2, tolerance = 1e-10)
})

test_that("att_o_agg uses the (T+2)/2 * M half-width", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 0:5, pre_periods = -7:-1, M = 0.2, B = 0)

  expect_true(!is.null(fit$att_o_agg))
  expect_true(is.data.frame(fit$att_o_agg))
  expect_true(all(c("Lambda_agg", "n_periods", "M",
                     "att_o_agg_lower", "att_o_agg_upper")
                    %in% names(fit$att_o_agg)))
  expect_equal(fit$att_o_agg$n_periods[1], 6)

  # For contiguous horizons 0..5 (T = 5): half-width = (T+2)/2 * M = 3.5 * M
  width <- fit$att_o_agg$att_o_agg_upper - fit$att_o_agg$att_o_agg_lower
  expect_equal(width, 2 * 3.5 * 0.2, tolerance = 1e-10)

  # Center is the time average of the per-period estimands
  expect_equal(fit$att_o_agg$Lambda_agg, mean(fit$att_o$att_o_bin),
               tolerance = 1e-10)
})

test_that("att_o_agg is NULL for single post-period", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, M = 0.2, B = 0.1)

  expect_null(fit$att_o_agg)
})

test_that("M and B sensitivity vectors are supported", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             M = c(0, 0.1, 0.2), B = c(0, 0.05))

  expect_equal(sort(unique(fit$datt$B)), c(0, 0.05))
  expect_equal(sort(unique(fit$att$M)), c(0, 0.1, 0.2))
  expect_equal(fit$M_hat, 0.2)
  expect_equal(fit$B_hat, 0.05)
})

test_that("mixed specification works: numeric B with calibrated M", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1,
             M = "calibrate", B = 0)

  expect_equal(fit$B_hat, 0)
  expect_equal(fit$B_source, "user-supplied")
  expect_gt(fit$M_hat, 0)
  expect_equal(fit$M_source, "calibrated")
  expect_true(!is.null(fit$calibration))
})

test_that("horizon column exists in output data frames", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(0, 5), pre_periods = -7:-1, M = 0, B = 0)

  expect_true("horizon" %in% names(fit$datt))
  expect_true("horizon" %in% names(fit$att))
  expect_true("horizon" %in% names(fit$att_o))
  expect_equal(fit$datt$horizon[fit$datt$period == 0][1], 0L)
  expect_equal(fit$datt$horizon[fit$datt$period == 5][1], 5L)
})

test_that("default evaluation grid lies in the treated dose support", {
  data(sru, package = "lpt")
  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 5, pre_periods = -7:-1, M = 0, B = 0)

  d_pos <- sru$dose[sru$year == -1 & sru$dose > 0]
  expect_true(all(fit$datt$d >= stats::quantile(d_pos, 0.05) - 1e-10))
  expect_true(all(fit$datt$d <= stats::quantile(d_pos, 0.95) + 1e-10))
  expect_true(all(fit$datt$d > 0))
})
