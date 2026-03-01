source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("plot.lpt works with time_restriction='rm'", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_lpt(n = 300, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 60)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = attr(sim, "true_B"),
             time_restriction = "rm", M_bar = 1)

  expect_no_error(plot(fit, type = "att"))
  expect_no_error(plot(fit, type = "datt"))
  expect_no_error(plot(fit, type = "sensitivity", estimand = "att_o"))
})

test_that("sensitivity plot with rm uses M_bar axis", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_lpt(n = 300, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 61)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = 4L, pre_periods = 1:3,
             B = 0,
             time_restriction = "rm", M_bar = 1)

  p <- plot(fit, type = "sensitivity", estimand = "att_o")
  # The x-axis label should reference M_bar or M[bar]
  x_label <- as.character(p$labels$x)
  expect_true(grepl("M_bar|M\\[bar\\]", x_label, ignore.case = FALSE) ||
              grepl("M_bar", x_label, fixed = TRUE),
              label = paste("x-axis references M_bar, got:", x_label))
})
