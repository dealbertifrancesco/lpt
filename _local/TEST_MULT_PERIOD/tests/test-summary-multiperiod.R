source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("summary reports time_restriction and t multipliers", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 50)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = attr(sim, "true_B"),
             time_restriction = "rm", M_bar = 1)

  # Summary output should mention time_restriction
  out <- capture.output(summary(fit))
  expect_true(any(grepl("rm", out, ignore.case = TRUE)),
              label = "summary mentions 'rm'")
  expect_true(any(grepl("t\\s*=\\s*[12]", out)),
              label = "summary mentions t = 1 or t = 2")
})

test_that("summary mentions binding restriction when applicable", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 51)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = 4L, pre_periods = 1:3,
             B = attr(sim, "true_B"),
             time_restriction = "rm", M_bar = 1)

  out <- capture.output(summary(fit))
  # Should mention binding_lower or binding_upper somewhere
  expect_true(any(grepl("binding|bind", out, ignore.case = TRUE)),
              label = "summary mentions binding restriction")
})
