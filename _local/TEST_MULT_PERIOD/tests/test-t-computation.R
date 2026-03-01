source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("t computation is correct for annual calendar-year data", {
  # Use simulate_lpt as proxy: n_pre=3, n_post=5 -> time indices 1..8
  # ref_period = 3 (last pre-period). Post-periods: 4,5,6,7,8.
  sim <- simulate_lpt(n = 200, n_pre = 3, n_post = 5, dgp = "linear_selection", seed = 1)
  all_times <- sort(unique(sim$time))
  ref_period <- 3L

  # t for post_period=4 should be 1; for post_period=6 should be 3
  t4 <- sum(all_times > ref_period & all_times <= 4L)
  t6 <- sum(all_times > ref_period & all_times <= 6L)

  expect_equal(t4, 1L)
  expect_equal(t6, 3L)
})

test_that("t computation is stored in lpt output", {
  sim <- simulate_lpt(n = 300, n_pre = 2, n_post = 3, dgp = "linear_selection", seed = 42)
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L),
             pre_periods = 1:2,
             B = attr(sim, "true_B"))

  # t values should be stored in fit$t_values
  expect_true(!is.null(fit$t_values))
  # ref_period = 2, so t for post 3 = 1, post 4 = 2, post 5 = 3
  expect_equal(fit$t_values[["3"]], 1L)
  expect_equal(fit$t_values[["4"]], 2L)
  expect_equal(fit$t_values[["5"]], 3L)
})
