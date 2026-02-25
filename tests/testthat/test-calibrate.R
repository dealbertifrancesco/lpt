test_that("calibrate_B recovers true B under linear_selection", {
  dat <- simulate_lpt(n = 5000, n_pre = 4, dgp = "linear_selection",
                       rho = 0.3, seed = 42)
  result <- calibrate_B(dat, "id", "time", "outcome", "dose",
                         pre_periods = 1:4)
  expect_gt(result$B_hat, 0.15)
  expect_lt(result$B_hat, 0.55)
})

test_that("calibrate_B returns near-zero under no_selection", {
  dat <- simulate_lpt(n = 3000, n_pre = 4, dgp = "no_selection", seed = 42)
  result <- calibrate_B(dat, "id", "time", "outcome", "dose",
                         pre_periods = 1:4)
  expect_lt(result$B_hat, 0.2)
})

test_that("calibrate_B errors with fewer than 2 pre-periods", {
  dat <- simulate_lpt(n = 100, n_pre = 1, dgp = "no_selection", seed = 1)
  expect_error(
    calibrate_B(dat, "id", "time", "outcome", "dose", pre_periods = 1),
    "at least 2"
  )
})

test_that("calibrate_B returns correct structure", {
  dat <- simulate_lpt(n = 500, n_pre = 4, dgp = "linear_selection", seed = 1)
  result <- calibrate_B(dat, "id", "time", "outcome", "dose",
                         pre_periods = 1:4)

  expect_true(is.list(result))
  expect_true(all(c("B_hat", "pre_slopes", "sup_by_period") %in% names(result)))
  expect_true(is.numeric(result$B_hat))
  expect_true(is.data.frame(result$pre_slopes))
  expect_equal(length(result$sup_by_period), 3)
  expect_true(all(c("period_pair", "d", "mu_prime_d", "se") %in%
                    names(result$pre_slopes)))
})
