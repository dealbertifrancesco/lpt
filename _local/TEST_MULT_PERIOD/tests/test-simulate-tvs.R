source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("simulate_lpt supports time_varying_selection DGP", {
  expect_no_error(
    simulate_lpt(n = 200, n_pre = 3, n_post = 2,
                 dgp = "time_varying_selection", rho = 0.2, rho2 = 0.05, seed = 40)
  )
})

test_that("time_varying_selection has correct true_B and true_delta_tilde_0", {
  sim <- simulate_lpt(n = 200, n_pre = 3, n_post = 2,
                      dgp = "time_varying_selection", rho = 0.2, rho2 = 0.05, seed = 41)

  # true_B should be abs(rho) (Lipschitz bound on cross-sectional slope)
  expect_equal(attr(sim, "true_B"), abs(0.2))

  # true_delta_tilde_0 should be stored
  expect_true(!is.null(attr(sim, "true_delta_tilde_0")))

  # For n_pre=3: true_delta_tilde_0(d) = rho*d + rho2*d*(2*3-1) = 0.2*d + 0.05*d*5
  true_dt0 <- attr(sim, "true_delta_tilde_0")
  test_d <- c(0, 1, 2, 3)
  expected <- 0.2 * test_d + 0.05 * test_d * (2 * 3 - 1)
  expect_equal(true_dt0(test_d), expected, tolerance = 1e-10)
})

test_that("linear_selection also carries true_delta_tilde_0 attribute", {
  sim <- simulate_lpt(n = 200, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 42)

  # true_delta_tilde_0(d) = rho * d = 0.3 * d for linear_selection
  true_dt0 <- attr(sim, "true_delta_tilde_0")
  expect_true(!is.null(true_dt0))
  test_d <- c(0, 1, 2)
  expect_equal(true_dt0(test_d), 0.3 * test_d, tolerance = 1e-10)
})
