source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("calibrate_B returns delta_tilde_s levels", {
  # Use n=2000 so GAM has enough data for accurate level estimation
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 10)

  cal <- calibrate_B(sim, "id", "time", "outcome", "dose",
                     pre_periods = 1:3)

  # Should have pre_slopes with delta_tilde_d column
  expect_true("delta_tilde_d" %in% names(cal$pre_slopes))

  # For linear_selection DGP, delta_tilde_s(d) = rho * d at any d
  # So delta_tilde_d should be approximately 0.3 * d at each eval point
  pair1_dt <- cal$pre_slopes[cal$pre_slopes$period_pair == "1-2", ]
  expected <- 0.3 * pair1_dt$d
  # GAM level estimation has more noise than derivative; use relative tol = 0.30
  expect_equal(pair1_dt$delta_tilde_d, expected, tolerance = 0.30)
})

test_that("calibrate_B returns delta_star_pre (max |delta_tilde_s| per d)", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 11)

  cal <- calibrate_B(sim, "id", "time", "outcome", "dose",
                     pre_periods = 1:3)

  expect_true(!is.null(cal$delta_star_pre))
  expect_true(all(cal$delta_star_pre$delta_star_pre >= 0))
  # delta_star_pre should equal max |delta_tilde_s(d)| over pairs
  # In this DGP all pairs have same delta_tilde ~= 0.3*d
  expect_equal(cal$delta_star_pre$delta_star_pre,
               abs(cal$delta_star_pre$d * 0.3),
               tolerance = 0.30)
})

test_that("calibrate_B returns delta_tilde_0 (last pre-period increment)", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 12)

  cal <- calibrate_B(sim, "id", "time", "outcome", "dose",
                     pre_periods = 1:3)

  # delta_tilde_0 is the delta_tilde for the last pre-period pair (2-3)
  expect_true(!is.null(cal$delta_tilde_0))
  expect_equal(cal$delta_tilde_0$delta_tilde_0, 0.3 * cal$delta_tilde_0$d,
               tolerance = 0.30)
})
