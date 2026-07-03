test_that("calibrate_bounds returns correct structure on sru data", {
  data(sru, package = "lpt")
  result <- calibrate_bounds(sru, "commune", "year", "outcome", "dose",
                             pre_periods = -7:-1)

  expect_true(is.list(result))
  expect_true(all(c("M_hat", "B_hat", "pre_deviations", "pre_slopes",
                     "sup_dev_by_period", "sup_slope_by_period")
                    %in% names(result)))
  expect_true(is.numeric(result$M_hat))
  expect_true(is.numeric(result$B_hat))
  expect_gt(result$M_hat, 0)
  expect_gt(result$B_hat, 0)

  # 7 pre-periods -> 6 consecutive pairs
  expect_equal(length(result$sup_slope_by_period), 6)
  expect_equal(length(result$sup_dev_by_period), 6)

  expect_true(is.data.frame(result$pre_slopes))
  expect_true(all(c("period_pair", "d", "mu_prime_d") %in%
                    names(result$pre_slopes)))
  expect_true(is.data.frame(result$pre_deviations))
  expect_true(all(c("period_pair", "d", "deviation") %in%
                    names(result$pre_deviations)))

  # Maxima are consistent with the per-pair suprema
  expect_equal(result$M_hat, max(result$sup_dev_by_period))
  expect_equal(result$B_hat, max(result$sup_slope_by_period))
})

test_that("calibrate_bounds errors with fewer than 2 pre-periods", {
  data(sru, package = "lpt")
  expect_error(
    calibrate_bounds(sru, "commune", "year", "outcome", "dose",
                     pre_periods = -1),
    "at least 2"
  )
})

test_that("calibrate_bounds returns NA M_hat without untreated units", {
  data(sru, package = "lpt")
  sru_treated <- sru[sru$dose > 0, ]
  result <- calibrate_bounds(sru_treated, "commune", "year", "outcome",
                             "dose", pre_periods = -7:-1)

  expect_true(is.na(result$M_hat))
  expect_null(result$pre_deviations)
  expect_null(result$sup_dev_by_period)
  expect_gt(result$B_hat, 0)
})
