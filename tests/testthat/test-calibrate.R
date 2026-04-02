test_that("calibrate_B returns correct structure on sru data", {
  data(sru, package = "lpt")
  result <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                        pre_periods = 1993:1999)

  expect_true(is.list(result))
  expect_true(all(c("B_hat", "pre_slopes", "sup_by_period") %in% names(result)))
  expect_true(is.numeric(result$B_hat))
  expect_gt(result$B_hat, 0)
  expect_true(is.data.frame(result$pre_slopes))
  expect_equal(length(result$sup_by_period), 6)  # 7 years -> 6 pairs
  expect_true(all(c("period_pair", "d", "mu_prime_d") %in%
                    names(result$pre_slopes)))
})

test_that("calibrate_B errors with fewer than 2 pre-periods", {
  data(sru, package = "lpt")
  expect_error(
    calibrate_B(sru, "commune", "year", "outcome", "dose", pre_periods = 1999),
    "at least 2"
  )
})
