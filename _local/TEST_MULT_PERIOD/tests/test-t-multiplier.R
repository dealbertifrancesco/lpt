source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("ATT bounds half-width scales with t (LPT alone)", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 4, dgp = "linear_selection",
                      beta = c(1, 0), rho = 0.3, seed = 1)
  B_true <- attr(sim, "true_B")  # = 0.3

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L, 6L),
             pre_periods = 1:2,
             B = B_true)

  # For a dose d, ATT half-width at post-period t should be t * B * d
  # Check at median dose
  ep <- fit$att$d[fit$att$period == 3L]
  d_med <- ep[which.min(abs(ep - stats::median(ep)))]

  for (t_idx in 1:4) {
    pp_t <- 2L + t_idx  # post-periods are 3,4,5,6
    att_pp <- fit$att[fit$att$period == pp_t & fit$att$B == B_true, ]
    row <- att_pp[which.min(abs(att_pp$d - d_med)), ]
    expected_hw <- t_idx * B_true * row$d
    actual_hw   <- (row$att_upper - row$att_lower) / 2
    expect_equal(actual_hw, expected_hw, tolerance = 1e-10,
                 label = paste("half-width at t =", t_idx))
  }
})

test_that("ATT bounds are wider for later periods", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 3, dgp = "linear_selection", seed = 2)
  B_true <- attr(sim, "true_B")
  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L), pre_periods = 1:2, B = B_true)

  # Use a non-zero dose: first positive eval point to avoid d=0 giving hw=0
  ep_vals <- unique(fit$att$d[fit$att$period == 3L])
  ep <- ep_vals[ep_vals > 0][1]

  hw_t1 <- with(fit$att[fit$att$period == 3L & fit$att$B == B_true, ],
                (att_upper - att_lower)[which.min(abs(d - ep))])
  hw_t2 <- with(fit$att[fit$att$period == 4L & fit$att$B == B_true, ],
                (att_upper - att_lower)[which.min(abs(d - ep))])
  hw_t3 <- with(fit$att[fit$att$period == 5L & fit$att$B == B_true, ],
                (att_upper - att_lower)[which.min(abs(d - ep))])

  expect_lt(hw_t1, hw_t2)
  expect_lt(hw_t2, hw_t3)
})
