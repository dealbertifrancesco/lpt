source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("dATT half-width scales with t", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 3, dgp = "linear_selection",
                      beta = c(1, 0), rho = 0.3, seed = 3)
  B_true <- attr(sim, "true_B")  # 0.3

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L), pre_periods = 1:2, B = B_true)

  # dATT half-width at each t should be t * B
  for (t_idx in 1:3) {
    pp_t  <- 2L + t_idx
    datt_pp <- fit$datt[fit$datt$period == pp_t & fit$datt$B == B_true, ]
    hw <- (datt_pp$datt_upper - datt_pp$datt_lower) / 2
    # hw is constant across eval points (lambda_d cancels); compare all values
    expect_equal(hw, rep(t_idx * B_true, length(hw)), tolerance = 1e-10,
                 label = paste("dATT half-width at t =", t_idx))
  }
})

test_that("ATT^o half-width scales with t", {
  sim <- simulate_lpt(n = 500, n_pre = 2, n_post = 3, dgp = "linear_selection",
                      beta = c(1, 0), rho = 0.3, seed = 4)
  B_true <- attr(sim, "true_B")  # 0.3

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(3L, 4L, 5L), pre_periods = 1:2, B = B_true)

  D_bar <- fit$att_o$D_bar[fit$att_o$period == 3L & fit$att_o$B == B_true]

  for (t_idx in 1:3) {
    pp_t  <- 2L + t_idx
    atto_row <- fit$att_o[fit$att_o$period == pp_t & fit$att_o$B == B_true, ]
    hw <- (atto_row$att_o_upper - atto_row$att_o_lower) / 2
    expect_equal(hw, t_idx * B_true * atto_row$D_bar, tolerance = 1e-10,
                 label = paste("ATT^o half-width at t =", t_idx))
  }
})
