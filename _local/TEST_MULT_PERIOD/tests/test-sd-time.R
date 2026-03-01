source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("SD-Time with B_t=0 gives point identification (linear_selection DGP)", {
  # For linear_selection: delta_tilde_s(d) = rho * d for all s (constant in time)
  # delta_tilde_0(d) = rho * d, so center = Lambda(d,t) - t * rho * d
  # With true ATT(d) = tau(d), Lambda(d,t) = tau(d) + delta_t(d) = tau(d) + t*rho*d
  # Center = tau(d) + t*rho*d - t*rho*d = tau(d) <- point identification
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3,
                      beta = c(1, 0), seed = 30)
  true_att <- attr(sim, "true_att")

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = 0, time_restriction = "sd", B_t = 0)

  for (pp in c(4L, 5L)) {
    att_pp <- fit$att[fit$att$period == pp & fit$att$B == 0, ]
    # Half-width should be zero
    hw <- (att_pp$att_upper - att_pp$att_lower) / 2
    expect_equal(hw, rep(0, length(hw)), tolerance = 1e-10)
    # Center should equal tau(d) up to estimation error
    expect_equal(att_pp$att_lower, true_att(att_pp$d), tolerance = 0.3)
  }
})

test_that("SD-Time half-width = B_t * t*(t+1)/2", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 31)
  B_t_val <- 0.1

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = 0, time_restriction = "sd", B_t = B_t_val)

  # t=1: hw = 0.1 * 1 * 2 / 2 = 0.1 (constant across all d)
  att_t1 <- fit$att[fit$att$period == 4L & fit$att$B == 0, ]
  hw_t1  <- (att_t1$att_upper - att_t1$att_lower) / 2
  expect_equal(hw_t1, rep(B_t_val * 1 * 2 / 2, length(hw_t1)), tolerance = 1e-10)

  # t=2: hw = 0.1 * 2 * 3 / 2 = 0.3 (constant across all d)
  att_t2 <- fit$att[fit$att$period == 5L & fit$att$B == 0, ]
  hw_t2  <- (att_t2$att_upper - att_t2$att_lower) / 2
  expect_equal(hw_t2, rep(B_t_val * 2 * 3 / 2, length(hw_t2)), tolerance = 1e-10)
})

test_that("SD-Time: true ATT inside identified set", {
  sim <- simulate_lpt(n = 1500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 32)
  true_att <- attr(sim, "true_att")
  B_t_val  <- 0.5

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = 0, time_restriction = "sd", B_t = B_t_val)

  for (pp in c(4L, 5L)) {
    att_pp <- fit$att[fit$att$period == pp & fit$att$B == 0, ]
    true_vals <- true_att(att_pp$d)
    expect_true(all(att_pp$att_lower - 0.5 <= true_vals &
                      true_vals <= att_pp$att_upper + 0.5),
                label = paste("True ATT inside SD-Time IS at period", pp))
  }
})
