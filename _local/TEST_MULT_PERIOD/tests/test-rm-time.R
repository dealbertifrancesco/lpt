source("_local/TEST_MULT_PERIOD/tests/helpers.R")

test_that("lpt accepts time_restriction = 'rm' argument without error", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 20)
  expect_no_error(
    lpt(sim, "id", "time", "outcome", "dose",
        post_period = c(4L, 5L), pre_periods = 1:3,
        B = attr(sim, "true_B"),
        time_restriction = "rm", M_bar = 1)
  )
})

test_that("RM-Time ATT half-width = t * M_bar * delta_star_pre(d)", {
  sim <- simulate_lpt(n = 1000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 21)
  B_true <- attr(sim, "true_B")

  fit_rm <- lpt(sim, "id", "time", "outcome", "dose",
                post_period = c(4L, 5L), pre_periods = 1:3,
                B = 0,  # no LPT, pure RM-Time
                time_restriction = "rm", M_bar = 1)

  # For linear_selection DGP: delta_star_pre(d) ~= rho * d = 0.3 * d
  # RM-Time half-width at t=1 (post=4): ~= 1 * 1 * 0.3 * d = 0.3 * d
  att_t1 <- fit_rm$att[fit_rm$att$period == 4L & fit_rm$att$B == 0, ]
  hw_t1  <- (att_t1$att_upper - att_t1$att_lower) / 2
  expect_equal(hw_t1, fit_rm$calibration$delta_star_pre$delta_star_pre,
               tolerance = 0.2, label = "RM half-width at t=1 ~= delta_star_pre")

  # At t=2, half-width should be double t=1 (exclude d=0 where hw=0 → NaN ratio)
  att_t2    <- fit_rm$att[fit_rm$att$period == 5L & fit_rm$att$B == 0, ]
  hw_t2     <- (att_t2$att_upper - att_t2$att_lower) / 2
  pos_idx   <- hw_t1 > 0
  expect_equal(hw_t2[pos_idx] / hw_t1[pos_idx], rep(2, sum(pos_idx)),
               tolerance = 0.01, label = "RM half-width doubles from t=1 to t=2")
})

test_that("RM-Time joint with LPT produces intersection (tighter bounds)", {
  sim <- simulate_lpt(n = 500, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 22)
  B_true <- attr(sim, "true_B")

  fit_lpt <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = 4L, pre_periods = 1:3, B = B_true)
  fit_rm  <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = 4L, pre_periods = 1:3, B = 0,
                 time_restriction = "rm", M_bar = 1)
  fit_joint <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3, B = B_true,
                   time_restriction = "rm", M_bar = 1)

  # Joint set should be at least as tight as LPT or RM alone at each point
  hw_lpt   <- (fit_lpt$att$att_upper - fit_lpt$att$att_lower) / 2
  hw_rm    <- (fit_rm$att$att_upper - fit_rm$att$att_lower) / 2
  hw_joint <- (fit_joint$att$att_upper - fit_joint$att$att_lower) / 2

  expect_true(all(hw_joint <= hw_lpt + 1e-10),
              label = "Joint no wider than LPT alone")
  expect_true(all(hw_joint <= hw_rm + 1e-10),
              label = "Joint no wider than RM alone")
})

test_that("RM-Time: true ATT is inside the identified set", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 23)
  B_true   <- attr(sim, "true_B")
  true_att <- attr(sim, "true_att")

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = c(4L, 5L), pre_periods = 1:3,
             B = B_true, time_restriction = "rm", M_bar = 1)

  for (pp in c(4L, 5L)) {
    att_pp <- fit$att[fit$att$period == pp & fit$att$B == B_true, ]
    true_vals <- true_att(att_pp$d)
    expect_true(all(att_pp$att_lower <= true_vals + 0.5 &
                      true_vals - 0.5 <= att_pp$att_upper),
                label = paste("True ATT inside IS at period", pp))
  }
})
