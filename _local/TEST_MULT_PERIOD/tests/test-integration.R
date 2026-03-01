source("_local/TEST_MULT_PERIOD/tests/helpers.R")

# --- Test 1: LPT t multiplier ---
test_that("Test 1 - LPT t multiplier: IS contains true ATT for all (d,t)", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 5,
                      dgp = "linear_selection", seed = 42)
  B_true   <- attr(sim, "true_B")
  true_att <- attr(sim, "true_att")

  fit <- lpt(sim, "id", "time", "outcome", "dose",
             post_period = 4:8, pre_periods = 1:3, B = B_true)

  for (pp in 4:8) {
    att_pp    <- fit$att[fit$att$period == pp & fit$att$B == B_true, ]
    t_val     <- fit$t_values[[as.character(pp)]]
    true_vals <- true_att(att_pp$d)

    # True ATT inside IS (with small tolerance for finite sample)
    expect_true(all(att_pp$att_lower <= true_vals + 0.5),
                label = paste("lower bound, period", pp))
    expect_true(all(true_vals <= att_pp$att_upper + 0.5),
                label = paste("upper bound, period", pp))

    # Half-width = t * B * d
    hw_expected <- t_val * B_true * att_pp$d
    hw_actual   <- (att_pp$att_upper - att_pp$att_lower) / 2
    expect_equal(hw_actual, hw_expected, tolerance = 1e-10,
                 label = paste("half-width formula at period", pp))
  }

  # Verify t values are 1..5 (names differ, so use unname)
  t_vals_all <- unlist(fit$t_values[as.character(4:8)])
  expect_equal(unname(t_vals_all), 1:5, label = "t values 1..5 for post-periods 4..8")
})

# --- Test 2: RM-Time ---
test_that("Test 2 - RM-Time: delta_star_pre ~= |rho| * d, IS contains true ATT", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3, seed = 43)
  B_true   <- attr(sim, "true_B")  # 0.3
  true_att <- attr(sim, "true_att")

  fit_rm <- lpt(sim, "id", "time", "outcome", "dose",
                post_period = c(4L, 5L), pre_periods = 1:3,
                B = B_true, time_restriction = "rm", M_bar = 1)

  # delta_star_pre(d) ~= |rho| * d = 0.3 * d
  dsp <- fit_rm$calibration$delta_star_pre
  expect_equal(dsp$delta_star_pre, abs(0.3 * dsp$d), tolerance = 0.20,
               label = "delta_star_pre ~= 0.3 * d")

  # True ATT inside IS
  for (pp in c(4L, 5L)) {
    att_pp    <- fit_rm$att[fit_rm$att$period == pp & fit_rm$att$B == B_true, ]
    true_vals <- true_att(att_pp$d)
    expect_true(all(att_pp$att_lower <= true_vals + 0.5),
                label = paste("lower, period", pp))
    expect_true(all(true_vals <= att_pp$att_upper + 0.5),
                label = paste("upper, period", pp))
  }

  # Joint set is tighter than LPT alone (or equal)
  fit_lpt <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = c(4L, 5L), pre_periods = 1:3, B = B_true)
  hw_lpt  <- (fit_lpt$att$att_upper - fit_lpt$att$att_lower) / 2
  hw_rm   <- (fit_rm$att$att_upper - fit_rm$att$att_lower) / 2
  expect_true(all(hw_rm <= hw_lpt + 1e-10),
              label = "Joint RM+LPT no wider than LPT alone")
})

# --- Test 3: SD-Time ---
test_that("Test 3 - SD-Time B_t=0: point identification in linear_selection DGP", {
  sim <- simulate_lpt(n = 2000, n_pre = 3, n_post = 2,
                      dgp = "linear_selection", rho = 0.3,
                      beta = c(1, 0), seed = 44)
  true_att <- attr(sim, "true_att")

  fit_sd <- lpt(sim, "id", "time", "outcome", "dose",
                post_period = c(4L, 5L), pre_periods = 1:3,
                B = 0, time_restriction = "sd", B_t = 0)

  for (pp in c(4L, 5L)) {
    att_pp <- fit_sd$att[fit_sd$att$period == pp & fit_sd$att$B == 0, ]
    hw <- (att_pp$att_upper - att_pp$att_lower) / 2
    expect_equal(hw, rep(0, length(hw)), tolerance = 1e-10,
                 label = paste("zero half-width at period", pp))

    # Center ~= true ATT (up to finite-sample error)
    expect_equal(att_pp$att_lower, true_att(att_pp$d), tolerance = 0.4,
                 label = paste("center ~= true ATT at period", pp))
  }

  # With B_t > 0, IS widens at t=2
  fit_sd2 <- lpt(sim, "id", "time", "outcome", "dose",
                 post_period = c(4L, 5L), pre_periods = 1:3,
                 B = 0, time_restriction = "sd", B_t = 0.2)
  att_t2_sd2 <- fit_sd2$att[fit_sd2$att$period == 5L & fit_sd2$att$B == 0, ]
  hw_t2 <- (att_t2_sd2$att_upper - att_t2_sd2$att_lower) / 2
  expect_equal(hw_t2, rep(0.2 * 2 * 3 / 2, length(hw_t2)), tolerance = 1e-10,
               label = "SD-Time hw at t=2 = B_t * 2 * 3/2 = 0.6")
})

# --- Test 4: Joint restriction tightening ---
test_that("Test 4 - Joint LPT + RM is tighter than either alone", {
  sim <- simulate_lpt(n = 1000, n_pre = 3, n_post = 1,
                      dgp = "linear_selection", rho = 0.3, seed = 99)
  B_true <- attr(sim, "true_B")

  fit_lpt   <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3, B = B_true)
  fit_rm    <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3,
                   B = 0, time_restriction = "rm", M_bar = 1)
  fit_joint <- lpt(sim, "id", "time", "outcome", "dose",
                   post_period = 4L, pre_periods = 1:3,
                   B = B_true, time_restriction = "rm", M_bar = 1)

  hw_lpt   <- (fit_lpt$att$att_upper - fit_lpt$att$att_lower) / 2
  hw_rm    <- (fit_rm$att$att_upper - fit_rm$att$att_lower) / 2
  hw_joint <- (fit_joint$att$att_upper - fit_joint$att$att_lower) / 2

  expect_true(all(hw_joint <= hw_lpt + 1e-10),
              label = "Joint no wider than LPT alone")
  expect_true(all(hw_joint <= hw_rm + 1e-10),
              label = "Joint no wider than RM alone")
})

# --- Test 5: Real SRU data ---
test_that("Test 5 - SRU data: later periods have wider IS and t values correct", {
  skip_if(!file.exists("data/sru.rda"), "sru data not found")
  data(sru)

  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2005, 2010, 2015, 2019),
             pre_periods = 1993:1999,
             B = "calibrate",
             time_restriction = "rm", M_bar = 1)

  # t values should be 6, 11, 16, 20 (years after 1999 ref period)
  expect_equal(fit$t_values[["2005"]], 6L,  label = "t for 2005")
  expect_equal(fit$t_values[["2010"]], 11L, label = "t for 2010")
  expect_equal(fit$t_values[["2015"]], 16L, label = "t for 2015")
  expect_equal(fit$t_values[["2019"]], 20L, label = "t for 2019")

  # Later post-periods should have wider ATT^o IS
  get_hw_atto <- function(pp) {
    row <- fit$att_o[fit$att_o$period == pp & fit$att_o$B == fit$B_hat, ]
    if (nrow(row) == 0) return(NA_real_)
    (row$att_o_upper - row$att_o_lower) / 2
  }

  hw_2005 <- get_hw_atto(2005)
  hw_2019 <- get_hw_atto(2019)
  expect_gt(hw_2019, hw_2005, label = "2019 IS wider than 2005 IS")
})

# --- Test 6: Backward compatibility ---
test_that("Test 6 - Backward compat: single period, time_restriction='none'", {
  skip_if(!file.exists("data/sru.rda"), "sru data not found")
  data(sru)

  fit <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0.1, time_restriction = "none")

  expect_s3_class(fit, "lpt")
  expect_equal(fit$t_values[["2019"]], 20L, label = "t for 2019 = 20")

  # dATT half-width should be t * B = 20 * 0.1 = 2.0 (constant across d)
  datt_row <- fit$datt[fit$datt$B == 0.1, ]
  hw_datt  <- (datt_row$datt_upper - datt_row$datt_lower) / 2
  expect_equal(hw_datt, rep(20 * 0.1, length(hw_datt)), tolerance = 1e-10,
               label = "dATT hw = t * B = 2.0")

  # ATT^o half-width should be t * B * D_bar = 20 * 0.1 * D_bar
  atto_row <- fit$att_o[fit$att_o$B == 0.1, ]
  hw_atto  <- (atto_row$att_o_upper - atto_row$att_o_lower) / 2
  expect_equal(hw_atto, 20 * 0.1 * atto_row$D_bar, tolerance = 1e-10,
               label = "ATT^o hw = t * B * D_bar")
})
