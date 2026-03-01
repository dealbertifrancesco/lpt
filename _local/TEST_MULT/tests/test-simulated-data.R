# test-simulated-data.R
cat("Testing with simulated data...\n")

# --- no_selection DGP: B=0 should recover truth ---
sim_ns <- simulate_lpt(n = 2000, n_pre = 4, n_post = 3,
                         dgp = "no_selection", seed = 42)

fit_ns <- lpt(sim_ns, "id", "time", "outcome", "dose",
              post_period = 5:7, pre_periods = 1:4,
              B = 0, lpt_type = "a")

# With no selection, lambda(d) should approximate tau'(d)
true_datt <- attr(sim_ns, "true_datt")
datt_pp5 <- fit_ns$datt[fit_ns$datt$period == 5, ]
true_vals <- true_datt(datt_pp5$d)
# Check correlation is high (not exact due to noise)
cor_val <- cor(datt_pp5$lambda_d, true_vals)
stopifnot(cor_val > 0.9)
cat(sprintf("  no_selection: cor(lambda, tau') = %.3f. OK\n", cor_val))

# --- linear_selection DGP: LPT-a and LPT-b should both work ---
sim_ls <- simulate_lpt(n = 2000, n_pre = 4, n_post = 3,
                         dgp = "linear_selection", seed = 42)

fit_ls_a <- lpt(sim_ls, "id", "time", "outcome", "dose",
                post_period = 5:7, pre_periods = 1:4,
                B = "calibrate", lpt_type = "a")
fit_ls_b <- lpt(sim_ls, "id", "time", "outcome", "dose",
                post_period = 5:7, pre_periods = 1:4,
                B = "calibrate", lpt_type = "b")

# Both should have positive B
stopifnot(fit_ls_a$B_hat > 0)
stopifnot(fit_ls_b$B_hat > 0)
# LPT-a uses cumulative calibration
stopifnot(fit_ls_a$calibration$type == "cumulative")
stopifnot(fit_ls_b$calibration$type == "first_diff")
cat(sprintf("  linear_selection: B_a=%.4f, B_b=%.4f. OK\n",
            fit_ls_a$B_hat, fit_ls_b$B_hat))

# --- LPT-b bounds should grow with t ---
datt_b_5 <- fit_ls_b$datt[fit_ls_b$datt$period == 5, ]
datt_b_7 <- fit_ls_b$datt[fit_ls_b$datt$period == 7, ]
w5 <- datt_b_5$datt_upper[1] - datt_b_5$datt_lower[1]
w7 <- datt_b_7$datt_upper[1] - datt_b_7$datt_lower[1]
stopifnot(w7 > w5)
cat(sprintf("  LPT-b: width at t=1 is %.3f, at t=3 is %.3f (grows). OK\n", w5, w7))

# --- LPT-a bounds should be constant ---
datt_a_5 <- fit_ls_a$datt[fit_ls_a$datt$period == 5, ]
datt_a_7 <- fit_ls_a$datt[fit_ls_a$datt$period == 7, ]
w5a <- datt_a_5$datt_upper[1] - datt_a_5$datt_lower[1]
w7a <- datt_a_7$datt_upper[1] - datt_a_7$datt_lower[1]
stopifnot(abs(w5a - w7a) < 1e-10)
cat(sprintf("  LPT-a: width constant at %.3f across periods. OK\n", w5a))

cat("Simulated data tests PASSED.\n")
