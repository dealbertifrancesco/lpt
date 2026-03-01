# test-simulated-data.R
cat("Testing with simulated data...\n")

# --- no_selection DGP: B=0 should recover truth ---
sim_ns <- simulate_lpt(n = 2000, n_pre = 4, n_post = 3,
                         dgp = "no_selection", seed = 42)

fit_ns <- lpt(sim_ns, "id", "time", "outcome", "dose",
              post_period = 5:7, pre_periods = 1:4,
              B = 0, lpt_type = "C")

true_datt <- attr(sim_ns, "true_datt")
datt_pp5 <- fit_ns$datt[fit_ns$datt$period == 5, ]
true_vals <- true_datt(datt_pp5$d)
cor_val <- cor(datt_pp5$lambda_d, true_vals)
stopifnot(cor_val > 0.9)
cat(sprintf("  no_selection: cor(lambda, tau') = %.3f. OK\n", cor_val))

# --- linear_selection DGP: LPT-C and LPT-P should both work ---
sim_ls <- simulate_lpt(n = 2000, n_pre = 4, n_post = 3,
                         dgp = "linear_selection", seed = 42)

fit_ls_c <- lpt(sim_ls, "id", "time", "outcome", "dose",
                post_period = 5:7, pre_periods = 1:4,
                B = "calibrate", lpt_type = "C")
fit_ls_p <- lpt(sim_ls, "id", "time", "outcome", "dose",
                post_period = 5:7, pre_periods = 1:4,
                B = "calibrate", lpt_type = "P")

stopifnot(fit_ls_c$B_hat > 0)
stopifnot(fit_ls_p$B_hat > 0)
stopifnot(fit_ls_c$calibration$type == "cumulative")
stopifnot(fit_ls_p$calibration$type == "first_diff")
cat(sprintf("  linear_selection: B_C=%.4f, B_P=%.4f. OK\n",
            fit_ls_c$B_hat, fit_ls_p$B_hat))

# --- LPT-P bounds should grow with t ---
datt_p_5 <- fit_ls_p$datt[fit_ls_p$datt$period == 5, ]
datt_p_7 <- fit_ls_p$datt[fit_ls_p$datt$period == 7, ]
w5 <- datt_p_5$datt_upper[1] - datt_p_5$datt_lower[1]
w7 <- datt_p_7$datt_upper[1] - datt_p_7$datt_lower[1]
stopifnot(w7 > w5)
cat(sprintf("  LPT-P: width at t=1 is %.3f, at t=3 is %.3f (grows). OK\n", w5, w7))

# --- LPT-C bounds should be constant ---
datt_c_5 <- fit_ls_c$datt[fit_ls_c$datt$period == 5, ]
datt_c_7 <- fit_ls_c$datt[fit_ls_c$datt$period == 7, ]
w5c <- datt_c_5$datt_upper[1] - datt_c_5$datt_lower[1]
w7c <- datt_c_7$datt_upper[1] - datt_c_7$datt_lower[1]
stopifnot(abs(w5c - w7c) < 1e-10)
cat(sprintf("  LPT-C: width constant at %.3f across periods. OK\n", w5c))

cat("Simulated data tests PASSED.\n")
