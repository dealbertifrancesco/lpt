# test-lpt-b.R
cat("Testing lpt() with lpt_type = 'b'...\n")

# --- Test 1: LPT-b with B=0 (same as LPT-a with B=0) ---
fit_b0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = 0, lpt_type = "b")

stopifnot(inherits(fit_b0, "lpt"))
stopifnot(fit_b0$lpt_type == "b")
datt0 <- fit_b0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
cat("  LPT-b, B=0: point identification. OK\n")

# --- Test 2: LPT-b with B>0 — bounds GROW with t ---
fit_b_b1 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = 0.5, lpt_type = "b")

datt_2018 <- fit_b_b1$datt[fit_b_b1$datt$period == 2018, ]
datt_2019 <- fit_b_b1$datt[fit_b_b1$datt$period == 2019, ]
t_2018 <- unique(datt_2018$t_index)
t_2019 <- unique(datt_2019$t_index)

# Under LPT-b, datt bounds width = 2 * t * B_d
width_2018 <- datt_2018$datt_upper[1] - datt_2018$datt_lower[1]
width_2019 <- datt_2019$datt_upper[1] - datt_2019$datt_lower[1]
expected_2018 <- 2 * t_2018 * 0.5
expected_2019 <- 2 * t_2019 * 0.5
stopifnot(abs(width_2018 - expected_2018) < 1e-10)
stopifnot(abs(width_2019 - expected_2019) < 1e-10)
stopifnot(width_2019 > width_2018)  # later period is wider!
cat(sprintf("  LPT-b, B=0.5: dATT width grows with t. t(2018)=%d->w=%.1f, t(2019)=%d->w=%.1f. OK\n",
            t_2018, width_2018, t_2019, width_2019))

# --- Test 3: LPT-b calibration uses first_diff ---
fit_b_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "b")

stopifnot(fit_b_cal$calibration$type == "first_diff")
cat(sprintf("  LPT-b, calibrated: B_hat=%.4f (first_diff). OK\n", fit_b_cal$B_hat))

# --- Test 4: ATT width grows with t under LPT-b ---
att_2018 <- fit_b_b1$att[fit_b_b1$att$period == 2018, ]
att_2019 <- fit_b_b1$att[fit_b_b1$att$period == 2019, ]
# Width at same dose d: 2*t*B_d*d
# Pick a dose point present in both
d_val <- att_2018$d[25]
w_2018 <- att_2018$att_upper[25] - att_2018$att_lower[25]
w_2019 <- att_2019$att_upper[25] - att_2019$att_lower[25]
expected_ratio <- t_2019 / t_2018
actual_ratio <- w_2019 / w_2018
stopifnot(abs(actual_ratio - expected_ratio) < 1e-6)
cat(sprintf("  LPT-b: ATT width ratio = %.2f (expected t2/t1 = %.2f). OK\n",
            actual_ratio, expected_ratio))

# --- Test 5: ATT^o width grows with t under LPT-b ---
atto_2018 <- fit_b_b1$att_o[fit_b_b1$att_o$period == 2018, ]
atto_2019 <- fit_b_b1$att_o[fit_b_b1$att_o$period == 2019, ]
w_o_2018 <- atto_2018$att_o_upper - atto_2018$att_o_lower
w_o_2019 <- atto_2019$att_o_upper - atto_2019$att_o_lower
stopifnot(w_o_2019 > w_o_2018)
cat("  LPT-b: ATT^o width grows with t. OK\n")

cat("lpt_type='b' tests PASSED.\n")
