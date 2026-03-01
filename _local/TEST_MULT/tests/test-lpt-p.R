# test-lpt-p.R
cat("Testing lpt() with lpt_type = 'P'...\n")

fit_p0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = 0, lpt_type = "P")
stopifnot(inherits(fit_p0, "lpt"))
stopifnot(fit_p0$lpt_type == "P")
datt0 <- fit_p0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
cat("  LPT-P, B=0: point identification. OK\n")

# LPT-P: dATT width = 2*t*B, grows linearly with t
fit_p_b1 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = 0.5, lpt_type = "P")
datt_2018 <- fit_p_b1$datt[fit_p_b1$datt$period == 2018, ]
datt_2019 <- fit_p_b1$datt[fit_p_b1$datt$period == 2019, ]
t_2018    <- unique(datt_2018$t_index)
t_2019    <- unique(datt_2019$t_index)
width_2018 <- datt_2018$datt_upper[1] - datt_2018$datt_lower[1]
width_2019 <- datt_2019$datt_upper[1] - datt_2019$datt_lower[1]
stopifnot(abs(width_2018 - 2 * t_2018 * 0.5) < 1e-10)
stopifnot(abs(width_2019 - 2 * t_2019 * 0.5) < 1e-10)
stopifnot(width_2019 > width_2018)
cat(sprintf("  LPT-P, B=0.5: width grows. t(2018)=%d->w=%.1f, t(2019)=%d->w=%.1f. OK\n",
            t_2018, width_2018, t_2019, width_2019))

fit_p_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "P")
stopifnot(fit_p_cal$calibration$type == "first_diff")
cat(sprintf("  LPT-P, calibrated: B_hat=%.4f (first_diff). OK\n", fit_p_cal$B_hat))

att_2018 <- fit_p_b1$att[fit_p_b1$att$period == 2018, ]
att_2019 <- fit_p_b1$att[fit_p_b1$att$period == 2019, ]
expected_ratio <- t_2019 / t_2018
actual_ratio   <- (att_2019$att_upper[25] - att_2019$att_lower[25]) /
                  (att_2018$att_upper[25] - att_2018$att_lower[25])
stopifnot(abs(actual_ratio - expected_ratio) < 1e-6)
cat(sprintf("  LPT-P: ATT width ratio = %.2f (expected t2/t1 = %.2f). OK\n",
            actual_ratio, expected_ratio))

atto_2018 <- fit_p_b1$att_o[fit_p_b1$att_o$period == 2018, ]
atto_2019 <- fit_p_b1$att_o[fit_p_b1$att_o$period == 2019, ]
stopifnot((atto_2019$att_o_upper - atto_2019$att_o_lower) >
          (atto_2018$att_o_upper - atto_2018$att_o_lower))
cat("  LPT-P: ATT^o width grows with t. OK\n")

cat("lpt_type='P' tests PASSED.\n")
