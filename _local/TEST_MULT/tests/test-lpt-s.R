# test-lpt-s.R
cat("Testing lpt() with lpt_type = 'S'...\n")

# --- Test 1: LPT-S with b0=0, C=0 collapses to point identification ---
fit_s0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = c(b0 = 0, C = 0), lpt_type = "S")
stopifnot(inherits(fit_s0, "lpt"))
stopifnot(fit_s0$lpt_type == "S")
datt0 <- fit_s0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
stopifnot(all(abs(datt0$datt_upper - datt0$lambda_d) < 1e-10))
cat("  LPT-S, b0=0, C=0: point identification. OK\n")

# --- Test 2: LPT-S IS width = 2*Phi(t) where Phi(t) = t*b0 + C*t*(t+1)/2 ---
fit_s_multi <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = c(2018, 2019), pre_periods = 1993:1999,
                    B = c(b0 = 0.1, C = 0.2), lpt_type = "S")

datt_2018 <- fit_s_multi$datt[fit_s_multi$datt$period == 2018, ]
datt_2019 <- fit_s_multi$datt[fit_s_multi$datt$period == 2019, ]
t1 <- unique(datt_2018$t_index)
t2 <- unique(datt_2019$t_index)

phi_t1 <- t1 * 0.1 + 0.2 * t1 * (t1 + 1) / 2
phi_t2 <- t2 * 0.1 + 0.2 * t2 * (t2 + 1) / 2

width_2018 <- datt_2018$datt_upper[1] - datt_2018$datt_lower[1]
width_2019 <- datt_2019$datt_upper[1] - datt_2019$datt_lower[1]

stopifnot(abs(width_2018 - 2 * phi_t1) < 1e-10)
stopifnot(abs(width_2019 - 2 * phi_t2) < 1e-10)
stopifnot(width_2019 > width_2018)
cat(sprintf("  LPT-S: Phi(t1=%d)=%.3f->w=%.3f, Phi(t2=%d)=%.3f->w=%.3f. OK\n",
            t1, phi_t1, width_2018, t2, phi_t2, width_2019))

# --- Test 3: LPT-S with calibration uses smooth mode ---
fit_s_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "S")
stopifnot(fit_s_cal$calibration$type == "smooth")
stopifnot(!is.null(fit_s_cal$b0))
stopifnot(!is.null(fit_s_cal$C_hat))
stopifnot(fit_s_cal$B_hat == fit_s_cal$b0)
cat(sprintf("  LPT-S calibrated: b0=%.4f, C_hat=%.4f. OK\n",
            fit_s_cal$b0, fit_s_cal$C_hat))

# --- Test 4: ATT width under LPT-S = 2*Phi(t)*d ---
stopifnot(!is.null(fit_s_multi$att))
att_2018 <- fit_s_multi$att[fit_s_multi$att$period == 2018, ]
att_2019 <- fit_s_multi$att[fit_s_multi$att$period == 2019, ]
# width ratio at same dose index: Phi(t2)/Phi(t1)
idx <- 25
w18 <- att_2018$att_upper[idx] - att_2018$att_lower[idx]
w19 <- att_2019$att_upper[idx] - att_2019$att_lower[idx]
expected_ratio <- phi_t2 / phi_t1
stopifnot(abs(w19 / w18 - expected_ratio) < 1e-6)
cat(sprintf("  LPT-S: ATT width ratio = %.3f (expected %.3f). OK\n",
            w19 / w18, expected_ratio))

# --- Test 5: ATT^o width grows with t ---
stopifnot(!is.null(fit_s_multi$att_o))
atto_2018 <- fit_s_multi$att_o[fit_s_multi$att_o$period == 2018, ]
atto_2019 <- fit_s_multi$att_o[fit_s_multi$att_o$period == 2019, ]
w_o18 <- atto_2018$att_o_upper - atto_2018$att_o_lower
w_o19 <- atto_2019$att_o_upper - atto_2019$att_o_lower
stopifnot(w_o19 > w_o18)
cat("  LPT-S: ATT^o width grows with t. OK\n")

# --- Test 6: LPT-S with C=0 gives width 2*t*b0 per period ---
fit_s_c0 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = c(b0 = 0.3, C = 0), lpt_type = "S")
datt_s_2018 <- fit_s_c0$datt[fit_s_c0$datt$period == 2018, ]
phi_t1_c0   <- t1 * 0.3  # C=0, so Phi(t) = t*b0
w_s         <- datt_s_2018$datt_upper[1] - datt_s_2018$datt_lower[1]
stopifnot(abs(w_s - 2 * phi_t1_c0) < 1e-10)
cat(sprintf("  LPT-S with C=0: width = 2*t*b0 = %.3f. OK\n", 2 * phi_t1_c0))

# --- Test 7: b0 and C_hat stored on lpt object; NULL for C/P ---
stopifnot(!is.null(fit_s_multi$b0))
stopifnot(!is.null(fit_s_multi$C_hat))
stopifnot(fit_s_multi$b0 == 0.1)
stopifnot(fit_s_multi$C_hat == 0.2)
fit_c_check <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = 2019, pre_periods = 1993:1999,
                    B = 0.1, lpt_type = "C")
stopifnot(is.null(fit_c_check$b0))
stopifnot(is.null(fit_c_check$C_hat))
cat("  b0/C_hat present for S, NULL for C. OK\n")

cat("lpt_type='S' tests PASSED.\n")
