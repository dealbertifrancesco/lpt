# test-lpt-a.R
cat("Testing lpt() with lpt_type = 'a'...\n")

# --- Test 1: LPT-a with single post-period, B=0 ---
fit_a0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = 0, lpt_type = "a")

stopifnot(inherits(fit_a0, "lpt"))
stopifnot(fit_a0$lpt_type == "a")
stopifnot(fit_a0$B_hat == 0)
stopifnot("t_index" %in% names(fit_a0$datt))
# With B=0, datt bounds should equal lambda_d
datt0 <- fit_a0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
stopifnot(all(abs(datt0$datt_upper - datt0$lambda_d) < 1e-10))
cat("  LPT-a, B=0, single period: OK\n")

# --- Test 2: LPT-a with multiple post-periods, B=0 ---
fit_a_multi <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = c(2018, 2019), pre_periods = 1993:1999,
                    B = 0, lpt_type = "a")

stopifnot(length(fit_a_multi$post_periods) == 2)
stopifnot(all(c(2018, 2019) %in% fit_a_multi$datt$period))
# Check t_index values
t_for_2018 <- unique(fit_a_multi$datt$t_index[fit_a_multi$datt$period == 2018])
t_for_2019 <- unique(fit_a_multi$datt$t_index[fit_a_multi$datt$period == 2019])
stopifnot(length(t_for_2018) == 1, length(t_for_2019) == 1)
stopifnot(t_for_2019 > t_for_2018)  # 2019 is further from ref
cat(sprintf("  LPT-a, multi-period: t(2018)=%d, t(2019)=%d. OK\n", t_for_2018, t_for_2019))

# --- Test 3: LPT-a with B>0 — bounds are CONSTANT width across periods ---
fit_a_b1 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = 0.5, lpt_type = "a")

datt_2018 <- fit_a_b1$datt[fit_a_b1$datt$period == 2018, ]
datt_2019 <- fit_a_b1$datt[fit_a_b1$datt$period == 2019, ]
# Under LPT-a, datt bounds width = 2*B regardless of t
width_2018 <- datt_2018$datt_upper - datt_2018$datt_lower
width_2019 <- datt_2019$datt_upper - datt_2019$datt_lower
stopifnot(all(abs(width_2018 - 1.0) < 1e-10))  # 2*0.5 = 1.0
stopifnot(all(abs(width_2019 - 1.0) < 1e-10))  # same width!
cat("  LPT-a, B=0.5: constant dATT width across periods. OK\n")

# --- Test 4: LPT-a with calibration uses cumulative diffs ---
fit_a_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "a")

stopifnot(!is.null(fit_a_cal$calibration))
stopifnot(fit_a_cal$calibration$type == "cumulative")
stopifnot(fit_a_cal$B_hat > 0)
cat(sprintf("  LPT-a, calibrated: B_hat=%.4f (cumulative). OK\n", fit_a_cal$B_hat))

# --- Test 5: ATT and ATT^o bounds exist when untreated present ---
stopifnot(!is.null(fit_a_b1$att))
stopifnot(!is.null(fit_a_b1$att_o))
# ATT width under LPT-a: 2*B*d (constant across periods)
att_2018 <- fit_a_b1$att[fit_a_b1$att$period == 2018, ]
att_2019 <- fit_a_b1$att[fit_a_b1$att$period == 2019, ]
att_width_2018 <- att_2018$att_upper - att_2018$att_lower  # should be 2*B*d = 2*0.5*d = d
att_width_2019 <- att_2019$att_upper - att_2019$att_lower
# Width should be same for same dose point
stopifnot(all(abs(att_width_2018 - att_width_2019) < 1e-10))
cat("  LPT-a: constant ATT width across periods. OK\n")

cat("lpt_type='a' tests PASSED.\n")
