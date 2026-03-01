# test-lpt-c.R
cat("Testing lpt() with lpt_type = 'C'...\n")

fit_c0 <- lpt(sru, "commune", "year", "outcome", "dose",
              post_period = 2019, pre_periods = 1993:1999,
              B = 0, lpt_type = "C")
stopifnot(inherits(fit_c0, "lpt"))
stopifnot(fit_c0$lpt_type == "C")
stopifnot(fit_c0$B_hat == 0)
stopifnot("t_index" %in% names(fit_c0$datt))
datt0 <- fit_c0$datt
stopifnot(all(abs(datt0$datt_lower - datt0$lambda_d) < 1e-10))
stopifnot(all(abs(datt0$datt_upper - datt0$lambda_d) < 1e-10))
cat("  LPT-C, B=0, single period: OK\n")

fit_c_multi <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = c(2018, 2019), pre_periods = 1993:1999,
                    B = 0, lpt_type = "C")
stopifnot(length(fit_c_multi$post_periods) == 2)
stopifnot(all(c(2018, 2019) %in% fit_c_multi$datt$period))
t_for_2018 <- unique(fit_c_multi$datt$t_index[fit_c_multi$datt$period == 2018])
t_for_2019 <- unique(fit_c_multi$datt$t_index[fit_c_multi$datt$period == 2019])
stopifnot(length(t_for_2018) == 1, length(t_for_2019) == 1)
stopifnot(t_for_2019 > t_for_2018)
cat(sprintf("  LPT-C, multi-period: t(2018)=%d, t(2019)=%d. OK\n", t_for_2018, t_for_2019))

# LPT-C: dATT width = 2*B, constant across all periods
fit_c_b1 <- lpt(sru, "commune", "year", "outcome", "dose",
                 post_period = c(2018, 2019), pre_periods = 1993:1999,
                 B = 0.5, lpt_type = "C")
datt_2018 <- fit_c_b1$datt[fit_c_b1$datt$period == 2018, ]
datt_2019 <- fit_c_b1$datt[fit_c_b1$datt$period == 2019, ]
width_2018 <- datt_2018$datt_upper - datt_2018$datt_lower
width_2019 <- datt_2019$datt_upper - datt_2019$datt_lower
stopifnot(all(abs(width_2018 - 1.0) < 1e-10))
stopifnot(all(abs(width_2019 - 1.0) < 1e-10))
cat("  LPT-C, B=0.5: constant dATT width across periods. OK\n")

fit_c_cal <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = "calibrate", lpt_type = "C")
stopifnot(!is.null(fit_c_cal$calibration))
stopifnot(fit_c_cal$calibration$type == "cumulative")
stopifnot(fit_c_cal$B_hat > 0)
cat(sprintf("  LPT-C, calibrated: B_hat=%.4f (cumulative). OK\n", fit_c_cal$B_hat))

stopifnot(!is.null(fit_c_b1$att))
stopifnot(!is.null(fit_c_b1$att_o))
att_2018 <- fit_c_b1$att[fit_c_b1$att$period == 2018, ]
att_2019 <- fit_c_b1$att[fit_c_b1$att$period == 2019, ]
att_width_2018 <- att_2018$att_upper - att_2018$att_lower
att_width_2019 <- att_2019$att_upper - att_2019$att_lower
stopifnot(all(abs(att_width_2018 - att_width_2019) < 1e-10))
cat("  LPT-C: constant ATT width across periods. OK\n")

cat("lpt_type='C' tests PASSED.\n")
