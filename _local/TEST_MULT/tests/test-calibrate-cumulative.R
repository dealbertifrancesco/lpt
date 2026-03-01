# test-calibrate-cumulative.R
cat("Testing calibrate_B cumulative mode...\n")

# --- Test 1: first_diff mode — consecutive pairs only ---
cal_fd <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                       pre_periods = 1993:1999, type = "first_diff")
stopifnot(is.numeric(cal_fd$B_hat), cal_fd$B_hat > 0)
stopifnot(cal_fd$type == "first_diff")
stopifnot(is.data.frame(cal_fd$pre_slopes))
stopifnot("period_pair" %in% names(cal_fd$pre_slopes))
# 6 consecutive pairs for 7 pre-periods
n_pairs_fd <- length(unique(cal_fd$pre_slopes$period_pair))
stopifnot(n_pairs_fd == 6)
cat(sprintf("  first_diff: B_hat = %.4f, %d pairs. OK\n", cal_fd$B_hat, n_pairs_fd))

# --- Test 2: cumulative mode uses ALL ordered pairs ---
cal_cum <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                        pre_periods = 1993:1999, type = "cumulative")
stopifnot(is.numeric(cal_cum$B_hat), cal_cum$B_hat > 0)
stopifnot(cal_cum$type == "cumulative")
stopifnot(is.data.frame(cal_cum$pre_slopes))
# C(7,2) = 21 pairs for 7 pre-periods (every ordered pair t0 < t1)
n_windows <- length(unique(cal_cum$pre_slopes$period_pair))
stopifnot(n_windows == 21)
# Must have window_length column (1 = adjacent, 2 = two steps apart, etc.)
stopifnot("window_length" %in% names(cal_cum$pre_slopes))
# window_length ranges from 1 to 6
stopifnot(min(cal_cum$pre_slopes$window_length) == 1)
stopifnot(max(cal_cum$pre_slopes$window_length) == 6)
cat(sprintf("  cumulative: B_hat = %.4f, %d pairs. OK\n", cal_cum$B_hat, n_windows))

# --- Test 3: B comparison (informational) ---
cat(sprintf("  B_hat comparison: cumulative=%.4f vs first_diff=%.4f\n",
            cal_cum$B_hat, cal_fd$B_hat))

# --- Test 4: default type is "cumulative" ---
cal_default <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                            pre_periods = 1993:1999)
stopifnot(cal_default$type == "cumulative")
cat("  default type is cumulative. OK\n")

# --- Test 5: with only 2 pre-periods, both modes give 1 pair and agree ---
cal_2_fd <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                          pre_periods = c(1998, 1999), type = "first_diff")
cal_2_cum <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                           pre_periods = c(1998, 1999), type = "cumulative")
# C(2,2) = 1 pair; 1 consecutive pair
stopifnot(length(unique(cal_2_fd$pre_slopes$period_pair)) == 1)
stopifnot(length(unique(cal_2_cum$pre_slopes$period_pair)) == 1)
# With 2 periods, cumulative(k=1) IS the first diff
stopifnot(abs(cal_2_fd$B_hat - cal_2_cum$B_hat) < 0.01)
cat("  2-period case: both modes agree. OK\n")

# --- Test 6: with 3 pre-periods, cumulative gives C(3,2)=3 pairs ---
cal_3_cum <- calibrate_B(sru, "commune", "year", "outcome", "dose",
                           pre_periods = c(1997, 1998, 1999), type = "cumulative")
stopifnot(length(unique(cal_3_cum$pre_slopes$period_pair)) == 3)
cat("  3-period case: C(3,2)=3 cumulative pairs. OK\n")

cat("calibrate_B cumulative tests PASSED.\n")
