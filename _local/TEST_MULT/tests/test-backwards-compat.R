# test-backwards-compat.R
cat("Testing backwards compatibility...\n")

# --- Test 1: B=0 produces identical results for C and P ---
fit_c <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0, lpt_type = "C")
fit_p <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0, lpt_type = "P")
stopifnot(all(abs(fit_c$datt$datt_lower - fit_p$datt$datt_lower) < 1e-10))
stopifnot(all(abs(fit_c$datt$datt_upper - fit_p$datt$datt_upper) < 1e-10))
cat("  B=0: LPT-C == LPT-P. OK\n")

# --- Test 2: Default lpt_type is "C" ---
fit_default <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = 2019, pre_periods = 1993:1999, B = 0)
stopifnot(fit_default$lpt_type == "C")
cat("  Default lpt_type is 'C'. OK\n")

# --- Test 3: t_index >= 1 for single post-period after ref ---
t_val <- fit_c$t_index_map[["2019"]]
stopifnot(t_val >= 1)
cat(sprintf("  t_index for 2019 (ref=1999): %d. OK\n", t_val))

# --- Test 4: Numeric B works with C and P ---
fit_c_num <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = 0.1, lpt_type = "C")
fit_p_num <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = 0.1, lpt_type = "P")
stopifnot(fit_c_num$B_hat == 0.1)
stopifnot(fit_p_num$B_hat == 0.1)
cat("  Numeric B works with C and P. OK\n")

# --- Test 5: lpt object has all expected fields ---
expected_fields <- c("datt", "att", "att_o", "B_hat", "B_values",
                      "b0", "C_hat", "calibration", "slopes", "call",
                      "n", "has_untreated", "post_periods", "ref_period",
                      "lpt_type", "t_index_map", "specifications")
missing <- setdiff(expected_fields, names(fit_c))
stopifnot(length(missing) == 0)
cat("  All expected fields present. OK\n")

# --- Test 6: b0 and C_hat are NULL for C and P ---
stopifnot(is.null(fit_c$b0))
stopifnot(is.null(fit_c$C_hat))
stopifnot(is.null(fit_p$b0))
stopifnot(is.null(fit_p$C_hat))
cat("  b0/C_hat NULL for C and P. OK\n")

# --- Test 7: datt has t_index column ---
stopifnot("t_index" %in% names(fit_c$datt))
cat("  datt has t_index column. OK\n")

cat("Backwards compatibility tests PASSED.\n")
