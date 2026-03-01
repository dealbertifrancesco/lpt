# test-backwards-compat.R
cat("Testing backwards compatibility...\n")

# --- Test 1: Single post-period with B=0 produces same results for both types ---
fit_a <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0, lpt_type = "a")
fit_b <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = 2019, pre_periods = 1993:1999,
             B = 0, lpt_type = "b")

# With B=0, LPT-a and LPT-b are identical (no bias)
stopifnot(all(abs(fit_a$datt$datt_lower - fit_b$datt$datt_lower) < 1e-10))
stopifnot(all(abs(fit_a$datt$datt_upper - fit_b$datt$datt_upper) < 1e-10))
cat("  B=0: LPT-a == LPT-b. OK\n")

# --- Test 2: Default lpt_type is "a" ---
fit_default <- lpt(sru, "commune", "year", "outcome", "dose",
                    post_period = 2019, pre_periods = 1993:1999, B = 0)
stopifnot(fit_default$lpt_type == "a")
cat("  Default lpt_type is 'a'. OK\n")

# --- Test 3: t_index >= 1 for single-period-after-ref case ---
t_val <- fit_a$t_index_map[["2019"]]
stopifnot(t_val >= 1)
cat(sprintf("  t_index for 2019 (ref=1999): %d. OK\n", t_val))

# --- Test 4: Numeric B works with both types ---
fit_a_num <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = 0.1, lpt_type = "a")
fit_b_num <- lpt(sru, "commune", "year", "outcome", "dose",
                  post_period = 2019, pre_periods = 1993:1999,
                  B = 0.1, lpt_type = "b")
stopifnot(fit_a_num$B_hat == 0.1)
stopifnot(fit_b_num$B_hat == 0.1)
cat("  Numeric B works with both types. OK\n")

# --- Test 5: lpt object has all expected fields ---
expected_fields <- c("datt", "att", "att_o", "B_hat", "B_values",
                      "calibration", "slopes", "call", "n", "has_untreated",
                      "post_periods", "ref_period", "lpt_type", "t_index_map",
                      "specifications")
missing <- setdiff(expected_fields, names(fit_a))
stopifnot(length(missing) == 0)
cat("  All expected fields present. OK\n")

# --- Test 6: datt has t_index column ---
stopifnot("t_index" %in% names(fit_a$datt))
cat("  datt has t_index column. OK\n")

cat("Backwards compatibility tests PASSED.\n")
