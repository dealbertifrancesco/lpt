# test-summary-plot.R
cat("Testing summary and plot methods...\n")

# --- Test 1: print works for both types ---
fit_a <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "a")
fit_b <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "b")

cat("\n--- print LPT-a ---\n")
print(fit_a)
cat("\n--- print LPT-b ---\n")
print(fit_b)
cat("  print works for both types. OK\n")

# --- Test 2: summary works for both types ---
cat("\n--- summary LPT-a ---\n")
summary(fit_a)
cat("\n--- summary LPT-b ---\n")
summary(fit_b)
cat("  summary works for both types. OK\n")

# --- Test 3: plots work (if ggplot2 available) ---
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # datt plot
  p1 <- plot(fit_a, type = "datt")
  p2 <- plot(fit_b, type = "datt")
  cat("  plot(type='datt') works. OK\n")

  # att plot
  p3 <- plot(fit_a, type = "att")
  p4 <- plot(fit_b, type = "att")
  cat("  plot(type='att') works. OK\n")

  # pretrends plot
  p5 <- plot(fit_a, type = "pretrends")
  p6 <- plot(fit_b, type = "pretrends")
  cat("  plot(type='pretrends') works. OK\n")

  # sensitivity plot
  p7 <- plot(fit_a, type = "sensitivity", estimand = "att_o")
  p8 <- plot(fit_b, type = "sensitivity", estimand = "att_o")
  cat("  plot(type='sensitivity', estimand='att_o') works. OK\n")

  d_median <- stats::median(fit_a$datt$d)
  p9  <- plot(fit_a, type = "sensitivity", estimand = "datt", d0 = d_median)
  p10 <- plot(fit_b, type = "sensitivity", estimand = "datt", d0 = d_median)
  cat("  plot(type='sensitivity', estimand='datt') works. OK\n")
} else {
  cat("  ggplot2 not available, skipping plot tests.\n")
}

cat("Summary and plot tests PASSED.\n")
