# test-summary-plot.R
cat("Testing summary and plot methods...\n")

fit_c <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "C")
fit_p <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "P")
fit_s <- lpt(sru, "commune", "year", "outcome", "dose",
             post_period = c(2018, 2019), pre_periods = 1993:1999,
             B = "calibrate", lpt_type = "S")

cat("\n--- print LPT-C ---\n"); print(fit_c)
cat("\n--- print LPT-P ---\n"); print(fit_p)
cat("\n--- print LPT-S ---\n"); print(fit_s)
cat("  print works for all three types. OK\n")

cat("\n--- summary LPT-C ---\n"); summary(fit_c)
cat("\n--- summary LPT-P ---\n"); summary(fit_p)
cat("\n--- summary LPT-S ---\n"); summary(fit_s)
cat("  summary works for all three types. OK\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p1 <- plot(fit_c, type = "datt")
  p2 <- plot(fit_p, type = "datt")
  p3 <- plot(fit_s, type = "datt")
  cat("  plot(type='datt') works for C/P/S. OK\n")

  p4 <- plot(fit_c, type = "att")
  p5 <- plot(fit_p, type = "att")
  p6 <- plot(fit_s, type = "att")
  cat("  plot(type='att') works for C/P/S. OK\n")

  p7 <- plot(fit_c, type = "pretrends")
  p8 <- plot(fit_p, type = "pretrends")
  p9 <- plot(fit_s, type = "pretrends")
  cat("  plot(type='pretrends') works for C/P/S. OK\n")

  p10 <- plot(fit_c, type = "sensitivity", estimand = "att_o")
  p11 <- plot(fit_p, type = "sensitivity", estimand = "att_o")
  p12 <- plot(fit_s, type = "sensitivity", estimand = "att_o")
  cat("  plot(type='sensitivity', estimand='att_o') works for C/P/S. OK\n")

  d_median <- stats::median(fit_c$datt$d)
  p13 <- plot(fit_c, type = "sensitivity", estimand = "datt", d0 = d_median)
  p14 <- plot(fit_p, type = "sensitivity", estimand = "datt", d0 = d_median)
  p15 <- plot(fit_s, type = "sensitivity", estimand = "datt", d0 = d_median)
  cat("  plot(type='sensitivity', estimand='datt') works for C/P/S. OK\n")
} else {
  cat("  ggplot2 not available, skipping plot tests.\n")
}

cat("Summary and plot tests PASSED.\n")
