#' Print an lpt object
#'
#' @param x An object of class \code{"lpt"}.
#' @param ... Additional arguments (unused).
#'
#' @return \code{x}, invisibly.
#'
#' @method print lpt
#' @export
print.lpt <- function(x, ...) {
  b_source <- if (!is.null(x$calibration)) "calibrated" else "user-supplied"
  n_post <- length(x$post_periods)
  post_str <- if (n_post == 1) {
    as.character(x$post_periods)
  } else {
    sprintf("%d periods (%s)", n_post,
            paste(x$post_periods, collapse = ", "))
  }
  cat(sprintf(
    "Lipschitz Parallel Trends (lpt) | n = %d | B = %.4f (%s) | Post: %s\n",
    x$n, x$B_hat, b_source, post_str
  ))
  invisible(x)
}


#' Summarize an lpt object
#'
#' @param object An object of class \code{"lpt"}.
#' @param ... Additional arguments (unused).
#'
#' @return \code{object}, invisibly.
#'
#' @method summary lpt
#' @export
summary.lpt <- function(object, ...) {
  rule <- function(char = "=", width = 55) {
    cat(paste(rep(char, width), collapse = ""), "\n")
  }

  rule()
  cat("  Lipschitz Parallel Trends Summary\n")
  rule()

  # --- B information ---
  if (!is.null(object$calibration)) {
    cat(sprintf("\n  Sensitivity parameter B: %.4f (calibrated)\n",
                object$B_hat))
    cat("  Pre-period sup|mu'(d)| by pair:\n")
    for (nm in names(object$calibration$sup_by_period)) {
      cat(sprintf("    %s: %.4f\n", nm,
                  object$calibration$sup_by_period[nm]))
    }
  } else {
    cat(sprintf("\n  Sensitivity parameter B: %.4f (user-supplied)\n",
                object$B_hat))
  }
  cat("  B = 0 is standard parallel trends (point identification).\n")

  # --- ATT^o summary (Corollary 1) ---
  if (!is.null(object$att_o)) {
    for (pp in object$post_periods) {
      atto_pp <- object$att_o[object$att_o$period == pp &
                                object$att_o$B == object$B_hat, ]
      if (nrow(atto_pp) == 0) next

      cat(sprintf("\n  --- ATT^o overall summary (period %s) ---\n", pp))
      cat(sprintf("  ATT^o_bin (binarized DiD): %.4f\n", atto_pp$att_o_bin))
      cat(sprintf("  Mean dose among treated:   %.4f\n", atto_pp$D_bar))
      cat(sprintf("  IS_{ATT^o}(B = %.4f):      [%.4f, %.4f]\n",
                  object$B_hat, atto_pp$att_o_lower, atto_pp$att_o_upper))
    }
  }

  # --- dATT results at selected quantiles, per period ---
  for (pp in object$post_periods) {
    datt_pp <- object$datt[object$datt$period == pp &
                             object$datt$B == object$B_hat, ]
    if (nrow(datt_pp) == 0) next

    cat(sprintf("\n  --- dATT identified sets (period %s) ---\n", pp))
    dose_vals <- datt_pp$d
    qtiles <- stats::quantile(dose_vals, probs = c(0.25, 0.5, 0.75))

    cat(sprintf("  %-10s %-12s %-22s %-10s\n",
                "Dose", "lambda(d)", "[dATT lower, upper]", "Excl. 0?"))
    for (qt in qtiles) {
      idx <- which.min(abs(dose_vals - qt))
      row <- datt_pp[idx, ]
      excludes_zero <- (row$datt_lower > 0) | (row$datt_upper < 0)
      cat(sprintf("  %-10.3f %-12.4f [%-8.4f, %-8.4f]  %-10s\n",
                  row$d, row$lambda_d, row$datt_lower, row$datt_upper,
                  ifelse(excludes_zero, "Yes", "No")))
    }
  }

  # --- ATT identified sets at dose quartiles ---
  if (!is.null(object$att)) {
    for (pp in object$post_periods) {
      att_pp <- object$att[object$att$period == pp &
                             object$att$B == object$B_hat, ]
      if (nrow(att_pp) == 0) next

      cat(sprintf("\n  --- ATT identified sets (period %s) ---\n", pp))
      dose_vals <- att_pp$d
      qtiles <- stats::quantile(dose_vals, probs = c(0.25, 0.5, 0.75))

      cat(sprintf("  %-10s %-12s %-22s %-10s\n",
                  "Dose", "Lambda(d)", "[ATT lower, upper]", "Excl. 0?"))
      for (qt in qtiles) {
        idx <- which.min(abs(dose_vals - qt))
        row <- att_pp[idx, ]
        excludes_zero <- (row$att_lower > 0) | (row$att_upper < 0)
        cat(sprintf("  %-10.3f %-12.4f [%-8.4f, %-8.4f]  %-10s\n",
                    row$d, row$Lambda_d, row$att_lower, row$att_upper,
                    ifelse(excludes_zero, "Yes", "No")))
      }
    }
  }

  # --- Footer ---
  cat(sprintf(
    "\n  n = %d | Spline: %s (k = %d) | alpha = %.2f\n",
    object$n,
    object$specifications$spline_bs,
    object$specifications$k,
    object$specifications$alpha
  ))
  rule()

  invisible(object)
}
