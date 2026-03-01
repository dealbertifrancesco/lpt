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
  lpt_label <- if (!is.null(x$lpt_type)) {
    sprintf("LPT-%s", x$lpt_type)
  } else {
    "LPT"
  }
  n_post <- length(x$post_periods)
  post_str <- if (n_post == 1) {
    as.character(x$post_periods)
  } else {
    sprintf("%d periods (%s)", n_post,
            paste(x$post_periods, collapse = ", "))
  }
  cat(sprintf(
    "%s | n = %d | B = %.4f (%s) | Post: %s\n",
    lpt_label, x$n, x$B_hat, b_source, post_str
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
  rule <- function(char = "=", width = 60) {
    cat(paste(rep(char, width), collapse = ""), "\n")
  }

  rule()
  lpt_label <- if (!is.null(object$lpt_type)) {
    sprintf("Local Parallel Trends Summary (LPT-%s)", object$lpt_type)
  } else {
    "Local Parallel Trends Summary"
  }
  cat(sprintf("  %s\n", lpt_label))
  rule()

  # --- LPT type explanation ---
  if (!is.null(object$lpt_type)) {
    if (object$lpt_type == "a") {
      cat("\n  LPT-a: smoothness on cumulative trend.")
      cat("\n  Identified set width is CONSTANT across post-periods.\n")
    } else {
      cat("\n  LPT-b: smoothness on period-specific increments.")
      cat("\n  Identified set width GROWS linearly with horizon t.\n")
    }
  }

  # --- B information ---
  if (!is.null(object$calibration)) {
    cal_label <- if (!is.null(object$calibration$type)) {
      sprintf("calibrated, %s", object$calibration$type)
    } else {
      "calibrated"
    }
    cat(sprintf("\n  Sensitivity parameter B: %.4f (%s)\n",
                object$B_hat, cal_label))
    cat("  Pre-period sup|mu'(d)| by pair/window:\n")
    for (nm in names(object$calibration$sup_by_period)) {
      cat(sprintf("    %s: %.4f\n", nm,
                  object$calibration$sup_by_period[nm]))
    }
  } else {
    cat(sprintf("\n  Sensitivity parameter B: %.4f (user-supplied)\n",
                object$B_hat))
  }
  cat("  B = 0 is standard parallel trends (point identification).\n")

  # --- t_index per post-period ---
  if (!is.null(object$t_index_map)) {
    cat("\n  Post-period horizons:\n")
    for (nm in names(object$t_index_map)) {
      t_val <- object$t_index_map[[nm]]
      if (!is.null(object$lpt_type) && object$lpt_type == "b") {
        cat(sprintf("    period %s: t = %d (bias multiplier = %d)\n", nm, t_val, t_val))
      } else {
        cat(sprintf("    period %s: t = %d\n", nm, t_val))
      }
    }
  }

  # --- ATT^o summary ---
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

  # --- dATT results at selected quantiles ---
  for (pp in object$post_periods) {
    datt_pp <- object$datt[object$datt$period == pp &
                             object$datt$B == object$B_hat, ]
    if (nrow(datt_pp) == 0) next

    t_str <- if ("t_index" %in% names(datt_pp)) {
      sprintf(", t=%d", datt_pp$t_index[1])
    } else {
      ""
    }
    cat(sprintf("\n  --- dATT identified sets (period %s%s) ---\n", pp, t_str))
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

  # --- ATT identified sets ---
  if (!is.null(object$att)) {
    for (pp in object$post_periods) {
      att_pp <- object$att[object$att$period == pp &
                             object$att$B == object$B_hat, ]
      if (nrow(att_pp) == 0) next

      t_str <- if ("t_index" %in% names(att_pp)) {
        sprintf(", t=%d", att_pp$t_index[1])
      } else {
        ""
      }
      cat(sprintf("\n  --- ATT identified sets (period %s%s) ---\n", pp, t_str))
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
