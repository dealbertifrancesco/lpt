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
  n_post <- length(x$post_periods)
  post_str <- if (n_post == 1) {
    as.character(x$post_periods)
  } else {
    sprintf("%d periods (%s)", n_post,
            paste(x$post_periods, collapse = ", "))
  }
  m_str <- if (is.finite(x$M_hat)) {
    sprintf("%.4f (%s)", x$M_hat, x$M_source)
  } else {
    "NA (no untreated units)"
  }
  cat(sprintf(
    "Local Parallel Trends (lpt) | n = %d | M = %s | B = %.4f (%s) | Post: %s\n",
    x$n, m_str, x$B_hat, x$B_source, post_str
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
  cat("  Local Parallel Trends Summary\n")
  rule()

  # --- Sensitivity parameters ---
  cat("\n  Sensitivity parameters:\n")
  if (is.finite(object$M_hat)) {
    cat(sprintf("    M (level bound, |mu_t(d) - mu_t(0)| <= M): %.4f (%s)\n",
                object$M_hat, object$M_source))
  } else {
    cat("    M (level bound): NA (no untreated units)\n")
  }
  cat(sprintf("    B (slope bound, |mu_t'(d)| <= B):          %.4f (%s)\n",
              object$B_hat, object$B_source))
  cat("  M = 0 (resp. B = 0) is standard parallel trends for levels\n")
  cat("  (resp. slopes): point identification.\n")

  # --- Calibration details ---
  if (!is.null(object$calibration)) {
    cal <- object$calibration
    if (!is.null(cal$sup_dev_by_period)) {
      cat("\n  Pre-period sup|mu_s(d) - mu_s(0)| by pair (calibrates M):\n")
      for (nm in names(cal$sup_dev_by_period)) {
        cat(sprintf("    %s: %.4f\n", nm, cal$sup_dev_by_period[nm]))
      }
    }
    cat("\n  Pre-period sup|mu_s'(d)| by pair (calibrates B):\n")
    for (nm in names(cal$sup_slope_by_period)) {
      cat(sprintf("    %s: %.4f\n", nm, cal$sup_slope_by_period[nm]))
    }
  }

  # --- ATT^o summary per period ---
  if (!is.null(object$att_o) && is.finite(object$M_hat)) {
    for (pp in object$post_periods) {
      atto_pp <- object$att_o[object$att_o$period == pp &
                                object$att_o$M == object$M_hat, ]
      if (nrow(atto_pp) == 0) next
      h <- atto_pp$horizon[1]

      cat(sprintf("\n  --- ATT^o (period %s, horizon %d, drift mult = %d) ---\n",
                  pp, h, h + 1L))
      cat(sprintf("  ATT^o_t (binary DiD):   %.4f\n", atto_pp$att_o_bin))
      cat(sprintf("  IS(M = %.4f):           [%.4f, %.4f]\n",
                  object$M_hat, atto_pp$att_o_lower, atto_pp$att_o_upper))
    }
  }

  # --- ATT^o aggregated ---
  if (!is.null(object$att_o_agg) && is.finite(object$M_hat)) {
    agg <- object$att_o_agg[object$att_o_agg$M == object$M_hat, ]
    if (nrow(agg) > 0) {
      cat(sprintf("\n  --- ATT^o time-aggregated (%d periods) ---\n",
                  agg$n_periods[1]))
      cat(sprintf("  Lambda-bar^agg: %.4f\n", agg$Lambda_agg))
      cat(sprintf("  IS(M = %.4f):   [%.4f, %.4f]\n",
                  object$M_hat, agg$att_o_agg_lower, agg$att_o_agg_upper))
    }
  }

  # --- dATT results at selected quantiles, per period ---
  for (pp in object$post_periods) {
    datt_pp <- object$datt[object$datt$period == pp &
                             object$datt$B == object$B_hat, ]
    if (nrow(datt_pp) == 0) next
    h <- datt_pp$horizon[1]

    cat(sprintf("\n  --- dATT identified sets (period %s, horizon %d; B) ---\n",
                pp, h))
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
  if (!is.null(object$att) && is.finite(object$M_hat)) {
    for (pp in object$post_periods) {
      att_pp <- object$att[object$att$period == pp &
                             object$att$M == object$M_hat, ]
      if (nrow(att_pp) == 0) next
      h <- att_pp$horizon[1]

      cat(sprintf("\n  --- ATT identified sets (period %s, horizon %d; M) ---\n",
                  pp, h))
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
  meth <- object$specifications$method
  if (is.null(meth) || meth == "gam") {
    cat(sprintf(
      "\n  n = %d | Method: GAM | Spline: %s (k = %d)\n",
      object$n,
      object$specifications$spline_bs,
      object$specifications$k
    ))
  } else if (meth == "contdid") {
    cd_args <- object$specifications$contdid_args
    cat(sprintf(
      "\n  n = %d | Method: contdid | knots = %s, degree = %s\n",
      object$n,
      if (!is.null(cd_args$num_knots)) cd_args$num_knots else "1",
      if (!is.null(cd_args$degree)) cd_args$degree else "3"
    ))
  } else if (meth == "npiv") {
    np_args <- object$specifications$npiv_args
    j_seg <- if (!is.null(np_args$J.x.segments)) np_args$J.x.segments else "2"
    j_deg <- if (!is.null(np_args$J.x.degree)) np_args$J.x.degree else "3"
    cat(sprintf(
      "\n  n = %d | Method: npiv | J.segments = %s, degree = %s\n",
      object$n, j_seg, j_deg
    ))
  } else {
    kr_args <- object$specifications$kernel_args
    bw_str <- if (!is.null(kr_args$bw) && is.numeric(kr_args$bw)) {
      sprintf("%.4f", kr_args$bw)
    } else {
      if (!is.null(kr_args$bw)) kr_args$bw else "cv.ls"
    }
    reg <- if (!is.null(kr_args$regtype)) kr_args$regtype else "ll"
    ker <- if (!is.null(kr_args$ckertype)) kr_args$ckertype else "gaussian"
    # Show actual selected bandwidths if available
    bw_actual <- NULL
    if (!is.null(object$kernel_fits)) {
      bws <- sapply(object$kernel_fits, function(x) x$bw)
      bw_actual <- sprintf("%.4f", mean(bws))
    }
    bw_display <- if (!is.null(bw_actual)) bw_actual else bw_str
    cat(sprintf(
      "\n  n = %d | Method: kernel (%s) | bw = %s, kernel = %s\n",
      object$n, reg, bw_display, ker
    ))
  }
  rule()

  invisible(object)
}
