#' Compute ATT level bounds (multi-period, with optional time restrictions)
#'
#' Computes the identified set for \eqn{ATT_t(d|d)} under LPT with optional
#' time-direction restrictions (RM-Time or SD-Time). Applies the t multiplier.
#'
#' @param slope_result Output from \code{\link{estimate_dose_slope}}.
#' @param B_values Numeric vector of sensitivity parameter values.
#' @param dose Numeric vector of observed doses.
#' @param period Scalar. The post-period label for this estimate.
#' @param t Integer. Number of time increments from ref to this post-period.
#' @param time_restriction Character. One of "none", "rm", "sd".
#' @param M_bar Numeric. RM-Time parameter (required if time_restriction = "rm").
#' @param delta_star_pre_d Numeric vector. Max |delta_tilde_s(d)| over pre-pairs,
#'   evaluated at eval_points (required for time_restriction = "rm").
#' @param B_t Numeric. SD-Time Lipschitz parameter (required if time_restriction = "sd").
#' @param delta_tilde_0_d Numeric vector. Last pre-period increment delta_tilde_0(d),
#'   evaluated at eval_points (required for time_restriction = "sd").
#'
#' @return A data frame with columns: \code{period}, \code{d}, \code{Lambda_d},
#'   \code{B}, \code{t}, \code{att_lower}, \code{att_upper},
#'   \code{binding_lower}, \code{binding_upper}.
#'
#' @keywords internal
compute_att_bounds <- function(slope_result, B_values, dose, period = NA, t = 1L,
                               time_restriction = "none",
                               M_bar = NULL, delta_star_pre_d = NULL,
                               B_t = NULL, delta_tilde_0_d = NULL) {
  ep       <- slope_result$eval_points
  gam_fit  <- slope_result$gam_fit

  # Lambda(d) = E[DeltaY | D=d] - E[DeltaY | D=0]
  cm_at_ep   <- as.numeric(stats::predict(gam_fit, newdata = data.frame(dose = ep)))
  cm_at_zero <- as.numeric(stats::predict(gam_fit, newdata = data.frame(dose = 0)))
  Lambda_d   <- cm_at_ep - cm_at_zero

  att_list <- lapply(B_values, function(b) {
    # --- LPT bounds (with t multiplier) ---
    lpt_hw     <- t * b * ep  # half-width under LPT
    lpt_lower  <- Lambda_d - lpt_hw
    lpt_upper  <- Lambda_d + lpt_hw

    # --- Time-restriction bounds ---
    if (time_restriction == "rm") {
      if (is.null(delta_star_pre_d))
        stop("delta_star_pre_d required for time_restriction = 'rm'.")
      rm_hw    <- t * M_bar * delta_star_pre_d
      rm_lower <- Lambda_d - rm_hw
      rm_upper <- Lambda_d + rm_hw
    } else if (time_restriction == "sd") {
      if (is.null(delta_tilde_0_d) || is.null(B_t))
        stop("delta_tilde_0_d and B_t required for time_restriction = 'sd'.")
      sd_hw     <- B_t * t * (t + 1) / 2
      sd_center <- Lambda_d - t * delta_tilde_0_d
      sd_lower  <- sd_center - sd_hw
      sd_upper  <- sd_center + sd_hw
    }

    # --- Determine final bounds and which restriction binds ---
    if (time_restriction == "none") {
      final_lower <- lpt_lower
      final_upper <- lpt_upper
      bind_lower  <- rep("lpt", length(ep))
      bind_upper  <- rep("lpt", length(ep))
    } else if (time_restriction == "rm") {
      if (b > 0) {
        final_lower <- pmax(lpt_lower, rm_lower)
        final_upper <- pmin(lpt_upper, rm_upper)
      } else {
        final_lower <- rm_lower
        final_upper <- rm_upper
      }
      bind_lower <- ifelse(b > 0 & lpt_lower >= rm_lower, "lpt",
                    ifelse(b > 0 & lpt_lower < rm_lower, "rm", "rm"))
      bind_upper <- ifelse(b > 0 & lpt_upper <= rm_upper, "lpt",
                    ifelse(b > 0 & lpt_upper > rm_upper, "rm", "rm"))
    } else {  # "sd"
      if (b > 0) {
        final_lower <- pmax(lpt_lower, sd_lower)
        final_upper <- pmin(lpt_upper, sd_upper)
      } else {
        final_lower <- sd_lower
        final_upper <- sd_upper
      }
      bind_lower <- ifelse(b > 0 & lpt_lower >= sd_lower, "lpt", "sd")
      bind_upper <- ifelse(b > 0 & lpt_upper <= sd_upper, "lpt", "sd")
    }

    # Check for empty intersection
    empty_idx <- final_lower > final_upper
    if (any(empty_idx)) {
      warning(sprintf("%d dose point(s) have empty intersection of LPT and time bounds.",
                      sum(empty_idx)))
      final_lower[empty_idx] <- NA_real_
      final_upper[empty_idx] <- NA_real_
    }

    data.frame(
      period        = period,
      d             = ep,
      Lambda_d      = Lambda_d,
      B             = b,
      t             = t,
      att_lower     = final_lower,
      att_upper     = final_upper,
      binding_lower = bind_lower,
      binding_upper = bind_upper
    )
  })
  do.call(rbind, att_list)
}
