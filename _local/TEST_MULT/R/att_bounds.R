#' Compute ATT level bounds
#'
#' Computes the identified set for \eqn{ATT_t(d|d)} under Local Parallel Trends:
#' \itemize{
#'   \item LPT-a: \eqn{IS_{ATT}(d; B) = [\Lambda(d,t) - Bd, \Lambda(d,t) + Bd]}
#'   \item LPT-b: \eqn{IS_{ATT}(d; B_d) = [\Lambda(d,t) - tB_d d, \Lambda(d,t) + tB_d d]}
#' }
#'
#' @param slope_result Output from \code{\link{estimate_dose_slope}}.
#' @param B_values Numeric vector of sensitivity parameter values.
#' @param dose Numeric vector of observed doses.
#' @param period Scalar. The post-period label for this estimate.
#' @param t_multiplier Numeric. Multiplier for the bias term.
#'   For LPT-a: 1 (constant across horizons).
#'   For LPT-b: t (the number of periods from baseline).
#'
#' @return A data frame with columns: \code{period}, \code{d},
#'   \code{Lambda_d}, \code{B}, \code{att_lower}, \code{att_upper}.
#'
#' @keywords internal
compute_att_bounds <- function(slope_result, B_values, dose, period = NA,
                                t_multiplier = 1) {
  ep       <- slope_result$eval_points
  gam_fit  <- slope_result$gam_fit

  # Lambda(d, t) = E[DeltaY | D=d] - E[DeltaY | D=0]
  cm_at_ep   <- as.numeric(stats::predict(gam_fit,
                                           newdata = data.frame(dose = ep)))
  cm_at_zero <- as.numeric(stats::predict(gam_fit,
                                           newdata = data.frame(dose = 0)))
  Lambda_d   <- cm_at_ep - cm_at_zero

  att_list <- lapply(B_values, function(b) {
    # LPT-a: bias_hw = b * d          (t_multiplier = 1)
    # LPT-b: bias_hw = t * b * d      (t_multiplier = t)
    bias_hw <- t_multiplier * b * ep
    data.frame(
      period    = period,
      d         = ep,
      Lambda_d  = Lambda_d,
      B         = b,
      att_lower = Lambda_d - bias_hw,
      att_upper = Lambda_d + bias_hw
    )
  })
  do.call(rbind, att_list)
}
