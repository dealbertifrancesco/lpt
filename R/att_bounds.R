#' Compute ATT level bounds
#'
#' Computes the identified set for \eqn{ATT(d,t|d)} under Local Parallel Trends:
#' \eqn{IS_{ATT}(d, t; B) = [\Lambda_t(d) - (t+1)Bd, \Lambda_t(d) + (t+1)Bd]}.
#'
#' @param slope_result Output from \code{\link{estimate_dose_slope}}.
#' @param B_values Numeric vector of sensitivity parameter values.
#' @param dose Numeric vector of observed doses.
#' @param period Scalar. The post-period label for this estimate.
#' @param horizon Integer. Horizon index (0 = first post-period). The identified
#'   set width is multiplied by \code{(horizon + 1)} to account for accumulated
#'   selection drift across post-treatment periods.
#'
#' @return A data frame with columns: \code{period}, \code{horizon}, \code{d},
#'   \code{Lambda_d}, \code{B}, \code{att_lower}, \code{att_upper}.
#'
#' @keywords internal
compute_att_bounds <- function(slope_result, B_values, dose, period = NA,
                               horizon = 0L) {
  ep <- slope_result$eval_points
  gam_fit <- slope_result$gam_fit

  cm_at_ep <- as.numeric(stats::predict(gam_fit,
                                         newdata = data.frame(dose = ep)))
  cm_at_zero <- as.numeric(stats::predict(gam_fit,
                                           newdata = data.frame(dose = 0)))
  Lambda_d <- cm_at_ep - cm_at_zero

  mult <- horizon + 1L

  att_list <- lapply(B_values, function(b) {
    data.frame(
      period = period,
      horizon = horizon,
      d = ep,
      Lambda_d = Lambda_d,
      B = b,
      att_lower = Lambda_d - mult * b * ep,
      att_upper = Lambda_d + mult * b * ep
    )
  })
  do.call(rbind, att_list)
}
