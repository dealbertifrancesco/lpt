#' Compute ATT level bounds
#'
#' Computes the identified set for \eqn{ATT(d,t|d)} under the level
#' restriction (bounded trend deviations, \eqn{|\mu_t(d) - \mu_t(0)| \le M}):
#' \eqn{IS_{ATT}(d, t; M) = [\Lambda_t(d) - (t+1)M, \Lambda_t(d) + (t+1)M]}.
#' The set widens linearly in the horizon: each post period contributes
#' another M of counterfactual drift.
#'
#' @param eval_points Numeric vector. Dose grid.
#' @param Lambda_d Numeric vector. DiD estimand \eqn{\hat\Lambda_t(d) =
#'   \hat{f}_t(d) - \hat{f}_t(0)} at each grid point, where
#'   \eqn{\hat{f}_t(0)} is the untreated sample mean of the long difference.
#' @param M_values Numeric vector of level-bound values.
#' @param period Scalar. The post-period label for this estimate.
#' @param horizon Integer. Horizon index (0 = first post-period).
#'
#' @return A data frame with columns: \code{period}, \code{horizon}, \code{d},
#'   \code{Lambda_d}, \code{M}, \code{att_lower}, \code{att_upper}.
#'
#' @keywords internal
compute_att_bounds <- function(eval_points, Lambda_d, M_values, period = NA,
                               horizon = 0L) {
  mult <- horizon + 1L

  att_list <- lapply(M_values, function(m) {
    data.frame(
      period = period,
      horizon = horizon,
      d = eval_points,
      Lambda_d = Lambda_d,
      M = m,
      att_lower = Lambda_d - mult * m,
      att_upper = Lambda_d + mult * m
    )
  })
  do.call(rbind, att_list)
}


#' Resolve the sensitivity parameters M and B
#'
#' Turns user input (\code{"calibrate"} or numeric) into concrete values,
#' drawing on a calibration result when requested. The level bound M is only
#' meaningful when untreated units exist (level targets need the untreated
#' reference); without them \code{M_hat} is \code{NA} and \code{M_values}
#' is empty.
#'
#' @param M,B User-supplied specification: \code{"calibrate"} or numeric
#'   scalar/vector.
#' @param calibration Output of \code{\link{calibrate_bounds}} (or an
#'   equivalent backend calibration), or NULL when neither parameter is
#'   calibrated.
#' @param has_untreated Logical. Whether untreated (dose = 0) units exist.
#'
#' @return A list with \code{M_hat}, \code{M_values}, \code{M_source},
#'   \code{B_hat}, \code{B_values}, \code{B_source}.
#'
#' @keywords internal
resolve_bounds <- function(M, B, calibration, has_untreated) {

  if (is.character(B) && length(B) == 1 && B == "calibrate") {
    if (is.null(calibration)) {
      stop("Internal error: B = 'calibrate' but no calibration was run.")
    }
    B_hat <- calibration$B_hat
    B_values <- B_hat
    B_source <- "calibrated"
  } else if (is.numeric(B)) {
    B_values <- as.numeric(B)
    if (any(!is.finite(B_values)) || any(B_values < 0)) {
      stop("B must be non-negative finite numeric values.")
    }
    B_hat <- max(B_values)
    B_source <- "user-supplied"
  } else {
    stop("B must be numeric or 'calibrate'.")
  }

  if (is.character(M) && length(M) == 1 && M == "calibrate") {
    if (!has_untreated) {
      M_hat <- NA_real_
      M_values <- numeric(0)
    } else {
      if (is.null(calibration)) {
        stop("Internal error: M = 'calibrate' but no calibration was run.")
      }
      M_hat <- calibration$M_hat
      M_values <- M_hat
    }
    M_source <- "calibrated"
  } else if (is.numeric(M)) {
    M_values <- as.numeric(M)
    if (any(!is.finite(M_values)) || any(M_values < 0)) {
      stop("M must be non-negative finite numeric values.")
    }
    M_hat <- max(M_values)
    M_source <- "user-supplied"
  } else {
    stop("M must be numeric or 'calibrate'.")
  }

  list(M_hat = M_hat, M_values = M_values, M_source = M_source,
       B_hat = B_hat, B_values = B_values, B_source = B_source)
}
