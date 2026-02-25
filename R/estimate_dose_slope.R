#' Estimate the dose-response slope via penalized splines
#'
#' Fits a GAM to estimate \eqn{E[\Delta Y \mid D = d]} and its derivative
#' \eqn{\lambda(d) = \frac{\partial}{\partial d} E[\Delta Y \mid D = d]}
#' using cubic regression splines with REML smoothing parameter selection.
#'
#' @param delta_y Numeric vector. First differences \eqn{\Delta Y_i}.
#' @param dose Numeric vector. Doses \eqn{D_i}. Same length as \code{delta_y}.
#' @param eval_points Numeric vector or NULL. Dose values at which to evaluate.
#'   Default: 50 evenly spaced points from 5th to 95th percentile of dose.
#' @param k Integer. Number of spline basis functions (default: 20).
#'   Reduced automatically if too large for the sample size.
#' @param spline_bs Character. Basis type: \code{"cr"} (cubic regression, default)
#'   or \code{"tp"} (thin plate).
#'
#' @return A list of class \code{"dose_slope"} with components:
#'   \describe{
#'     \item{eval_points}{Dose grid used for evaluation.}
#'     \item{conditional_mean}{\eqn{\hat{E}[\Delta Y \mid D = d]} at each evaluation point.}
#'     \item{lambda_d}{\eqn{\hat{\lambda}(d)} — estimated derivative at each point.}
#'     \item{se_lambda}{Pointwise standard error of the derivative.}
#'     \item{gam_fit}{The fitted \code{mgcv::gam} object.}
#'   }
#'
#' @details
#' The derivative is computed via finite differencing of the prediction matrix
#' (lpmatrix). Standard errors use the delta method. See the package vignette
#' for theoretical details.
#'
#' @examples
#' data(sru)
#' ref <- sru[sru$year == 1999, ]
#' post <- sru[sru$year == 2019, ]
#' m <- merge(ref[, c("commune", "outcome", "dose")],
#'            post[, c("commune", "outcome")],
#'            by = "commune", suffixes = c("_ref", "_post"))
#' m$delta_y <- m$outcome_post - m$outcome_ref
#' result <- estimate_dose_slope(m$delta_y, m$dose)
#' plot(result$eval_points, result$lambda_d, type = "l")
#'
#' @export
estimate_dose_slope <- function(delta_y, dose, eval_points = NULL,
                                 k = 20, spline_bs = "cr") {
  # --- Input validation ---
  if (!is.numeric(delta_y)) stop("delta_y must be numeric.")
  if (!is.numeric(dose)) stop("dose must be numeric.")
  if (length(delta_y) != length(dose)) {
    stop("delta_y and dose must have the same length.")
  }

  # Remove NAs
  nas <- is.na(delta_y) | is.na(dose)
  if (any(nas)) {
    warning(sprintf("Removing %d observations with NA values.", sum(nas)))
    delta_y <- delta_y[!nas]
    dose <- dose[!nas]
  }
  n <- length(delta_y)
  if (n < 10) stop("Need at least 10 observations.")

  # Validate spline basis
  spline_bs <- match.arg(spline_bs, c("cr", "tp"))

  # Adjust k if too large relative to n
  if (k >= n / 5) {
    k_new <- max(4L, as.integer(floor(n / 5)))
    warning(sprintf("k = %d too large for n = %d. Reducing to k = %d.", k, n, k_new))
    k <- k_new
  }

  # --- Default eval_points ---
  if (is.null(eval_points)) {
    dose_range <- stats::quantile(dose, probs = c(0.05, 0.95))
    eval_points <- seq(dose_range[[1]], dose_range[[2]], length.out = 50)
  }

  # --- Fit GAM ---
  fit_data <- data.frame(delta_y = delta_y, dose = dose)
  gam_fit <- mgcv::gam(delta_y ~ s(dose, bs = spline_bs, k = k),
                        data = fit_data, method = "REML")

  # --- Check EDF ---
  edf <- summary(gam_fit)$s.table[, "edf"]
  if (edf > k - 2) {
    warning(sprintf(
      "Effective degrees of freedom (%.1f) close to k-1 = %d. Consider increasing k.",
      edf, k - 1
    ))
  }

  # --- Conditional mean ---
  newdata <- data.frame(dose = eval_points)
  conditional_mean <- as.numeric(stats::predict(gam_fit, newdata = newdata,
                                                 type = "response"))

  # --- Derivative via finite differencing of lpmatrix ---
  deriv_result <- get_derivative(gam_fit, eval_points)

  structure(
    list(
      eval_points = eval_points,
      conditional_mean = conditional_mean,
      lambda_d = deriv_result$derivative,
      se_lambda = deriv_result$se,
      gam_fit = gam_fit
    ),
    class = "dose_slope"
  )
}


#' Compute derivative of a GAM fit via finite differencing of the lpmatrix
#'
#' @param gam_fit A fitted \code{mgcv::gam} object.
#' @param eval_points Numeric vector. Points at which to evaluate the derivative.
#' @param eps Numeric. Step size for finite differencing (default: 1e-7).
#'
#' @return A list with \code{derivative} and \code{se} (both numeric vectors).
#'
#' @keywords internal
get_derivative <- function(gam_fit, eval_points, eps = 1e-7) {
  newdata0 <- data.frame(dose = eval_points)
  newdata1 <- data.frame(dose = eval_points + eps)

  X0 <- stats::predict(gam_fit, newdata = newdata0, type = "lpmatrix")
  X1 <- stats::predict(gam_fit, newdata = newdata1, type = "lpmatrix")

  Xd <- (X1 - X0) / eps

  beta_hat <- stats::coef(gam_fit)
  V <- stats::vcov(gam_fit)

  deriv <- as.numeric(Xd %*% beta_hat)
  se_deriv <- sqrt(pmax(rowSums((Xd %*% V) * Xd), 0))

  list(derivative = deriv, se = se_deriv)
}
