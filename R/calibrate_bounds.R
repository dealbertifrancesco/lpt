#' Calibrate the sensitivity parameters M and B from pre-treatment periods
#'
#' Estimates the level bound \eqn{\hat{M}} and the slope bound \eqn{\hat{B}}
#' from pre-treatment periods. Under no anticipation, pre-period first
#' differences identify the selection function directly:
#' \eqn{E[\Delta Y_s \mid D = d] = \mu_s(d)} for \eqn{s < 0}, and likewise
#' \eqn{\mu'_s(d)} on the interior of the dose support. Under stable
#' selection (post-treatment selection is no worse than the worst
#' pre-treatment selection), the natural data-driven candidates are
#' \deqn{\hat{M} = \max_{s<0} \sup_{d \in D_+} |\hat{\mu}_s(d) - \hat{\mu}_s(0)|,
#'   \qquad
#'   \hat{B} = \max_{s<0} \sup_{d \in (d_l, d_u)} |\hat{\mu}'_s(d)|,}
#' with suprema taken over the evaluation grid.
#'
#' @param data Data frame in long format.
#' @param id_col Character. Unit identifier column.
#' @param time_col Character. Time period column.
#' @param outcome_col Character. Outcome column.
#' @param dose_col Character. Dose column.
#' @param pre_periods Vector. Pre-treatment period identifiers (in order).
#'   Must have at least 2 elements.
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#'   If NULL, 50 points over the 5th-95th percentile of positive doses
#'   (the trimmed interior of the treated support).
#' @param k Integer. Spline basis dimension. Default: 5.
#' @param spline_bs Character. Spline basis type. Default: \code{"cr"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{M_hat}{Numeric. Calibrated level bound
#'       \eqn{\hat{M} = \max_s \sup_d |\hat{\mu}_s(d) - \hat{\mu}_s(0)|}.
#'       \code{NA} if no untreated units exist.}
#'     \item{B_hat}{Numeric. Calibrated slope bound
#'       \eqn{\hat{B} = \max_s \sup_d |\hat{\mu}'_s(d)|}.}
#'     \item{pre_deviations}{Data frame with columns \code{period_pair},
#'       \code{d}, \code{deviation} (\eqn{\hat{\mu}_s(d) - \hat{\mu}_s(0)}).
#'       NULL if no untreated units.}
#'     \item{pre_slopes}{Data frame with columns \code{period_pair},
#'       \code{d}, \code{mu_prime_d}.}
#'     \item{sup_dev_by_period}{Named numeric vector of
#'       \eqn{\sup_d |\hat{\mu}_s(d) - \hat{\mu}_s(0)|} per pre-period pair.
#'       NULL if no untreated units.}
#'     \item{sup_slope_by_period}{Named numeric vector of
#'       \eqn{\sup_d |\hat{\mu}'_s(d)|} per pre-period pair.}
#'     \item{eval_points}{The dose grid used.}
#'   }
#'
#' @details
#' Each consecutive pair of pre-periods contributes one first difference
#' \eqn{\Delta Y_s}. The selection function \eqn{\hat{\mu}_s} is fit by
#' penalized regression spline on \strong{treated units only}: the dose
#' support has a hole at \eqn{(0, d_l)}, so no smooth fit should bridge it.
#' \eqn{\hat{\mu}_s(0)} is the sample mean of \eqn{\Delta Y_s} among
#' untreated units.
#'
#' Note the practical asymmetry: \eqn{\hat{M}} requires only conditional
#' means of pre-period first differences, whereas \eqn{\hat{B}} requires
#' their derivatives, a harder and noisier estimation problem.
#'
#' @examples
#' data(sru)
#' cal <- calibrate_bounds(sru, "commune", "year", "outcome", "dose",
#'                         pre_periods = -7:-1)
#' cal$M_hat
#' cal$B_hat
#'
#' @export
calibrate_bounds <- function(data, id_col, time_col, outcome_col, dose_col,
                             pre_periods, eval_points = NULL,
                             k = 5, spline_bs = "cr") {

  fit_fun <- function(delta_y, dose, ep) {
    sr <- estimate_dose_slope(delta_y, dose, eval_points = ep,
                              k = k, spline_bs = spline_bs)
    list(mean = sr$conditional_mean, deriv = sr$lambda_d)
  }

  calibrate_bounds_engine(
    data = data, id_col = id_col, time_col = time_col,
    outcome_col = outcome_col, dose_col = dose_col,
    pre_periods = pre_periods, eval_points = eval_points,
    fit_fun = fit_fun
  )
}


#' Calibration engine shared across estimation backends
#'
#' Iterates over consecutive pre-period pairs, fits the selection function
#' \eqn{\hat{\mu}_s} on treated units via the supplied fitter, and computes
#' the level deviations (against the untreated sample mean) and slopes that
#' calibrate M and B.
#'
#' @param data,id_col,time_col,outcome_col,dose_col As in
#'   \code{\link{calibrate_bounds}}.
#' @param pre_periods Sorted vector of pre-period identifiers (>= 2).
#' @param eval_points Numeric vector or NULL. If NULL, a default grid over
#'   the 5th-95th percentile of positive doses in the first pair is used.
#' @param fit_fun Function \code{(delta_y, dose, eval_points)} returning a
#'   list with elements \code{mean} (fitted conditional mean on the grid)
#'   and \code{deriv} (its derivative on the grid).
#'
#' @return Same structure as \code{\link{calibrate_bounds}}.
#'
#' @keywords internal
calibrate_bounds_engine <- function(data, id_col, time_col, outcome_col,
                                    dose_col, pre_periods, eval_points,
                                    fit_fun) {

  pre_periods <- sort(pre_periods)
  if (length(pre_periods) < 2) {
    stop("Need at least 2 pre-treatment periods to calibrate M and B.")
  }

  slope_rows <- list()
  dev_rows <- list()
  sup_slope <- numeric(0)
  sup_dev <- numeric(0)
  pair_labels <- character(0)
  any_untreated <- FALSE

  for (j in seq_len(length(pre_periods) - 1)) {
    t0 <- pre_periods[j]
    t1 <- pre_periods[j + 1]
    pair_label <- paste0(t0, "-", t1)
    pair_labels <- c(pair_labels, pair_label)

    dat_t0 <- data[data[[time_col]] == t0, ]
    dat_t1 <- data[data[[time_col]] == t1, ]

    merged <- merge(
      dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
      dat_t1[, c(id_col, outcome_col), drop = FALSE],
      by = id_col, suffixes = c("_0", "_1")
    )

    delta_y <- merged[[paste0(outcome_col, "_1")]] -
               merged[[paste0(outcome_col, "_0")]]
    dose_vec <- merged[[dose_col]]
    treated <- dose_vec > 0
    untreated <- dose_vec == 0

    if (is.null(eval_points)) {
      d_pos <- dose_vec[treated]
      rng <- stats::quantile(d_pos, probs = c(0.05, 0.95))
      eval_points <- seq(rng[[1]], rng[[2]], length.out = 50)
    }

    ft <- fit_fun(delta_y[treated], dose_vec[treated], eval_points)

    slope_rows[[j]] <- data.frame(
      period_pair = pair_label,
      d = eval_points,
      mu_prime_d = as.numeric(ft$deriv)
    )
    sup_slope <- c(sup_slope, max(abs(ft$deriv)))

    if (any(untreated)) {
      any_untreated <- TRUE
      mu_at_zero <- mean(delta_y[untreated])
      deviation <- as.numeric(ft$mean) - mu_at_zero
      dev_rows[[j]] <- data.frame(
        period_pair = pair_label,
        d = eval_points,
        deviation = deviation
      )
      sup_dev <- c(sup_dev, max(abs(deviation)))
    } else {
      sup_dev <- c(sup_dev, NA_real_)
    }
  }

  names(sup_slope) <- pair_labels
  names(sup_dev) <- pair_labels

  pre_slopes <- do.call(rbind, slope_rows)
  rownames(pre_slopes) <- NULL
  pre_deviations <- NULL
  if (any_untreated) {
    pre_deviations <- do.call(rbind, dev_rows)
    rownames(pre_deviations) <- NULL
  }

  list(
    M_hat = if (any_untreated) max(sup_dev, na.rm = TRUE) else NA_real_,
    B_hat = max(sup_slope),
    pre_deviations = pre_deviations,
    pre_slopes = pre_slopes,
    sup_dev_by_period = if (any_untreated) sup_dev else NULL,
    sup_slope_by_period = sup_slope,
    eval_points = eval_points
  )
}
