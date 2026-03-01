#' Local Parallel Trends estimation
#'
#' Estimates identified sets for \eqn{\partial ATT(d|d)/\partial d} and
#' \eqn{ATT(d|d)} under the Local Parallel Trends assumption in
#' continuous difference-in-differences designs.
#'
#' @param data Data frame in long format with unit, time, outcome, and dose columns.
#' @param id_col Character. Column name for unit identifier.
#' @param time_col Character. Column name for time period.
#' @param outcome_col Character. Column name for outcome variable.
#' @param dose_col Character. Column name for dose/treatment intensity.
#' @param post_period Scalar or vector. Identifier(s) for post-treatment period(s).
#'   If a vector, estimation is done separately for each post-period using
#'   long differences from the reference period.
#' @param pre_periods Vector or NULL. Identifiers for pre-treatment periods
#'   (used for calibrating B). If NULL (default), all periods before the earliest
#'   post-period are used.
#' @param B Numeric scalar, numeric vector, or \code{"calibrate"}. Sensitivity
#'   parameter bounding the selection slope \eqn{|\mu'(d)| \leq B}. If
#'   \code{"calibrate"}, estimated from pre-treatment periods. If 0, standard
#'   parallel trends (point identification). If a numeric vector, bounds are
#'   computed for each value (sensitivity analysis). Default: \code{"calibrate"}.
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#'   Default: 50 points over 5th-95th percentile of dose.
#' @param k Integer. Spline basis dimension. Default: 20.
#' @param spline_bs Character. Spline basis type (\code{"cr"} or \code{"tp"}).
#'   Default: \code{"cr"}.
#' @param alpha Numeric. Significance level for pointwise CIs. Default: 0.05.
#'
#' @return An S3 object of class \code{"lpt"} containing:
#'   \describe{
#'     \item{datt}{Data frame with columns \code{period}, \code{d},
#'       \code{lambda_d}, \code{se_lambda}, \code{B}, \code{datt_lower},
#'       \code{datt_upper}, \code{ci_lower}, \code{ci_upper}.
#'       Identified set for \eqn{\partial ATT(d|d)/\partial d}.}
#'     \item{att}{Data frame with ATT level bounds (NULL if no untreated units).
#'       Identified set for \eqn{ATT(d|d)}.}
#'     \item{att_o}{List with overall ATT summary (Corollary 1):
#'       \code{att_o_bin} (binarized DiD), \code{D_bar} (mean dose among treated),
#'       and for each B value: \code{att_o_lower}, \code{att_o_upper}.
#'       NULL if no untreated units.}
#'     \item{B_hat}{The primary sensitivity parameter used.}
#'     \item{B_values}{All B values computed.}
#'     \item{calibration}{Output from \code{\link{calibrate_B}} if calibration
#'       was used.}
#'     \item{slopes}{Named list of \code{\link{estimate_dose_slope}} results,
#'       one per post-period.}
#'     \item{call}{The matched call.}
#'     \item{n}{Number of units in estimation sample.}
#'     \item{has_untreated}{Logical. Whether untreated units (D=0) exist.}
#'     \item{post_periods}{The post-period(s) estimated.}
#'     \item{ref_period}{The reference (last pre-) period used for differencing.}
#'     \item{specifications}{List of all estimation settings.}
#'   }
#'
#' @details
#' The key decomposition is:
#' \deqn{\lambda(d) = \partial ATT(d|d)/\partial d + \mu'(d)}
#' where \eqn{\lambda(d)} is the observable dose-slope and \eqn{\mu'(d)} is
#' the unobservable selection slope bounded by B.
#'
#' Under LPT (\eqn{|\mu'(d)| \leq B}):
#' \itemize{
#'   \item \eqn{IS_{\partial ATT}(d; B) = [\lambda(d) - B, \lambda(d) + B]}
#'   \item \eqn{IS_{ATT}(d; B) = [\Lambda(d) - Bd, \Lambda(d) + Bd]}
#'   \item \eqn{IS_{ATT^o}(B) = [ATT^o_{bin} - B \bar{D}_+, ATT^o_{bin} + B \bar{D}_+]}
#' }
#' where \eqn{\Lambda(d) = E[\Delta Y | D=d] - E[\Delta Y | D=0]} and
#' \eqn{ATT^o_{bin} = E[\Delta Y | D>0] - E[\Delta Y | D=0]}.
#'
#' @examples
#' data(sru)
#' fit <- lpt(sru, "commune", "year", "outcome", "dose",
#'            post_period = 2019, pre_periods = 1993:1999,
#'            B = "calibrate")
#' fit
#'
#' @export
lpt <- function(data, id_col, time_col, outcome_col, dose_col,
                post_period, pre_periods = NULL,
                B = "calibrate", eval_points = NULL,
                k = 20, spline_bs = "cr", alpha = 0.05) {

  # --- Input validation ---
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  required_cols <- c(id_col, time_col, outcome_col, dose_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found in data: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  all_times <- sort(unique(data[[time_col]]))
  bad_periods <- setdiff(post_period, all_times)
  if (length(bad_periods) > 0) {
    stop(sprintf("post_period value(s) not found in '%s': %s",
                 time_col, paste(bad_periods, collapse = ", ")))
  }

  the_call <- match.call()

  # --- Period structure ---
  post_periods <- sort(post_period)

  if (!is.null(pre_periods)) {
    pre_period_set <- sort(pre_periods)
    bad_pre <- setdiff(pre_period_set, all_times)
    if (length(bad_pre) > 0) {
      stop(sprintf("pre_periods value(s) not found in '%s': %s",
                   time_col, paste(bad_pre, collapse = ", ")))
    }
    ref_period <- max(pre_period_set)
  } else {
    pre_period_set <- all_times[all_times < min(post_periods)]
    ref_period <- max(pre_period_set)
  }

  # Compute t (# time increments from ref to each post-period) from data
  all_times_sorted <- sort(unique(data[[time_col]]))
  t_values <- vapply(post_periods, function(pp) {
    sum(all_times_sorted > ref_period & all_times_sorted <= pp)
  }, integer(1L))
  names(t_values) <- as.character(post_periods)

  if (length(pre_period_set) < 1) {
    stop("Need at least one pre-period before the earliest post-period.")
  }

  # --- Extract reference-period data ---
  ref_data <- data[data[[time_col]] == ref_period, ]
  dose_vec <- ref_data[[dose_col]]
  has_untreated <- any(dose_vec == 0)

  # --- Estimate dose slope for each post-period ---
  slopes <- list()
  datt_all <- list()
  att_all <- list()
  att_o_all <- list()

  for (pp in post_periods) {
    post_data <- data[data[[time_col]] == pp, ]

    merged <- merge(
      ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
      post_data[, c(id_col, outcome_col), drop = FALSE],
      by = id_col, suffixes = c("_ref", "_post")
    )

    outcome_ref_col <- paste0(outcome_col, "_ref")
    outcome_post_col <- paste0(outcome_col, "_post")
    delta_y <- merged[[outcome_post_col]] - merged[[outcome_ref_col]]
    wd <- merged[[dose_col]]

    slope_result <- estimate_dose_slope(
      delta_y = delta_y,
      dose = wd,
      eval_points = eval_points,
      k = k,
      spline_bs = spline_bs
    )

    # Lock in eval_points for consistency across periods
    if (is.null(eval_points) && length(slopes) == 0) {
      eval_points <- slope_result$eval_points
    }

    slopes[[as.character(pp)]] <- slope_result
  }

  # --- Handle B ---
  calibration_result <- NULL

  if (is.character(B) && B == "calibrate") {
    if (length(pre_period_set) < 2) {
      stop("B = 'calibrate' requires at least 2 pre-treatment periods. ",
           "Supply B as a numeric value, or provide data with more pre-periods.")
    }

    first_slope <- slopes[[1]]

    calibration_result <- calibrate_B(
      data = data,
      id_col = id_col,
      time_col = time_col,
      outcome_col = outcome_col,
      dose_col = dose_col,
      pre_periods = pre_period_set,
      eval_points = first_slope$eval_points,
      k = k,
      spline_bs = spline_bs
    )
    B_hat <- calibration_result$B_hat
    B_values <- B_hat
    message(sprintf("Calibrated B = %.4f from %d pre-period pair(s).",
                    B_hat, length(pre_period_set) - 1))
  } else if (is.numeric(B)) {
    B_values <- as.numeric(B)
    if (any(!is.finite(B_values)) || any(B_values < 0)) {
      stop("B must be non-negative finite numeric values.")
    }
    B_hat <- max(B_values)
  } else {
    stop("B must be numeric or 'calibrate'.")
  }

  # --- Construct identified sets per period ---
  z_alpha <- stats::qnorm(1 - alpha / 2)

  for (pp in post_periods) {
    pp_char <- as.character(pp)
    sr <- slopes[[pp_char]]
    ep <- sr$eval_points
    lambda_d <- sr$lambda_d
    se_lambda <- sr$se_lambda

    # IS_{dATT}(d; B) = [lambda(d) - B, lambda(d) + B]
    for (b in B_values) {
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period = pp,
        d = ep,
        lambda_d = lambda_d,
        se_lambda = se_lambda,
        B = b,
        datt_lower = lambda_d - b,
        datt_upper = lambda_d + b,
        ci_lower = lambda_d - b - z_alpha * se_lambda,
        ci_upper = lambda_d + b + z_alpha * se_lambda
      )
    }

    # IS_{ATT}(d; B) = [Lambda(d,t) - t*B*d, Lambda(d,t) + t*B*d]
    if (has_untreated) {
      t_val  <- t_values[[pp_char]]
      att_pp <- compute_att_bounds(sr, B_values, dose_vec, period = pp, t = t_val)
      att_all[[length(att_all) + 1]] <- att_pp
    }

    # IS_{ATT^o}(B) = [ATT^o_bin - B * D_bar, ATT^o_bin + B * D_bar]
    # (Corollary 1: overall ATT summary)
    if (has_untreated) {
      # Recompute delta_y for this period to get ATT^o_bin
      post_data <- data[data[[time_col]] == pp, ]
      merged_atto <- merge(
        ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
        post_data[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_ref", "_post")
      )
      dy <- merged_atto[[paste0(outcome_col, "_post")]] -
            merged_atto[[paste0(outcome_col, "_ref")]]
      d_merged <- merged_atto[[dose_col]]

      treated_idx <- d_merged > 0
      untreated_idx <- d_merged == 0
      att_o_bin <- mean(dy[treated_idx]) - mean(dy[untreated_idx])
      D_bar <- mean(d_merged[treated_idx])

      for (b in B_values) {
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period = pp,
          att_o_bin = att_o_bin,
          D_bar = D_bar,
          B = b,
          att_o_lower = att_o_bin - b * D_bar,
          att_o_upper = att_o_bin + b * D_bar
        )
      }
    }
  }

  datt <- do.call(rbind, datt_all)
  att <- if (length(att_all) > 0) do.call(rbind, att_all) else NULL
  att_o <- if (length(att_o_all) > 0) do.call(rbind, att_o_all) else NULL
  rownames(datt) <- NULL
  if (!is.null(att)) rownames(att) <- NULL
  if (!is.null(att_o)) rownames(att_o) <- NULL

  if (!has_untreated) {
    message("No untreated units (dose = 0). ATT level bounds not computed.")
  }

  # --- Assemble output ---
  structure(
    list(
      datt = datt,
      att = att,
      att_o = att_o,
      B_hat = B_hat,
      B_values = B_values,
      calibration = calibration_result,
      slopes = slopes,
      call = the_call,
      n = nrow(ref_data),
      has_untreated = has_untreated,
      post_periods = post_periods,
      ref_period = ref_period,
      t_values = t_values,          # NEW: t multiplier per post-period
      specifications = list(
        k = k, spline_bs = spline_bs, alpha = alpha,
        post_periods = post_periods, ref_period = ref_period
      )
    ),
    class = "lpt"
  )
}
