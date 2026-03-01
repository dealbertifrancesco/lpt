#' Calibrate the sensitivity parameter B from pre-treatment periods
#'
#' Estimates the sensitivity parameter B by computing dose-response slopes
#' from pre-treatment data. Two modes:
#' \itemize{
#'   \item \code{"cumulative"} (for LPT-a): uses cumulative differences
#'     \eqn{Y_0 - Y_{-k}} for all window lengths \eqn{k = 1, \ldots, T_{pre}}.
#'   \item \code{"first_diff"} (for LPT-b): uses first differences
#'     \eqn{Y_s - Y_{s-1}} for consecutive pre-period pairs.
#' }
#'
#' @param data Data frame in long format.
#' @param id_col Character. Unit identifier column.
#' @param time_col Character. Time period column.
#' @param outcome_col Character. Outcome column.
#' @param dose_col Character. Dose column.
#' @param pre_periods Vector. Pre-treatment period identifiers (in order).
#'   Must have at least 2 elements.
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#' @param k Integer. Spline basis dimension. Default: 20.
#' @param spline_bs Character. Spline basis type. Default: \code{"cr"}.
#' @param type Character. Calibration type: \code{"cumulative"} (default,
#'   for LPT-a) or \code{"first_diff"} (for LPT-b).
#'
#' @return A list with components:
#'   \describe{
#'     \item{B_hat}{Numeric. Calibrated sensitivity parameter.}
#'     \item{pre_slopes}{Data frame with columns \code{period_pair},
#'       \code{window_length}, \code{d}, \code{mu_prime_d}, \code{se}.}
#'     \item{sup_by_period}{Named numeric vector of suprema.}
#'     \item{type}{Character. The calibration type used.}
#'   }
#'
#' @examples
#' data(sru)
#' cal <- calibrate_B(sru, "commune", "year", "outcome", "dose",
#'                     pre_periods = 1993:1999)
#' cal$B_hat
#'
#' @export
calibrate_B <- function(data, id_col, time_col, outcome_col, dose_col,
                         pre_periods, eval_points = NULL,
                         k = 20, spline_bs = "cr",
                         type = c("cumulative", "first_diff")) {

  type <- match.arg(type)
  pre_periods <- sort(pre_periods)
  if (length(pre_periods) < 2) {
    stop("Need at least 2 pre-treatment periods to calibrate B.")
  }

  all_slopes <- list()
  sup_vals <- numeric(0)
  pair_labels <- character(0)

  ref <- max(pre_periods)  # baseline = last pre-period

  if (type == "cumulative") {
    # LPT-a calibration: cumulative differences Y_ref - Y_{earlier}
    # for each window length k = 1, ..., T_pre
    earlier_periods <- pre_periods[pre_periods < ref]

    for (t0 in earlier_periods) {
      window_len <- sum(pre_periods > t0 & pre_periods <= ref)
      pair_label <- paste0("k=", window_len, " (", ref, "-", t0, ")")
      pair_labels <- c(pair_labels, pair_label)

      dat_ref <- data[data[[time_col]] == ref, ]
      dat_t0  <- data[data[[time_col]] == t0, ]

      merged <- merge(
        dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
        dat_ref[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_early", "_ref")
      )

      # Cumulative diff: Y_ref - Y_t0
      delta_y  <- merged[[paste0(outcome_col, "_ref")]] -
                   merged[[paste0(outcome_col, "_early")]]
      dose_vec <- merged[[dose_col]]

      slope_result <- estimate_dose_slope(delta_y, dose_vec,
                                           eval_points = eval_points,
                                           k = k, spline_bs = spline_bs)

      # Lock in eval_points after first fit
      if (is.null(eval_points)) eval_points <- slope_result$eval_points

      slopes_df <- data.frame(
        period_pair   = pair_label,
        window_length = window_len,
        d             = slope_result$eval_points,
        mu_prime_d    = slope_result$lambda_d,
        se            = slope_result$se_lambda
      )
      all_slopes[[length(all_slopes) + 1]] <- slopes_df

      sup_val <- max(abs(slope_result$lambda_d))
      sup_vals <- c(sup_vals, sup_val)
    }

  } else {
    # LPT-b calibration: first differences Y_{s} - Y_{s-1} for consecutive pairs
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

      outcome_0_col <- paste0(outcome_col, "_0")
      outcome_1_col <- paste0(outcome_col, "_1")
      delta_y  <- merged[[outcome_1_col]] - merged[[outcome_0_col]]
      dose_vec <- merged[[dose_col]]

      slope_result <- estimate_dose_slope(delta_y, dose_vec,
                                           eval_points = eval_points,
                                           k = k, spline_bs = spline_bs)

      if (is.null(eval_points)) eval_points <- slope_result$eval_points

      slopes_df <- data.frame(
        period_pair   = pair_label,
        window_length = 1L,
        d             = slope_result$eval_points,
        mu_prime_d    = slope_result$lambda_d,
        se            = slope_result$se_lambda
      )
      all_slopes[[length(all_slopes) + 1]] <- slopes_df

      sup_val <- max(abs(slope_result$lambda_d))
      sup_vals <- c(sup_vals, sup_val)
    }
  }

  pre_slopes <- do.call(rbind, all_slopes)
  rownames(pre_slopes) <- NULL
  names(sup_vals) <- pair_labels

  list(
    B_hat         = max(sup_vals),
    pre_slopes    = pre_slopes,
    sup_by_period = sup_vals,
    type          = type
  )
}
