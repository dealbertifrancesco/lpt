#' Calibrate the sensitivity parameter B from pre-treatment periods
#'
#' Estimates B by fitting dose-response slopes across pre-treatment differences.
#' Two modes:
#' \itemize{
#'   \item \code{"cumulative"} (for LPT-a): uses ALL cumulative differences
#'     \eqn{Y_{t_1} - Y_{t_0}} for every ordered pair \eqn{t_0 < t_1} in
#'     pre-treatment. This is C(n_pre, 2) pairs total.
#'   \item \code{"first_diff"} (for LPT-b): uses consecutive first differences
#'     \eqn{Y_s - Y_{s-1}} for adjacent pre-period pairs only.
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
#' @param k Integer. Spline basis dimension. Default: 5.
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
#' @export
calibrate_B <- function(data, id_col, time_col, outcome_col, dose_col,
                         pre_periods, eval_points = NULL,
                         k = 5, spline_bs = "cr",
                         type = c("cumulative", "first_diff")) {

  type <- match.arg(type)
  pre_periods <- sort(pre_periods)
  if (length(pre_periods) < 2) {
    stop("Need at least 2 pre-treatment periods to calibrate B.")
  }

  all_slopes <- list()
  sup_vals <- numeric(0)
  pair_labels <- character(0)

  if (type == "cumulative") {
    # LPT-a calibration: ALL ordered pairs (t0, t1) with t0 < t1 in pre_periods.
    # For n_pre periods this gives C(n_pre, 2) pairs.
    # Each pair estimates the slope of the cumulative trend mu(d, t1-t0),
    # and we take the supremum over all pairs.
    n <- length(pre_periods)
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1L, n)) {
        t0 <- pre_periods[i]
        t1 <- pre_periods[j]
        window_len <- j - i   # number of period steps between t0 and t1
        pair_label <- paste0(t0, "-", t1)
        pair_labels <- c(pair_labels, pair_label)

        dat_t0 <- data[data[[time_col]] == t0, ]
        dat_t1 <- data[data[[time_col]] == t1, ]

        merged <- merge(
          dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
          dat_t1[, c(id_col, outcome_col), drop = FALSE],
          by = id_col, suffixes = c("_early", "_late")
        )

        # Cumulative diff: Y_t1 - Y_t0
        delta_y  <- merged[[paste0(outcome_col, "_late")]] -
                    merged[[paste0(outcome_col, "_early")]]
        dose_vec <- merged[[dose_col]]

        slope_result <- estimate_dose_slope(delta_y, dose_vec,
                                             eval_points = eval_points,
                                             k = k, spline_bs = spline_bs)

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

      delta_y  <- merged[[paste0(outcome_col, "_1")]] -
                  merged[[paste0(outcome_col, "_0")]]
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
