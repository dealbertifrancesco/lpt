#' Calibrate the sensitivity parameter B from pre-treatment periods
#'
#' Two modes for LPT-C and LPT-P; one mode for LPT-S.
#' \itemize{
#'   \item \code{"cumulative"} (LPT-C): ALL C(n,2) ordered pairs.
#'   \item \code{"first_diff"} (LPT-P): consecutive first differences only.
#'   \item \code{"smooth"} (LPT-S): same fits as first_diff, plus extracts
#'     b_0 (curvature at last pre-period) and C_hat (max change in curvature).
#' }
#'
#' @param data Data frame in long format.
#' @param id_col Character. Unit identifier column.
#' @param time_col Character. Time period column.
#' @param outcome_col Character. Outcome column.
#' @param dose_col Character. Dose column.
#' @param pre_periods Vector. Pre-treatment period identifiers. At least 2.
#' @param eval_points Numeric vector or NULL.
#' @param k Integer. Spline basis dimension. Default: 5.
#' @param spline_bs Character. Spline basis type. Default: \code{"cr"}.
#' @param type Character. One of \code{"cumulative"} (default, LPT-C),
#'   \code{"first_diff"} (LPT-P), \code{"smooth"} (LPT-S).
#'
#' @return A list with:
#'   \describe{
#'     \item{B_hat}{Calibrated B (or b_0 for smooth).}
#'     \item{C_hat}{NULL for cumulative/first_diff; Ĉ for smooth.}
#'     \item{b_s_sequence}{NULL for cumulative/first_diff; named b_s vector for smooth.}
#'     \item{pre_slopes}{Data frame of fitted slopes per pair.}
#'     \item{sup_by_period}{Named vector of sup|mu'(d)| per pair.}
#'     \item{type}{The calibration type used.}
#'   }
#' @export
calibrate_B <- function(data, id_col, time_col, outcome_col, dose_col,
                         pre_periods, eval_points = NULL,
                         k = 5, spline_bs = "cr",
                         type = c("cumulative", "first_diff", "smooth")) {

  type <- match.arg(type)
  pre_periods <- sort(pre_periods)
  if (length(pre_periods) < 2) {
    stop("Need at least 2 pre-treatment periods to calibrate B.")
  }

  all_slopes  <- list()
  sup_vals    <- numeric(0)
  pair_labels <- character(0)

  if (type == "cumulative") {
    # LPT-C: ALL C(n,2) ordered pairs — estimates slope of cumulative mu(d,k)
    n <- length(pre_periods)
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1L, n)) {
        t0 <- pre_periods[i]; t1 <- pre_periods[j]
        window_len <- j - i
        pair_label <- paste0(t0, "-", t1)
        pair_labels <- c(pair_labels, pair_label)

        dat_t0 <- data[data[[time_col]] == t0, ]
        dat_t1 <- data[data[[time_col]] == t1, ]
        merged <- merge(
          dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
          dat_t1[, c(id_col, outcome_col), drop = FALSE],
          by = id_col, suffixes = c("_early", "_late")
        )
        delta_y  <- merged[[paste0(outcome_col, "_late")]] -
                    merged[[paste0(outcome_col, "_early")]]
        dose_vec <- merged[[dose_col]]

        sr <- estimate_dose_slope(delta_y, dose_vec,
                                  eval_points = eval_points,
                                  k = k, spline_bs = spline_bs)
        if (is.null(eval_points)) eval_points <- sr$eval_points

        all_slopes[[length(all_slopes) + 1]] <- data.frame(
          period_pair = pair_label, window_length = window_len,
          d = sr$eval_points, mu_prime_d = sr$lambda_d, se = sr$se_lambda
        )
        sup_vals <- c(sup_vals, max(abs(sr$lambda_d)))
      }
    }
    names(sup_vals) <- pair_labels
    pre_slopes <- do.call(rbind, all_slopes); rownames(pre_slopes) <- NULL
    return(list(B_hat = max(sup_vals), C_hat = NULL, b_s_sequence = NULL,
                pre_slopes = pre_slopes, sup_by_period = sup_vals, type = type))

  } else {
    # "first_diff" and "smooth": consecutive first-difference pairs
    b_s_vec <- numeric(0)

    for (j in seq_len(length(pre_periods) - 1)) {
      t0 <- pre_periods[j]; t1 <- pre_periods[j + 1]
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

      sr <- estimate_dose_slope(delta_y, dose_vec,
                                eval_points = eval_points,
                                k = k, spline_bs = spline_bs)
      if (is.null(eval_points)) eval_points <- sr$eval_points

      b_s <- max(abs(sr$lambda_d))
      b_s_vec <- c(b_s_vec, b_s)

      all_slopes[[length(all_slopes) + 1]] <- data.frame(
        period_pair = pair_label, window_length = 1L,
        d = sr$eval_points, mu_prime_d = sr$lambda_d, se = sr$se_lambda
      )
      sup_vals <- c(sup_vals, b_s)
    }
    names(sup_vals) <- pair_labels
    names(b_s_vec)  <- pair_labels
    pre_slopes <- do.call(rbind, all_slopes); rownames(pre_slopes) <- NULL

    if (type == "first_diff") {
      return(list(B_hat = max(sup_vals), C_hat = NULL, b_s_sequence = NULL,
                  pre_slopes = pre_slopes, sup_by_period = sup_vals, type = type))
    }

    # type == "smooth" (LPT-S)
    # b0 = b_s at the last pre-period pair (s=0, closest to treatment)
    b0    <- b_s_vec[length(b_s_vec)]
    # C_hat = max change in max-curvature across consecutive pre-periods
    C_hat <- if (length(b_s_vec) >= 2) max(abs(diff(b_s_vec))) else 0

    return(list(B_hat = b0, C_hat = C_hat, b_s_sequence = b_s_vec,
                pre_slopes = pre_slopes, sup_by_period = sup_vals, type = type))
  }
}
