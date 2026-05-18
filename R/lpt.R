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
#' @param method Character. Estimation method: \code{"gam"} (default, penalized
#'   splines via \code{mgcv}) or \code{"contdid"} (B-splines via the
#'   \code{contdid} package; Callaway, Goodman-Bacon & Sant'Anna 2024).
#' @param contdid_args Named list. Additional arguments passed to
#'   \code{contdid::cont_did()} when \code{method = "contdid"}. Common options:
#'   \code{num_knots} (default 1), \code{degree} (default 3),
#'   \code{biters} (bootstrap iterations, default 500),
#'   \code{control_group} (default \code{"notyettreated"}).
#'   Ignored when \code{method = "gam"}.
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#'   Default: 50 points over 5th-95th percentile of dose. Ignored when
#'   \code{method = "contdid"} (determined by contdid internally).
#' @param k Integer. Spline basis dimension for GAM method. Default: 5.
#'   Ignored when \code{method = "contdid"}.
#' @param spline_bs Character. Spline basis type (\code{"cr"} or \code{"tp"}).
#'   Default: \code{"cr"}. Ignored when \code{method = "contdid"}.
#' @return An S3 object of class \code{"lpt"} containing:
#'   \describe{
#'     \item{datt}{Data frame with columns \code{period}, \code{horizon},
#'       \code{d}, \code{lambda_d}, \code{B}, \code{datt_lower},
#'       \code{datt_upper}. Identified set width is \code{2*(horizon+1)*B}.}
#'     \item{att}{Data frame with ATT level bounds (NULL if no untreated units).
#'       Includes \code{horizon} column; width is \code{2*(horizon+1)*B*d}.}
#'     \item{att_o}{Data frame with per-period overall ATT summary.
#'       Includes \code{horizon} column; width is
#'       \code{2*(horizon+1)*B*D_bar}. NULL if no untreated units.}
#'     \item{att_o_agg}{Data frame with time-aggregated overall ATT across
#'       all requested post-periods. NULL if single post-period or no
#'       untreated units.}
#'     \item{pre_att_o}{Data frame with pre-period binary DiD estimates
#'       (columns \code{period}, \code{att_o_bin}). Used by the event study
#'       plot. NULL if no untreated units or fewer than 2 pre-periods.}
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
#'     \item{contdid_fit}{Raw \code{pte_results} object from \code{contdid}
#'       (NULL when \code{method = "gam"}).}
#'     \item{specifications}{List of all estimation settings.}
#'   }
#'
#' @references
#' Callaway B, Goodman-Bacon A, Sant'Anna PHC (2024).
#' \dQuote{Difference-in-differences with a continuous treatment.}
#' \emph{National Bureau of Economic Research}.
#'
#' Wood SN (2017).
#' \emph{Generalized Additive Models: An Introduction with R} (2nd ed.).
#' Chapman and Hall/CRC.
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
#'            post_period = 0:5, pre_periods = -7:-1,
#'            B = "calibrate")
#' fit
#'
#' @export
lpt <- function(data, id_col, time_col, outcome_col, dose_col,
                post_period, pre_periods = NULL,
                B = "calibrate", method = c("gam", "contdid"),
                contdid_args = list(), eval_points = NULL,
                k = 5, spline_bs = "cr") {

  method <- match.arg(method)

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
  min_post <- min(post_periods)

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

  if (length(pre_period_set) < 1) {
    stop("Need at least one pre-period before the earliest post-period.")
  }

  # --- Extract reference-period data ---
  ref_data <- data[data[[time_col]] == ref_period, ]
  dose_vec <- ref_data[[dose_col]]
  has_untreated <- any(dose_vec == 0)

  contdid_fit <- NULL

  # ==================================================================
  #  Estimation: branch on method
  # ==================================================================

  if (method == "gam") {

    # --- GAM: estimate dose slope for each post-period ---
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
    for (pp in post_periods) {
      pp_char <- as.character(pp)
      sr <- slopes[[pp_char]]
      ep <- sr$eval_points
      lambda_d <- sr$lambda_d
      horizon_t <- as.integer(pp - min_post)
      mult <- horizon_t + 1L

      # IS_{dATT}(d, t; B) = [lambda_t(d) - (t+1)B, lambda_t(d) + (t+1)B]
      for (b in B_values) {
        datt_all[[length(datt_all) + 1]] <- data.frame(
          period = pp,
          horizon = horizon_t,
          d = ep,
          lambda_d = lambda_d,
          B = b,
          datt_lower = lambda_d - mult * b,
          datt_upper = lambda_d + mult * b
        )
      }

      # IS_{ATT}(d, t; B) = [Lambda_t(d) - (t+1)Bd, Lambda_t(d) + (t+1)Bd]
      if (has_untreated) {
        att_pp <- compute_att_bounds(sr, B_values, dose_vec,
                                      period = pp, horizon = horizon_t)
        att_all[[length(att_all) + 1]] <- att_pp
      }

      # IS_{ATT^o_t}(B) = [att_o_bin - (t+1)*B*D_bar, att_o_bin + (t+1)*B*D_bar]
      if (has_untreated) {
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
            horizon = horizon_t,
            att_o_bin = att_o_bin,
            D_bar = D_bar,
            B = b,
            att_o_lower = att_o_bin - mult * b * D_bar,
            att_o_upper = att_o_bin + mult * b * D_bar
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

  } else {

    # --- contdid: delegate to backend ---
    cd <- run_contdid(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      post_periods = post_periods, pre_period_set = pre_period_set,
      min_post = min_post, has_untreated = has_untreated,
      dose_vec = dose_vec, B = B, contdid_args = contdid_args
    )
    slopes             <- cd$slopes
    calibration_result <- cd$calibration
    B_hat              <- cd$B_hat
    B_values           <- cd$B_values
    datt               <- cd$datt
    att                <- cd$att
    att_o              <- cd$att_o
    contdid_fit        <- cd$contdid_fit
  }

  # ==================================================================
  #  Common: ATT^o aggregation, pre-period event study, output
  # ==================================================================

  # --- Time-aggregated ATT^0 (Eq. 31 generalized) ---
  att_o_agg <- NULL
  if (has_untreated && length(post_periods) > 1 && !is.null(att_o)) {
    horizons <- as.integer(post_periods - min_post)
    K <- length(post_periods)
    avg_mult <- mean(horizons + 1L)

    for (b in B_values) {
      att_o_b <- att_o[att_o$B == b, ]
      Lambda_agg <- mean(att_o_b$att_o_bin)
      D_bar_agg <- att_o_b$D_bar[1]
      half_width <- avg_mult * b * D_bar_agg

      att_o_agg <- rbind(att_o_agg, data.frame(
        Lambda_agg = Lambda_agg,
        D_bar = D_bar_agg,
        n_periods = K,
        B = b,
        att_o_agg_lower = Lambda_agg - half_width,
        att_o_agg_upper = Lambda_agg + half_width
      ))
    }
    rownames(att_o_agg) <- NULL
  }

  if (!has_untreated) {
    message("No untreated units (dose = 0). ATT level bounds not computed.")
  }

  # --- Pre-period binary DiD (for event study plot) ---
  pre_att_o <- NULL
  if (has_untreated && length(pre_period_set) >= 2) {
    pre_att_o_list <- list()
    pre_no_ref <- pre_period_set[pre_period_set != ref_period]

    for (s in pre_no_ref) {
      s_data <- data[data[[time_col]] == s, ]
      merged_pre <- merge(
        ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
        s_data[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_ref", "_s")
      )
      dy_pre <- merged_pre[[paste0(outcome_col, "_s")]] -
                merged_pre[[paste0(outcome_col, "_ref")]]
      d_pre <- merged_pre[[dose_col]]
      treated_pre <- d_pre > 0
      untreated_pre <- d_pre == 0

      if (sum(treated_pre) > 0 && sum(untreated_pre) > 0) {
        beta_s <- mean(dy_pre[treated_pre]) - mean(dy_pre[untreated_pre])
        pre_att_o_list[[length(pre_att_o_list) + 1]] <- data.frame(
          period = s,
          att_o_bin = beta_s
        )
      }
    }

    if (length(pre_att_o_list) > 0) {
      pre_att_o <- do.call(rbind, pre_att_o_list)
      rownames(pre_att_o) <- NULL
    }
  }

  # --- Assemble output ---
  structure(
    list(
      datt = datt,
      att = att,
      att_o = att_o,
      att_o_agg = att_o_agg,
      pre_att_o = pre_att_o,
      B_hat = B_hat,
      B_values = B_values,
      calibration = calibration_result,
      slopes = slopes,
      contdid_fit = contdid_fit,
      call = the_call,
      n = nrow(ref_data),
      has_untreated = has_untreated,
      post_periods = post_periods,
      ref_period = ref_period,
      specifications = list(
        method = method,
        k = if (method == "gam") k else NULL,
        spline_bs = if (method == "gam") spline_bs else NULL,
        contdid_args = if (method == "contdid") contdid_args else NULL,
        post_periods = post_periods, ref_period = ref_period
      )
    ),
    class = "lpt"
  )
}
