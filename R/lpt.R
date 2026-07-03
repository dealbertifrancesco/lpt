#' Local Parallel Trends estimation
#'
#' Estimates identified sets for treatment-effect levels and dose-response
#' slopes in continuous difference-in-differences designs under local
#' restrictions on counterfactual trends. A level restriction (bound M on
#' trend deviations) delivers identified sets for \eqn{ATT(d,t|d)} and the
#' overall summary \eqn{ATT^o}; a separate slope restriction (bound B on
#' selection slopes) delivers identified sets for the dose-response
#' derivative \eqn{\partial ATT(d,t|d)/\partial d}.
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
#'   (used for calibrating M and B). If NULL (default), all periods before the
#'   earliest post-period are used.
#' @param M Numeric scalar, numeric vector, or \code{"calibrate"}. Level
#'   sensitivity parameter bounding trend deviations,
#'   \eqn{|\mu_t(d) - \mu_t(0)| \leq M} (Assumption L). Governs identified
#'   sets for \eqn{ATT(d,t|d)} and \eqn{ATT^o}. If \code{"calibrate"},
#'   estimated from pre-treatment periods. If 0, standard parallel trends
#'   (point identification of levels). If a numeric vector, bounds are
#'   computed for each value (sensitivity analysis). Requires untreated
#'   units; ignored otherwise. Default: \code{"calibrate"}.
#' @param B Numeric scalar, numeric vector, or \code{"calibrate"}. Slope
#'   sensitivity parameter bounding the selection slope
#'   \eqn{|\mu'_t(d)| \leq B} on the interior of the treated dose support
#'   (Assumption S). Governs identified sets for the dose-response
#'   derivative only. If \code{"calibrate"}, estimated from pre-treatment
#'   periods. Default: \code{"calibrate"}.
#' @param method Character. Estimation method: \code{"gam"} (default, penalized
#'   splines via \code{mgcv}), \code{"contdid"} (B-splines via the
#'   \code{contdid} package; Callaway, Goodman-Bacon & Sant'Anna 2024),
#'   \code{"npiv"} (B-spline sieve regression via the \code{npiv} package;
#'   Chen, Christensen & Kankanala 2024), or \code{"kernel"} (local polynomial
#'   kernel regression via the \code{np} package; Hayfield & Racine 2008).
#' @param contdid_args Named list. Additional arguments passed to
#'   \code{contdid::cont_did()} when \code{method = "contdid"}. Common options:
#'   \code{num_knots} (default 1), \code{degree} (default 3),
#'   \code{biters} (bootstrap iterations, default 500),
#'   \code{control_group} (default \code{"notyettreated"}).
#'   Ignored when \code{method = "gam"} or \code{"npiv"}.
#' @param npiv_args Named list. Additional arguments passed to
#'   \code{npiv::npiv()} when \code{method = "npiv"}. Common options:
#'   \code{J.x.segments} (default 2), \code{J.x.degree} (default 3),
#'   \code{knots} (default \code{"uniform"}).
#'   Ignored when \code{method != "npiv"}.
#' @param kernel_args Named list. Additional arguments when
#'   \code{method = "kernel"}. \code{bw}: bandwidth — \code{"cv.ls"}
#'   (default, least-squares CV), \code{"cv.aic"}, or a positive numeric
#'   scalar. \code{regtype}: \code{"ll"} (local linear, default) or
#'   \code{"lc"} (local constant). \code{ckertype}: \code{"gaussian"}
#'   (default), \code{"epanechnikov"}, or \code{"uniform"}.
#'   \code{bw_inflate}: logical (default \code{TRUE}), inflate CV bandwidth
#'   by \eqn{n^{2/35}} for derivative estimation (Fan & Gijbels, 1996).
#'   Ignored when \code{method != "kernel"}.
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#'   Default: 50 points over the 5th-95th percentile of \emph{positive}
#'   doses — a trimmed interior of the treated support, where spline
#'   derivative estimates are reliable. Ignored when
#'   \code{method = "contdid"} (determined by contdid internally).
#'   Used by \code{"gam"}, \code{"npiv"}, and \code{"kernel"} methods.
#' @param k Integer. Spline basis dimension for GAM method. Default: 5.
#'   Ignored when \code{method != "gam"}.
#' @param spline_bs Character. Spline basis type (\code{"cr"} or \code{"tp"}).
#'   Default: \code{"cr"}. Ignored when \code{method != "gam"}.
#' @return An S3 object of class \code{"lpt"} containing:
#'   \describe{
#'     \item{datt}{Data frame with columns \code{period}, \code{horizon},
#'       \code{d}, \code{lambda_d}, \code{B}, \code{datt_lower},
#'       \code{datt_upper}. Identified set width is \code{2*(horizon+1)*B}.}
#'     \item{att}{Data frame with ATT level bounds (NULL if no untreated
#'       units). Columns \code{period}, \code{horizon}, \code{d},
#'       \code{Lambda_d}, \code{M}, \code{att_lower}, \code{att_upper};
#'       width is \code{2*(horizon+1)*M}, constant across doses.}
#'     \item{att_o}{Data frame with per-period overall ATT summary. Columns
#'       \code{period}, \code{horizon}, \code{att_o_bin}, \code{M},
#'       \code{att_o_lower}, \code{att_o_upper}; width is
#'       \code{2*(horizon+1)*M}. NULL if no untreated units.}
#'     \item{att_o_agg}{Data frame with the time-aggregated overall ATT
#'       across all requested post-periods; half-width is
#'       \code{mean(horizon+1)*M} (equal to \code{(T+2)/2*M} for contiguous
#'       horizons 0..T). NULL if single post-period or no untreated units.}
#'     \item{pre_att_o}{Data frame with pre-period binary DiD estimates
#'       (columns \code{period}, \code{att_o_bin}). Used by the event study
#'       plot. NULL if no untreated units or fewer than 2 pre-periods.}
#'     \item{M_hat}{The primary level bound used (NA if no untreated units).}
#'     \item{M_values}{All M values computed.}
#'     \item{B_hat}{The primary slope bound used.}
#'     \item{B_values}{All B values computed.}
#'     \item{M_source, B_source}{\code{"calibrated"} or \code{"user-supplied"}.}
#'     \item{calibration}{Output from \code{\link{calibrate_bounds}} if
#'       calibration was used.}
#'     \item{slopes}{Named list of \code{\link{estimate_dose_slope}} results,
#'       one per post-period (fit on treated units only).}
#'     \item{call}{The matched call.}
#'     \item{n}{Number of units in estimation sample.}
#'     \item{has_untreated}{Logical. Whether untreated units (D=0) exist.}
#'     \item{post_periods}{The post-period(s) estimated.}
#'     \item{ref_period}{The reference (last pre-) period used for differencing.}
#'     \item{contdid_fit}{Raw \code{pte_results} object from \code{contdid}
#'       (NULL when \code{method != "contdid"}).}
#'     \item{npiv_fits}{Named list of raw \code{npiv} objects per post-period
#'       (NULL when \code{method != "npiv"}).}
#'     \item{kernel_fits}{Named list of kernel fit objects per post-period,
#'       each with \code{fit} and \code{bw} elements
#'       (NULL when \code{method != "kernel"}).}
#'     \item{specifications}{List of all estimation settings.}
#'   }
#'
#' @references
#' Callaway B, Goodman-Bacon A, Sant'Anna PHC (2024).
#' \dQuote{Difference-in-differences with a continuous treatment.}
#' \emph{National Bureau of Economic Research}.
#'
#' Rambachan A, Roth J (2023).
#' \dQuote{A More Credible Approach to Parallel Trends.}
#' \emph{Review of Economic Studies}, \strong{90}(5), 2555--2591.
#'
#' Wood SN (2017).
#' \emph{Generalized Additive Models: An Introduction with R} (2nd ed.).
#' Chapman and Hall/CRC.
#'
#' Chen X, Christensen T, Kankanala S (2024).
#' \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric
#' Structural Functions and Elasticities.}
#' \emph{Review of Economic Studies}, \strong{91}(6), 3337--3369.
#'
#' Hayfield T, Racine JS (2008).
#' \dQuote{Nonparametric Econometrics: The np Package.}
#' \emph{Journal of Statistical Software}, \strong{27}(5).
#'
#' @details
#' Let \eqn{\mu_t(d) = E[\Delta Y_t(0) | D = d]} denote the per-period
#' untreated trend conditional on dose (the selection function), and let
#' \eqn{\Lambda_t(d) = E[Y_t - Y_{ref} | D=d] - E[Y_t - Y_{ref} | D=0]} be
#' the horizon-t DiD estimand with \eqn{\lambda_t(d)} its dose derivative.
#' Two logically independent assumptions bound the selection function:
#' \itemize{
#'   \item \strong{Level (Assumption L):} \eqn{|\mu_t(d) - \mu_t(0)| \leq M}
#'     for all \eqn{d} in the treated support. Then
#'     \deqn{IS_{ATT}(d, t; M) = [\Lambda_t(d) - (t+1)M,\
#'       \Lambda_t(d) + (t+1)M].}
#'   \item \strong{Slope (Assumption S):} \eqn{|\mu'_t(d)| \leq B} on the
#'     interior of the treated support. Then
#'     \deqn{IS_{\partial ATT}(d, t; B) = [\lambda_t(d) - (t+1)B,\
#'       \lambda_t(d) + (t+1)B].}
#' }
#' Identified sets widen linearly in the horizon because each post period
#' contributes another unit of counterfactual drift. The time-aggregated
#' summary \eqn{ATT^o} has identified set
#' \eqn{[\bar\Lambda^{agg} \mp \frac{T+2}{2} M]} for contiguous horizons
#' 0..T. Setting \eqn{M = 0} (resp. \eqn{B = 0}) recovers point
#' identification under standard parallel trends for the corresponding
#' target. Neither assumption implies the other: a slope bound does not
#' control the deviation between treated doses and the untreated group,
#' because the dose support has a hole at \eqn{(0, d_l)} which no
#' within-support restriction can bridge.
#'
#' Estimation follows the pipeline of the accompanying paper: conditional
#' means and derivatives are fit on \strong{treated units only} (no smooth
#' fit bridges the support gap), \eqn{\hat{f}_t(0)} is the sample mean of
#' the long difference among untreated units, and all evaluations are on a
#' grid in a trimmed interior of the treated dose support.
#'
#' @examples
#' data(sru)
#' fit <- lpt(sru, "commune", "year", "outcome", "dose",
#'            post_period = 0:5, pre_periods = -7:-1,
#'            M = "calibrate", B = "calibrate")
#' fit
#'
#' @export
lpt <- function(data, id_col, time_col, outcome_col, dose_col,
                post_period, pre_periods = NULL,
                M = "calibrate", B = "calibrate",
                method = c("gam", "contdid", "npiv", "kernel"),
                contdid_args = list(), npiv_args = list(),
                kernel_args = list(),
                eval_points = NULL,
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
  npiv_fits <- NULL
  kernel_fits <- NULL

  # ==================================================================
  #  Estimation: branch on method
  # ==================================================================

  if (method == "gam") {

    # --- Evaluation grid: trimmed interior of the treated dose support ---
    if (is.null(eval_points)) {
      d_pos <- dose_vec[dose_vec > 0]
      rng <- stats::quantile(d_pos, probs = c(0.05, 0.95))
      eval_points <- seq(rng[[1]], rng[[2]], length.out = 50)
    }

    # --- Per post-period: fit long difference on treated units only ---
    slopes <- list()
    Lambda_by_period <- list()
    att_o_bin_by_period <- list()

    for (pp in post_periods) {
      post_data <- data[data[[time_col]] == pp, ]

      merged <- merge(
        ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
        post_data[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_ref", "_post")
      )

      delta_y <- merged[[paste0(outcome_col, "_post")]] -
                 merged[[paste0(outcome_col, "_ref")]]
      wd <- merged[[dose_col]]
      treated_idx <- wd > 0
      untreated_idx <- wd == 0

      slope_result <- estimate_dose_slope(
        delta_y = delta_y[treated_idx],
        dose = wd[treated_idx],
        eval_points = eval_points,
        k = k,
        spline_bs = spline_bs
      )
      slopes[[as.character(pp)]] <- slope_result

      # Lambda_t(d) = f_t(d) - f_t(0), with f_t(0) the untreated sample mean
      if (has_untreated) {
        f_at_zero <- mean(delta_y[untreated_idx])
        Lambda_by_period[[as.character(pp)]] <-
          slope_result$conditional_mean - f_at_zero
        att_o_bin_by_period[[as.character(pp)]] <-
          mean(delta_y[treated_idx]) - f_at_zero
      }
    }

    # --- Calibrate / resolve M and B ---
    calibration_result <- NULL
    if (wants_calibration(M, B, has_untreated)) {
      if (length(pre_period_set) < 2) {
        stop("Calibration requires at least 2 pre-treatment periods. ",
             "Supply numeric M and B, or provide data with more pre-periods.")
      }
      calibration_result <- calibrate_bounds(
        data = data,
        id_col = id_col,
        time_col = time_col,
        outcome_col = outcome_col,
        dose_col = dose_col,
        pre_periods = pre_period_set,
        eval_points = eval_points,
        k = k,
        spline_bs = spline_bs
      )
    }
    rb <- resolve_bounds(M, B, calibration_result, has_untreated)
    announce_calibration(rb, length(pre_period_set) - 1)

    # --- Construct identified sets per period ---
    datt_all <- list()
    att_all <- list()
    att_o_all <- list()

    for (pp in post_periods) {
      pp_char <- as.character(pp)
      sr <- slopes[[pp_char]]
      lambda_d <- sr$lambda_d
      horizon_t <- as.integer(pp - min_post)
      mult <- horizon_t + 1L

      # IS_{dATT}(d, t; B) = [lambda_t(d) - (t+1)B, lambda_t(d) + (t+1)B]
      for (b in rb$B_values) {
        datt_all[[length(datt_all) + 1]] <- data.frame(
          period = pp,
          horizon = horizon_t,
          d = eval_points,
          lambda_d = lambda_d,
          B = b,
          datt_lower = lambda_d - mult * b,
          datt_upper = lambda_d + mult * b
        )
      }

      if (has_untreated && length(rb$M_values) > 0) {
        # IS_{ATT}(d, t; M) = [Lambda_t(d) - (t+1)M, Lambda_t(d) + (t+1)M]
        att_all[[length(att_all) + 1]] <- compute_att_bounds(
          eval_points = eval_points,
          Lambda_d = Lambda_by_period[[pp_char]],
          M_values = rb$M_values,
          period = pp, horizon = horizon_t
        )

        # IS_{ATT^o_t}(M) = [Lambda-bar_t - (t+1)M, Lambda-bar_t + (t+1)M]
        att_o_bin <- att_o_bin_by_period[[pp_char]]
        for (m in rb$M_values) {
          att_o_all[[length(att_o_all) + 1]] <- data.frame(
            period = pp,
            horizon = horizon_t,
            att_o_bin = att_o_bin,
            M = m,
            att_o_lower = att_o_bin - mult * m,
            att_o_upper = att_o_bin + mult * m
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

  } else if (method == "contdid") {

    # --- contdid: delegate to backend ---
    cd <- run_contdid(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      post_periods = post_periods, pre_period_set = pre_period_set,
      min_post = min_post, has_untreated = has_untreated,
      dose_vec = dose_vec, M = M, B = B, contdid_args = contdid_args
    )
    slopes             <- cd$slopes
    calibration_result <- cd$calibration
    rb                 <- cd$bounds
    datt               <- cd$datt
    att                <- cd$att
    att_o              <- cd$att_o
    contdid_fit        <- cd$contdid_fit

  } else if (method == "npiv") {

    # --- npiv: delegate to backend ---
    np <- run_npiv(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      post_periods = post_periods, pre_period_set = pre_period_set,
      min_post = min_post, ref_period = ref_period,
      has_untreated = has_untreated,
      dose_vec = dose_vec, M = M, B = B,
      eval_points = eval_points, npiv_args = npiv_args
    )
    slopes             <- np$slopes
    calibration_result <- np$calibration
    rb                 <- np$bounds
    datt               <- np$datt
    att                <- np$att
    att_o              <- np$att_o
    npiv_fits          <- np$npiv_fits

  } else {

    # --- kernel: delegate to backend ---
    kr <- run_kernel(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      post_periods = post_periods, pre_period_set = pre_period_set,
      min_post = min_post, ref_period = ref_period,
      has_untreated = has_untreated,
      dose_vec = dose_vec, M = M, B = B,
      eval_points = eval_points, kernel_args = kernel_args
    )
    slopes             <- kr$slopes
    calibration_result <- kr$calibration
    rb                 <- kr$bounds
    datt               <- kr$datt
    att                <- kr$att
    att_o              <- kr$att_o
    kernel_fits        <- kr$kernel_fits
  }

  # ==================================================================
  #  Common: ATT^o aggregation, pre-period event study, output
  # ==================================================================

  # --- Time-aggregated ATT^o (Eq. 4.19): half-width mean(t+1)*M ---
  att_o_agg <- NULL
  if (has_untreated && length(post_periods) > 1 && !is.null(att_o)) {
    horizons <- as.integer(post_periods - min_post)
    K <- length(post_periods)
    # mean(horizons + 1) equals (T+2)/2 for contiguous horizons 0..T
    avg_mult <- mean(horizons + 1L)

    for (m in rb$M_values) {
      att_o_m <- att_o[att_o$M == m, ]
      Lambda_agg <- mean(att_o_m$att_o_bin)
      half_width <- avg_mult * m

      att_o_agg <- rbind(att_o_agg, data.frame(
        Lambda_agg = Lambda_agg,
        n_periods = K,
        M = m,
        att_o_agg_lower = Lambda_agg - half_width,
        att_o_agg_upper = Lambda_agg + half_width
      ))
    }
    rownames(att_o_agg) <- NULL
  }

  if (!has_untreated) {
    message("No untreated units (dose = 0). ATT level bounds and the ",
            "level bound M are not computed.")
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
      M_hat = rb$M_hat,
      M_values = rb$M_values,
      B_hat = rb$B_hat,
      B_values = rb$B_values,
      M_source = rb$M_source,
      B_source = rb$B_source,
      calibration = calibration_result,
      slopes = slopes,
      contdid_fit = contdid_fit,
      npiv_fits = npiv_fits,
      kernel_fits = kernel_fits,
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
        npiv_args = if (method == "npiv") npiv_args else NULL,
        kernel_args = if (method == "kernel") kernel_args else NULL,
        post_periods = post_periods, ref_period = ref_period
      )
    ),
    class = "lpt"
  )
}


#' Does the M/B specification require pre-period calibration?
#'
#' @param M,B User-supplied specifications.
#' @param has_untreated Logical. Calibrating M requires untreated units.
#'
#' @keywords internal
wants_calibration <- function(M, B, has_untreated) {
  (is.character(B) && length(B) == 1 && B == "calibrate") ||
    (is.character(M) && length(M) == 1 && M == "calibrate" && has_untreated)
}


#' Message the user about calibrated bounds
#'
#' @param rb Output of \code{\link{resolve_bounds}}.
#' @param n_pairs Number of pre-period pairs used.
#'
#' @keywords internal
announce_calibration <- function(rb, n_pairs) {
  parts <- character(0)
  if (rb$M_source == "calibrated" && is.finite(rb$M_hat)) {
    parts <- c(parts, sprintf("M = %.4f", rb$M_hat))
  }
  if (rb$B_source == "calibrated") {
    parts <- c(parts, sprintf("B = %.4f", rb$B_hat))
  }
  if (length(parts) > 0) {
    message(sprintf("Calibrated %s from %d pre-period pair(s).",
                    paste(parts, collapse = ", "), n_pairs))
  }
  invisible(NULL)
}
