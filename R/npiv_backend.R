#' Run npiv estimation backend
#'
#' Internal function that fits nonparametric B-spline sieve regressions
#' on long differences per period, following the same approach as the GAM
#' backend but using \code{npiv::npiv()} as the estimator. Fits are on
#' treated units only; the untreated reference \eqn{\hat{f}_t(0)} is the
#' sample mean of the long difference among untreated units.
#'
#' @param data Data frame in long format.
#' @param id_col,time_col,outcome_col,dose_col Column name strings.
#' @param post_periods Sorted numeric vector of post-period identifiers.
#' @param pre_period_set Sorted numeric vector of pre-period identifiers.
#' @param min_post Minimum post-period value.
#' @param ref_period Scalar. Reference (last pre-) period for differencing.
#' @param has_untreated Logical. Whether untreated units exist.
#' @param dose_vec Numeric vector of doses from the reference period.
#' @param M Level bound specification ("calibrate" or numeric).
#' @param B Slope bound specification ("calibrate" or numeric).
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#' @param npiv_args List of additional arguments for \code{npiv()}.
#'
#' @return A list with components: \code{datt}, \code{att}, \code{att_o},
#'   \code{slopes}, \code{calibration}, \code{bounds}, \code{npiv_fits}.
#'
#' @references
#' Chen X, Christensen T, Kankanala S (2024).
#' \dQuote{Adaptive Estimation and Uniform Confidence Bands for Nonparametric
#' Structural Functions and Elasticities.}
#' \emph{Review of Economic Studies}, \strong{91}(6), 3337--3369.
#'
#' @keywords internal
run_npiv <- function(data, id_col, time_col, outcome_col, dose_col,
                     post_periods, pre_period_set, min_post, ref_period,
                     has_untreated, dose_vec, M, B, eval_points, npiv_args) {

  if (!requireNamespace("npiv", quietly = TRUE)) {
    stop("Package 'npiv' required for method = 'npiv'. ",
         "Install with: install.packages('npiv')")
  }

  ref_data <- data[data[[time_col]] == ref_period, ]

  # --- Evaluation grid: trimmed interior of the treated dose support ---
  if (is.null(eval_points)) {
    d_pos <- dose_vec[dose_vec > 0]
    dose_range <- stats::quantile(d_pos, probs = c(0.05, 0.95))
    eval_points <- seq(dose_range[[1]], dose_range[[2]], length.out = 50)
  }

  # --- npiv call defaults ---
  npiv_defaults <- list(
    J.x.segments = 2,
    J.x.degree   = 3,
    knots        = "uniform",
    progress     = FALSE,
    ucb.h        = FALSE,
    ucb.deriv    = FALSE
  )

  # --- Helper: fit npiv on one (delta_y, dose) pair ---
  fit_one_npiv <- function(delta_y, dose, ep) {
    args <- npiv_defaults
    # User overrides (except Y, X, W, X.eval which we control)
    for (nm in names(npiv_args)) {
      args[[nm]] <- npiv_args[[nm]]
    }
    # K.w.segments default: 2 * J.x.segments (following K.w.smooth = 2)
    if (is.null(args[["K.w.segments"]]) && !is.null(args[["J.x.segments"]])) {
      args[["K.w.segments"]] <- 2L * as.integer(args[["J.x.segments"]])
    }
    # Core args (always set internally)
    args[["Y"]]      <- delta_y
    args[["X"]]      <- dose
    args[["W"]]      <- dose          # W = X for regression (not IV)
    args[["X.eval"]] <- matrix(ep, ncol = 1)
    do.call(npiv::npiv, args)
  }

  # --- Calibrate / resolve M and B ---
  calibration_result <- NULL
  if (wants_calibration(M, B, has_untreated)) {
    if (length(pre_period_set) < 2) {
      stop("Calibration requires at least 2 pre-treatment periods. ",
           "Supply numeric M and B, or provide data with more pre-periods.")
    }
    calibration_result <- calibrate_bounds_engine(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      pre_periods = pre_period_set, eval_points = eval_points,
      fit_fun = function(dy, dd, ep) {
        f <- fit_one_npiv(dy, dd, ep)
        list(mean = as.numeric(f$h), deriv = as.numeric(f$deriv))
      }
    )
  }
  rb <- resolve_bounds(M, B, calibration_result, has_untreated)
  announce_calibration(rb, length(pre_period_set) - 1)

  # --- Post-period estimation ---
  slopes    <- list()
  npiv_fits <- list()
  datt_all  <- list()
  att_all   <- list()
  att_o_all <- list()

  for (pp in post_periods) {
    pp_char <- as.character(pp)
    post_data <- data[data[[time_col]] == pp, ]

    merged <- merge(
      ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
      post_data[, c(id_col, outcome_col), drop = FALSE],
      by = id_col, suffixes = c("_ref", "_post")
    )

    dy <- merged[[paste0(outcome_col, "_post")]] - merged[[paste0(outcome_col, "_ref")]]
    wd <- merged[[dose_col]]
    treated_idx <- wd > 0
    untreated_idx <- wd == 0

    # Fit on treated units only; untreated enter via their sample mean
    fit_pp <- fit_one_npiv(dy[treated_idx], wd[treated_idx], eval_points)

    h_at_eval <- as.numeric(fit_pp$h)
    lambda_d  <- as.numeric(fit_pp$deriv)

    npiv_fits[[pp_char]] <- fit_pp

    # Build dose_slope-compatible object for plot/summary
    slopes[[pp_char]] <- structure(
      list(
        eval_points      = eval_points,
        conditional_mean = h_at_eval,
        lambda_d         = lambda_d,
        gam_fit          = NULL
      ),
      class = "dose_slope"
    )

    horizon_t <- as.integer(pp - min_post)
    mult <- horizon_t + 1L

    # IS_{dATT}(d, t; B) = [lambda_t(d) - (t+1)B, lambda_t(d) + (t+1)B]
    for (b in rb$B_values) {
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period     = pp,
        horizon    = horizon_t,
        d          = eval_points,
        lambda_d   = lambda_d,
        B          = b,
        datt_lower = lambda_d - mult * b,
        datt_upper = lambda_d + mult * b
      )
    }

    if (has_untreated && length(rb$M_values) > 0) {
      f_at_zero <- mean(dy[untreated_idx])
      Lambda_d <- h_at_eval - f_at_zero

      # IS_{ATT}(d, t; M) = [Lambda_t(d) - (t+1)M, Lambda_t(d) + (t+1)M]
      att_all[[length(att_all) + 1]] <- compute_att_bounds(
        eval_points = eval_points, Lambda_d = Lambda_d,
        M_values = rb$M_values, period = pp, horizon = horizon_t
      )

      # IS_{ATT^o_t}(M)
      att_o_bin <- mean(dy[treated_idx]) - f_at_zero
      for (m in rb$M_values) {
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period      = pp,
          horizon     = horizon_t,
          att_o_bin   = att_o_bin,
          M           = m,
          att_o_lower = att_o_bin - mult * m,
          att_o_upper = att_o_bin + mult * m
        )
      }
    }
  }

  datt  <- do.call(rbind, datt_all)
  att   <- if (length(att_all)   > 0) do.call(rbind, att_all)   else NULL
  att_o <- if (length(att_o_all) > 0) do.call(rbind, att_o_all) else NULL
  rownames(datt) <- NULL
  if (!is.null(att))   rownames(att)   <- NULL
  if (!is.null(att_o)) rownames(att_o) <- NULL

  list(
    datt        = datt,
    att         = att,
    att_o       = att_o,
    slopes      = slopes,
    calibration = calibration_result,
    bounds      = rb,
    npiv_fits   = npiv_fits
  )
}
