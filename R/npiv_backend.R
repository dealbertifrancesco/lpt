#' Run npiv estimation backend
#'
#' Internal function that fits nonparametric B-spline sieve regressions
#' on pairwise first differences per period, following the same approach
#' as the GAM backend but using \code{npiv::npiv()} as the estimator.
#'
#' @param data Data frame in long format.
#' @param id_col,time_col,outcome_col,dose_col Column name strings.
#' @param post_periods Sorted numeric vector of post-period identifiers.
#' @param pre_period_set Sorted numeric vector of pre-period identifiers.
#' @param min_post Minimum post-period value.
#' @param ref_period Scalar. Reference (last pre-) period for differencing.
#' @param has_untreated Logical. Whether untreated units exist.
#' @param dose_vec Numeric vector of doses from the reference period.
#' @param B Sensitivity parameter specification ("calibrate" or numeric).
#' @param eval_points Numeric vector or NULL. Dose grid for evaluation.
#' @param npiv_args List of additional arguments for \code{npiv()}.
#'
#' @return A list with components: \code{datt}, \code{att}, \code{att_o},
#'   \code{slopes}, \code{calibration}, \code{B_hat}, \code{B_values},
#'   \code{npiv_fits}.
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
                     has_untreated, dose_vec, B, eval_points, npiv_args) {

  if (!requireNamespace("npiv", quietly = TRUE)) {
    stop("Package 'npiv' required for method = 'npiv'. ",
         "Install with: install.packages('npiv')")
  }

  ref_data <- data[data[[time_col]] == ref_period, ]

  # --- Default eval_points ---
  if (is.null(eval_points)) {
    dose_range <- stats::quantile(dose_vec, probs = c(0.05, 0.95))
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

  # --- B calibration from pre-period pairs ---
  calibration_result <- NULL
  B_hat <- NULL
  B_values <- NULL

  if (is.character(B) && B == "calibrate") {
    if (length(pre_period_set) < 2) {
      stop("B = 'calibrate' requires at least 2 pre-treatment periods. ",
           "Supply B as a numeric value, or provide data with more pre-periods.")
    }

    pre_slopes_list <- list()
    sup_vals <- numeric(0)
    pair_labels <- character(0)

    for (j in seq_len(length(pre_period_set) - 1)) {
      t0 <- pre_period_set[j]
      t1 <- pre_period_set[j + 1]
      pair_label <- paste0(t0, "-", t1)
      pair_labels <- c(pair_labels, pair_label)

      dat_t0 <- data[data[[time_col]] == t0, ]
      dat_t1 <- data[data[[time_col]] == t1, ]

      merged <- merge(
        dat_t0[, c(id_col, outcome_col, dose_col), drop = FALSE],
        dat_t1[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_0", "_1")
      )

      dy <- merged[[paste0(outcome_col, "_1")]] - merged[[paste0(outcome_col, "_0")]]
      d_pre <- merged[[dose_col]]

      fit_pre <- fit_one_npiv(dy, d_pre, eval_points)

      pre_slopes_list[[j]] <- data.frame(
        period_pair = pair_label,
        d           = eval_points,
        mu_prime_d  = as.numeric(fit_pre$deriv)
      )
      sup_vals <- c(sup_vals, max(abs(fit_pre$deriv)))
    }

    names(sup_vals) <- pair_labels
    pre_slopes <- do.call(rbind, pre_slopes_list)
    rownames(pre_slopes) <- NULL

    B_hat <- max(sup_vals)
    B_values <- B_hat

    calibration_result <- list(
      B_hat         = B_hat,
      pre_slopes    = pre_slopes,
      sup_by_period = sup_vals
    )

    message(sprintf("Calibrated B = %.4f from %d pre-period pair(s) [npiv].",
                    B_hat, length(sup_vals)))

  } else if (is.numeric(B)) {
    B_values <- as.numeric(B)
    if (any(!is.finite(B_values)) || any(B_values < 0)) {
      stop("B must be non-negative finite numeric values.")
    }
    B_hat <- max(B_values)
  } else {
    stop("B must be numeric or 'calibrate'.")
  }

  # --- Post-period estimation ---
  D_bar <- mean(dose_vec[dose_vec > 0])
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

    # Fit npiv on all observations; prepend 0 to eval grid for Lambda_d
    ep_with_zero <- c(0, eval_points)
    fit_pp <- fit_one_npiv(dy, wd, ep_with_zero)

    h_all     <- as.numeric(fit_pp$h)
    deriv_all <- as.numeric(fit_pp$deriv)
    h_at_zero <- h_all[1]
    h_at_eval <- h_all[-1]
    lambda_d  <- deriv_all[-1]
    Lambda_d  <- h_at_eval - h_at_zero

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
    for (b in B_values) {
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

    # IS_{ATT}(d, t; B) = [Lambda_t(d) - (t+1)Bd, Lambda_t(d) + (t+1)Bd]
    if (has_untreated) {
      for (b in B_values) {
        att_all[[length(att_all) + 1]] <- data.frame(
          period    = pp,
          horizon   = horizon_t,
          d         = eval_points,
          Lambda_d  = Lambda_d,
          B         = b,
          att_lower = Lambda_d - mult * b * eval_points,
          att_upper = Lambda_d + mult * b * eval_points
        )
      }
    }

    # IS_{ATT^o_t}(B)
    if (has_untreated) {
      treated_idx   <- wd > 0
      untreated_idx <- wd == 0
      att_o_bin <- mean(dy[treated_idx]) - mean(dy[untreated_idx])

      for (b in B_values) {
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period      = pp,
          horizon     = horizon_t,
          att_o_bin   = att_o_bin,
          D_bar       = D_bar,
          B           = b,
          att_o_lower = att_o_bin - mult * b * D_bar,
          att_o_upper = att_o_bin + mult * b * D_bar
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
    B_hat       = B_hat,
    B_values    = B_values,
    npiv_fits   = npiv_fits
  )
}
