#' Run kernel estimation backend
#'
#' Internal function that fits local polynomial kernel regressions
#' on long differences per period, following the same approach as the GAM
#' backend but using \code{np::npreg()} as the estimator. Fits are on
#' treated units only; the untreated reference \eqn{\hat{f}_t(0)} is the
#' sample mean of the long difference among untreated units.
#'
#' Local linear regression (\code{regtype = "ll"}) is the default and
#' automatically corrects for boundary bias (Fan & Gijbels, 1996).
#' Bandwidth can be selected via least-squares cross-validation or
#' supplied manually.
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
#' @param kernel_args List of additional arguments. Recognized entries:
#'   \describe{
#'     \item{bw}{Bandwidth: \code{"cv.ls"} (default, least-squares CV),
#'       \code{"cv.aic"} (AIC-based CV), or a positive numeric scalar.
#'       When using CV, the selected bandwidth is inflated for derivative
#'       estimation (see \code{bw_inflate}).}
#'     \item{regtype}{\code{"ll"} (local linear, default) or
#'       \code{"lc"} (local constant / Nadaraya-Watson).}
#'     \item{ckertype}{Kernel function: \code{"gaussian"} (default),
#'       \code{"epanechnikov"}, or \code{"uniform"}.}
#'     \item{bw_inflate}{Logical. If \code{TRUE} (default), inflate
#'       CV-selected bandwidths by \eqn{n^{2/35}} for derivative estimation
#'       (Fan & Gijbels, 1996, Theorem 3.1). CV optimizes for level
#'       estimation; derivatives need a larger bandwidth. Ignored when
#'       \code{bw} is numeric. Set to \code{FALSE} to use the raw CV
#'       bandwidth.}
#'   }
#'
#' @return A list with components: \code{datt}, \code{att}, \code{att_o},
#'   \code{slopes}, \code{calibration}, \code{bounds}, \code{kernel_fits}.
#'
#' @references
#' Hayfield T, Racine JS (2008).
#' \dQuote{Nonparametric Econometrics: The np Package.}
#' \emph{Journal of Statistical Software}, \strong{27}(5).
#'
#' Fan J, Gijbels I (1996).
#' \emph{Local Polynomial Modelling and Its Applications}.
#' Chapman and Hall/CRC.
#'
#' @keywords internal
run_kernel <- function(data, id_col, time_col, outcome_col, dose_col,
                       post_periods, pre_period_set, min_post, ref_period,
                       has_untreated, dose_vec, M, B, eval_points, kernel_args) {

  if (!requireNamespace("np", quietly = TRUE)) {
    stop("Package 'np' required for method = 'kernel'. ",
         "Install with: install.packages('np')")
  }

  ref_data <- data[data[[time_col]] == ref_period, ]

  # --- Evaluation grid: trimmed interior of the treated dose support ---
  if (is.null(eval_points)) {
    d_pos <- dose_vec[dose_vec > 0]
    dose_range <- stats::quantile(d_pos, probs = c(0.05, 0.95))
    eval_points <- seq(dose_range[[1]], dose_range[[2]], length.out = 50)
  }

  # --- Parse kernel_args ---
  bw_spec  <- if (!is.null(kernel_args[["bw"]]))       kernel_args[["bw"]]       else "cv.ls"
  regtype  <- if (!is.null(kernel_args[["regtype"]]))   kernel_args[["regtype"]]   else "ll"
  ckertype <- if (!is.null(kernel_args[["ckertype"]])) kernel_args[["ckertype"]]  else "gaussian"
  # Derivative bandwidth inflation: CV optimizes for h(d), not h'(d).
  # The optimal bandwidth for first-derivative estimation is larger by a
  # factor of n^{2/35} (Fan & Gijbels 1996, Theorem 3.1). Since lpt always
  # needs derivatives (lambda_d, calibration slopes), we inflate by default.
  bw_inflate <- if (!is.null(kernel_args[["bw_inflate"]])) kernel_args[["bw_inflate"]] else TRUE

  if (!regtype %in% c("ll", "lc")) {
    stop("kernel_args$regtype must be 'll' (local linear) or 'lc' (local constant).")
  }
  if (!ckertype %in% c("gaussian", "epanechnikov", "uniform")) {
    stop("kernel_args$ckertype must be 'gaussian', 'epanechnikov', or 'uniform'.")
  }

  # --- Helper: fit kernel regression on one (delta_y, dose) pair ---
  fit_one_kernel <- function(delta_y, dose, ep) {
    if (is.numeric(bw_spec)) {
      # Manual bandwidth — skip CV, no inflation
      bw_use <- bw_spec
    } else {
      # Data-driven bandwidth selection
      bw_method <- match.arg(bw_spec, c("cv.ls", "cv.aic"))
      bw_obj <- np::npregbw(
        ydat = delta_y, xdat = dose,
        regtype = regtype, ckertype = ckertype,
        bwmethod = bw_method
      )
      bw_use <- bw_obj$bw
      # Inflate for derivative estimation (Fan & Gijbels 1996)
      if (isTRUE(bw_inflate)) {
        n <- length(delta_y)
        bw_use <- bw_use * n^(2/35)
      }
    }
    # Fit with final bandwidth
    bw_final <- np::npregbw(
      ydat = delta_y, xdat = dose,
      regtype = regtype, ckertype = ckertype,
      bws = bw_use, bandwidth.compute = FALSE
    )
    fit <- np::npreg(bw_final, exdat = ep, gradients = TRUE)
    list(
      h     = as.numeric(stats::fitted(fit)),
      deriv = as.numeric(np::gradients(fit)),
      bw    = bw_use,
      fit   = fit
    )
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
        f <- fit_one_kernel(dy, dd, ep)
        list(mean = f$h, deriv = f$deriv)
      }
    )
  }
  rb <- resolve_bounds(M, B, calibration_result, has_untreated)
  announce_calibration(rb, length(pre_period_set) - 1)

  # --- Post-period estimation ---
  slopes      <- list()
  kernel_fits <- list()
  datt_all    <- list()
  att_all     <- list()
  att_o_all   <- list()

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
    kfit <- fit_one_kernel(dy[treated_idx], wd[treated_idx], eval_points)

    h_at_eval <- kfit$h
    lambda_d  <- kfit$deriv

    kernel_fits[[pp_char]] <- list(fit = kfit$fit, bw = kfit$bw)

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
    kernel_fits = kernel_fits
  )
}
