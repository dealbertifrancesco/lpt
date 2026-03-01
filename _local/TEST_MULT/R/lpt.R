#' Local Parallel Trends estimation
#'
#' Supports three multi-period LPT variants:
#' \itemize{
#'   \item \code{"C"} (LPT-C): smoothness on the cumulative trend; IS width constant.
#'   \item \code{"P"} (LPT-P): smoothness on period-specific increments; IS width grows linearly.
#'   \item \code{"S"} (LPT-S): smooth curvature path; IS width grows via Phi(t).
#' }
#'
#' @export
lpt <- function(data, id_col, time_col, outcome_col, dose_col,
                post_period, pre_periods = NULL,
                B = "calibrate", eval_points = NULL,
                k = 5, spline_bs = "cr", alpha = 0.05,
                lpt_type = c("C", "P", "S")) {

  lpt_type <- match.arg(lpt_type)

  # --- Input validation ---
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  required_cols <- c(id_col, time_col, outcome_col, dose_col)
  missing_cols  <- setdiff(required_cols, names(data))
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
    ref_period     <- max(pre_period_set)
  }

  if (length(pre_period_set) < 1) {
    stop("Need at least one pre-period before the earliest post-period.")
  }

  # --- Compute t_index for each post-period ---
  ref_rank    <- which(all_times == ref_period)
  t_index_map <- vapply(post_periods, function(pp) {
    as.integer(which(all_times == pp) - ref_rank)
  }, integer(1L))
  names(t_index_map) <- as.character(post_periods)

  # --- Reference-period data ---
  ref_data      <- data[data[[time_col]] == ref_period, ]
  dose_vec      <- ref_data[[dose_col]]
  has_untreated <- any(dose_vec == 0)

  # --- Estimate dose slope for each post-period ---
  slopes <- list()
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
    sr <- estimate_dose_slope(delta_y = delta_y, dose = wd,
                              eval_points = eval_points,
                              k = k, spline_bs = spline_bs)
    if (is.null(eval_points) && length(slopes) == 0) eval_points <- sr$eval_points
    slopes[[as.character(pp)]] <- sr
  }

  # --- Handle B ---
  calibration_result <- NULL
  b0_val <- NULL   # only populated for LPT-S
  C_hat  <- NULL   # only populated for LPT-S

  if (is.character(B) && B == "calibrate") {
    if (length(pre_period_set) < 2) {
      stop("B = 'calibrate' requires at least 2 pre-treatment periods.")
    }
    cal_type <- switch(lpt_type,
      "C" = "cumulative",
      "P" = "first_diff",
      "S" = "smooth"
    )
    calibration_result <- calibrate_B(
      data = data, id_col = id_col, time_col = time_col,
      outcome_col = outcome_col, dose_col = dose_col,
      pre_periods = pre_period_set,
      eval_points = slopes[[1]]$eval_points,
      k = k, spline_bs = spline_bs, type = cal_type
    )
    if (lpt_type == "S") {
      b0_val   <- calibration_result$B_hat  # B_hat stores b0 for smooth type
      C_hat    <- calibration_result$C_hat
      B_hat    <- b0_val
      B_values <- b0_val
      message(sprintf("Calibrated b0 = %.4f, C = %.4f (%s) from %d pre-period(s).",
                      b0_val, C_hat, cal_type, length(pre_period_set) - 1))
    } else {
      B_hat    <- calibration_result$B_hat
      B_values <- B_hat
      message(sprintf("Calibrated B = %.4f (%s) from %d pre-period(s).",
                      B_hat, cal_type, length(pre_period_set) - 1))
    }

  } else if (is.numeric(B)) {
    if (lpt_type == "S") {
      # Accept c(b0=..., C=...) or length-2 unnamed, or scalar (C defaults to 0)
      if (!is.null(names(B)) && all(c("b0", "C") %in% names(B))) {
        b0_val <- B[["b0"]]; C_hat <- B[["C"]]
      } else if (length(B) == 2) {
        b0_val <- B[1]; C_hat <- B[2]
      } else if (length(B) == 1) {
        b0_val <- B[1]; C_hat <- 0
      } else {
        stop("For lpt_type='S', B must be c(b0=..., C=...) or length 1-2 numeric.")
      }
      if (b0_val < 0 || C_hat < 0) stop("b0 and C must be non-negative.")
      B_hat    <- b0_val
      B_values <- b0_val
    } else {
      B_values <- as.numeric(B)
      if (any(!is.finite(B_values)) || any(B_values < 0)) {
        stop("B must be non-negative finite numeric values.")
      }
      B_hat <- max(B_values)
    }
  } else {
    stop("B must be numeric or 'calibrate'.")
  }

  # --- Construct identified sets per period ---
  z_alpha   <- stats::qnorm(1 - alpha / 2)
  datt_all  <- list(); att_all <- list(); att_o_all <- list()

  for (pp in post_periods) {
    pp_char   <- as.character(pp)
    sr        <- slopes[[pp_char]]
    ep        <- sr$eval_points
    lambda_d  <- sr$lambda_d
    se_lambda <- sr$se_lambda
    t_idx     <- t_index_map[[pp_char]]

    # Per-period multiplier
    # LPT-C: t_mult = 1        (IS width constant in t)
    # LPT-P: t_mult = t        (IS width linear in t)
    # LPT-S: t_mult = Phi(t)   (IS width ~quadratic in t)
    if (lpt_type == "C") {
      t_mult <- 1L
    } else if (lpt_type == "P") {
      t_mult <- t_idx
    } else {
      t_mult <- t_idx * b0_val + C_hat * t_idx * (t_idx + 1) / 2
    }

    # For LPT-S, the full bias is Phi(t) = t_mult, independent of B.
    # Pass b_eff=1 so that bias_margin = t_mult * 1 = Phi(t).
    b_loop  <- if (lpt_type == "S") 1 else B_values

    for (b in b_loop) {
      bias_margin <- t_mult * b
      b_store     <- if (lpt_type == "S") B_hat else b
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period     = pp,        t_index   = t_idx,
        d          = ep,        lambda_d  = lambda_d,
        se_lambda  = se_lambda, B         = b_store,
        datt_lower = lambda_d - bias_margin,
        datt_upper = lambda_d + bias_margin,
        ci_lower   = lambda_d - bias_margin - z_alpha * se_lambda,
        ci_upper   = lambda_d + bias_margin + z_alpha * se_lambda
      )
    }

    if (has_untreated) {
      b_for_att <- if (lpt_type == "S") 1 else B_values
      att_pp    <- compute_att_bounds(sr, b_for_att, dose_vec,
                                       period = pp, t_multiplier = t_mult)
      if (lpt_type == "S") att_pp$B <- B_hat
      att_pp$t_index <- t_idx
      att_all[[length(att_all) + 1]] <- att_pp
    }

    if (has_untreated) {
      post_data   <- data[data[[time_col]] == pp, ]
      merged_atto <- merge(
        ref_data[, c(id_col, outcome_col, dose_col), drop = FALSE],
        post_data[, c(id_col, outcome_col), drop = FALSE],
        by = id_col, suffixes = c("_ref", "_post")
      )
      dy       <- merged_atto[[paste0(outcome_col, "_post")]] -
                  merged_atto[[paste0(outcome_col, "_ref")]]
      d_merged <- merged_atto[[dose_col]]

      treated_idx   <- d_merged > 0
      untreated_idx <- d_merged == 0
      att_o_bin     <- mean(dy[treated_idx]) - mean(dy[untreated_idx])
      D_bar         <- mean(d_merged[treated_idx])

      b_for_atto <- if (lpt_type == "S") 1 else B_values
      for (b in b_for_atto) {
        bias_margin_o <- t_mult * b * D_bar
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period      = pp,        t_index    = t_idx,
          att_o_bin   = att_o_bin, D_bar      = D_bar,
          B           = if (lpt_type == "S") B_hat else b,
          att_o_lower = att_o_bin - bias_margin_o,
          att_o_upper = att_o_bin + bias_margin_o
        )
      }
    }
  }

  datt  <- do.call(rbind, datt_all); rownames(datt) <- NULL
  att   <- if (length(att_all)   > 0) { r <- do.call(rbind, att_all);   rownames(r) <- NULL; r } else NULL
  att_o <- if (length(att_o_all) > 0) { r <- do.call(rbind, att_o_all); rownames(r) <- NULL; r } else NULL

  if (!has_untreated) {
    message("No untreated units (dose = 0). ATT level bounds not computed.")
  }

  structure(
    list(
      datt          = datt,
      att           = att,
      att_o         = att_o,
      B_hat         = B_hat,
      B_values      = B_values,
      b0            = b0_val,
      C_hat         = C_hat,
      calibration   = calibration_result,
      slopes        = slopes,
      call          = the_call,
      n             = nrow(ref_data),
      has_untreated = has_untreated,
      post_periods  = post_periods,
      ref_period    = ref_period,
      lpt_type      = lpt_type,
      t_index_map   = t_index_map,
      specifications = list(
        k = k, spline_bs = spline_bs, alpha = alpha,
        post_periods = post_periods, ref_period = ref_period,
        lpt_type = lpt_type
      )
    ),
    class = "lpt"
  )
}
