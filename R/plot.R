#' Plot LPT results
#'
#' Creates visualizations for identified sets, pre-trend diagnostics,
#' and sensitivity analysis.
#'
#' @param x An object of class \code{"lpt"}.
#' @param type Character. Plot type: \code{"datt"} (default), \code{"att"},
#'   \code{"pretrends"}, \code{"sensitivity"}, \code{"eventstudy"}.
#' @param period Scalar or NULL. Which post-period to plot. If NULL and
#'   multiple post-periods exist, plots all periods faceted. Default: NULL.
#' @param d0 Numeric. Required dose value when \code{estimand} is \code{"datt"}
#'   or \code{"att"}. The identified set for these estimands depends on the
#'   chosen dose level; no default is provided. Ignored for \code{"att_o"}.
#' @param B_grid Numeric vector or NULL. Grid of B values for the
#'   \code{"datt"} sensitivity plot.
#'   Default: \code{seq(0, 3 * B_hat, length.out = 50)}.
#' @param M_grid Numeric vector or NULL. Grid of M values for the
#'   \code{"att"} and \code{"att_o"} sensitivity plots.
#'   Default: \code{seq(0, 3 * M_hat, length.out = 50)}.
#' @param estimand Character. Estimand for the sensitivity plot:
#'   \code{"att_o"} (default, overall ATT \eqn{ATT^o}; no \code{d0} needed;
#'   traces M), \code{"datt"} (dose-response slope
#'   \eqn{\partial ATT/\partial d} at \code{d0}; traces B), or \code{"att"}
#'   (ATT level \eqn{ATT(d|d)} at \code{d0}; traces M).
#'   Only used when \code{type = "sensitivity"}.
#' @param true_curve Function or NULL. If provided, a function \code{f(d)} that
#'   returns the true value at dose \code{d}, overlaid as a dashed black line.
#'   Applies to \code{"datt"} and \code{"att"} plot types only.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' data(sru)
#' fit <- lpt(sru, "commune", "year", "outcome", "dose",
#'            post_period = 5, pre_periods = -7:-1, M = 0, B = 0)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(fit, type = "datt")
#' }
#'
#' @method plot lpt
#' @export
plot.lpt <- function(x, type = c("datt", "att", "pretrends", "sensitivity",
                                  "eventstudy"),
                      period = NULL, d0 = NULL, B_grid = NULL, M_grid = NULL,
                      estimand = c("att_o", "datt", "att"),
                      true_curve = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.lpt(). ",
         "Install with install.packages('ggplot2').")
  }

  type <- match.arg(type)
  estimand <- match.arg(estimand)

  if (!is.null(true_curve) && !is.function(true_curve)) {
    stop("true_curve must be a function(d) or NULL.")
  }

  # Okabe-Ito palette (colorblind-friendly)
  col_band <- "#0072B2"     # blue
  col_line <- "#D55E00"     # vermilion
  col_pretrend <- "#009E73" # green
  col_marker <- "#CC79A7"   # pink

  p <- switch(type,
    "datt" = plot_datt(x, period, col_band, col_line, true_curve),
    "att" = plot_att(x, period, col_band, col_line, true_curve),
    "pretrends" = plot_pretrends(x, col_pretrend, col_marker),
    "sensitivity" = plot_sensitivity(x, d0, B_grid, M_grid, col_band,
                                      col_line, col_marker,
                                      estimand = estimand, period = period),
    "eventstudy" = plot_eventstudy(x, col_band, col_pretrend)
  )

  print(p)
  invisible(p)
}


# --- Internal plot helpers ---

plot_datt <- function(x, period, col_band, col_line, true_curve = NULL) {
  df <- x$datt
  if (!is.null(period)) {
    df <- df[df$period %in% period, ]
  }

  multi_B <- length(unique(df$B)) > 1
  multi_period <- length(unique(df$period)) > 1

  if (multi_B) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$datt_lower, ymax = .data$datt_upper,
                     group = factor(.data$B), alpha = factor(.data$B)),
        fill = col_band
      ) +
      ggplot2::scale_alpha_discrete(
        name = "B",
        range = c(0.15, 0.4)
      )
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$datt_lower, ymax = .data$datt_upper),
        fill = col_band, alpha = 0.35
      )
  }

  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$lambda_d),
                        color = col_line, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      x = "Dose (d)",
      y = expression(partialdiff * ATT(d) / partialdiff * d),
      title = "Identified Set for Dose-Response Slope (slope bound B)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (multi_period) {
    p <- p + ggplot2::facet_wrap(~ period, labeller = ggplot2::label_both)
  }

  if (!is.null(true_curve)) {
    truth_df <- data.frame(d = unique(df$d))
    truth_df$true_val <- true_curve(truth_df$d)
    p <- p +
      ggplot2::geom_line(
        data = truth_df,
        ggplot2::aes(x = .data$d, y = .data$true_val),
        linetype = "dashed", color = "black", linewidth = 0.7
      )
  }

  p
}


plot_att <- function(x, period, col_band, col_line, true_curve = NULL) {
  if (is.null(x$att)) {
    stop("ATT level bounds not available (no untreated units in data).")
  }

  df <- x$att
  if (!is.null(period)) {
    df <- df[df$period %in% period, ]
  }

  multi_M <- length(unique(df$M)) > 1
  multi_period <- length(unique(df$period)) > 1

  if (multi_M) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$att_lower, ymax = .data$att_upper,
                     group = factor(.data$M), alpha = factor(.data$M)),
        fill = col_band
      ) +
      ggplot2::scale_alpha_discrete(
        name = "M",
        range = c(0.15, 0.4)
      )
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$att_lower, ymax = .data$att_upper),
        fill = col_band, alpha = 0.3
      )
  }

  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$Lambda_d),
                        color = col_line, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      x = "Dose (d)",
      y = expression(ATT(d)),
      title = "Identified Set for ATT(d|d) (level bound M)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (multi_period) {
    p <- p + ggplot2::facet_wrap(~ period, labeller = ggplot2::label_both)
  }

  if (!is.null(true_curve)) {
    truth_df <- data.frame(d = unique(df$d))
    truth_df$true_val <- true_curve(truth_df$d)
    p <- p +
      ggplot2::geom_line(
        data = truth_df,
        ggplot2::aes(x = .data$d, y = .data$true_val),
        linetype = "dashed", color = "black", linewidth = 0.7
      )
  }

  p
}


plot_pretrends <- function(x, col_pretrend, col_marker) {
  cal <- x$calibration
  if (is.null(cal)) {
    stop("Pre-trends plot requires calibration ",
         "(M = 'calibrate' or B = 'calibrate').")
  }

  panel_dev <- "Trend deviations (calibrate M)"
  panel_slope <- "Selection slopes (calibrate B)"

  dfs <- list()
  hlines <- list()

  if (!is.null(cal$pre_deviations) && is.finite(cal$M_hat)) {
    dfs[[length(dfs) + 1]] <- data.frame(
      panel = panel_dev,
      period_pair = cal$pre_deviations$period_pair,
      d = cal$pre_deviations$d,
      value = cal$pre_deviations$deviation
    )
    hlines[[length(hlines) + 1]] <- data.frame(
      panel = panel_dev,
      bound = c(-cal$M_hat, cal$M_hat)
    )
  }

  dfs[[length(dfs) + 1]] <- data.frame(
    panel = panel_slope,
    period_pair = cal$pre_slopes$period_pair,
    d = cal$pre_slopes$d,
    value = cal$pre_slopes$mu_prime_d
  )
  hlines[[length(hlines) + 1]] <- data.frame(
    panel = panel_slope,
    bound = c(-cal$B_hat, cal$B_hat)
  )

  df <- do.call(rbind, dfs)
  hl <- do.call(rbind, hlines)
  panel_levels <- intersect(c(panel_dev, panel_slope), unique(df$panel))
  df$panel <- factor(df$panel, levels = panel_levels)
  hl$panel <- factor(hl$panel, levels = panel_levels)

  subtitle <- if (is.finite(cal$M_hat)) {
    sprintf("M-hat = %.3f  |  B-hat = %.3f", cal$M_hat, cal$B_hat)
  } else {
    sprintf("B-hat = %.3f", cal$B_hat)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$d, y = .data$value,
                                    color = .data$period_pair)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(
      data = hl,
      ggplot2::aes(yintercept = .data$bound),
      linetype = "dashed", color = col_marker, linewidth = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "grey40",
                         linewidth = 0.3) +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::labs(
      x = "Dose (d)",
      y = NULL,
      title = "Pre-Period Selection Diagnostics",
      subtitle = subtitle,
      color = "Period pair"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


plot_sensitivity <- function(x, d0, B_grid, M_grid, col_band, col_line,
                              col_marker, estimand = "att_o", period = NULL) {

  # Select period and determine horizon
  pp <- if (is.null(period)) x$post_periods[1] else period
  pp_char <- as.character(pp)
  sr <- x$slopes[[pp_char]]
  if (is.null(sr)) {
    stop("No results for period ", pp, ".")
  }
  ep <- sr$eval_points

  horizon_t <- if ("horizon" %in% names(x$datt)) {
    h_rows <- x$datt[x$datt$period == pp, ]
    if (nrow(h_rows) > 0) h_rows$horizon[1] else 0L
  } else {
    0L
  }
  mult <- horizon_t + 1L

  if (estimand %in% c("datt", "att") && is.null(d0)) {
    stop("'d0' must be specified when estimand = '", estimand, "'. ",
         "Supply a dose value, e.g. d0 = ", round(stats::median(ep), 3), ".")
  }

  # The slope estimand traces B; the level estimands trace M
  if (estimand == "datt") {
    param_name <- "B"
    param_hat <- x$B_hat
    grid <- B_grid
  } else {
    param_name <- "M"
    param_hat <- x$M_hat
    grid <- M_grid
  }

  if (is.null(grid)) {
    if (is.finite(param_hat) && param_hat > 0) {
      grid <- seq(0, 3, length.out = 50) * param_hat
    } else {
      grid <- seq(0, 1, length.out = 50)
    }
  }
  use_ratio <- is.finite(param_hat) && param_hat > 0
  x_var <- if (use_ratio) "ratio" else "param"

  if (estimand == "datt") {
    # IS_dATT(d0, t; B) = [lambda_t(d0) - (t+1)B, lambda_t(d0) + (t+1)B]
    idx        <- which.min(abs(ep - d0))
    center_val <- sr$lambda_d[idx]
    d0_actual  <- ep[idx]
    y_label    <- expression(partialdiff * ATT(d) / partialdiff * d ~ "at" ~ d[0])
    plot_title <- sprintf("Sensitivity in B: dATT at d = %.2f (horizon %d)",
                          d0_actual, horizon_t)

  } else if (estimand == "att") {
    # IS_ATT(d0, t; M) = [Lambda_t(d0) - (t+1)M, Lambda_t(d0) + (t+1)M]
    if (is.null(x$att)) {
      stop("ATT level bounds not available (no untreated units in data). ",
           "Use estimand = 'datt'.")
    }
    att_pp <- x$att[x$att$period == pp, ]
    if (nrow(att_pp) == 0) {
      stop("No ATT data for period ", pp, ".")
    }
    # Lambda_d is constant across M; use any M value
    att_sub    <- att_pp[att_pp$M == att_pp$M[1], ]
    idx        <- which.min(abs(att_sub$d - d0))
    center_val <- att_sub$Lambda_d[idx]
    d0_actual  <- att_sub$d[idx]
    y_label    <- expression(ATT(d[0]))
    plot_title <- sprintf("Sensitivity in M: ATT(d|d) at d = %.2f (horizon %d)",
                          d0_actual, horizon_t)

  } else {
    # estimand == "att_o"
    # IS_{ATT^o_t}(M) = [Lambda-bar_t - (t+1)M, Lambda-bar_t + (t+1)M]
    if (is.null(x$att_o)) {
      stop("Overall ATT bounds not available (no untreated units in data). ",
           "Use estimand = 'datt'.")
    }
    att_o_pp <- x$att_o[x$att_o$period == pp, ]
    if (nrow(att_o_pp) == 0) {
      stop("No ATT^o data for period ", pp, ".")
    }
    center_val <- att_o_pp$att_o_bin[1]
    y_label    <- expression(ATT^o)
    plot_title <- sprintf("Sensitivity in M: ATT^o (horizon %d)", horizon_t)
  }

  sens_df <- data.frame(
    param    = grid,
    ratio    = if (use_ratio) grid / param_hat else grid,
    center   = center_val,
    is_lower = center_val - mult * grid,
    is_upper = center_val + mult * grid
  )

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = .data[[x_var]]))

  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$is_lower, ymax = .data$is_upper),
      fill = col_band, alpha = 0.35
    ) +
    ggplot2::geom_hline(yintercept = sens_df$center[1], color = col_line,
                         linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

  if (use_ratio) {
    x_label <- if (param_name == "B") {
      expression(B / hat(B))
    } else {
      expression(M / hat(M))
    }
    p <- p +
      ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                           color = col_marker, linewidth = 0.7) +
      ggplot2::annotate("text", x = 1,
                         y = max(sens_df$is_upper, na.rm = TRUE),
                         label = sprintf("hat(%s)", param_name),
                         parse = TRUE,
                         vjust = -0.5, hjust = -0.1, color = col_marker,
                         size = 3.5) +
      ggplot2::labs(x = x_label, y = y_label, title = plot_title)
  } else {
    p <- p +
      ggplot2::labs(x = sprintf("Sensitivity parameter %s", param_name),
                    y = y_label, title = plot_title)
  }

  p + ggplot2::theme_minimal(base_size = 12)
}


plot_eventstudy <- function(x, col_band, col_pretrend) {
  if (!x$has_untreated) {
    stop("Event study plot requires untreated units (dose = 0).")
  }
  if (is.null(x$att_o)) {
    stop("Event study plot requires att_o (run with untreated units).")
  }

  # --- Post-period data (point estimate + IS at primary M) ---
  att_o_M <- x$att_o[x$att_o$M == x$M_hat, ]
  post_df <- data.frame(
    period   = att_o_M$period,
    estimate = att_o_M$att_o_bin,
    lower    = att_o_M$att_o_lower,
    upper    = att_o_M$att_o_upper,
    group    = "Post"
  )

  # --- Pre-period data (point estimates only) ---
  pre_df <- NULL
  if (!is.null(x$pre_att_o)) {
    pre_df <- data.frame(
      period   = x$pre_att_o$period,
      estimate = x$pre_att_o$att_o_bin,
      lower    = NA_real_,
      upper    = NA_real_,
      group    = "Pre"
    )
  }

  # --- Reference period (normalized to zero) ---
  ref_df <- data.frame(
    period   = x$ref_period,
    estimate = 0,
    lower    = NA_real_,
    upper    = NA_real_,
    group    = "Reference"
  )

  plot_df <- rbind(pre_df, ref_df, post_df)
  plot_df$group <- factor(plot_df$group, levels = c("Pre", "Reference", "Post"))

  # --- Build plot ---
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$period, y = .data$estimate))

  # IS bars for post-periods only
  post_rows <- plot_df$group == "Post"
  if (any(post_rows)) {
    p <- p +
      ggplot2::geom_errorbar(
        data = plot_df[post_rows, ],
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        width = 0.2, color = col_band, linewidth = 0.7
      )
  }

  p <- p +
    # Points colored by group
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$group, color = .data$group),
      size = 2.5
    ) +
    # Zero reference line
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    # Treatment onset line
    ggplot2::geom_vline(
      xintercept = min(x$post_periods) - 0.5,
      linetype = "dotted", color = "grey40", linewidth = 0.6
    ) +
    ggplot2::scale_color_manual(
      values = c("Pre" = col_pretrend, "Reference" = "grey40", "Post" = col_band),
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Pre" = 16, "Reference" = 1, "Post" = 16),
      guide = "none"
    ) +
    ggplot2::labs(
      x = "Period",
      y = expression(ATT^o ~ (binary ~ DiD)),
      title = "Event Study",
      subtitle = sprintf("Identified sets under level bound M = %.3f", x$M_hat)
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$period))) +
    ggplot2::theme_minimal(base_size = 12)

  p
}
