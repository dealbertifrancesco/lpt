#' Plot LPT results
#'
#' Creates visualizations for identified sets, pre-trend diagnostics,
#' and sensitivity analysis.
#'
#' @param x An object of class \code{"lpt"}.
#' @param type Character. Plot type: \code{"datt"} (default), \code{"att"},
#'   \code{"pretrends"}, \code{"sensitivity"}.
#' @param period Scalar or NULL. Which post-period to plot. If NULL and
#'   multiple post-periods exist, plots all periods faceted. Default: NULL.
#' @param d0 Numeric. Required dose value when \code{estimand} is \code{"datt"}
#'   or \code{"att"}. The identified set for these estimands depends on the
#'   chosen dose level; no default is provided. Ignored for \code{"att_o"}.
#' @param B_grid Numeric vector or NULL. Grid of B values for sensitivity plot.
#'   Default: \code{seq(0, 3 * B_hat, length.out = 50)}.
#' @param estimand Character. Estimand for the sensitivity plot:
#'   \code{"att_o"} (default, overall ATT \eqn{ATT^o}; no \code{d0} needed),
#'   \code{"datt"} (dose-response slope \eqn{\partial ATT/\partial d} at
#'   \code{d0}), or \code{"att"} (ATT level \eqn{ATT(d|d)} at \code{d0}).
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
#'            post_period = 2019, pre_periods = 1993:1999, B = 0)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(fit, type = "datt")
#' }
#'
#' @method plot lpt
#' @export
plot.lpt <- function(x, type = c("datt", "att", "pretrends", "sensitivity"),
                      period = NULL, d0 = NULL, B_grid = NULL,
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
    "sensitivity" = plot_sensitivity(x, d0, B_grid, col_band, col_line, col_marker,
                                      estimand = estimand, period = period)
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
      title = "Identified Set for Dose-Response Slope"
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
    stop("ATT bounds not available (no untreated units in data).")
  }

  df <- x$att
  if (!is.null(period)) {
    df <- df[df$period %in% period, ]
  }

  multi_period <- length(unique(df$period)) > 1

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$att_lower, ymax = .data$att_upper,
                   group = factor(.data$B)),
      fill = col_band, alpha = 0.3
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$Lambda_d),
                        color = col_line, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      x = "Dose (d)",
      y = expression(ATT(d)),
      title = "Identified Set for ATT(d|d)"
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
  if (is.null(x$calibration)) {
    stop("Pre-trends plot requires calibration (B = 'calibrate').")
  }

  df <- x$calibration$pre_slopes
  B_hat <- x$B_hat

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d, y = .data$mu_prime_d,
                                          color = .data$period_pair)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$mu_prime_d - 1.96 * .data$se,
                   ymax = .data$mu_prime_d + 1.96 * .data$se,
                   fill = .data$period_pair),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_hline(yintercept = c(-B_hat, B_hat),
                         linetype = "dashed", color = col_marker, linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "grey40",
                         linewidth = 0.3) +
    ggplot2::annotate("text", x = min(df$d), y = B_hat,
                       label = sprintf("B = %.3f", B_hat),
                       vjust = -0.5, hjust = 0, color = col_marker, size = 3.5) +
    ggplot2::labs(
      x = "Dose (d)",
      y = expression(hat(mu) * "'" * (d)),
      title = "Pre-Period Selection Slopes",
      color = "Period pair",
      fill = "Period pair"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  p
}


plot_sensitivity_rm <- function(x, d0, B_grid, col_band, col_line, col_marker,
                                 estimand = "att_o", period = NULL) {
  # RM-Time sensitivity: x-axis is M_bar; IS = [att_o_bin +/- t * M * delta_bar_star]
  if (estimand != "att_o") {
    stop("For pure RM-Time sensitivity (B=0), only estimand = 'att_o' is supported.")
  }
  if (is.null(x$att_o)) {
    stop("ATT^o bounds not available (no untreated units in data).")
  }

  pp    <- if (is.null(period)) x$post_periods[1] else period
  t_val <- x$t_values[[as.character(pp)]]
  M_hat <- x$specifications$M_bar

  att_o_pp <- x$att_o[x$att_o$period == pp & x$att_o$B == 0, ]
  if (nrow(att_o_pp) == 0) {
    stop("No ATT^o data for period ", pp, " with B = 0.")
  }
  center_val <- att_o_pp$att_o_bin[1]

  # Recover delta_bar_star_pre from stored bounds (valid since B=0 means RM only)
  hw_stored      <- (att_o_pp$att_o_upper[1] - att_o_pp$att_o_lower[1]) / 2
  if (t_val > 0 && M_hat > 0) {
    delta_bar_star <- hw_stored / (t_val * M_hat)
  } else {
    delta_bar_star <- 0
  }

  M_grid  <- seq(0, 3 * max(M_hat, 1e-6), length.out = 50)
  sens_df <- data.frame(
    M_bar    = M_grid,
    is_lower = center_val - t_val * M_grid * delta_bar_star,
    is_upper = center_val + t_val * M_grid * delta_bar_star
  )

  ggplot2::ggplot(sens_df, ggplot2::aes(x = .data$M_bar)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$is_lower, ymax = .data$is_upper),
      fill = col_band, alpha = 0.35
    ) +
    ggplot2::geom_hline(yintercept = center_val, color = col_line, linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = M_hat, linetype = "dotted",
                         color = col_marker, linewidth = 0.7) +
    ggplot2::annotate("text", x = M_hat,
                       y = max(sens_df$is_upper, na.rm = TRUE),
                       label = expression(hat(M)[bar]),
                       vjust = -0.5, hjust = -0.1, color = col_marker, size = 3.5) +
    ggplot2::labs(
      x = "M_bar",
      y = expression(ATT^o),
      title = sprintf("Sensitivity: ATT^o (RM-Time, t = %d)", t_val)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


plot_sensitivity <- function(x, d0, B_grid, col_band, col_line, col_marker,
                              estimand = "datt", period = NULL) {

  # Dispatch to RM-Time helper when appropriate
  tr <- x$specifications$time_restriction
  if (!is.null(tr) && tr == "rm" && x$B_hat == 0) {
    return(plot_sensitivity_rm(x, d0, B_grid, col_band, col_line, col_marker,
                               estimand = estimand, period = period))
  }

  # Select period
  pp <- if (is.null(period)) x$post_periods[1] else period
  pp_char <- as.character(pp)
  sr <- x$slopes[[pp_char]]
  ep <- sr$eval_points

  if (estimand %in% c("datt", "att") && is.null(d0)) {
    stop("'d0' must be specified when estimand = '", estimand, "'. ",
         "Supply a dose value, e.g. d0 = ", round(stats::median(ep), 3), ".")
  }

  # Build B_grid (absolute values; optionally expressed as B/B_hat on axis)
  if (is.null(B_grid)) {
    if (x$B_hat > 0) {
      B_grid <- seq(0, 3, length.out = 50) * x$B_hat
    } else {
      B_grid <- seq(0, 1, length.out = 50)
    }
  }

  use_ratio <- x$B_hat > 0
  z_alpha   <- stats::qnorm(1 - x$specifications$alpha / 2)
  x_var     <- if (use_ratio) "B_ratio" else "B"

  if (estimand == "datt") {
    # IS_dATT(d0; B) = [lambda(d0) - B, lambda(d0) + B]
    idx          <- which.min(abs(ep - d0))
    center_val   <- sr$lambda_d[idx]
    se_val       <- sr$se_lambda[idx]
    d0_actual    <- ep[idx]

    sens_df <- data.frame(
      B       = B_grid,
      B_ratio = if (use_ratio) B_grid / x$B_hat else B_grid,
      center  = center_val,
      is_lower = center_val - B_grid,
      is_upper = center_val + B_grid,
      ci_lower = center_val - B_grid - z_alpha * se_val,
      ci_upper = center_val + B_grid + z_alpha * se_val
    )
    y_label    <- expression(partialdiff * ATT(d) / partialdiff * d ~ "at" ~ d[0])
    plot_title <- sprintf("Sensitivity: Dose-Response Slope at d = %.2f", d0_actual)
    show_ci    <- TRUE

  } else if (estimand == "att") {
    # IS_ATT(d0; B) = [Lambda(d0) - B*d0, Lambda(d0) + B*d0]
    if (is.null(x$att)) {
      stop("ATT level bounds not available (no untreated units in data). ",
           "Use estimand = 'datt'.")
    }
    att_pp <- x$att[x$att$period == pp, ]
    if (nrow(att_pp) == 0) {
      stop("No ATT data for period ", pp, ".")
    }
    # Lambda_d is constant across B; use any B value
    att_sub      <- att_pp[att_pp$B == att_pp$B[1], ]
    idx          <- which.min(abs(att_sub$d - d0))
    center_val   <- att_sub$Lambda_d[idx]
    d0_actual    <- att_sub$d[idx]

    sens_df <- data.frame(
      B        = B_grid,
      B_ratio  = if (use_ratio) B_grid / x$B_hat else B_grid,
      center   = center_val,
      is_lower = center_val - B_grid * d0_actual,
      is_upper = center_val + B_grid * d0_actual,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
    y_label    <- expression(ATT(d[0]))
    plot_title <- sprintf("Sensitivity: ATT(d|d) at d = %.2f", d0_actual)
    show_ci    <- FALSE

  } else {
    # estimand == "att_o"
    # IS_{ATT^o}(B) = [att_o_bin - B*D_bar, att_o_bin + B*D_bar]
    if (is.null(x$att_o)) {
      stop("Overall ATT bounds not available (no untreated units in data). ",
           "Use estimand = 'datt'.")
    }
    att_o_pp <- x$att_o[x$att_o$period == pp, ]
    if (nrow(att_o_pp) == 0) {
      stop("No ATT^o data for period ", pp, ".")
    }
    center_val <- att_o_pp$att_o_bin[1]
    D_bar      <- att_o_pp$D_bar[1]

    sens_df <- data.frame(
      B        = B_grid,
      B_ratio  = if (use_ratio) B_grid / x$B_hat else B_grid,
      center   = center_val,
      is_lower = center_val - B_grid * D_bar,
      is_upper = center_val + B_grid * D_bar,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
    y_label    <- expression(ATT^o)
    plot_title <- "Sensitivity: Overall ATT"
    show_ci    <- FALSE
  }

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = .data[[x_var]]))

  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = col_band, alpha = 0.2
    )
  }

  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$is_lower, ymax = .data$is_upper),
      fill = col_band, alpha = 0.35
    ) +
    ggplot2::geom_hline(yintercept = sens_df$center[1], color = col_line,
                         linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

  if (use_ratio) {
    p <- p +
      ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                           color = col_marker, linewidth = 0.7) +
      ggplot2::annotate("text", x = 1,
                         y = max(sens_df$is_upper, na.rm = TRUE),
                         label = expression(hat(B)),
                         vjust = -0.5, hjust = -0.1, color = col_marker,
                         size = 3.5) +
      ggplot2::labs(x = expression(B / hat(B)), y = y_label, title = plot_title)
  } else {
    p <- p +
      ggplot2::labs(x = "Sensitivity parameter B", y = y_label,
                    title = plot_title)
  }

  p + ggplot2::theme_minimal(base_size = 12)
}
