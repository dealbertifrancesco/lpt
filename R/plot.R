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
#' @param d0 Numeric or NULL. Dose value for the sensitivity plot.
#'   Default: median of evaluation grid.
#' @param B_grid Numeric vector or NULL. Grid of B values for sensitivity plot.
#'   Default: \code{seq(0, 2 * B_hat, length.out = 50)}.
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
                      true_curve = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.lpt(). ",
         "Install with install.packages('ggplot2').")
  }

  type <- match.arg(type)

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
    "sensitivity" = plot_sensitivity(x, d0, B_grid, col_band, col_line, col_marker)
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


plot_sensitivity <- function(x, d0, B_grid, col_band, col_line, col_marker) {
  ep <- x$slopes[[1]]$eval_points
  if (is.null(d0)) {
    d0 <- stats::median(ep)
  }

  idx <- which.min(abs(ep - d0))
  lambda_at_d0 <- x$slopes[[1]]$lambda_d[idx]
  se_at_d0 <- x$slopes[[1]]$se_lambda[idx]
  z_alpha <- stats::qnorm(1 - x$specifications$alpha / 2)

  # B_grid can be supplied in absolute terms; default is multiples of B_hat
  if (is.null(B_grid)) {
    if (x$B_hat > 0) {
      # Default: 0 to 3x calibrated B, in multiples
      mult_grid <- seq(0, 3, length.out = 50)
      B_grid <- mult_grid * x$B_hat
    } else {
      # B_hat = 0 (user-supplied): use absolute scale, no meaningful ratio
      B_grid <- seq(0, 1, length.out = 50)
    }
  }

  # Determine if we can express as B/B_hat multiples
  use_ratio <- x$B_hat > 0

  sens_df <- data.frame(
    B = B_grid,
    B_ratio = if (use_ratio) B_grid / x$B_hat else B_grid,
    datt_lower = lambda_at_d0 - B_grid,
    datt_upper = lambda_at_d0 + B_grid,
    ci_lower = lambda_at_d0 - B_grid - z_alpha * se_at_d0,
    ci_upper = lambda_at_d0 + B_grid + z_alpha * se_at_d0
  )

  x_var <- if (use_ratio) "B_ratio" else "B"

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = .data[[x_var]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = col_band, alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$datt_lower, ymax = .data$datt_upper),
      fill = col_band, alpha = 0.35
    ) +
    ggplot2::geom_hline(yintercept = lambda_at_d0, color = col_line,
                         linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")

  if (use_ratio) {
    # Mark B_hat at ratio = 1
    p <- p +
      ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                           color = col_marker, linewidth = 0.7) +
      ggplot2::annotate("text", x = 1, y = max(sens_df$ci_upper),
                         label = expression(hat(B)),
                         vjust = -0.5, hjust = -0.1, color = col_marker,
                         size = 3.5) +
      ggplot2::labs(
        x = expression(B / hat(B)),
        y = expression(partialdiff * ATT / partialdiff * d ~ "at" ~ d[0]),
        title = sprintf("Sensitivity Analysis at d = %.2f", d0)
      )
  } else {
    p <- p +
      ggplot2::labs(
        x = "Sensitivity parameter B",
        y = expression(partialdiff * ATT / partialdiff * d ~ "at" ~ d[0]),
        title = sprintf("Sensitivity Analysis at d = %.2f", d0)
      )
  }

  p + ggplot2::theme_minimal(base_size = 12)
}
