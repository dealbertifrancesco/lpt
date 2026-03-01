#' Simulate data for the LPT framework (internal)
#'
#' Generates panel data with known ground truth for testing.
#'
#' @param n Integer. Number of units (default: 1000).
#' @param n_pre Integer. Number of pre-treatment periods (default: 3).
#' @param n_post Integer. Number of post-treatment periods (default: 1).
#' @param dgp Character. One of \code{"no_selection"}, \code{"linear_selection"},
#'   \code{"time_varying_selection"}.
#' @param beta Numeric vector of length 2. Dose-response coefficients.
#' @param rho Numeric. Selection strength (default: 0.3).
#' @param rho2 Numeric. Quadratic time coefficient for \code{"time_varying_selection"}
#'   (default: 0.05). Ignored for other DGPs.
#' @param sigma_eps Numeric. Noise SD (default: 1).
#' @param seed Integer or NULL. Random seed.
#'
#' @return A data frame with columns: \code{id}, \code{time}, \code{outcome},
#'   \code{dose}, \code{treated}. Attributes: \code{true_datt}, \code{true_att}
#'   (functions), \code{true_B} (scalar), \code{true_delta_tilde_0} (function,
#'   for time-varying DGPs), \code{dgp} (character).
#'
#' @keywords internal
simulate_lpt <- function(n = 1000, n_pre = 3, n_post = 1,
                          dgp = c("no_selection", "linear_selection",
                                  "time_varying_selection"),
                          beta = c(1, -0.5), rho = 0.3, rho2 = 0.05,
                          sigma_eps = 1, seed = NULL) {
  dgp <- match.arg(dgp)
  stopifnot(is.numeric(n), length(n) == 1, n >= 10, n == floor(n))
  stopifnot(is.numeric(n_pre), length(n_pre) == 1, n_pre >= 1, n_pre == floor(n_pre))
  stopifnot(is.numeric(n_post), length(n_post) == 1, n_post >= 1, n_post == floor(n_post))
  stopifnot(is.numeric(beta), length(beta) == 2)
  stopifnot(is.numeric(rho), length(rho) == 1, is.finite(rho))
  stopifnot(is.numeric(sigma_eps), length(sigma_eps) == 1, sigma_eps > 0)

  if (!is.null(seed)) set.seed(seed)

  n_periods <- n_pre + n_post

  # Unit-level
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  alpha_i <- stats::rnorm(n)

  # Dose: exponential with point mass at zero
  latent_dose <- exp(0.5 * x1 + 0.3 * x2 + stats::rnorm(n, sd = 0.5))
  untreated <- stats::rbinom(n, 1, 0.2) == 1
  dose <- ifelse(untreated, 0, latent_dose)

  if (any(dose > 0)) {
    d_max <- stats::quantile(dose[dose > 0], 0.95)
    dose <- pmin(dose, d_max)
  }

  # Selection function: E[Y_0(t) | D=d] - E[Y_0(t) | D=0] = selection_fn(d, t)
  selection_fn <- switch(dgp,
    "no_selection"           = function(d, t) rep(0, length(d)),
    "linear_selection"       = function(d, t) rho * d * t,
    "time_varying_selection" = function(d, t) rho * d * t + rho2 * d * t^2
  )

  # Treatment effect
  tau <- function(d) beta[1] * d + beta[2] * d^2
  tau_prime <- function(d) beta[1] + 2 * beta[2] * d

  # LPT bound: max derivative of selection_fn w.r.t. d
  # For linear_selection: d/dd [rho*d*t] = rho  => true_B = |rho|
  # For time_varying_selection: d/dd [...] = rho*t + rho2*t^2 -- but B is the
  #   Lipschitz bound on the *cross-sectional slope* not time trends; use |rho|
  true_B <- switch(dgp,
    "no_selection"           = 0,
    "linear_selection"       = abs(rho),
    "time_varying_selection" = abs(rho)
  )

  # Generate panel
  gamma <- 1
  records <- vector("list", n_periods)

  for (t_idx in seq_len(n_periods)) {
    eps_it <- stats::rnorm(n, sd = sigma_eps)
    y0 <- alpha_i + gamma * t_idx + selection_fn(dose, t_idx) + eps_it

    if (t_idx > n_pre) {
      y_obs <- y0 + tau(dose)
    } else {
      y_obs <- y0
    }

    records[[t_idx]] <- data.frame(
      id = seq_len(n),
      time = t_idx,
      outcome = y_obs,
      dose = dose,
      treated = dose > 0
    )
  }

  result <- do.call(rbind, records)
  rownames(result) <- NULL

  attr(result, "true_datt") <- tau_prime
  attr(result, "true_att")  <- tau
  attr(result, "true_B")    <- true_B
  attr(result, "dgp")       <- dgp

  # true_delta_tilde_0(d): increment in selection from period (n_pre-1) to n_pre
  # delta_tilde_s(d) = selection_fn(d, s) - selection_fn(d, s-1) for consecutive s
  # Last pre-period pair: s = n_pre
  #   linear_selection:       rho*d*n_pre - rho*d*(n_pre-1)     = rho * d
  #   time_varying_selection: [rho*d*n_pre + rho2*d*n_pre^2]
  #                         - [rho*d*(n_pre-1) + rho2*d*(n_pre-1)^2]
  #                         = rho*d + rho2*d*(2*n_pre - 1)
  true_delta_tilde_0 <- if (dgp == "time_varying_selection") {
    function(d) rho * d + rho2 * d * (2 * n_pre - 1)
  } else {
    function(d) rho * d   # also correct for no_selection (rho=0)
  }
  attr(result, "true_delta_tilde_0") <- true_delta_tilde_0

  result
}
