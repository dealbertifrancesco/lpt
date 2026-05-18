#' French Social Housing Data (SRU Law)
#'
#' Panel data on French communes affected by the Solidarity and Urban
#' Renewal (SRU) law, which imposed a 20% social housing target.
#' Communes below the threshold faced continuous treatment intensity
#' proportional to their gap from the target. Periods have been
#' re-indexed so that treatment onset corresponds to period 0.
#'
#' @format A data frame with 12,883 rows and 4 columns:
#' \describe{
#'   \item{commune}{Character. Commune geographic identifier.}
#'   \item{year}{Integer. Period index (-7 to 5). Pre-treatment:
#'     -7 to -1; post-treatment: 0 to 5. Treatment onset at period 0.}
#'   \item{outcome}{Numeric. Social housing share, standardized to the
#'     reference period (period -1).}
#'   \item{dose}{Numeric. Treatment intensity: 0 for never-treated communes,
#'     gap from 20\% target (normalized to [0,1]) for treated communes.
#'     Time-invariant.}
#' }
#'
#' @details
#' The sample contains 991 communes: 358 treated and 633 never-treated
#' controls. Pre-treatment periods are -7 to -1 (7 periods, 6 consecutive
#' pairs for B calibration); post-treatment periods are 0 to 5 (6 periods,
#' horizons 0 through 5).
#'
#' Recommended usage with \code{\link{lpt}}:
#' \preformatted{
#' data(sru)
#' fit <- lpt(sru, "commune", "year", "outcome", "dose",
#'            post_period = 0:5, pre_periods = -7:-1,
#'            B = "calibrate")
#' }
#'
#' @source Inspired by Dealberti (2026), Local Parallel Trends for
#'   Continuous Difference-in-Differences.
#'
#' @name sru
#' @docType data
#' @usage data(sru)
NULL
