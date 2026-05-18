#' Run contdid estimation backend
#'
#' Internal function that calls \code{contdid::cont_did()} and extracts
#' per-period dose-response estimates into lpt's data structures.
#'
#' @param data Data frame in long format.
#' @param id_col,time_col,outcome_col,dose_col Column name strings.
#' @param post_periods Sorted numeric vector of post-period identifiers.
#' @param pre_period_set Sorted numeric vector of pre-period identifiers.
#' @param min_post Minimum post-period value.
#' @param has_untreated Logical. Whether untreated units exist.
#' @param dose_vec Numeric vector of doses from the reference period.
#' @param B Sensitivity parameter specification ("calibrate" or numeric).
#' @param contdid_args List of additional arguments for \code{cont_did()}.
#'
#' @return A list with components: \code{datt}, \code{att}, \code{att_o},
#'   \code{slopes}, \code{calibration}, \code{B_hat}, \code{B_values},
#'   \code{contdid_fit}.
#'
#' @references
#' Callaway B, Goodman-Bacon A, Sant'Anna PHC (2024).
#' \dQuote{Difference-in-differences with a continuous treatment.}
#' \emph{National Bureau of Economic Research}.
#'
#' @keywords internal
run_contdid <- function(data, id_col, time_col, outcome_col, dose_col,
                        post_periods, pre_period_set, min_post,
                        has_untreated, dose_vec, B, contdid_args) {

  if (!requireNamespace("contdid", quietly = TRUE)) {
    stop("Package 'contdid' required for method = 'contdid'. ",
         "Install with: remotes::install_github('bcallaway11/contdid')")
  }

  all_periods <- sort(unique(data[[time_col]]))

  # --- Data prep: contdid needs positive integer periods + gname column ---
  shift <- 1L - as.integer(min(all_periods))
  min_post_shifted <- as.integer(min_post) + shift

  cd_data <- data.frame(data, check.names = FALSE)
  cd_data[[".cd_period"]] <- as.integer(data[[time_col]]) + shift
  cd_data[[".cd_G2p"]] <- ifelse(data[[dose_col]] > 0, min_post_shifted, 0L)

  # --- Build cont_did() arguments ---
  cd_args <- list(
    yname            = outcome_col,
    tname            = ".cd_period",
    idname           = id_col,
    dname            = dose_col,
    gname            = ".cd_G2p",
    data             = cd_data,
    target_parameter = "slope",
    aggregation      = "eventstudy",
    treatment_type   = "continuous",
    control_group    = "notyettreated"
  )

  # Defaults for contdid-specific args
  cd_defaults <- list(biters = 500, cband = TRUE, num_knots = 1, degree = 3)
  for (nm in names(cd_defaults)) {
    cd_args[[nm]] <- if (!is.null(contdid_args[[nm]])) {
      contdid_args[[nm]]
    } else {
      cd_defaults[[nm]]
    }
  }
  # Forward any extra contdid args the user supplied
  for (nm in setdiff(names(contdid_args), names(cd_defaults))) {
    cd_args[[nm]] <- contdid_args[[nm]]
  }

  # --- Call cont_did ---
  cd_fit <- do.call(contdid::cont_did, cd_args)

  # --- Extract per-period data from extra_gt_returns ---
  att_gt <- cd_fit$att_gt
  n_gt <- length(att_gt$group)
  eval_points <- as.numeric(cd_fit$ptep$dvals)

  period_data <- list()          # keyed by original lpt period
  pre_slopes_list <- list()      # for calibration / pretrends plot
  sup_vals <- numeric(0)
  pre_labels <- character(0)

  for (i in seq_len(n_gt)) {
    orig_period <- as.numeric(att_gt$t[i]) - shift
    extra <- att_gt$extra_gt_returns[[i]]$extra_gt_returns

    if (is.null(extra)) {
      warning(sprintf("No dose-resolved data for period %s (contdid returned NULL). Skipping.",
                      orig_period))
      next
    }

    period_data[[as.character(orig_period)]] <- list(
      acrt_d       = as.numeric(extra$acrt.d),
      att_d        = as.numeric(extra$att.d),
      att_overall  = as.numeric(extra$att.overall),
      acrt_overall = as.numeric(extra$acrt.overall)
    )

    # Pre-period: collect for B calibration
    if (orig_period < min_post) {
      lbl <- as.character(orig_period)
      pre_labels <- c(pre_labels, lbl)
      sup_val <- max(abs(extra$acrt.d))
      sup_vals <- c(sup_vals, sup_val)

      pre_slopes_list[[length(pre_slopes_list) + 1]] <- data.frame(
        period_pair = lbl,
        d           = eval_points,
        mu_prime_d  = as.numeric(extra$acrt.d)
      )
    }
  }

  if (length(period_data) == 0) {
    stop("contdid returned no dose-resolved data for any period.")
  }

  # --- Handle B ---
  calibration_result <- NULL

  if (is.character(B) && B == "calibrate") {
    if (length(sup_vals) < 1) {
      stop("B = 'calibrate' with method = 'contdid' requires pre-period data. ",
           "Ensure pre_periods covers at least two pre-treatment periods.")
    }
    names(sup_vals) <- pre_labels

    pre_slopes <- do.call(rbind, pre_slopes_list)
    rownames(pre_slopes) <- NULL

    B_hat <- max(sup_vals)
    B_values <- B_hat

    calibration_result <- list(
      B_hat          = B_hat,
      pre_slopes     = pre_slopes,
      sup_by_period  = sup_vals
    )

    message(sprintf("Calibrated B = %.4f from %d pre-period(s) [contdid].",
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

  # --- Construct identified sets per post-period ---
  D_bar <- mean(dose_vec[dose_vec > 0])
  datt_all <- list()
  att_all  <- list()
  att_o_all <- list()

  for (pp in post_periods) {
    pd <- period_data[[as.character(pp)]]
    if (is.null(pd)) next

    horizon_t <- as.integer(pp - min_post)
    mult <- horizon_t + 1L

    # IS_{dATT}
    for (b in B_values) {
      datt_all[[length(datt_all) + 1]] <- data.frame(
        period     = pp,
        horizon    = horizon_t,
        d          = eval_points,
        lambda_d   = pd$acrt_d,
        B          = b,
        datt_lower = pd$acrt_d - mult * b,
        datt_upper = pd$acrt_d + mult * b
      )
    }

    # IS_{ATT}  (Lambda_d comes directly from contdid's att.d)
    if (has_untreated) {
      for (b in B_values) {
        att_all[[length(att_all) + 1]] <- data.frame(
          period    = pp,
          horizon   = horizon_t,
          d         = eval_points,
          Lambda_d  = pd$att_d,
          B         = b,
          att_lower = pd$att_d - mult * b * eval_points,
          att_upper = pd$att_d + mult * b * eval_points
        )
      }
    }

    # IS_{ATT^o}
    if (has_untreated) {
      for (b in B_values) {
        att_o_all[[length(att_o_all) + 1]] <- data.frame(
          period     = pp,
          horizon    = horizon_t,
          att_o_bin  = pd$att_overall,
          D_bar      = D_bar,
          B          = b,
          att_o_lower = pd$att_overall - mult * b * D_bar,
          att_o_upper = pd$att_overall + mult * b * D_bar
        )
      }
    }
  }

  datt <- do.call(rbind, datt_all)
  att  <- if (length(att_all)  > 0) do.call(rbind, att_all)  else NULL
  att_o <- if (length(att_o_all) > 0) do.call(rbind, att_o_all) else NULL
  rownames(datt) <- NULL
  if (!is.null(att))  rownames(att) <- NULL
  if (!is.null(att_o)) rownames(att_o) <- NULL

  # --- Build slopes list (plot/summary compatibility) ---
  slopes <- list()
  for (pp in post_periods) {
    pd <- period_data[[as.character(pp)]]
    if (is.null(pd)) next
    slopes[[as.character(pp)]] <- structure(
      list(
        eval_points      = eval_points,
        conditional_mean = pd$att_d,
        lambda_d         = pd$acrt_d,
        gam_fit          = NULL
      ),
      class = "dose_slope"
    )
  }

  list(
    datt        = datt,
    att         = att,
    att_o       = att_o,
    slopes      = slopes,
    calibration = calibration_result,
    B_hat       = B_hat,
    B_values    = B_values,
    contdid_fit = cd_fit
  )
}
