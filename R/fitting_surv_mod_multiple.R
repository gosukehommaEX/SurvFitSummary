#' Fit Multiple Parametric Survival Models with Time Scale Conversion
#'
#' This function fits multiple parametric survival models to time-to-event data
#' using various distributions. It supports both dependent and independent models
#' simultaneously, converts time scales, and automatically generates landmark
#' survival probabilities at 4-week intervals. Cox PH model is always fitted
#' regardless of dependent parameter value.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distribution Character vector specifying parametric distributions.
#'   Must be a subset of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param dependent Logical vector or character string. If c(TRUE, FALSE) or "both",
#'   fits both dependent and independent models. If single logical value (TRUE or FALSE),
#'   fits only that model type. Default is c(TRUE, FALSE)
#' @param control_arm Character string specifying which ARM level should be used as
#'   the reference (control) group. If NULL (default), uses the first factor level
#' @param time_scale Character string specifying time scale of input data.
#'   Must be one of: "day", "week", "month", "year". Default is "month".
#'   Data is assumed to be in this unit and will be converted to weeks for analysis
#' @param time_horizon Numeric value specifying the time horizon (in original time_scale units)
#'   for landmark survival probabilities and RMST calculation. If NULL (default),
#'   uses maximum observed survival time. Automatically converted to weeks
#'
#' @return A list containing seven data frames:
#'   \describe{
#'     \item{coef_estimates}{Parameter estimates with distribution and model type}
#'     \item{varcovmat}{Variance-covariance matrices in long format}
#'     \item{cholesky}{Cholesky decompositions in long format}
#'     \item{gof_statistics}{Goodness-of-fit statistics (AIC, BIC)}
#'     \item{landmark_survival}{Landmark survival probabilities at 4-week intervals}
#'     \item{survival_metrics}{Median, mean, and restricted mean survival times (in weeks)}
#'     \item{cox_ph_results}{Cox proportional hazards model results with HR and 95% CI}
#'   }
#'
#' @details
#' Time scale conversion: Data is assumed to be in the specified time_scale units.
#' All analysis is performed in weeks (4-week landmark intervals). Landmark survival
#' probabilities are automatically generated at 4-week intervals up to time_horizon.
#'
#' Cox PH model is always fitted and returned, regardless of the dependent parameter.
#' This allows comparison of parametric and semi-parametric estimates of treatment effect.
#'
#' @importFrom dplyr bind_rows mutate arrange select if_else
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom tidyr everything
#' @export
#'
#' @examples
#' \dontrun{
#' dataset <- generate_dummy_survival_data(
#'   arm = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazards = log(2) / c(20, 15),
#'   dropout_per_year = 0.1,
#'   seed = 123
#' )
#'
#' dataset_processed <- processing_dataset(
#'   dataset = dataset,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = NULL
#' )
#'
#' results <- fitting_surv_mod_multiple(
#'   dataset = dataset_processed,
#'   distribution = c("exp", "weibull", "gompertz"),
#'   dependent = "both",
#'   time_scale = "month",
#'   time_horizon = 60
#' )
#'
#' print(results$cox_ph_results)
#' }

fitting_surv_mod_multiple <- function(dataset,
                                      distribution,
                                      dependent = c(TRUE, FALSE),
                                      control_arm = NULL,
                                      time_scale = "month",
                                      time_horizon = NULL) {

  # Validate distribution parameter
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!all(distribution %in% valid_distributions)) {
    invalid <- distribution[!distribution %in% valid_distributions]
    stop("Invalid distribution(s): ", paste(invalid, collapse = ", "),
         "\nMust be subset of: ", paste(valid_distributions, collapse = ", "))
  }

  # Validate and convert time_scale
  scale_factors <- c(
    "day" = 1 / 7,
    "week" = 1,
    "month" = 30.44 / 7,
    "year" = 365.25 / 7
  )

  if (!time_scale %in% names(scale_factors)) {
    stop("time_scale must be one of: ", paste(names(scale_factors), collapse = ", "))
  }

  scale_factor <- scale_factors[[time_scale]]

  # Convert dataset to weeks
  dataset_weeks <- dataset %>%
    dplyr::mutate(SURVTIME = SURVTIME * scale_factor)

  # Determine time_horizon in weeks
  if (is.null(time_horizon)) {
    time_horizon_weeks <- max(dataset_weeks$SURVTIME, na.rm = TRUE)
  } else {
    time_horizon_weeks <- time_horizon * scale_factor
  }

  # Handle "both" option for dependent
  if (length(dependent) == 1 && is.character(dependent) && dependent == "both") {
    dependent <- c(TRUE, FALSE)
  }

  # Ensure dependent is logical vector
  if (!is.logical(dependent)) {
    stop("dependent must be logical vector or 'both'")
  }

  # Remove FALSE if ARM has only one level
  arm_levels <- unique(dataset_weeks$ARM)
  if (length(arm_levels) == 1 && TRUE %in% dependent) {
    message("ARM has only one level. Removing dependent = TRUE from analysis.")
    dependent <- FALSE
  }

  # Auto-generate landmark_times in weeks (4-week intervals)
  landmark_times <- seq(4, floor(time_horizon_weeks / 4) * 4, by = 4)
  if (length(landmark_times) == 0) {
    landmark_times <- c(4)
  }
  message("Time scale: ", time_scale, " | Data converted to weeks for analysis")
  message("Time horizon: ", time_horizon_weeks, " weeks (from ", time_horizon %||% "max observed", " ", time_scale, ")")
  message("Landmark survival probabilities at 4-week intervals: ", paste(landmark_times, collapse = ", "), " weeks")

  # Create grid of parameters to iterate over
  param_grid <- expand.grid(
    distribution = distribution,
    dependent = dependent,
    stringsAsFactors = FALSE
  )

  # Initialize result lists
  all_coef_estimates <- list()
  all_varcovmat <- list()
  all_cholesky <- list()
  all_gof_stats <- list()
  all_landmark <- list()
  all_metrics <- list()

  # Iterate over each parameter combination
  for (i in seq_len(nrow(param_grid))) {
    dist <- param_grid$distribution[i]
    dep <- param_grid$dependent[i]

    message("Fitting distribution: ", dist, ", dependent: ", dep)

    # Fit model
    fit_result <- tryCatch({
      fitting_surv_mod(
        dataset = dataset_weeks,
        distribution = dist,
        dependent = dep,
        control_arm = control_arm,
        landmark_times = landmark_times,
        time_horizon = time_horizon_weeks
      )
    }, error = function(e) {
      stop("FITTING ERROR - distribution: ", dist, ", dependent: ", dep,
           "\nError message: ", e$message, call. = FALSE)
    })

    # Add model_type column
    model_type <- ifelse(dep, "dependent", "independent")

    # Helper function to add model_type and ARM columns safely
    add_model_arm_columns <- function(df, model_type) {
      df <- df %>% dplyr::mutate(model_type = model_type, .after = distribution)

      if ("ARM" %in% colnames(df)) {
        df <- df %>%
          dplyr::mutate(ARM = dplyr::if_else(is.na(ARM), "Overall", as.character(ARM)))
      } else {
        df <- df %>%
          dplyr::mutate(ARM = "Overall", .after = model_type)
      }
      return(df)
    }

    # Store coefficient estimates
    coef_est <- add_model_arm_columns(fit_result$coef_estimates, model_type)
    all_coef_estimates[[i]] <- coef_est

    # Store variance-covariance matrix
    varcov <- add_model_arm_columns(fit_result$varcovmat, model_type)
    all_varcovmat[[i]] <- varcov

    # Store Cholesky decomposition
    chol <- add_model_arm_columns(fit_result$cholesky, model_type)
    all_cholesky[[i]] <- chol

    # Store goodness-of-fit statistics
    gof <- add_model_arm_columns(fit_result$gof_statistics, model_type)
    all_gof_stats[[i]] <- gof

    # Store landmark survival (if available)
    if (!is.null(fit_result$landmark_survival)) {
      lm <- add_model_arm_columns(fit_result$landmark_survival, model_type)
      all_landmark[[i]] <- lm
    }

    # Store survival metrics
    metrics <- add_model_arm_columns(fit_result$survival_metrics, model_type)
    all_metrics[[i]] <- metrics

    message("Successfully fitted: ", dist, " (", model_type, ")")
  }

  # Combine results - safe arrange and select
  coef_estimates_combined <- dplyr::bind_rows(all_coef_estimates) %>%
    dplyr::arrange(distribution, model_type, ARM) %>%
    dplyr::select(distribution, model_type, ARM, tidyr::everything())

  varcovmat_combined <- dplyr::bind_rows(all_varcovmat) %>%
    dplyr::arrange(distribution, model_type, ARM) %>%
    dplyr::select(distribution, model_type, ARM, tidyr::everything())

  cholesky_combined <- dplyr::bind_rows(all_cholesky) %>%
    dplyr::arrange(distribution, model_type, ARM) %>%
    dplyr::select(distribution, model_type, ARM, tidyr::everything())

  gof_stats_combined <- dplyr::bind_rows(all_gof_stats) %>%
    dplyr::arrange(distribution, model_type, ARM) %>%
    dplyr::select(distribution, model_type, ARM, tidyr::everything())

  landmark_combined <- dplyr::bind_rows(all_landmark) %>%
    dplyr::arrange(distribution, model_type, ARM) %>%
    dplyr::select(distribution, model_type, ARM, tidyr::everything())

  metrics_combined <- dplyr::bind_rows(all_metrics) %>%
    dplyr::arrange(distribution, model_type, ARM) %>%
    dplyr::select(distribution, model_type, ARM, tidyr::everything())

  # ========== Cox PH model fitting (embedded) ==========
  message("Fitting Cox proportional hazards model...")

  # Ensure ARM is a factor with proper reference level
  if (!is.factor(dataset_weeks$ARM)) {
    dataset_weeks$ARM <- as.factor(dataset_weeks$ARM)
  }

  # Fit Cox model
  cox_fit <- tryCatch({
    survival::coxph(
      survival::Surv(SURVTIME, EVENT) ~ ARM,
      data = dataset_weeks
    )
  }, error = function(e) {
    warning("Cox PH fitting failed: ", e$message)
    return(NULL)
  })

  cox_ph_result <- NULL

  if (!is.null(cox_fit)) {
    # Extract results with confidence intervals
    cox_ph_result <- tryCatch({
      broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
        dplyr::mutate(
          HR = round(estimate, 3),
          `L95%` = round(conf.low, 3),
          `U95%` = round(conf.high, 3),
          `p.value` = round(p.value, 4)
        ) %>%
        dplyr::select(term, HR, `L95%`, `U95%`, `p.value`) %>%
        dplyr::mutate(
          event_count = sum(dataset$EVENT, na.rm = TRUE),
          sample_size = nrow(dataset)
        )
    }, error = function(e) {
      # Fallback: extract from summary.coxph
      cox_summary <- summary(cox_fit)
      cox_ci <- cox_summary$conf.int

      cox_ph_result <- data.frame(
        term = rownames(cox_ci),
        HR = round(cox_ci[, 1], 3),
        `L95%` = round(cox_ci[, 3], 3),
        `U95%` = round(cox_ci[, 4], 3),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(
          `p.value` = round(cox_summary$coefficients[, "Pr(>|z|)"], 4)
        ) %>%
        dplyr::select(term, HR, `L95%`, `U95%`, `p.value`) %>%
        dplyr::mutate(
          event_count = sum(dataset$EVENT, na.rm = TRUE),
          sample_size = nrow(dataset)
        )
      return(cox_ph_result)
    })

    message("Successfully fitted: Cox PH model")
  } else {
    warning("Cox PH model could not be fitted")
  }

  # Return combined results
  return(list(
    coef_estimates = coef_estimates_combined,
    varcovmat = varcovmat_combined,
    cholesky = cholesky_combined,
    gof_statistics = gof_stats_combined,
    landmark_survival = landmark_combined,
    survival_metrics = metrics_combined,
    cox_ph_results = cox_ph_result
  ))
}


# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
