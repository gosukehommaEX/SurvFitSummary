#' Fit Multiple Parametric Survival Models
#'
#' This function fits multiple parametric survival models to time-to-event data
#' using various distributions. It supports both dependent and independent models
#' simultaneously. Cox PH model is always fitted regardless of dependent parameter value.
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
#' @param landmark_increment Numeric value specifying the time interval for landmark survival
#'   probability estimation. If NULL (default), landmark survival is not calculated.
#'   Results will be calculated at times 0, landmark_increment, 2*landmark_increment, etc.
#' @param time_horizon Numeric value specifying the time horizon for restricted mean
#'   survival time (RMST) calculation. If NULL (default), RMST is not calculated
#'
#' @return A list containing seven data frames:
#'   \describe{
#'     \item{coef_estimates}{Parameter estimates with distribution and model type}
#'     \item{varcovmat}{Variance-covariance matrices in long format}
#'     \item{cholesky}{Cholesky decompositions in long format}
#'     \item{gof_statistics}{Goodness-of-fit statistics (AIC, BIC)}
#'     \item{landmark_survival}{Landmark survival probabilities with columns:
#'       model_type, ARM, Time, KM, distribution1, distribution2, ...}
#'     \item{survival_metrics}{Median, mean, and restricted mean survival times}
#'     \item{cox_ph_results}{Cox proportional hazards model results with HR and 95% CI}
#'   }
#'
#' @details
#' Landmark survival probabilities are calculated for each distribution at specified time intervals
#' for each ARM. The Kaplan-Meier estimate is included for each ARM as a reference.
#' Results are combined across distributions using model_type, ARM, and Time as join keys.
#'
#' Cox PH model is always fitted and returned, regardless of the dependent parameter.
#' This allows comparison of parametric and semi-parametric estimates of treatment effect.
#'
#' @importFrom dplyr bind_rows mutate arrange select if_else left_join group_by summarise across
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom tidyr everything
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage with multiple distributions
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
#'   landmark_increment = 6,
#'   time_horizon = 36
#' )
#'
#' # Landmark survival with multiple distributions and ARMs
#' print(results$landmark_survival)
#'
#' # Survival metrics
#' print(results$survival_metrics)
#'
#' # Cox PH results
#' print(results$cox_ph_results)
#'
#' # Example 2: Multiple distributions with dependent models only
#' results2 <- fitting_surv_mod_multiple(
#'   dataset = dataset_processed,
#'   distribution = c("weibull", "lnorm", "gengamma"),
#'   dependent = TRUE,
#'   control_arm = "Control",
#'   landmark_increment = 3,
#'   time_horizon = 36
#' )
#'
#' print(results2$landmark_survival)
#' print(results2$survival_metrics)
#'
#' # Example 3: Single distribution
#' results3 <- fitting_surv_mod_multiple(
#'   dataset = dataset_processed,
#'   distribution = "exp",
#'   dependent = c(TRUE, FALSE),
#'   landmark_increment = 6,
#'   time_horizon = 36
#' )
#'
#' print(results3$landmark_survival)
#' }
fitting_surv_mod_multiple <- function(dataset,
                                      distribution,
                                      dependent = c(TRUE, FALSE),
                                      control_arm = NULL,
                                      landmark_increment = NULL,
                                      time_horizon = NULL) {

  # Validate distribution parameter
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!all(distribution %in% valid_distributions)) {
    invalid <- distribution[!distribution %in% valid_distributions]
    stop("Invalid distribution(s): ", paste(invalid, collapse = ", "),
         "\nMust be subset of: ", paste(valid_distributions, collapse = ", "))
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
  arm_levels <- unique(dataset$ARM)
  if (length(arm_levels) == 1 && TRUE %in% dependent) {
    message("ARM has only one level. Removing dependent = TRUE from analysis.")
    dependent <- FALSE
  }

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
        dataset = dataset,
        distribution = dist,
        dependent = dep,
        control_arm = control_arm,
        landmark_increment = landmark_increment,
        time_horizon = time_horizon
      )
    }, error = function(e) {
      stop("FITTING ERROR - distribution: ", dist, ", dependent: ", dep,
           "\nError message: ", e$message, call.=F)
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
      all_landmark[[i]] <- fit_result$landmark_survival
    }

    # Store survival metrics
    metrics <- add_model_arm_columns(fit_result$survival_metrics, model_type)
    all_metrics[[i]] <- metrics

    message("Successfully fitted: ", dist, " (", model_type, ")")
  }

  # Combine results
  coef_estimates_combined <- dplyr::bind_rows(all_coef_estimates)
  varcovmat_combined <- dplyr::bind_rows(all_varcovmat)
  cholesky_combined <- dplyr::bind_rows(all_cholesky)
  gof_stats_combined <- dplyr::bind_rows(all_gof_stats)
  metrics_combined <- dplyr::bind_rows(all_metrics)

  # Combine landmark survival results
  landmark_survival_combined <- NULL
  if (length(all_landmark) > 0) {
    landmark_survival_combined <- dplyr::bind_rows(all_landmark) %>%
      dplyr::group_by(model_type, ARM, Time) %>%
      dplyr::summarise(
        across(everything(), ~dplyr::first(na.omit(.))),
        .groups = 'drop'
      ) %>%
      dplyr::arrange(model_type, ARM, Time)
  }

  # Fit Cox PH model
  cox_fit <- survival::coxph(
    survival::Surv(SURVTIME, EVENT) ~ ARM,
    data = dataset
  )
  cox_summary <- summary(cox_fit)
  cox_ci <- cox_summary$conf.int
  cox_results <- data.frame(
    term = rownames(cox_ci),
    HR = round(cox_ci[, 1], 3),
    check.names = FALSE
  )
  cox_results[["L95%"]] <- round(cox_ci[, 3], 3)
  cox_results[["U95%"]] <- round(cox_ci[, 4], 3)
  cox_results[["p.value"]] <- round(cox_summary$coefficients[, "Pr(>|z|)"], 4)
  cox_results <- cox_results %>%
    dplyr::select(term, HR, `L95%`, `U95%`, `p.value`) %>%
    dplyr::mutate(
      event_count = sum(dataset$EVENT, na.rm = TRUE),
      sample_size = nrow(dataset)
    )

  # Return results as list
  list(
    coef_estimates = coef_estimates_combined,
    varcovmat = varcovmat_combined,
    cholesky = cholesky_combined,
    gof_statistics = gof_stats_combined,
    landmark_survival = landmark_survival_combined,
    survival_metrics = metrics_combined,
    cox_ph_results = cox_results
  )
}
