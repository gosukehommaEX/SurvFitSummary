#' Fit Parametric Survival Models with Multiple Distributions
#'
#' This function fits parametric survival models to time-to-event data using
#' various distributions. It supports both dependent (jointly fitted with ARM
#' as covariate) and independent models.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distribution Character string specifying the parametric distribution.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param dependent Logical indicating whether to fit dependent model (ARM as covariate).
#'   If TRUE, fits model with ~ ARM formula. If FALSE, fits separate models for each ARM
#'   with ~ 1 formula. Default is TRUE
#' @param control_arm Character string specifying which ARM level should be used as
#'   the reference (control) group. If NULL (default), uses the first factor level
#' @param landmark_increment Numeric value specifying the time interval for landmark survival
#'   probability estimation. If NULL (default), landmark survival is not calculated.
#'   Results will be calculated at times 0, landmark_increment, 2*landmark_increment, etc.
#' @param time_horizon Numeric value specifying the time horizon for restricted mean
#'   survival time (RMST) calculation. If NULL (default), RMST is not calculated
#'
#' @return A list containing six components:
#'   \describe{
#'     \item{coef_estimates}{Data frame with parameter estimates, standard errors,
#'       confidence intervals, and p-values. For independent models, statistic and
#'       p.value columns are excluded}
#'     \item{varcovmat}{Data frame in long format with columns: distribution, ARM,
#'       row_param, col_param, value}
#'     \item{cholesky}{Data frame in long format with columns: distribution, ARM,
#'       row_param, col_param, value}
#'     \item{gof_statistics}{Data frame with goodness-of-fit statistics (AIC, BIC)}
#'     \item{landmark_survival}{Data frame with columns: model_type, ARM, Time, KM (Kaplan-Meier),
#'       and distribution name, at intervals specified by landmark_increment (NULL if not specified)}
#'     \item{survival_metrics}{Data frame with median, mean, and RMST}
#'   }
#'
#' @importFrom dplyr as_tibble select mutate bind_cols group_by reframe filter rename tibble
#' @importFrom purrr map_dfr
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv survfit
#' @importFrom stats vcov AIC BIC
#' @importFrom broom tidy
#' @importFrom tidyr pivot_longer expand_grid
#' @importFrom tibble rownames_to_column
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Dependent model (ARM as covariate)
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
#' result_dep <- fitting_surv_mod(
#'   dataset = dataset_processed,
#'   distribution = "weibull",
#'   dependent = TRUE,
#'   control_arm = "Control",
#'   landmark_increment = 6,
#'   time_horizon = 36
#' )
#'
#' print(result_dep$coef_estimates)
#' print(result_dep$landmark_survival)
#' print(result_dep$survival_metrics)
#'
#' # Example 2: Independent models (separate fits per ARM)
#' result_indep <- fitting_surv_mod(
#'   dataset = dataset_processed,
#'   distribution = "weibull",
#'   dependent = FALSE,
#'   landmark_increment = 6,
#'   time_horizon = 36
#' )
#'
#' print(result_indep$coef_estimates)
#' print(result_indep$gof_statistics)
#'
#' # Example 3: With different increment
#' result_scaled <- fitting_surv_mod(
#'   dataset = dataset_processed,
#'   distribution = "exp",
#'   dependent = TRUE,
#'   landmark_increment = 4,
#'   time_horizon = 52
#' )
#'
#' print(result_scaled$survival_metrics)
#' }
fitting_surv_mod <- function(dataset,
                             distribution,
                             dependent = TRUE,
                             control_arm = NULL,
                             landmark_increment = NULL,
                             time_horizon = NULL) {

  # Validate distribution parameter
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!distribution %in% valid_distributions) {
    stop("Distribution must be one of: ", paste(valid_distributions, collapse = ", "))
  }

  # Check ARM levels and adjust dependent parameter if needed
  arm_levels <- unique(dataset$ARM)
  if (length(arm_levels) == 1) {
    if (dependent) {
      message("ARM has only one level. Setting dependent = FALSE automatically.")
      dependent <- FALSE
    }
  }

  # Handle control_arm specification
  if (!is.null(control_arm)) {
    if (!control_arm %in% arm_levels) {
      stop("Specified control_arm '", control_arm, "' not found in ARM levels: ",
           paste(arm_levels, collapse = ", "))
    }
    dataset$ARM <- relevel(factor(dataset$ARM), ref = control_arm)
    arm_levels <- levels(dataset$ARM)
    message("ARM releveled with '", control_arm, "' as reference group")
  } else {
    if (!is.factor(dataset$ARM)) {
      dataset$ARM <- as.factor(dataset$ARM)
    }
  }

  # Set seed for reproducibility
  set.seed(1)

  # Model fitting with conditional formula
  model_fit <- dataset %>%
    {if(dependent) .
      else dplyr::group_by(., ARM)} %>%
    dplyr::reframe(
      distribution = distribution,
      fit_result = list(
        if(dependent) {
          flexsurv::flexsurvreg(
            survival::Surv(SURVTIME, EVENT) ~ ARM,
            data = dplyr::pick(everything()),
            dist = distribution
          )
        } else {
          flexsurv::flexsurvreg(
            survival::Surv(SURVTIME, EVENT) ~ 1,
            data = dplyr::pick(everything()),
            dist = distribution
          )
        }
      )
    )

  # Extract fit result and ARM levels
  arm_levels_subset <- unique(dataset[['ARM']])

  # Coefficient estimates
  coef_estimates <- if(dependent) {
    fit_result <- model_fit$fit_result[[1]]

    broom::tidy(fit_result) %>%
      dplyr::bind_cols(
        fit_result[['res']] %>%
          dplyr::as_tibble() %>%
          dplyr::select('L95%', 'U95%')
      ) %>%
      dplyr::select(term, estimate, std.error, 'L95%', 'U95%', statistic, p.value) %>%
      dplyr::mutate(distribution = distribution, .before = 1)
  } else {
    purrr::map_dfr(1:nrow(model_fit), function(i) {
      broom::tidy(model_fit$fit_result[[i]]) %>%
        dplyr::bind_cols(
          model_fit$fit_result[[i]][['res']] %>%
            dplyr::as_tibble() %>%
            dplyr::select('L95%', 'U95%')
        ) %>%
        dplyr::select(term, estimate, std.error, 'L95%', 'U95%') %>%
        dplyr::mutate(
          distribution = distribution,
          ARM = model_fit$ARM[i],
          .before = 1
        )
    })
  }

  # Helper function to convert variance-covariance matrix to long format
  vcov_to_long <- function(vcov_mat, distribution, arm = NULL, is_dependent = TRUE) {
    # Special handling for exponential distribution with independent models
    # (single parameter "rate" case)
    if (distribution == "exp" && !is_dependent) {
      result <- dplyr::tibble(
        row_param = "rate",
        col_param = "rate",
        value = vcov_mat[1, 1]
      )
    } else {
      result <- vcov_mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column("row_param") %>%
        tidyr::pivot_longer(
          cols = -row_param,
          names_to = "col_param",
          values_to = "value"
        )
    }

    result <- result %>%
      dplyr::mutate(distribution = distribution, .before = 1)

    if (!is.null(arm)) {
      result <- result %>%
        dplyr::mutate(ARM = arm, .after = distribution)
    }

    return(result)
  }

  # Variance-covariance matrix (long format)
  varcovmat <- if(dependent) {
    fit_result <- model_fit$fit_result[[1]]
    vcov_mat <- stats::vcov(fit_result)
    vcov_to_long(vcov_mat, distribution, is_dependent = TRUE)
  } else {
    purrr::map_dfr(1:nrow(model_fit), function(i) {
      vcov_mat <- stats::vcov(model_fit$fit_result[[i]])
      vcov_to_long(vcov_mat, distribution, arm = model_fit$ARM[i], is_dependent = FALSE)
    })
  }

  # Cholesky decomposition (long format)
  cholesky <- if(dependent) {
    fit_result <- model_fit$fit_result[[1]]
    vcov_mat <- stats::vcov(fit_result)
    chol_mat <- chol(vcov_mat)
    vcov_to_long(chol_mat, distribution, is_dependent = TRUE)
  } else {
    purrr::map_dfr(1:nrow(model_fit), function(i) {
      vcov_mat <- stats::vcov(model_fit$fit_result[[i]])
      chol_mat <- chol(vcov_mat)
      vcov_to_long(chol_mat, distribution, arm = model_fit$ARM[i], is_dependent = FALSE)
    })
  }

  # Goodness of fit statistics
  gof_stats <- if(dependent) {
    fit_result <- model_fit$fit_result[[1]]

    dplyr::tibble(
      distribution = distribution,
      AIC = stats::AIC(fit_result),
      BIC = stats::BIC(fit_result)
    )
  } else {
    purrr::map_dfr(1:nrow(model_fit), function(i) {
      dplyr::tibble(
        distribution = distribution,
        ARM = model_fit$ARM[i],
        AIC = stats::AIC(model_fit$fit_result[[i]]),
        BIC = stats::BIC(model_fit$fit_result[[i]])
      )
    })
  }

  # Landmark survival probabilities (ARM-specific)
  landmark_survival <- NULL
  if(!is.null(landmark_increment)) {
    # Calculate maximum time from dataset
    max_time <- max(dataset$SURVTIME, na.rm = TRUE)

    # Generate time sequence from 0 with specified increment
    landmark_times <- seq(0, max_time, by = landmark_increment)

    # Calculate ARM-specific landmark survival
    landmark_results <- purrr::map_dfr(arm_levels_subset, function(arm_val) {
      # Fit Kaplan-Meier for this ARM
      km_fit <- survival::survfit(
        survival::Surv(SURVTIME, EVENT) ~ 1,
        data = dataset %>% dplyr::filter(ARM == arm_val)
      )

      # Get K-M survival probabilities at landmark times
      km_surv <- summary(km_fit, times = landmark_times, extend = TRUE)
      km_probs <- km_surv$surv

      # Get parametric model predictions for this ARM
      if(dependent) {
        fit_result <- model_fit$fit_result[[1]]
        newdata <- data.frame(ARM = arm_val)
        param_summ <- summary(fit_result, newdata = newdata, t = landmark_times, type = 'survival')
        param_probs <- param_summ[[1]]$est
      } else {
        # Find the index for this ARM in model_fit
        arm_idx <- which(model_fit$ARM == arm_val)
        summ <- summary(model_fit$fit_result[[arm_idx]], t = landmark_times, type = 'survival')
        param_probs <- summ[[1]]$est
      }

      # Create output data frame for this ARM
      result_df <- dplyr::tibble(
        model_type = ifelse(dependent, "dependent", "independent"),
        ARM = arm_val,
        Time = landmark_times,
        KM = km_probs
      )
      result_df[[distribution]] <- param_probs
      result_df
    })

    landmark_survival <- landmark_results
  }

  # Survival metrics
  survival_metrics <- if(dependent) {
    fit_result <- model_fit$fit_result[[1]]
    purrr::map_dfr(arm_levels_subset, function(arm_val) {
      newdata <- data.frame(ARM = arm_val)
      dplyr::tibble(
        distribution = distribution,
        ARM = arm_val,
        'Median Survival Time' = tryCatch(
          summary(fit_result, newdata = newdata, type = 'median')[[1]][['est']][1],
          error = function(e) NA_real_
        ),
        'Mean Survival Time' = tryCatch(
          summary(fit_result, newdata = newdata, type = 'mean')[[1]][['est']][1],
          error = function(e) NA_real_
        ),
        'Restricted Mean Survival Time' = if(!is.null(time_horizon)) {
          tryCatch(
            summary(fit_result, newdata = newdata, t = time_horizon, type = 'rmst')[[1]][['est']][1],
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
      )
    })
  } else {
    purrr::map_dfr(1:nrow(model_fit), function(i) {
      dplyr::tibble(
        distribution = distribution,
        ARM = model_fit$ARM[i],
        'Median Survival Time' = tryCatch(
          summary(model_fit$fit_result[[i]], type = 'median')[[1]][['est']][1],
          error = function(e) NA_real_
        ),
        'Mean Survival Time' = tryCatch(
          summary(model_fit$fit_result[[i]], type = 'mean')[[1]][['est']][1],
          error = function(e) NA_real_
        ),
        'Restricted Mean Survival Time' = if(!is.null(time_horizon)) {
          tryCatch(
            summary(model_fit$fit_result[[i]], t = time_horizon, type = 'rmst')[[1]][['est']][1],
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
      )
    })
  }

  # Return results as list
  list(
    coef_estimates = coef_estimates,
    varcovmat = varcovmat,
    cholesky = cholesky,
    gof_statistics = gof_stats,
    landmark_survival = landmark_survival,
    survival_metrics = survival_metrics
  )
}
