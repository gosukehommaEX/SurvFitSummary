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
#' @param landmark_times Numeric vector specifying time points for landmark survival
#'   probability estimation. If NULL (default), landmark survival is not calculated
#' @param time_horizon Numeric value specifying the time horizon for restricted mean
#'   survival time (RMST) calculation. If NULL (default), RMST is not calculated
#' @param original_time_scale Numeric value for converting results back to original time scale.
#'   If NULL (default), results are in the same units as input data. When provided,
#'   survival metrics are divided by this value to convert back to original scale
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
#'     \item{landmark_survival}{Data frame with landmark survival probabilities (NULL if not specified)}
#'     \item{survival_metrics}{Data frame with median, mean, and RMST}
#'   }
#'
#' @importFrom dplyr as_tibble select mutate bind_cols group_by reframe filter rename tibble
#' @importFrom purrr map_dfr
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv
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
#'   landmark_times = c(6, 12, 18, 24),
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
#'   landmark_times = c(6, 12, 18, 24),
#'   time_horizon = 36
#' )
#'
#' print(result_indep$coef_estimates)
#' print(result_indep$gof_statistics)
#'
#' # Example 3: With time scale conversion
#' # Data in weeks, convert survival metrics back to months
#' result_scaled <- fitting_surv_mod(
#'   dataset = dataset_processed,
#'   distribution = "exp",
#'   dependent = TRUE,
#'   landmark_times = seq(4, 52, by = 4),
#'   time_horizon = 52,
#'   original_time_scale = 30.44 / 7  # Convert weeks to months
#' )
#'
#' print(result_scaled$survival_metrics)
#' }
fitting_surv_mod <- function(dataset,
                             distribution,
                             dependent = TRUE,
                             control_arm = NULL,
                             landmark_times = NULL,
                             time_horizon = NULL,
                             original_time_scale = NULL) {

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
    {if(dependent) . else dplyr::group_by(., ARM)} %>%
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

  # Landmark survival probabilities
  landmark_survival <- NULL
  if(!is.null(landmark_times)) {
    landmark_survival <- if(dependent) {
      fit_result <- model_fit$fit_result[[1]]
      purrr::map_dfr(arm_levels_subset, function(arm_val) {
        newdata <- data.frame(ARM = arm_val)
        summ <- summary(fit_result, newdata = newdata, t = landmark_times, type = 'survival')
        dplyr::tibble(ARM = arm_val, summ[[1]])
      })
    } else {
      purrr::map_dfr(1:nrow(model_fit), function(i) {
        summ <- summary(model_fit$fit_result[[i]], t = landmark_times, type = 'survival')[[1]]
        dplyr::tibble(ARM = model_fit$ARM[i], summ)
      })
    }
    landmark_survival <- landmark_survival %>%
      dplyr::rename(
        'Survival Probability' = est,
        'L95%' = lcl,
        'U95%' = ucl
      ) %>%
      dplyr::mutate(distribution = distribution, .before = ARM)
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

  # Convert survival metrics back to original time scale if specified
  if (!is.null(original_time_scale)) {
    survival_metrics <- survival_metrics %>%
      dplyr::mutate(
        `Median Survival Time` = `Median Survival Time` / original_time_scale,
        `Mean Survival Time` = `Mean Survival Time` / original_time_scale,
        `Restricted Mean Survival Time` = `Restricted Mean Survival Time` / original_time_scale
      )
  }

  # Return results
  return(list(
    coef_estimates = coef_estimates,
    varcovmat = varcovmat,
    cholesky = cholesky,
    gof_statistics = gof_stats,
    landmark_survival = landmark_survival,
    survival_metrics = survival_metrics
  ))
}
