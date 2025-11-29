#' Plot Kaplan-Meier Curves with Multiple Parametric Distributions
#'
#' This function creates multiple visualizations of Kaplan-Meier survival curves,
#' each overlaid with a different parametric distribution fit. It wraps
#' \code{plot_km_and_parametric()} and executes it for multiple distributions.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distributions Character vector specifying parametric distributions to plot.
#'   Each element must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz",
#'   "gengamma", "gamma". Example: c("weibull", "gompertz", "lnorm")
#' @param conf_int Logical indicating whether to display 95% confidence intervals on
#'   Kaplan-Meier curves. Default is TRUE. Applied to all distributions
#' @param time_scale Character string specifying the unit of time scale in the input dataset.
#'   Must be one of: "day", "week", "month", "year". Default is "week"
#' @param time_horizon Numeric value specifying the maximum time point in the specified
#'   time_scale unit. If NULL (default), uses the maximum observed survival time in the
#'   dataset
#'
#' @return A list containing:
#'   \describe{
#'     \item{plots}{Named list of ggplot objects, one for each distribution.
#'       Names correspond to distribution names (e.g., "weibull", "gompertz")}
#'     \item{risktable}{Data frame showing number at risk at 4-week intervals.
#'       This is shared across all distributions as it depends only on the data,
#'       not the fitted distribution}
#'   }
#'
#' @details
#' The function iterates through each distribution and calls \code{plot_km_and_parametric()}
#' internally. The risk table is calculated once from the raw data and is identical
#' for all distributions, as it only depends on the observed survival times, not the
#' fitted parametric model.
#'
#' @importFrom dplyr filter select mutate arrange
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Compare multiple distributions
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
#' result <- plot_km_and_parametric_multiple(
#'   dataset = dataset_processed,
#'   distributions = c("weibull", "gompertz", "lnorm"),
#'   conf_int = TRUE,
#'   time_scale = "week",
#'   time_horizon = 100
#' )
#'
#' # Access individual plots
#' print(result$plots$weibull)
#' print(result$plots$gompertz)
#' print(result$plots$lnorm)
#'
#' # Access shared risk table
#' print(result$risktable)
#' }

plot_km_and_parametric_multiple <- function(dataset,
                                            distributions,
                                            conf_int = TRUE,
                                            time_scale = "week",
                                            time_horizon = NULL) {

  # Validate distributions parameter
  if (!is.character(distributions)) {
    stop("distributions must be a character vector")
  }

  if (length(distributions) == 0) {
    stop("distributions must contain at least one distribution")
  }

  # Valid distributions
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!all(distributions %in% valid_distributions)) {
    invalid <- distributions[!distributions %in% valid_distributions]
    stop("Invalid distribution(s): ", paste(invalid, collapse = ", "))
  }

  # Initialize list to store plots and risktable
  plots_list <- list()
  risktable <- NULL

  # Execute plot_km_and_parametric for each distribution
  for (i in seq_along(distributions)) {
    dist <- distributions[i]

    # Call plot_km_and_parametric
    result <- plot_km_and_parametric(
      dataset = dataset,
      distribution = dist,
      conf_int = conf_int,
      time_scale = time_scale,
      time_horizon = time_horizon
    )

    # Store plot with distribution name as key
    plots_list[[dist]] <- result$plot

    # Store risktable from first iteration only (shared across all distributions)
    if (is.null(risktable)) {
      risktable <- result$risktable
    }
  }

  # Return list with plots and shared risktable
  return(list(
    plots = plots_list,
    risktable = risktable
  ))
}
