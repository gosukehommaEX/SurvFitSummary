#' Plot Smoothed Empirical Hazard with Multiple Parametric Distributions
#'
#' This function creates multiple visualizations of smoothed empirical hazard functions,
#' each overlaid with a different parametric distribution fit. Each plot retains its own
#' independent y-axis scale (not normalized across distributions).
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distributions Character vector specifying parametric distributions to plot.
#'   Each element must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz",
#'   "gengamma", "gamma". Example: c("weibull", "gompertz", "lnorm")
#' @param conf_int Logical indicating whether to display 95% confidence intervals on
#'   hazard curves. Default is TRUE. Applied to all distributions
#'
#' @return A list containing:
#'   \describe{
#'     \item{plots}{Named list of ggplot objects, one for each distribution.
#'       Names correspond to distribution names (e.g., "weibull", "gompertz").
#'       Each plot has its own independent y-axis scale}
#'   }
#'
#' @details
#' The function calls \code{plot_hazard_and_parametric()} for each distribution.
#' Unlike the previous version, each distribution plot retains its own y-axis scale,
#' which is particularly important when distributions show very different hazard patterns
#' (e.g., exponential vs generalized gamma).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Compare multiple distributions with independent y-axes
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
#' result <- plot_hazard_and_parametric_multiple(
#'   dataset = dataset_processed,
#'   distributions = c("weibull", "gompertz", "lnorm"),
#'   conf_int = TRUE
#' )
#'
#' # Access individual plots
#' print(result$plots$weibull)
#' print(result$plots$gompertz)
#' print(result$plots$lnorm)
#' }

plot_hazard_and_parametric_multiple <- function(dataset,
                                                distributions,
                                                conf_int = TRUE) {

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

  # Initialize list to store plots
  message("Creating hazard plots for all distributions with independent y-axis scaling...")
  plots_list <- list()

  # Execute plot_hazard_and_parametric for each distribution
  for (i in seq_along(distributions)) {
    dist <- distributions[i]

    tryCatch({
      result <- plot_hazard_and_parametric(
        dataset = dataset,
        distribution = dist,
        conf_int = conf_int
      )

      # Store plot with distribution name as key
      plots_list[[dist]] <- result$plot

      message("  Successfully created plot for: ", dist)

    }, error = function(e) {
      warning("Failed to process distribution '", dist, "': ", e$message)
    })
  }

  # Check if any plots were created
  if (length(plots_list) == 0) {
    stop("No hazard plots could be created from any distribution")
  }

  message("All hazard plots completed successfully!")

  # Return list with plots (each with independent y-axis)
  return(list(
    plots = plots_list
  ))
}
