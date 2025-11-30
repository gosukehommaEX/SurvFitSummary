#' Plot Smoothed Empirical Hazard with Parametric Distribution Overlay (Facet Grid Version)
#'
#' This function creates visualizations of smoothed empirical hazard functions
#' overlaid with parametric distribution hazard curves using facet_grid.
#' Rows represent treatment arms, columns represent time scales (Time vs log(Time)).
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distribution Character string specifying the parametric distribution.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param conf_int Logical indicating whether to display 95% confidence intervals on
#'   hazard curves. Default is TRUE
#'
#' @return A list containing:
#'   \describe{
#'     \item{plot}{A ggplot object with facet_grid layout. Rows = ARMs, Columns = Time scales}
#'     \item{hazard_data}{A named list of data frames, one for each treatment arm.
#'       Each data frame contains columns: time, hazard, lower_ci, upper_ci, type, ARM}
#'   }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Calculates smoothed empirical hazard using \code{bshazard::bshazard()}
#'   \item Fits parametric survival models independently for each ARM
#'   \item Parametric hazard is calculated over the same time range as empirical hazard
#'   \item Creates facet_grid with ARM as rows and time scale as columns
#'   \item Each ARM has independent x-axis range based on its data
#' }
#'
#' @importFrom survival Surv
#' @importFrom bshazard bshazard
#' @importFrom flexsurv flexsurvreg
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_color_manual scale_fill_manual
#'   scale_x_continuous labs theme_bw theme element_text element_blank
#'   element_rect element_line guides guide_legend facet_grid vars
#' @importFrom dplyr filter mutate bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Single-arm trial
#' dataset1 <- generate_dummy_survival_data(
#'   arm = c("Treatment"),
#'   n = c(150),
#'   hazards = log(2) / 20,
#'   dropout_per_year = 0.1,
#'   seed = 123
#' )
#'
#' dataset1_processed <- processing_dataset(
#'   dataset = dataset1,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = NULL
#' )
#'
#' p1 <- plot_hazard_and_parametric(
#'   dataset = dataset1_processed,
#'   distribution = "weibull",
#'   conf_int = TRUE
#' )
#' print(p1$plot)
#'
#' # Example 2: Two-arm trial
#' dataset2 <- generate_dummy_survival_data(
#'   arm = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazards = log(2) / c(20, 15),
#'   dropout_per_year = 0.1,
#'   seed = 456
#' )
#'
#' dataset2_processed <- processing_dataset(
#'   dataset = dataset2,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = NULL
#' )
#'
#' p2 <- plot_hazard_and_parametric(
#'   dataset = dataset2_processed,
#'   distribution = "exp",
#'   conf_int = TRUE
#' )
#' print(p2$plot)
#'
#' # Example 3: Three-arm trial
#' dataset3 <- generate_dummy_survival_data(
#'   arm = c("Treatment A", "Treatment B", "Control"),
#'   n = c(80, 80, 80),
#'   hazards = log(2) / c(22, 18, 15),
#'   dropout_per_year = 0.1,
#'   seed = 789
#' )
#'
#' dataset3_processed <- processing_dataset(
#'   dataset = dataset3,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = NULL
#' )
#'
#' p3 <- plot_hazard_and_parametric(
#'   dataset = dataset3_processed,
#'   distribution = "gompertz",
#'   conf_int = TRUE
#' )
#' print(p3$plot)
#' }

plot_hazard_and_parametric <- function(dataset,
                                       distribution,
                                       conf_int = TRUE) {

  # Load required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!requireNamespace("bshazard", quietly = TRUE)) {
    stop("Package 'bshazard' is required but not installed.")
  }
  if (!requireNamespace("flexsurv", quietly = TRUE)) {
    stop("Package 'flexsurv' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }

  # Validate distribution parameter
  valid_distributions <- c("exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma")
  if (!distribution %in% valid_distributions) {
    stop("Distribution must be one of: ", paste(valid_distributions, collapse = ", "))
  }

  # Get ARM levels
  arm_levels <- unique(dataset$ARM)
  if (!is.factor(dataset$ARM)) {
    dataset$ARM <- as.factor(dataset$ARM)
    arm_levels <- levels(dataset$ARM)
  }
  n_arms <- length(arm_levels)

  # Distribution names
  dist_names <- c(
    "exp" = "Exponential",
    "weibull" = "Weibull",
    "lnorm" = "Log-Normal",
    "llogis" = "Log-Logistic",
    "gompertz" = "Gompertz",
    "gengamma" = "Generalized Gamma",
    "gamma" = "Gamma"
  )

  # Colors
  color_smoothed <- '#004C97'
  color_parametric <- '#F0B323'

  # Collect all data for faceting
  all_hazard_data <- list()
  hazard_data_list <- list()

  for (i in seq_along(arm_levels)) {
    arm_name <- arm_levels[i]
    arm_data <- dataset[dataset$ARM == arm_name, ]

    # Calculate smoothed empirical hazard
    tryCatch({
      bsh_fit <- bshazard::bshazard(
        survival::Surv(SURVTIME, EVENT) ~ 1,
        data = arm_data,
        verbose = FALSE
      )

      hazard_empirical <- data.frame(
        time = bsh_fit$time,
        hazard = bsh_fit$hazard,
        lower_ci = bsh_fit$lower.ci,
        upper_ci = bsh_fit$upper.ci,
        type = 'Smoothed Empirical',
        ARM = arm_name
      ) %>%
        dplyr::filter(time > 0, hazard > 0, is.finite(hazard))

      # Fit parametric distribution
      param_fit <- flexsurv::flexsurvreg(
        survival::Surv(SURVTIME, EVENT) ~ 1,
        data = arm_data,
        dist = distribution
      )

      # Get parametric hazard - match the range of smoothed empirical
      max_empirical_time <- max(hazard_empirical$time, na.rm = TRUE)
      time_seq <- seq(0.01, max_empirical_time, length.out = 200)
      param_pred <- summary(param_fit, t = time_seq, type = 'hazard')[[1]]

      hazard_parametric <- data.frame(
        time = time_seq,
        hazard = param_pred$est,
        lower_ci = param_pred$lcl,
        upper_ci = param_pred$ucl,
        type = 'Parametric',
        ARM = arm_name
      ) %>%
        dplyr::filter(time > 0, hazard > 0, is.finite(hazard))

      # Combine data for this ARM
      hazard_combined <- rbind(hazard_empirical, hazard_parametric)

      # Store hazard data for external use
      hazard_data_list[[as.character(arm_name)]] <- hazard_combined

      # Add to list for faceting
      all_hazard_data[[i]] <- hazard_combined

    }, error = function(e) {
      warning(paste('Failed to process ARM', arm_name, ':', e$message))
    })
  }

  # Check if any data was created
  if (length(all_hazard_data) == 0) {
    stop("No valid hazard data could be created for any ARM")
  }

  # Combine all ARM data
  combined_data <- dplyr::bind_rows(all_hazard_data)

  # Create data for Time scale
  data_time <- combined_data %>%
    dplyr::mutate(
      time_value = time,
      time_scale_label = "Time"
    )

  # Create data for log(Time) scale
  data_logtime <- combined_data %>%
    dplyr::mutate(
      time_value = log(time),
      time_scale_label = "log(Time)"
    )

  # Combine both scales
  plot_data <- dplyr::bind_rows(data_time, data_logtime)

  # Ensure ARM and time_scale_label are factors with correct order
  plot_data$ARM <- factor(plot_data$ARM, levels = arm_levels)
  plot_data$time_scale_label <- factor(plot_data$time_scale_label, levels = c("Time", "log(Time)"))

  # Create faceted plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time_value, y = hazard, color = type, linetype = type)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(
      values = c('Smoothed Empirical' = color_smoothed, 'Parametric' = color_parametric),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c('Smoothed Empirical' = 'solid', 'Parametric' = 'dashed'),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(ARM),
      cols = ggplot2::vars(time_scale_label),
      scales = "free_x"
    ) +
    ggplot2::labs(
      title = 'Smoothed Empirical Hazard with Parametric Distribution Overlay',
      subtitle = paste('Distribution:', dist_names[distribution]),
      x = 'Time',
      y = 'Hazard h(Time)',
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, face = 'bold', hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 16, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 16, face = 'bold'),
      strip.background = ggplot2::element_rect(fill = '#E8E8E8', color = 'black'),
      legend.position = 'bottom',
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.key.width = ggplot2::unit(2, "cm"),
      panel.grid.major = ggplot2::element_line(color = 'gray90'),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add confidence intervals if requested
  if (conf_int) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_ci, ymax = upper_ci, fill = type),
        alpha = 0.2,
        color = NA
      ) +
      ggplot2::scale_fill_manual(
        values = c('Smoothed Empirical' = color_smoothed, 'Parametric' = color_parametric),
        name = NULL
      )
  }

  # Return list with plot and hazard data
  return(list(
    plot = p,
    hazard_data = hazard_data_list
  ))
}
