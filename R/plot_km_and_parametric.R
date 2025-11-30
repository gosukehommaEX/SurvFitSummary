#' Plot Kaplan-Meier Curves with Parametric Distribution Overlay
#'
#' This function creates a combined visualization of Kaplan-Meier survival curves
#' with overlaid parametric distribution fits. It supports multiple treatment arms
#' and displays confidence intervals optionally.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distribution Character string specifying the parametric distribution.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param conf_int Logical indicating whether to display 95% confidence intervals on
#'   Kaplan-Meier curves. Default is TRUE
#'
#' @return A list containing:
#'   \describe{
#'     \item{plot}{ggplot object with Kaplan-Meier curves and parametric overlays}
#'     \item{risktable}{Data frame showing number at risk at regular intervals}
#'   }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Fits Kaplan-Meier curves using \code{survival::survfit()}
#'   \item Fits parametric survival models independently for each ARM
#'   \item Calculates number at risk at regular intervals
#'   \item Axes begin at (0, 0) for improved visualization
#'   \item X-axis extends to the longest ARM's last event time
#'   \item Handles single and multiple ARM datasets
#' }
#'
#' @importFrom survival Surv survfit
#' @importFrom flexsurv flexsurvreg
#' @importFrom ggplot2 ggplot aes geom_line geom_step geom_ribbon scale_color_manual
#'   scale_fill_manual scale_x_continuous scale_y_continuous labs theme element_text
#'   element_blank element_line annotate expansion
#' @importFrom dplyr group_by summarise n mutate select filter arrange
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_wider
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Single arm analysis
#' dataset1 <- generate_dummy_survival_data(
#'   arm = c("Treatment"),
#'   n = c(100),
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
#' result <- plot_km_and_parametric(
#'   dataset = dataset1_processed,
#'   distribution = "weibull",
#'   conf_int = TRUE
#' )
#'
#' print(result$plot)
#' print(result$risktable)
#'
#' # Example 2: Two-arm trial
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
#' result <- plot_km_and_parametric(
#'   dataset = dataset_processed,
#'   distribution = "weibull",
#'   conf_int = TRUE
#' )
#'
#' print(result$plot)
#' print(result$risktable)
#'
#' # Example 3: Three-arm trial
#' dataset3 <- generate_dummy_survival_data(
#'   arm = c("Treatment A", "Treatment B", "Control"),
#'   n = c(80, 90, 100),
#'   hazards = log(2) / c(25, 20, 15),
#'   dropout_per_year = 0.1,
#'   seed = 456
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
#' result3 <- plot_km_and_parametric(
#'   dataset = dataset3_processed,
#'   distribution = "weibull",
#'   conf_int = TRUE
#' )
#'
#' print(result3$plot)
#' print(result3$risktable)
#' }

plot_km_and_parametric <- function(dataset,
                                   distribution,
                                   conf_int = TRUE) {

  # Load required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
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

  # Define color palette
  color_palette <- c('#004C97', '#F0B323', '#658D1B', '#A62B4E', '#D91E49', '#939597')

  # Get ARM levels and assign colors
  arm_levels <- unique(dataset$ARM)
  n_arms <- length(arm_levels)

  # Assign colors based on number of ARMs
  colors <- color_palette[1:min(n_arms, length(color_palette))]
  if (n_arms > length(color_palette)) {
    colors <- rep(color_palette, length.out = n_arms)
  }
  names(colors) <- as.character(arm_levels)

  # Distribution name mapping
  dist_names <- c(
    'exp' = 'Exponential',
    'weibull' = 'Weibull',
    'lnorm' = 'Lognormal',
    'llogis' = 'Loglogistic',
    'gompertz' = 'Gompertz',
    'gengamma' = 'Generalized Gamma',
    'gamma' = 'Gamma'
  )

  # Fit Kaplan-Meier model
  km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ ARM, data = dataset)

  # Extract KM curve data (handle both single and multiple ARMs)
  if (n_arms == 1) {
    # Single ARM case
    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      ARM = arm_levels[1],
      lower = km_fit$lower,
      upper = km_fit$upper,
      stringsAsFactors = FALSE
    )
  } else {
    # Multiple ARMs case
    # Extract ARM names from strata names (remove "ARM=" prefix)
    strata_names <- names(km_fit$strata)
    arm_names <- gsub("ARM=", "", strata_names)

    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      ARM = rep(arm_names, km_fit$strata),
      lower = km_fit$lower,
      upper = km_fit$upper,
      stringsAsFactors = FALSE
    )
  }

  # Add point at time = 0 with survival = 1 for each ARM
  km_data_start <- data.frame(
    time = 0,
    surv = 1,
    ARM = arm_levels,
    lower = 1,
    upper = 1,
    stringsAsFactors = FALSE
  )

  km_data <- rbind(km_data_start, km_data)

  # Sort by ARM and time
  km_data <- km_data %>%
    dplyr::arrange(ARM, time)

  # Determine the plot time horizon (longest ARM's last event time)
  max_time_per_arm <- km_data %>%
    dplyr::group_by(ARM) %>%
    dplyr::summarise(max_time = max(time), .groups = 'drop')

  plot_time_horizon <- max(max_time_per_arm$max_time)

  # Extend KM data to plot_time_horizon for each ARM
  km_data_extended <- data.frame()

  for (arm in arm_levels) {
    km_data_arm <- km_data %>% dplyr::filter(ARM == arm)

    # Get last observation for this ARM
    last_obs <- km_data_arm[nrow(km_data_arm), ]

    # If last time < plot_time_horizon, add endpoint
    if (last_obs$time < plot_time_horizon) {
      end_point <- data.frame(
        time = plot_time_horizon,
        surv = last_obs$surv,
        ARM = arm,
        lower = last_obs$lower,
        upper = last_obs$upper,
        stringsAsFactors = FALSE
      )
      km_data_arm <- rbind(km_data_arm, end_point)
    }

    km_data_extended <- rbind(km_data_extended, km_data_arm)
  }

  km_data <- km_data_extended

  # Fit parametric survival models for each ARM
  parametric_data <- data.frame()

  for (arm in arm_levels) {
    # Subset data for current ARM
    data_arm <- dataset %>% dplyr::filter(ARM == arm)

    # Fit parametric model
    fit <- flexsurv::flexsurvreg(
      survival::Surv(SURVTIME, EVENT) ~ 1,
      data = data_arm,
      dist = distribution
    )

    # Generate prediction times (from 0 to plot_time_horizon)
    pred_times <- seq(0, plot_time_horizon, length.out = 200)

    # Predict survival probabilities
    surv_pred <- summary(fit, t = pred_times, type = "survival", ci = FALSE)

    # Extract survival probabilities
    surv_values <- surv_pred[[1]]$est

    # Create dataframe
    temp_df <- data.frame(
      time = pred_times,
      surv = surv_values,
      ARM = arm,
      stringsAsFactors = FALSE
    )

    parametric_data <- rbind(parametric_data, temp_df)
  }

  # Calculate number at risk at regular intervals (10% of plot_time_horizon)
  interval <- max(1, round(plot_time_horizon / 10))
  risk_times <- seq(0, plot_time_horizon, by = interval)

  risk_table <- data.frame()

  for (t in risk_times) {
    at_risk <- dataset %>%
      dplyr::filter(SURVTIME >= t) %>%
      dplyr::group_by(ARM) %>%
      dplyr::summarise(n = dplyr::n(), .groups = 'drop') %>%
      dplyr::mutate(Time = t)

    risk_table <- rbind(risk_table, at_risk)
  }

  # Reshape risk table for output
  risk_table_wide <- risk_table %>%
    tidyr::pivot_wider(names_from = ARM, values_from = n, values_fill = 0) %>%
    dplyr::select(Time, dplyr::everything())

  # Create base plot
  p <- ggplot2::ggplot(km_data, ggplot2::aes(x = time, y = surv, color = ARM)) +
    ggplot2::scale_x_continuous(
      limits = c(0, plot_time_horizon),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = 'Kaplan-Meier Curves with Parametric Distribution Overlay',
      subtitle = paste('Distribution:', dist_names[distribution]),
      x = 'Time',
      y = 'Survival Probability',
      color = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 26, face = 'bold', hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 24, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 22),
      axis.title.y = ggplot2::element_text(size = 22),
      axis.text.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 20),
      legend.position = 'bottom',
      legend.text = ggplot2::element_text(size = 18),
      legend.title = ggplot2::element_text(size = 18),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.key.width = ggplot2::unit(2, "cm"),
      panel.grid.major = ggplot2::element_line(color = 'gray90'),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add confidence intervals if requested
  if (conf_int) {
    # For step plots, we need to duplicate points to create step effect in ribbon
    km_data_ribbon <- data.frame()

    for (arm in arm_levels) {
      km_data_arm <- km_data %>% dplyr::filter(ARM == arm)

      # Create step data for ribbon
      n_points <- nrow(km_data_arm)
      step_data <- data.frame()

      for (i in 1:(n_points - 1)) {
        # Add current point
        step_data <- rbind(step_data, km_data_arm[i, ])
        # Add horizontal segment endpoint (next time point with current values)
        next_point <- km_data_arm[i, ]
        next_point$time <- km_data_arm[i + 1, "time"]
        step_data <- rbind(step_data, next_point)
      }
      # Add last point
      step_data <- rbind(step_data, km_data_arm[n_points, ])

      km_data_ribbon <- rbind(km_data_ribbon, step_data)
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = km_data_ribbon,
        ggplot2::aes(ymin = lower, ymax = upper, fill = ARM),
        alpha = 0.2,
        color = NA
      ) +
      ggplot2::scale_fill_manual(values = colors, guide = "none")
  }

  # Add KM step lines
  p <- p +
    ggplot2::geom_step(linewidth = 1)

  # Add parametric overlay
  p <- p +
    ggplot2::geom_line(
      data = parametric_data,
      ggplot2::aes(x = time, y = surv, color = ARM),
      linetype = 'dashed',
      linewidth = 1
    )

  # Add legend for ARM info
  legend_lines <- c()
  for (i in seq_along(arm_levels)) {
    arm <- arm_levels[i]
    n_arm <- sum(dataset$ARM == arm)
    events_arm <- sum(dataset$ARM == arm & dataset$EVENT == 1)
    legend_lines <- c(legend_lines, sprintf("%s: N=%d, Events=%d", arm, n_arm, events_arm))
  }

  legend_text <- paste(legend_lines, collapse = "\n")

  p <- p +
    ggplot2::annotate(
      'text',
      x = plot_time_horizon * 0.65,
      y = 0.95,
      label = legend_text,
      hjust = 0,
      vjust = 1,
      size = 8,
      fontface = 'bold'
    )

  # Return list with plot and risk table
  return(list(
    plot = p,
    risktable = risk_table_wide
  ))
}
