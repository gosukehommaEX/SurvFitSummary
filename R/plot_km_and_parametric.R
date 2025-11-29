#' Plot Kaplan-Meier Curves with Parametric Distribution Overlay (Fixed for Single ARM)
#'
#' This function creates a combined visualization of Kaplan-Meier survival curves
#' with overlaid parametric distribution fits. It supports multiple treatment arms,
#' displays confidence intervals optionally, and handles multiple time scales.
#' Fixed version handles single ARM datasets correctly.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distribution Character string specifying the parametric distribution.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param conf_int Logical indicating whether to display 95% confidence intervals on
#'   Kaplan-Meier curves. Default is TRUE
#' @param time_scale Character string specifying the unit of time scale in the input dataset.
#'   Must be one of: "day", "week", "month", "year". Default is "week".
#'   Data will be converted internally to weeks for analysis
#' @param time_horizon Numeric value specifying the maximum time point in the specified
#'   time_scale unit. If NULL (default), uses the maximum observed survival time in the
#'   dataset. All outputs (plots and tables) are displayed in weeks
#'
#' @return A list containing:
#'   \describe{
#'     \item{plot}{ggplot object with Kaplan-Meier curves and parametric overlays,
#'       with x-axis in weeks}
#'     \item{risktable}{Data frame showing number at risk at 4-week intervals,
#'       with Time column in weeks}
#'   }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Converts input SURVTIME from time_scale to weeks internally
#'   \item Fits Kaplan-Meier curves using \code{survival::survfit()}
#'   \item Fits parametric survival models independently for each ARM
#'   \item Calculates number at risk at 4-week intervals
#'   \item Returns plot and risk table with all time values in weeks
#'   \item Axes begin at (0, 0) for improved visualization
#'   \item Handles both single and multiple ARM datasets
#' }
#'
#' @importFrom survival Surv survfit
#' @importFrom flexsurv flexsurvreg
#' @importFrom ggplot2 ggplot aes geom_line geom_step geom_ribbon scale_color_manual
#'   scale_fill_manual scale_x_continuous scale_y_continuous labs theme element_text
#'   element_blank element_line coord_cartesian annotate geom_hline
#' @importFrom dplyr group_by summarise n mutate select
#' @importFrom tibble as_tibble
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
#'   conf_int = TRUE,
#'   time_scale = "week",
#'   time_horizon = 100
#' )
#'
#' print(result$plot)
#' print(result$risktable)
#'
#' #' # Example 2: Two-arm trial with input data in weeks
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
#'   conf_int = TRUE,
#'   time_scale = "week",
#'   time_horizon = 100
#' )
#'
#' print(result$plot)
#' print(result$risktable)
#' }

plot_km_and_parametric <- function(dataset,
                                   distribution,
                                   conf_int = TRUE,
                                   time_scale = "week",
                                   time_horizon = NULL) {

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

  # Validate and set time scale conversion factor
  scale_factors <- c(
    "day" = 1 / 7,
    "week" = 1,
    "month" = 30.44 / 7,
    "year" = 365.25 / 7
  )

  if (!time_scale %in% names(scale_factors)) {
    stop("time_scale must be one of: ", paste(names(scale_factors), collapse = ", "))
  }

  scale_factor <- scale_factors[time_scale]

  # Convert SURVTIME to weeks
  dataset <- dataset %>%
    dplyr::mutate(SURVTIME = SURVTIME * scale_factor)

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

  # Set time_horizon in weeks
  if (is.null(time_horizon)) {
    time_horizon_weeks <- max(dataset$SURVTIME, na.rm = TRUE)
  } else {
    time_horizon_weeks <- time_horizon * scale_factor
  }

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

  # Calculate sample size and event counts for each ARM
  arm_summary <- dataset %>%
    dplyr::group_by(ARM) %>%
    dplyr::summarise(
      N = dplyr::n(),
      Events = sum(EVENT),
      .groups = 'drop'
    )

  # Fit Kaplan-Meier curves - FIXED for single ARM
  if (n_arms == 1) {
    # Single ARM: use ~ 1 formula
    km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = dataset)
    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      upper = km_fit$upper,
      lower = km_fit$lower,
      ARM = as.character(arm_levels[1])
    )
  } else {
    # Multiple ARMs: use ~ ARM formula
    km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ ARM, data = dataset)
    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      upper = km_fit$upper,
      lower = km_fit$lower,
      ARM = gsub('ARM=', '', rep(names(km_fit$strata), km_fit$strata))
    )
  }

  # Calculate risk table (number at risk at each 4-week time point)
  time_points_weeks <- seq(0, time_horizon_weeks, by = 4)
  risktable_list <- list()

  for (time_weeks in time_points_weeks) {
    row_data <- data.frame(Time = time_weeks)

    for (arm in arm_levels) {
      arm_data <- dataset[dataset$ARM == arm, ]
      # Count subjects still at risk (SURVTIME >= current time)
      nrisk <- sum(arm_data$SURVTIME >= time_weeks)
      row_data[[as.character(arm)]] <- nrisk
    }

    risktable_list[[length(risktable_list) + 1]] <- row_data
  }

  # Combine risk table rows into data frame
  risktable_df <- do.call(rbind, risktable_list)
  rownames(risktable_df) <- NULL

  # Create base plot
  p <- suppressWarnings({
    ggplot2::ggplot(km_data, ggplot2::aes(x = time, y = surv, color = ARM)) +
      ggplot2::geom_step(linewidth = 1) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(
        title = 'Kaplan-Meier Curves with Parametric Distribution Overlay',
        subtitle = paste('Distribution:', dist_names[distribution]),
        x = 'Time (weeks)',
        y = 'Survival Probability',
        color = NULL,
        fill = NULL
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, time_horizon_weeks),
        ylim = c(0, 1)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 26, face = 'bold', hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 20, hjust = 0.5),
        axis.title.x = ggplot2::element_text(size = 22),
        axis.title.y = ggplot2::element_text(size = 22),
        axis.text.x = ggplot2::element_text(size = 18),
        axis.text.y = ggplot2::element_text(size = 18),
        legend.position = 'bottom',
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16),
        legend.key = ggplot2::element_rect(colour = NA, fill = NA),
        legend.key.width = ggplot2::unit(2, "cm"),
        panel.grid.major = ggplot2::element_line(color = 'gray90'),
        panel.grid.minor = ggplot2::element_blank()
      )
  })

  # Add confidence intervals if requested
  if (conf_int) {
    p <- suppressWarnings({
      p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper, fill = ARM),
          alpha = 0.2,
          color = NA
        ) +
        ggplot2::scale_fill_manual(values = colors)
    })
  }

  # Fit parametric distributions and add overlay for each ARM
  time_seq <- seq(0, time_horizon_weeks, length.out = 200)

  for (i in seq_along(arm_levels)) {
    arm_name <- arm_levels[i]
    arm_subset <- dataset[dataset$ARM == arm_name, ]

    # Fit parametric model independently for this ARM
    param_fit <- flexsurv::flexsurvreg(
      survival::Surv(SURVTIME, EVENT) ~ 1,
      data = arm_subset,
      dist = distribution
    )

    # Get survival predictions
    param_pred <- summary(param_fit, t = time_seq, type = 'survival')[[1]]

    param_data <- data.frame(
      time = time_seq,
      surv = param_pred$est,
      ARM = as.character(arm_name)
    )

    # Add parametric curve overlay
    p <- p +
      ggplot2::geom_line(
        data = param_data,
        ggplot2::aes(x = time, y = surv),
        color = colors[i],
        linetype = 'dashed',
        linewidth = 1.2,
        inherit.aes = FALSE
      )
  }

  # Create annotation text for sample sizes and events
  annotation_lines <- paste0(
    arm_summary$ARM, ': N=', arm_summary$N, ', Events=', arm_summary$Events
  )
  annotation_text <- paste(annotation_lines, collapse = '\n')

  # Add annotation in top-right corner
  p <- p +
    ggplot2::annotate(
      'text',
      x = time_horizon_weeks * 0.98,
      y = 0.98,
      label = annotation_text,
      hjust = 1,
      vjust = 1,
      size = 5
    )

  # Return list with plot and risk table
  return(list(
    plot = p,
    risktable = risktable_df
  ))
}
