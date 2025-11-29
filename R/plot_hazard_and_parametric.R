#' Plot Smoothed Empirical Hazard with Parametric Distribution Overlay
#'
#' This function creates visualizations of smoothed empirical hazard functions
#' overlaid with parametric distribution hazard curves. For each treatment arm,
#' it generates two plots side-by-side: hazard versus time (natural scale) and
#' hazard versus log(time). Separate panels are created for each treatment arm
#' to assess the visual fit of the parametric model to the observed data.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param distribution Character string specifying the parametric distribution.
#'   Must be one of: "exp", "weibull", "lnorm", "llogis", "gompertz", "gengamma", "gamma"
#' @param conf_int Logical indicating whether to display 95% confidence intervals on
#'   hazard curves. Default is TRUE
#' @param time_scale Character string specifying the unit of time scale in the input dataset.
#'   Must be one of: "day", "week", "month", "year". Default is "week"
#' @param time_horizon Numeric value specifying the maximum time point in the specified
#'   time_scale unit. If NULL (default), uses the maximum observed survival time in the
#'   dataset
#'
#' @return A list containing:
#'   \describe{
#'     \item{plot}{A ggplot object containing multiple panels arranged by ARM.
#'       Each ARM has side-by-side plots for time and log(time) axes with
#'       smoothed empirical hazard (solid line) overlaid with parametric hazard (dashed line)}
#'     \item{hazard_data}{A named list of data frames, one for each treatment arm.
#'       Each data frame contains columns: time, hazard, lower_ci, upper_ci, type, ARM}
#'   }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Converts input SURVTIME from time_scale to weeks internally
#'   \item Calculates smoothed empirical hazard using \code{bshazard::bshazard()}
#'   \item Fits parametric survival models independently for each ARM
#'   \item Creates side-by-side plots (time vs log(time)) for each ARM
#'   \item All time values are displayed in the input time_scale unit
#' }
#'
#' Interpretation guidelines:
#' \itemize{
#'   \item Close alignment between empirical and parametric curves suggests good model fit
#'   \item Large deviations indicate the chosen distribution may not adequately
#'         represent the data for that ARM
#'   \item log(time) axis plots help identify distributional shape violations
#' }
#'
#' @importFrom survival Surv
#' @importFrom bshazard bshazard
#' @importFrom flexsurv flexsurvreg
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_color_manual scale_fill_manual
#'   scale_x_continuous scale_y_continuous labs theme_bw theme element_text element_blank
#'   element_rect element_line guides guide_legend
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @importFrom dplyr filter mutate
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Two-arm trial
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
#' # Plot for Weibull distribution
#' p <- plot_hazard_and_parametric(
#'   dataset = dataset_processed,
#'   distribution = "weibull",
#'   conf_int = TRUE,
#'   time_scale = "week",
#'   time_horizon = 100
#' )
#'
#' print(p$plot)
#' # Access hazard data for each ARM
#' head(p$hazard_data$Treatment)
#' }

plot_hazard_and_parametric <- function(dataset,
                                       distribution,
                                       conf_int = TRUE,
                                       time_scale = "week",
                                       time_horizon = NULL) {

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
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required but not installed.")
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

  # Get ARM levels
  arm_levels <- unique(dataset$ARM)
  n_arms <- length(arm_levels)

  # Color palette for emphasis
  color_smoothed <- '#004C97'  # Blue for smoothed empirical
  color_parametric <- '#F0B323'  # Gold for parametric

  # Base theme for all plots
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 26, face = 'bold', hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 20, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 22),
      axis.title.y = ggplot2::element_text(size = 22),
      axis.text.x = ggplot2::element_text(size = 18),
      axis.text.y = ggplot2::element_text(size = 18),
      legend.position = 'bottom',
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.key.width = ggplot2::unit(2, "cm"),
      panel.grid.major = ggplot2::element_line(color = 'gray90'),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Create individual plots for each ARM
  plot_list <- list()
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

      # Get parametric hazard
      time_seq <- seq(0.01, time_horizon_weeks, length.out = 200)
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

      # Combine data
      hazard_combined <- rbind(hazard_empirical, hazard_parametric)

      # Store hazard data for external use
      hazard_data_list[[as.character(arm_name)]] <- hazard_combined

      # Get y-axis ranges
      y_max <- max(hazard_combined$upper_ci, na.rm = TRUE) * 1.1

      # ===== PLOT 1: Time axis (with legend only for first ARM) =====
      p_time <- suppressWarnings({
        base_plot <- ggplot2::ggplot(hazard_combined, ggplot2::aes(x = time, y = hazard, color = type, linetype = type)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::scale_color_manual(
            values = c('Smoothed Empirical' = color_smoothed, 'Parametric' = color_parametric),
            name = NULL
          ) +
          ggplot2::scale_linetype_manual(
            values = c('Smoothed Empirical' = 'solid', 'Parametric' = 'dashed'),
            name = NULL
          ) +
          ggplot2::scale_x_continuous(limits = c(0, time_horizon_weeks)) +
          ggplot2::scale_y_continuous(limits = c(0, y_max)) +
          ggplot2::labs(
            title = paste(arm_name, '(Time Scale)'),
            x = paste('Time (', time_scale, 's)', sep = ''),
            y = 'Hazard h(t)',
            color = NULL,
            linetype = NULL
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 26, face = 'bold', hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 20, hjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 22),
            axis.title.y = ggplot2::element_text(size = 22),
            axis.text.x = ggplot2::element_text(size = 18),
            axis.text.y = ggplot2::element_text(size = 18),
            legend.position = if (i == 1) 'bottom' else 'none',
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_blank(),
            legend.key = ggplot2::element_rect(colour = NA, fill = NA),
            legend.key.width = ggplot2::unit(2, "cm"),
            panel.grid.major = ggplot2::element_line(color = 'gray90'),
            panel.grid.minor = ggplot2::element_blank()
          )

        # Add confidence intervals if requested
        if (conf_int) {
          base_plot <- base_plot +
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

        base_plot
      })

      # ===== PLOT 2: Log-time axis (without legend) =====
      hazard_combined_log <- hazard_combined %>%
        dplyr::mutate(time_log = log(time))

      p_log_time <- suppressWarnings({
        base_plot <- ggplot2::ggplot(hazard_combined_log, ggplot2::aes(x = time_log, y = hazard, color = type, linetype = type)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::scale_color_manual(
            values = c('Smoothed Empirical' = color_smoothed, 'Parametric' = color_parametric),
            name = NULL
          ) +
          ggplot2::scale_linetype_manual(
            values = c('Smoothed Empirical' = 'solid', 'Parametric' = 'dashed'),
            name = NULL
          ) +
          ggplot2::scale_y_continuous(limits = c(0, y_max)) +
          ggplot2::labs(
            title = paste(arm_name, '(Log-Time Scale)'),
            x = paste('log(Time)'),
            y = 'Hazard h(t)',
            color = NULL,
            linetype = NULL
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 26, face = 'bold', hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 20, hjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 22),
            axis.title.y = ggplot2::element_text(size = 22),
            axis.text.x = ggplot2::element_text(size = 18),
            axis.text.y = ggplot2::element_text(size = 18),
            legend.position = 'none'
          ) +
          ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = 'gray90'),
            panel.grid.minor = ggplot2::element_blank()
          )

        # Add confidence intervals if requested
        if (conf_int) {
          base_plot <- base_plot +
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

        base_plot
      })

      # Combine time and log-time plots side-by-side for this ARM
      # p_log_time has legend.position = 'none', so only p_time legend is shown
      arm_combined_plots <- p_time + p_log_time +
        patchwork::plot_layout(ncol = 2)

      plot_list[[as.character(arm_name)]] <- arm_combined_plots

    }, error = function(e) {
      warning(paste('Failed to process ARM', arm_name, ':', e$message))
    })
  }

  # Check if any plots were created
  if (length(plot_list) == 0) {
    stop("No valid plots could be created for any ARM")
  }

  # Combine all ARM plots vertically with legend collection
  combined_plot <- patchwork::wrap_plots(
    plot_list,
    ncol = 1,
    nrow = n_arms
  ) +
    patchwork::plot_layout(guides = 'collect') +
    patchwork::plot_annotation(
      title = 'Smoothed Empirical Hazard with Parametric Distribution Overlay',
      subtitle = paste('Distribution:', dist_names[distribution]),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 20,
          face = 'bold',
          color = '#666666',
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 16,
          color = '#666666',
          hjust = 0.5
        ),
        legend.position = 'bottom'
      )
    )

  # Return list with plot and hazard data
  return(list(
    plot = combined_plot,
    hazard_data = hazard_data_list
  ))
}
