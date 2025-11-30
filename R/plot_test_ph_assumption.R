#' Plot Diagnostic Plots for Parametric Distribution Assumption Testing
#'
#' This function creates diagnostic plots to assess the appropriateness of
#' parametric survival distributions. It generates three plots: Weibull/Exponential
#' diagnostic (log(-log(S(t))) vs log(t)), Gompertz diagnostic (log(-log(S(t))) vs t),
#' and Log-normal diagnostic (Î¦^-1(1-S(t)) vs log(t)).
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param add_regression_line Logical indicating whether to add linear regression
#'   lines with confidence intervals to the plots. Default is FALSE
#'
#' @return A ggplot object containing a single row of 3 diagnostic plots with shared legend at bottom
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Fits Kaplan-Meier curves using \code{survival::survfit()}
#'   \item Transforms survival estimates for diagnostic plotting
#'   \item Creates three diagnostic plots to assess distributional assumptions
#' }
#'
#' Interpretation guidelines:
#' \itemize{
#'   \item \strong{Weibull/Exponential}: If log(-log(S(t))) vs log(t) is approximately
#'         linear, Weibull (or Exponential if parallel) is appropriate
#'   \item \strong{Gompertz}: If log(-log(S(t))) vs t is approximately linear,
#'         Gompertz distribution is appropriate
#'   \item \strong{Log-normal}: If Î¦^-1(1-S(t)) vs log(t) is approximately linear,
#'         Log-normal distribution is appropriate
#' }
#'
#' Colors are automatically assigned based on number of ARM levels:
#' \itemize{
#'   \item 1 ARM: #004C97 (blue)
#'   \item 2 ARMs: #004C97 (blue), #F0B323 (gold)
#'   \item 3 ARMs: #004C97 (blue), #F0B323 (gold), #658D1B (green)
#'   \item 4 ARMs: #004C97, #F0B323, #658D1B, #A62B4E
#'   \item 5 ARMs: #004C97, #F0B323, #658D1B, #A62B4E, #D91E49
#'   \item 6+ ARMs: #004C97, #F0B323, #658D1B, #A62B4E, #D91E49, #939597
#' }
#'
#' @importFrom survival survfit Surv
#' @importFrom ggplot2 ggplot aes geom_step geom_smooth scale_color_manual scale_fill_manual
#'   scale_x_continuous scale_y_continuous labs theme_minimal theme element_text
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @importFrom dplyr filter mutate
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Single ARM analysis
#' dataset1 <- generate_dummy_survival_data(
#'   arm = c("DrugA"),
#'   n = c(200),
#'   hazards = log(2) / 18,
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
#' # Plot
#' p1 <- plot_test_ph_assumption(dataset1_processed)
#' print(p1)
#'
#' # Example 2: Two-arm analysis
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
#' # Plot without regression lines
#' p2 <- plot_test_ph_assumption(dataset2_processed)
#' print(p2)
#'
#' # Plot with regression lines
#' p2_reg <- plot_test_ph_assumption(dataset2_processed, add_regression_line = TRUE)
#' print(p2_reg)
#'
#' # Example 3: Three-arm trial
#' dataset3 <- generate_dummy_survival_data(
#'   arm = c("Treatment1", "Treatment2", "Control"),
#'   n = c(100, 100, 100),
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
#' # Plot
#' p3 <- plot_test_ph_assumption(dataset3_processed)
#' print(p3)
#'
#' # Example 4: Four-arm trial with regression lines
#' dataset4 <- generate_dummy_survival_data(
#'   arm = c("TrtA", "TrtB", "TrtC", "Control"),
#'   n = c(75, 75, 75, 75),
#'   hazards = log(2) / c(24, 20, 18, 15),
#'   dropout_per_year = 0.1,
#'   seed = 101
#' )
#'
#' dataset4_processed <- processing_dataset(
#'   dataset = dataset4,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = NULL
#' )
#'
#' # Plot
#' p4 <- plot_test_ph_assumption(dataset4_processed, add_regression_line = TRUE)
#' print(p4)
#' }
plot_test_ph_assumption <- function(dataset,
                                    add_regression_line = FALSE) {

  # Load required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required but not installed.")
  }

  # Define color palette
  color_palette <- c('#004C97', '#F0B323', '#658D1B', '#A62B4E', '#D91E49', '#939597')

  # Get ARM levels
  arm_levels <- unique(dataset$ARM)
  n_arms <- length(arm_levels)

  # Fit Kaplan-Meier
  if (n_arms > 1) {
    km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ ARM, data = dataset)

    # Extract and transform Kaplan-Meier estimates
    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      arm = gsub('ARM=', '', rep(names(km_fit$strata), km_fit$strata))
    )
  } else {
    # Single ARM case: no stratification in survfit
    km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = dataset)

    # Extract and transform Kaplan-Meier estimates
    km_data <- data.frame(
      time = km_fit$time,
      surv = km_fit$surv,
      arm = as.character(arm_levels[1])
    )
  }

  # Filter and transform data
  km_data <- km_data %>%
    dplyr::filter(
      time > 0,      # Required for log(time)
      surv > 0,      # Required for log(-log(surv))
      surv < 1       # Required for qnorm(1 - surv)
    ) %>%
    dplyr::mutate(
      log_time = log(time),
      cloglog = log(-log(surv)),
      probit = stats::qnorm(1 - surv)
    )

  # Check if we have valid data
  if (nrow(km_data) == 0) {
    stop(
      'No valid data points after transformation. ',
      'Check survival times and event indicators.',
      call. = FALSE
    )
  }

  # Assign colors
  colors <- color_palette[1:n_arms]
  names(colors) <- as.character(arm_levels)

  # Get x-axis and y-axis ranges for alignment
  log_time_range <- range(km_data$log_time, na.rm = TRUE)
  time_range <- range(km_data$time, na.rm = TRUE)
  cloglog_range <- range(km_data$cloglog, na.rm = TRUE)
  probit_range <- range(km_data$probit, na.rm = TRUE)

  # Add buffer to y-axis ranges
  cloglog_buffer <- diff(cloglog_range) * 0.1
  cloglog_limits <- c(cloglog_range[1] - cloglog_buffer, cloglog_range[2] + cloglog_buffer)

  probit_buffer <- diff(probit_range) * 0.1
  probit_limits <- c(probit_range[1] - probit_buffer, probit_range[2] + probit_buffer)

  # Base theme for all plots
  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, face = 'bold', hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 18, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      legend.position = 'bottom',
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.key.width = ggplot2::unit(2, "cm"),
      panel.grid.major = ggplot2::element_line(color = 'gray90'),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Plot 1: Weibull/Exponential diagnostic
  p1 <- ggplot2::ggplot(km_data, ggplot2::aes(x = log_time, y = cloglog, color = arm, fill = arm)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(limits = log_time_range) +
    ggplot2::scale_y_continuous(limits = cloglog_limits) +
    ggplot2::labs(
      title = 'Weibull/Exponential',
      x = 'log(Time)',
      y = 'log(-log(S(Time)))',
      color = 'Treatment',
      fill = 'Treatment'
    ) +
    base_theme

  # Add regression line if requested
  if (add_regression_line) {
    p1 <- p1 + ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = TRUE,
      linetype = 'dashed',
      linewidth = 0.6,
      alpha = 0.2
    )
  }

  # Plot 2: Gompertz diagnostic
  p2 <- ggplot2::ggplot(km_data, ggplot2::aes(x = time, y = cloglog, color = arm, fill = arm)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(limits = time_range) +
    ggplot2::scale_y_continuous(limits = cloglog_limits) +
    ggplot2::labs(
      title = 'Gompertz',
      x = 'Time',
      y = 'log(-log(S(Time)))',
      color = 'Treatment',
      fill = 'Treatment'
    ) +
    base_theme

  # Add regression line if requested
  if (add_regression_line) {
    p2 <- p2 + ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = TRUE,
      linetype = 'dashed',
      linewidth = 0.6,
      alpha = 0.2
    )
  }

  # Plot 3: Log-normal diagnostic
  p3 <- ggplot2::ggplot(km_data, ggplot2::aes(x = log_time, y = probit, color = arm, fill = arm)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(limits = log_time_range) +
    ggplot2::scale_y_continuous(limits = probit_limits) +
    ggplot2::labs(
      title = 'Log-normal',
      x = 'log(Time)',
      y = expression(Phi^{-1}*'(1 - S(Time))'),
      color = 'Treatment',
      fill = 'Treatment'
    ) +
    base_theme

  # Add regression line if requested
  if (add_regression_line) {
    p3 <- p3 + ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = TRUE,
      linetype = 'dashed',
      linewidth = 0.6,
      alpha = 0.2
    )
  }

  # Arrange 3 plots horizontally
  combined_plot <- patchwork::wrap_plots(
    list(p1, p2, p3),
    ncol = 3,
    nrow = 1
  ) +
    patchwork::plot_layout(guides = 'collect') +
    patchwork::plot_annotation(
      title = 'Parametric Distribution Diagnostic Plots',
      subtitle = 'Assess linearity for distributional assumptions',
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 20,
          face = 'bold',
          color = '#666666',
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          size = 20,
          color = '#666666',
          hjust = 0.5
        ),
        legend.position = 'bottom'
      )
    )

  return(combined_plot)
}
