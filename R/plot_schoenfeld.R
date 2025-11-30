#' Plot Scaled Schoenfeld Residuals for Proportional Hazards Assumption Testing
#'
#' This function creates scaled Schoenfeld residual plots to test the proportional
#' hazards assumption in Cox regression models. It supports multiple treatment arms.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param control_arm Character string specifying which ARM level should be used as
#'   the reference (control) group. If NULL (default), uses the first factor level
#'
#' @return A ggplot object:
#'   \itemize{
#'     \item Single plot without facet if ARM has 2 levels
#'     \item Horizontal facet panels if ARM has 3+ levels with shared y-axis
#'   }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Fits Cox proportional hazards model using \code{survival::coxph()}
#'   \item Tests proportional hazards assumption using \code{survival::cox.zph()}
#'   \item Creates scaled Schoenfeld residual plots for each treatment comparison
#'   \item Adds LOESS smoothing curve to visualize time-varying effects
#'   \item Displays beta coefficient and p-value on each plot
#' }
#'
#' The proportional hazards assumption is satisfied if:
#' \itemize{
#'   \item The LOESS curve is approximately horizontal (parallel to y = 0)
#'   \item The p-value from Grambsch-Therneau test is > 0.05
#' }
#'
#' Colors used in plots:
#' \itemize{
#'   \item Points: #004C97 (blue)
#'   \item LOESS curve: #F0B323 (gold)
#'   \item Reference line (y = 0): #D91E49 (red)
#'   \item Green (#658D1B) for p >= 0.05, Red (#D91E49) for p < 0.05
#' }
#'
#' Diagnostic information is printed to console showing:
#' \itemize{
#'   \item ARM coefficient name
#'   \item Beta value
#'   \item P-value from proportional hazards test
#'   \item Number of events
#'   \item Residuals range
#' }
#'
#' @importFrom survival coxph Surv cox.zph
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_hline labs theme_bw theme element_text element_blank element_rect element_line annotate coord_cartesian scale_x_continuous scale_y_continuous scale_shape_manual scale_linetype_manual scale_color_manual guide_legend guides unit facet_grid vars
#' @importFrom dplyr bind_rows mutate group_by summarise filter
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic two-arm analysis
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
#' # Plot
#' p1 <- plot_schoenfeld(dataset = dataset_processed)
#' print(p1)
#'
#' # Example 2: Three-arm trial with specified control
#' dataset2 <- generate_dummy_survival_data(
#'   arm = c("Treatment1", "Treatment2", "Control"),
#'   n = c(100, 100, 100),
#'   hazards = log(2) / c(20, 18, 15),
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
#' # Plot with Control as reference (2 plots: Treatment1 vs Control, Treatment2 vs Control)
#' p2 <- plot_schoenfeld(
#'   dataset = dataset2_processed,
#'   control_arm = "Control"
#' )
#' print(p2)
#'
#' # Example 3: Four-arm trial
#' dataset3 <- generate_dummy_survival_data(
#'   arm = c("A", "B", "C", "Control"),
#'   n = c(100, 100, 100, 100),
#'   hazards = log(2) / c(22, 18, 16, 15),
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
#' # Creates 3 panels (A vs Control, B vs Control, C vs Control)
#' p3 <- plot_schoenfeld(
#'   dataset = dataset3_processed,
#'   control_arm = "Control"
#' )
#' print(p3)
#' }
plot_schoenfeld <- function(dataset,
                            control_arm = NULL) {

  # Load required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }

  # Check ARM levels
  arm_levels <- unique(dataset$ARM)
  if (length(arm_levels) < 2) {
    stop("ARM must have at least 2 levels for proportional hazards test.")
  }

  # Convert ARM to factor if not already
  dataset$ARM <- factor(dataset$ARM, levels = arm_levels)

  # Determine control group
  if (is.null(control_arm)) {
    control <- levels(dataset$ARM)[1]
    warning(paste0("control_arm not specified. Using '", control, "' as reference."))
  } else {
    if (!control_arm %in% arm_levels) {
      stop(paste0("control_arm '", control_arm, "' not found in ARM levels."))
    }
    control <- control_arm
  }

  # Set control as reference level
  dataset$ARM <- relevel(dataset$ARM, ref = control)

  # Get treatment arms (exclude control)
  treatment_arms <- setdiff(levels(dataset$ARM), control)

  # Store all data for combining
  all_schoenfeld_data <- list()

  # Process each treatment arm vs control comparison individually
  for (trt_arm in treatment_arms) {
    # Create subset with only control and this treatment arm
    subset_data <- dataset %>%
      dplyr::filter(ARM %in% c(control, trt_arm))

    # Re-level to ensure control is reference
    subset_data$ARM <- factor(subset_data$ARM, levels = c(control, trt_arm))

    # Fit Cox model for this comparison only
    cox_formula <- survival::Surv(SURVTIME, EVENT) ~ ARM
    fit_cox_subset <- survival::coxph(cox_formula, data = subset_data)

    # Test proportional hazards assumption
    ph_test_subset <- survival::cox.zph(fit_cox_subset)

    # Extract scaled residuals
    if (is.matrix(ph_test_subset$y)) {
      residuals <- ph_test_subset$y[, 1]
    } else {
      residuals <- ph_test_subset$y
    }

    schoenfeld_data <- data.frame(
      Time = ph_test_subset$time,
      Residuals = residuals
    )

    # Get beta coefficient and p-value
    beta <- coef(fit_cox_subset)[1]
    pval <- ph_test_subset$table["ARM", "p"]

    # Store data with comparison label
    schoenfeld_data$Comparison <- paste0(trt_arm, " vs ", control)

    # Store annotation info
    pval_text <- ifelse(pval < 0.001, 'p < 0.001', sprintf('p = %.3f', pval))
    annot_color <- ifelse(pval < 0.05, '#D91E49', '#658D1B')
    annot <- sprintf('%s\n\u03B2 = %.3f', pval_text, beta)

    schoenfeld_data$Annotation <- annot
    schoenfeld_data$Annot_Color <- annot_color

    all_schoenfeld_data[[trt_arm]] <- list(
      data = schoenfeld_data,
      beta = beta,
      pval = pval
    )

    # Print diagnostics
    cat(sprintf(
      '\n=== Diagnostic Information ===\nComparison: %s vs %s\nBeta: %.4f\nP-value: %.4f\nEvents: %d\nResiduals range: [%.4f, %.4f]\n',
      trt_arm,
      control,
      beta,
      pval,
      length(ph_test_subset$time),
      min(schoenfeld_data$Residuals, na.rm = TRUE),
      max(schoenfeld_data$Residuals, na.rm = TRUE)
    ))
  }

  # Combine all data
  combined_data <- dplyr::bind_rows(lapply(all_schoenfeld_data, function(x) x$data))

  # Convert Comparison to factor to control panel order
  combined_data$Comparison <- factor(combined_data$Comparison,
                                     levels = unique(combined_data$Comparison))

  # Calculate common y-axis range
  y_limits <- range(combined_data$Residuals, na.rm = TRUE)
  y_buffer <- diff(y_limits) * 0.1
  y_limits <- c(y_limits[1] - y_buffer, y_limits[2] + y_buffer)

  # Create annotation data frame for each panel
  annot_data <- combined_data %>%
    dplyr::group_by(Comparison, Annotation, Annot_Color) %>%
    dplyr::summarise(
      x_pos = max(Time, na.rm = TRUE) * 0.98,
      y_pos = y_limits[2] * 0.95,
      .groups = 'drop'
    )

  # Create base plot
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = Time, y = Residuals)) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = 0, linetype = 'Reference line (y = 0)'),
      color = '#D91E49',
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = 'Residuals'),
      alpha = 0.5,
      color = '#004C97',
      size = 2
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(color = 'LOESS smoothing (95% CI)'),
      method = 'loess',
      formula = y ~ x,
      se = TRUE,
      fill = '#F0B323',
      alpha = 0.2,
      linewidth = 1
    ) +
    ggplot2::geom_text(
      data = annot_data,
      ggplot2::aes(x = x_pos, y = y_pos, label = Annotation),
      color = annot_data$Annot_Color,
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 5,
      fontface = 'bold',
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c('LOESS smoothing (95% CI)' = '#F0B323'),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c('Residuals' = 16),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c('Reference line (y = 0)' = 'dashed'),
      name = NULL
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0)) +
    ggplot2::labs(
      title = 'Scaled Schoenfeld Residuals',
      subtitle = 'Test for PH assumption using Grambsch-Therneau method',
      x = 'Time',
      y = 'Scaled Schoenfeld Residuals'
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, face = 'bold', hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 18, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      strip.text = ggplot2::element_text(size = 16, face = 'bold'),
      strip.background = ggplot2::element_rect(fill = '#E8E8E8', color = 'black'),
      legend.position = 'bottom',
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.key.width = ggplot2::unit(2, "cm"),
      panel.grid.major = ggplot2::element_line(color = 'gray90'),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2),
      color = ggplot2::guide_legend(order = 3)
    )

  # Add facet_grid only if there are 2 or more treatment arms (3+ total arms)
  if (length(treatment_arms) >= 2) {
    p <- p + ggplot2::facet_grid(
      cols = ggplot2::vars(Comparison),
      scales = "free_x"
    )
  }

  return(p)
}
