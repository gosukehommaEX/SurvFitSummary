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
#' @return A single ggplot object containing:
#'   \itemize{
#'     \item Single plot if ARM has 2 levels
#'     \item Horizontal panel if ARM has 3+ levels with aligned y-axes
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
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_hline labs theme_bw theme element_text element_blank element_rect element_line annotate coord_cartesian scale_x_continuous scale_y_continuous scale_shape_manual scale_linetype_manual scale_color_manual guide_legend guides unit
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
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
#' # Creates 3 plots (A vs Control, B vs Control, C vs Control)
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
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required but not installed.")
  }

  # Check ARM levels
  arm_levels <- unique(dataset$ARM)
  if (length(arm_levels) < 2) {
    stop("ARM must have at least 2 levels for proportional hazards test. Found: ",
         length(arm_levels), " level(s)")
  }

  # Handle control_arm specification
  if (!is.null(control_arm)) {
    if (!control_arm %in% arm_levels) {
      stop("Specified control_arm '", control_arm, "' not found in ARM levels: ",
           paste(arm_levels, collapse = ", "))
    }
    # Relevel ARM to set control_arm as reference
    dataset$ARM <- relevel(factor(dataset$ARM), ref = control_arm)
    arm_levels <- levels(dataset$ARM)
    cat("ARM releveled with '", control_arm, "' as reference group\n", sep = "")
  } else {
    # Ensure ARM is a factor
    if (!is.factor(dataset$ARM)) {
      dataset$ARM <- as.factor(dataset$ARM)
      arm_levels <- levels(dataset$ARM)
    }
  }

  # Control arm is the reference (first level)
  control <- arm_levels[1]
  treatment_arms <- arm_levels[-1]

  # Store all schoenfeld data to calculate common y-axis range
  all_schoenfeld_data <- list()

  # First pass: collect all data
  for (trt_arm in treatment_arms) {
    # Create binary comparison dataset: control vs current treatment
    binary_data <- dataset[dataset$ARM %in% c(control, trt_arm), ]
    binary_data$ARM <- droplevels(binary_data$ARM)
    binary_data$ARM <- relevel(binary_data$ARM, ref = control)

    # Fit Cox model and test PH assumption
    cox_fit <- survival::coxph(survival::Surv(SURVTIME, EVENT) ~ ARM, data = binary_data)
    ph_test <- survival::cox.zph(cox_fit)

    # Get ARM coefficient name
    arm_coef <- grep('^ARM', names(coef(cox_fit)), value = TRUE)[1]

    # Extract data
    beta <- coef(cox_fit)[arm_coef]
    pval <- ph_test$table["ARM", 'p']

    # Get residuals (should be a single column for binary comparison)
    if (is.matrix(ph_test$y)) {
      residuals <- ph_test$y[, 1]
    } else {
      residuals <- ph_test$y
    }

    schoenfeld_data <- data.frame(
      Time = ph_test$time,
      Residuals = residuals,
      Beta_t = beta + residuals
    )

    all_schoenfeld_data[[trt_arm]] <- list(
      data = schoenfeld_data,
      beta = beta,
      pval = pval,
      ph_test = ph_test
    )
  }

  # Calculate common y-axis range across all plots
  all_residuals <- unlist(lapply(all_schoenfeld_data, function(x) x$data$Residuals))
  y_limits <- range(all_residuals, na.rm = TRUE)
  y_buffer <- diff(y_limits) * 0.1
  y_limits <- c(y_limits[1] - y_buffer, y_limits[2] + y_buffer)

  # Calculate common x-axis range across all plots
  all_times <- unlist(lapply(all_schoenfeld_data, function(x) x$data$Time))
  x_max <- max(all_times, na.rm = TRUE)

  # Create plot list for each treatment vs control comparison
  plot_list <- list()

  for (trt_arm in treatment_arms) {
    # Extract stored data
    schoenfeld_data <- all_schoenfeld_data[[trt_arm]]$data
    beta <- all_schoenfeld_data[[trt_arm]]$beta
    pval <- all_schoenfeld_data[[trt_arm]]$pval
    ph_test <- all_schoenfeld_data[[trt_arm]]$ph_test

    # Create plot title (without subtitle at individual plot level)
    plot_title <- paste0(trt_arm, " vs ", control)

    # Create annotation text (only beta and p-value)
    pval_text <- ifelse(pval < 0.001, 'p < 0.001', sprintf('p = %.3f', pval))
    annot_color <- ifelse(pval < 0.05, '#D91E49', '#658D1B')
    annot <- sprintf('%s\nβ = %.3f', pval_text, beta)

    # Create plot
    p <- ggplot2::ggplot(schoenfeld_data, ggplot2::aes(x = Time, y = Residuals)) +
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
      ggplot2::scale_shape_manual(
        name = '',
        values = c('Residuals' = 16)
      ) +
      ggplot2::scale_linetype_manual(
        name = '',
        values = c('Reference line (y = 0)' = 'dashed')
      ) +
      ggplot2::scale_color_manual(
        name = '',
        values = c('LOESS smoothing (95% CI)' = '#F0B323')
      ) +
      ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = y_limits) +
      ggplot2::labs(
        title = plot_title,
        x = 'Time',
        y = 'Scaled Schoenfeld Residuals'
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
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
      ggplot2::annotate(
        'text',
        x = max(schoenfeld_data$Time) * 0.7,
        y = y_limits[2] * 0.9,
        label = annot,
        hjust = 0,
        vjust = 1,
        size = 6,
        fontface = 'bold',
        color = annot_color
      )

    # Store plot
    plot_list[[trt_arm]] <- p

    # Print diagnostics
    cat(sprintf(
      '\n=== Diagnostic Information ===\nComparison: %s vs %s\nBeta: %.4f\nP-value: %.4f\nEvents: %d\nResiduals range: [%.4f, %.4f]\n',
      trt_arm,
      control,
      beta,
      pval,
      length(ph_test$time),
      min(schoenfeld_data$Residuals, na.rm = TRUE),
      max(schoenfeld_data$Residuals, na.rm = TRUE)
    ))
  }

  # Arrange plots horizontally if multiple ARM comparisons
  if (length(plot_list) > 1) {
    combined_plot <- patchwork::wrap_plots(
      plot_list,
      ncol = length(plot_list),
      nrow = 1
    ) +
      patchwork::plot_layout(guides = 'collect') +
      patchwork::plot_annotation(
        title = 'Scaled Schoenfeld Residuals',
        subtitle = 'Test for proportional hazards assumption using Grambsch-Therneau method',
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
  } else {
    combined_plot <- plot_list[[1]] +
      patchwork::plot_annotation(
        title = 'Scaled Schoenfeld Residuals',
        subtitle = 'Test for proportional hazards assumption using Grambsch-Therneau method',
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
  }

  return(combined_plot)
}
