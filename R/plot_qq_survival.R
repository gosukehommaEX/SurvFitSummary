#' Plot Q-Q Plots for Survival Time Distributions
#'
#' This function creates quantile-quantile (Q-Q) plots to compare survival time
#' distributions between treatment arms. For 2-arm trials, it creates a single
#' Q-Q plot. For 3+ arms, it creates faceted Q-Q plots comparing each treatment
#' arm against the control arm.
#'
#' @param dataset A data frame containing survival analysis data with columns:
#'   ARM (treatment arm), SURVTIME (survival time), EVENT (event indicator)
#' @param control_arm Character string specifying which ARM level should be used as
#'   the reference (control) group for comparison. If NULL (default), uses the
#'   first factor level. Required when ARM has 3+ levels
#'
#' @return A ggplot object containing:
#'   \itemize{
#'     \item Single Q-Q plot if ARM has 2 levels
#'     \item Faceted Q-Q plots if ARM has 3+ levels, each comparing
#'           a treatment arm against the control arm
#'     \item NULL if ARM has only 1 level (with informative message)
#'   }
#'   Each plot includes:
#'   \itemize{
#'     \item Observed quantiles (red points)
#'     \item Identity line y = x (dashed black line)
#'     \item Linear regression fit with 95% CI (gray solid line with shaded area)
#'   }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Checks if ARM has at least 2 levels; returns NULL with message if not
#'   \item Fits Kaplan-Meier curves for each ARM using \code{survival::survfit()}
#'   \item Extracts survival quantiles at percentiles: 0.05, 0.10, ..., 0.95
#'   \item Creates Q-Q plots comparing quantiles between arms
#'   \item Adds linear regression fit to assess distributional similarity
#' }
#'
#' Interpretation guidelines:
#' \itemize{
#'   \item \strong{Points near identity line}: Survival distributions are similar
#'   \item \strong{Points above identity line}: Y-arm has longer survival times
#'   \item \strong{Points below identity line}: Y-arm has shorter survival times
#'   \item \strong{Non-linear pattern}: Different distributional shapes
#' }
#'
#' Colors and styling:
#' \itemize{
#'   \item Observed quantiles: #D91E49 (red points)
#'   \item Identity line: Black dashed
#'   \item Regression fit: #939597 (gray) with 95% CI
#' }
#'
#' @importFrom survival survfit Surv
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_smooth geom_segment labs theme_bw
#'   theme element_text scale_shape_manual scale_linetype_manual guide_legend guides
#'   facet_grid vars element_rect element_blank element_line unit scale_x_continuous scale_y_continuous
#' @importFrom dplyr filter bind_rows
#' @importFrom purrr map_dfr
#' @importFrom tidyr pivot_wider
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Two-arm trial
#' dataset1 <- generate_dummy_survival_data(
#'   arm = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazards = log(2) / c(20, 15),
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
#' # Single Q-Q plot
#' p1 <- plot_qq_survival(dataset1_processed)
#' print(p1)
#'
#' # Example 2: Three-arm trial with specified control
#' dataset2 <- generate_dummy_survival_data(
#'   arm = c("Treatment1", "Treatment2", "Control"),
#'   n = c(100, 100, 100),
#'   hazards = log(2) / c(22, 18, 15),
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
#' # Faceted Q-Q plots (Treatment1 vs Control, Treatment2 vs Control)
#' p2 <- plot_qq_survival(
#'   dataset = dataset2_processed,
#'   control_arm = "Control"
#' )
#' print(p2)
#'
#' # Example 3: Four-arm trial
#' dataset3 <- generate_dummy_survival_data(
#'   arm = c("TrtA", "TrtB", "TrtC", "Control"),
#'   n = c(80, 80, 80, 80),
#'   hazards = log(2) / c(24, 20, 18, 15),
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
#' # Three faceted Q-Q plots (TrtA vs Control, TrtB vs Control, TrtC vs Control)
#' p3 <- plot_qq_survival(
#'   dataset = dataset3_processed,
#'   control_arm = "Control"
#' )
#' print(p3)
#'
#' # Example 4: Single ARM (returns NULL with message)
#' dataset4 <- generate_dummy_survival_data(
#'   arm = c("DrugA"),
#'   n = c(200),
#'   hazards = log(2) / 18,
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
#' # Returns NULL with message
#' p4 <- plot_qq_survival(dataset4_processed)
#' # Message: "Q-Q plot requires at least 2 ARM levels for comparison.
#' #          Found: 1 level(s). Returning NULL."
#' }
plot_qq_survival <- function(dataset,
                             control_arm = NULL) {

  # Load required packages
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }

  # Check ARM levels
  arm_levels <- unique(dataset$ARM)
  if (length(arm_levels) < 2) {
    message("Q-Q plot requires at least 2 ARM levels for comparison. Found: ",
            length(arm_levels), " level(s). Returning NULL.")
    return(NULL)
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

  # Define percentiles for quantile extraction
  probs <- seq(0.05, 0.95, by = 0.05)

  # Collect all quantile data for all comparisons
  all_qq_data <- dplyr::bind_rows(
    lapply(treatment_arms, function(trt_arm) {
      # Extract quantiles for control and treatment arms
      qq_data <- purrr::map_dfr(c(control, trt_arm), function(arm) {
        arm_data <- dataset %>% dplyr::filter(ARM == arm)
        km_fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = arm_data)
        quantiles <- stats::quantile(km_fit, probs = probs)$quantile
        data.frame(
          ARM = arm,
          percentile = probs,
          quantile = quantiles,
          stringsAsFactors = FALSE
        )
      }) %>%
        tidyr::pivot_wider(names_from = ARM, values_from = quantile)

      # Add comparison label for faceting
      qq_data$Comparison <- paste0(trt_arm, ' vs ', control)
      qq_data$Control_Quantile <- qq_data[[control]]
      qq_data$Treatment_Quantile <- qq_data[[trt_arm]]

      return(qq_data[, c("Comparison", "percentile", "Control_Quantile", "Treatment_Quantile")])
    })
  )

  # Determine axis limits (same for all facets)
  all_values <- c(all_qq_data$Control_Quantile, all_qq_data$Treatment_Quantile)
  axis_max <- max(all_values, na.rm = TRUE)

  # Create combined plot with faceting
  p <- ggplot2::ggplot(all_qq_data, ggplot2::aes(x = Control_Quantile, y = Treatment_Quantile)) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = 'black',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(linetype = 'Regression fit'),
      method = 'lm',
      se = TRUE,
      color = '#939597',
      fill = '#939597',
      alpha = 0.2,
      linewidth = 0.8,
      formula = y ~ x,
      fullrange = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = 'Observed quantiles'),
      size = 3,
      color = '#D91E49'
    ) +
    ggplot2::geom_segment(
      data = data.frame(
        x = -Inf,
        xend = -Inf,
        y = -Inf,
        yend = -Inf,
        Comparison = unique(all_qq_data$Comparison)[1]
      ),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, linetype = 'Identity line (y = x)'),
      color = 'black',
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = '',
      values = c('Observed quantiles' = 16)
    ) +
    ggplot2::scale_linetype_manual(
      name = '',
      values = c(
        'Identity line (y = x)' = 'dashed',
        'Regression fit' = 'solid'
      )
    ) +
    ggplot2::scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ggplot2::labs(
      title = 'Q-Q Plots of Survival Percentiles',
      subtitle = 'Comparing survival time distributions across treatment arms',
      x = paste0('Quantiles (', control, ')'),
      y = 'Quantiles (Treatment)'
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
      linetype = ggplot2::guide_legend(
        order = 1,
        override.aes = list(linewidth = 0.8)
      ),
      shape = ggplot2::guide_legend(order = 2)
    )

  # Add faceting only if there are multiple comparisons
  if (length(treatment_arms) > 1) {
    p <- p +
      ggplot2::facet_grid(
        cols = ggplot2::vars(Comparison)
      )
  }

  return(p)
}
