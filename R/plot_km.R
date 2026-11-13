#' Kaplan-Meier Curves
#'
#' @description
#' Plots the Kaplan-Meier curve of each selected arm, with or without
#' pointwise confidence bands, and the numbers at risk below the plot. Each
#' curve ends at the largest observed time of its arm, and the axes start at
#' (0, 0).
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param arms \code{NULL} (the default, all arms in one plot) or the arms
#'   to plot, for example a single arm.
#' @param conf_int Logical; \code{TRUE} adds the pointwise confidence bands
#'   (log transformation, as in \code{\link{km_table}}). The default is
#'   \code{FALSE}.
#' @param conf_level Confidence level, a single number in (0, 1). The
#'   default is 0.95.
#' @param risk_table Logical; \code{TRUE} (the default) adds the numbers at
#'   risk at the breaks of the time axis.
#' @param time_unit Time unit of the data, for the label of the time axis
#'   and the spacing of its breaks: \code{"months"} (the default),
#'   \code{"days"}, \code{"weeks"} or \code{"years"}.
#' @param y_label Label of the y axis, for example \code{"Overall
#'   survival"}. The default is \code{"Survival probability"}.
#' @param break_by \code{NULL} (the default, chosen from usual intervals of
#'   the time unit) or the spacing of the breaks of the time axis.
#' @param legend_position Position of the legend: \code{"bottom"} (the
#'   default), \code{"right"}, \code{"top"}, \code{"left"} or \code{"none"}.
#'
#' @return A ggplot object, or with \code{risk_table = TRUE} a patchwork
#'   object whose first element is the plot and second the table.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' plot_km(dat)
#' plot_km(dat, arms = "Treatment", conf_int = TRUE)
#'
#' @export
plot_km <- function(data,
                    arms = NULL,
                    conf_int = FALSE,
                    conf_level = 0.95,
                    risk_table = TRUE,
                    time_unit = c("months", "days", "weeks", "years"),
                    y_label = "Survival probability",
                    break_by = NULL,
                    legend_position = c("bottom", "right", "top",
                                        "left", "none")) {
  .check_survival_data(data)
  time_unit <- match.arg(time_unit)
  legend_position <- match.arg(legend_position)
  all_arms <- levels(droplevels(data$ARM))
  if (is.null(arms)) arms <- all_arms
  if (!is.character(arms) || length(arms) == 0L ||
      !all(arms %in% all_arms)) {
    stop("'arms' must be NULL or a subset of: ",
         paste(all_arms, collapse = ", "), call. = FALSE)
  }
  arms <- all_arms[all_arms %in% arms]
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  data <- data[data$ARM %in% arms, , drop = FALSE]
  max_time <- max(data$SURVTIME)
  axis <- .time_axis(max_time, time_unit, break_by)
  colours <- .plot_colours(arms = all_arms)[arms]
  km <- .km_curve(data, arms, conf_level)
  km$arm <- factor(km$arm, levels = arms)
  p <- ggplot2::ggplot()
  if (isTRUE(conf_int)) {
    band <- .step_ribbon(km)
    band$arm <- factor(band$arm, levels = arms)
    p <- p + ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = .data$time, ymin = .data$lcl, ymax = .data$ucl,
                   fill = .data$arm),
      alpha = 0.2
    )
  }
  p <- p +
    ggplot2::geom_step(
      data = km,
      ggplot2::aes(x = .data$time, y = .data$surv, colour = .data$arm),
      linewidth = 0.8
    ) +
    ggplot2::scale_colour_manual(values = colours, breaks = arms) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::scale_x_continuous(breaks = axis$breaks) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    ggplot2::coord_cartesian(xlim = c(0, max_time), ylim = c(0, 1),
                             expand = FALSE) +
    ggplot2::labs(x = axis$label, y = y_label) +
    .plot_theme(legend_position = legend_position)
  if (isTRUE(risk_table)) {
    p <- .add_risk_table(p, data, arms, axis$breaks, max_time, colours)
  }
  p
}
