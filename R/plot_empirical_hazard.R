#' Smoothed Empirical Hazard
#'
#' @description
#' Plots the smoothed empirical hazard of each selected arm
#' (\code{\link{estimate_empirical_hazard}}), without fitted models, with or
#' without pointwise confidence bands. The axes start at (0, 0).
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param arms \code{NULL} (the default, all arms in one plot) or the arms
#'   to plot.
#' @param conf_int Logical; \code{TRUE} adds the pointwise confidence bands
#'   of \code{bshazard::bshazard()}. The default is \code{FALSE}.
#' @param conf_level Confidence level, a single number in (0, 1). The
#'   default is 0.95.
#' @param time_unit Time unit of the data: \code{"months"} (the default),
#'   \code{"days"}, \code{"weeks"} or \code{"years"}.
#' @param y_label Label of the y axis. The default is \code{"Hazard rate"}.
#' @param break_by \code{NULL} (the default) or the spacing of the breaks
#'   of the time axis.
#' @param legend_position Position of the legend: \code{"bottom"} (the
#'   default), \code{"right"}, \code{"top"}, \code{"left"} or \code{"none"}.
#' @param colours \code{NULL} (the default colours of the arms, the same in
#'   all plots) or colours of the arms: a vector named by arm, for some or
#'   all of the plotted arms, or an unnamed vector with one colour per
#'   plotted arm, in the order of the arms. The confidence bands take the
#'   colours of their arms.
#' @param linetypes \code{NULL} (the default, solid lines) or line types of
#'   the arms, given as for \code{colours} (see \code{\link{plot_km}}).
#' @inheritParams plot_km
#'
#' @return A ggplot object.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' plot_empirical_hazard(dat, conf_int = TRUE)
#'
#' @export
plot_empirical_hazard <- function(data,
                                  arms = NULL,
                                  conf_int = FALSE,
                                  conf_level = 0.95,
                                  time_unit = c("months", "days", "weeks",
                                                "years"),
                                  y_label = "Hazard rate",
                                  break_by = NULL,
                                  legend_position = c("bottom", "right", "top",
                                                      "left", "none"),
                                  colours = NULL,
                                  linetypes = NULL,
                                  text_sizes = NULL) {
  .check_survival_data(data)
  time_unit <- match.arg(time_unit)
  legend_position <- match.arg(legend_position)
  text_sizes <- .check_text_sizes(text_sizes)
  all_arms <- levels(droplevels(data$ARM))
  if (is.null(arms)) arms <- all_arms
  if (!is.character(arms) || length(arms) == 0L ||
      !all(arms %in% all_arms)) {
    stop("'arms' must be NULL or a subset of: ",
         paste(all_arms, collapse = ", "), call. = FALSE)
  }
  arms <- all_arms[all_arms %in% arms]
  colours <- .plot_values(colours, arms, .plot_colours(arms = all_arms)[arms],
                          "colour", "colours")
  linetypes <- .plot_values(linetypes, arms, rep("solid", length(arms)),
                            "linetype", "linetypes")
  data <- data[data$ARM %in% arms, , drop = FALSE]
  hazard <- estimate_empirical_hazard(data, conf_level = conf_level)
  hazard$ARM <- factor(hazard$ARM, levels = arms)
  max_time <- max(data$SURVTIME)
  y_max <- 1.05 * max(if (isTRUE(conf_int)) hazard$ucl else hazard$hazard)
  axis <- .time_axis(max_time, time_unit, break_by)
  p <- ggplot2::ggplot()
  if (isTRUE(conf_int)) {
    p <- p + ggplot2::geom_ribbon(
      data = hazard,
      ggplot2::aes(x = .data$time, ymin = .data$lcl, ymax = .data$ucl,
                   fill = .data$ARM),
      alpha = 0.2
    )
  }
  p +
    ggplot2::geom_line(
      data = hazard,
      ggplot2::aes(x = .data$time, y = .data$hazard, colour = .data$ARM,
                   linetype = .data$ARM),
      linewidth = 0.8
    ) +
    ggplot2::scale_colour_manual(values = colours, breaks = arms) +
    ggplot2::scale_linetype_manual(values = linetypes, breaks = arms) +
    ggplot2::scale_fill_manual(values = colours, guide = "none") +
    ggplot2::scale_x_continuous(breaks = axis$breaks) +
    ggplot2::coord_cartesian(xlim = c(0, max_time), ylim = c(0, y_max),
                             expand = FALSE) +
    ggplot2::labs(x = axis$label, y = y_label) +
    .plot_theme(legend_position = legend_position, text_sizes = text_sizes)
}
