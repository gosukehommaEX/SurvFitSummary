#' Observed and Fitted Survival or Hazard
#'
#' @description
#' Internal helper of \code{\link{plot_survival_fits}} and
#' \code{\link{plot_hazard_fits}}. With one arm, the Kaplan-Meier curve (or
#' smoothed empirical hazard) is black and each model has its colour, with
#' its AIC in the legend; with several arms and one model, each arm has its
#' colour, the observed curve is solid and the model dashed. The observed
#' curves end at the largest observed time of each arm; the models are
#' extrapolated to \code{time_horizon}.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param sel List returned by \code{.select_models()}.
#' @param quantity \code{"survival"} or \code{"hazard"}.
#' @param conf_int Logical; confidence bands of the observed curves, and of
#'   the model when one model is plotted.
#' @param time_horizon \code{NULL} or the largest time of the plot.
#' @param y_label Label of the y axis.
#' @param y_max \code{NULL} or the upper limit of the y axis (hazard only).
#' @param break_by \code{NULL} or the spacing of the breaks of the time
#'   axis.
#' @param n_points Number of times at which the models are evaluated.
#' @param n_sim,seed Simulation of the confidence band of a model.
#' @param risk_table Logical; numbers at risk below the plot (survival
#'   only).
#' @param legend_position Position of the legend (validated by the caller).
#'
#' @return A ggplot or patchwork object.
#'
#' @keywords internal
#' @noRd
.plot_fits <- function(x, sel, quantity, conf_int, time_horizon, y_label,
                       y_max, break_by, n_points, n_sim, seed, risk_table,
                       legend_position) {
  if (!is.logical(conf_int) || length(conf_int) != 1L || is.na(conf_int)) {
    stop("'conf_int' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(time_horizon) &&
      (!is.numeric(time_horizon) || length(time_horizon) != 1L ||
       !is.finite(time_horizon) || time_horizon <= 0)) {
    stop("'time_horizon' must be NULL or a single positive number.",
         call. = FALSE)
  }
  if (!is.numeric(n_points) || length(n_points) != 1L ||
      !is.finite(n_points) || n_points < 2 || n_points != round(n_points)) {
    stop("'n_points' must be a single whole number of at least 2.",
         call. = FALSE)
  }
  if (!is.null(y_max) &&
      (!is.numeric(y_max) || length(y_max) != 1L || !is.finite(y_max) ||
       y_max <= 0)) {
    stop("'y_max' must be NULL or a single positive number.", call. = FALSE)
  }
  survival_plot <- quantity == "survival"
  data <- x$data[x$data$ARM %in% sel$arms, , drop = FALSE]
  if (is.null(time_horizon)) time_horizon <- max(data$SURVTIME)
  axis <- .time_axis(time_horizon, x$time_unit, break_by)
  model_ci <- conf_int && length(sel$models) == 1L

  if (survival_plot) {
    times <- seq(0, time_horizon, length.out = n_points)
    km <- .km_curve(data, sel$arms, x$conf_level)
    observed <- data.frame(arm = km$arm, time = km$time, value = km$surv,
                           lcl = km$lcl, ucl = km$ucl,
                           stringsAsFactors = FALSE)
    observed_name <- "Kaplan-Meier"
  } else {
    times <- seq(0, time_horizon, length.out = n_points + 1L)[-1L]
    hz <- estimate_empirical_hazard(data, conf_level = x$conf_level)
    observed <- data.frame(arm = hz$ARM, time = hz$time, value = hz$hazard,
                           lcl = hz$lcl, ucl = hz$ucl,
                           stringsAsFactors = FALSE)
    observed_name <- "Smoothed empirical hazard"
  }
  pred <- predict_survival(sel$x, times, type = quantity,
                           conf_int = model_ci, n_sim = n_sim, seed = seed)
  pred <- pred[pred$arm %in% sel$arms, , drop = FALSE]
  if (survival_plot) {
    y_max <- 1
  } else if (is.null(y_max)) {
    # Model hazards very close to time 0 can be very large; they do not set
    # the scale
    values <- c(observed$value, if (conf_int) observed$ucl,
                pred$est[pred$time >= 0.05 * time_horizon])
    y_max <- 1.2 * max(values[is.finite(values)])
  }
  band_of <- function(d) if (survival_plot) .step_ribbon(d) else d
  observed_geom <- if (survival_plot) ggplot2::geom_step else ggplot2::geom_line

  p <- ggplot2::ggplot()
  if (length(sel$arms) == 1L) {
    gof <- gof_table(sel$x)
    aic <- gof$AIC[match(sel$models, gof$model)]
    labels <- paste0(sel$models, " (AIC = ",
                     formatC(aic, format = "f", digits = 1), ")")
    observed_label <- paste0(observed_name, " (", sel$arms, ")")
    curves <- c(observed_label, labels)
    colours <- c(stats::setNames("black", observed_label),
                 stats::setNames(unname(.plot_colours(models = sel$rows)),
                                 labels))
    observed$curve <- factor(observed_label, levels = curves)
    pred$curve <- factor(labels[match(pred$model, sel$models)],
                         levels = curves)
    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(
        data = band_of(observed),
        ggplot2::aes(x = .data$time, ymin = .data$lcl, ymax = .data$ucl),
        fill = "grey60", alpha = 0.3
      )
      if (model_ci) {
        p <- p + ggplot2::geom_ribbon(
          data = pred,
          ggplot2::aes(x = .data$time, ymin = .data$lcl, ymax = .data$ucl),
          fill = colours[[labels]], alpha = 0.2
        )
      }
    }
    p <- p +
      observed_geom(
        data = observed,
        ggplot2::aes(x = .data$time, y = .data$value, colour = .data$curve),
        linewidth = 0.8
      ) +
      ggplot2::geom_line(
        data = pred,
        ggplot2::aes(x = .data$time, y = .data$est, colour = .data$curve),
        linewidth = 0.8
      ) +
      ggplot2::scale_colour_manual(values = colours, breaks = curves) +
      # Below or above the plot, the long labels of the models in two
      # columns, so that the legend is not wider than the plot
      ggplot2::guides(colour = ggplot2::guide_legend(
        ncol = if (legend_position %in% c("bottom", "top")) 2L else 1L,
        byrow = TRUE
      ))
  } else {
    types <- c(observed_name, sel$models)
    colours <- .plot_colours(arms = x$arms)[sel$arms]
    observed$arm <- factor(observed$arm, levels = sel$arms)
    observed$type <- factor(observed_name, levels = types)
    pred$arm <- factor(pred$arm, levels = sel$arms)
    pred$type <- factor(sel$models, levels = types)
    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(
        data = band_of(observed),
        ggplot2::aes(x = .data$time, ymin = .data$lcl, ymax = .data$ucl,
                     fill = .data$arm),
        alpha = 0.15
      )
      if (model_ci) {
        p <- p + ggplot2::geom_ribbon(
          data = pred,
          ggplot2::aes(x = .data$time, ymin = .data$lcl, ymax = .data$ucl,
                       fill = .data$arm),
          alpha = 0.15
        )
      }
    }
    p <- p +
      observed_geom(
        data = observed,
        ggplot2::aes(x = .data$time, y = .data$value, colour = .data$arm,
                     linetype = .data$type),
        linewidth = 0.8
      ) +
      ggplot2::geom_line(
        data = pred,
        ggplot2::aes(x = .data$time, y = .data$est, colour = .data$arm,
                     linetype = .data$type),
        linewidth = 0.8
      ) +
      ggplot2::scale_colour_manual(values = colours, breaks = sel$arms) +
      ggplot2::scale_fill_manual(values = colours, guide = "none") +
      ggplot2::scale_linetype_manual(
        values = stats::setNames(c("solid", "dashed"), types), breaks = types
      )
  }
  p <- p +
    ggplot2::scale_x_continuous(breaks = axis$breaks) +
    ggplot2::coord_cartesian(xlim = c(0, time_horizon), ylim = c(0, y_max),
                             expand = FALSE) +
    ggplot2::labs(x = axis$label, y = y_label) +
    .plot_theme(legend_position = legend_position)
  if (survival_plot) {
    p <- p + ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.2))
    if (isTRUE(risk_table)) {
      # Black numbers with the black Kaplan-Meier curve of one arm
      table_colours <- if (length(sel$arms) == 1L) {
        stats::setNames("black", sel$arms)
      } else {
        .plot_colours(arms = x$arms)[sel$arms]
      }
      p <- .add_risk_table(p, data, sel$arms, axis$breaks, time_horizon,
                           table_colours)
    }
  }
  p
}
