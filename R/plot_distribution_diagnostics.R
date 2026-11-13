#' Diagnostic Plots of the Parametric Distributions
#'
#' @description
#' Plots transformations of the Kaplan-Meier estimate (or of the smoothed
#' empirical hazard) of each arm that are straight lines when the data
#' follow a distribution:
#' \itemize{
#'   \item \code{"weibull"}: \code{log(-log(S(t)))} against \code{log(t)},
#'     straight for the Weibull (slope = shape; slope 1 for the
#'     exponential). Parallel lines also indicate proportional hazards, so
#'     this is the log cumulative hazard plot of the proportional hazards
#'     assessment.
#'   \item \code{"loglogistic"}: \code{log((1 - S(t)) / S(t))} against
#'     \code{log(t)}, straight for the log-logistic.
#'   \item \code{"lognormal"}: \code{qnorm(1 - S(t))} against
#'     \code{log(t)}, straight for the lognormal.
#'   \item \code{"gompertz"}: \code{log(h(t))} of the smoothed empirical
#'     hazard against \code{t}, straight for the Gompertz.
#' }
#' The Kaplan-Meier estimate is used at the event times where it is
#' strictly between 0 and 1. With \code{regression = TRUE} the least-squares
#' line of each arm is added (dashed).
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param which Character vector with the plots, a subset of
#'   \code{c("weibull", "loglogistic", "lognormal", "gompertz")} (the
#'   default, all four).
#' @param regression Logical; \code{TRUE} (the default) adds the
#'   least-squares line of each arm.
#' @param time_unit Time unit of the data: \code{"months"} (the default),
#'   \code{"days"}, \code{"weeks"} or \code{"years"}.
#' @param legend_position Position of the legend: \code{"bottom"} (the
#'   default), \code{"right"}, \code{"top"}, \code{"left"} or \code{"none"}.
#' @param colours \code{NULL} (the default colours of the arms, the same in
#'   all plots) or colours of the arms: a vector named by arm, for some or
#'   all of the arms, or an unnamed vector with one colour per arm, in the
#'   order of the arms.
#' @param linetypes \code{NULL} (the default, solid lines) or line types of
#'   the curves of the arms, given as for \code{colours} (see
#'   \code{\link{plot_km}}); the least-squares lines are dashed.
#' @inheritParams plot_km
#'
#' @return A ggplot object for a single plot, otherwise a patchwork object
#'   with the plots in the order of \code{which}, two per row, and one
#'   legend.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' plot_distribution_diagnostics(dat)
#' plot_distribution_diagnostics(dat, which = "weibull")
#'
#' @export
plot_distribution_diagnostics <- function(data,
                                          which = c("weibull",
                                                    "loglogistic",
                                                    "lognormal",
                                                    "gompertz"),
                                          regression = TRUE,
                                          time_unit = c("months", "days",
                                                        "weeks", "years"),
                                          legend_position = c("bottom",
                                                              "right",
                                                              "top", "left",
                                                              "none"),
                                          colours = NULL,
                                          linetypes = NULL,
                                          text_sizes = NULL) {
  .check_survival_data(data)
  which <- match.arg(which, several.ok = TRUE)
  time_unit <- match.arg(time_unit)
  legend_position <- match.arg(legend_position)
  text_sizes <- .check_text_sizes(text_sizes)
  data$ARM <- droplevels(data$ARM)
  arms <- levels(data$ARM)
  colours <- .plot_values(colours, arms, .plot_colours(arms = arms),
                          "colour", "colours")
  linetypes <- .plot_values(linetypes, arms, rep("solid", length(arms)),
                            "linetype", "linetypes")
  time_label <- .time_axis(1, time_unit)$label
  log_time_label <- paste0("log(", time_label, ")")

  # Kaplan-Meier estimate at the event times, strictly between 0 and 1
  km <- do.call(rbind, lapply(arms, function(a) {
    d <- data[data$ARM == a, , drop = FALSE]
    fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = d)
    keep <- fit$n.event > 0 & fit$surv > 0 & fit$surv < 1
    data.frame(arm = a, time = fit$time[keep], surv = fit$surv[keep],
               stringsAsFactors = FALSE)
  }))
  panels <- list(
    weibull = list(title = "Weibull / exponential", x = log(km$time),
                   y = log(-log(km$surv)), arm = km$arm,
                   x_label = log_time_label, y_label = "log(-log(S(t)))"),
    loglogistic = list(title = "Log-logistic", x = log(km$time),
                       y = log((1 - km$surv) / km$surv), arm = km$arm,
                       x_label = log_time_label,
                       y_label = "log((1 - S(t)) / S(t))"),
    lognormal = list(title = "Lognormal", x = log(km$time),
                     y = stats::qnorm(1 - km$surv), arm = km$arm,
                     x_label = log_time_label, y_label = "qnorm(1 - S(t))")
  )
  if ("gompertz" %in% which) {
    hz <- estimate_empirical_hazard(data)
    hz <- hz[hz$hazard > 0, , drop = FALSE]
    panels$gompertz <- list(title = "Gompertz", x = hz$time,
                            y = log(hz$hazard), arm = hz$ARM,
                            x_label = time_label, y_label = "log(h(t))",
                            breaks = .time_axis(max(hz$time),
                                                time_unit)$breaks)
  }

  plots <- lapply(which, function(w) {
    panel <- panels[[w]]
    points <- data.frame(arm = factor(panel$arm, levels = arms),
                         x = panel$x, y = panel$y)
    geom <- if (w == "gompertz") ggplot2::geom_line else ggplot2::geom_step
    p <- ggplot2::ggplot(
      points, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$arm)
    ) +
      geom(ggplot2::aes(linetype = .data$arm), linewidth = 0.8)
    if (isTRUE(regression)) {
      lines <- do.call(rbind, lapply(arms, function(a) {
        d <- points[points$arm == a, , drop = FALSE]
        if (nrow(d) < 2L) return(NULL)
        coefs <- stats::coef(stats::lm(y ~ x, data = d))
        ends <- range(d$x)
        data.frame(arm = factor(a, levels = arms), x = ends,
                   y = unname(coefs[1L] + coefs[2L] * ends))
      }))
      p <- p + ggplot2::geom_line(
        data = lines,
        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$arm),
        linetype = "dashed", linewidth = 0.6
      )
    }
    if (!is.null(panel$breaks)) {
      p <- p + ggplot2::scale_x_continuous(breaks = panel$breaks)
    }
    p +
      ggplot2::scale_colour_manual(values = colours, breaks = arms) +
      ggplot2::scale_linetype_manual(values = linetypes, breaks = arms) +
      ggplot2::labs(title = panel$title, x = panel$x_label,
                    y = panel$y_label) +
      .plot_theme(legend_position = legend_position, text_sizes = text_sizes)
  })
  if (length(plots) == 1L) return(plots[[1L]])
  # The collected legend is placed by the theme of the patchwork
  patchwork::wrap_plots(plots, ncol = 2, guides = "collect") &
    ggplot2::theme(legend.position = legend_position)
}
