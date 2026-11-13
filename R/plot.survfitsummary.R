#' Plot the Survival Models Fitted by fit_survival_models()
#'
#' @description
#' The \code{plot()} method of the object returned by
#' \code{\link{fit_survival_models}}: plots the fitted survival curves with
#' \code{\link{plot_survival_fits}} (\code{type = "survival"}, the default)
#' or the fitted hazards with \code{\link{plot_hazard_fits}}
#' (\code{type = "hazard"}). These functions plot either one arm with one
#' or more models or one model with one or more arms, so with neither
#' \code{arms} nor a single model given, the models of the control arm (the
#' first arm) are plotted.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param type \code{"survival"} (the default) or \code{"hazard"}.
#' @param models \code{NULL} (the default, all fitted models of
#'   \code{fit_type}) or model labels as in \code{x$models$model}.
#' @param arms \code{NULL} (the default: all arms with a single model,
#'   otherwise the control arm) or the arms to plot.
#' @param ... Other arguments of \code{\link{plot_survival_fits}} or
#'   \code{\link{plot_hazard_fits}}, for example \code{fit_type},
#'   \code{conf_int}, \code{time_horizon}, \code{colours} or
#'   \code{legend_position}.
#'
#' @return The plot returned by \code{\link{plot_survival_fits}} or
#'   \code{\link{plot_hazard_fits}} (a ggplot or patchwork object), drawn
#'   when it is printed.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull",
#'                                                    "llogis"),
#'                             spline_scales = NULL)
#' plot(fits, time_horizon = 120)
#' plot(fits, type = "hazard", models = "Weibull")
#'
#' @export
plot.survfitsummary <- function(x, type = c("survival", "hazard"),
                                models = NULL, arms = NULL, ...) {
  .check_survfitsummary(x)
  type <- match.arg(type)
  # One arm with its models unless the arms or a single model are given
  if (is.null(arms) && length(models) != 1L) arms <- x$arms[1L]
  fun <- if (type == "survival") plot_survival_fits else plot_hazard_fits
  fun(x, models = models, arms = arms, ...)
}
