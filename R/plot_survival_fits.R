#' Kaplan-Meier Curves with Fitted Survival Curves
#'
#' @description
#' Plots the Kaplan-Meier curve with the survival curves of fitted
#' parametric or spline models, extrapolated to a time horizon, with the
#' numbers at risk below the plot. Two layouts are available:
#' \itemize{
#'   \item one arm with one or more models: the Kaplan-Meier curve is black
#'     and each model has its colour, with its AIC in the legend (for
#'     example \code{"Log-logistic (AIC = 1300.7)"});
#'   \item one model with one or more arms: each arm has its colour, the
#'     Kaplan-Meier curve is solid and the model dashed.
#' }
#' The Kaplan-Meier curves end at the largest observed time of each arm,
#' and the axes start at (0, 0).
#'
#' @details
#' With \code{conf_int = TRUE} the pointwise confidence bands of the
#' Kaplan-Meier curves are drawn, and, when one model is plotted, the
#' confidence band of the model simulated as in
#' \code{\link{predict_survival}}. The plot is drawn from the fitted object,
#' so plots with and without confidence bands do not refit the models.
#'
#' With \code{range_times}, a vertical line at each of these times marks
#' the range of the survival of the models plotted, from the lowest to the
#' highest, labelled in percent (for example "5-7\%"), as in the plots of
#' the spline models in the Word report (at 5 and 10 years).
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param models \code{NULL} (the default, all fitted models of
#'   \code{fit_type}) or model labels as in \code{x$models$model}, for
#'   example \code{c("Weibull", "Hazards, 2 knots")}.
#' @param arms \code{NULL} (the default, all arms) or the arms to plot.
#' @param fit_type \code{"independent"} (the default, the model fitted to
#'   each arm) or \code{"dependent"} (the model fitted to both arms with the
#'   arm as a covariate).
#' @param conf_int Logical; \code{TRUE} adds confidence bands (see
#'   Details). The default is \code{FALSE}.
#' @param time_horizon \code{NULL} (the default, the largest observed time
#'   of the plotted arms) or the largest time of the plot, to show the
#'   extrapolation.
#' @param risk_table Logical; \code{TRUE} (the default) adds the numbers at
#'   risk at the breaks of the time axis.
#' @param y_label Label of the y axis, for example \code{"Overall
#'   survival"}. The default is \code{"Survival probability"}.
#' @param break_by \code{NULL} (the default) or the spacing of the breaks
#'   of the time axis.
#' @param n_points Number of times from 0 to the time horizon at which the
#'   models are evaluated. The default is 200.
#' @param n_sim Number of simulated parameter vectors for the confidence
#'   band of a model. The default is 1000.
#' @param seed Seed of the simulation. The default is 1.
#' @param legend_position Position of the legend: \code{"bottom"} (the
#'   default), \code{"right"}, \code{"top"}, \code{"left"} or \code{"none"}.
#'   Below or above the plot, the curves of one arm with several models
#'   are listed in two columns.
#' @param colours \code{NULL} (the default colours) or colours of the
#'   curves. With one arm, the curves are the Kaplan-Meier curve (named
#'   \code{"Kaplan-Meier"}, black by default) and the models (named by
#'   their labels, as in \code{models}); with several arms, the curves are
#'   named by arm. A named vector sets the colours of the curves it names;
#'   an unnamed vector gives one colour per curve, in the order of the
#'   legend. With one arm, the confidence band of the model and the numbers
#'   at risk take the colours of the model and of the Kaplan-Meier curve;
#'   with several arms, the bands and the numbers at risk take the colours
#'   of the arms.
#' @param linetypes \code{NULL} (the default) or line types of the
#'   Kaplan-Meier curve (\code{"Kaplan-Meier"}) and of the models (named by
#'   their labels), in both layouts, given as for \code{colours} (see
#'   \code{\link{plot_km}} for the line types). By default the lines are
#'   solid, except that of the model with several arms (dashed).
#' @inheritParams plot_km
#' @param range_times \code{NULL} (the default) or times at which the range
#'   of the survival of the models is marked (see Details). Times after the
#'   time horizon are left out.
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
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull",
#'                                                    "llogis"),
#'                             spline_scales = NULL)
#' plot_survival_fits(fits, arms = "Control", time_horizon = 120,
#'                    range_times = c(60, 120))
#' plot_survival_fits(fits, models = "Weibull", fit_type = "dependent",
#'                    conf_int = TRUE)
#'
#' @export
plot_survival_fits <- function(x,
                               models = NULL,
                               arms = NULL,
                               fit_type = c("independent", "dependent"),
                               conf_int = FALSE,
                               time_horizon = NULL,
                               risk_table = TRUE,
                               y_label = "Survival probability",
                               break_by = NULL,
                               n_points = 200,
                               n_sim = 1000,
                               seed = 1,
                               legend_position = c("bottom", "right", "top",
                                                   "left", "none"),
                               colours = NULL,
                               linetypes = NULL,
                               text_sizes = NULL,
                               range_times = NULL) {
  .check_survfitsummary(x)
  fit_type <- match.arg(fit_type)
  legend_position <- match.arg(legend_position)
  text_sizes <- .check_text_sizes(text_sizes)
  sel <- .select_models(x, models, arms, fit_type)
  .plot_fits(x, sel, quantity = "survival", conf_int = conf_int,
             time_horizon = time_horizon, y_label = y_label, y_max = NULL,
             break_by = break_by, n_points = n_points, n_sim = n_sim,
             seed = seed, risk_table = risk_table,
             legend_position = legend_position, colours = colours,
             linetypes = linetypes, text_sizes = text_sizes,
             range_times = range_times)
}
