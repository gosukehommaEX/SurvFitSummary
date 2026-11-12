#' Smoothed Empirical Hazard with Fitted Hazards
#'
#' @description
#' Plots the smoothed empirical hazard (\code{\link{estimate_empirical_hazard}})
#' with the hazards of fitted parametric or spline models, extrapolated to a
#' time horizon. As in \code{\link{plot_survival_fits}}, either one arm with
#' one or more models (the empirical hazard black, each model with its
#' colour and its AIC in the legend) or one model with one or more arms (each
#' arm with its colour, the empirical hazard solid and the model dashed).
#' The empirical hazard ends at the largest observed time of each arm, and
#' the axes start at (0, 0). For the empirical hazard alone, use
#' \code{\link{plot_empirical_hazard}}.
#'
#' @details
#' The hazards of the models are evaluated at \code{n_points} times from
#' \code{time_horizon / n_points} to \code{time_horizon} (not at time 0,
#' where some hazards are infinite). By default the y axis goes up to 1.2
#' times the largest of the empirical hazard (and its upper confidence
#' limit with \code{conf_int = TRUE}) and of the model hazards after 5\% of
#' the time horizon, so that the very large hazards of some models just
#' after time 0 do not set the scale. With \code{conf_int = TRUE} the
#' confidence bands of the empirical hazards are drawn, and, when one model
#' is plotted, the confidence band of the model simulated as in
#' \code{\link{predict_survival}}.
#'
#' @inheritParams plot_survival_fits
#' @param y_max \code{NULL} (the default, see Details) or the upper limit
#'   of the y axis.
#' @param y_label Label of the y axis. The default is \code{"Hazard rate"}.
#'
#' @return A ggplot object.
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
#' plot_hazard_fits(fits, arms = "Control", conf_int = TRUE)
#'
#' @export
plot_hazard_fits <- function(x,
                             models = NULL,
                             arms = NULL,
                             fit_type = c("independent", "dependent"),
                             conf_int = FALSE,
                             time_horizon = NULL,
                             y_max = NULL,
                             y_label = "Hazard rate",
                             break_by = NULL,
                             n_points = 200,
                             n_sim = 1000,
                             seed = 1) {
  .check_survfitsummary(x)
  fit_type <- match.arg(fit_type)
  sel <- .select_models(x, models, arms, fit_type)
  .plot_fits(x, sel, quantity = "hazard", conf_int = conf_int,
             time_horizon = time_horizon, y_label = y_label, y_max = y_max,
             break_by = break_by, n_points = n_points, n_sim = n_sim,
             seed = seed, risk_table = FALSE)
}
