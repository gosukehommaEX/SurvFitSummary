#' Contents of the Word Report
#'
#' @description
#' Internal helper of \code{\link{write_survival_report}}: the elements of
#' the report in the order of the report (title and overview,
#' methods, data overview, diagnostics, standard parametric models,
#' spline-based extrapolation and notes on model fitting). Figures are not
#' drawn here.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param overview List of the overview items (see
#'   \code{.report_front()}).
#' @param conf_int Logical vector: versions of the figures without
#'   (\code{FALSE}) and/or with (\code{TRUE}) confidence intervals.
#' @param n_sim,seed Simulation of the confidence intervals of the models.
#'
#' @return A list of elements made by \code{.report_block()}.
#'
#' @keywords internal
#' @noRd
.report_blocks <- function(x, s, overview, conf_int, n_sim, seed) {
  c(.report_front(x, s, overview),
    .report_methods(x, s),
    .report_data_results(x, s, conf_int),
    .report_diagnostics(x, s, conf_int),
    .report_parametric(x, s, conf_int, n_sim, seed),
    .report_splines(x, s, conf_int, n_sim, seed),
    .report_problems(x))
}
