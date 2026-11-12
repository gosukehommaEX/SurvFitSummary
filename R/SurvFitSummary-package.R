#' SurvFitSummary: Survival Extrapolation with Standard Parametric and Spline Models
#'
#' @description
#' Prepares time-to-event data with \code{\link{processing_dataset}},
#' summarises them with Kaplan-Meier estimates (\code{\link{summarise_data}},
#' \code{\link{km_table}}), Cox proportional hazards models
#' (\code{\link{estimate_hazard_ratio}}), the Grambsch-Therneau test
#' (\code{\link{test_proportional_hazards}}) and smoothed empirical hazards
#' (\code{\link{estimate_empirical_hazard}}). Test data can be simulated with
#' \code{\link{generate_dummy_survival_data}}.
#'
#' Standard parametric and Royston-Parmar spline models are fitted with
#' \code{\link{fit_survival_models}}, and their results are reported by
#' \code{\link{parameter_table}}, \code{\link{covariance_table}},
#' \code{\link{gof_table}}, \code{\link{predict_survival}},
#' \code{\link{survival_summary}} and \code{\link{spline_knot_table}}.
#'
#' The plots are \code{\link{plot_km}}, \code{\link{plot_empirical_hazard}},
#' \code{\link{plot_survival_fits}}, \code{\link{plot_hazard_fits}},
#' \code{\link{plot_distribution_diagnostics}},
#' \code{\link{plot_schoenfeld}} and \code{\link{plot_qq_survival}}.
#'
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"
