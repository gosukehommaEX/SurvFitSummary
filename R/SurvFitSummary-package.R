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
#' @keywords internal
"_PACKAGE"
