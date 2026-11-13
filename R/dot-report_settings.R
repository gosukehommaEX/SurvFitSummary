#' Settings Shared by the Excel Workbook and the Word Report
#'
#' @description
#' Internal helper that checks the arguments shared by
#' \code{\link{write_survival_workbook}} and
#' \code{\link{write_survival_report}} and fills in their defaults: the
#' reporting times with \code{\link{reporting_times}} (the reporting
#' interval of the Kaplan-Meier tables, the landmark times of the predicted
#' survival and the time horizon) and the hazard ratios (Cox models without
#' adjustment when none are given).
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param interval \code{NULL} or the spacing of the reporting times.
#' @param landmarks \code{NULL} or the landmark times.
#' @param time_horizon \code{NULL} or the time horizon.
#' @param hazard_ratios \code{NULL} or a data frame returned by
#'   \code{\link{estimate_hazard_ratio}}.
#'
#' @return A list with \code{interval}, \code{landmarks} (a data frame with
#'   the columns \code{time} and \code{label}), \code{time_horizon},
#'   \code{unit} (for example \code{"months"}) and \code{hazard_ratios}
#'   (\code{NULL} with one arm).
#'
#' @keywords internal
#' @noRd
.report_settings <- function(x, interval, landmarks, time_horizon,
                             hazard_ratios) {
  times <- reporting_times(x, interval = interval, landmarks = landmarks,
                           time_horizon = time_horizon)
  columns <- c("model", "covariates", "strata", "comparison", "hr", "lcl",
               "ucl", "p_value")
  if (length(x$arms) < 2L) {
    hazard_ratios <- NULL
  } else if (is.null(hazard_ratios)) {
    hazard_ratios <- estimate_hazard_ratio(x$data,
                                           conf_level = x$conf_level)
  } else if (!is.data.frame(hazard_ratios) ||
             !all(columns %in% names(hazard_ratios))) {
    stop("'hazard_ratios' must be NULL or a data frame returned by ",
         "estimate_hazard_ratio().", call. = FALSE)
  }
  list(interval = times$interval, landmarks = times$landmarks,
       time_horizon = times$time_horizon, unit = times$time_unit,
       hazard_ratios = hazard_ratios)
}
