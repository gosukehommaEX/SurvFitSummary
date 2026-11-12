#' Settings Shared by the Excel Workbook and the Word Report
#'
#' @description
#' Internal helper that checks the arguments shared by
#' \code{\link{write_survival_workbook}} and
#' \code{\link{write_survival_report}} and fills in their defaults: the
#' reporting interval of the Kaplan-Meier tables (the spacing of the time
#' axis of the plots of the data), the landmark times of the predicted
#' survival (6 months and 1, 2, 5 and 10 years, as in the report template),
#' the time horizon (the last landmark time) and the hazard ratios (Cox
#' models without adjustment when none are given).
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
  positive_number <- function(v) {
    is.numeric(v) && length(v) == 1L && is.finite(v) && v > 0
  }
  unit <- x$time_unit
  if (is.null(interval)) {
    interval <- .time_axis(max(x$data$SURVTIME), unit)$by
  } else if (!positive_number(interval)) {
    stop("'interval' must be NULL or a single positive number.",
         call. = FALSE)
  }
  if (is.null(landmarks)) {
    landmarks <- .landmark_times(unit)
  } else {
    if (!is.numeric(landmarks) || length(landmarks) == 0L ||
        any(!is.finite(landmarks)) || any(landmarks <= 0)) {
      stop("'landmarks' must be NULL or a vector of positive numbers.",
           call. = FALSE)
    }
    landmarks <- sort(unique(landmarks))
    landmarks <- data.frame(time = landmarks,
                            label = paste(format(landmarks, trim = TRUE),
                                          unit),
                            stringsAsFactors = FALSE)
  }
  if (is.null(time_horizon)) {
    time_horizon <- max(landmarks$time)
  } else if (!positive_number(time_horizon)) {
    stop("'time_horizon' must be NULL or a single positive number.",
         call. = FALSE)
  }
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
  list(interval = interval, landmarks = landmarks,
       time_horizon = time_horizon, unit = unit,
       hazard_ratios = hazard_ratios)
}
