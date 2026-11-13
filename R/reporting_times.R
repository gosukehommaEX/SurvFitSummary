#' Reporting Times of the Tables, Predictions and Plots
#'
#' @description
#' Returns the reporting times that \code{\link{write_survival_workbook}}
#' and \code{\link{write_survival_report}} use, with their defaults filled
#' in, so that tables, predictions and plots made with the other functions
#' of the package can use the same times as the workbook and the report:
#' \itemize{
#'   \item the spacing of the reporting times of the Kaplan-Meier tables, by
#'     default the spacing of the time axis of the plots of the data (the
#'     smallest of the usual intervals of the time unit, for months 1, 2,
#'     3, 6, 12, 24, 60, 120 or 240, that gives at most 10 intervals up to
#'     the largest observed time);
#'   \item the landmark times of the predicted survival, by default 6 months
#'     and 1, 2, 5 and 10 years in the time unit of the data (1 year =
#'     365.25 days = 12 months), labelled \code{"6 months"}, \code{"1 year"},
#'     \code{"2 years"}, \code{"5 years"} and \code{"10 years"};
#'   \item the time horizon of the extrapolation and of the restricted mean,
#'     by default the last landmark time.
#' }
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param interval \code{NULL} (the default) or the spacing of the
#'   reporting times, a positive number in the time unit of the data.
#' @param landmarks \code{NULL} (the default) or the landmark times, a
#'   vector of positive numbers in the time unit of the data. They are
#'   sorted and labelled with the time unit, for example \code{"36 months"}.
#' @param time_horizon \code{NULL} (the default) or the time horizon, a
#'   positive number in the time unit of the data.
#'
#' @return A list with \code{interval} (number), \code{landmarks} (data
#'   frame with the columns \code{time} and \code{label}),
#'   \code{time_horizon} (number) and \code{time_unit} (the time unit of
#'   \code{x}, for example \code{"months"}).
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = NULL)
#' times <- reporting_times(fits)
#' times
#' km_table(dat, interval = times$interval)
#' predict_survival(fits, times = times$landmarks$time)
#' reporting_times(fits, interval = 6, landmarks = c(24, 36),
#'                 time_horizon = 60)
#'
#' @export
reporting_times <- function(x,
                            interval = NULL,
                            landmarks = NULL,
                            time_horizon = NULL) {
  .check_survfitsummary(x)
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
  list(interval = interval, landmarks = landmarks,
       time_horizon = time_horizon, time_unit = unit)
}
