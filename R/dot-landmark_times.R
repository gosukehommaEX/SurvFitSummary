#' Default Landmark Times of the Predicted Survival
#'
#' @description
#' Internal helper that returns the default times at which the predicted
#' survival is reported (6 months and 1, 2, 5 and 10 years), in the time
#' unit of the data (1 year = 365.25 days = 12 months).
#'
#' @param time_unit \code{"days"}, \code{"weeks"}, \code{"months"} or
#'   \code{"years"}.
#'
#' @return A data frame with the columns \code{time} and \code{label}.
#'
#' @keywords internal
#' @noRd
.landmark_times <- function(time_unit) {
  per_year <- switch(time_unit, days = 365.25, weeks = 365.25 / 7,
                     months = 12, years = 1)
  data.frame(time = c(0.5, 1, 2, 5, 10) * per_year,
             label = c("6 months", "1 year", "2 years", "5 years",
                       "10 years"),
             stringsAsFactors = FALSE)
}
