#' Breaks and Label of the Time Axis
#'
#' @description
#' Internal helper that chooses the spacing of the breaks of the time axis
#' from usual intervals of the time unit (for months 1, 2, 3, 6, 12, 24, 60,
#' 120 or 240), the smallest one that gives at most 10 intervals up to
#' \code{max_time}, unless \code{break_by} is given. The numbers at risk of
#' the plots are reported at these breaks.
#'
#' @param max_time Largest time of the axis.
#' @param time_unit \code{"days"}, \code{"weeks"}, \code{"months"} or
#'   \code{"years"}.
#' @param break_by \code{NULL} or the spacing of the breaks.
#'
#' @return A list with \code{breaks} (from 0 up to \code{max_time}) and
#'   \code{label} (for example \code{"Time (Months)"}).
#'
#' @keywords internal
#' @noRd
.time_axis <- function(max_time, time_unit, break_by = NULL) {
  if (is.null(break_by)) {
    candidates <- switch(
      time_unit,
      days = c(7, 14, 30, 60, 90, 180, 365, 730, 1825, 3650),
      weeks = c(1, 2, 4, 8, 13, 26, 52, 104, 260, 520),
      months = c(1, 2, 3, 6, 12, 24, 60, 120, 240),
      years = c(0.25, 0.5, 1, 2, 5, 10, 20)
    )
    suitable <- candidates[max_time / candidates <= 10]
    break_by <- if (length(suitable) > 0L) {
      suitable[1L]
    } else {
      candidates[length(candidates)]
    }
  } else if (!is.numeric(break_by) || length(break_by) != 1L ||
             !is.finite(break_by) || break_by <= 0) {
    stop("'break_by' must be NULL or a single positive number.",
         call. = FALSE)
  }
  list(
    breaks = seq(0, max_time * (1 + 1e-10), by = break_by),
    label = paste0("Time (", toupper(substr(time_unit, 1L, 1L)),
                   substring(time_unit, 2L), ")")
  )
}
