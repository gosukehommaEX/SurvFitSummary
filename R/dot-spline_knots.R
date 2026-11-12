#' Knots of a Royston-Parmar Spline Model
#'
#' @description
#' Internal helper that returns the knots of a spline model on the log time
#' scale. The boundary knots are the smallest and largest log event times.
#' Without \code{knot_times}, the internal knots are the equally spaced
#' quantiles (type 7) of the log event times, as in
#' \code{flexsurv::flexsurvspline()}: the median for one knot, the 1/3 and
#' 2/3 quantiles for two knots, and the quartiles for three knots.
#'
#' @param time Observed times.
#' @param event Event indicator (1 = event, 0 = censored).
#' @param n_knots Number of internal knots, used when \code{knot_times} is
#'   \code{NULL}.
#' @param knot_times \code{NULL} or the internal knots on the time scale.
#'
#' @return A list with \code{internal} and \code{boundary} (knots on the log
#'   time scale) and \code{quantile} (the probability of each internal knot;
#'   for knots given in \code{knot_times}, the proportion of log event times
#'   at or below the knot).
#'
#' @keywords internal
#' @noRd
.spline_knots <- function(time, event, n_knots, knot_times = NULL) {
  log_events <- log(time[event == 1L])
  if (length(unique(log_events)) < 2L) {
    stop("At least two distinct event times are needed for spline models.",
         call. = FALSE)
  }
  boundary <- range(log_events)
  if (is.null(knot_times)) {
    probs <- seq(0, 1, length.out = n_knots + 2L)[-c(1L, n_knots + 2L)]
    internal <- unname(stats::quantile(log_events, probs, type = 7))
    quantile <- probs
  } else {
    internal <- log(knot_times)
    quantile <- vapply(internal, function(k) mean(log_events <= k),
                       numeric(1))
  }
  list(internal = internal, boundary = boundary, quantile = quantile)
}
