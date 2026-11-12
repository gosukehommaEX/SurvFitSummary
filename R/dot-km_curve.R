#' Kaplan-Meier Curves for Plotting
#'
#' @description
#' Internal helper that returns the Kaplan-Meier estimate of each arm at
#' time 0 and at every observed time (\code{survival::survfit()}), so that
#' each curve ends at the largest observed time of its arm, with pointwise
#' confidence limits (log transformation, the default of \code{survfit()}).
#'
#' @param data Data frame with the columns \code{ARM}, \code{SURVTIME} and
#'   \code{EVENT}.
#' @param arms Arms to include, in this order.
#' @param conf_level Confidence level.
#'
#' @return A data frame with the columns \code{arm}, \code{time},
#'   \code{surv}, \code{lcl} and \code{ucl}.
#'
#' @keywords internal
#' @noRd
.km_curve <- function(data, arms, conf_level = 0.95) {
  rows <- lapply(arms, function(a) {
    d <- data[data$ARM == a, , drop = FALSE]
    fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1, data = d,
                             conf.int = conf_level)
    data.frame(arm = a, time = c(0, fit$time), surv = c(1, fit$surv),
               lcl = c(1, fit$lower), ucl = c(1, fit$upper),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$lcl <- pmax(out$lcl, 0)
  out$ucl <- pmin(out$ucl, 1)
  rownames(out) <- NULL
  out
}
