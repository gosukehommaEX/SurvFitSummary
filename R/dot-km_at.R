#' Kaplan-Meier Estimate at Given Times
#'
#' @description
#' Internal helper that evaluates the right-continuous Kaplan-Meier
#' estimate, its Greenwood standard error and the number at risk at the
#' times \code{x}. With \code{log_scale = TRUE}, \code{x} is on the log time
#' scale and the comparisons with the observed times are made on that
#' scale, so that a knot placed at an observed log time is counted
#' consistently (the event at that time is included in the estimate and the
#' patient in the number at risk).
#'
#' @param time Observed times.
#' @param event Event indicator (1 = event, 0 = censored).
#' @param x Times (or log times) at which the estimate is evaluated.
#' @param log_scale Logical; \code{TRUE} when \code{x} is on the log scale.
#'
#' @return A data frame with the columns \code{n_risk} (number with an
#'   observed time at or after \code{x}), \code{surv} and \code{se}.
#'
#' @keywords internal
#' @noRd
.km_at <- function(time, event, x, log_scale = FALSE) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  grid <- if (log_scale) log(fit$time) else fit$time
  observed <- if (log_scale) log(time) else time
  idx <- findInterval(x, grid) + 1L
  surv <- c(1, fit$surv)[idx]
  se <- surv * c(0, fit$std.err)[idx]
  se[!is.finite(se)] <- NA_real_
  data.frame(
    n_risk = vapply(x, function(v) sum(observed >= v), integer(1)),
    surv = surv,
    se = se
  )
}
