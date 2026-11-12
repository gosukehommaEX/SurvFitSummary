#' Predicted Survival, Hazard or Cumulative Hazard of the Fitted Models
#'
#' @description
#' Predicts the survival probability, the hazard or the cumulative hazard
#' of each fitted model at the given times, for each arm (the arm of an
#' independent fit, every arm for a dependent fit), with the Kaplan-Meier
#' estimate of the arm for comparison. Used for the landmark survival
#' probabilities (for example every 4 weeks up to the time horizon) and for
#' the curves of the plots.
#'
#' @details
#' The predictions are those of \code{summary.flexsurvreg()}. The
#' confidence limits, when requested, are the quantiles of the predictions
#' for \code{n_sim} parameter vectors drawn from the asymptotic normal
#' distribution of the estimates on the transformed scale (flexsurv's
#' \code{normboot.flexsurvreg()}); the random number generator is seeded
#' with \code{seed} for each model and arm, so that the limits do not
#' depend on the other models fitted, and the state of the generator of the
#' session is restored afterwards. At time 0 the survival is set to 1, the
#' cumulative hazard to 0 and the hazard to \code{NA} (it is infinite or
#' undefined for some models), without calling flexsurv.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param times Non-negative times in the time unit of the data.
#' @param type \code{"survival"} (the default), \code{"hazard"} or
#'   \code{"cumhaz"}.
#' @param conf_int Logical; \code{TRUE} to add simulated confidence limits.
#'   The default is \code{FALSE}.
#' @param n_sim Number of simulated parameter vectors for the confidence
#'   limits. The default is 1000.
#' @param seed Seed of the simulation, a single whole number. The default
#'   is 1.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm} (the
#'   arm predicted), \code{model}, \code{time}, \code{est}, \code{lcl} and
#'   \code{ucl} (\code{NA} without \code{conf_int}) and \code{km} (for
#'   \code{type = "survival"}, the Kaplan-Meier estimate of the arm, and
#'   \code{NA} after its largest observed time; \code{NA} for the other
#'   types).
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = NULL)
#' predict_survival(fits, times = seq(0, 60, by = 12))
#'
#' @export
predict_survival <- function(x,
                             times,
                             type = c("survival", "hazard", "cumhaz"),
                             conf_int = FALSE,
                             n_sim = 1000,
                             seed = 1) {
  .check_survfitsummary(x)
  type <- match.arg(type)
  if (!is.numeric(times) || length(times) == 0L || any(!is.finite(times)) ||
      any(times < 0)) {
    stop("'times' must be a vector of non-negative numbers.", call. = FALSE)
  }
  if (!is.logical(conf_int) || length(conf_int) != 1L || is.na(conf_int)) {
    stop("'conf_int' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) ||
      n_sim < 2 || n_sim != round(n_sim)) {
    stop("'n_sim' must be a single whole number of at least 2.",
         call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed)) {
    stop("'seed' must be a single whole number.", call. = FALSE)
  }
  positive <- times > 0
  value_at_zero <- switch(type, survival = 1, cumhaz = 0, hazard = NA_real_)
  rows <- list()
  for (i in which(!vapply(x$fits, is.null, logical(1)))) {
    for (pa in .model_arms(x, i)) {
      pred <- list(est = rep(value_at_zero, length(times)),
                   lcl = rep(NA_real_, length(times)),
                   ucl = rep(NA_real_, length(times)))
      if (conf_int) {
        pred$lcl[!positive] <- value_at_zero
        pred$ucl[!positive] <- value_at_zero
      }
      if (any(positive)) {
        fitted <- .summary_one(x$fits[[i]], pa$newdata, type = type,
                               t = times[positive], conf_int = conf_int,
                               n_sim = n_sim, seed = seed,
                               conf_level = x$conf_level)
        pred$est[positive] <- fitted$est
        pred$lcl[positive] <- fitted$lcl
        pred$ucl[positive] <- fitted$ucl
      }
      km <- rep(NA_real_, length(times))
      if (type == "survival") {
        d <- x$data[x$data$ARM == pa$arm, , drop = FALSE]
        km <- .km_at(d$SURVTIME, d$EVENT, times)$surv
        km[times > max(d$SURVTIME)] <- NA_real_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        fit_type = x$models$fit_type[i], arm = pa$arm,
        model = x$models$model[i], time = times, est = pred$est,
        lcl = pred$lcl, ucl = pred$ucl, km = km, stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(fit_type = character(0), arm = character(0),
                      model = character(0), time = numeric(0),
                      est = numeric(0), lcl = numeric(0), ucl = numeric(0),
                      km = numeric(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
