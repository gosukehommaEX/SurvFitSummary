#' Fit One Parametric or Spline Model
#'
#' @description
#' Internal helper of \code{\link{fit_survival_models}}. Fits one model with
#' \code{flexsurv::flexsurvreg()} or \code{flexsurv::flexsurvspline()},
#' collects the warnings instead of printing them, and returns the error
#' message instead of stopping when the fit fails. The Rayleigh model is the
#' Weibull model with the shape fixed at 2 (\code{fixedpars = 1}); its
#' starting value of the scale is the maximum likelihood estimate without
#' covariates, \code{sqrt(sum(time^2) / events)}.
#'
#' flexsurv has no analytic gradient for the spline models on the normal
#' scale, so \code{optim()} uses central differences with the step
#' \code{ndeps}. Unless \code{control} sets \code{ndeps}, the step is
#' \code{1e-5} for these models instead of the default \code{1e-3} of
#' \code{optim()}, with which the fits of
#' \code{inst/extdata/fit_example.csv} stopped 2e-5 to 6e-5 short of the
#' maximum of the log-likelihood (2026-11-13, found by
#' \code{inst/validation/validate_splines_rstpm2.R}).
#'
#' @param data Data frame with the columns \code{ARM}, \code{SURVTIME} and
#'   \code{EVENT}.
#' @param dependent Logical; \code{TRUE} puts the arm in the model as a
#'   covariate of the location parameter (\code{gamma0} for splines).
#' @param distribution Distribution code of \code{.distribution_info()}, or
#'   \code{NA} for a spline model.
#' @param scale \code{NULL} for a parametric model, or the spline scale
#'   (\code{"hazard"}, \code{"odds"} or \code{"normal"}).
#' @param knots \code{NULL}, or the list returned by \code{.spline_knots()}.
#' @param conf_level Confidence level of the parameter estimates.
#' @param control List passed to \code{stats::optim()} (with \code{ndeps}
#'   added for the normal scale, see Description).
#'
#' @return A list with \code{fit} (the flexsurvreg object, or \code{NULL}
#'   if the fit failed), \code{status} (\code{"ok"}, \code{"warning"} or
#'   \code{"error"}) and \code{message} (warnings or error, \code{NA} if
#'   none).
#'
#' @keywords internal
#' @noRd
.fit_one_model <- function(data, dependent, distribution, scale = NULL,
                           knots = NULL, conf_level, control) {
  model_formula <- if (dependent) {
    survival::Surv(SURVTIME, EVENT) ~ ARM
  } else {
    survival::Surv(SURVTIME, EVENT) ~ 1
  }
  if (identical(scale, "normal") && is.null(control$ndeps)) {
    # One step per parameter: gamma0, ..., gamma(k + 1) and the arm effects
    n_cov <- if (dependent) nlevels(droplevels(data$ARM)) - 1L else 0L
    control$ndeps <- rep(1e-5, length(knots$internal) + 2L + n_cov)
  }
  collected <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      if (!is.null(scale)) {
        flexsurv::flexsurvspline(model_formula, data = data,
                                 knots = knots$internal,
                                 bknots = knots$boundary, scale = scale,
                                 cl = conf_level, control = control)
      } else if (distribution == "rayleigh") {
        n_cov <- if (dependent) nlevels(droplevels(data$ARM)) - 1L else 0L
        inits <- c(2, sqrt(sum(data$SURVTIME^2) / sum(data$EVENT)),
                   rep(0, n_cov))
        flexsurv::flexsurvreg(model_formula, data = data, dist = "weibull",
                              inits = inits, fixedpars = 1L,
                              cl = conf_level, control = control)
      } else {
        flexsurv::flexsurvreg(model_formula, data = data,
                              dist = distribution, cl = conf_level,
                              control = control)
      },
      warning = function(w) {
        collected <<- c(collected, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(fit = NULL, status = "error",
                message = conditionMessage(fit)))
  }
  notes <- unique(collected)
  convergence <- fit$opt$convergence
  if (!is.null(convergence) && convergence != 0) {
    notes <- c(notes, paste0("optim() did not report convergence (code ",
                             convergence, ")."))
  }
  # Only the estimated parameters (not the fixed Rayleigh shape)
  if (anyNA(.model_parameters(fit, distribution)$cov)) {
    notes <- c(notes, "The covariance matrix is not available.")
  }
  list(fit = fit,
       status = if (length(notes) > 0L) "warning" else "ok",
       message = if (length(notes) > 0L) {
         paste(notes, collapse = " ")
       } else {
         NA_character_
       })
}
