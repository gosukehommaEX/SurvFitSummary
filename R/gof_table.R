#' Goodness-of-Fit Statistics of the Fitted Models
#'
#' @description
#' Reports the log-likelihood, the Akaike information criterion (AIC) and
#' the Bayesian information criterion (BIC) of each fitted model, sorted by
#' AIC in ascending order within each fit (each arm of the independent
#' fits, then the dependent fit).
#'
#' @details
#' \code{AIC = -2 * loglik + 2 * npars} and
#' \code{BIC = -2 * loglik + log(n) * npars}, where \code{npars} is the
#' number of estimated parameters (the fixed shape of the Rayleigh model is
#' not counted) and \code{n} is the number of patients in the data fitted,
#' censored or not (the default of \code{flexsurv::BIC.flexsurvreg()}). The
#' log-likelihood is that of flexsurv, on the time scale (it includes the
#' log density of the event times).
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm},
#'   \code{model}, \code{n}, \code{events}, \code{npars}, \code{loglik},
#'   \code{AIC} and \code{BIC}. Models that failed are left out.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull",
#'                                                    "llogis"),
#'                             spline_scales = NULL)
#' gof_table(fits)
#'
#' @export
gof_table <- function(x) {
  .check_survfitsummary(x)
  rows <- list()
  for (i in which(!vapply(x$fits, is.null, logical(1)))) {
    mp <- .model_parameters(x$fits[[i]], x$models$distribution[i])
    npars <- length(mp$est)
    loglik <- x$fits[[i]]$loglik
    n <- x$models$n[i]
    rows[[length(rows) + 1L]] <- data.frame(
      fit_type = x$models$fit_type[i], arm = x$models$arm[i],
      model = x$models$model[i], n = n, events = x$models$events[i],
      npars = npars, loglik = loglik,
      AIC = -2 * loglik + 2 * npars,
      BIC = -2 * loglik + log(n) * npars,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(fit_type = character(0), arm = character(0),
                      model = character(0), n = integer(0),
                      events = integer(0), npars = integer(0),
                      loglik = numeric(0), AIC = numeric(0),
                      BIC = numeric(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  group <- paste(out$fit_type, out$arm, sep = "\t")
  out <- out[order(match(group, unique(group)), out$AIC), , drop = FALSE]
  rownames(out) <- NULL
  out
}
