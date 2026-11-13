#' Summary of the Survival Models Fitted by fit_survival_models()
#'
#' @description
#' The \code{summary()} method of the object returned by
#' \code{\link{fit_survival_models}}: the fit statistics of
#' \code{\link{gof_table}} for all the models, with their status, including
#' the models that failed (whose statistics are \code{NA}). The warnings and
#' errors themselves are printed by \code{\link{print.survfitsummary}}.
#'
#' @param object Object returned by \code{\link{fit_survival_models}}.
#' @param ... Not used.
#'
#' @return A data frame with one row per model and the columns
#'   \code{fit_type}, \code{arm}, \code{model}, \code{status}
#'   (\code{"ok"}, \code{"warning"} or \code{"error"}), \code{n},
#'   \code{events}, \code{npars}, \code{loglik}, \code{AIC} and \code{BIC},
#'   sorted by AIC within each fit as in \code{\link{gof_table}}, with the
#'   models that failed last.
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
#' summary(fits)
#'
#' @export
summary.survfitsummary <- function(object, ...) {
  if (!inherits(object, "survfitsummary")) {
    stop("'object' must be an object returned by fit_survival_models().",
         call. = FALSE)
  }
  m <- object$models
  gof <- gof_table(object)
  key <- function(d) paste(d$fit_type, d$arm, d$model, sep = "\t")
  k <- match(key(m), key(gof))
  out <- data.frame(fit_type = m$fit_type, arm = m$arm, model = m$model,
                    status = m$status, n = m$n, events = m$events,
                    npars = gof$npars[k], loglik = gof$loglik[k],
                    AIC = gof$AIC[k], BIC = gof$BIC[k],
                    stringsAsFactors = FALSE)
  # The fits in the order of the models, each sorted by AIC (failed models,
  # without AIC, last)
  group <- paste(out$fit_type, out$arm, sep = "\t")
  out <- out[order(match(group, unique(group)), out$AIC), , drop = FALSE]
  rownames(out) <- NULL
  out
}
