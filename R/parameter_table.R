#' Parameter Estimates of the Fitted Models
#'
#' @description
#' Reports the maximum likelihood estimates of the parameters of each
#' fitted model with their standard errors and Wald confidence intervals,
#' on the natural scale of the distribution and/or on the transformed
#' (real-line) scale on which flexsurv estimates them.
#'
#' @details
#' On the transformed scale, the positive parameters (for example the shape
#' and scale of the Weibull, see \code{\link{fit_survival_models}}) are on
#' the log scale and the others, the arm effects and all spline parameters
#' are untransformed; the confidence interval is
#' \code{est +/- z * se} with \code{z} the normal quantile of the confidence
#' level of the fit. On the natural scale, the estimate and the confidence
#' limits of a log-transformed parameter are the exponentials of those on
#' the transformed scale (as in \code{flexsurvreg()$res}) and the standard
#' error is obtained by the delta method (\code{est * se} of the log). Arm
#' effects in dependent fits are effects on the transformed location
#' parameter (for example the log of the scale of the Weibull) on both
#' scales. The shape of the Rayleigh model, fixed at 2, is not reported.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param scale Character vector with \code{"natural"} and/or
#'   \code{"transformed"}. The default is both.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm},
#'   \code{model}, \code{parameter} (name used by flexsurv, for example
#'   \code{"shape"}, \code{"gamma1"} or \code{"ARMTreatment"}),
#'   \code{scale}, \code{est}, \code{lcl}, \code{ucl} and \code{se}, with
#'   one row per parameter and scale of each fitted model (models that
#'   failed are left out).
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = NULL)
#' parameter_table(fits)
#'
#' @export
parameter_table <- function(x, scale = c("natural", "transformed")) {
  .check_survfitsummary(x)
  scale <- match.arg(scale, several.ok = TRUE)
  z <- stats::qnorm(1 - (1 - x$conf_level) / 2)
  rows <- list()
  for (i in which(!vapply(x$fits, is.null, logical(1)))) {
    mp <- .model_parameters(x$fits[[i]], x$models$distribution[i])
    se_t <- sqrt(diag(mp$cov))
    lcl_t <- mp$est - z * se_t
    ucl_t <- mp$est + z * se_t
    for (s in scale) {
      if (s == "natural") {
        est <- ifelse(mp$log_scale, exp(mp$est), mp$est)
        lcl <- ifelse(mp$log_scale, exp(lcl_t), lcl_t)
        ucl <- ifelse(mp$log_scale, exp(ucl_t), ucl_t)
        se <- ifelse(mp$log_scale, est * se_t, se_t)
      } else {
        est <- mp$est
        lcl <- lcl_t
        ucl <- ucl_t
        se <- se_t
      }
      rows[[length(rows) + 1L]] <- data.frame(
        fit_type = x$models$fit_type[i], arm = x$models$arm[i],
        model = x$models$model[i], parameter = names(mp$est), scale = s,
        est = unname(est), lcl = unname(lcl), ucl = unname(ucl),
        se = unname(se), stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(fit_type = character(0), arm = character(0),
                      model = character(0), parameter = character(0),
                      scale = character(0), est = numeric(0),
                      lcl = numeric(0), ucl = numeric(0), se = numeric(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
