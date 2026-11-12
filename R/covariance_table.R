#' Covariance Matrices and Cholesky Factors of the Fitted Models
#'
#' @description
#' Reports, for each fitted model, the covariance matrix of the parameter
#' estimates or its Cholesky factor, on the transformed (real-line) scale on
#' which flexsurv estimates the parameters and/or on the natural scale.
#'
#' @details
#' The covariance matrix on the transformed scale is the inverse of the
#' Hessian of minus the log-likelihood computed by flexsurv. On the natural
#' scale it is \code{J V J} (delta method), with \code{J} the diagonal
#' matrix of the derivatives of the natural parameters with respect to the
#' transformed ones (the estimate itself for a log-transformed parameter, 1
#' otherwise). The Cholesky factor is the lower triangular matrix \code{L}
#' with \code{L \%*\% t(L)} equal to the covariance matrix; correlated
#' normal draws for probabilistic sensitivity analysis are \code{est + L z}
#' with \code{z} independent standard normal. Draws should be made on the
#' transformed scale, where the asymptotic normal distribution is used by
#' flexsurv and positive parameters stay positive after back-transformation.
#' The parameters are those of \code{\link{parameter_table}}.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param type \code{"covariance"} (the default) or \code{"cholesky"}.
#' @param scale Character vector with \code{"transformed"} and/or
#'   \code{"natural"}. The default is both.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm},
#'   \code{model}, \code{scale}, \code{row}, \code{col} (parameter names)
#'   and \code{value}, with one row per element of the matrix of each
#'   fitted model, row by row (the elements above the diagonal of the
#'   Cholesky factor are 0). Values are \code{NA} when the covariance matrix
#'   is not available or not positive definite.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = "weibull",
#'                             spline_scales = NULL)
#' covariance_table(fits)
#' covariance_table(fits, type = "cholesky", scale = "transformed")
#'
#' @export
covariance_table <- function(x,
                             type = c("covariance", "cholesky"),
                             scale = c("transformed", "natural")) {
  .check_survfitsummary(x)
  type <- match.arg(type)
  scale <- match.arg(scale, several.ok = TRUE)
  rows <- list()
  for (i in which(!vapply(x$fits, is.null, logical(1)))) {
    mp <- .model_parameters(x$fits[[i]], x$models$distribution[i])
    p <- length(mp$est)
    for (s in scale) {
      v <- mp$cov
      if (s == "natural") {
        jac <- ifelse(mp$log_scale, exp(mp$est), 1)
        v <- v * outer(jac, jac)
      }
      m <- v
      if (type == "cholesky") {
        m <- tryCatch(t(chol(v)), error = function(e) {
          matrix(NA_real_, p, p)
        })
      }
      rows[[length(rows) + 1L]] <- data.frame(
        fit_type = x$models$fit_type[i], arm = x$models$arm[i],
        model = x$models$model[i], scale = s,
        row = rep(names(mp$est), each = p),
        col = rep(names(mp$est), times = p),
        value = as.vector(t(m)), stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(fit_type = character(0), arm = character(0),
                      model = character(0), scale = character(0),
                      row = character(0), col = character(0),
                      value = numeric(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
