#' Arms Predicted by a Fitted Model
#'
#' @description
#' Internal helper that returns the treatment arms for which a fitted model
#' gives predictions: the arm of an independent fit, or every arm of a
#' dependent fit, with the \code{newdata} to pass to
#' \code{summary.flexsurvreg()}.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param i Index of the model in \code{x$models}.
#'
#' @return A list of lists with \code{arm} (character) and \code{newdata}
#'   (\code{NULL} for an independent fit, otherwise a data frame with the
#'   factor \code{ARM}).
#'
#' @keywords internal
#' @noRd
.model_arms <- function(x, i) {
  if (x$models$fit_type[i] == "dependent") {
    lapply(x$arms, function(a) {
      list(arm = a, newdata = data.frame(ARM = factor(a, levels = x$arms)))
    })
  } else {
    list(list(arm = x$models$arm[i], newdata = NULL))
  }
}
