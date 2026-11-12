#' Check an Object Returned by fit_survival_models()
#'
#' @description
#' Internal helper that stops with an informative message when \code{x} is
#' not of class \code{survfitsummary}.
#'
#' @param x Object to check.
#'
#' @return \code{TRUE}, invisibly.
#'
#' @keywords internal
#' @noRd
.check_survfitsummary <- function(x) {
  if (!inherits(x, "survfitsummary")) {
    stop("'x' must be an object returned by fit_survival_models().",
         call. = FALSE)
  }
  invisible(TRUE)
}
