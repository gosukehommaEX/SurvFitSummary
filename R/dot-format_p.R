#' Format P-Values for the Report
#'
#' @description
#' Internal helper: p-values with 3 decimals, and \code{"<0.001"} below
#' 0.001.
#'
#' @param p Numeric vector of p-values.
#'
#' @return A character vector.
#'
#' @keywords internal
#' @noRd
.format_p <- function(p) {
  out <- .format_number(p, 3L)
  out[!is.na(p) & p < 0.001] <- "<0.001"
  out
}
