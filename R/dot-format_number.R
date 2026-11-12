#' Format Numbers for the Tables of the Report
#'
#' @description
#' Internal helper that rounds numbers to a fixed number of decimals for
#' the Word report, without a minus sign on values rounded to 0.
#'
#' @param x Numeric vector.
#' @param digits Number of decimals.
#' @param na Text of missing values (\code{"NE"}, not estimable, by
#'   default).
#'
#' @return A character vector.
#'
#' @keywords internal
#' @noRd
.format_number <- function(x, digits, na = "NE") {
  x <- round(x, digits)
  x[!is.na(x) & x == 0] <- 0
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- na
  out
}
