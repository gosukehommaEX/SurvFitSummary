#' Check the Text Sizes Given to a Plot
#'
#' @description
#' Internal helper of the plots: checks \code{text_sizes} and returns the
#' sizes of the five kinds of text of the plots, \code{NA} for those not
#' given, which keep their default sizes.
#'
#' @param text_sizes \code{NULL} or a named vector of font sizes in points,
#'   with some or all of the names \code{"title"}, \code{"axis_title"},
#'   \code{"axis_text"}, \code{"legend"} and \code{"risk_table"}.
#'
#' @return A numeric vector with the five names above.
#'
#' @keywords internal
#' @noRd
.check_text_sizes <- function(text_sizes) {
  texts <- c("title", "axis_title", "axis_text", "legend", "risk_table")
  out <- stats::setNames(rep(NA_real_, length(texts)), texts)
  if (is.null(text_sizes)) return(out)
  given <- names(text_sizes)
  if (!is.numeric(text_sizes) || length(text_sizes) == 0L ||
      is.null(given) || anyNA(given) || any(given == "") ||
      anyDuplicated(given) > 0L || !all(given %in% texts) ||
      any(!is.finite(text_sizes)) || any(text_sizes <= 0)) {
    stop("'text_sizes' must be NULL or a named vector of positive numbers ",
         "with distinct names among: ", paste(texts, collapse = ", "),
         call. = FALSE)
  }
  out[given] <- unname(text_sizes)
  out
}
