#' Join Words into an English List
#'
#' @description
#' Internal helper of the report text: \code{c("a", "b", "c")} becomes
#' \code{"a, b and c"}.
#'
#' @param words Character vector.
#'
#' @return A single string.
#'
#' @keywords internal
#' @noRd
.and_list <- function(words) {
  n <- length(words)
  if (n <= 1L) return(paste(words, collapse = ""))
  paste(paste(words[-n], collapse = ", "), "and", words[n])
}
