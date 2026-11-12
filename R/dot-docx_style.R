#' Choose a Style of the Word Document
#'
#' @description
#' Internal helper of the Word report: returns the first of the candidate
#' styles that the document defines (names compared without regard to
#' case), so that the report can be written with the default document of
#' officer or with a document given by the user.
#'
#' @param styles Data frame returned by \code{officer::styles_info()}.
#' @param candidates Names of the styles, in order of preference.
#' @param type \code{"paragraph"} or \code{"table"}.
#'
#' @return The name of the style as defined in the document, or
#'   \code{NULL} (the default style of the document) if none is defined.
#'
#' @keywords internal
#' @noRd
.docx_style <- function(styles, candidates, type = "paragraph") {
  available <- styles$style_name[styles$style_type == type]
  for (candidate in candidates) {
    hit <- available[tolower(available) == tolower(candidate)]
    if (length(hit) > 0L) return(hit[1L])
  }
  NULL
}
