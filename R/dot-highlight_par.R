#' Paragraph of the Word Report with a Part Highlighted
#'
#' @description
#' Internal helper of \code{.render_report()}: a paragraph whose part
#' \code{part} (for example a placeholder in square brackets for the user
#' to complete) is highlighted in yellow (black text on a yellow
#' background), as in the template of the report; the rest of the text
#' keeps the formatting of the paragraph style. officer writes a shading
#' colour only together with a text colour, so the colour of the part is
#' set to black.
#'
#' @param text Text of the paragraph.
#' @param part Part of \code{text} to highlight. When \code{text} does not
#'   contain it, nothing is highlighted.
#'
#' @return An \code{fpar} object of officer.
#'
#' @keywords internal
#' @noRd
.highlight_par <- function(text, part) {
  at <- regexpr(part, text, fixed = TRUE)
  if (at < 0L) return(officer::fpar(officer::ftext(text)))
  before <- substr(text, 1L, at - 1L)
  after <- substr(text, at + nchar(part), nchar(text))
  pieces <- list(
    if (nzchar(before)) officer::ftext(before),
    officer::ftext(part,
                   prop = officer::fp_text_lite(color = "black",
                                                shading.color = "yellow")),
    if (nzchar(after)) officer::ftext(after)
  )
  do.call(officer::fpar, pieces[!vapply(pieces, is.null, logical(1))])
}
