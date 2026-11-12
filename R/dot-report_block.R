#' One Element of the Word Report
#'
#' @description
#' Internal helper: the report is first built as a list of elements
#' (title, headings, paragraphs, tables and figures), which
#' \code{.render_report()} then writes to the Word document. Figures are
#' kept as a plotting function and its arguments, and drawn only when the
#' document is written.
#'
#' @param type \code{"title"}, \code{"subtitle"}, \code{"heading"},
#'   \code{"paragraph"}, \code{"table"} or \code{"figure"}.
#' @param text Text of the element, or the caption of a table or figure
#'   (without its number).
#' @param level Level of a heading (1 to 4).
#' @param table Data frame of character columns (tables).
#' @param plot List with \code{fun} (a plotting function of the package)
#'   and \code{args} (its arguments) (figures).
#' @param width,height Size of the figure when drawn, in inches.
#' @param display_width Width of the figure in the document, in inches (the
#'   text of an A4 page with 2.5 cm margins is 6.3 inches wide).
#'
#' @return A list with the arguments as elements.
#'
#' @keywords internal
#' @noRd
.report_block <- function(type, text = NULL, level = NULL, table = NULL,
                          plot = NULL, width = 9, height = 6,
                          display_width = 6.2) {
  list(type = type, text = text, level = level, table = table, plot = plot,
       width = width, height = height, display_width = display_width)
}
