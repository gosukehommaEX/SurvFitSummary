#' Text of an Item of a Bulleted List in the Word Report
#'
#' @description
#' Internal helper of \code{.render_report()}: the text of an item of a
#' bulleted list, preceded by a bullet when the document has no style
#' "List Bullet" to draw it.
#'
#' @param text Text of the item.
#' @param styled Logical; \code{TRUE} when the style "List Bullet" of the
#'   document draws the bullet.
#'
#' @return The text, preceded by a bullet and a space when \code{styled} is
#'   \code{FALSE}.
#'
#' @keywords internal
#' @noRd
.bullet_text <- function(text, styled) {
  if (styled) text else paste0("\u2022 ", text)
}
