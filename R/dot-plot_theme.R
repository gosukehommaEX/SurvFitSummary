#' Theme of the Plots of the Package
#'
#' @description
#' Internal helper: a classic theme (axes without grid) with a base font
#' size large enough for figures in Word reports.
#'
#' @param legend_position Position of the legend, for example
#'   \code{"right"} or \code{"bottom"}.
#' @param base_size Base font size in points.
#' @param text_sizes \code{NULL} or the sizes returned by
#'   \code{.check_text_sizes()}: those that are not \code{NA} replace the
#'   sizes of the theme (\code{"title"} sets the titles of the plot and of
#'   its panels).
#'
#' @return A ggplot2 theme.
#'
#' @keywords internal
#' @noRd
.plot_theme <- function(legend_position = "right", base_size = 14,
                        text_sizes = NULL) {
  theme <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 16, 8, 8)
    )
  # Elements of the theme set by each name of text_sizes
  elements <- c(plot.title = "title", strip.text = "title",
                axis.title = "axis_title", axis.text = "axis_text",
                legend.text = "legend")
  for (e in names(elements)) {
    size <- unname(text_sizes[elements[[e]]])
    if (length(size) == 1L && !is.na(size)) {
      theme <- theme + do.call(
        ggplot2::theme,
        stats::setNames(list(ggplot2::element_text(size = size)), e)
      )
    }
  }
  theme
}
