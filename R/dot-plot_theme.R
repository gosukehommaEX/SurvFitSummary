#' Theme of the Plots of the Package
#'
#' @description
#' Internal helper: a classic theme (axes without grid) with a base font
#' size large enough for figures in Word reports.
#'
#' @param legend_position Position of the legend, for example
#'   \code{"right"} or \code{"bottom"}.
#' @param base_size Base font size in points.
#'
#' @return A ggplot2 theme.
#'
#' @keywords internal
#' @noRd
.plot_theme <- function(legend_position = "right", base_size = 14) {
  ggplot2::theme_classic(base_size = base_size) +
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
}
