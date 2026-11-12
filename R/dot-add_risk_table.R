#' Add a Table of Numbers at Risk Below a Plot
#'
#' @description
#' Internal helper that places the numbers at risk of each arm at the
#' breaks of the time axis below the plot (patchwork), with the same time
#' axis. A number is shown only up to the largest observed time of the arm.
#'
#' @param p ggplot object with the time on the x axis from 0 to
#'   \code{max_time}.
#' @param data Data frame with the columns \code{ARM} and \code{SURVTIME}.
#' @param arms Arms, in the order of the rows of the table.
#' @param breaks Times at which the numbers are reported.
#' @param max_time Largest time of the x axis.
#' @param colours Named colours of the arms.
#'
#' @return A patchwork object: the plot above the table.
#'
#' @keywords internal
#' @noRd
.add_risk_table <- function(p, data, arms, breaks, max_time, colours) {
  rows <- lapply(arms, function(a) {
    t <- data$SURVTIME[data$ARM == a]
    shown <- breaks[breaks <= max(t)]
    data.frame(arm = a, time = shown,
               n = vapply(shown, function(b) sum(t >= b), integer(1)),
               stringsAsFactors = FALSE)
  })
  table_data <- do.call(rbind, rows)
  table_data$arm <- factor(table_data$arm, levels = rev(arms))
  base_size <- 14
  table_plot <- ggplot2::ggplot(
    table_data,
    ggplot2::aes(x = .data$time, y = .data$arm, label = .data$n,
                 colour = .data$arm)
  ) +
    ggplot2::geom_text(size = 0.8 * base_size / ggplot2::.pt,
                       show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_x_continuous(breaks = breaks,
                                expand = ggplot2::expansion(0, 0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::coord_cartesian(xlim = c(0, max_time), clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, title = "N at risk") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 0.9 * base_size),
      plot.margin = ggplot2::margin(0, 16, 4, 8)
    )
  patchwork::wrap_plots(p, table_plot, ncol = 1,
                        heights = c(1, 0.1 + 0.06 * length(arms)))
}
