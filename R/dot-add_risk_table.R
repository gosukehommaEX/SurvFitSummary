#' Add a Table of Numbers at Risk Below a Plot
#'
#' @description
#' Internal helper that places the numbers at risk of each arm at the
#' breaks of the time axis below the plot (patchwork), with the same time
#' axis. A number is shown only up to the largest observed time of the arm.
#' The numbers at time 0 are centred on the left edge of the panel, so the
#' names of the arms are moved to the left by half the width of the widest
#' of these numbers (a digit is about 0.6 times the font size wide).
#'
#' @param p ggplot object with the time on the x axis from 0 to
#'   \code{max_time}.
#' @param data Data frame with the columns \code{ARM} and \code{SURVTIME}.
#' @param arms Arms, in the order of the rows of the table.
#' @param breaks Times at which the numbers are reported.
#' @param max_time Largest time of the x axis.
#' @param colours Named colours of the numbers of each arm.
#' @param text_size \code{NULL} or \code{NA} (the default, 0.8 times the
#'   base size of 14 points) or the font size in points of the numbers and
#'   of the names of the arms; the title of the table is 1.125 times as
#'   large (0.9 times the base size by default).
#'
#' @return A patchwork object: the plot above the table.
#'
#' @keywords internal
#' @noRd
.add_risk_table <- function(p, data, arms, breaks, max_time, colours,
                            text_size = NA) {
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
  if (is.null(text_size) || is.na(text_size)) text_size <- 0.8 * base_size
  first <- table_data$n[table_data$time == min(table_data$time)]
  label_gap <- 0.5 * max(nchar(first)) * 0.6 * text_size + 4
  table_plot <- ggplot2::ggplot(
    table_data,
    ggplot2::aes(x = .data$time, y = .data$arm, label = .data$n,
                 colour = .data$arm)
  ) +
    ggplot2::geom_text(size = text_size / ggplot2::.pt,
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
      axis.text.y = ggplot2::element_text(
        colour = "black", size = text_size,
        margin = ggplot2::margin(r = label_gap)
      ),
      plot.title = ggplot2::element_text(size = 1.125 * text_size),
      plot.margin = ggplot2::margin(0, 16, 4, 8)
    )
  patchwork::wrap_plots(p, table_plot, ncol = 1,
                        heights = c(1, 0.1 + 0.06 * length(arms)))
}
