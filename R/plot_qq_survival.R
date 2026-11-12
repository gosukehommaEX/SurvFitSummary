#' Quantile-Quantile Plot of Survival Times
#'
#' @description
#' For each arm other than the control arm, plots the quantiles of the
#' survival time of the arm against those of the control arm, both from the
#' Kaplan-Meier estimates (\code{quantile()} of \code{survival::survfit()}),
#' with the identity line (dashed) and the least-squares line through the
#' points (solid). Points on a straight line through the origin indicate
#' that the survival times of the arm are a multiple of those of the control
#' arm (an accelerated failure time effect); points on the identity line
#' indicate no difference. The axes go from 0 to the largest quantile,
#' with the same scale on both axes.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}, with
#'   two or more treatment arms.
#' @param probs Probabilities of the quantiles, numbers in (0, 1). The
#'   default is \code{seq(0.05, 0.95, by = 0.05)}. Quantiles that are not
#'   reached by a Kaplan-Meier curve are left out.
#' @param time_unit Time unit of the data: \code{"months"} (the default),
#'   \code{"days"}, \code{"weeks"} or \code{"years"}.
#'
#' @return A ggplot object (one panel per comparison with the control arm).
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' plot_qq_survival(dat)
#'
#' @export
plot_qq_survival <- function(data,
                             probs = seq(0.05, 0.95, by = 0.05),
                             time_unit = c("months", "days", "weeks",
                                           "years")) {
  .check_survival_data(data, min_arms = 2L)
  time_unit <- match.arg(time_unit)
  if (!is.numeric(probs) || length(probs) == 0L || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1)) {
    stop("'probs' must contain numbers in (0, 1).", call. = FALSE)
  }
  data$ARM <- droplevels(data$ARM)
  arms <- levels(data$ARM)
  quantiles_of <- function(a) {
    fit <- survival::survfit(survival::Surv(SURVTIME, EVENT) ~ 1,
                             data = data[data$ARM == a, , drop = FALSE])
    q <- stats::quantile(fit, probs = probs, conf.int = FALSE)
    if (is.list(q)) q <- q$quantile
    unname(as.vector(q))
  }
  control <- quantiles_of(arms[1L])
  comparisons <- paste(arms[-1L], "vs", arms[1L])
  points <- do.call(rbind, lapply(seq_along(arms)[-1L], function(k) {
    out <- data.frame(comparison = comparisons[k - 1L], prob = probs,
                      control = control, arm = quantiles_of(arms[k]),
                      stringsAsFactors = FALSE)
    out[!is.na(out$control) & !is.na(out$arm), , drop = FALSE]
  }))
  if (nrow(points) == 0L) {
    stop("No quantile is reached by the Kaplan-Meier curves of both arms.",
         call. = FALSE)
  }
  lines <- do.call(rbind, lapply(comparisons, function(cmp) {
    d <- points[points$comparison == cmp, , drop = FALSE]
    coefs <- if (nrow(d) >= 2L) {
      unname(stats::coef(stats::lm(arm ~ control, data = d)))
    } else {
      c(NA_real_, NA_real_)
    }
    data.frame(comparison = cmp,
               type = c("Identity line (y = x)", "Least-squares line"),
               intercept = c(0, coefs[1L]), slope = c(1, coefs[2L]),
               stringsAsFactors = FALSE)
  }))
  lines <- lines[!is.na(lines$slope), , drop = FALSE]
  points$comparison <- factor(points$comparison, levels = comparisons)
  lines$comparison <- factor(lines$comparison, levels = comparisons)
  limit <- 1.05 * max(c(points$control, points$arm))
  unit_label <- substring(.time_axis(1, time_unit)$label, 6L)
  point_label <- paste0("Quantiles (", format(100 * min(probs)), "% to ",
                        format(100 * max(probs)), "%)")
  p <- ggplot2::ggplot() +
    ggplot2::geom_abline(
      data = lines,
      ggplot2::aes(intercept = .data$intercept, slope = .data$slope,
                   linetype = .data$type),
      colour = "grey30", linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = points,
      ggplot2::aes(x = .data$control, y = .data$arm, shape = point_label),
      colour = "#D91E49", size = 2.5
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Identity line (y = x)" = "dashed",
                 "Least-squares line" = "solid")
    ) +
    ggplot2::scale_shape_manual(values = stats::setNames(16, point_label)) +
    ggplot2::coord_fixed(xlim = c(0, limit), ylim = c(0, limit),
                         expand = FALSE) +
    ggplot2::labs(
      x = paste("Quantiles of", arms[1L], unit_label),
      y = if (length(comparisons) == 1L) {
        paste("Quantiles of", arms[2L], unit_label)
      } else {
        paste("Quantiles of the arm", unit_label)
      }
    ) +
    .plot_theme(legend_position = "bottom")
  if (length(comparisons) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$comparison))
  }
  p
}
