#' Scaled Schoenfeld Residuals of the Treatment Effect
#'
#' @description
#' For each arm other than the control arm, fits a Cox model to the data of
#' the two arms and plots the scaled Schoenfeld residuals of the treatment
#' effect (\code{survival::cox.zph()}) against time with:
#' \itemize{
#'   \item the estimated log hazard ratio of the Cox model (pink line);
#'   \item the regression of the residuals on time (black line), a natural
#'     cubic spline with 4 degrees of freedom as in
#'     \code{survival::plot.cox.zph()}, with its pointwise confidence
#'     interval (dotted lines);
#'   \item the log hazard ratio and the p-value of the Grambsch-Therneau
#'     test of \code{\link{test_proportional_hazards}}, in the title of the
#'     panel.
#' }
#' A flat regression line close to the log hazard ratio indicates that the
#' treatment effect does not depend on time, that is, proportional hazards.
#'
#' @details
#' The scaled residuals are centred on the log hazard ratio (they estimate
#' the time-varying coefficient beta(t)), so the reference line is the
#' estimated log hazard ratio and not 0. The confidence interval of the
#' regression uses the variance returned by \code{cox.zph()}, as in
#' \code{plot.cox.zph()}, with the normal quantile of \code{conf_level}.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}, with
#'   two or more treatment arms.
#' @param conf_level Confidence level of the interval of the regression
#'   line, a single number in (0, 1). The default is 0.95.
#' @param time_unit Time unit of the data: \code{"months"} (the default),
#'   \code{"days"}, \code{"weeks"} or \code{"years"}.
#'
#' @return A ggplot object (one panel per comparison with the control arm,
#'   titled with the comparison, the log hazard ratio and the p-value).
#'
#' @references
#' Grambsch PM, Therneau TM (1994). Proportional hazards tests and
#' diagnostics based on weighted residuals. Biometrika 81(3):515-526.
#' \doi{10.1093/biomet/81.3.515}
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' plot_schoenfeld(dat)
#'
#' @export
plot_schoenfeld <- function(data,
                            conf_level = 0.95,
                            time_unit = c("months", "days", "weeks",
                                          "years")) {
  .check_survival_data(data, min_arms = 2L)
  time_unit <- match.arg(time_unit)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  data$ARM <- droplevels(data$ARM)
  arms <- levels(data$ARM)
  tests <- test_proportional_hazards(data)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  residuals <- list()
  smooth <- list()
  reference <- list()
  for (k in seq_along(arms)[-1L]) {
    d <- data[data$ARM %in% arms[c(1L, k)], , drop = FALSE]
    d$ARM <- droplevels(d$ARM)
    fit <- survival::coxph(survival::Surv(SURVTIME, EVENT) ~ ARM, data = d)
    zph <- survival::cox.zph(fit, transform = "identity")
    comparison <- paste(arms[k], "vs", arms[1L])
    time <- unname(zph$time)
    y <- unname(zph$y[, 1L])
    # Natural spline regression on time, as in survival::plot.cox.zph()
    grid <- seq(min(time), max(time), length.out = 100L)
    basis <- splines::ns(c(grid, time), df = 4, intercept = TRUE)
    pmat <- basis[seq_along(grid), , drop = FALSE]
    xmat <- basis[-seq_along(grid), , drop = FALSE]
    qmat <- qr(xmat)
    fitted <- drop(pmat %*% qr.coef(qmat, y))
    bk <- backsolve(qmat$qr[1:4, 1:4], diag(4))
    seval <- drop(((pmat %*% (bk %*% t(bk))) * pmat) %*% rep(1, 4))
    half <- z * sqrt(zph$var[1L, 1L] * seval)
    test <- tests[tests$comparison == comparison, ]
    residuals[[comparison]] <- data.frame(comparison = comparison,
                                          time = time, y = y)
    smooth[[comparison]] <- data.frame(comparison = comparison, time = grid,
                                       fit = fitted, lcl = fitted - half,
                                       ucl = fitted + half)
    p_text <- if (test$p_value < 0.001) {
      "p < 0.001"
    } else {
      sprintf("p = %.3f", test$p_value)
    }
    reference[[comparison]] <- data.frame(
      comparison = comparison, log_hr = unname(stats::coef(fit)[1L]),
      label = sprintf("%s (log HR = %.3f; Grambsch-Therneau test %s)",
                      comparison, test$log_hr, p_text)
    )
  }
  residuals <- do.call(rbind, residuals)
  smooth <- do.call(rbind, smooth)
  reference <- do.call(rbind, reference)
  comparisons <- paste(arms[-1L], "vs", arms[1L])
  residuals$comparison <- factor(residuals$comparison, levels = comparisons)
  smooth$comparison <- factor(smooth$comparison, levels = comparisons)
  reference$comparison <- factor(reference$comparison, levels = comparisons)
  # Panel titles with the log hazard ratio and the p-value of the test
  panels <- reference$label[match(comparisons, reference$comparison)]
  residuals$panel <- factor(panels[match(residuals$comparison, comparisons)],
                            levels = panels)
  smooth$panel <- factor(panels[match(smooth$comparison, comparisons)],
                         levels = panels)
  reference$panel <- factor(reference$label, levels = panels)
  axis <- .time_axis(max(residuals$time), time_unit)

  lines <- c("Estimated log hazard ratio", "Regression on time",
             paste0(round(100 * conf_level), "% confidence interval"))
  colours <- stats::setNames(c("#E7298A", "black", "black"), lines)
  types <- stats::setNames(c("solid", "solid", "dotted"), lines)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = residuals,
      ggplot2::aes(x = .data$time, y = .data$y,
                   shape = "Scaled Schoenfeld residual"),
      colour = "#004C97", alpha = 0.5, size = 1.8
    ) +
    ggplot2::geom_hline(
      data = reference,
      ggplot2::aes(yintercept = .data$log_hr, colour = lines[1L],
                   linetype = lines[1L]),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = smooth,
      ggplot2::aes(x = .data$time, y = .data$fit, colour = lines[2L],
                   linetype = lines[2L]),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = smooth,
      ggplot2::aes(x = .data$time, y = .data$lcl, colour = lines[3L],
                   linetype = lines[3L]),
      linewidth = 0.6
    ) +
    ggplot2::geom_line(
      data = smooth,
      ggplot2::aes(x = .data$time, y = .data$ucl, colour = lines[3L],
                   linetype = lines[3L]),
      linewidth = 0.6
    ) +
    ggplot2::scale_colour_manual(values = colours, breaks = lines) +
    ggplot2::scale_linetype_manual(values = types, breaks = lines) +
    ggplot2::scale_shape_manual(values = c("Scaled Schoenfeld residual" = 16)) +
    ggplot2::scale_x_continuous(breaks = axis$breaks) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_x") +
    ggplot2::labs(x = axis$label,
                  y = "Scaled Schoenfeld residual, beta(t)") +
    .plot_theme(legend_position = "bottom") +
    # The legend of the points above that of the lines, so that the legend
    # is not wider than the plot
    ggplot2::theme(legend.box = "vertical") +
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1),
                    colour = ggplot2::guide_legend(order = 2),
                    linetype = ggplot2::guide_legend(order = 2))
  p
}
