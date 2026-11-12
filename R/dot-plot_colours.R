#' Colours of the Arms and of the Models in the Plots
#'
#' @description
#' Internal helper. Arms take the colours of version 1.0.0 in the order of
#' the arms. Each parametric distribution has a fixed colour, so that a
#' distribution has the same colour in every plot; spline models take a
#' colour family by scale (hazard orange, odds purple, normal green), darker
#' with more knots.
#'
#' @param arms Character vector of arms, or \code{NULL}.
#' @param models \code{NULL}, or rows of \code{x$models} (columns
#'   \code{distribution}, \code{scale}, \code{n_knots} and \code{model}).
#'
#' @return A named character vector of colours: names are the arms, or the
#'   model labels when \code{models} is given.
#'
#' @keywords internal
#' @noRd
.plot_colours <- function(arms = NULL, models = NULL) {
  if (is.null(models)) {
    palette <- c("#004C97", "#F0B323", "#658D1B", "#A62B4E", "#D91E49",
                 "#939597")
    return(stats::setNames(rep(palette, length.out = length(arms)), arms))
  }
  parametric <- c(exp = "#984EA3", weibull = "#E41A1C",
                  weibullPH = "#99000D", rayleigh = "#FB9A99",
                  lnorm = "#A65628", llogis = "#377EB8",
                  gompertz = "#4DAF4A", gamma = "#F781BF",
                  gengamma = "#FF7F00", gengamma.orig = "#B15928",
                  genf = "#1B9E77", genf.orig = "#666666")
  spline <- list(hazard = c("#FDAE6B", "#F16913", "#A63603"),
                 odds = c("#BCBDDC", "#807DBA", "#54278F"),
                 normal = c("#A1D99B", "#41AB5D", "#006D2C"))
  colours <- vapply(seq_len(nrow(models)), function(i) {
    if (is.na(models$scale[i])) {
      unname(parametric[models$distribution[i]])
    } else {
      spline[[models$scale[i]]][min(models$n_knots[i], 3L)]
    }
  }, character(1))
  stats::setNames(colours, models$model)
}
