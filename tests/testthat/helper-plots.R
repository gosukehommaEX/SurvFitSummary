# Access to the layers of the plots in the tests: the data given to a layer
# (with the columns of the package, for example the model labels) and the
# ranges of the axes of the built plot.

layer_index <- function(p, geom) {
  which(vapply(p$layers, function(l) class(l$geom)[1L], character(1)) == geom)
}

layer_input <- function(p, geom, k = 1L) {
  p$layers[[layer_index(p, geom)[k]]]$data
}

panel_ranges <- function(p) {
  params <- ggplot2::ggplot_build(p)$layout$panel_params[[1L]]
  list(x = params$x.range, y = params$y.range)
}

# Data of a single arm with event times drawn from a distribution, without
# censoring, for the checks of the diagnostic plots
single_arm <- function(time) {
  data.frame(ARM = factor(rep("A", length(time))), SURVTIME = time,
             EVENT = rep(1L, length(time)))
}

# Slope of the regression line drawn in a diagnostic plot
drawn_slope <- function(p) {
  line <- layer_input(p, "GeomLine", length(layer_index(p, "GeomLine")))
  diff(line$y) / diff(line$x)
}
