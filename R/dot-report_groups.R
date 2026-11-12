#' Groups of Fitted Models in the Report
#'
#' @description
#' Internal helper: the report presents the models of each arm fitted
#' separately and then the models of both arms fitted together, as in the
#' report template. Models that failed are left out.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param family \code{"parametric"} or \code{"spline"}.
#'
#' @return A list with one element per group: \code{fit_type}, \code{arm}
#'   (the arm, or \code{"Both arms"}), \code{models} (labels, in the order
#'   of fitting) and \code{arms} (the arms predicted).
#'
#' @keywords internal
#' @noRd
.report_groups <- function(x, family) {
  m <- x$models
  ok <- !vapply(x$fits, is.null, logical(1)) & m$family == family
  groups <- list()
  for (a in x$arms) {
    models <- unique(m$model[ok & m$fit_type == "independent" &
                               m$arm == a])
    if (length(models) > 0L) {
      groups[[length(groups) + 1L]] <- list(fit_type = "independent",
                                            arm = a, models = models,
                                            arms = a)
    }
  }
  models <- unique(m$model[ok & m$fit_type == "dependent"])
  if (length(models) > 0L) {
    groups[[length(groups) + 1L]] <- list(fit_type = "dependent",
                                          arm = "Both arms",
                                          models = models, arms = x$arms)
  }
  groups
}
