#' Select Fitted Models and Arms for a Plot
#'
#' @description
#' Internal helper of \code{\link{plot_survival_fits}} and
#' \code{\link{plot_hazard_fits}}. Keeps the fitted models of one type
#' (independent or dependent) whose labels are in \code{models}, for the
#' arms in \code{arms}, and checks that the plot shows either one arm with
#' one or more models or one model with one or more arms.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param models \code{NULL} (all fitted models of the type) or model
#'   labels.
#' @param arms \code{NULL} (all arms) or arms.
#' @param fit_type \code{"independent"} or \code{"dependent"}.
#'
#' @return A list with \code{x} (the object reduced to the selected fits),
#'   \code{models} (labels, in the order of \code{x$models}), \code{rows}
#'   (one row of \code{x$models} per label) and \code{arms}.
#'
#' @keywords internal
#' @noRd
.select_models <- function(x, models, arms, fit_type) {
  if (is.null(arms)) arms <- x$arms
  if (!is.character(arms) || length(arms) == 0L || !all(arms %in% x$arms)) {
    stop("'arms' must be NULL or a subset of: ",
         paste(x$arms, collapse = ", "), call. = FALSE)
  }
  arms <- x$arms[x$arms %in% arms]
  keep <- x$models$fit_type == fit_type &
    !vapply(x$fits, is.null, logical(1))
  if (fit_type == "independent") keep <- keep & x$models$arm %in% arms
  available <- unique(x$models$model[keep])
  if (length(available) == 0L) {
    stop("No ", fit_type, " model was fitted for the selected arms.",
         call. = FALSE)
  }
  if (is.null(models)) models <- available
  if (!is.character(models) || length(models) == 0L ||
      !all(models %in% available)) {
    stop("'models' must be NULL or a subset of: ",
         paste(available, collapse = "; "), call. = FALSE)
  }
  models <- available[available %in% models]
  if (length(arms) > 1L && length(models) > 1L) {
    stop("Plot either one arm with one or more models, or one model with ",
         "one or more arms.", call. = FALSE)
  }
  keep <- keep & x$models$model %in% models
  y <- x
  y$models <- x$models[keep, , drop = FALSE]
  y$fits <- x$fits[keep]
  rownames(y$models) <- NULL
  list(x = y, models = models,
       rows = y$models[match(models, y$models$model), , drop = FALSE],
       arms = arms)
}
