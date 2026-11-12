#' Step-Shaped Confidence Band of Kaplan-Meier Curves
#'
#' @description
#' Internal helper that duplicates the points of the confidence limits so
#' that \code{geom_ribbon()} draws them as step functions, like the curves
#' drawn by \code{geom_step()}. Rows with missing limits (after the curve
#' reaches 0) are left out.
#'
#' @param curve Data frame returned by \code{.km_curve()}.
#'
#' @return A data frame with the columns \code{arm}, \code{time},
#'   \code{lcl} and \code{ucl}.
#'
#' @keywords internal
#' @noRd
.step_ribbon <- function(curve) {
  pieces <- lapply(unique(curve$arm), function(a) {
    d <- curve[curve$arm == a & !is.na(curve$lcl) & !is.na(curve$ucl), ,
               drop = FALSE]
    n <- nrow(d)
    if (n < 2L) return(NULL)
    data.frame(arm = a,
               time = rep(d$time, each = 2L)[-1L],
               lcl = rep(d$lcl, each = 2L)[-(2L * n)],
               ucl = rep(d$ucl, each = 2L)[-(2L * n)],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}
