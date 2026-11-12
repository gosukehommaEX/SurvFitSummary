#' Check a Dataset Returned by processing_dataset()
#'
#' @description
#' Internal helper that stops with an informative message when \code{data}
#' is not a data frame with the columns \code{ARM} (factor), \code{SURVTIME}
#' and \code{EVENT}, or has fewer treatment arms than \code{min_arms}.
#'
#' @param data Object to check.
#' @param min_arms Smallest number of treatment arms that is allowed.
#'
#' @return \code{TRUE}, invisibly.
#'
#' @keywords internal
#' @noRd
.check_survival_data <- function(data, min_arms = 1L) {
  if (!is.data.frame(data) ||
      !all(c("ARM", "SURVTIME", "EVENT") %in% names(data))) {
    stop("'data' must be a data frame with the columns ARM, SURVTIME and ",
         "EVENT, as returned by processing_dataset().", call. = FALSE)
  }
  if (!is.factor(data$ARM)) {
    stop("The column ARM of 'data' must be a factor; prepare the data ",
         "with processing_dataset().", call. = FALSE)
  }
  if (nlevels(droplevels(data$ARM)) < min_arms) {
    stop("At least ", min_arms, " treatment arms are needed.", call. = FALSE)
  }
  invisible(TRUE)
}
