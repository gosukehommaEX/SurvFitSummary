#' Models with Fitting Problems in the Word Report
#'
#' @description
#' Internal helper: a table of the models fitted with warnings or that
#' could not be fitted, with their messages, at the end of the Results
#' section. Models that could not be fitted are left out of the other
#' tables and figures.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#'
#' @return A list of elements made by \code{.report_block()}; an empty list
#'   when all models were fitted without problems.
#'
#' @keywords internal
#' @noRd
.report_problems <- function(x) {
  m <- x$models[x$models$status != "ok", , drop = FALSE]
  if (nrow(m) == 0L) return(list())
  table <- data.frame(
    "Fit" = ifelse(m$fit_type == "dependent", "Both arms", "Each arm"),
    "Arm" = m$arm,
    "Model" = m$model,
    "Status" = ifelse(m$status == "error", "Not fitted",
                      "Fitted with a warning"),
    "Message" = m$message,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  list(
    .report_block("heading", "Notes on model fitting", level = 2L),
    .report_block("paragraph", paste(
      "The following models were fitted with a warning or could not be",
      "fitted. Models that could not be fitted are left out of the tables",
      "and figures."
    )),
    .report_block("table", "Models with warnings or errors", table = table)
  )
}
