#' Print the Survival Models Fitted by fit_survival_models()
#'
#' @description
#' Prints the time unit, the fits (arms), the numbers of models fitted
#' without problems, with warnings and failed, and the warnings and errors
#' of the models that had any.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = "exp",
#'                             spline_scales = NULL)
#' print(fits)
#'
#' @export
print.survfitsummary <- function(x, ...) {
  m <- x$models
  counts <- table(factor(m$status, levels = c("ok", "warning", "error")))
  fitted <- unique(paste0(m$arm, " (", m$fit_type, ")"))
  writeLines(c(
    paste0("Survival models fitted with SurvFitSummary (time unit: ",
           x$time_unit, ")"),
    paste0("Fits: ", paste(fitted, collapse = ", ")),
    paste0(nrow(m), " models: ", counts[["ok"]], " without problems, ",
           counts[["warning"]], " with warnings, ", counts[["error"]],
           " failed")
  ))
  problems <- m[m$status != "ok", c("fit_type", "arm", "model", "status",
                                    "message"), drop = FALSE]
  if (nrow(problems) > 0L) {
    writeLines("Models with warnings or errors:")
    print(problems, row.names = FALSE)
  }
  invisible(x)
}
