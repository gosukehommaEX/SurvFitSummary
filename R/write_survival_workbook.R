#' Write the Results to an Excel Workbook
#'
#' @description
#' Writes the tables of the analysis to an Excel workbook (with the writexl
#' package), one table per sheet with bold column names, with unrounded
#' numbers (empty cells for missing values) and the print-ready names of
#' the models (for example \code{"Log-logistic"} or
#' \code{"Hazards, 2 knots"}). The sheets are:
#' \describe{
#'   \item{Settings}{Time unit, confidence level, reporting times and the
#'     versions of R and of the packages.}
#'   \item{Data summary}{\code{\link{summarise_data}}.}
#'   \item{KM table}{\code{\link{km_table}} at \code{interval}.}
#'   \item{Hazard ratios, PH test}{\code{hazard_ratios} and
#'     \code{\link{test_proportional_hazards}} (two or more arms).}
#'   \item{Model status}{Each model with its status (\code{"ok"},
#'     \code{"warning"} or \code{"error"}) and message.}
#'   \item{Fit statistics}{\code{\link{gof_table}}.}
#'   \item{Parameters}{Parameters on the natural scale, with the rows of
#'     their variance-covariance matrix in the columns \code{cov_1},
#'     \code{cov_2}, ... (column k is the covariance with the k-th parameter
#'     of the model).}
#'   \item{Parameters (transformed)}{The same on the scale on which the
#'     parameters are estimated (for probabilistic sensitivity analysis).}
#'   \item{Covariance, Cholesky}{\code{\link{covariance_table}} on both
#'     scales.}
#'   \item{Survival summary}{\code{\link{survival_summary}} with the
#'     restricted mean up to \code{time_horizon}.}
#'   \item{Landmark survival}{\code{\link{predict_survival}} at
#'     \code{0, interval, 2 * interval, ...} up to \code{time_horizon} and
#'     at the landmark times.}
#'   \item{Spline knots}{\code{\link{spline_knot_table}} (spline models
#'     only).}
#' }
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param file Path of the \code{.xlsx} file to write (replaced if it
#'   exists).
#' @param interval Spacing of the reporting times of the Kaplan-Meier table
#'   and of the predictions. The default \code{NULL} uses the spacing of the
#'   time axis of the plots of the data (for example 12 months).
#' @param landmarks Times at which the survival is predicted. The default
#'   \code{NULL} uses 6 months and 1, 2, 5 and
#'   10 years, in the time unit of the data.
#' @param time_horizon Time horizon of the predictions and of the
#'   restricted mean. The default \code{NULL} uses the last landmark time.
#' @param hazard_ratios Data frame returned by
#'   \code{\link{estimate_hazard_ratio}} (for example with the covariates
#'   and strata chosen for the analysis). The default \code{NULL} uses the
#'   Cox model of the arm without adjustment.
#' @param conf_int Logical; \code{TRUE} (the default) adds the simulated
#'   confidence intervals of the predictions and of the median, mean and
#'   restricted mean, which takes longer.
#' @param n_sim,seed Number of simulations and seed of the confidence
#'   intervals (see \code{\link{predict_survival}}).
#'
#' @return The path \code{file}, invisibly.
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = NULL)
#' path <- tempfile(fileext = ".xlsx")
#' write_survival_workbook(fits, path, conf_int = FALSE)
#' unlink(path)
#'
#' @export
write_survival_workbook <- function(x,
                                    file,
                                    interval = NULL,
                                    landmarks = NULL,
                                    time_horizon = NULL,
                                    hazard_ratios = NULL,
                                    conf_int = TRUE,
                                    n_sim = 1000,
                                    seed = 1) {
  .check_survfitsummary(x)
  .check_output_file(file, "xlsx")
  if (!is.logical(conf_int) || length(conf_int) != 1L || is.na(conf_int)) {
    stop("'conf_int' must be TRUE or FALSE.", call. = FALSE)
  }
  s <- .report_settings(x, interval, landmarks, time_horizon,
                        hazard_ratios)
  times <- sort(unique(c(seq(0, s$time_horizon, by = s$interval),
                         s$landmarks$time)))
  versions <- vapply(c("SurvFitSummary", "flexsurv", "survival",
                       "bshazard"),
                     function(p) as.character(utils::packageVersion(p)),
                     character(1))
  settings <- data.frame(
    item = c("Created", "Time unit", "Confidence level",
             "Reporting interval", "Landmark times", "Time horizon",
             "Confidence intervals of the predictions",
             "Simulations of the confidence intervals", "Seed",
             "R version", paste(names(versions), "version")),
    value = c(format(Sys.Date()), s$unit, format(x$conf_level),
              format(s$interval),
              paste(format(s$landmarks$time, trim = TRUE),
                    collapse = ", "),
              format(s$time_horizon), if (conf_int) "Yes" else "No",
              format(n_sim), format(seed),
              paste(R.version$major, R.version$minor, sep = "."),
              unname(versions)),
    stringsAsFactors = FALSE
  )
  two_arms <- length(x$arms) >= 2L
  status <- x$models[c("fit_type", "arm", "model", "family", "distribution",
                       "scale", "n_knots", "n", "events", "status",
                       "message")]
  sheets <- list(
    "Settings" = settings,
    "Data summary" = summarise_data(x$data, conf_level = x$conf_level),
    "KM table" = km_table(x$data, interval = s$interval,
                          conf_level = x$conf_level),
    "Hazard ratios" = s$hazard_ratios,
    "PH test" = if (two_arms) test_proportional_hazards(x$data),
    "Model status" = status,
    "Fit statistics" = gof_table(x),
    "Parameters" = .parameter_sheet(x, "natural"),
    "Parameters (transformed)" = .parameter_sheet(x, "transformed"),
    "Covariance" = covariance_table(x, type = "covariance",
                                    scale = c("transformed", "natural")),
    "Cholesky" = covariance_table(x, type = "cholesky",
                                  scale = c("transformed", "natural")),
    "Survival summary" = survival_summary(x, time_horizon = s$time_horizon,
                                          conf_int = conf_int,
                                          n_sim = n_sim, seed = seed),
    "Landmark survival" = predict_survival(x, times, conf_int = conf_int,
                                           n_sim = n_sim, seed = seed),
    "Spline knots" = spline_knot_table(x)
  )
  keep <- vapply(sheets, function(d) !is.null(d) && nrow(d) > 0L,
                 logical(1))
  keep[["Model status"]] <- TRUE
  sheets <- lapply(sheets[keep], function(d) {
    d[] <- lapply(d, function(v) if (is.factor(v)) as.character(v) else v)
    d
  })
  writexl::write_xlsx(sheets, path = file, col_names = TRUE,
                      format_headers = TRUE)
  invisible(file)
}
